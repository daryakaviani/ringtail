package sign

import (
	"bytes"
	"fmt"
	"log"
	"math"
	"math/big"

	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/gonum/mat"

	"lattice-threshold-signature/primitives"
	"lattice-threshold-signature/utils"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

// Party struct holds all state and methods for a party in the protocol
type Party struct {
	ID             int
	Ring           *ring.Ring
	RingXi         *ring.Ring
	RingNu         *ring.Ring
	UniformSampler *ring.UniformSampler
	SkShare        structs.Vector[ring.Poly]
	Seed           map[int][][]byte
	R              structs.Matrix[ring.Poly]
	C              ring.Poly
	H              structs.Vector[ring.Poly]
	Lambda         ring.Poly
	D              structs.Matrix[ring.Poly]
	MACKeys        map[int][]byte // Pairwise MAC keys
	MACs           map[int][]byte // MACs generated for each party
}

// NewParty initializes a new Party instance
func NewParty(id int, r *ring.Ring, r_xi *ring.Ring, r_nu *ring.Ring, sampler *ring.UniformSampler) *Party {
	return &Party{
		ID:             id,
		Ring:           r,
		RingXi:         r_xi,
		RingNu:         r_nu,
		UniformSampler: sampler,
		MACKeys:        make(map[int][]byte),
		MACs:           make(map[int][]byte),
	}
}

// Gen generates the secret shares, seeds, MAC keys, and the public parameter b
func Gen(r *ring.Ring, r_xi *ring.Ring, uniformSampler *ring.UniformSampler, trustedDealerKey []byte, lagrangeCoefficients structs.Vector[ring.Poly]) (structs.Matrix[ring.Poly], map[int]structs.Vector[ring.Poly], map[int][][]byte, map[int]map[int][]byte, structs.Vector[ring.Poly]) {
	A := utils.SamplePolyMatrix(r, M, N, uniformSampler, true, true)

	precomputeSize := (K * K * KeySize) + (r.N() * N * (K - 1) * len(r.Modulus().Bytes())) + (K * (K - 1) * KeySize)
	utils.PrecomputeRandomness(precomputeSize, trustedDealerKey)

	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: SigmaE, Bound: BoundE}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := utils.SamplePolyVector(r, N, gaussianSampler, false, false)
	skShares := primitives.ShamirSecretSharing(r, s, Threshold, lagrangeCoefficients)

	for _, skShare := range skShares {
		utils.ConvertVectorToNTT(r, skShare)
	}
	utils.ConvertVectorToNTT(r, s)

	e := utils.SamplePolyVector(r, M, gaussianSampler, true, true)
	b := utils.InitializeVector(r, M)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	// Round b
	utils.ConvertVectorFromNTT(r, b)
	bTilde := utils.RoundVector(r, r_xi, b, Xi)

	seeds := make(map[int][][]byte)
	MACKeys := make(map[int]map[int][]byte)
	MACKeys[0] = make(map[int][]byte)

	for i := 0; i < K; i++ {
		seeds[i] = make([][]byte, K)
		for j := 0; j < K; j++ {
			seeds[i][j] = utils.GetRandomBytes(KeySize)
			if i != j {
				if MACKeys[j] == nil {
					MACKeys[j] = make(map[int][]byte)
				}
				if MACKeys[i][j] == nil && MACKeys[j][i] == nil {
					MACKeys[i][j] = utils.GetRandomBytes(KeySize)
					MACKeys[j][i] = MACKeys[i][j]
				}
			}
		}
	}

	return A, skShares, seeds, MACKeys, bTilde
}

// SignRound1 performs the first round of signing
func (party *Party) SignRound1(A structs.Matrix[ring.Poly], sid int, PRFKey []byte, T []int) (structs.Matrix[ring.Poly], map[int][]byte) {
	r := party.Ring

	// Initialize r_star and e_star
	skHash := primitives.PRNGKey(party.SkShare)
	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: SigmaStar, Bound: BoundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)
	r_star := utils.SamplePolyVector(r, N, gaussianSampler, true, true)
	e_star := utils.SamplePolyVector(r, M, gaussianSampler, true, true)

	// Initialize R_i and E_i
	gaussianParams = ring.DiscreteGaussian{Sigma: SigmaE, Bound: BoundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	R_i := utils.SamplePolyMatrix(r, N, Dbar, gaussianSampler, true, true)
	E_i := utils.SamplePolyMatrix(r, M, Dbar, gaussianSampler, true, true)

	concatenatedR := utils.InitializeMatrix(r, N, Dbar+1)
	for i := range concatenatedR {
		concatenatedR[i] = append([]ring.Poly{r_star[i]}, R_i[i]...)
	}
	party.R = concatenatedR

	// Ensure concatenatedE is properly initialized
	concatenatedE := utils.InitializeMatrix(r, M, Dbar+1)
	for i := range concatenatedE {
		concatenatedE[i] = append([]ring.Poly{e_star[i]}, E_i[i]...)
	}

	D := utils.InitializeMatrix(r, M, Dbar+1)

	utils.MatrixMatrixMul(r, A, concatenatedR, D)
	utils.MatrixAdd(r, concatenatedE, D, D)

	party.D = D

	// Generate MACs for each party
	MACs := make(map[int][]byte)
	for _, j := range T {
		if j != party.ID {
			MACs[j] = primitives.GenerateMAC(D, party.MACKeys[j], party.ID, sid, T, j, false)
		}
	}

	return D, MACs
}

// SignRound2Preprocess verifies the MACs received in round 1 and performs the minimum eigenvalue check
func (party *Party) SignRound2Preprocess(A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], D map[int]structs.Matrix[ring.Poly], MACs map[int]map[int][]byte, sid int, T []int) (bool, structs.Matrix[ring.Poly], []byte) {
	hash := primitives.Hash(A, b, D, sid, T)

	for _, j := range T {
		if j != party.ID {
			MAC := MACs[j][party.ID]
			expectedMAC := primitives.GenerateMAC(D[j], party.MACKeys[j], party.ID, sid, T, j, true)
			if !bytes.Equal(MAC, expectedMAC) {
				return false, nil, nil
			}
		}
	}

	DSum := utils.InitializeMatrix(party.Ring, M, Dbar+1)
	for _, D_j := range D {
		utils.MatrixAdd(party.Ring, D_j, DSum, DSum)
	}

	DBar := make(structs.Matrix[ring.Poly], len(DSum))
	for i := range DSum {
		DBar[i] = make([]ring.Poly, len(DSum[i])-1)
		for j := 1; j < len(DSum[i]); j++ {
			DBar[i][j-1] = *DSum[i][j].CopyNew()
		}
	}

	utils.ConvertMatrixFromNTT(party.Ring, DBar)

	if !CheckMinEigenvalue(party.Ring, DBar) {
		log.Fatalf("Minimum eigenvalue check failed! Aborting.")
		return false, nil, nil
	}

	if !CheckCoprimeEntries(party.Ring, DBar) {
		log.Fatalf("Coprime check failed! Aborting.")
		return false, nil, nil
	}

	return true, DSum, hash
}

// SignRound2 performs the second round of signing
func (party *Party) SignRound2(A structs.Matrix[ring.Poly], bTilde structs.Vector[ring.Poly], DSum structs.Matrix[ring.Poly], sid int, mu string, T []int, PRFKey []byte, hash []byte) (structs.Vector[ring.Poly], structs.Vector[ring.Poly]) {
	r := party.Ring
	r_nu := party.RingNu
	partyID := party.ID
	concatR := party.R
	seeds := party.Seed

	s_i := party.SkShare
	lambda := party.Lambda

	onePoly := r.NewMonomialXi(0)
	r.NTT(onePoly, onePoly)
	r.MForm(onePoly, onePoly)

	u := structs.Vector[ring.Poly]{}
	oneSlice := structs.Vector[ring.Poly]{onePoly}
	if Dbar > 0 {
		h_u := primitives.GaussianHash(r, hash, mu, SigmaU, BoundU, Dbar)
		u = append(oneSlice, h_u...)
	}

	h := utils.InitializeVector(r, M)
	utils.MatrixVectorMul(r, DSum, u, h)

	utils.ConvertVectorFromNTT(r, h)
	roundedH := utils.RoundVector(r, r_nu, h, Nu)
	party.H = roundedH

	c := primitives.LowNormHash(r, A, bTilde, roundedH, mu, Kappa)
	party.C = c

	seed_i := party.Seed[party.ID]
	mask := utils.InitializeVector(r, N)
	for _, j := range T {
		mask_j := primitives.PRF(r, seed_i[j], PRFKey, mu, hash, N)
		utils.VectorAdd(r, mask, mask_j, mask)
	}

	maskPrime := utils.InitializeVector(r, N)
	for _, j := range T {
		mask_j := primitives.PRF(r, seeds[j][partyID], PRFKey, mu, hash, N)
		utils.VectorAdd(r, maskPrime, mask_j, maskPrime)
	}

	z_i := utils.InitializeVector(r, N)

	utils.MatrixVectorMul(r, concatR, u, z_i)

	utils.VectorAdd(r, z_i, maskPrime, z_i)

	s_c_lambda := utils.InitializeVector(r, N)

	utils.VectorPolyMul(r, s_i, lambda, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, c, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	return z_i, mask
}

// SignFinalize finalizes the signature
func (party *Party) SignFinalize(z map[int]structs.Vector[ring.Poly], masks map[int]structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], bTilde structs.Vector[ring.Poly]) (ring.Poly, structs.Vector[ring.Poly], structs.Vector[ring.Poly]) {
	r := party.Ring
	r_xi := party.RingXi
	r_nu := party.RingNu
	c := party.C
	h := party.H

	z_sum := utils.InitializeVector(r, N)

	for j, z_j := range z {
		mask_j := masks[j]
		for i := 0; i < N; i++ {
			if len(z_j) <= i {
				fmt.Printf("Error: z_j[%d] is out of range\n", i)
			}
			if len(mask_j) <= i {
				fmt.Printf("Error: mask_j[%d] is out of range\n", i)
			}
			r.Sub(z_j[i], mask_j[i], z_j[i])
			r.Add(z_sum[i], z_j[i], z_sum[i])
		}
	}

	Az_bc := utils.InitializeVector(r, M)
	utils.MatrixVectorMul(r, A, z_sum, Az_bc)
	bc := utils.InitializeVector(r, M)

	b := utils.RestoreVector(r, r_xi, bTilde, Xi)
	utils.ConvertVectorToNTT(r, b)

	utils.VectorPolyMul(r, b, c, bc)
	utils.VectorSub(r, Az_bc, bc, Az_bc)

	utils.ConvertVectorFromNTT(r, Az_bc)
	roundedAz_bc := utils.RoundVector(r, r_nu, Az_bc, Nu)

	Delta := utils.InitializeVector(r_nu, M)
	utils.VectorSub(r_nu, h, roundedAz_bc, Delta)

	return party.C, z_sum, Delta
}

// Verify verifies the correctness of the signature
func Verify(r *ring.Ring, r_xi *ring.Ring, r_nu *ring.Ring, z structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], mu string, bTilde structs.Vector[ring.Poly], c ring.Poly, roundedDelta structs.Vector[ring.Poly]) bool {
	Az_bc := utils.InitializeVector(r, M)
	utils.MatrixVectorMul(r, A, z, Az_bc)
	bc := utils.InitializeVector(r, M)

	b := utils.RestoreVector(r, r_xi, bTilde, Xi)
	utils.ConvertVectorToNTT(r, b)

	utils.VectorPolyMul(r, b, c, bc)
	utils.VectorSub(r, Az_bc, bc, Az_bc)

	utils.ConvertVectorFromNTT(r, Az_bc)
	roundedAz_bc := utils.RoundVector(r, r_nu, Az_bc, Nu)

	Az_bc_Delta := utils.InitializeVector(r_nu, M)
	utils.VectorAdd(r_nu, roundedAz_bc, roundedDelta, Az_bc_Delta)

	computedC := primitives.LowNormHash(r, A, bTilde, Az_bc_Delta, mu, Kappa)
	if !r.Equal(c, computedC) {
		return false
	}

	Delta := utils.RestoreVector(r, r_nu, roundedDelta, Nu)
	utils.ConvertVectorFromNTT(r, z)

	return CheckL2Norm(r, Delta, z)
}

// CheckL2Norm checks if the L2 norm of the vector of Delta is less than or equal to Bsquare
func CheckL2Norm(r *ring.Ring, Delta structs.Vector[ring.Poly], z structs.Vector[ring.Poly]) bool {
	sumSquares := big.NewInt(0)
	qBig := new(big.Int).SetUint64(Q)
	halfQ := new(big.Int).Div(qBig, big.NewInt(2))

	DeltaCoeffsBigInt := make(structs.Vector[[]*big.Int], r.N())
	for i, polyCoeffs := range Delta {
		DeltaCoeffsBigInt[i] = make([]*big.Int, r.N())
		r.PolyToBigint(polyCoeffs, 1, DeltaCoeffsBigInt[i])
	}

	for _, polyCoeffs := range DeltaCoeffsBigInt {
		for _, coeff := range polyCoeffs {
			if coeff.Cmp(halfQ) > 0 {
				coeff.Sub(coeff, qBig)
			}
			coeffSquare := new(big.Int).Mul(coeff, coeff)
			sumSquares.Add(sumSquares, coeffSquare)
		}
	}

	zCoeffsBigInt := make(structs.Vector[[]*big.Int], r.N())
	for i, polyCoeffs := range z {
		zCoeffsBigInt[i] = make([]*big.Int, r.N())
		r.PolyToBigint(polyCoeffs, 1, zCoeffsBigInt[i])
	}

	for _, polyCoeffs := range zCoeffsBigInt {
		for _, coeff := range polyCoeffs {
			if coeff.Cmp(halfQ) > 0 {
				coeff.Sub(coeff, qBig)
			}
			coeffSquare := new(big.Int).Mul(coeff, coeff)
			sumSquares.Add(sumSquares, coeffSquare)
		}
	}

	log.Println("Sum of Squares:", sumSquares)
	log.Println("Bsquare:", Bsquare)

	Bsquare, _ := new(big.Int).SetString(Bsquare, 10)
	return sumSquares.Cmp(Bsquare) <= 0
}

// Eigenvalue check on DSum matrix
func CheckMinEigenvalue(r *ring.Ring, DBar structs.Matrix[ring.Poly]) bool {
	phi := r.N()

	DBarFFT := make(structs.Matrix[[]complex128], len(DBar))
	for i := range DBar {
		DBarFFT[i] = make([][]complex128, len(DBar[i]))
		for j := range DBar[i] {
			poly := DBar[i][j]
			coeffs := make([]*big.Int, phi)
			r.PolyToBigint(poly, 1, coeffs)
			utils.SignedRepresentation(coeffs, Q)

			complexValues := make([]complex128, phi)
			for k, coeff := range coeffs {
				realPart, _ := new(big.Float).SetInt(coeff).Float64()
				complexValues[k] = complex(realPart, 0)
			}

			fft := fourier.NewCmplxFFT(len(complexValues))
			fftResult := fft.Coefficients(nil, complexValues)

			primitiveRootFFT := make([]complex128, phi/2)
			for k := 0; k < phi; k += 2 {
				index := k + 1
				if index < phi {
					primitiveRootFFT[k/2] = fftResult[index]
				}
			}
			DBarFFT[i][j] = primitiveRootFFT
		}
	}

	submatrices := make([]structs.Matrix[complex128], phi/2)
	for i := 0; i < phi/2; i++ {
		submatrix := make(structs.Matrix[complex128], len(DBarFFT))
		for j := range DBarFFT {
			submatrix[j] = make([]complex128, len(DBarFFT[j]))
			for k := range DBarFFT[j] {
				submatrix[j][k] = DBarFFT[j][k][i]
			}
		}
		submatrices[i] = submatrix
	}

	for i, submatrix := range submatrices {
		submatrices[i] = utils.MultiplyByConjugateTranspose(submatrix)
	}

	minEigenvalue := math.Inf(1)
	for _, submatrix := range submatrices {
		realMatrix := utils.ConvertComplexToRealMatrix(submatrix)

		var eig mat.Eigen
		ok := eig.Factorize(realMatrix, mat.EigenRight)
		if !ok {
			log.Fatal("Eigendecomposition failed")
		}

		eigenvalues := eig.Values(nil)
		for _, ev := range eigenvalues {
			realPart := real(ev)
			if realPart < minEigenvalue {
				minEigenvalue = realPart
			}
		}
	}

	sqrtMinEigenvalue := math.Sqrt(minEigenvalue) / float64(phi)
	EtaEpsilonQ := new(big.Float).Mul(big.NewFloat(EtaEpsilon), new(big.Float).SetInt64(Q))
	EtaEpsilonQSigmaU := new(big.Float).Quo(EtaEpsilonQ, new(big.Float).SetFloat64(SigmaU))
	EtaEpsilonQSigmaUFloat, _ := EtaEpsilonQSigmaU.Float64()

	return sqrtMinEigenvalue > EtaEpsilonQSigmaUFloat
}

// CheckCoprimeEntries checks if the entries of all rows in DBar are setwise coprime
func CheckCoprimeEntries(r *ring.Ring, DBar structs.Matrix[ring.Poly]) bool {
	for _, row := range DBar {
		if !utils.CheckRowCoprime(r, row) {
			return false
		}
	}
	return true
}
