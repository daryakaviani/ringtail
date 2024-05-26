package main

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"fmt"
	"lattice-threshold-signature/utils"
	"log"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
	"golang.org/x/crypto/sha3"
)

// PARAMETERS
const (
	m         = 10 // TESTING WITH SMALL PARAMS, TODO: UPDATE
	n         = 9  // TESTING WITH SMALL PARAMS, TODO: UPDATE
	d         = 91 // TESTING WITH SMALL PARAMS, TODO: UPDATE
	dbar      = d - 1
	p         = 1 << 30
	t         = 1 // Active threshold, TESTING WITH SMALL PARAMS, TODO: UPDATE
	k         = 1 // Total number of parties, TESTING WITH SMALL PARAMS, TODO: UPDATE
	ell       = 5
	B         = 655327113807872.217683738776675 // 2^49.219208552119575
	kappa     = 23
	logN      = 8
	sigma_e   = 3.6649122421707636
	bound_e   = sigma_e * 2
	sigmaE    = sigma_e
	boundE    = sigmaE * 2
	sigmaStar = 1022622718524.366463783348896 // 2^39.8954111209675
	boundStar = sigmaStar * 2
	sigmaU    = 3808175.1950851358
	boundU    = sigmaU * 2
	keySize   = 30
	q         = 0x8000000001c01 // 51-bit NTT-friendly prime
)

// Party struct holds all state and methods for a party in the protocol
type Party struct {
	ID             int
	Ring           *ring.Ring
	UniformSampler *ring.UniformSampler
	SkShare        structs.Vector[ring.Poly]
	Seed           [][]byte
	Mask           structs.Vector[ring.Poly]
	R              structs.Matrix[ring.Poly]
	C              ring.Poly
	H              structs.Vector[[]*big.Int]
	Lambda         ring.Poly
	D              structs.Matrix[ring.Poly]
}

// NewParty initializes a new Party instance
func NewParty(id int, r *ring.Ring, sampler *ring.UniformSampler) *Party {
	return &Party{
		ID:             id,
		Ring:           r,
		UniformSampler: sampler,
	}
}

// main function orchestrates the threshold signature protocol
func main() {
	var setupDuration, genDuration, signRound1Duration, signRound2Duration, finalizeDuration, verifyDuration time.Duration

	randomKey := make([]byte, keySize)

	r, _ := ring.NewRing(1<<logN, []uint64{q})

	prng, _ := sampling.NewKeyedPRNG(randomKey)
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := randomKey

	// SETUP phase
	start := time.Now()
	A := Setup(uniformSampler)
	utils.PrintMatrix("Matrix A:", A)
	setupDuration = time.Since(start)

	parties := make([]*Party, k)
	for i := range parties {
		parties[i] = NewParty(i, r, uniformSampler)
	}

	// GEN: Generate secret shares and seeds
	start = time.Now()
	skShares, seeds, b := Gen(r, A, uniformSampler, trustedDealerKey)
	utils.PrintVector("Vector b:", b)
	genDuration = time.Since(start)

	for partyID := 0; partyID < k; partyID++ {
		parties[partyID].SkShare = skShares[partyID]
		parties[partyID].Seed = seeds[partyID]
		utils.PrintVector(fmt.Sprintf("Secret share for party %d:", partyID), skShares[partyID])
	}

	// SIGNATURE ROUND 1
	mu := "Message"
	sid := 1
	PRFKey := "PRF Key" // PRF key with fixed randomness for testing, TODO: Update
	T := []int{0}       // Active parties
	start = time.Now()
	lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(q)))
	D := make(map[int]structs.Matrix[ring.Poly])
	masks := make(map[int]structs.Vector[ring.Poly])
	for _, partyID := range T {
		parties[partyID].Lambda = lagrangeCoeffs[partyID]
		parties[partyID].Seed = seeds[partyID]
		D[partyID], masks[partyID] = parties[partyID].SignRound1(A, sid, []byte(PRFKey), T)
		utils.PrintMatrix(fmt.Sprintf("Matrix D for party %d:", partyID), D[partyID])
		utils.PrintVector(fmt.Sprintf("Mask for party %d:", partyID), masks[partyID])
	}
	signRound1Duration = time.Since(start)

	// SIGNATURE ROUND 2
	start = time.Now()
	z := make(map[int]structs.Vector[ring.Poly])
	for _, partyID := range T {

		// Create new maps excluding the party's own D and mask
		DExcludingParty := make(map[int]structs.Matrix[ring.Poly])
		masksExcludingParty := make(map[int]structs.Vector[ring.Poly])
		for key, value := range D {
			if key != partyID {
				DExcludingParty[key] = value
			}
		}
		for key, value := range masks {
			if key != partyID {
				masksExcludingParty[key] = value
			}
		}

		z[partyID] = parties[partyID].SignRound2(A, b, DExcludingParty, masksExcludingParty, sid, mu, T, []byte(PRFKey), seeds)
		utils.PrintVector(fmt.Sprintf("Vector z for party %d:", partyID), z[partyID])
	}
	signRound2Duration = time.Since(start)

	// SIGNATURE FINALIZE
	start = time.Now()
	finalParty := parties[0]
	_, sig, Delta := finalParty.SignFinalize(z, masks, A, b)
	finalizeDuration = time.Since(start)

	// Verify the signature
	betaDelta := calculateBetaDelta()
	start = time.Now()
	valid := Verify(r, sig, A, mu, b, finalParty.C, Delta, betaDelta)
	verifyDuration = time.Since(start)
	fmt.Printf("Signature Verification Result: %v\n", valid)

	// Print all durations
	fmt.Println("Setup duration:", setupDuration)
	fmt.Println("Gen duration:", genDuration)
	fmt.Println("Signature Round 1 duration:", signRound1Duration)
	fmt.Println("Signature Round 2 duration:", signRound2Duration)
	fmt.Println("Finalize duration:", finalizeDuration)
	fmt.Println("Verify duration:", verifyDuration)
}

// Setup generates the public parameters
func Setup(uniformSampler *ring.UniformSampler) structs.Matrix[ring.Poly] {
	matrix := utils.SamplePolyMatrix(m, n, uniformSampler)
	utils.PrintMatrix("Setup Matrix:", matrix)
	return matrix
}

// Gen generates the secret shares, seeds, and the public parameter b
func Gen(r *ring.Ring, A structs.Matrix[ring.Poly], uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (map[int]structs.Vector[ring.Poly], map[int][][]byte, structs.Vector[ring.Poly]) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: bound_e}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := utils.SamplePolyVector(n, uniformSampler)
	utils.PrintVector("Secret Vector s:", s)
	skShares := ShamirSecretSharing(r, s, t, k)
	e := utils.SamplePolyVector(m, gaussianSampler)
	utils.PrintVector("Error Vector e:", e)
	b := utils.InitializeVector(r, m)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	seeds := make(map[int][][]byte)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			seeds[i][j] = generateRandomSeed()
		}
	}

	return skShares, seeds, b
}

// generateRandomSeed generates a random seed of length ell
func generateRandomSeed() []byte {
	sd := make([]byte, ell)
	_, err := rand.Read(sd)
	if err != nil {
		log.Fatalf("Error generating random seed %v\n", err)
	}
	return sd
}

// SignRound1 performs the first round of signing
func (party *Party) SignRound1(A structs.Matrix[ring.Poly], sid int, PRFKey []byte, T []int) (structs.Matrix[ring.Poly], structs.Vector[ring.Poly]) {
	r := party.Ring
	uniformSampler := party.UniformSampler
	seeds := party.Seed

	// Initialize the mask vector using the helper function
	mask := utils.InitializeVector(r, n)

	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j], PRFKey)
		utils.VectorAdd(r, mask, mask_j, mask)
	}

	// Initialize r_star
	r_star := utils.SamplePolyVector(n, uniformSampler)

	// Initialize e_star
	skHash := party.PRNGKey()
	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaStar, Bound: boundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)
	e_star := utils.SamplePolyVector(m, gaussianSampler)

	// Initialize R_i
	R_i := utils.SamplePolyMatrix(n, dbar, uniformSampler)

	// Initialize E_i
	gaussianParams = ring.DiscreteGaussian{Sigma: sigmaE, Bound: boundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	E_i := utils.SamplePolyMatrix(m, dbar, gaussianSampler)

	// Ensure concatenatedR is properly initialized
	concatenatedR := utils.InitializeMatrix(r, n, d)
	for i := range concatenatedR {
		concatenatedR[i] = append([]ring.Poly{r_star[i]}, R_i[i]...)
	}
	party.R = concatenatedR

	// Ensure concatenatedE is properly initialized
	concatenatedE := utils.InitializeMatrix(r, m, d)
	for i := range concatenatedE {
		concatenatedE[i] = append([]ring.Poly{e_star[i]}, E_i[i]...)
	}

	// utils.PrintMatrix("Try:", concatenatedE)

	utils.PrintMatrix("ConcatE", concatenatedE)
	utils.PrintMatrix("ConcatR", concatenatedR)

	D := utils.InitializeMatrix(r, m, n)
	utils.MatrixMatrixMul(r, A, concatenatedR, D)
	utils.MatrixAdd(r, concatenatedE, D, D)

	utils.PrintMatrix("Matrix D:", D)
	utils.PrintVector("Mask vector:", mask)

	party.Mask = mask
	party.D = D

	return D, mask
}

// SignRound2 performs the second round of signing
func (party *Party) SignRound2(A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], DExcludingParty map[int]structs.Matrix[ring.Poly], masksExcludingParty map[int]structs.Vector[ring.Poly], sid int, mu string, T []int, PRFKey []byte, seeds map[int][][]byte) structs.Vector[ring.Poly] {
	r := party.Ring
	partyID := party.ID
	concatR := party.R
	s_i := party.SkShare
	lambda := party.Lambda

	u := make(map[int]structs.Vector[ring.Poly], d)
	onePoly := r.NewMonomialXi(0)

	D := utils.CopyMatrixMap(DExcludingParty)
	D[partyID] = party.D
	masks := utils.CopyVectorMap(masksExcludingParty)
	masks[partyID] = party.Mask

	for _, j := range T {
		oneSlice := structs.Vector[ring.Poly]{onePoly}
		u[j] = oneSlice
		if d > 1 {
			u_j := H_u(r, A, b, sid, T, j, D, mu, masks)
			u[j] = append(oneSlice, u_j...)
		}
	}

	h := utils.InitializeVector(r, m)

	for j, D_j := range D {
		D_j_u_j := utils.InitializeVector(r, m)
		utils.MatrixVectorMul(r, D_j, u[j], D_j_u_j)
		utils.VectorAdd(r, h, D_j_u_j, h)
	}

	roundedH := utils.RoundCoeffsToNearestMultiple(r, h, p, q)
	party.H = roundedH

	c := H_c(r, A, b, roundedH, mu)
	party.C = c

	mPrime := utils.InitializeVector(r, n)

	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j][partyID], PRFKey)
		utils.VectorAdd(r, mPrime, mask_j, mPrime)
	}

	utils.PrintVector("Mask prime vector:", mPrime)

	z_i := utils.InitializeVector(r, n)
	utils.MatrixVectorMul(r, concatR, u[partyID], z_i)
	utils.VectorAdd(r, z_i, mPrime, z_i)
	s_c_lambda := utils.InitializeVector(r, n)
	utils.VectorPolyMul(r, s_i, c, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, lambda, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	utils.PrintVector(fmt.Sprintf("Vector z for party %d:", partyID), z_i)

	return z_i
}

// SignFinalize finalizes the signature
func (party *Party) SignFinalize(z map[int]structs.Vector[ring.Poly], masks map[int]structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly]) (ring.Poly, structs.Vector[ring.Poly], structs.Vector[[]*big.Int]) {
	r := party.Ring
	c := party.C
	h := party.H

	z_sum := utils.InitializeVector(r, n)

	for j, z_j := range z {
		mask_j := masks[j]
		for i := 0; i < n; i++ {
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

	utils.PrintVector("z_sum", z_sum)

	Az := utils.InitializeVector(r, m)
	utils.MatrixVectorMul(r, A, z_sum, Az)

	utils.PrintVector("Az", Az)

	bc := utils.InitializeVector(r, m)
	utils.VectorPolyMul(r, b, c, bc)

	utils.PrintVector("bc", bc)

	Az_bc := utils.InitializeVector(r, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	utils.PrintVector("Az-bc before rounding", Az_bc)

	Az_bc_bigint := utils.RoundCoeffsToNearestMultiple(r, Az_bc, p, q)
	utils.PrintBigIntVector("Az-bc after rounding", Az_bc_bigint)

	Delta := make(structs.Vector[[]*big.Int], m)
	for i := range Delta {
		Delta[i] = make([]*big.Int, r.N())
		for j := range Delta[i] {
			Delta[i][j] = new(big.Int).Sub(h[i][j], Az_bc_bigint[i][j])
		}
	}

	utils.PrintBigIntVector("Vector Delta:", Delta)

	return party.C, z_sum, Delta
}

// Verify verifies the correctness of the signature
func Verify(r *ring.Ring, z structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], mu string, b structs.Vector[ring.Poly], c ring.Poly, Delta structs.Vector[[]*big.Int], betaDelta *big.Int) bool {
	Az := utils.InitializeVector(r, m)
	utils.MatrixVectorMul(r, A, z, Az)

	bc := utils.InitializeVector(r, m)
	utils.VectorPolyMul(r, b, c, bc)

	Az_bc := utils.InitializeVector(r, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	Az_bc_bigint := utils.RoundCoeffsToNearestMultiple(r, Az_bc, p, q)

	Az_bc_Delta := make(structs.Vector[[]*big.Int], m)
	for i := range Az_bc_Delta {
		Az_bc_Delta[i] = make([]*big.Int, r.N())
		for j := range Az_bc_Delta[i] {
			Az_bc_Delta[i][j] = new(big.Int).Add(Az_bc_bigint[i][j], Delta[i][j])
		}
	}

	//utils.PrintBigIntVector("Az_bc_Delta:", Az_bc_Delta)

	computedC := H_c(r, A, b, Az_bc_Delta, mu)
	if !r.Equal(c, computedC) {
		return false
	}

	return CheckInfinityNorm(Delta, betaDelta, new(big.Int).SetUint64(p))
}

// CheckInfinityNorm checks if the infinity norm of the vector of *big.Int Delta is less than or equal to betaDelta
func CheckInfinityNorm(Delta structs.Vector[[]*big.Int], betaDelta *big.Int, modulus *big.Int) bool {
	maxValue := big.NewInt(0)

	for _, polyCoeffs := range Delta {
		for _, coeff := range polyCoeffs {
			absCoeff := new(big.Int).Abs(coeff)
			p_absCoeff := new(big.Int).Sub(modulus, absCoeff)
			infNorm := big.NewInt(0)

			if absCoeff.Cmp(modulus) > 0 {
				return false
			}
			// |a|_inf := min{a, p-a} for a in Z_p
			if absCoeff.Cmp(p_absCoeff) > 0 {
				infNorm = p_absCoeff
			} else {
				infNorm = absCoeff
			}
			if infNorm.Cmp(maxValue) > 0 {
				maxValue = infNorm
			}
		}
	}

	log.Println("Delta Norm:", maxValue)
	log.Println("Beta Delta:", betaDelta)

	return maxValue.Cmp(betaDelta) <= 0
}

// calculateBetaDelta computes ((B * p) - q) / (2 * q) as a big.Int
func calculateBetaDelta() *big.Int {
	// Convert constants to big.Int
	pInt := new(big.Int).SetUint64(p)
	BInt := new(big.Int)
	BInt.SetString(fmt.Sprintf("%.0f", B), 10) // Convert float64 to string, then to big.Int
	qInt := new(big.Int).SetUint64(q)

	// B * p
	Bp := new(big.Int).Mul(BInt, pInt)
	// B * p - q
	BpSubQ := new(big.Int).Sub(Bp, qInt)
	// 2 * q
	twoQ := new(big.Int).Mul(big.NewInt(2), qInt)
	// ((B * p) - q) / (2 * q)
	betaDelta := new(big.Int).Div(BpSubQ, twoQ)
	return betaDelta
}

// HASHES & PRF

// H_u hashes parameters to a Gaussian distribution
func H_u(r *ring.Ring, A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], sid int, T []int, j int, D map[int]structs.Matrix[ring.Poly], mu string, masks map[int]structs.Vector[ring.Poly]) structs.Vector[ring.Poly] {
	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	binary.Write(buf, binary.BigEndian, int64(sid))
	binary.Write(buf, binary.BigEndian, T)
	binary.Write(buf, binary.BigEndian, int64(j))

	for _, D_h := range D {
		if _, err := D_h.WriteTo(buf); err != nil {
			log.Fatalf("Error writing matrix D_h: %v\n", err)
		}
	}

	buf.WriteString(mu)

	for _, mask := range masks {
		if _, err := mask.WriteTo(buf); err != nil {
			log.Fatalf("Error writing vector mask: %v\n", err)
		}
	}

	_, err := hasher.Write(buf.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	//gaussianParams := ring.DiscreteGaussian{}
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaU, Bound: boundU}
	hashGaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	u_j := utils.SamplePolyVector(dbar, hashGaussianSampler)

	return u_j
}

// PRF generates pseudorandom ring elements
func PRF(r *ring.Ring, sid int, sd_ij []byte, PRFKey []byte) structs.Vector[ring.Poly] {
	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)

	binary.Write(buf, binary.BigEndian, PRFKey)
	binary.Write(buf, binary.BigEndian, sd_ij)
	binary.Write(buf, binary.BigEndian, int64(sid))

	_, err := hasher.Write(buf.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	PRFUniformSampler := ring.NewUniformSampler(prng, r)
	mask := utils.SamplePolyVector(n, PRFUniformSampler)
	return mask
}

// H_c hashes to low norm ring elements
func H_c(r *ring.Ring, A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], h structs.Vector[[]*big.Int], mu string) ring.Poly {
	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	for _, coeffs := range h {
		for _, coeff := range coeffs {
			binary.Write(buf, binary.BigEndian, coeff.Bytes())
		}
	}

	binary.Write(buf, binary.BigEndian, []byte(mu))

	_, err := hasher.Write(buf.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	ternaryParams := ring.Ternary{H: kappa}
	ternarySampler, err := ring.NewTernarySampler(prng, r, ternaryParams, false)
	if err != nil {
		log.Fatalf("Error creating ternary sampler: %v", err)
	}
	c := ternarySampler.ReadNew()

	utils.PrintPolynomial("Polynomial c:", c)

	return c
}

// ShamirSecretSharing shares each coefficient of a vector of ring.Poly across k parties using (t, k)-threshold Shamir secret sharing.
func ShamirSecretSharing(r *ring.Ring, s structs.Vector[ring.Poly], t, k int) map[int]structs.Vector[ring.Poly] {
	degree := r.N()

	shares := make(map[int]structs.Vector[ring.Poly], k)
	for i := 0; i < k; i++ {
		shares[i] = structs.Vector[ring.Poly](make([]ring.Poly, len(s)))
		for j := range shares[i] {
			newPoly := r.NewPoly()
			shares[i][j] = newPoly
		}
	}

	for polyIndex, poly := range s {
		coeffs := make([]*big.Int, degree)
		r.PolyToBigint(poly, 1, coeffs)

		for coeffIndex, secret := range coeffs {
			polyCoeffs := make([]*big.Int, t)
			polyCoeffs[0] = secret
			for i := 1; i < t; i++ {
				randomCoeff := big.NewInt(15) // TODO: Fixed randomness for testing, make it random later
				polyCoeffs[i] = randomCoeff
			}

			for i := 1; i <= k; i++ {
				x := big.NewInt(int64(i))
				shareValue := big.NewInt(0)
				xPow := big.NewInt(1)

				for _, coeff := range polyCoeffs {
					term := new(big.Int).Mul(coeff, xPow)
					shareValue.Add(shareValue, term)
					shareValue.Mod(shareValue, big.NewInt(int64(q)))
					xPow.Mul(xPow, x)
				}

				if shares[i-1][polyIndex].Coeffs[0] == nil {
					shares[i-1][polyIndex].Coeffs[0] = make([]uint64, degree)
				}

				shares[i-1][polyIndex].Coeffs[0][coeffIndex] = shareValue.Uint64()
			}
		}
	}

	return shares
}

// ComputeLagrangeCoefficients computes the Lagrange coefficients for interpolation based on the indices of available shares.
func ComputeLagrangeCoefficients(r *ring.Ring, T []int, modulus *big.Int) map[int]ring.Poly {
	lagrangeCoefficients := make(map[int]ring.Poly)
	for i := 0; i < len(T); i++ {
		xi := big.NewInt(int64(T[i] + 1))
		numerator := big.NewInt(1)
		denominator := big.NewInt(1)
		for j := 0; j < len(T); j++ {
			if i != j {
				xj := big.NewInt(int64(T[j] + 1))
				numerator.Mul(numerator, new(big.Int).Neg(xj))
				numerator.Mod(numerator, modulus)
				temp := new(big.Int).Sub(xi, xj)
				denominator.Mul(denominator, temp)
				denominator.Mod(denominator, modulus)
			}
		}
		denomInv := new(big.Int).ModInverse(denominator, modulus)
		coeff := new(big.Int).Mul(numerator, denomInv)
		coeff.Mod(coeff, modulus)
		lagrangePoly := r.NewPoly()
		r.SetCoefficientsBigint([]*big.Int{coeff}, lagrangePoly)
		lagrangeCoefficients[T[i]] = lagrangePoly
		utils.PrintPolynomial(fmt.Sprintf("Lagrange Coefficient for party %d:", T[i]), lagrangePoly)
	}
	return lagrangeCoefficients
}

// PRNGKey generates a key for PRNG using the secret key share
func (party *Party) PRNGKey() []byte {
	skShare := party.SkShare

	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)
	for _, poly := range skShare {
		if _, err := poly.WriteTo(buf); err != nil {
			log.Fatalf("Error writing poly to buffer: %v\n", err)
		}
	}

	hasher.Write(buf.Bytes())

	hashOutputLength := keySize
	skHash := make([]byte, hashOutputLength)
	_, err := hasher.Read(skHash)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}
	return skHash
}
