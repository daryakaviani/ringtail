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
	m = 11
	n = 9
	d = 10 // Length of joint noise vector

	p         = 2 ^ 30
	t         = 64   // Active threshold
	k         = 1024 // Total number of parties
	ell       = 1
	beta      = 10
	betaDelta = 100000000000
	kappa     = 23
	logN      = 8
	bound_e   = 1000
	sigma_e   = 32
	boundE    = 1000
	sigmaE    = 32
	boundStar = 24296004000 * 2
	sigmaStar = 24296004000
	boundU    = 0
	sigmaU    = 0
	keySize   = 30
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
	H              structs.Vector[ring.Poly]
	Lambda         ring.Poly
}

var q = ring.Qi60[0]

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
	setupDuration = time.Since(start)

	parties := make([]*Party, k)
	for i := range parties {
		parties[i] = NewParty(i, r, uniformSampler)
	}

	// GEN: Generate secret shares and seeds
	start = time.Now()
	skShares, seeds, b := Gen(r, A, uniformSampler, trustedDealerKey)
	genDuration = time.Since(start)

	for partyID := 0; partyID < k; partyID++ {
		parties[partyID].SkShare = skShares[partyID]
		parties[partyID].Seed = seeds[partyID]
	}

	// SIGNATURE ROUND 1
	mu := "Hello, Threshold Signature!"
	sid := 1
	PRFKey := "PRF Key"
	T := []int{0, 1} // Active parties
	start = time.Now()
	lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(q)))
	D := make(map[int]structs.Matrix[ring.Poly])
	masks := make(map[int]structs.Vector[ring.Poly])
	for _, partyID := range T {
		parties[partyID].Lambda = lagrangeCoeffs[partyID]
		parties[partyID].Seed = seeds[partyID]
		D[partyID], masks[partyID] = parties[partyID].SignRound1(A, sid, []byte(PRFKey), T)
	}
	signRound1Duration = time.Since(start)

	// SIGNATURE ROUND 2
	start = time.Now()
	z := make(map[int]structs.Vector[ring.Poly])
	for _, partyID := range T {
		z[partyID] = parties[partyID].SignRound2(A, b, D, masks, sid, mu, T, []byte(PRFKey), seeds)
	}
	signRound2Duration = time.Since(start)

	// SIGNATURE FINALIZE
	start = time.Now()
	finalParty := parties[0]
	_, sig, Delta := finalParty.SignFinalize(z, masks, A, b)
	finalizeDuration = time.Since(start)

	// Verify the signature
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
	return utils.SamplePolyMatrix(m, n, uniformSampler)
}

// Gen generates the secret shares, seeds, and the public parameter b
func Gen(r *ring.Ring, A structs.Matrix[ring.Poly], uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (map[int]structs.Vector[ring.Poly], map[int][][]byte, structs.Vector[ring.Poly]) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: bound_e}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := utils.SamplePolyVector(n, uniformSampler)
	skShares := ShamirSecretSharing(r, s, t, k)
	e := utils.SamplePolyVector(m, gaussianSampler)
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
	R_i := utils.SamplePolyMatrix(n, d-1, uniformSampler)

	// Initialize E_i
	gaussianParams = ring.DiscreteGaussian{Sigma: sigmaE, Bound: boundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	E_i := utils.SamplePolyMatrix(m, d-1, gaussianSampler)

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

	D := utils.InitializeMatrix(r, m, n)
	utils.MatrixMatrixMul(r, A, concatenatedR, D)
	utils.MatrixAdd(r, concatenatedE, D, D)

	return D, mask
}

// SignRound2 performs the second round of signing
func (party *Party) SignRound2(A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], D map[int]structs.Matrix[ring.Poly], masks map[int]structs.Vector[ring.Poly], sid int, mu string, T []int, PRFKey []byte, seeds map[int][][]byte) structs.Vector[ring.Poly] {
	r := party.Ring
	partyID := party.ID
	concatR := party.R
	s_i := party.SkShare
	lambda := party.Lambda

	u := make(map[int]structs.Vector[ring.Poly], d)
	onePoly := r.NewMonomialXi(0)

	for _, j := range T {
		oneSlice := structs.Vector[ring.Poly]{onePoly}
		u[j] = oneSlice
		if d > 1 {
			u_j := H_u(r, A, b, sid, j, D, mu, masks)
			u[j] = append(oneSlice, u_j...)
		}
	}

	h := utils.InitializeVector(r, m)

	for j, D_j := range D {
		D_j_u_j := utils.InitializeVector(r, m)
		utils.MatrixVectorMul(r, D_j, u[j], D_j_u_j)
		utils.VectorAdd(r, h, D_j_u_j, h)
	}

	for _, poly := range h {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	party.H = h

	c := H_c(r, A, b, h, mu)
	party.C = c

	mPrime := utils.InitializeVector(r, n)

	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j][partyID], PRFKey)
		utils.VectorAdd(r, mPrime, mask_j, mPrime)
	}

	z_i := utils.InitializeVector(r, n)
	utils.MatrixVectorMul(r, concatR, u[partyID], z_i)
	utils.VectorAdd(r, z_i, mPrime, z_i)
	s_c_lambda := utils.InitializeVector(r, n)
	utils.VectorPolyMul(r, s_i, c, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, lambda, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	return z_i
}

// SignFinalize finalizes the signature
func (party *Party) SignFinalize(z map[int]structs.Vector[ring.Poly], masks map[int]structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly]) (ring.Poly, structs.Vector[ring.Poly], structs.Vector[ring.Poly]) {
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

	Az := utils.InitializeVector(r, m)
	utils.MatrixVectorMul(r, A, z_sum, Az)

	bc := utils.InitializeVector(r, m)
	utils.VectorPolyMul(r, b, c, bc)

	Az_bc := utils.InitializeVector(r, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}

	Delta := utils.InitializeVector(r, m)
	utils.VectorSub(r, h, Az_bc, Delta)

	return party.C, z_sum, Delta
}

// Verify verifies the correctness of the signature
func Verify(r *ring.Ring, z structs.Vector[ring.Poly], A structs.Matrix[ring.Poly], mu string, b structs.Vector[ring.Poly], c ring.Poly, Delta structs.Vector[ring.Poly], betaDelta uint64) bool {
	Az := utils.InitializeVector(r, m)
	utils.MatrixVectorMul(r, A, z, Az)

	bc := utils.InitializeVector(r, m)
	utils.VectorPolyMul(r, b, c, bc)

	Az_bc := utils.InitializeVector(r, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}

	Az_bc_Delta := utils.InitializeVector(r, m)
	utils.VectorAdd(r, Az_bc, Delta, Az_bc_Delta)

	computedC := H_c(r, A, b, Az_bc_Delta, mu)
	if !r.Equal(c, computedC) {
		return false
	}

	return CheckInfinityNorm(r, Delta, betaDelta)
}

// CheckInfinityNorm checks if the infinity norm of the vector of polynomial Delta is less than or equal to betaDelta
func CheckInfinityNorm(r *ring.Ring, Delta structs.Vector[ring.Poly], betaDelta uint64) bool {
	maxValue := big.NewInt(0)

	for _, poly := range Delta {
		coeffsBigint := make([]*big.Int, r.N())
		for i := range coeffsBigint {
			coeffsBigint[i] = big.NewInt(0)
		}
		r.PolyToBigintCentered(poly, 1, coeffsBigint)

		for _, coeff := range coeffsBigint {
			absCoeff := new(big.Int).Abs(coeff)
			if absCoeff.Cmp(maxValue) == 1 {
				maxValue.Set(absCoeff)
			}
		}
	}

	betaDeltaBig := new(big.Int).SetUint64(betaDelta)

	log.Print("DELTA NORM:", maxValue)

	return new(big.Int).Mod(maxValue, new(big.Int).SetUint64(q)).Cmp(betaDeltaBig) <= 0
}

// HASHES & PRF

// H_u hashes parameters to a Gaussian distribution
func H_u(r *ring.Ring, A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], sid int, j int, D map[int]structs.Matrix[ring.Poly], mu string, masks map[int]structs.Vector[ring.Poly]) structs.Vector[ring.Poly] {
	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	binary.Write(buf, binary.BigEndian, int64(sid))
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
	gaussianParams := ring.DiscreteGaussian{}
	hashGaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	u_j := utils.SamplePolyVector(d-1, hashGaussianSampler)

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
func H_c(r *ring.Ring, A structs.Matrix[ring.Poly], b structs.Vector[ring.Poly], h structs.Vector[ring.Poly], mu string) ring.Poly {
	hasher := sha3.NewShake128()
	buf := new(bytes.Buffer)

	if _, err := A.WriteTo(buf); err != nil {
		log.Fatalf("Error writing matrix A: %v\n", err)
	}

	if _, err := b.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector b: %v\n", err)
	}

	if _, err := h.WriteTo(buf); err != nil {
		log.Fatalf("Error writing vector h: %v\n", err)
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
				randomCoeff := big.NewInt(20)
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
