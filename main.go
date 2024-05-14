// lattice/main.go

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
	"golang.org/x/crypto/sha3"
)

// PARAMETERS
const (
	m = 11
	n = 9
	d = 10 // Length of joint noise vector

	p         = 2 ^ 30
	t         = 2 // Active threshold
	k         = 3 // Total number of parties
	ell       = 1
	beta      = 10
	betaDelta = 100000000000
	kappa     = 23
	logN      = 8
	bound_e   = 1000
	sigma_e   = 32
	boundE    = 1000
	sigmaE    = 32
	boundStar = 1000
	sigmaStar = 32
	boundU    = 0
	sigmaU    = 0
	keySize   = 30
)

// Party struct holds all state and methods for a party in the protocol
type Party struct {
	ID             int
	Ring           *ring.Ring
	UniformSampler *ring.UniformSampler
	SkShare        []*ring.Poly
	Seed           [][]byte
	Mask           []*ring.Poly
	R              *[][]*ring.Poly
	C              *ring.Poly
	H              []*ring.Poly
	Lambda         *ring.Poly
}

var q = ring.Qi60[0]

func NewParty(id int, r *ring.Ring, sampler *ring.UniformSampler) *Party {
	return &Party{
		ID:             id,
		Ring:           r,
		UniformSampler: sampler,
	}
}

func main() {
	// Prepare time measurement variables
	var setupDuration, genDuration, signRound1Duration, signRound2Duration, finalizeDuration, verifyDuration time.Duration

	randomKey := make([]byte, keySize)

	r, _ := ring.NewRing(1<<logN, []uint64{q})

	prng, _ := sampling.NewKeyedPRNG(randomKey) // Consider using a truly random key for the PRNG seeding
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := randomKey

	// SETUP
	start := time.Now()
	A := Setup(uniformSampler)
	setupDuration = time.Since(start)
	// utils.printMatrix("A: ", A)

	parties := make([]*Party, k)
	for i := range parties {
		parties[i] = NewParty(i, r, uniformSampler)
	}

	// GEN: Generate secret shares and seeds
	start = time.Now()
	skShares, seeds, b := Gen(r, A, uniformSampler, trustedDealerKey)
	genDuration = time.Since(start)
	// utils.printVector("b: ", b)
	for partyID := 0; partyID < k; partyID++ {
		parties[partyID].SkShare = (*skShares)[partyID]
		parties[partyID].Seed = (*seeds)[partyID]
	}

	// SIGNATURE ROUND 1
	mu := "Hello, Threshold Signature!"
	sid := 1
	PRFKey := "PRF Key"
	T := []int{0, 1} // Active parties
	start = time.Now()
	lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(q)))
	D := make(map[int]*[][]*ring.Poly)
	masks := make(map[int]*[]*ring.Poly)
	for _, partyID := range T {
		parties[partyID].Lambda = lagrangeCoeffs[partyID]
		parties[partyID].Seed = (*seeds)[partyID]
		D[partyID], masks[partyID] = parties[partyID].SignRound1(A, sid, mu, []byte(PRFKey), T)
	}
	signRound1Duration = time.Since(start)
	// utils.printPolynomial("Lagrange:", lagrangeCoeffs[0])
	// utils.printVector("sk:", (*skShares)[0])

	// SIGNATURE ROUND 2
	start = time.Now()
	z := make(map[int]*[]*ring.Poly)
	for _, partyID := range T {
		z[partyID] = parties[partyID].SignRound2(A, b, &D, &masks, sid, mu, T, []byte(PRFKey), seeds)
	}
	signRound2Duration = time.Since(start)

	// SIGNATURE FINALIZE
	start = time.Now()
	finalParty := parties[0] // Example: Let the first party finalize the signature
	_, sig, Delta := finalParty.SignFinalize(z, masks, A, b)
	finalizeDuration = time.Since(start)
	// utils.printVector("Signature: ", *sig)
	// utils.printPolynomial("c: ", c)
	// utils.printVector("Delta: ", *Delta)

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

// Generate the public parameters
func Setup(uniformSampler *ring.UniformSampler) *[][]*ring.Poly {
	A := utils.SamplePolyMatrix(m, n, uniformSampler)
	return &A
}

func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (*map[int][]*ring.Poly, *map[int][][]byte, []*ring.Poly) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: bound_e}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := utils.SamplePolyVector(n, uniformSampler)
	// utils.printVector("s:", s)

	// Share the secret key vector
	skShares := ShamirSecretSharing(r, s, t, k) // Shares the secret 's' across 'k' parties with threshold 't'

	e := utils.SamplePolyVector(m, gaussianSampler)
	// utils.printVector("e: ", e)

	// Compute b = As + e mod q
	b := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	// Generate random seeds for all possible i, j
	seeds := make(map[int][][]byte)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			seeds[i][j] = []byte{0, 1, 2, 3, 4, 5, 6, 7, 8, 9} // TODO: change to random seed again
		}
	}

	return &skShares, &seeds, b
}

func generateRandomSeed() []byte {
	// Generate a random binary string of length ell
	sd := make([]byte, ell)
	_, err := rand.Read(sd)
	if err != nil {
		log.Fatalf("Error generating random seed %v\n", err)
	}
	return sd
}

func (party *Party) SignRound1(A *[][]*ring.Poly, sid int, mu string, PRFKey []byte, T []int) (*[][]*ring.Poly, *[]*ring.Poly) {
	r := party.Ring
	uniformSampler := party.UniformSampler
	seeds := party.Seed

	// Generate the row-wise mask from PRFs
	mask := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j], PRFKey)
		utils.VectorAdd(r, mask, mask_j, mask)
	}

	// Sample r*
	r_star := utils.SamplePolyVector(n, uniformSampler)

	// Sample e*
	skHash := party.PRNGKey()
	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaStar, Bound: boundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)
	e_star := utils.SamplePolyVector(m, gaussianSampler)

	// Sample the R_i matrix
	R_i := utils.SamplePolyMatrix(n, d-1, uniformSampler)

	// Sample the E_i matrix
	gaussianParams = ring.DiscreteGaussian{Sigma: sigmaE, Bound: boundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	E_i := utils.SamplePolyMatrix(m, d-1, gaussianSampler)

	// Concatenate r_i_star with R_i, and e_i_star with E_i
	concatenatedR := make([][]*ring.Poly, n)
	concatenatedE := make([][]*ring.Poly, m)

	for i := 0; i < n; i++ {
		concatenatedR[i] = append([]*ring.Poly{r_star[i]}, R_i[i]...)
	}
	party.R = &concatenatedR

	for i := 0; i < m; i++ {
		concatenatedE[i] = append([]*ring.Poly{e_star[i]}, E_i[i]...)
	}
	// utils.printMatrix("concatE: ", &concatenatedE)

	// Compute D_i = A(concatenatedR) + concatenatedE
	D := make([][]*ring.Poly, m)
	utils.MatrixMatrixMul(r, A, &concatenatedR, &D)
	utils.MatrixAdd(r, &concatenatedE, &D, &D)

	// Broacasted values
	return &D, &mask
}

func (party *Party) SignRound2(A *[][]*ring.Poly, b []*ring.Poly, D *map[int]*[][]*ring.Poly, masks *map[int]*[]*ring.Poly, sid int, mu string, T []int, PRFKey []byte, seeds *map[int][][]byte) *[]*ring.Poly {
	r := party.Ring
	partyID := party.ID
	concatR := party.R
	s_i := party.SkShare
	lambda := party.Lambda

	// Create the noise vectors u
	u := make(map[int][]*ring.Poly, d)
	onePoly := r.NewMonomialXi(0)

	for _, j := range T {
		oneSlice := []*ring.Poly{onePoly.CopyNew()}
		u[j] = oneSlice
		if d > 1 {
			u_j := H_u(r, A, b, sid, j, D, mu, masks)
			u[j] = append(oneSlice, u_j...)
		}
	}

	h := make([]*ring.Poly, m)

	for j, D_j := range *D {
		D_j_u_j := make([]*ring.Poly, m)
		// utils.printMatrix("D_j", D_j)
		// utils.printVector("u[j]", u[j])

		utils.MatrixVectorMul(r, D_j, u[j], D_j_u_j)

		// utils.printVector("D_j_u_j", D_j_u_j)

		utils.VectorAdd(r, h, D_j_u_j, h)
		// utils.printVector("h rn", h)

	}

	// Round h to the nearest multiple of p
	for _, poly := range h {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	party.H = h

	// c = H_c
	c := H_c(r, A, b, h, mu)
	party.C = c

	// Generate the column-wise mask from PRFs
	mPrime := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, (*seeds)[j][partyID], PRFKey)
		utils.VectorAdd(r, mPrime, mask_j, mPrime)
	}

	// Compute z_i as a vector of ring.Poly
	z_i := make([]*ring.Poly, n)
	utils.MatrixVectorMul(r, concatR, u[partyID], z_i)
	// utils.printVector("z after R*u", z_i)
	utils.VectorAdd(r, z_i, mPrime, z_i)
	// utils.printVector("z after R*u", z_i)
	s_c_lambda := make([]*ring.Poly, n)
	utils.VectorPolyMul(r, s_i, c, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, lambda, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	// Broadcast z_i
	return &z_i
}

func (party *Party) SignFinalize(z map[int]*[]*ring.Poly, masks map[int]*[]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly) (*ring.Poly, *[]*ring.Poly, *[]*ring.Poly) {
	r := party.Ring
	c := party.C
	h := party.H

	// Initialize z_sum as an array of zero polynomials
	z_sum := make([]*ring.Poly, n)

	// Compute z_sum = Σ(z_j - m_j) mod q
	for j, z_j := range z { // Looping through each z_j
		mask_j := masks[j]
		for i := 0; i < n; i++ {
			if z_sum[i] == nil {
				newPoly := r.NewPoly()
				z_sum[i] = &newPoly
			}
			r.Sub(*(*z_j)[i], *(*mask_j)[i], *(*z_j)[i]) // z_j[i] - m_j[i] mod q
			r.Add(*z_sum[i], *(*z_j)[i], *z_sum[i])
		}
	}

	// utils.printVector("z_sum", z_sum)

	// Compute Az using MatrixVectorMul
	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, z_sum, Az)

	// Compute bc using VectorPolyMul
	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	// Subtract bc from Az to get Az_bc
	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	// utils.printVector("Az - bc: ", Az_bc)

	// Round Az_bc to the nearest multiple of p
	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	// utils.printVector("Rounded Az_bc: ", Az_bc)

	// Compute Δ = [h_p] - [Az - bc]_p
	Delta := make([]*ring.Poly, m)
	utils.VectorSub(r, h, Az_bc, Delta)

	return party.C, &z_sum, &Delta
}

func Verify(r *ring.Ring, z *[]*ring.Poly, A *[][]*ring.Poly, mu string, b []*ring.Poly, c *ring.Poly, Delta *[]*ring.Poly, betaDelta uint64) bool {
	// Compute Az using MatrixVectorMul
	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, *z, Az)

	// Compute bc using VectorPolyMul
	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	// Subtract bc from Az to get Az_bc
	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	// utils.printVector("Az - bc verification: ", Az_bc)

	// Round Az_bc to the nearest multiple of p
	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	// utils.printVector("Rounded Az_bc: ", Az_bc)

	Az_bc_Delta := make([]*ring.Poly, m)
	utils.VectorAdd(r, Az_bc, *Delta, Az_bc_Delta)

	// utils.printVector("Rounded Az_bc_delta: ", Az_bc_Delta)

	// Verify that c equals H_c([Az - bc]_p + Delta, mu)
	computedC := H_c(r, A, b, Az_bc_Delta, mu)
	if !r.Equal(*c, *computedC) {
		return false
	}

	// utils.printPolynomial("Computed c: %v", computedC)

	// Verify that ||Delta||_inf <= betaDelta
	return CheckInfinityNorm(r, Delta, betaDelta)
}

// CheckInfinityNorm checks if the infinity norm of the vector of polynomial Delta is less than or equal to betaDelta
func CheckInfinityNorm(r *ring.Ring, Delta *[]*ring.Poly, betaDelta uint64) bool {
	maxValue := big.NewInt(0) // Temporary variable to store the maximum value found

	for _, poly := range *Delta {
		coeffsBigint := make([]*big.Int, r.N())
		for i := range coeffsBigint {
			coeffsBigint[i] = big.NewInt(0)
		}
		r.PolyToBigintCentered(*poly, 1, coeffsBigint)

		for _, coeff := range coeffsBigint {
			absCoeff := new(big.Int).Abs(coeff) // Get the absolute value of the coefficient
			if absCoeff.Cmp(maxValue) == 1 {    // Compare absCoeff with maxValue
				maxValue.Set(absCoeff) // Update maxValue if absCoeff is greater
			}
		}
	}

	betaDeltaBig := new(big.Int).SetUint64(betaDelta)

	log.Print("DELTA NORM:", maxValue)

	// Check if the maximum absolute value is less than or equal to betaDelta
	return new(big.Int).Mod(maxValue, new(big.Int).SetUint64(q)).Cmp(betaDeltaBig) <= 0
}

// HASHES & PRF

// Hash parameters to a Gaussian distribution
func H_u(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, sid int, j int, D *map[int]*[][]*ring.Poly, mu string, masks *map[int]*[]*ring.Poly) []*ring.Poly {
	hasher := sha3.NewShake128()

	// Buffer to store all concatenated bytes
	var buffer bytes.Buffer

	// Handle matrix A
	for _, row := range *A {
		for _, poly := range row {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	// Handle slice b
	for _, poly := range b {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	// Handle integer sid
	binary.Write(&buffer, binary.BigEndian, int64(sid))

	// Handle integer j
	binary.Write(&buffer, binary.BigEndian, int64(j))

	// Handle map of matrices D
	for _, D_h := range *D {
		for _, row := range *D_h {
			for _, poly := range row {
				data, err := poly.MarshalBinary()
				if err != nil {
					log.Fatalf("Error marshalling poly: %v\n", err)
				}
				buffer.Write(data)
			}
		}
	}

	// Handle string mu
	buffer.WriteString(mu)

	// Handle map of slices m
	for _, mask := range *masks {
		for _, poly := range *mask {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	// Write the final concatenated data to the hasher
	_, err := hasher.Write(buffer.Bytes())
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

	// utils.printVector("u_j", u_j)

	return u_j
}

// PRF to ring elements
func PRF(r *ring.Ring, sid int, sd_ij []byte, PRFKey []byte) []*ring.Poly {
	hasher := sha3.NewShake128()

	// Buffer to store all concatenated bytes
	var buffer bytes.Buffer

	// Handle PRF key
	binary.Write(&buffer, binary.BigEndian, PRFKey)

	// Handle seed
	binary.Write(&buffer, binary.BigEndian, sd_ij)

	// Handle integer sid
	binary.Write(&buffer, binary.BigEndian, int64(sid))

	// Write the final concatenated data to the hasher
	_, err := hasher.Write(buffer.Bytes())
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

// Hash to low norm ring elements
func H_c(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, h []*ring.Poly, mu string) *ring.Poly {
	hasher := sha3.NewShake128()

	// Buffer to store all concatenated bytes
	var buffer bytes.Buffer

	// Handle matrix A
	for _, row := range *A {
		for _, poly := range row {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	// Handle slice b
	for _, poly := range b {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	// Handle slice h
	for _, poly := range h {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	// Handle string mu
	binary.Write(&buffer, binary.BigEndian, []byte(mu))

	// Write the final concatenated data to the hasher
	_, err := hasher.Write(buffer.Bytes())
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

	return &c
}

// ShamirSecretSharing shares each coefficient of a vector of ring.Poly across k parties using (t, k)-threshold Shamir secret sharing.
func ShamirSecretSharing(r *ring.Ring, s []*ring.Poly, t, k int) map[int][]*ring.Poly {
	degree := r.N() // Number of coefficients in each ring.Poly

	// Initialize shares for each party
	shares := make(map[int][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		shares[i] = make([]*ring.Poly, len(s))
		for j := range shares[i] {
			newPoly := r.NewPoly()
			shares[i][j] = &newPoly
		}
	}

	// Perform secret sharing for each coefficient in each polynomial
	for polyIndex, poly := range s {
		coeffs := make([]*big.Int, degree)
		r.PolyToBigint(*poly, 1, coeffs) // Extract bigint coefficients from the polynomial

		for coeffIndex, secret := range coeffs {
			// Create a polynomial for sharing this coefficient
			polyCoeffs := make([]*big.Int, t)
			polyCoeffs[0] = secret // the secret is the coefficient itself
			for i := 1; i < t; i++ {
				// randomCoeff, _ := rand.Int(rand.Reader, big.NewInt(int64(q))) // TODO make random
				randomCoeff := big.NewInt(20)
				polyCoeffs[i] = randomCoeff
			}

			// Evaluate this polynomial at points 1 to k
			for i := 1; i <= k; i++ {
				x := big.NewInt(int64(i))
				shareValue := big.NewInt(0)
				xPow := big.NewInt(1) // x^0

				for _, coeff := range polyCoeffs {
					term := new(big.Int).Mul(coeff, xPow)
					shareValue.Add(shareValue, term)
					shareValue.Mod(shareValue, big.NewInt(int64(q)))
					xPow.Mul(xPow, x)
				}

				// Update the coefficient in the share polynomial
				if shares[i-1][polyIndex].Coeffs[0] == nil {
					shares[i-1][polyIndex].Coeffs[0] = make([]uint64, degree)
				}

				// Set the calculated share value to the specific coefficient
				shares[i-1][polyIndex].Coeffs[0][coeffIndex] = shareValue.Uint64()
			}
		}
	}

	return shares
}

// ComputeLagrangeCoefficients computes the Lagrange coefficients for interpolation based on the indices of available shares.
// It returns a map where keys are party IDs from the list T, and values are pointers to ring.Poly objects containing the coefficients.
func ComputeLagrangeCoefficients(r *ring.Ring, T []int, modulus *big.Int) map[int]*ring.Poly {
	lagrangeCoefficients := make(map[int]*ring.Poly)

	for i := 0; i < len(T); i++ {
		xi := big.NewInt(int64(T[i] + 1)) // Convert party ID to x value (x values are 1-based in Shamir's scheme)
		numerator := big.NewInt(1)
		denominator := big.NewInt(1)

		for j := 0; j < len(T); j++ {
			if i != j {
				xj := big.NewInt(int64(T[j] + 1))
				// numerator *= (0 - xj)
				numerator.Mul(numerator, new(big.Int).Neg(xj))
				numerator.Mod(numerator, modulus)
				// denominator *= (xi - xj)
				temp := new(big.Int).Sub(xi, xj)
				denominator.Mul(denominator, temp)
				denominator.Mod(denominator, modulus)
			}
		}

		// Calculate the modular inverse of the denominator
		denomInv := new(big.Int).ModInverse(denominator, modulus)
		coeff := new(big.Int).Mul(numerator, denomInv)
		coeff.Mod(coeff, modulus)

		// Create a polynomial with this coefficient
		lagrangePoly := r.NewPoly()
		r.SetCoefficientsBigint([]*big.Int{coeff}, lagrangePoly)
		lagrangeCoefficients[T[i]] = &lagrangePoly
	}

	return lagrangeCoefficients
}

func (party *Party) PRNGKey() []byte {
	skShare := party.SkShare

	hasher := sha3.NewShake128()
	// Buffer to store all concatenated bytes
	var buffer bytes.Buffer
	// Handle slice sk
	for _, poly := range skShare {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	hasher.Write(buffer.Bytes())

	hashOutputLength := keySize
	skHash := make([]byte, hashOutputLength)
	_, err := hasher.Read(skHash)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}
	return skHash
}
