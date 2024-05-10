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

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

// PARAMETERS
const (
	m         = 10
	n         = 10
	p         = 10  // TODO: Experiment
	t         = 4   // Active threshold
	k         = 6   // Total number of parties
	d         = 100 // Length of joint noise vector
	ell       = 1
	beta      = 10
	betaDelta = 3291932728626317921
	kappa     = 10
	logN      = 3
	e_bound   = 0.5 // Keep these as 0 for now!
	sigma_e   = 2   // Standard deviation for the error distribution
	bound_c   = 3
	sigma_c   = 6 // Standard deviation for the hash output distribution
	boundE    = 0.5
	sigmaE    = 4
	boundStar = 0.5
	sigmaStar = 2
	boundU    = 0.5
	sigmaU    = 2
	keySize   = 30
)

// var q = uint64(61)
var q = ring.Qi60[0]

func main() {
	randomKey := make([]byte, keySize)
	if _, err := rand.Read(randomKey); err != nil {
		log.Fatalf("failed to generate random key: %v", err)
	}
	fmt.Printf("Random PRNG Key: %v", randomKey)

	r, _ := ring.NewRing(1<<logN, []uint64{q})
	prng, _ := sampling.NewKeyedPRNG(randomKey) // TODO: Change to random key for the PRNG seeding
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := randomKey

	// Setup
	A := Setup(uniformSampler)
	utils.PrintMatrix("A: ", A)

	// Gen
	skShares, seeds, b := Gen(r, A, uniformSampler, []byte(trustedDealerKey))

	utils.PrintVector("b: ", b)

	// ROUND 1
	mu := "Hello, Threshold Signature!"
	sid := 1
	D := make(map[int]*[][]*ring.Poly)
	concatR := make(map[int]*[][]*ring.Poly)
	mask := make(map[int][]*ring.Poly)
	PRFKey := "PRF Key"
	T := []int{0, 1, 2, 3}
	lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(q)))

	for partyInt := 0; partyInt < len(T); partyInt++ {
		D[partyInt], mask[partyInt], concatR[partyInt] = SignRound1(r, uniformSampler, A, partyInt, sid, (*skShares)[partyInt], mu, []byte(PRFKey), (*seeds)[partyInt], T)
	}

	for key, matrix := range D {
		utils.PrintMatrix("D_"+fmt.Sprint(key), matrix)
	}
	for key, vec := range mask {
		utils.PrintVector("Mask_m_"+fmt.Sprint(key), vec)
	}
	for key, matrix := range concatR {
		utils.PrintMatrix("concatR_"+fmt.Sprint(key), matrix)
	}

	// ROUND 2
	z := make(map[int][]*ring.Poly)
	c := make(map[int]*ring.Poly)
	h := make(map[int][]*ring.Poly)

	for partyInt := 0; partyInt < len(T); partyInt++ {
		z[partyInt], c[partyInt], h[partyInt] = SignRound2(r, partyInt, D, mask, A, b, (*skShares)[partyInt], sid, mu, T, []byte(PRFKey), seeds, concatR[partyInt], lagrangeCoeffs[partyInt])
	}

	// After computation in SignRound2
	for i, polys := range z {
		utils.PrintVector("z_"+fmt.Sprint(i), polys)
	}

	for i, poly := range c {
		utils.PrintPolynomial("c_"+fmt.Sprint(i), poly)
	}

	for i, poly := range h {
		utils.PrintVector("h_"+fmt.Sprint(i), poly)
	}

	// Aggregate the signature
	// TODO: Update to take in only the local user's h value, this part is not broadcasted
	Delta, sig := SignFinalize(r, z, mask, A, b, c[0], h[0])
	utils.PrintVector("Delta: ", Delta)
	utils.PrintVector("Signature: ", sig)

	// Verify the signature
	valid := Verify(r, sig, A, mu, b, c[0], Delta, betaDelta)

	fmt.Printf("Signature Verification Result: %v\n", valid)
}

// Generate the public parameters
func Setup(uniformSampler *ring.UniformSampler) *[][]*ring.Poly {
	// Create a new matrix
	A := make([][]*ring.Poly, m)

	// Generate elements for each row and column
	for i := 0; i < m; i++ {
		A[i] = make([]*ring.Poly, n)
		for j := 0; j < n; j++ {
			// Generate a new ring element
			element := uniformSampler.ReadNew()
			A[i][j] = &element
		}
	}
	return &A
}

func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (*map[int][]*ring.Poly, *map[int][][]byte, []*ring.Poly) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: e_bound}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		element := uniformSampler.ReadNew()
		s[i] = &element
	}

	// Share the secret key vector
	skShares := ShamirSecretSharing(r, s, t, k) // Shares the secret 's' across 'k' parties with threshold 't'

	e := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		element := gaussianSampler.ReadNew()
		e[i] = &element
	}
	utils.PrintVector("e: ", e)

	// Compute b = As + e mod q
	b := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	// Generate random seeds for all possible i, j
	seeds := make(map[int][][]byte)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			seeds[i][j] = generateRandomSeed()
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

// Sign function signs a message using the secret key and returns the signature
func SignRound1(r *ring.Ring, uniformSampler *ring.UniformSampler, A *[][]*ring.Poly, partyInt int, sid int, skShare []*ring.Poly, mu string, PRFKey []byte, seeds [][]byte, T []int) (*[][]*ring.Poly, []*ring.Poly, *[][]*ring.Poly) {
	// Generate the row-wise mask from PRFs
	mask := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j], PRFKey)
		utils.VectorAdd(r, mask, mask_j, mask)
	}

	// Sample r*
	r_star := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		element := uniformSampler.ReadNew()
		r_star[i] = &element
	}

	// Sample e*
	skHash := userPRNGKey(skShare)
	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaStar, Bound: boundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)
	e_star := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		element := gaussianSampler.ReadNew()
		e_star[i] = &element
	}

	// Sample the R_i matrix
	R_i := make([][]*ring.Poly, n)
	for i := 0; i < n; i++ {
		R_i[i] = make([]*ring.Poly, d-1)
		for j := 0; j < d-1; j++ {
			element := uniformSampler.ReadNew()
			R_i[i][j] = &element
		}
	}

	// Sample the E_i matrix
	gaussianParams = ring.DiscreteGaussian{Sigma: sigmaE, Bound: boundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	E_i := make([][]*ring.Poly, m)
	for i := 0; i < m; i++ {
		E_i[i] = make([]*ring.Poly, d-1)
		for j := 0; j < d-1; j++ {
			element := gaussianSampler.ReadNew()
			E_i[i][j] = &element
		}
	}
	// TODO: Make a util function for sampling a new matrix or vector

	// Concatenate r_i_star with R_i, and e_i_star with E_i
	concatenatedR := make([][]*ring.Poly, n)
	concatenatedE := make([][]*ring.Poly, m)

	for i := 0; i < n; i++ {
		concatenatedR[i] = append([]*ring.Poly{r_star[i]}, R_i[i]...)
	}
	for i := 0; i < m; i++ {
		concatenatedE[i] = append([]*ring.Poly{e_star[i]}, E_i[i]...)
	}
	utils.PrintMatrix("concatE: ", &concatenatedE)

	// Compute D_i = A(concatenatedR) + concatenatedE
	D := make([][]*ring.Poly, m)
	utils.MatrixMatrixMul(r, A, &concatenatedR, &D)
	utils.MatrixAdd(r, &concatenatedE, &D, &D)

	return &D, mask, &concatenatedR
}

// Sign function signs a message using the secret key and returns the signature
func SignRound2(r *ring.Ring, partyInt int, DMap map[int]*[][]*ring.Poly, mMap map[int][]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly, s_i []*ring.Poly, sid int, mu string, T []int, PRFKey []byte, seeds *map[int][][]byte, concatR *[][]*ring.Poly, lambda *ring.Poly) ([]*ring.Poly, *ring.Poly, []*ring.Poly) {
	// Create the noise vectors u
	u := make(map[int][]*ring.Poly, d)
	onePoly := r.NewMonomialXi(0)

	for _, j := range T {
		oneSlice := []*ring.Poly{onePoly.CopyNew()}
		if d > 1 {
			u_j := H_u(r, A, b, sid, j, DMap, mu, mMap)
			u[j] = append(oneSlice, u_j...)
		}
	}

	h := make([]*ring.Poly, m)

	for j, D_j := range DMap {
		D_j_u_j := make([]*ring.Poly, m)
		utils.MatrixVectorMul(r, D_j, u[j], D_j_u_j)
		utils.VectorAdd(r, h, D_j_u_j, h)
	}

	// Round h to the nearest multiple of p
	for _, poly := range h {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}

	// c = H_c
	c := H_c(r, A, b, h, mu)

	// Generate the column-wise mask from PRFs
	mPrime := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, (*seeds)[j][partyInt], PRFKey)
		utils.VectorAdd(r, mPrime, mask_j, mPrime)
	}

	// Compute z_i as a vector of ring.Poly
	z_i := make([]*ring.Poly, n)
	utils.MatrixVectorMul(r, concatR, u[partyInt], z_i)
	utils.VectorAdd(r, z_i, mPrime, z_i)
	s_c_lambda := make([]*ring.Poly, n)
	utils.VectorPolyMul(r, s_i, c, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, lambda, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	return z_i, c, h
}

func SignFinalize(r *ring.Ring, z map[int][]*ring.Poly, mMap map[int][]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly, c *ring.Poly, h []*ring.Poly) ([]*ring.Poly, []*ring.Poly) {
	// Initialize z_sum as an array of zero polynomials
	z_sum := make([]*ring.Poly, n)

	// Compute z_sum = Σ(z_j - m_j) mod q
	for j, z_j := range z { // Looping through each z_j
		mask_j := mMap[j]
		for i := 0; i < n; i++ {
			newPoly := r.NewPoly()
			z_sum[i] = &newPoly
			r.Sub(*z_j[i], *mask_j[i], *z_j[i]) // z_j[i] - m_j[i] mod q
			r.Add(*z_sum[i], *z_j[i], *z_sum[i])
		}
	}

	// Compute Az using MatrixVectorMul
	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, z_sum, Az)

	// Compute bc using VectorPolyMul
	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	// Subtract bc from Az to get Az_bc
	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	utils.PrintVector("Az - bc: ", Az_bc)

	// Round Az_bc to the nearest multiple of p
	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	utils.PrintVector("Rounded Az_bc: ", Az_bc)

	// Compute Δ = [h_p] - [Az - bc]_p
	Delta := make([]*ring.Poly, m)
	utils.VectorSub(r, h, Az_bc, Delta)

	return Delta, z_sum
}

func Verify(r *ring.Ring, z []*ring.Poly, A *[][]*ring.Poly, mu string, b []*ring.Poly, c *ring.Poly, Delta []*ring.Poly, betaDelta uint64) bool {
	// Compute Az using MatrixVectorMul
	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, z, Az)

	// Compute bc using VectorPolyMul
	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	// Subtract bc from Az to get Az_bc
	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	utils.PrintVector("Az - bc: ", Az_bc)

	// Round Az_bc to the nearest multiple of p
	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	utils.PrintVector("Rounded Az_bc: ", Az_bc)

	Az_bc_Delta := make([]*ring.Poly, m)
	utils.VectorAdd(r, Az_bc, Delta, Az_bc_Delta)

	utils.PrintVector("Rounded Az_bc_delta: ", Az_bc_Delta)

	// Verify that c equals H_c([Az - bc]_p + Delta, mu)
	computedC := H_c(r, A, b, Az_bc_Delta, mu)
	if !r.Equal(*c, *computedC) {
		return false
	}

	utils.PrintPolynomial("computed c: %v", computedC)

	// Verify that ||Delta||_inf <= betaDelta
	if !checkInfinityNorm(r, Delta, betaDelta) {
		return false
	}

	return true
}

// checkInfinityNorm checks if the infinity norm of the vector of polynomial Delta is less than or equal to betaDelta
func checkInfinityNorm(r *ring.Ring, Delta []*ring.Poly, betaDelta uint64) bool {
	maxValue := big.NewInt(0) // Temporary variable to store the maximum value found

	for _, poly := range Delta {
		coeffsBigint := make([]*big.Int, poly.N())
		r.PolyToBigint(*poly, 1, coeffsBigint)

		for _, coeff := range coeffsBigint {
			absCoeff := new(big.Int).Abs(coeff) // Get the absolute value of the coefficient
			if absCoeff.Cmp(maxValue) == 1 {    // Compare absCoeff with maxValue
				maxValue.Set(absCoeff) // Update maxValue if absCoeff is greater
			}
		}
	}

	// Convert betaDelta to big.Int for comparison
	betaDeltaBig := new(big.Int).SetUint64(betaDelta)

	log.Printf("Norm: %v", maxValue)

	// Check if the maximum absolute value is less than or equal to betaDelta
	return maxValue.Cmp(betaDeltaBig) <= 0
}

func userPRNGKey(skShare []*ring.Poly) []byte {
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

// HASHES & PRF

// Hash parameters to a Gaussian distribution
func H_u(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, sid int, j int, DMap map[int]*[][]*ring.Poly, mu string, mMap map[int][]*ring.Poly) []*ring.Poly {
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
	binary.Write(&buffer, binary.BigEndian, sid)

	// Handle integer j
	binary.Write(&buffer, binary.BigEndian, j)

	// Handle map of matrices D
	for _, D := range DMap {
		for _, row := range *D {
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
	for _, m := range mMap {
		for _, poly := range m {
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
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaU, Bound: boundU}
	hashGaussiamSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	u_j := make([]*ring.Poly, d-1)
	for i := 0; i < d-1; i++ {
		element := hashGaussiamSampler.ReadNew()
		u_j[i] = &element
	}

	utils.PrintVector("u_j", u_j)

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
	binary.Write(&buffer, binary.BigEndian, sid)

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

	mask := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		element := PRFUniformSampler.ReadNew()
		mask[i] = &element
	}
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
	binary.Write(&buffer, binary.BigEndian, mu)

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
	modulus := q    // Assuming using the first modulus for simplicity
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
				randomCoeff, _ := rand.Int(rand.Reader, big.NewInt(int64(modulus)))
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
					shareValue.Mod(shareValue, big.NewInt(int64(modulus)))
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
// It returns a map where keys are indices, and values are pointers to ring.Poly objects containing the coefficients.
func ComputeLagrangeCoefficients(r *ring.Ring, indices []int, modulus *big.Int) map[int]*ring.Poly {
	k := len(indices)
	lagrangeCoefficients := make(map[int]*ring.Poly)

	for i := 0; i < k; i++ {
		xi := big.NewInt(int64(indices[i] + 1)) // Convert index to x value (x values are 1-based in Shamir's scheme)
		numerator := big.NewInt(1)
		denominator := big.NewInt(1)

		for j := 0; j < k; j++ {
			if i != j {
				xj := big.NewInt(int64(indices[j] + 1))
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
		lagrangeCoefficients[i] = &lagrangePoly
	}

	return lagrangeCoefficients
}
