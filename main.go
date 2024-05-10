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

const (
	m         = 1
	n         = 1
	p         = 2 // TODO: Experiment
	t         = 1 // Active threshold
	k         = 1 // Total number of parties
	d         = 1 // Length of joint noise vector
	ell       = 1
	beta      = 10
	betaDelta = 3291932728626317921
	kappa     = 10
	logN      = 3
	e_bound   = 0.5 // Keep these as 0 for now!
	sigma_e   = 2   // Standard deviation for the error distribution
	bound_c   = 3
	sigma_c   = 2 // Standard deviation for the hash output distribution
	boundE    = 0.5
	sigmaE    = 2
	boundStar = 0.5
	sigmaStar = 2
	boundU    = 0.5
	sigmaU    = 2
)

var q = uint64(61)

func main() {
	r, _ := ring.NewRing(1<<logN, []uint64{q})
	prng, _ := sampling.NewKeyedPRNG([]byte("0")) // TODO: Change to random key for the PRNG seeding
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := "Trusted dealer key"

	// Setup
	A := Setup(uniformSampler)

	utils.PrintMatrix("A: ", A)

	// Gen
	skShare, seeds, b := Gen(r, A, uniformSampler, []byte(trustedDealerKey))

	// fmt.Println("Secret Key Shares (skShares):")
	// for i, shareRow := range *skShares {
	// 	for j, poly := range shareRow {
	// 		fmt.Printf("skShares[%d][%d]: %v\n", i, j, poly.Coeffs[0])
	// 	}
	// }
	skShares := make([][]*ring.Poly, 1)
	skShares[0] = skShare

	fmt.Println("Seeds:")
	for i, seedRow := range seeds {
		for j, seed := range seedRow {
			fmt.Printf("Seed[%d][%d]: %x\n", i, j, seed) // Printing bytes as hex
		}
	}

	utils.PrintVector("b: ", b)

	// Test signing round 1 by simulating each of the parties
	mu := "Hello, Threshold Signature!"
	sid := 1
	D := make(map[int]*[][]*ring.Poly)
	concatR := make(map[int]*[][]*ring.Poly)
	mask := make(map[int][]*ring.Poly)
	PRFKey := "PRF Key"
	T := []int{0}
	// lagrangeCoeffs := GenLagrangeCoefficients(r, T)

	for i := 0; i < m; i++ {
		// TODO: Change this to be a different seed per user pair
		D[i], mask[i], concatR[i] = SignRound1(r, uniformSampler, A, i, sid, skShares[0], mu, []byte(PRFKey), seeds[0], T)
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

	// Testing signing round 2 by simulating each of the parties
	z := make(map[int][]*ring.Poly)
	c := make(map[int]*ring.Poly)
	h := make(map[int][]*ring.Poly)

	for i := 0; i < m; i++ {
		z[i], c[i], h[i] = SignRound2(r, i, D, nil, A, b, skShares[0], sid, mu, T, []byte(PRFKey), seeds, concatR[i], nil)
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

// Function to generate the secret-shared polynomials
func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, trustedDealerKey []byte) ([]*ring.Poly, [][][]byte, []*ring.Poly) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: e_bound}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	// Sample the secret key from the ring
	s := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		// Generate a new ring element
		element := uniformSampler.ReadNew()
		s[i] = &element
	}

	utils.PrintVector("s: ", s)

	// Sample the error from the Gaussian distribution
	e := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		// Generate a new ring element from the Gaussian distribution
		element := gaussianSampler.ReadNew()
		e[i] = &element
	}

	for i, poly := range e {
		fmt.Printf("e[%d]: %v\n", i, poly.Coeffs[0]) // Accessing the underlying slice
	}

	// Initialize result vector
	b := make([]*ring.Poly, m)

	// Compute b = As + e mod q

	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	// for i := 0; i < m; i++ {
	// 	// Initialize b[i] as a new zero polynomial by taking the address of the result from r.NewPoly()
	// 	newPoly := r.NewPoly()
	// 	b[i] = &newPoly

	// 	// Compute A*s for row i
	// 	for j := 0; j < n; j++ {
	// 		// Compute A[i][j] * s[j] mod q and add it to b[i]

	// 		newPoly := r.NewPoly()
	// 		utils.MulPoly(r, (*A)[i][j], s[j], &newPoly)
	// 		r.Add(*b[i], newPoly, *b[i])
	// 	}

	// 	// Add e[i] to b[i] mod q
	// 	r.Add(*b[i], *e[i], *b[i])
	// }

	// Secret-share s
	// sharesCoeffs := make([][][]*big.Int, k)
	// for l := 0; l < k; l++ {
	// 	sharesCoeffs[l] = make([][]*big.Int, n) // Initialize each ring element vector in the share
	// 	for i := 0; i < n; i++ {
	// 		sharesCoeffs[l][i] = make([]*big.Int, r.N()) // Initialize each coefficient in the ring element vector
	// 	}
	// }

	// Indexing over all ring elements of the s vector
	for i := 0; i < n; i++ {
		// // These are the coefficients of the ring element s[i]
		// coeffsBigint := make([]*big.Int, r.N())
		// r.PolyToBigint(*s[i], 1, coeffsBigint)

		// // For each coefficient of the ring element polynomial
		// for j, coeff := range coeffsBigint {
		// 	coeffSharesBytes, err := shamir.Split(coeff.Bytes(), k, t)
		// 	if err != nil {
		// 		log.Fatalf("Error splitting coefficient: %v\n", err)
		// 	}

		// 	// For each share of this coefficient of the ring element polynomial
		// 	for l := 0; l < k; l++ {
		// 		coeffShare := new(big.Int)
		// 		coeffShare.SetBytes(coeffSharesBytes[l])

		// 		// Set the coefficient share at the k-th secret share, i-th secret vector element, j-th coefficient
		// 		sharesCoeffs[l][i][j] = coeffShare
		// 	}
		// }
	}

	// // Produce secret-shared elements
	// skShares := make([][]*ring.Poly, k)
	// for l := 0; l < k; l++ {
	// 	skShares[l] = make([]*ring.Poly, n)
	// 	for i := 0; i < n; i++ {
	// 		elem := r.NewPoly()
	// 		r.SetCoefficientsBigint(sharesCoeffs[l][i], elem)
	// 		skShares[l][i] = &elem
	// 	}
	// }

	// Generate random seeds for i, j in [d]
	seeds := make([][][]byte, k)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			// Generate random seed sd
			seeds[i][j] = generateRandomSeed()
		}
	}

	return s, seeds, b
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
func SignRound1(r *ring.Ring, uniformSampler *ring.UniformSampler, A *[][]*ring.Poly, partyInt int, sid int, skShare []*ring.Poly, mu string, PRFKey []byte, seeds_i [][]byte, T []int) (*[][]*ring.Poly, []*ring.Poly, *[][]*ring.Poly) {
	// Generate the row-wise mask from PRFs
	mask := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := r.NewPoly()
		mask[i] = &newPoly
	}

	for _, j := range T {
		sd_ij := seeds_i[j]
		mask_j := PRF(r, sid, sd_ij, PRFKey)
		for k := 0; k < n; k++ {
			r.Add(*mask[k], *mask_j[k], *mask[k])
		}
	}

	// Sample the error from the Gaussian distribution
	r_star := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		// Generate a new ring element from the Gaussian distribution
		element := uniformSampler.ReadNew()
		r_star[i] = &element
	}

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

	hashOutputLength := n // TODO: Check
	skHash := make([]byte, hashOutputLength)
	_, err := hasher.Read(skHash)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaStar, Bound: boundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	// Sample the error from the Gaussian distribution
	e_star := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		// Generate a new ring element from the Gaussian distribution
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

	// Initialize D_i as m x d matrix of polynomials
	D_i := make([][]*ring.Poly, m)
	for i := range D_i {
		D_i[i] = make([]*ring.Poly, d)
		for j := range D_i[i] {
			newPoly := r.NewPoly()
			D_i[i][j] = &newPoly
		}
	}

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

	for i := 0; i < m; i++ {
		for j := 0; j < d; j++ {
			newPoly := r.NewPoly()

			for k := 0; k < len(concatenatedR); k++ {
				// Multiply A[i][k] with concatenatedR[k][j] and accumulate
				tempPoly := r.NewPoly()
				utils.MulPoly(r, (*A)[i][k], concatenatedR[k][j], &tempPoly)
				r.Add(newPoly, tempPoly, newPoly)
			}
			// Add the result of A * concatenatedR with concatenatedE
			r.Add(newPoly, *concatenatedE[i][j], *D_i[i][j])
		}
	}

	//TODO: fix

	return &D_i, mask, &concatenatedR
}

// Sign function signs a message using the secret key and returns the signature
func SignRound2(r *ring.Ring, partyInt int, DMap map[int]*[][]*ring.Poly, mMap map[int][]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly, s_i []*ring.Poly, sid int, mu string, T []int, PRFKey []byte, seeds [][][]byte, concatR_i *[][]*ring.Poly, lambda_T_i *ring.Poly) ([]*ring.Poly, *ring.Poly, []*ring.Poly) {
	// // Create the noise matrix u
	// u := make([][]*ring.Poly, m)
	// onePoly := r.NewMonomialXi(0)

	// for j := 0; j < m; j++ {
	// 	oneSlice := []*ring.Poly{onePoly.CopyNew()}
	// 	u_j := H_u(r, A, b, sid, j, DMap, mu, mMap)
	// 	u[j] = append(oneSlice, u_j...)
	// }
	// Create the noise matrix u, now simplified to use onePoly for all elements
	u := make([][]*ring.Poly, m)
	onePoly := r.NewMonomialXi(0) // Assumes onePoly is defined as a monomial of degree 0 (constant polynomial of 1)

	for j := 0; j < m; j++ {
		u[j] = make([]*ring.Poly, d) // Ensure u[j] has d entries
		for k := 0; k < d; k++ {
			u[j][k] = onePoly.CopyNew() // Each entry of u is just a copy of onePoly
		}
	}

	h := make([]*ring.Poly, m)

	for j, D_j := range DMap {
		for i := 0; i < d; i++ {
			newPoly := r.NewPoly()
			h[j] = &newPoly

			// Compute D_j * u_j for row i
			for k := 0; k < m; k++ {
				// Compute D_j[i][k] * u_j[k] mod q and add it to h[j][i]
				tempPoly := r.NewPoly()

				utils.MulPoly(r, (*D_j)[i][k], u[j][k], &tempPoly)
				r.Add(*h[j], tempPoly, *h[j])
			}
		}
	}

	// Round h to the nearest multiple of p
	for _, poly := range h {
		RoundCoeffsToNearestMultiple(r, poly, p)
	}

	// c = H_c
	c := H_c(r, A, b, h, mu)

	// // Compute the column-wise mask from PRFs
	// mask_prime := make([]*ring.Poly, n)
	// for i := 0; i < n; i++ {
	// 	newPoly := r.NewPoly()
	// 	mask_prime[i] = &newPoly
	// }
	// for _, j := range T {
	// 	sd_ji := seeds[j][partyInt]
	// 	mask_j_prime := PRF(r, sid, sd_ji, PRFKey)
	// 	for k := 0; k < n; k++ {
	// 		r.Add(*mask_prime[k], *mask_j_prime[k], *mask_prime[k])
	// 	}
	// }

	// Compute z_i as a vector of ring.Poly
	z_i := make([]*ring.Poly, len(s_i))
	for index, poly := range s_i {
		finalPoly := r.NewPoly()
		utils.MulPoly(r, poly, c, &finalPoly) // s_i_elem * c
		// r.Mul(finalPoly, *lambda_T_i, finalPoly) // (s_i_elem * c) * lambda_T_i TODO, add back when adding lagrange

		// Calculate (concatR_i * u_i) mod q for this poly
		weightedSum := r.NewPoly()
		for j := 0; j < len((*concatR_i)[partyInt]); j++ {
			temp := r.NewPoly()

			utils.MulPoly(r, (*concatR_i)[partyInt][j], u[partyInt][j], &temp) // concatR_i_elem * u_i_elem
			r.Add(weightedSum, temp, weightedSum)                              // Accumulate
		}

		// Combine all terms
		r.Add(finalPoly, weightedSum, finalPoly) // Add weightedSum to (s_i_elem * c) * lambda_T_i
		// r.Add(finalPoly, *mask_prime[partyInt], finalPoly) // Add mask_prime[partyInt] to the finalPoly

		z_i[index] = &finalPoly
	}

	return z_i, c, h
}

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

	hashOutputLength := d - 1
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
		fmt.Printf("H_u element: %v\n", element)
		u_j[i] = &element
	}

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

	hashOutputLength := n
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	prfUniformSampler := ring.NewUniformSampler(prng, r)

	mask := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		element := prfUniformSampler.ReadNew()
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

	hashOutputLength := n
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

	// fmt.Printf("H_c: %v\n", c)

	return &c
}

// RoundPolyCoefficientsToNearestMultiple rounds the coefficients of a polynomial to the nearest multiple of p
// and updates the polynomial using the SetCoefficientsBigint method.
func RoundCoeffsToNearestMultiple(r *ring.Ring, poly *ring.Poly, p uint64) {
	pBig := new(big.Int).SetUint64(p)
	halfP := new(big.Int).Div(new(big.Int).SetUint64(p), big.NewInt(2)) // p/2 for rounding calculation
	roundedCoeffs := make([]*big.Int, poly.N())
	coeffsBigint := make([]*big.Int, poly.N())
	r.PolyToBigint(*poly, 1, coeffsBigint)

	// Calculate rounded coefficients
	for i, coeff := range coeffsBigint {
		// Initialize if nil
		if roundedCoeffs[i] == nil {
			roundedCoeffs[i] = new(big.Int)
		}

		// Perform rounding
		mod := new(big.Int).Mod(coeff, pBig)
		if mod.Cmp(halfP) > 0 {
			coeff.Add(coeff, pBig)
			coeff.Sub(coeff, mod)
		} else {
			coeff.Sub(coeff, mod)
		}
		roundedCoeffs[i].Set(coeff)
	}

	r.SetCoefficientsBigint(roundedCoeffs, *poly)
}

func GenLagrangeCoefficients(r *ring.Ring, T []int) []*ring.Poly {
	lagrangePolynomials := make([]*ring.Poly, len(T)) // TODO: confirm len T

	// For each index in T, we need to create a polynomial l_j(x) that is 1 at x_j and 0 at all other x_m in T
	for j, xj := range T {
		// Start with l_j(x) = 1
		l_j := r.NewPoly()
		r.AddScalar(l_j, 1, l_j)

		// Create the Lagrange polynomial for each xj
		for _, xm := range T {
			if xm != xj {
				// Compute (x - xm)
				tempPoly := r.NewMonomialXi(1)
				r.SubScalar(tempPoly, uint64(xm), tempPoly) // tempPoly = x - xm

				// Multiply l_j by (x - xm)
				resultPoly := r.NewPoly()
				utils.MulPoly(r, &l_j, &tempPoly, &resultPoly)

				l_j = resultPoly
			}
		}

		// Normalize l_j(x) to make l_j(x_j) = 1
		// Compute the denominator product (x_j - xm) for m != j
		denom := uint64(1)
		for _, xm := range T {
			if xm != xj {
				denom *= uint64(xj - xm)
				denom %= q // Perform modulo operation to keep within field limits
			}
		}

		// Compute multiplicative inverse of denom mod q
		denomInv := new(big.Int).ModInverse(big.NewInt(int64(denom)), big.NewInt(int64(q)))

		// Multiply each coefficient of l_j by denomInv
		r.MulScalarBigint(l_j, denomInv, l_j)

		lagrangePolynomials[j] = &l_j
	}

	return lagrangePolynomials
}

func SignFinalize(r *ring.Ring, z map[int][]*ring.Poly, m map[int][]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly, c *ring.Poly, h []*ring.Poly) ([]*ring.Poly, []*ring.Poly) {
	// Initialize z_sum as an array of zero polynomials
	z_sum := make([]*ring.Poly, n)
	for i := range z_sum {
		newPoly := r.NewPoly()
		z_sum[i] = &newPoly
	}

	// Compute z_sum = Σ(z_j - m_j) mod q
	for _, z_j := range z { // Looping through each z_j
		// mask_j := m[j]
		for i := 0; i < n; i++ {
			// r.Sub(*z_j[i], *mask_j[i], *z_j[i])  // z_j[i] - m_j[i] mod q
			r.Add(*z_sum[i], *z_j[i], *z_sum[i]) // Accumulate the result
		}
	}

	// Compute Az
	Az := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := r.NewPoly()
		Az[i] = &newPoly
		for k := 0; k < n; k++ {
			temp := r.NewPoly()
			utils.MulPoly(r, (*A)[i][k], z_sum[k], &temp)
			r.Add(*Az[i], temp, *Az[i])
		}
	}

	fmt.Printf("Az: %v\n", Az[0].Coeffs[0])

	// Compute Az - bc
	bc := make([]*ring.Poly, n)
	Az_bc := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := r.NewPoly()
		bc[i] = &newPoly
		newPoly = r.NewPoly()
		Az_bc[i] = &newPoly
		utils.MulPoly(r, b[i], c, bc[i])

		r.Sub(*Az[i], *bc[i], *Az_bc[i]) // Az[i] - bc[i]
	}

	fmt.Printf("Az - bc: %v\n", Az_bc[0].Coeffs[0])

	// Round Az - bc to the nearest multiple of p
	for _, poly := range Az_bc {
		RoundCoeffsToNearestMultiple(r, poly, p)
	}

	// Compute Δ = [h_p] - [Az - bc]_p
	Delta := make([]*ring.Poly, len(h)) // Assume hLength is the new defined length of vector h
	for i := range Delta {
		newPoly := r.NewPoly()
		Delta[i] = &newPoly
		r.Sub(*h[i], *Az_bc[i], *Delta[i])
	}

	return Delta, z_sum
}

func Verify(r *ring.Ring, z []*ring.Poly, A *[][]*ring.Poly, mu string, b []*ring.Poly, c *ring.Poly, Delta []*ring.Poly, betaDelta uint64) bool {
	// Calculate [Az - bc]_p + Delta
	// computedH := make([]*ring.Poly, len(z))
	// for i := range computedH {
	// 	newPoly := r.NewPoly()
	// 	computedH[i] = &newPoly
	// 	for k := range z {
	// 		temp := r.NewPoly()
	// 		MulPoly(r, (*A)[i][k], z[k], &temp)
	// 		r.Add(*computedH[i], temp, *computedH[i])
	// 	}
	// }

	// bc := make([]*ring.Poly, len(b))
	// for i, bi := range b {
	// 	newPoly := r.NewPoly()
	// 	bc[i] = &newPoly
	// 	MulPoly(r, bi, c, bc[i])                         // Compute b_i * c
	// 	r.Sub(*computedH[i], *bc[i], *computedH[i])      // Compute Az_i - b_i*c
	// 	RoundCoeffsToNearestMultiple(r, computedH[i], p) // Round to nearest multiple of p
	// }

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
