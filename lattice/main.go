// lattice/main.go

package main

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"fmt"
	"log"
	"math/big"

	"github.com/hashicorp/vault/shamir"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

const m = 5
const n = 5
const q = 5 // Prime
const p = 4
const t = 3 // Active threshold
const k = 5 // Total number of parties
const d = 5 // Length of joint noise vector
const ell = 5
const beta = 10

func main() {
	LogN := 10
	r, _ := ring.NewRing(1<<LogN, ring.Qi60[:4])
	prng, _ := sampling.NewPRNG()
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := "Trusted dealer key"

	// Setup
	A := Setup(uniformSampler)

	fmt.Println("A: ", A)

	// Gen
	skShares, seeds, b := Gen(r, A, uniformSampler, []byte(trustedDealerKey))

	fmt.Println("Shares: ", skShares)
	fmt.Println("Seeds: ", seeds)

	// Test signing round 1 by simulating each of the parties
	mu := "Hello, Threshold Signature!"
	sid := 1
	numActiveParties := 4
	D := make(map[int]*[][]*ring.Poly)
	concatR := make(map[int]*[][]*ring.Poly)
	m := make(map[int][]*ring.Poly)
	PRFKey := "PRF Key"
	T := []int{1, 2, 3, 4}

	for i := 0; i < numActiveParties; i++ {
		D[i], m[i], concatR[i] = SignRound1(r, uniformSampler, A, i, sid, (*skShares)[i], mu, []byte(PRFKey), seeds[i], T)
	}

	fmt.Println("D: ", D)
	fmt.Println("m: ", m)

	// Testing signing round 2 by simulating each of the parties
	c := make(map[int]*ring.Poly)
	z := make(map[int]*ring.Poly)
	for i := 0; i < numActiveParties; i++ {
		c[i], z[i] = SignRound2(r, i, D, m, A, b, (*skShares)[i], sid, mu, T, []byte(PRFKey), seeds, concatR[i])
	}

	// // Aggregate the signature
	// sig := SignFinalize(z)

	// // Verify the signature
	// valid := Verify(sig, A, mu, b, c, z)
	// fmt.Printf("Signature Verification Result: %v\n", valid)
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
func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (*[][]*ring.Poly, [][][]byte, []*ring.Poly) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

	// Sample the secret key from the ring
	s := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		// Generate a new ring element
		element := uniformSampler.ReadNew()
		s[i] = &element
	}

	fmt.Println("s: ", s)

	// Sample the error from the Gaussian distribution
	e := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		// Generate a new ring element from the Gaussian distribution
		element := gaussianSampler.ReadNew()
		e[i] = &element
	}

	fmt.Println("e: ", e)

	// Initialize result vector
	b := make([]*ring.Poly, m)

	// Compute b = As + e mod q
	for i := 0; i < m; i++ {
		// Initialize b[i] as a new zero polynomial by taking the address of the result from r.NewPoly()
		newPoly := r.NewPoly()
		b[i] = &newPoly

		// Compute A*s for row i
		for j := 0; j < n; j++ {
			// Compute A[i][j] * s[j] mod q and add it to b[i]
			r.MulCoeffsMontgomeryThenAdd(*(*A)[i][j], *s[j], *b[i])
		}

		// Add e[i] to b[i] mod q
		r.Add(*b[i], *e[i], *b[i])
	}

	fmt.Println("b: ", b)

	// Secret-share s
	sharesCoeffs := make([][][]*big.Int, k)
	for l := 0; l < k; l++ {
		sharesCoeffs[l] = make([][]*big.Int, n) // Initialize each ring element vector in the share
		for i := 0; i < n; i++ {
			sharesCoeffs[l][i] = make([]*big.Int, r.N()) // Initialize each coefficient in the ring element vector
		}
	}

	// Indexing over all ring elements of the s vector
	for i := 0; i < n; i++ {
		// These are the coefficients of the ring element s[i]
		coeffsBigint := make([]*big.Int, r.N())
		r.PolyToBigint(*s[i], 1, coeffsBigint)

		// For each coefficient of the ring element polynomial
		for j, coeff := range coeffsBigint {
			coeffSharesBytes, err := shamir.Split(coeff.Bytes(), k, t)
			if err != nil {
				log.Fatalf("Error splitting coefficient: %v\n", err)
			}

			// For each share of this coefficient of the ring element polynomial
			for l := 0; l < k; l++ {
				coeffShare := new(big.Int)
				coeffShare.SetBytes(coeffSharesBytes[l])

				// Set the coefficient share at the k-th secret share, i-th secret vector element, j-th coefficient
				sharesCoeffs[l][i][j] = coeffShare
			}
		}
	}

	fmt.Println("sharesCoeffs: ", sharesCoeffs)

	// Produce secret-shared elements
	skShares := make([][]*ring.Poly, k)
	for l := 0; l < k; l++ {
		skShares[l] = make([]*ring.Poly, n)
		for i := 0; i < n; i++ {
			elem := r.NewPoly()
			r.SetCoefficientsBigint(sharesCoeffs[l][i], elem)
			skShares[l][i] = &elem
		}
	}

	// Generate random seeds for i, j in [d]
	seeds := make([][][]byte, k)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			// Generate random seed sd
			seeds[i][j] = generateRandomSeed()
		}
	}

	return &skShares, seeds, b
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
	m_i := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := r.NewPoly()
		m_i[i] = &newPoly
	}

	for _, j := range T {
		sd_ij := seeds_i[j]
		m_ij := PRF(r, sid, sd_ij, PRFKey)
		for k := 0; k < n; k++ {
			r.Add(*m_i[k], *m_ij[k], *m_i[k])
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
	gaussianParams := ring.DiscreteGaussian{}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

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
	gaussianParams = ring.DiscreteGaussian{}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, true)
	E_i := make([][]*ring.Poly, m)
	for i := 0; i < m; i++ {
		E_i[i] = make([]*ring.Poly, d-1)
		for j := 0; j < d-1; j++ {
			element := gaussianSampler.ReadNew()
			E_i[i][j] = &element
		}
	}

	// Compute D_i here
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

	// Compute D_i = A(concatenatedR) + concatenatedE
	for i := 0; i < m; i++ {
		for j := 0; j < d; j++ {
			newPoly := r.NewPoly()

			for k := 0; k < len(concatenatedR); k++ {
				// Multiply A[i][k] with concatenatedR[k][j] and accumulate in tempPoly
				r.MulCoeffsMontgomeryThenAdd(*(*A)[i][k], *concatenatedR[k][j], newPoly)
			}
			// Add the result of A * concatenatedR with concatenatedE
			r.Add(newPoly, *concatenatedE[i][j], *D_i[i][j])
		}
	}

	return &D_i, m_i, &concatenatedR
}

// Sign function signs a message using the secret key and returns the signature
func SignRound2(r *ring.Ring, partyInt int, DMap map[int]*[][]*ring.Poly, mMap map[int][]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly, s_i []*ring.Poly, sid int, mu string, T []int, PRFKey []byte, seeds [][][]byte, concatR_i *[][]*ring.Poly) (*ring.Poly, *ring.Poly) {
	// Create the noise matrix u
	u := make([][]*ring.Poly, len(T))
	onePoly := r.NewMonomialXi(0)

	for _, j := range T {
		oneSlice := []*ring.Poly{onePoly.CopyNew()}
		u_j := H_u(r, A, b, sid, j, DMap, mu, mMap)
		u[j] = append(oneSlice, u_j...)
	}

	h := make([]*ring.Poly, len(T))

	for j, D_j := range DMap {
		for i := 0; i < d; i++ {
			newPoly := r.NewPoly()
			h[j] = &newPoly

			// Compute D_j * u_j for row i
			for k := 0; k < m; k++ {
				// Compute D_j[i][k] * u_j[k] mod q and add it to h[j][i]
				r.MulCoeffsMontgomeryThenAdd(*(*D_j)[i][k], *u[j][k], *h[j])
			}
		}
	}

	h_sum := r.NewPoly()
	for _, j := range T {
		r.Add(h_sum, *h[j], h_sum)
	}

	// Round h to the nearest multiple of p
	RoundCoeffsToNearestMultiple(r, &h_sum, p)

	// c = H_c
	c := H_c(r, A, b, &h_sum, mu)

	// Compute the column-wise mask from PRFs
	m_i_prime := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := r.NewPoly()
		m_i_prime[i] = &newPoly
	}
	for _, j := range T {
		sd_ji := seeds[j][partyInt]
		m_ij := PRF(r, sid, sd_ji, PRFKey)
		for k := 0; k < n; k++ {
			r.Add(*m_i_prime[k], *m_ij[k], *m_i_prime[k])
		}
	}

	// First term: s_i * c * lambda_T_i
	interpolator, err := ring.NewInterpolator(r.N(), uint64(len(T)))
	if err != nil {
		log.Fatalf("Cannot create interpolator: %v", err)
	}

	var concatSkElems []uint64
	for _, s_iElem := range s_i {
		for i := 0; i < len(s_iElem.Coeffs); i++ {
			concatSkElems = append(concatSkElems, s_iElem.Coeffs[i]...)
		}
	}
	lagrangeCoeffs, err := interpolator.Lagrange([]uint64{uint64(partyInt)}, concatSkElems)
	if err != nil {
		log.Fatalf("Cannot get Lagrange coefficients: %v", err)
	}
	lambda_T_i := lagrangeCoeffs[partyInt] // This is a scalar value

	// First term: Iterate over each polynomial in s_i to compute (s_i * c * lambda_T_i)
	scalerProductSum := r.NewPoly()
	for _, poly := range s_i {
		scalerProduct := r.NewPoly()
		r.MulCoeffsMontgomery(*poly, *c, scalerProduct) // s_i_elem * c
		scaledProduct := r.NewPoly()
		r.MulScalar(scalerProduct, lambda_T_i, scaledProduct)    // (s_i_elem * c) * lambda_T_i
		r.Add(scalerProductSum, scaledProduct, scalerProductSum) // Accumulate results
	}

	// Second term: (r_i_star, R_i) * (1, u_i) mod q
	weightedSum := r.NewPoly()
	r.Add(weightedSum, *(*concatR_i)[partyInt][0], weightedSum) // Start with r_i_star
	for j := 1; j < len((*concatR_i)[partyInt]); j++ {
		temp := r.NewPoly()
		r.MulCoeffsMontgomery(*(*concatR_i)[partyInt][j], *u[partyInt][j-1], temp) // R_i[j-1] * u_i[j-1]
		r.Add(weightedSum, temp, weightedSum)                                      // Accumulate
	}

	// Third term: m_i_prime
	z_i := r.NewPoly()
	r.Add(scalerProductSum, weightedSum, z_i) // Add first two terms (scalerProductSum and weightedSum)
	r.Add(z_i, *m_i_prime[partyInt], z_i)     // Add m_i_prime

	return &z_i, c
}

// An arbitrary party aggregates the signatures
func SignFinalize(z []*ring.Poly) *ring.Poly {
	return nil
}

// An verify
func Verify(sig *ring.Poly, A *[][]*ring.Poly, mu string, b []*ring.Poly, c []*ring.Poly, z []*ring.Poly) bool {
	return false
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

	// Print the hash as a hexadecimal string
	fmt.Printf("SHAKE128 Hash: %x\n", hashOutput)

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	gaussianParams := ring.DiscreteGaussian{}
	hashGaussiamSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

	u_j := make([]*ring.Poly, d-1)
	for i := 0; i < d-1; i++ {
		element := hashGaussiamSampler.ReadNew()
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

	// Print the hash as a hexadecimal string
	fmt.Printf("SHAKE128 Hash: %x\n", hashOutput)

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	prfUniformSampler := ring.NewUniformSampler(prng, r)

	m_i := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		element := prfUniformSampler.ReadNew()
		m_i[i] = &element
	}
	return m_i
}

// Hash to low norm ring elements
func H_c(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, h *ring.Poly, mu string) *ring.Poly {
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
	data, err := h.MarshalBinary()
	if err != nil {
		log.Fatalf("Error marshalling poly: %v\n", err)
	}
	buffer.Write(data)

	// Handle string mu
	binary.Write(&buffer, binary.BigEndian, mu)

	// Write the final concatenated data to the hasher
	_, err = hasher.Write(buffer.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := n
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	// Print the hash as a hexadecimal string
	fmt.Printf("SHAKE128 Hash: %x\n", hashOutput)

	prng, _ := sampling.NewKeyedPRNG(hashOutput)

	lowNormSampler := newLowNormSampler(r)
	poly := lowNormSampler.baseRing.NewPoly()

	for i := range lowNormSampler.coeffs {
		lowNormSampler.coeffs[i] = bignum.RandInt(prng, big.NewInt(beta))
	}

	lowNormSampler.baseRing.AtLevel(poly.Level()).SetCoefficientsBigint(lowNormSampler.coeffs, poly)

	return &poly
}

// Low norm sampler (taken from the VOLE example in LattiGo)
type lowNormSampler struct {
	baseRing *ring.Ring
	coeffs   []*big.Int
}

func newLowNormSampler(baseRing *ring.Ring) (lns *lowNormSampler) {
	lns = new(lowNormSampler)
	lns.baseRing = baseRing
	lns.coeffs = make([]*big.Int, baseRing.N())
	return
}

// Samples a uniform polynomial in Z_{norm}/(X^N + 1)
func (lns *lowNormSampler) newPolyLowNorm(norm *big.Int) (pol ring.Poly) {

	pol = lns.baseRing.NewPoly()

	prng, _ := sampling.NewPRNG()

	for i := range lns.coeffs {
		lns.coeffs[i] = bignum.RandInt(prng, norm)
	}

	lns.baseRing.AtLevel(pol.Level()).SetCoefficientsBigint(lns.coeffs, pol)

	return
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
