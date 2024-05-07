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
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

const m = 5
const n = 5
const q = 5 // Prime
const t = 3 // Active threshold
const k = 5 // Total number of parties
const d = 5 // Length of joint noise vector
const ell = 5

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
	D := make([]*ring.Poly, numActiveParties)
	m := make([]*ring.Poly, numActiveParties)
	PRFKey := "PRF Key"
	T := []int{1, 2, 3, 4}

	for i := 0; i < numActiveParties; i++ {
		D[i], m[i] = SignRound1(r, uniformSampler, i, sid, skShares[i], mu, []byte(PRFKey), seeds[i], T)
	}

	// Testing signing round 2 by simulating each of the parties
	c := make([]*ring.Poly, numActiveParties)
	z := make([]*ring.Poly, numActiveParties)
	for i := 0; i < numActiveParties; i++ {
		c[i], z[i] = SignRound2(i, D, m)
	}

	// Aggregate the signature
	sig := SignFinalize(z)

	// Verify the signature
	valid := Verify(sig, A, mu, b, c, z)
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
	fmt.Println("sharesCoeffs: ", sharesCoeffs)

	// Produce secret-shared elements
	skShares := make([]*ring.Poly, k)
	for l := 0; l < k; l++ {
		for i := 0; i < n; i++ {
			ringElem := r.NewPoly()
			r.SetCoefficientsBigint(sharesCoeffs[l][i], ringElem)
			skShares[l] = &ringElem
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

	return skShares, seeds, b
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
func SignRound1(r *ring.Ring, uniformSampler *ring.UniformSampler, partyInt int, sid int, skShare *ring.Poly, mu string, PRFKey []byte, seeds_i [][]byte, T []int) (*ring.Poly, *ring.Poly) {
	// Generate the row-wise mask from PRFs
	m_i := make([]*ring.Poly, n)
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

	// Hash the secret key to seed the Gaussian
	hasher := sha3.NewShake128()
	skMarshalled, err := skShare.MarshalBinary()
	if err != nil {
		log.Fatalf("Error marshalling poly: %v\n", err)
	}
	hasher.Write(skMarshalled)
	hashOutputLength := n // TODO: Check
	skHash := make([]byte, hashOutputLength)
	_, err = hasher.Read(skHash)
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

	return nil, nil
}

// Sign function signs a message using the secret key and returns the signature
func SignRound2(partyInt int, D []*ring.Poly, m []*ring.Poly) (c_i *ring.Poly, z_i *ring.Poly) {
	return nil, nil
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
func H_u(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, sid int, j int, D []*ring.Poly, mu string, m []*ring.Poly) []*ring.Poly {
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

	// Handle slice D
	for _, poly := range D {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	// Handle string mu
	buffer.WriteString(mu)

	// Handle slice m
	for _, poly := range m {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
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
