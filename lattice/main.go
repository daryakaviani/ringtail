// lattice/main.go

package main

import (
	"crypto/rand"
	"fmt"
	"log"
	"math/big"

	"github.com/hashicorp/vault/shamir"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
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
	gaussianParams := ring.DiscreteGaussian{}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

	// Setup
	A := Setup(uniformSampler)

	fmt.Println("A: ", A)

	// Gen
	skShares, seeds, b := Gen(r, A, uniformSampler, gaussianSampler)

	fmt.Println("Shares: ", skShares)
	fmt.Println("Seeds: ", seeds)

	// Test signing round 1 by simulating each of the parties
	mu := "Hello, Threshold Signature!"
	sid := 1
	const numActiveParties = 4
	D := make([]*ring.Poly, numActiveParties)
	m := make([]*ring.Poly, numActiveParties)

	for i := 0; i < numActiveParties; i++ {
		D[i], m[i] = SignRound1(i, sid, skShares[i], mu)
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
func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, gaussianSampler *ring.GaussianSampler) ([]*ring.Poly, [][][]byte, []*ring.Poly) {
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
	seeds := make([][][]byte, d)
	for i := 0; i < d; i++ {
		seeds[i] = make([][]byte, d) // Initialize each second-level slice within the outer slice
		for j := 0; j < d; j++ {
			// Generate random seed sd
			seeds[i][j] = generateRandomSeed()
		}
	}

	return skShares, seeds, b
}

func generateRandomSeed() []byte {
	// Generate a random binary string of length ell
	sd := make([]byte, (ell+7)/8)
	_, err := rand.Read(sd)
	if err != nil {
		panic(err)
	}
	return sd
}

// Sign function signs a message using the secret key and returns the signature
func SignRound1(partyInt int, sid int, skShare *ring.Poly, mu string) (D_i *ring.Poly, m_i *ring.Poly) {
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

// // Verify function verifies the signature of a message using the public key
// func Verify(context *lattigo.Context, pk *lattigo.PublicKey, message string, signature *lattigo.Signature) bool {
// 	messageBytes := []byte(message)
// 	valid := pk.Verify(signature, messageBytes)
// 	return valid
// }
