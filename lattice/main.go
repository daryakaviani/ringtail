// lattice/main.go

package main

import (
	"crypto/rand"
	"math/big"

	"github.com/hashicorp/vault/shamir"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

const m = 5
const n = 5
const q = 5
const t = 3
const k = 5
const d = 5
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

	// Gen
	shares, seeds := Gen(r, A, uniformSampler, gaussianSampler)

	// Sign
	// message := "Hello, Threshold Signature!"
	// signature := Sign(context, sk, message)

	// Verify
	// valid := Verify(context, pk, message, signature)
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
func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, gaussianSampler *ring.GaussianSampler) ([]*ring.Poly, [][][]byte) {
	// Sample the secret key from the ring
	s := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		// Generate a new ring element
		element := uniformSampler.ReadNew()
		s[i] = &element
	}

	// Sample the error from the Gaussian distribution
	e := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		// Generate a new ring element from the Gaussian distribution
		element := gaussianSampler.ReadNew()
		e[i] = &element
	}

	// Initialize result vector
	b := make([]*ring.Poly, m)

	// Compute b = As + e mod q
	for i := 0; i < m; i++ {
		// Initialize b[i] as a zero polynomial
		*b[i] = r.NewPoly()

		// Compute A*s for row i
		for j := 0; j < n; j++ {
			// Compute A[i][j] * s[j] mod q and add it to b[i]
			r.MulCoeffsMontgomeryThenAdd(*(*A)[i][j], *s[j], *b[i])
		}

		// Add e[i] to b[i] mod q
		r.Add(*b[i], *e[i], *b[i])
	}

	// Secret-share s
	sharesCoeffs := make([][][]*big.Int, k)

	// Indexing over all ring elements of the s vector
	for i := 0; i < n; i++ {
		// These are the coefficients of the ring element s[i]
		coeffsBigint := make([]*big.Int, r.N())
		r.PolyToBigint(*s[i], 1, coeffsBigint)

		// For each coefficient of the ring element polynomial
		for j, coeff := range coeffsBigint {
			coeffSharesBytes, _ := shamir.Split(coeff.Bytes(), k, t)

			// For each share of this coefficient of the ring element polynomial
			for l := 0; l < k; l++ {
				coeffShare := new(big.Int)
				coeffShare.SetBytes(coeffSharesBytes[l])

				// Set the coefficient share at the k-th secret share, i-th secret vector element, j-th coefficient
				sharesCoeffs[l][i][j] = coeffShare
			}
		}
	}

	// Produce secret-shared elements
	secretKeyShares := make([]*ring.Poly, k)
	for l := 0; l < k; l++ {
		for i := 0; i < n; i++ {
			ringElem := r.NewPoly()
			r.SetCoefficientsBigint(sharesCoeffs[l][i], ringElem)
			secretKeyShares[l] = &ringElem
		}
	}

	// Generate random seeds for i, j in [d]
	seeds := make([][][]byte, d)
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			// Generate random seed sd
			seeds[i][j] = generateRandomSeed()
		}
	}

	return secretKeyShares, seeds
}

func generateRandomSeed() []byte {
	// Generate a random binary string of length ell
	sd := make([]byte, (ell+7)/8)
	_, err := rand.Read(sd)
	if err != nil {
		panic(err) // Handle error appropriately
	}
	return sd
}

// r.SetCoefficientsBigint()

// // Sign function signs a message using the secret key and returns the signature
// func Sign(context *lattigo.Context, sk *lattigo.SecretKey, message string) *lattigo.Signature {
// 	messageBytes := []byte(message)
// 	signature := sk.Sign(messageBytes)
// 	return signature
// }

// // Verify function verifies the signature of a message using the public key
// func Verify(context *lattigo.Context, pk *lattigo.PublicKey, message string, signature *lattigo.Signature) bool {
// 	messageBytes := []byte(message)
// 	valid := pk.Verify(signature, messageBytes)
// 	return valid
// }
