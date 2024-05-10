package main

import (
	"bytes"
	"encoding/binary"
	"log"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

const (
	m = 20
	n = 20
	// q       = uint64(61) // Modulo parameter
	p       = 5
	e_bound = 0 // Keep these as 0 for now!
	sigma_e = 0 // Standard deviation for the error distribution
	sigma_c = 2 // Standard deviation for the hash output distribution
	c_bound = 3
	logN    = 3
	kappa   = 10
)

var q = ring.Qi60[0]

func main() {
	// Initialize ring and samplers
	r, _ := ring.NewRing(1<<logN, []uint64{q})
	prng, _ := sampling.NewKeyedPRNG([]byte("0"))
	uniformSampler := ring.NewUniformSampler(prng, r)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: e_bound}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

	// Setup
	A := BCMSetup(r, uniformSampler)
	PrintMatrix("Matrix A:", A)

	// Key Generation
	s, b := BCMGen(r, A, uniformSampler, gaussianSampler)
	PrintVector("Vector s:", s)
	PrintVector("Vector b:", b)

	// Signing a message
	mu := "Hello, BCM!"
	c, z := BCMSign(r, uniformSampler, A, s, b, mu)
	PrintPolynomial("Polynomial c:", c)
	PrintVector("Vector z:", z)

	// Verify the signature
	valid := BCMVer(r, A, b, mu, c, z)
	if valid {
		log.Println("Signature verification successful.")
	} else {
		log.Println("Signature verification failed.")
	}
}

func PrintMatrix(label string, matrix *[][]*ring.Poly) {
	log.Println(label)
	for i, row := range *matrix {
		for j, poly := range row {
			log.Printf("A[%d][%d]: %v\n", i, j, poly.Coeffs[0])
		}
	}
}

func PrintVector(label string, vector []*ring.Poly) {
	log.Println(label)
	for i, poly := range vector {
		log.Printf("[%d]: %v\n", i, poly.Coeffs[0])
	}
}

func PrintPolynomial(label string, poly *ring.Poly) {
	log.Println(label)
	log.Printf("%v\n", poly.Coeffs[0])
}

// BCMSetup generates public parameters.
func BCMSetup(r *ring.Ring, uniformSampler *ring.UniformSampler) *[][]*ring.Poly {
	A := make([][]*ring.Poly, m)
	for i := 0; i < m; i++ {
		A[i] = make([]*ring.Poly, n)
		for j := 0; j < n; j++ {
			newPoly := uniformSampler.ReadNew()
			A[i][j] = &newPoly
		}
	}
	return &A
}

// BCMGen generates a secret key and computes the public vector b.
func BCMGen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, gaussianSampler *ring.GaussianSampler) ([]*ring.Poly, []*ring.Poly) {
	s := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		newPoly := uniformSampler.ReadNew()
		s[i] = &newPoly
	}

	e := make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		newPoly := gaussianSampler.ReadNew()
		e[i] = &newPoly
	}

	b := make([]*ring.Poly, m)
	MatrixVectorMul(r, A, s, b)
	VectorAdd(r, b, e, b)

	return s, b
}

// BCMSign generates a signature for the message mu.
func BCMSign(r *ring.Ring, uniformSampler *ring.UniformSampler, A *[][]*ring.Poly, s []*ring.Poly, b []*ring.Poly, mu string) (*ring.Poly, []*ring.Poly) {
	rVec := make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		rVec[i] = r.NewPoly().CopyNew()
		*rVec[i] = uniformSampler.ReadNew()
	}

	// Compute h using the new helper function
	h := make([]*ring.Poly, m)
	MatrixVectorMul(r, A, rVec, h)
	PrintVector("Original h: ", h)

	// Round h to the nearest multiple of p
	for _, poly := range h {
		RoundCoeffsToNearestMultiple(r, poly, p)
	}
	PrintVector("Rounded original h: ", h)

	// Hash h and mu to obtain c
	c := H(r, A, b, h, mu)

	// Compute z using the new helper function
	z := make([]*ring.Poly, n)
	VectorPolyMul(r, s, c, z)
	VectorAdd(r, z, rVec, z)

	return c, z
}

// Hash to Gaussian
func H(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, h []*ring.Poly, mu string) *ring.Poly {
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
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_c, Bound: c_bound}
	hashGaussiamSampler := ring.NewGaussianSampler(prng, r, gaussianParams, true)

	c := hashGaussiamSampler.ReadNew()

	return &c
}

// BCMVer verifies the signature using public parameters and the message.
func BCMVer(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, mu string, c *ring.Poly, z []*ring.Poly) bool {
	// Compute Az using MatrixVectorMul
	Az := make([]*ring.Poly, m)
	MatrixVectorMul(r, A, z, Az)

	// Compute bc using VectorPolyMul
	bc := make([]*ring.Poly, m)
	VectorPolyMul(r, b, c, bc)

	// Subtract bc from Az to get Az_bc
	Az_bc := make([]*ring.Poly, m)
	VectorSub(r, Az, bc, Az_bc)

	PrintVector("Az - bc: ", Az_bc)

	// Round Az_bc to the nearest multiple of p
	for _, poly := range Az_bc {
		RoundCoeffsToNearestMultiple(r, poly, p)
	}
	PrintVector("Rounded Az_bc: ", Az_bc)

	// Hash Az_bc and mu to obtain a hash value
	computedC := H(r, A, b, Az_bc, mu)

	PrintPolynomial("Computed c: ", computedC)

	// Compare computed c with the provided c
	return r.Equal(*computedC, *c)
}

// TODO: Make this general to other rings which are larger
// polyMultInCyclotomicRing multiplies two polynomials within the cyclotomic ring x^8 + 1 and returns the result.
// Each polynomial is represented as a slice of big.Int pointers, sorted from least to most significant.
func MulPoly(r *ring.Ring, p1 *ring.Poly, p2 *ring.Poly, p3 *ring.Poly) {
	degree := 1 << logN // Since we are in a ring modulo x^8 + 1
	// Initialize result slice with big.Ints set to zero
	result := make([]*big.Int, degree)
	for i := range result {
		result[i] = big.NewInt(0)
	}

	p1Coeffs := make([]*big.Int, degree)
	r.PolyToBigint(*p1, 1, p1Coeffs)

	p2Coeffs := make([]*big.Int, degree)
	r.PolyToBigint(*p2, 1, p2Coeffs)

	// Polynomial multiplication (convolution)
	for i := range p1Coeffs {
		for j := range p2Coeffs {
			if i+j < degree {
				// Multiply coefficients and add to the right place
				temp := new(big.Int).Mul(p1Coeffs[i], p2Coeffs[j])
				result[i+j].Add(result[i+j], temp)
			} else {
				if (((i+j)-((i+j)%degree))/8)%2 == 0 {
					// Wrap around due to cyclotomic ring, i+j >= degree
					temp := new(big.Int).Mul(p1Coeffs[i], p2Coeffs[j])
					result[(i+j)%degree].Add(result[(i+j)%degree], temp)
				} else {
					// Wrap around due to cyclotomic ring, i+j >= degree
					temp := new(big.Int).Mul(p1Coeffs[i], p2Coeffs[j])
					result[(i+j)%degree].Sub(result[(i+j)%degree], temp) // subtracting because of x^8 + 1
				}
			}
		}
	}

	// Reduce each coefficient modulo q
	for i := range result {
		result[i].Mod(result[i], big.NewInt(int64(q)))
	}

	r.SetCoefficientsBigint(result, *p3)
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

// MatrixVectorMul performs matrix-vector multiplication.
// It takes a matrix of ring.Poly pointers, a vector of ring.Poly pointers, and outputs the result in a given result vector.
func MatrixVectorMul(r *ring.Ring, M *[][]*ring.Poly, vec []*ring.Poly, result []*ring.Poly) {
	for i := range *M {
		result[i] = r.NewPoly().CopyNew()
		for j := range (*M)[i] {
			temp := r.NewPoly()
			MulPoly(r, (*M)[i][j], vec[j], &temp)
			r.Add(*result[i], temp, *result[i]) // Accumulate the result
		}
	}
}

// VectorPolyMul performs element-wise multiplication of a vector by a polynomial.
// It takes a vector of ring.Poly pointers, a single ring.Poly pointer, and outputs the result in a given result vector.
func VectorPolyMul(r *ring.Ring, vec []*ring.Poly, poly *ring.Poly, result []*ring.Poly) {
	for i := range vec {
		result[i] = r.NewPoly().CopyNew()   // Initialize a new polynomial for each result entry
		MulPoly(r, vec[i], poly, result[i]) // Multiply each vector element by the polynomial
	}
}

// VectorAdd adds two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorAdd(r *ring.Ring, v1, v2, result []*ring.Poly) {
	for i := range v1 {
		if result[i] == nil {
			newPoly := r.NewPoly()
			result[i] = &newPoly
		}

		r.Add(*v1[i], *v2[i], *result[i])
	}
}

// VectorAdd subtracts two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorSub(r *ring.Ring, v1, v2, result []*ring.Poly) {
	for i := range v1 {
		if result[i] == nil {
			newPoly := r.NewPoly()
			result[i] = &newPoly
		}

		r.Sub(*v1[i], *v2[i], *result[i])
	}
}
