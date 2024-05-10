package main

import (
	"bytes"
	"encoding/binary"
	"lattice-threshold-signature/utils"
	"log"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

const (
	m = 5
	n = 5
	// q       = uint64(61) // Modulo parameter
	p       = 100
	e_bound = 0.5 // Keep these as 0 for now!
	sigma_e = 2   // Standard deviation for the error distribution
	sigma_c = 2   // Standard deviation for the hash output distribution
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
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	// Setup
	A := BCMSetup(r, uniformSampler)
	utils.PrintMatrix("Matrix A:", A)

	// Key Generation
	s, b := BCMGen(r, A, uniformSampler, gaussianSampler)
	utils.PrintVector("Vector s:", s)
	utils.PrintVector("Vector b:", b)

	// Signing a message
	mu := "Hello, BCM!"
	c, z := BCMSign(r, uniformSampler, A, s, b, mu)
	utils.PrintPolynomial("Polynomial c:", c)
	utils.PrintVector("Vector z:", z)

	// Verify the signature
	valid := BCMVer(r, A, b, mu, c, z)
	if valid {
		log.Println("Signature verification successful.")
	} else {
		log.Println("Signature verification failed.")
	}
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

	utils.PrintVector("e: ", e)

	b := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

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
	utils.MatrixVectorMul(r, A, rVec, h)
	utils.PrintVector("Original h: ", h)

	// Round h to the nearest multiple of p
	for _, poly := range h {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	utils.PrintVector("Rounded original h: ", h)

	// Hash h and mu to obtain c
	c := H(r, A, b, h, mu)

	// Compute z using the new helper function
	z := make([]*ring.Poly, n)
	utils.VectorPolyMul(r, s, c, z)
	utils.VectorAdd(r, z, rVec, z)

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
	hashGaussiamSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	c := hashGaussiamSampler.ReadNew()

	return &c
}

// BCMVer verifies the signature using public parameters and the message.
func BCMVer(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, mu string, c *ring.Poly, z []*ring.Poly) bool {
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

	// Hash Az_bc and mu to obtain a hash value
	computedC := H(r, A, b, Az_bc, mu)

	utils.PrintPolynomial("Computed c: ", computedC)

	// Compare computed c with the provided c
	return r.Equal(*computedC, *c)
}
