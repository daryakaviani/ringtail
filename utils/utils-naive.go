package utils

import (
	"log"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

const logN = 8

func MulPolyNaive(r *ring.Ring, p1 ring.Poly, p2 ring.Poly, p3 ring.Poly) {
	degree := 1 << logN // Since we are in a ring modulo x^256 + 1

	q := r.Modulus()
	// Initialize result slice with big.Ints set to zero
	result := make([]*big.Int, degree)
	for i := range result {
		result[i] = big.NewInt(0)
	}

	p1Coeffs := make([]*big.Int, degree)
	r.PolyToBigint(p1, 1, p1Coeffs)

	p2Coeffs := make([]*big.Int, degree)
	r.PolyToBigint(p2, 1, p2Coeffs)

	// Polynomial multiplication (convolution)
	for i := range p1Coeffs {
		for j := range p2Coeffs {
			coeff1 := p1Coeffs[i]
			coeff2 := p2Coeffs[j]
			temp := new(big.Int).Mul(coeff1, coeff2)
			if i+j < degree {
				// Multiply coefficients and add to the right place
				result[i+j].Add(result[i+j], temp)
			} else {
				index := (i + j) % degree
				if (((i+j)-((i+j)%degree))/(2<<logN))%2 == 0 {
					// Wrap around due to cyclotomic ring, i+j >= degree
					result[index].Add(result[index], temp)
				} else {
					// Wrap around due to cyclotomic ring, i+j >= degree
					result[index].Sub(result[index], temp) // subtracting because of x^256 + 1
				}
			}
		}
	}

	// Reduce each coefficient modulo q
	for i := range result {
		result[i].Mod(result[i], q)
	}

	r.SetCoefficientsBigint(result, p3)
}

// Sets p3 = p1 * p2 by first converting into NTT, multiplying coefficient-wise, and then converting out of NTT with INTT
func MulPolyNTT(r *ring.Ring, p1 ring.Poly, p2 ring.Poly, p3 ring.Poly) {
	// Transform p1 and p2 to the NTT domain
	r.NTT(p1, p1)
	r.NTT(p2, p2)

	p1Coeffs := make([]*big.Int, r.N())
	p2Coeffs := make([]*big.Int, r.N())
	p3Coeffs := make([]*big.Int, r.N())

	r.PolyToBigint(p1, 1, p1Coeffs)
	r.PolyToBigint(p2, 1, p2Coeffs)

	// Perform coefficient-wise multiplication in the NTT domain
	for i := 0; i < r.N(); i++ {
		p3Coeffs[i] = new(big.Int).Mul(p1Coeffs[i], p2Coeffs[i])
		p3Coeffs[i] = new(big.Int).Mod(p3Coeffs[i], r.Modulus())
	}

	r.SetCoefficientsBigint(p3Coeffs, p3)

	// Transform the result back to the standard domain using INTT
	r.INTT(p3, p3)
	r.INTT(p1, p1)
	r.INTT(p2, p2)
}

// MulCoeffsNTT performs coefficient-wise multiplication of two polynomials in the NTT domain.
func MulCoeffsNTT(r *ring.Ring, p1, p2, result ring.Poly) {
	p1Coeffs := make([]*big.Int, r.N())
	p2Coeffs := make([]*big.Int, r.N())
	p3Coeffs := make([]*big.Int, r.N())

	r.PolyToBigint(p1, 1, p1Coeffs)
	r.PolyToBigint(p2, 1, p2Coeffs)

	for i := 0; i < r.N(); i++ {
		p3Coeffs[i] = new(big.Int).Mul(p1Coeffs[i], p2Coeffs[i])
		p3Coeffs[i] = new(big.Int).Mod(p3Coeffs[i], r.Modulus())
	}

	r.SetCoefficientsBigint(p3Coeffs, result)
}

// NAIVE MULTIPLICATIONS

// MatrixVectorMul performs matrix-vector multiplication.
// It takes a matrix of ring.Poly pointers, a vector of ring.Poly pointers, and outputs the result in a given result vector.
func MatrixVectorMulNaive(r *ring.Ring, M structs.Matrix[ring.Poly], vec structs.Vector[ring.Poly], result structs.Vector[ring.Poly]) {
	for i := range M {
		result[i] = r.NewPoly()
		for j := range M[i] {
			temp := r.NewPoly()
			MulPolyNaive(r, M[i][j], vec[j], temp)
			r.Add(result[i], temp, result[i])
		}
	}
}

// MatrixMatrixMul performs matrix-matrix multiplication.
// It takes two matrices of ring.Poly pointers, M1 of dimensions m x p and M2 of dimensions p x n,
// and outputs the result in a given result matrix of dimensions m x n.
func MatrixMatrixMulNaive(r *ring.Ring, M1, M2 structs.Matrix[ring.Poly], result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len((M1)[0]) != len(M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(M1)
	p := len(M1[0])
	n := len(M2[0])

	// Initialize the result matrix with zeros
	for i := 0; i < m; i++ {
		result[i] = make([]ring.Poly, n)
		for j := 0; j < n; j++ {
			result[i][j] = r.NewPoly()
		}
	}

	// Perform matrix multiplication
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				temp := r.NewPoly()
				MulPolyNaive(r, M1[i][k], M2[k][j], temp)
				r.Add(result[i][j], temp, result[i][j])
			}
		}
	}
}

// VectorPolyMulNaive performs element-wise multiplication of a vector by a polynomial.
// It takes a vector of ring.Poly pointers, a single ring.Poly pointer, and outputs the result in a given result vector.
func VectorPolyMulNaive(r *ring.Ring, vec structs.Vector[ring.Poly], poly ring.Poly, result structs.Vector[ring.Poly]) {
	for i := range vec {
		MulPolyNaive(r, vec[i], poly, result[i]) // Multiply each vector element by the polynomial
	}
}
