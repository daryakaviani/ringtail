package utils

import (
	"math/big"
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
	"gonum.org/v1/gonum/mat"
)

const Q = 0x1000000004A01

// Multiply a matrix by its conjugate transpose
func MultiplyByConjugateTranspose(matrix structs.Matrix[complex128]) structs.Matrix[complex128] {
	rows := len(matrix)
	cols := len(matrix[0])

	conjTranspose := make(structs.Matrix[complex128], cols)
	for i := 0; i < cols; i++ {
		conjTranspose[i] = make([]complex128, rows)
		for j := 0; j < rows; j++ {
			conjTranspose[i][j] = cmplx.Conj(matrix[j][i])
		}
	}

	result := make(structs.Matrix[complex128], rows)
	for i := 0; i < rows; i++ {
		result[i] = make([]complex128, rows)
		for j := 0; j < rows; j++ {
			var sum complex128
			for k := 0; k < cols; k++ {
				sum += matrix[i][k] * conjTranspose[k][j]
			}
			result[i][j] = sum
		}
	}

	return result
}

// Convert a complex matrix to a real matrix
func ConvertComplexToRealMatrix(complexMatrix structs.Matrix[complex128]) *mat.Dense {
	rows := len(complexMatrix)
	cols := len(complexMatrix[0])
	realMatrix := mat.NewDense(2*rows, 2*cols, nil)

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			c := complexMatrix[i][j]
			realMatrix.Set(2*i, 2*j, real(c))
			realMatrix.Set(2*i, 2*j+1, -imag(c))
			realMatrix.Set(2*i+1, 2*j, imag(c))
			realMatrix.Set(2*i+1, 2*j+1, real(c))
		}
	}

	return realMatrix
}

// CheckRowCoprime checks if the entries of a row are setwise coprime
func CheckRowCoprime(r *ring.Ring, row []ring.Poly) bool {
	if len(row) == 0 {
		return true
	}
	gcd := CalculateAlgebraicNorm(r, row[0])
	for i := 1; i < len(row); i++ {
		norm := CalculateAlgebraicNorm(r, row[i])
		newGCD := new(big.Int).GCD(nil, nil, gcd, norm)
		if newGCD.Cmp(big.NewInt(1)) == 0 {
			return true
		}
		gcd = newGCD
	}
	return gcd.Cmp(big.NewInt(1)) == 0
}

// Calculate the algebraic norm using Karatsuba multiplication
func CalculateAlgebraicNorm(r *ring.Ring, poly ring.Poly) *big.Int {
	coeffs := make([]*big.Int, r.N())
	r.PolyToBigint(poly, 1, coeffs)
	SignedRepresentation(coeffs, Q)
	return KaratsubaNorm(coeffs)
}
