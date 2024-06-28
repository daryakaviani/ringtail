package utils

import (
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v5/utils/structs"
	"gonum.org/v1/gonum/mat"
)

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
