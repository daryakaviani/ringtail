package utils

import (
	"fmt"
	"log"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

// RoundCoeffsToNearestMultiple rounds the coefficients of a polynomial to the nearest multiple of p
func RoundCoeffsToNearestMultiple(r *ring.Ring, poly ring.Poly, p, q uint64) {
	qBig := new(big.Int).SetUint64(q)
	halfQ := new(big.Int).Div(qBig, big.NewInt(2))
	coeffsBigint := make([]*big.Int, r.N())
	roundedCoeffs := make([]*big.Int, r.N())

	r.PolyToBigint(poly, 1, coeffsBigint)

	for i, coeff := range coeffsBigint {
		if roundedCoeffs[i] == nil {
			roundedCoeffs[i] = new(big.Int)
		}

		// Convert to signed representation
		if coeff.Cmp(halfQ) > 0 {
			coeff.Sub(coeff, qBig)
		}

		// Scale the coefficient
		scaledCoeff := new(big.Float).Quo(new(big.Float).SetInt(coeff), new(big.Float).SetUint64(q))
		scaledCoeff.Mul(scaledCoeff, new(big.Float).SetUint64(p))

		// Round to the nearest integer
		roundedCoeff, _ := scaledCoeff.Int(nil)
		if scaledCoeff.Sign() >= 0 {
			scaledCoeff.Sub(scaledCoeff, new(big.Float).SetInt(roundedCoeff))
			if scaledCoeff.Cmp(new(big.Float).SetFloat64(0.5)) >= 0 {
				roundedCoeff.Add(roundedCoeff, big.NewInt(1))
			}
		} else {
			scaledCoeff.Sub(scaledCoeff, new(big.Float).SetInt(roundedCoeff))
			if scaledCoeff.Cmp(new(big.Float).SetFloat64(-0.5)) <= 0 {
				roundedCoeff.Sub(roundedCoeff, big.NewInt(1))
			}
		}

		// Map back to Z_q
		roundedCoeffs[i].Mod(roundedCoeff, qBig)
		if roundedCoeffs[i].Cmp(big.NewInt(0)) < 0 {
			roundedCoeffs[i].Add(roundedCoeffs[i], qBig)
		}
	}

	r.SetCoefficientsBigint(roundedCoeffs, poly)
}

// MatrixVectorMul performs matrix-vector multiplication.
func MatrixVectorMul(r *ring.Ring, M structs.Matrix[ring.Poly], vec structs.Vector[ring.Poly], result structs.Vector[ring.Poly]) {
	// Convert all elements of the matrix and the vector to the NTT domain
	ConvertMatrixToNTT(r, M)
	ConvertVectorToNTT(r, vec)

	// Perform the multiplications coefficient-wise
	temp := r.NewPoly()

	for i := range M {
		result[i] = r.NewPoly()
		for j := range (M)[i] {
			r.MForm(M[i][j], temp)
			r.MulCoeffsMontgomeryThenAdd(temp, vec[j], result[i])
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertMatrixFromNTT(r, M)
	ConvertVectorFromNTT(r, vec)
}

// MatrixMatrixMul performs matrix-matrix multiplication.
func MatrixMatrixMul(r *ring.Ring, M1, M2 structs.Matrix[ring.Poly], result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len((M1)[0]) != len(M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(M1)
	p := len(M1[0]) // Assuming all rows in M1 are of the same length
	n := len(M2[0]) // Assuming all rows in M2 are of the same length

	// Convert all elements of M1 and M2 to the NTT domain
	ConvertMatrixToNTT(r, M1)
	ConvertMatrixToNTT(r, M2)

	// Initialize the result matrix with zeros
	for i := 0; i < m; i++ {
		result[i] = make([]ring.Poly, n)
		for j := 0; j < n; j++ {
			result[i][j] = r.NewPoly()
		}
	}

	temp := r.NewPoly()

	// Perform matrix multiplication coefficient-wise
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				r.MForm(M1[i][k], temp)
				r.MulCoeffsMontgomeryThenAdd(temp, M2[k][j], result[i][j])
			}
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertMatrixFromNTT(r, result)
	ConvertMatrixFromNTT(r, M1)
	ConvertMatrixFromNTT(r, M2)
}

// VectorPolyMul performs element-wise multiplication of a vector by a polynomial.
func VectorPolyMul(r *ring.Ring, vec structs.Vector[ring.Poly], poly ring.Poly, result structs.Vector[ring.Poly]) {
	// Convert the polynomial to the NTT domain
	r.NTT(poly, poly)

	// Convert all elements of the vector to the NTT domain
	ConvertVectorToNTT(r, vec)
	ConvertVectorToNTT(r, result)
	temp := r.NewPoly()

	// Perform the multiplications coefficient-wise
	for i := range vec {
		r.MForm(vec[i], temp)
		r.MulCoeffsMontgomery(temp, poly, result[i])
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertVectorFromNTT(r, vec)
	r.INTT(poly, poly)
}

// MatrixAdd adds two matrices of ring.Poly element-wise and stores the result in a given result matrix.
func MatrixAdd(r *ring.Ring, M1, M2, result structs.Matrix[ring.Poly]) {
	if M1 == nil || M2 == nil || len(M1) == 0 || len(M2) == 0 || len(M1) != len(M2) || len((M1)[0]) != len((M2)[0]) {
		log.Fatalf("Matrix dimensions must match for element-wise addition.")
		return
	}

	m := len(M1)
	n := len(M1[0])

	for i := 0; i < m; i++ {
		if result[i] == nil {
			result[i] = make([]ring.Poly, n)
		}
		for j := 0; j < n; j++ {
			r.Add(M1[i][j], M2[i][j], result[i][j])
		}
	}
}

// VectorAdd adds two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorAdd(r *ring.Ring, v1, v2, result structs.Vector[ring.Poly]) {
	for i := range v1 {
		r.Add(v1[i], v2[i], result[i])
	}
}

// VectorSub subtracts two vectors of ring.Poly element-wise and stores the result in a result vector.
func VectorSub(r *ring.Ring, v1, v2, result structs.Vector[ring.Poly]) {
	for i := range v1 {
		r.Sub(v1[i], v2[i], result[i])
	}
}

// SAMPLER HELPERS

// SamplePolyVector samples a vector of polynomials of a given length using the provided sampler.
func SamplePolyVector(length int, sampler ring.Sampler) structs.Vector[ring.Poly] {
	vector := structs.Vector[ring.Poly](make([]ring.Poly, length))
	for i := 0; i < length; i++ {
		element := sampler.ReadNew()
		vector[i] = element
	}
	return vector
}

// SamplePolyMatrix samples a matrix of polynomials with given dimensions (rows and cols) using the provided sampler.
func SamplePolyMatrix(rows, cols int, sampler ring.Sampler) structs.Matrix[ring.Poly] {
	matrix := structs.Matrix[ring.Poly](make([][]ring.Poly, rows))
	for i := 0; i < rows; i++ {
		matrix[i] = make([]ring.Poly, cols)
		for j := 0; j < cols; j++ {
			element := sampler.ReadNew()
			matrix[i][j] = element
		}
	}
	return matrix
}

// PRINT FUNCTIONS

func PrintMatrix(label string, matrix structs.Matrix[ring.Poly]) {
	log.Println(label)
	for i, row := range matrix {
		for j, poly := range row {
			log.Printf("[%d][%d]: %s\n", i, j, formatCoeffs(poly.Coeffs[0]))
		}
	}
}

func PrintVector(label string, vector structs.Vector[ring.Poly]) {
	log.Println(label)
	for i, poly := range vector {
		log.Printf("[%d]: %s\n", i, formatCoeffs(poly.Coeffs[0]))
	}
}

func PrintPolynomial(label string, poly ring.Poly) {
	log.Println(label)
	log.Printf("%s\n", formatCoeffs(poly.Coeffs[0]))
}

func formatCoeffs(coeffs []uint64) string {
	var coeffStr []string
	for _, coeff := range coeffs {
		coeffStr = append(coeffStr, fmt.Sprintf("%v", coeff))
	}
	return strings.Join(coeffStr, ", ")
}

// NTT CONVERSION

// ConvertMatrixToNTT converts a matrix of polynomials to the NTT domain.
func ConvertMatrixToNTT(r *ring.Ring, M structs.Matrix[ring.Poly]) {
	for i := range M {
		for j := range M[i] {
			r.NTT(M[i][j], M[i][j])
		}
	}
}

// ConvertMatrixFromNTT converts a matrix of polynomials from the NTT domain back to the standard domain.
func ConvertMatrixFromNTT(r *ring.Ring, M structs.Matrix[ring.Poly]) {
	for i := range M {
		for j := range M[i] {
			r.INTT(M[i][j], M[i][j])
		}
	}
}

// ConvertVectorToNTT converts a vector of polynomials to the NTT domain.
func ConvertVectorToNTT(r *ring.Ring, vec structs.Vector[ring.Poly]) {
	for i := range vec {
		r.NTT(vec[i], vec[i])
	}
}

// ConvertVectorFromNTT converts a vector of polynomials from the NTT domain back to the standard domain.
func ConvertVectorFromNTT(r *ring.Ring, vec structs.Vector[ring.Poly]) {
	for i := range vec {
		r.INTT(vec[i], vec[i])
	}
}

// INITIALIZE HELPERS

// InitializeVector creates and returns a vector of the given length, initializing each element as a new polynomial.
func InitializeVector(r *ring.Ring, length int) structs.Vector[ring.Poly] {
	vector := make(structs.Vector[ring.Poly], length)
	for i := range vector {
		vector[i] = r.NewPoly()
	}
	return vector
}

// InitializeMatrix creates and returns a matrix of the given dimensions, initializing each element as a new polynomial.
func InitializeMatrix(r *ring.Ring, rows, cols int) structs.Matrix[ring.Poly] {
	matrix := make(structs.Matrix[ring.Poly], rows)
	for i := range matrix {
		matrix[i] = make([]ring.Poly, cols)
		for j := range matrix[i] {
			matrix[i][j] = r.NewPoly()
		}
	}
	return matrix
}

// Copy map helpers

func CopyMatrixMap(original map[int]structs.Matrix[ring.Poly]) map[int]structs.Matrix[ring.Poly] {
	copy := make(map[int]structs.Matrix[ring.Poly])
	for key, value := range original {
		copy[key] = value
	}
	return copy
}

func CopyVectorMap(original map[int]structs.Vector[ring.Poly]) map[int]structs.Vector[ring.Poly] {
	copy := make(map[int]structs.Vector[ring.Poly])
	for key, value := range original {
		copy[key] = value
	}
	return copy
}
