package utils

import (
	"log"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
)

const logN = 8

// TODO: Make this general to other rings which are larger
// polyMultInCyclotomicRing multiplies two polynomials within the cyclotomic ring x^8 + 1 and returns the result.
// Each polynomial is represented as a slice of big.Int pointers, sorted from least to most significant.
func MulPolyNaive(r *ring.Ring, p1 *ring.Poly, p2 *ring.Poly, p3 *ring.Poly) {
	degree := 1 << logN // Since we are in a ring modulo x^8 + 1

	q := r.Modulus()
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
				if (((i+j)-((i+j)%degree))/(2<<logN))%2 == 0 {
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
		result[i].Mod(result[i], q)
	}

	r.SetCoefficientsBigint(result, *p3)
}

// Sets p3 = p1 * p2 by first converting into NTT, multiplying coefficient-wise, and then converting out of NTT with INTT
func MulPolyNTT(r *ring.Ring, p1 *ring.Poly, p2 *ring.Poly, p3 *ring.Poly) {
	// Transform p1 and p2 to the NTT domain
	r.NTT(*p1, *p1)
	r.NTT(*p2, *p2)

	PrintPolynomial("p1:", p1)
	PrintPolynomial("p2:", p2)
	p1Coeffs := make([]*big.Int, r.N())
	p2Coeffs := make([]*big.Int, r.N())
	p3Coeffs := make([]*big.Int, r.N())

	r.PolyToBigint(*p1, 1, p1Coeffs)
	r.PolyToBigint(*p2, 1, p2Coeffs)

	// Perform coefficient-wise multiplication in the NTT domain
	for i := 0; i < r.N(); i++ {
		p3Coeffs[i] = new(big.Int).Mul(p1Coeffs[i], p2Coeffs[i])
		p3Coeffs[i] = new(big.Int).Mod(p3Coeffs[i], r.Modulus())
	}

	r.SetCoefficientsBigint(p3Coeffs, *p3)

	// Transform the result back to the standard domain using INTT
	r.INTT(*p3, *p3)
	r.INTT(*p1, *p1)
	r.INTT(*p2, *p2)
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
	// Convert all elements of the matrix and the vector to the NTT domain
	ConvertMatrixToNTT(r, M)
	ConvertVectorToNTT(r, vec)

	// Perform the multiplications coefficient-wise
	for i := range *M {
		result[i] = r.NewPoly().CopyNew()
		for j := range (*M)[i] {
			temp := r.NewPoly()
			p1Coeffs := make([]*big.Int, r.N())
			p2Coeffs := make([]*big.Int, r.N())
			p3Coeffs := make([]*big.Int, r.N())

			r.PolyToBigint(*(*M)[i][j], 1, p1Coeffs)
			r.PolyToBigint(*vec[j], 1, p2Coeffs)

			for k := 0; k < r.N(); k++ {
				p3Coeffs[k] = new(big.Int).Mul(p1Coeffs[k], p2Coeffs[k])
				p3Coeffs[k] = new(big.Int).Mod(p3Coeffs[k], r.Modulus())
			}

			r.SetCoefficientsBigint(p3Coeffs, temp)
			r.Add(*result[i], temp, *result[i])
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertMatrixFromNTT(r, M)
	ConvertVectorFromNTT(r, vec)
}

// MatrixMatrixMul performs matrix-matrix multiplication.
// It takes two matrices of ring.Poly pointers, M1 of dimensions m x p and M2 of dimensions p x n,
// and outputs the result in a given result matrix of dimensions m x n.
func MatrixMatrixMul(r *ring.Ring, M1, M2 *[][]*ring.Poly, result *[][]*ring.Poly) {
	if M1 == nil || M2 == nil || len(*M1) == 0 || len(*M2) == 0 || len((*M1)[0]) != len(*M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(*M1)
	p := len((*M1)[0]) // Assuming all rows in M1 are of the same length
	n := len((*M2)[0]) // Assuming all rows in M2 are of the same length

	// Convert all elements of M1 and M2 to the NTT domain
	// Convert all elements of M1 and M2 to the NTT domain
	ConvertMatrixToNTT(r, M1)
	ConvertMatrixToNTT(r, M2)

	// Initialize the result matrix with zeros
	for i := 0; i < m; i++ {
		(*result)[i] = make([]*ring.Poly, n)
		for j := 0; j < n; j++ {
			newPoly := r.NewPoly()
			(*result)[i][j] = &newPoly
		}
	}

	// Perform matrix multiplication coefficient-wise
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				temp := r.NewPoly()
				p1Coeffs := make([]*big.Int, r.N())
				p2Coeffs := make([]*big.Int, r.N())
				p3Coeffs := make([]*big.Int, r.N())

				r.PolyToBigint(*(*M1)[i][k], 1, p1Coeffs)
				r.PolyToBigint(*(*M2)[k][j], 1, p2Coeffs)

				for l := 0; l < r.N(); l++ {
					p3Coeffs[l] = new(big.Int).Mul(p1Coeffs[l], p2Coeffs[l])
					p3Coeffs[l] = new(big.Int).Mod(p3Coeffs[l], r.Modulus())
				}

				r.SetCoefficientsBigint(p3Coeffs, temp)
				r.Add(*(*result)[i][j], temp, *(*result)[i][j])
			}
		}
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertMatrixFromNTT(r, result)
	ConvertMatrixFromNTT(r, M1)
	ConvertMatrixFromNTT(r, M2)
}

// MatrixAdd adds two matrices of ring.Poly element-wise and stores the result in a given result matrix.
func MatrixAdd(r *ring.Ring, M1, M2, result *[][]*ring.Poly) {
	if M1 == nil || M2 == nil || len(*M1) == 0 || len(*M2) == 0 || len(*M1) != len(*M2) || len((*M1)[0]) != len((*M2)[0]) {
		log.Fatalf("Matrix dimensions must match for element-wise addition.")
		return
	}

	m := len(*M1)      // Number of rows
	n := len((*M1)[0]) // Number of columns

	// Ensure the result matrix is initialized
	for i := 0; i < m; i++ {
		if (*result)[i] == nil {
			(*result)[i] = make([]*ring.Poly, n)
		}
		for j := 0; j < n; j++ {
			if (*result)[i][j] == nil {
				newPoly := r.NewPoly()
				(*result)[i][j] = &newPoly
			}
			r.Add(*(*M1)[i][j], *(*M2)[i][j], *(*result)[i][j])
		}
	}
}

// VectorPolyMul performs element-wise multiplication of a vector by a polynomial.
// It takes a vector of ring.Poly pointers, a single ring.Poly pointer, and outputs the result in a given result vector.
func VectorPolyMul(r *ring.Ring, vec []*ring.Poly, poly *ring.Poly, result []*ring.Poly) {
	// Convert the polynomial to the NTT domain
	r.NTT(*poly, *poly)

	for i := range vec {
		if result[i] == nil {
			result[i] = r.NewPoly().CopyNew()
		}
	}

	// Convert all elements of the vector to the NTT domain
	ConvertVectorToNTT(r, vec)
	ConvertVectorToNTT(r, result)

	// Perform the multiplications coefficient-wise
	for i := range vec {
		p1Coeffs := make([]*big.Int, r.N())
		p2Coeffs := make([]*big.Int, r.N())
		p3Coeffs := make([]*big.Int, r.N())

		r.PolyToBigint(*vec[i], 1, p1Coeffs)
		r.PolyToBigint(*poly, 1, p2Coeffs)

		for j := 0; j < r.N(); j++ {
			p3Coeffs[j] = new(big.Int).Mul(p1Coeffs[j], p2Coeffs[j])
			p3Coeffs[j] = new(big.Int).Mod(p3Coeffs[j], r.Modulus())
		}

		r.SetCoefficientsBigint(p3Coeffs, *result[i])
	}

	// Convert the result and all other polynomials back to the original domain
	ConvertVectorFromNTT(r, result)
	ConvertVectorFromNTT(r, vec)
	r.INTT(*poly, *poly)
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

// NAIVE MULTIPLICATIONS

// MatrixVectorMul performs matrix-vector multiplication.
// It takes a matrix of ring.Poly pointers, a vector of ring.Poly pointers, and outputs the result in a given result vector.
func MatrixVectorMulNaive(r *ring.Ring, M *[][]*ring.Poly, vec []*ring.Poly, result []*ring.Poly) {
	for i := range *M {
		result[i] = r.NewPoly().CopyNew()
		for j := range (*M)[i] {
			temp := r.NewPoly()
			MulPolyNaive(r, (*M)[i][j], vec[j], &temp)
			r.Add(*result[i], temp, *result[i]) // Accumulate the result
		}
	}
}

// MatrixMatrixMul performs matrix-matrix multiplication.
// It takes two matrices of ring.Poly pointers, M1 of dimensions m x p and M2 of dimensions p x n,
// and outputs the result in a given result matrix of dimensions m x n.
func MatrixMatrixMulNaive(r *ring.Ring, M1, M2 *[][]*ring.Poly, result *[][]*ring.Poly) {
	if M1 == nil || M2 == nil || len(*M1) == 0 || len(*M2) == 0 || len((*M1)[0]) != len(*M2) {
		log.Fatalf("Matrix dimensions are not compatible for multiplication.")
		return
	}

	m := len(*M1)
	p := len((*M1)[0]) // Assuming all rows in M1 are of the same length
	n := len((*M2)[0]) // Assuming all rows in M2 are of the same length

	// Initialize the result matrix with zeros
	for i := 0; i < m; i++ {
		(*result)[i] = make([]*ring.Poly, n)
		for j := 0; j < n; j++ {
			newPoly := r.NewPoly()
			(*result)[i][j] = &newPoly
		}
	}

	// Perform matrix multiplication
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < p; k++ {
				temp := r.NewPoly()
				MulPolyNaive(r, (*M1)[i][k], (*M2)[k][j], &temp)
				r.Add(*(*result)[i][j], temp, *(*result)[i][j])
			}
		}
	}
}

// VectorPolyMulNaive performs element-wise multiplication of a vector by a polynomial.
// It takes a vector of ring.Poly pointers, a single ring.Poly pointer, and outputs the result in a given result vector.
func VectorPolyMulNaive(r *ring.Ring, vec []*ring.Poly, poly *ring.Poly, result []*ring.Poly) {
	for i := range vec {
		if result[i] == nil {
			newPoly := r.NewPoly()
			result[i] = &newPoly
		}
		MulPolyNaive(r, vec[i], poly, result[i]) // Multiply each vector element by the polynomial
	}
}

// SAMPLER FUNCTIONS

// SamplePolyVector samples a vector of polynomials of a given length using the provided sampler.
func SamplePolyVector(length int, sampler ring.Sampler) []*ring.Poly {
	vector := make([]*ring.Poly, length)
	for i := 0; i < length; i++ {
		element := sampler.ReadNew()
		vector[i] = &element
	}
	return vector
}

// SamplePolyMatrix samples a matrix of polynomials with given dimensions (rows and cols) using the provided sampler.
func SamplePolyMatrix(rows, cols int, sampler ring.Sampler) [][]*ring.Poly {
	matrix := make([][]*ring.Poly, rows)
	for i := 0; i < rows; i++ {
		matrix[i] = make([]*ring.Poly, cols)
		for j := 0; j < cols; j++ {
			element := sampler.ReadNew()
			matrix[i][j] = &element
		}
	}
	return matrix
}

// PRINT FUNCTIONS

func PrintMatrix(label string, matrix *[][]*ring.Poly) {
	log.Println(label)
	for i, row := range *matrix {
		for j, poly := range row {
			log.Printf("[%d][%d]: %v\n", i, j, poly.Coeffs[0])
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

// NTT CONVERSION

// ConvertMatrixToNTT converts a matrix of polynomials to the NTT domain.
func ConvertMatrixToNTT(r *ring.Ring, M *[][]*ring.Poly) {
	for i := range *M {
		for j := range (*M)[i] {
			r.NTT(*(*M)[i][j], *(*M)[i][j])
		}
	}
}

// ConvertMatrixFromNTT converts a matrix of polynomials from the NTT domain back to the standard domain.
func ConvertMatrixFromNTT(r *ring.Ring, M *[][]*ring.Poly) {
	for i := range *M {
		for j := range (*M)[i] {
			r.INTT(*(*M)[i][j], *(*M)[i][j])
		}
	}
}

// ConvertVectorToNTT converts a vector of polynomials to the NTT domain.
func ConvertVectorToNTT(r *ring.Ring, vec []*ring.Poly) {
	for i := range vec {
		r.NTT(*vec[i], *vec[i])
	}
}

// ConvertVectorFromNTT converts a vector of polynomials from the NTT domain back to the standard domain.
func ConvertVectorFromNTT(r *ring.Ring, vec []*ring.Poly) {
	for i := range vec {
		r.INTT(*vec[i], *vec[i])
	}
}
