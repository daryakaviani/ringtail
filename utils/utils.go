package utils

import (
	"log"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
)

var q = ring.Qi60[0]

const logN = 3

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

// PRINT FUNCTIONS

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
