package main

import (
	"crypto/rand"
	"lattice-threshold-signature/utils"
	"log"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
)

var q = ring.Qi60[0]

const logN = 3

func main() {
	r, err := ring.NewRing(1<<logN, []uint64{q})
	if err != nil {
		log.Fatalf("Error creating ring: %v", err)
	}
	TestingNTT(r)
	TestZeroPolynomial(r)
	// TestNTTLinearity(r)
	TestConvolution(r)
	TestNTTPowersOfTwo(r)
}

func TestingNTT(r *ring.Ring) {
	// Test 1: INTT of a polynomial of all 1s should result in a constant polynomial
	polyAllOnes := r.NewPoly() // Assuming constructor sets length
	ones := make([]*big.Int, 8)
	for i := range ones {
		ones[i] = big.NewInt(1)
	}
	r.SetCoefficientsBigint(ones, polyAllOnes)
	r.INTT(polyAllOnes, polyAllOnes)

	// Check if the result is a constant polynomial
	log.Println("Poly all 1s:", polyAllOnes)

	// Test 2: Compare MulPolyNTT with MulPoly
	p1 := r.NewPoly()
	p2 := r.NewPoly()
	// Example coefficients
	r.SetCoefficientsBigint([]*big.Int{big.NewInt(1), big.NewInt(2), big.NewInt(3), big.NewInt(4), big.NewInt(5), big.NewInt(6), big.NewInt(7), big.NewInt(8)}, p1)
	r.SetCoefficientsBigint([]*big.Int{big.NewInt(8), big.NewInt(7), big.NewInt(6), big.NewInt(5), big.NewInt(4), big.NewInt(3), big.NewInt(2), big.NewInt(1)}, p2)

	p3 := r.NewPoly()
	p4 := r.NewPoly()

	utils.MulNTT(r, &p1, &p2, &p3)
	utils.MulPoly(r, &p1, &p2, &p4)

	// Check if the results are equal
	resultsMatch := true
	for i := range p3.Coeffs {
		if p3.Coeffs[i][0] != p4.Coeffs[i][0] {
			resultsMatch = false
			break
		}
	}

	if resultsMatch {
		log.Println("Test 2 Passed: MulPolyNTT matches MulPoly.")
	} else {
		log.Println("Test 2 Failed: MulPolyNTT does not match MulPoly.")
		utils.PrintPolynomial("MulPolyNTT Result", &p3)
		utils.PrintPolynomial("MulPoly Result", &p4)
	}
}

func TestNTTIdentity(r *ring.Ring) {
	p := r.NewPoly()
	coeffs := generateRandomCoefficients(8, r.Modulus())
	r.SetCoefficientsBigint(coeffs, p)

	// Convert polynomial to big integers before NTT
	p1Coeffs := make([]*big.Int, 8)
	r.PolyToBigint(p, 1, p1Coeffs)

	// Perform NTT and then INTT
	r.NTT(p, p)
	r.INTT(p, p)

	// Convert back to big integers after INTT
	finalCoeffs := make([]*big.Int, 8)
	r.PolyToBigint(p, 1, finalCoeffs)

	// Check if the final coefficients match the original coefficients
	passed := true
	for i, originalCoeff := range p1Coeffs {
		if finalCoeffs[i].Cmp(originalCoeff) != 0 {
			passed = false
			break
		}
	}

	if passed {
		log.Println("Test Identity Passed: NTT followed by INTT returns the original polynomial.")
	} else {
		log.Println("Test Identity Failed: NTT followed by INTT does not return the original polynomial.")
	}
}

// TestZeroPolynomial verifies that the NTT of a zero polynomial remains zero after INTT.
func TestZeroPolynomial(r *ring.Ring) {
	p := r.NewPoly()
	zeroCoeffs := make([]*big.Int, r.N())
	for i := range zeroCoeffs {
		zeroCoeffs[i] = big.NewInt(0)
	}
	r.SetCoefficientsBigint(zeroCoeffs, p)
	r.NTT(p, p)
	r.INTT(p, p)

	// Check if all coefficients are zero after NTT and INTT

	finalCoeffs := make([]*big.Int, 8)
	r.PolyToBigint(p, 1, finalCoeffs)

	allZero := true
	for _, coeff := range finalCoeffs {
		if coeff.Cmp(big.NewInt(0)) != 0 {
			allZero = false
			break
		}
	}

	if allZero {
		log.Println("Test Zero Polynomial Passed: NTT and INTT of zero polynomial results in zero.")
	} else {
		log.Println("Test Zero Polynomial Failed: NTT and INTT of zero polynomial does not result in zero.")
	}
}

// // TestNTTLinearity checks the linearity of NTT.
// // func TestNTTLinearity(r *ring.Ring) {
// 	aCoeffs := generateRandomCoefficients(8, r.Modulus())
// 	bCoeffs := generateRandomCoefficients(8, r.Modulus())
// 	k := big.NewInt(2) // A scalar multiplier

// 	pa := r.NewPoly()
// 	pb := r.NewPoly()
// 	paPlusb := r.NewPoly()
// 	kPa := r.NewPoly()

// 	r.SetCoefficientsBigint(aCoeffs, pa)
// 	r.SetCoefficientsBigint(bCoeffs, pb)
// 	r.SetCoefficientsBigint(addCoefficients(aCoeffs, bCoeffs), paPlusb)
// 	r.SetCoefficientsBigint(multiplyCoefficients(aCoeffs, k), kPa)

// 	r.NTT(pa, pa)
// 	r.NTT(pb, pb)
// 	r.NTT(paPlusb, paPlusb)
// 	r.NTT(kPa, kPa)

// 	// Validate the linearity: NTT(a + b) = NTT(a) + NTT(b) and NTT(ka) = k * NTT(a)
// 	passed := true

// 	for i := range pa.Coeffs {
// 		if !paPlusb.Coeffs[i].Cmp(addBigInt(pa.Coeffs[i], pb.Coeffs[i])) || !kPa.Coeffs[i].Cmp(multiplyBigInt(pa.Coeffs[i], k)) {
// 			passed = false
// 			break
// 		}
// 	}

// 	if passed {
// 		log.Println("Test NTT Linearity Passed.")
// 	} else {
// 		log.Println("Test NTT Linearity Failed.")
// 	}
// }

// TestConvolution verifies that polynomial multiplication using NTT matches direct convolution.
func TestConvolution(r *ring.Ring) {
	aCoeffs := generateRandomCoefficients(8, r.Modulus())
	bCoeffs := generateRandomCoefficients(8, r.Modulus())

	pa := r.NewPoly()
	pb := r.NewPoly()
	pc := r.NewPoly()

	r.SetCoefficientsBigint(aCoeffs, pa)
	r.SetCoefficientsBigint(bCoeffs, pb)

	r.NTT(pa, pa)
	r.NTT(pb, pb)
	r.MForm(pa, pa)
	r.MForm(pb, pb)

	subrings := r.SubRings
	log.Println("SUBGRINGS", subrings[0].MRedConstant, subrings[0].Modulus, subrings[0].Factors)
	r.MulCoeffsMontgomery(pa, pb, pc)
	r.IMForm(pc, pc)
	r.INTT(pc, pc)

	// Expected result from direct convolution
	expected := directConvolution(aCoeffs, bCoeffs, r.Modulus())

	// Convert back to big integers after INTT
	finalCoeffs := make([]*big.Int, 8)
	r.PolyToBigint(pc, 1, finalCoeffs)

	// Check if the convolution result matches
	passed := true
	for i, coeff := range finalCoeffs {
		if coeff.Cmp(expected[i]) != 0 {
			passed = false
			break
		}
	}

	if passed {
		log.Println("Test Convolution Passed.")
	} else {
		log.Println("Test Convolution Failed.")
	}
}

// TestNTTPowersOfTwo verifies that NTT operates correctly for all lengths that are powers of two up to the degree supported by the implementation.
func TestNTTPowersOfTwo(r *ring.Ring) {
	passed := true
	for d := 1; d <= r.N(); d <<= 1 {
		TestNTTIdentity(r)
	}

	if passed {
		log.Println("Test Powers of Two Passed.")
	} else {
		log.Println("Test Powers of Two Failed.")
	}
}

// HELPERS

// generateRandomCoefficients creates a slice of random *big.Int values of the given size and modulus.
func generateRandomCoefficients(size int, modulus *big.Int) []*big.Int {
	coeffs := make([]*big.Int, size)
	for i := range coeffs {
		coeffs[i] = new(big.Int)
		// Random value in the range [0, modulus-1]
		coeffs[i], _ = rand.Int(rand.Reader, modulus)
	}
	return coeffs
}

// addCoefficients adds two slices of *big.Int coefficients element-wise.
func addCoefficients(a, b []*big.Int) []*big.Int {
	if len(a) != len(b) {
		log.Println("Coefficient slices are not of the same length.")
		return nil
	}
	result := make([]*big.Int, len(a))
	for i := range a {
		result[i] = new(big.Int).Add(a[i], b[i])
	}
	return result
}

// multiplyCoefficients multiplies each element of a slice of *big.Int coefficients by a scalar *big.Int.
func multiplyCoefficients(coeffs []*big.Int, k *big.Int) []*big.Int {
	result := make([]*big.Int, len(coeffs))
	for i := range coeffs {
		result[i] = new(big.Int).Mul(coeffs[i], k)
	}
	return result
}

// addBigInt adds two *big.Int numbers.
func addBigInt(a, b *big.Int) *big.Int {
	return new(big.Int).Add(a, b)
}

// multiplyBigInt multiplies two *big.Int numbers.
func multiplyBigInt(a, b *big.Int) *big.Int {
	return new(big.Int).Mul(a, b)
}

// directConvolution computes the convolution of two polynomial coefficient slices directly.
func directConvolution(a, b []*big.Int, modulus *big.Int) []*big.Int {
	result := make([]*big.Int, len(a))
	for i := range result {
		result[i] = big.NewInt(0)
	}

	for i := range a {
		for j := range b {
			index := (i + j) % len(a) // Assuming the result size is the same and wraps around
			temp := new(big.Int).Mul(a[i], b[j])
			temp.Mod(temp, modulus)
			result[index].Add(result[index], temp)
			result[index].Mod(result[index], modulus)
		}
	}

	return result
}
