package main

import (
	"fmt"
	"math/big"
)

// polyMultInCyclotomicRing multiplies two polynomials within the cyclotomic ring x^8 + 1 and returns the result.
// Each polynomial is represented as a slice of big.Int pointers, sorted from least to most significant.
func polyMultInCyclotomicRing(p1, p2 []*big.Int, q *big.Int) []*big.Int {
	degree := 8 // Since we are in a ring modulo x^8 + 1
	// Initialize result slice with big.Ints set to zero
	result := make([]*big.Int, degree)
	for i := range result {
		result[i] = big.NewInt(0)
	}

	// Polynomial multiplication (convolution)
	for i := range p1 {
		for j := range p2 {
			if i+j < degree {
				// Multiply coefficients and add to the right place
				temp := new(big.Int).Mul(p1[i], p2[j])
				result[i+j].Add(result[i+j], temp)
			} else {
				if (((i+j)-((i+j)%degree))/8)%2 == 0 {
					// Wrap around due to cyclotomic ring, i+j >= degree
					temp := new(big.Int).Mul(p1[i], p2[j])
					result[(i+j)%degree].Add(result[(i+j)%degree], temp)
				} else {
					// Wrap around due to cyclotomic ring, i+j >= degree
					temp := new(big.Int).Mul(p1[i], p2[j])
					result[(i+j)%degree].Sub(result[(i+j)%degree], temp) // subtracting because of x^8 + 1
				}
			}
		}
	}

	// Reduce each coefficient modulo q
	for i := range result {
		result[i].Mod(result[i], q)
	}

	return result
}

func main() {
	// Example usage:
	q := big.NewInt(61) // Modulus

	// Example polynomials
	p1 := []*big.Int{
		big.NewInt(34),
		big.NewInt(6),
		big.NewInt(30),
		big.NewInt(16),
		big.NewInt(9),
		big.NewInt(60),
		big.NewInt(5),
		big.NewInt(15),
	}
	p2 := []*big.Int{
		big.NewInt(28),
		big.NewInt(42),
		big.NewInt(26),
		big.NewInt(37),
		big.NewInt(21),
		big.NewInt(21),
		big.NewInt(10),
		big.NewInt(0),
	}

	result := polyMultInCyclotomicRing(p1, p2, q)

	fmt.Print(result)

	// correct := []*big.Int{
	// 	big.NewInt(1394366894685604571),
	// 	big.NewInt(1315606860825056612),
	// 	big.NewInt(204594538006257517),
	// 	big.NewInt(1703265942265006841),
	// 	big.NewInt(1069366183506519705),
	// 	big.NewInt(69891534085996005),
	// 	big.NewInt(1910119678996225288),
	// 	big.NewInt(1460472761970572730),
	// }

	// // Check and print results
	// allCorrect := true
	// for i, coef := range result {
	// 	if coef.Cmp(correct[i]) != 0 {
	// 		allCorrect = false
	// 		fmt.Printf("Mismatch at y^%d: got %s, expected %s\n", i, coef.String(), correct[i].String())
	// 	} else {
	// 		fmt.Printf("Coefficient of y^%d is correct: %s\n", i, coef.String())
	// 	}
	// }

	// if allCorrect {
	// 	fmt.Println("All coefficients are correct.")
	// } else {
	// 	fmt.Println("There are incorrect coefficients.")
	// }
}
