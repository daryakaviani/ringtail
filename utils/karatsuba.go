package utils

import (
	"fmt"
	"math/big"
	"math/rand"
	"time"
)

func karatsubaSquare(f []*big.Int) []*big.Int {
	n := len(f)
	if n == 1 {
		res := big.NewInt(0).Mul(f[0], f[0])
		return []*big.Int{res, big.NewInt(0)}
	}

	fe := make([]*big.Int, n/2)
	fo := make([]*big.Int, n/2)
	fs := make([]*big.Int, n/2)

	for i := 0; i < n/2; i++ {
		fe[i] = f[2*i]
		fo[i] = f[2*i+1]
		fs[i] = new(big.Int).Add(fe[i], fo[i])
	}

	fe = karatsubaSquare(fe)
	fo = karatsubaSquare(fo)
	fs = karatsubaSquare(fs)

	g := make([]*big.Int, 2*n)
	for i := 0; i < n; i++ {
		g[2*i] = new(big.Int).Add(fe[i], func() *big.Int {
			if i > 0 {
				return fo[i-1]
			}
			return big.NewInt(0)
		}())
		g[2*i+1] = new(big.Int).Sub(new(big.Int).Sub(fs[i], fe[i]), fo[i])
	}

	return g
}

func redKaratsubaSquare(f []*big.Int) []*big.Int {
	n := len(f)
	g := karatsubaSquare(f)
	res := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		res[i] = new(big.Int).Sub(g[i], g[i+n])
	}
	return res
}

func KaratsubaNorm(f []*big.Int) *big.Int {
	n := len(f)
	if n == 1 {
		return f[0]
	}

	fe := make([]*big.Int, n/2)
	fo := make([]*big.Int, n/2)

	for i := 0; i < n/2; i++ {
		fe[i] = f[2*i]
		fo[i] = f[2*i+1]
	}

	m := n / 2
	fe = redKaratsubaSquare(fe)
	fo = redKaratsubaSquare(fo)

	foShifted := make([]*big.Int, m)
	foShifted[0] = new(big.Int).Neg(fo[m-1])
	for i := 1; i < m; i++ {
		foShifted[i] = fo[i-1]
	}

	g := make([]*big.Int, m)
	for i := 0; i < m; i++ {
		g[i] = new(big.Int).Sub(fe[i], foShifted[i])
	}

	return KaratsubaNorm(g)
}

// For testing purposes
func TestKaratsubaNorm(n int, q *big.Int) ([]*big.Int, *big.Int) {
	f := make([]*big.Int, n)
	for i := range f {
		halfQ := new(big.Int).Div(q, big.NewInt(2))
		f[i] = new(big.Int).Sub(new(big.Int).Rand(rand.New(rand.NewSource(time.Now().UnixNano())), q), halfQ)
	}

	start := time.Now()
	r := KaratsubaNorm(f)
	elapsed := time.Since(start)

	fmt.Printf("Computed norm of a random element of R_q in %.2f ms.\n", float64(elapsed.Milliseconds()))
	fmt.Printf("The result is a %d-bit integer.\n", r.BitLen())

	return f, r
}
