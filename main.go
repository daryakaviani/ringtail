// lattice/main.go

package main

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"fmt"
	"lattice-threshold-signature/utils"
	"log"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"golang.org/x/crypto/sha3"
)

// PARAMETERS
const (
	m = 11
	n = 9
	d = 10 // Length of joint noise vector

	p         = 2 ^ 30
	t         = 2 // Active threshold
	k         = 3 // Total number of parties
	ell       = 1
	beta      = 10
	betaDelta = 100000000000
	kappa     = 23
	logN      = 8
	bound_e   = 1000
	sigma_e   = 32
	boundE    = 1000
	sigmaE    = 32
	boundStar = 24296004000 * 2
	sigmaStar = 24296004000
	boundU    = 0
	sigmaU    = 0
	keySize   = 30
)

// Party struct holds all state and methods for a party in the protocol
type Party struct {
	ID             int
	Ring           *ring.Ring
	UniformSampler *ring.UniformSampler
	SkShare        []*ring.Poly
	Seed           [][]byte
	Mask           []*ring.Poly
	R              *[][]*ring.Poly
	C              *ring.Poly
	H              []*ring.Poly
	Lambda         *ring.Poly
}

var q = ring.Qi60[0]

// NewParty initializes a new Party instance
func NewParty(id int, r *ring.Ring, sampler *ring.UniformSampler) *Party {
	return &Party{
		ID:             id,
		Ring:           r,
		UniformSampler: sampler,
	}
}

// main function orchestrates the threshold signature protocol
func main() {
	var setupDuration, genDuration, signRound1Duration, signRound2Duration, finalizeDuration, verifyDuration time.Duration

	randomKey := make([]byte, keySize)

	r, _ := ring.NewRing(1<<logN, []uint64{q})

	prng, _ := sampling.NewKeyedPRNG(randomKey)
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := randomKey

	// SETUP phase
	start := time.Now()
	A := Setup(uniformSampler)
	setupDuration = time.Since(start)

	parties := make([]*Party, k)
	for i := range parties {
		parties[i] = NewParty(i, r, uniformSampler)
	}

	// GEN: Generate secret shares and seeds
	start = time.Now()
	skShares, seeds, b := Gen(r, A, uniformSampler, trustedDealerKey)
	genDuration = time.Since(start)

	for partyID := 0; partyID < k; partyID++ {
		parties[partyID].SkShare = (*skShares)[partyID]
		parties[partyID].Seed = (*seeds)[partyID]
	}

	// SIGNATURE ROUND 1
	mu := "Hello, Threshold Signature!"
	sid := 1
	PRFKey := "PRF Key"
	T := []int{0, 1} // Active parties
	start = time.Now()
	lagrangeCoeffs := ComputeLagrangeCoefficients(r, T, big.NewInt(int64(q)))
	D := make(map[int]*[][]*ring.Poly)
	masks := make(map[int]*[]*ring.Poly)
	for _, partyID := range T {
		parties[partyID].Lambda = lagrangeCoeffs[partyID]
		parties[partyID].Seed = (*seeds)[partyID]
		D[partyID], masks[partyID] = parties[partyID].SignRound1(A, sid, []byte(PRFKey), T)
	}
	signRound1Duration = time.Since(start)

	// SIGNATURE ROUND 2
	start = time.Now()
	z := make(map[int]*[]*ring.Poly)
	for _, partyID := range T {
		z[partyID] = parties[partyID].SignRound2(A, b, &D, &masks, sid, mu, T, []byte(PRFKey), seeds)
	}
	signRound2Duration = time.Since(start)

	// SIGNATURE FINALIZE
	start = time.Now()
	finalParty := parties[0]
	_, sig, Delta := finalParty.SignFinalize(z, masks, A, b)
	finalizeDuration = time.Since(start)

	// Verify the signature
	start = time.Now()
	valid := Verify(r, sig, A, mu, b, finalParty.C, Delta, betaDelta)
	verifyDuration = time.Since(start)
	fmt.Printf("Signature Verification Result: %v\n", valid)

	// Print all durations
	fmt.Println("Setup duration:", setupDuration)
	fmt.Println("Gen duration:", genDuration)
	fmt.Println("Signature Round 1 duration:", signRound1Duration)
	fmt.Println("Signature Round 2 duration:", signRound2Duration)
	fmt.Println("Finalize duration:", finalizeDuration)
	fmt.Println("Verify duration:", verifyDuration)
}

// Setup generates the public parameters
func Setup(uniformSampler *ring.UniformSampler) *[][]*ring.Poly {
	A := utils.SamplePolyMatrix(m, n, uniformSampler)
	return &A
}

// Gen generates the secret shares, seeds, and the public parameter b
func Gen(r *ring.Ring, A *[][]*ring.Poly, uniformSampler *ring.UniformSampler, trustedDealerKey []byte) (*map[int][]*ring.Poly, *map[int][][]byte, []*ring.Poly) {
	prng, _ := sampling.NewKeyedPRNG(trustedDealerKey)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigma_e, Bound: bound_e}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	s := utils.SamplePolyVector(n, uniformSampler)

	skShares := ShamirSecretSharing(r, s, t, k)

	e := utils.SamplePolyVector(m, gaussianSampler)

	b := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, s, b)
	utils.VectorAdd(r, b, e, b)

	seeds := make(map[int][][]byte)
	for i := 0; i < k; i++ {
		seeds[i] = make([][]byte, k)
		for j := 0; j < k; j++ {
			seeds[i][j] = generateRandomSeed()
		}
	}

	return &skShares, &seeds, b
}

// generateRandomSeed generates a random seed of length ell
func generateRandomSeed() []byte {
	sd := make([]byte, ell)
	_, err := rand.Read(sd)
	if err != nil {
		log.Fatalf("Error generating random seed %v\n", err)
	}
	return sd
}

// SignRound1 performs the first round of signing
func (party *Party) SignRound1(A *[][]*ring.Poly, sid int, PRFKey []byte, T []int) (*[][]*ring.Poly, *[]*ring.Poly) {
	r := party.Ring
	uniformSampler := party.UniformSampler
	seeds := party.Seed

	mask := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, seeds[j], PRFKey)
		utils.VectorAdd(r, mask, mask_j, mask)
	}

	r_star := utils.SamplePolyVector(n, uniformSampler)

	skHash := party.PRNGKey()
	prng, _ := sampling.NewKeyedPRNG(skHash)
	gaussianParams := ring.DiscreteGaussian{Sigma: sigmaStar, Bound: boundStar}
	gaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)
	e_star := utils.SamplePolyVector(m, gaussianSampler)

	R_i := utils.SamplePolyMatrix(n, d-1, uniformSampler)

	gaussianParams = ring.DiscreteGaussian{Sigma: sigmaE, Bound: boundE}
	gaussianSampler = ring.NewGaussianSampler(prng, r, gaussianParams, false)
	E_i := utils.SamplePolyMatrix(m, d-1, gaussianSampler)

	concatenatedR := make([][]*ring.Poly, n)
	concatenatedE := make([][]*ring.Poly, m)

	for i := 0; i < n; i++ {
		concatenatedR[i] = append([]*ring.Poly{r_star[i]}, R_i[i]...)
	}
	party.R = &concatenatedR

	for i := 0; i < m; i++ {
		concatenatedE[i] = append([]*ring.Poly{e_star[i]}, E_i[i]...)
	}

	D := make([][]*ring.Poly, m)
	utils.MatrixMatrixMul(r, A, &concatenatedR, &D)
	utils.MatrixAdd(r, &concatenatedE, &D, &D)

	return &D, &mask
}

// SignRound2 performs the second round of signing
func (party *Party) SignRound2(A *[][]*ring.Poly, b []*ring.Poly, D *map[int]*[][]*ring.Poly, masks *map[int]*[]*ring.Poly, sid int, mu string, T []int, PRFKey []byte, seeds *map[int][][]byte) *[]*ring.Poly {
	r := party.Ring
	partyID := party.ID
	concatR := party.R
	s_i := party.SkShare
	lambda := party.Lambda

	u := make(map[int][]*ring.Poly, d)
	onePoly := r.NewMonomialXi(0)

	for _, j := range T {
		oneSlice := []*ring.Poly{onePoly.CopyNew()}
		u[j] = oneSlice
		if d > 1 {
			u_j := H_u(r, A, b, sid, j, D, mu, masks)
			u[j] = append(oneSlice, u_j...)
		}
	}

	h := make([]*ring.Poly, m)

	for j, D_j := range *D {
		D_j_u_j := make([]*ring.Poly, m)
		utils.MatrixVectorMul(r, D_j, u[j], D_j_u_j)
		utils.VectorAdd(r, h, D_j_u_j, h)
	}

	for _, poly := range h {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}
	party.H = h

	c := H_c(r, A, b, h, mu)
	party.C = c

	mPrime := make([]*ring.Poly, n)
	for _, j := range T {
		mask_j := PRF(r, sid, (*seeds)[j][partyID], PRFKey)
		utils.VectorAdd(r, mPrime, mask_j, mPrime)
	}

	z_i := make([]*ring.Poly, n)
	utils.MatrixVectorMul(r, concatR, u[partyID], z_i)
	utils.VectorAdd(r, z_i, mPrime, z_i)
	s_c_lambda := make([]*ring.Poly, n)
	utils.VectorPolyMul(r, s_i, c, s_c_lambda)
	utils.VectorPolyMul(r, s_c_lambda, lambda, s_c_lambda)
	utils.VectorAdd(r, z_i, s_c_lambda, z_i)

	return &z_i
}

// SignFinalize finalizes the signature
func (party *Party) SignFinalize(z map[int]*[]*ring.Poly, masks map[int]*[]*ring.Poly, A *[][]*ring.Poly, b []*ring.Poly) (*ring.Poly, *[]*ring.Poly, *[]*ring.Poly) {
	r := party.Ring
	c := party.C
	h := party.H

	z_sum := make([]*ring.Poly, n)

	for j, z_j := range z {
		mask_j := masks[j]
		for i := 0; i < n; i++ {
			if z_sum[i] == nil {
				z_sum[i] = r.NewPoly().CopyNew()
			}
			r.Sub(*(*z_j)[i], *(*mask_j)[i], *(*z_j)[i])
			r.Add(*z_sum[i], *(*z_j)[i], *z_sum[i])
		}
	}

	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, z_sum, Az)

	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}

	Delta := make([]*ring.Poly, m)
	utils.VectorSub(r, h, Az_bc, Delta)

	return party.C, &z_sum, &Delta
}

// Verify verifies the correctness of the signature
func Verify(r *ring.Ring, z *[]*ring.Poly, A *[][]*ring.Poly, mu string, b []*ring.Poly, c *ring.Poly, Delta *[]*ring.Poly, betaDelta uint64) bool {
	Az := make([]*ring.Poly, m)
	utils.MatrixVectorMul(r, A, *z, Az)

	bc := make([]*ring.Poly, m)
	utils.VectorPolyMul(r, b, c, bc)

	Az_bc := make([]*ring.Poly, m)
	utils.VectorSub(r, Az, bc, Az_bc)

	for _, poly := range Az_bc {
		utils.RoundCoeffsToNearestMultiple(r, poly, p)
	}

	Az_bc_Delta := make([]*ring.Poly, m)
	utils.VectorAdd(r, Az_bc, *Delta, Az_bc_Delta)

	computedC := H_c(r, A, b, Az_bc_Delta, mu)
	if !r.Equal(*c, *computedC) {
		return false
	}

	return CheckInfinityNorm(r, Delta, betaDelta)
}

// CheckInfinityNorm checks if the infinity norm of the vector of polynomial Delta is less than or equal to betaDelta
func CheckInfinityNorm(r *ring.Ring, Delta *[]*ring.Poly, betaDelta uint64) bool {
	maxValue := big.NewInt(0)

	for _, poly := range *Delta {
		coeffsBigint := make([]*big.Int, r.N())
		for i := range coeffsBigint {
			coeffsBigint[i] = big.NewInt(0)
		}
		r.PolyToBigintCentered(*poly, 1, coeffsBigint)

		for _, coeff := range coeffsBigint {
			absCoeff := new(big.Int).Abs(coeff)
			if absCoeff.Cmp(maxValue) == 1 {
				maxValue.Set(absCoeff)
			}
		}
	}

	betaDeltaBig := new(big.Int).SetUint64(betaDelta)

	log.Print("DELTA NORM:", maxValue)

	return new(big.Int).Mod(maxValue, new(big.Int).SetUint64(q)).Cmp(betaDeltaBig) <= 0
}

// HASHES & PRF

// H_u hashes parameters to a Gaussian distribution
func H_u(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, sid int, j int, D *map[int]*[][]*ring.Poly, mu string, masks *map[int]*[]*ring.Poly) []*ring.Poly {
	hasher := sha3.NewShake128()

	var buffer bytes.Buffer

	for _, row := range *A {
		for _, poly := range row {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	for _, poly := range b {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	binary.Write(&buffer, binary.BigEndian, int64(sid))

	binary.Write(&buffer, binary.BigEndian, int64(j))

	for _, D_h := range *D {
		for _, row := range *D_h {
			for _, poly := range row {
				data, err := poly.MarshalBinary()
				if err != nil {
					log.Fatalf("Error marshalling poly: %v\n", err)
				}
				buffer.Write(data)
			}
		}
	}

	buffer.WriteString(mu)

	for _, mask := range *masks {
		for _, poly := range *mask {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	_, err := hasher.Write(buffer.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	gaussianParams := ring.DiscreteGaussian{}
	hashGaussianSampler := ring.NewGaussianSampler(prng, r, gaussianParams, false)

	u_j := utils.SamplePolyVector(d-1, hashGaussianSampler)

	return u_j
}

// PRF generates pseudorandom ring elements
func PRF(r *ring.Ring, sid int, sd_ij []byte, PRFKey []byte) []*ring.Poly {
	hasher := sha3.NewShake128()

	var buffer bytes.Buffer

	binary.Write(&buffer, binary.BigEndian, PRFKey)
	binary.Write(&buffer, binary.BigEndian, sd_ij)
	binary.Write(&buffer, binary.BigEndian, int64(sid))

	_, err := hasher.Write(buffer.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	PRFUniformSampler := ring.NewUniformSampler(prng, r)
	mask := utils.SamplePolyVector(n, PRFUniformSampler)
	return mask
}

// H_c hashes to low norm ring elements
func H_c(r *ring.Ring, A *[][]*ring.Poly, b []*ring.Poly, h []*ring.Poly, mu string) *ring.Poly {
	hasher := sha3.NewShake128()

	var buffer bytes.Buffer

	for _, row := range *A {
		for _, poly := range row {
			data, err := poly.MarshalBinary()
			if err != nil {
				log.Fatalf("Error marshalling poly: %v\n", err)
			}
			buffer.Write(data)
		}
	}

	for _, poly := range b {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	for _, poly := range h {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	binary.Write(&buffer, binary.BigEndian, []byte(mu))

	_, err := hasher.Write(buffer.Bytes())
	if err != nil {
		log.Fatalf("Error writing hash: %v\n", err)
	}

	hashOutputLength := keySize
	hashOutput := make([]byte, hashOutputLength)
	_, err = hasher.Read(hashOutput)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}

	prng, _ := sampling.NewKeyedPRNG(hashOutput)
	ternaryParams := ring.Ternary{H: kappa}
	ternarySampler, err := ring.NewTernarySampler(prng, r, ternaryParams, false)
	if err != nil {
		log.Fatalf("Error creating ternary sampler: %v", err)
	}
	c := ternarySampler.ReadNew()

	return &c
}

// ShamirSecretSharing shares each coefficient of a vector of ring.Poly across k parties using (t, k)-threshold Shamir secret sharing.
func ShamirSecretSharing(r *ring.Ring, s []*ring.Poly, t, k int) map[int][]*ring.Poly {
	degree := r.N()

	shares := make(map[int][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		shares[i] = make([]*ring.Poly, len(s))
		for j := range shares[i] {
			newPoly := r.NewPoly()
			shares[i][j] = &newPoly
		}
	}

	for polyIndex, poly := range s {
		coeffs := make([]*big.Int, degree)
		r.PolyToBigint(*poly, 1, coeffs)

		for coeffIndex, secret := range coeffs {
			polyCoeffs := make([]*big.Int, t)
			polyCoeffs[0] = secret
			for i := 1; i < t; i++ {
				randomCoeff := big.NewInt(20)
				polyCoeffs[i] = randomCoeff
			}

			for i := 1; i <= k; i++ {
				x := big.NewInt(int64(i))
				shareValue := big.NewInt(0)
				xPow := big.NewInt(1)

				for _, coeff := range polyCoeffs {
					term := new(big.Int).Mul(coeff, xPow)
					shareValue.Add(shareValue, term)
					shareValue.Mod(shareValue, big.NewInt(int64(q)))
					xPow.Mul(xPow, x)
				}

				if shares[i-1][polyIndex].Coeffs[0] == nil {
					shares[i-1][polyIndex].Coeffs[0] = make([]uint64, degree)
				}

				shares[i-1][polyIndex].Coeffs[0][coeffIndex] = shareValue.Uint64()
			}
		}
	}

	return shares
}

// ComputeLagrangeCoefficients computes the Lagrange coefficients for interpolation based on the indices of available shares.
func ComputeLagrangeCoefficients(r *ring.Ring, T []int, modulus *big.Int) map[int]*ring.Poly {
	lagrangeCoefficients := make(map[int]*ring.Poly)

	for i := 0; i < len(T); i++ {
		xi := big.NewInt(int64(T[i] + 1))
		numerator := big.NewInt(1)
		denominator := big.NewInt(1)

		for j := 0; j < len(T); j++ {
			if i != j {
				xj := big.NewInt(int64(T[j] + 1))
				numerator.Mul(numerator, new(big.Int).Neg(xj))
				numerator.Mod(numerator, modulus)
				temp := new(big.Int).Sub(xi, xj)
				denominator.Mul(denominator, temp)
				denominator.Mod(denominator, modulus)
			}
		}

		denomInv := new(big.Int).ModInverse(denominator, modulus)
		coeff := new(big.Int).Mul(numerator, denomInv)
		coeff.Mod(coeff, modulus)

		lagrangePoly := r.NewPoly()
		r.SetCoefficientsBigint([]*big.Int{coeff}, lagrangePoly)
		lagrangeCoefficients[T[i]] = &lagrangePoly
	}

	return lagrangeCoefficients
}

// PRNGKey generates a key for PRNG using the secret key share
func (party *Party) PRNGKey() []byte {
	skShare := party.SkShare

	hasher := sha3.NewShake128()
	var buffer bytes.Buffer
	for _, poly := range skShare {
		data, err := poly.MarshalBinary()
		if err != nil {
			log.Fatalf("Error marshalling poly: %v\n", err)
		}
		buffer.Write(data)
	}

	hasher.Write(buffer.Bytes())

	hashOutputLength := keySize
	skHash := make([]byte, hashOutputLength)
	_, err := hasher.Read(skHash)
	if err != nil {
		log.Fatalf("Error reading hash: %v\n", err)
	}
	return skHash
}
