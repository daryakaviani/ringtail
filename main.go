package main

import (
	"bufio"
	"fmt"
	"lattice-threshold-signature/networking"
	"lattice-threshold-signature/primitives"
	"lattice-threshold-signature/sign"
	"log"
	"math/big"
	"net"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run main.go partyID")
		os.Exit(1)
	}

	partyIDString := os.Args[1]

	iters, err := strconv.Atoi(os.Args[2])
	if err != nil {
		fmt.Println("Error: Please enter a valid integer.")
		os.Exit(1)
	}

	parties, err := strconv.Atoi(os.Args[3])
	if err != nil {
		fmt.Println("Error: Please enter a valid integer.")
		os.Exit(1)
	}
	sign.K = parties
	sign.Threshold = parties

	if partyIDString == "l" {
		sign.LocalRun(iters)
		return
	}

	partyID, err := strconv.Atoi(partyIDString)
	if err != nil {
		fmt.Println("Error: Please enter a valid integer.")
		os.Exit(1)
	}

	// Initialize P2P communication
	comm := &networking.P2PComm{
		Socks: make(map[int]*net.Conn),
		Rank:  partyID,
	}

	// Establish connections
	var connWg sync.WaitGroup
	connWg.Add(1)
	go networking.EstablishConnections(&connWg, comm, partyID, sign.K)
	connWg.Wait()

	var setupDuration, genDuration, signRound1Duration, signRound2PreprocessDuration, signRound2Duration, finalizeDuration, verifyDuration time.Duration
	var genStart, genEnd, signRound1Start, signRound1End, signRound2Start, signRound2End, combinerReceiveEnd, combinerFinalizeEnd time.Time
	var A structs.Matrix[ring.Poly]
	var b structs.Vector[ring.Poly]
	randomKey := make([]byte, sign.KeySize)
	r, _ := ring.NewRing(1<<sign.LogN, []uint64{sign.Q})
	r_xi, _ := ring.NewRing(1<<sign.LogN, []uint64{sign.QXi})
	r_nu, _ := ring.NewRing(1<<sign.LogN, []uint64{sign.QNu})
	prng, _ := sampling.NewKeyedPRNG(randomKey)
	uniformSampler := ring.NewUniformSampler(prng, r)
	trustedDealerKey := randomKey
	// Create your own party
	party := sign.NewParty(partyID, r, r_xi, r_nu, uniformSampler)

	// Variables to be used in SIGNATURE ROUNDS
	D := make(map[int]structs.Matrix[ring.Poly])
	MACs := make(map[int]map[int][]byte)
	mu := "Message"
	T := make([]int, sign.K)
	for i := 0; i < sign.K; i++ {
		T[i] = i
	}
	lagrangeCoeffs := primitives.ComputeLagrangeCoefficients(r, T, big.NewInt(int64(sign.Q)))

	// Trusted dealer
	if partyID == sign.TrustedDealerID {
		// GEN: Generate secret shares, seeds, and MAC keys
		start := time.Now()
		AMat, skShares, seeds, MACKeys, bVec := sign.Gen(r, r_xi, uniformSampler, trustedDealerKey, lagrangeCoeffs)

		b = bVec
		A = AMat
		genDuration = time.Since(start)
		genStart = time.Now()

		// Give party 0 its own values
		party.SkShare = skShares[partyID]
		party.Seed = seeds
		party.MACKeys = MACKeys[partyID]

		// Send out public information & trusted dealer data
		var sendWg sync.WaitGroup
		for i := 0; i < sign.K; i++ {
			if i != sign.TrustedDealerID {
				sendWg.Add(1)
				go func(i int) {
					defer sendWg.Done()
					writer := bufio.NewWriter(*comm.GetSock(i))
					comm.SendVector(writer, i, b)
					comm.SendMatrix(writer, i, A)
					comm.SendVector(writer, i, skShares[i])
					comm.SendBytesSliceMap(writer, i, seeds)
					comm.SendBytesMap(writer, i, MACKeys[i])
				}(i)
			}
		}
		sendWg.Wait()
		genEnd = time.Now()
	} else {
		reader := bufio.NewReader(*comm.GetSock(sign.TrustedDealerID))
		b = comm.RecvVector(reader, sign.TrustedDealerID, sign.M)
		A = comm.RecvMatrix(reader, sign.TrustedDealerID, sign.M)
		party.SkShare = comm.RecvVector(reader, sign.TrustedDealerID, sign.N)
		party.Seed = comm.RecvBytesSliceMap(reader, sign.TrustedDealerID)
		party.MACKeys = comm.RecvBytesMap(reader, sign.TrustedDealerID)
	}

	time.Sleep(time.Second * 5)
	// SIGNATURE ROUND 1
	sid := 1
	PRFKey := "PRFKey"

	r.NTT(lagrangeCoeffs[partyID], lagrangeCoeffs[partyID])
	r.MForm(lagrangeCoeffs[partyID], lagrangeCoeffs[partyID])
	party.Lambda = lagrangeCoeffs[partyID]

	fmt.Printf("Timestamp before Sign Round 1 compute: %s\n", time.Now().Format("15:04:05.000000"))
	start := time.Now()
	D[partyID], MACs[partyID] = party.SignRound1(A, sid, []byte(PRFKey), T)
	signRound1Duration = time.Since(start)
	log.Println("Completed R1")

	signRound1Start = time.Now()
	// Concurrently send and receive data
	var round1Wg sync.WaitGroup
	for i := 0; i < sign.K; i++ {
		if i != partyID {
			round1Wg.Add(2)
			go func(i int) {
				defer round1Wg.Done()
				writer := bufio.NewWriter(*comm.GetSock(i))
				comm.SendMatrix(writer, i, D[partyID])
				comm.SendBytesMap(writer, i, MACs[partyID])
			}(i)

			go func(i int) {
				defer round1Wg.Done()
				reader := bufio.NewReader(*comm.GetSock(i))
				D[i] = comm.RecvMatrix(reader, i, sign.M)
				MACs[i] = comm.RecvBytesMap(reader, i)
			}(i)
		}
	}
	round1Wg.Wait()
	signRound1End = time.Now()

	// SIGN ROUND 2
	z := make(map[int]structs.Vector[ring.Poly])

	fmt.Printf("Timestamp before Sign Round 1 verify: %s\n", time.Now().Format("15:04:05.000000"))
	start = time.Now()
	valid, DSum, hash := party.SignRound2Preprocess(A, b, D, MACs, sid, T)
	if !valid {
		log.Fatalf("MAC verification failed for party %d", partyID)
	} else {
		log.Println("Verification passed, moving onto round 2")
	}
	signRound2PreprocessDuration = time.Since(start)

	fmt.Printf("Timestamp before Sign Round 2 compute: %s\n", time.Now().Format("15:04:05.000000"))
	start = time.Now()
	z[partyID] = party.SignRound2(A, b, DSum, sid, mu, T, []byte(PRFKey), hash)
	signRound2Duration = time.Since(start)

	signRound2Start = time.Now()
	if partyID != sign.CombinerID {
		writer := bufio.NewWriter(*comm.GetSock(sign.CombinerID))
		comm.SendVector(writer, sign.CombinerID, z[partyID])
		signRound2End = time.Now()
	} else {
		for i := 0; i < sign.K; i++ {
			if i != sign.CombinerID {
				reader := bufio.NewReader(*comm.GetSock(i))
				z[i] = comm.RecvVector(reader, i, sign.N)
			}
		}
		combinerReceiveEnd = time.Now()

		fmt.Printf("Timestamp before Finalize: %s\n", time.Now().Format("15:04:05.000000"))
		// SIGNATURE FINALIZE
		start := time.Now()
		_, sig, Delta := party.SignFinalize(z, A, b)
		finalizeDuration = time.Since(start)
		combinerFinalizeEnd = time.Now()
		// Verify the signature
		start = time.Now()
		verified := sign.Verify(r, r_xi, r_nu, sig, A, mu, b, party.C, Delta)
		verifyDuration = time.Since(start)
		fmt.Printf("Signature Verification Result: %v\n", verified)
	}

	// Print all durations
	fmt.Println("Setup duration:", setupDuration)
	fmt.Println("Gen duration:", genDuration)
	fmt.Println("Signature Round 1 duration:", signRound1Duration)
	fmt.Println("Signature Round 1 verify duration:", signRound2PreprocessDuration)
	fmt.Println("Signature Round 2 duration:", signRound2Duration)
	fmt.Println("Finalize duration:", finalizeDuration)
	fmt.Println("Verify duration:", verifyDuration)

	// Print timestamps for networking
	fmt.Printf("Gen start timestamp: %s\n", genStart.Format("15:04:05.000000"))
	fmt.Printf("Gen end timestamp: %s\n", genEnd.Format("15:04:05.000000"))

	fmt.Printf("Sign Round 1 sending/receiving start timestamp: %s\n", signRound1Start.Format("15:04:05.000000"))
	fmt.Printf("Sign Round 1 sending/receiving end timestamp: %s\n", signRound1End.Format("15:04:05.000000"))

	fmt.Printf("Sign Round 2 sending/receiving start timestamp: %s\n", signRound2Start.Format("15:04:05.000000"))
	fmt.Printf("Sign Round 2 sending end timestamp: %s\n", signRound2End.Format("15:04:05.000000"))

	if partyID == sign.CombinerID {
		fmt.Printf("Combiner receive end timestamp: %s\n", combinerReceiveEnd.Format("15:04:05.000000"))
		fmt.Printf("Combiner finalize end timestamp: %s\n", combinerFinalizeEnd.Format("15:04:05.000000"))
	}
}
