package networking

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"log"
	"net"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

type Communicator interface {
	Send(dst int, msg []byte) (int, error)
	Recv(src int) ([]byte, int, error)
	Close() error
}

type P2PComm struct {
	Socks map[int]*net.Conn
	Rank  int
	mu    sync.Mutex // Added mutex for safe concurrent access
}

func (comm *P2PComm) SetSock(key int, conn *net.Conn) {
	comm.mu.Lock()
	defer comm.mu.Unlock()
	comm.Socks[key] = conn
}

func (comm *P2PComm) GetSock(key int) *net.Conn {
	comm.mu.Lock()
	defer comm.mu.Unlock()
	return comm.Socks[key]
}

func (comm *P2PComm) SendBytes(writer *bufio.Writer, dst int, msg []byte) (int, error) {
	length := make([]byte, 4)
	binary.BigEndian.PutUint32(length, uint32(len(msg)))

	n, err := writer.Write(length)
	if err != nil {
		return 0, err
	}
	totalBytesSent := n

	for len(msg) > 0 {
		n, err = writer.Write(msg)
		if err != nil {
			return 0, err
		}
		totalBytesSent += n
		msg = msg[n:]
	}

	err = writer.Flush()
	if err != nil {
		return 0, err
	}

	return totalBytesSent, nil
}

func (comm *P2PComm) Recv(reader *bufio.Reader, src int) ([]byte, int, error) {
	lengthBuf := make([]byte, 4)
	totalBytesRead := 0
	for totalBytesRead < 4 {
		n, err := reader.Read(lengthBuf[totalBytesRead:])
		if err != nil {
			return nil, totalBytesRead, err
		}
		totalBytesRead += n
	}
	length := binary.BigEndian.Uint32(lengthBuf)

	data := make([]byte, length)
	bytesRead := 0
	for bytesRead < int(length) {
		n, err := reader.Read(data[bytesRead:])
		if err != nil {
			return nil, totalBytesRead, err
		}
		bytesRead += n
		totalBytesRead += n
	}

	return data, totalBytesRead, nil
}

func (comm *P2PComm) Close() error {
	comm.mu.Lock()
	defer comm.mu.Unlock()
	for _, sock := range comm.Socks {
		err := (*sock).Close()
		if err != nil {
			return err
		}
	}
	return nil
}

func (comm *P2PComm) SendVector(writer *bufio.Writer, dst int, msg structs.Vector[ring.Poly]) {
	if _, err := msg.WriteTo(writer); err != nil {
		log.Fatalf("Failed to write vector: %v", err)
	}

	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush writer: %v", err)
	}
}

func (comm *P2PComm) RecvVector(reader *bufio.Reader, src int, length int) structs.Vector[ring.Poly] {
	vec := make(structs.Vector[ring.Poly], length)
	if _, err := vec.ReadFrom(reader); err != nil {
		log.Fatalf("Failed to read vector: %v", err)
	}
	return vec
}

func (comm *P2PComm) SendMatrix(writer *bufio.Writer, dst int, msg structs.Matrix[ring.Poly]) {
	if _, err := msg.WriteTo(writer); err != nil {
		log.Fatalf("Error sending matrix: %v", err)
	}

	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush writer: %v", err)
	}
}

func (comm *P2PComm) RecvMatrix(reader *bufio.Reader, src int, length int) structs.Matrix[ring.Poly] {
	matrix := make(structs.Matrix[ring.Poly], length)
	if _, err := matrix.ReadFrom(reader); err != nil {
		log.Fatalf("Failed to read matrix: %v", err)
	}
	return matrix
}

func (comm *P2PComm) SendBytesSlice(writer *bufio.Writer, dst int, data [][]byte) {
	numSlices := uint32(len(data))
	if err := binary.Write(writer, binary.BigEndian, numSlices); err != nil {
		log.Fatalf("Failed to write number of slices: %v", err)
	}

	for _, slice := range data {
		length := uint32(len(slice))
		if err := binary.Write(writer, binary.BigEndian, length); err != nil {
			log.Fatalf("Failed to write slice length: %v", err)
		}

		for len(slice) > 0 {
			n, err := writer.Write(slice)
			if err != nil {
				log.Fatalf("Failed to write slice data: %v", err)
			}
			slice = slice[n:]
		}
	}

	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush writer: %v", err)
	}
}

func (comm *P2PComm) RecvBytesSlice(reader *bufio.Reader, src int) [][]byte {
	var numSlices uint32
	if err := binary.Read(reader, binary.BigEndian, &numSlices); err != nil {
		log.Fatalf("Failed to read number of slices: %v", err)
	}

	data := make([][]byte, numSlices)
	for i := uint32(0); i < numSlices; i++ {
		var length uint32
		if err := binary.Read(reader, binary.BigEndian, &length); err != nil {
			log.Fatalf("Failed to read slice length: %v", err)
		}

		slice := make([]byte, length)
		bytesRead := 0
		for bytesRead < int(length) {
			n, err := reader.Read(slice[bytesRead:])
			if err != nil {
				log.Fatalf("Failed to read slice data: %v", err)
			}
			bytesRead += n
		}

		data[i] = slice
	}

	return data
}

func (comm *P2PComm) SendBytesMap(writer *bufio.Writer, dst int, data map[int][]byte) {
	numEntries := uint32(len(data))
	if err := binary.Write(writer, binary.BigEndian, numEntries); err != nil {
		log.Fatalf("Failed to write number of map entries: %v", err)
	}

	for key, value := range data {
		if err := binary.Write(writer, binary.BigEndian, int32(key)); err != nil {
			log.Fatalf("Failed to write map key: %v", err)
		}

		length := uint32(len(value))
		if err := binary.Write(writer, binary.BigEndian, length); err != nil {
			log.Fatalf("Failed to write value length: %v", err)
		}

		for len(value) > 0 {
			n, err := writer.Write(value)
			if err != nil {
				log.Fatalf("Failed to write value data: %v", err)
			}
			value = value[n:]
		}
	}

	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush writer: %v", err)
	}
}

func (comm *P2PComm) RecvBytesMap(reader *bufio.Reader, src int) map[int][]byte {
	var numEntries uint32
	if err := binary.Read(reader, binary.BigEndian, &numEntries); err != nil {
		log.Fatalf("Failed to read number of map entries: %v", err)
	}

	data := make(map[int][]byte, numEntries)
	for i := uint32(0); i < numEntries; i++ {
		var key int32
		if err := binary.Read(reader, binary.BigEndian, &key); err != nil {
			log.Fatalf("Failed to read map key: %v", err)
		}

		var length uint32
		if err := binary.Read(reader, binary.BigEndian, &length); err != nil {
			log.Fatalf("Failed to read value length: %v", err)
		}

		value := make([]byte, length)
		bytesRead := 0
		for bytesRead < int(length) {
			n, err := reader.Read(value[bytesRead:])
			if err != nil {
				log.Fatalf("Failed to read value data: %v", err)
			}
			bytesRead += n
		}

		data[int(key)] = value
	}

	return data
}

func (comm *P2PComm) SendBytesSliceMap(writer *bufio.Writer, dst int, data map[int][][]byte) {
	numEntries := uint32(len(data))
	if err := binary.Write(writer, binary.BigEndian, numEntries); err != nil {
		log.Fatalf("Failed to write number of map entries: %v", err)
	}

	for key, value := range data {
		if err := binary.Write(writer, binary.BigEndian, int32(key)); err != nil {
			log.Fatalf("Failed to write map key: %v", err)
		}

		numSlices := uint32(len(value))
		if err := binary.Write(writer, binary.BigEndian, numSlices); err != nil {
			log.Fatalf("Failed to write number of slices: %v", err)
		}

		for _, slice := range value {
			length := uint32(len(slice))
			if err := binary.Write(writer, binary.BigEndian, length); err != nil {
				log.Fatalf("Failed to write slice length: %v", err)
			}

			for len(slice) > 0 {
				n, err := writer.Write(slice)
				if err != nil {
					log.Fatalf("Failed to write slice data: %v", err)
				}
				slice = slice[n:]
			}
		}
	}

	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush writer: %v", err)
	}
}

func (comm *P2PComm) RecvBytesSliceMap(reader *bufio.Reader, src int) map[int][][]byte {
	var numEntries uint32
	if err := binary.Read(reader, binary.BigEndian, &numEntries); err != nil {
		log.Fatalf("Failed to read number of map entries: %v", err)
	}

	data := make(map[int][][]byte, numEntries)
	for i := uint32(0); i < numEntries; i++ {
		var key int32
		if err := binary.Read(reader, binary.BigEndian, &key); err != nil {
			log.Fatalf("Failed to read map key: %v", err)
		}

		var numSlices uint32
		if err := binary.Read(reader, binary.BigEndian, &numSlices); err != nil {
			log.Fatalf("Failed to read number of slices: %v", err)
		}

		slices := make([][]byte, numSlices)
		for j := uint32(0); j < numSlices; j++ {
			var length uint32
			if err := binary.Read(reader, binary.BigEndian, &length); err != nil {
				log.Fatalf("Failed to read slice length: %v", err)
			}

			slice := make([]byte, length)
			bytesRead := 0
			for bytesRead < int(length) {
				n, err := reader.Read(slice[bytesRead:])
				if err != nil {
					log.Fatalf("Failed to read slice data: %v", err)
				}
				bytesRead += n
			}

			slices[j] = slice
		}

		data[int(key)] = slices
	}

	return data
}

func ListenTCP(comm *P2PComm, port string, src int) {
	l, err := net.Listen("tcp", "0.0.0.0:"+port)
	if err != nil {
		log.Fatal(err)
		return
	}
	defer l.Close()
	for {
		conn, err := l.Accept()
		if err != nil {
			log.Println(err)
			continue
		}

		comm.SetSock(src, &conn)
		break
	}
}

func DialTCP(comm *P2PComm, dst int, address string) {
	log.Println("Trying to dial", dst, "at address", address)

	// Check if the socket for the destination is already initialized
	sock := comm.GetSock(dst)
	if sock != nil && *sock != nil {
		log.Printf("Socket for destination %d is already initialized.", dst)
		return
	}

	var conn net.Conn
	var err error

	for i := 0; i < 10; i++ {
		conn, err = net.Dial("tcp", address)
		if err == nil {
			break
		}
		log.Printf("Failed to dial TCP to %s: %v, retrying...", address, err)
		time.Sleep(1 * time.Second)
	}

	if err != nil {
		log.Fatalf("Failed to dial TCP to %s after retries: %v", address, err)
		return
	}

	comm.SetSock(dst, &conn)
	log.Println("Successfully connected to", address)
}

func EstablishConnections(wg *sync.WaitGroup, comm *P2PComm, partyID int, totalParties int) {
	defer wg.Done()
	var localWg sync.WaitGroup

	// List of IP addresses of all parties
	ips := []string{"50.18.79.144", "34.226.247.5", "18.199.237.28", "34.244.57.20", "54.248.28.234", "47.129.58.89", "54.232.160.71", "3.107.98.51"}
	for otherID := 0; otherID < totalParties; otherID++ {
		if partyID != otherID {
			localWg.Add(1)
			go func(otherID int) {
				defer localWg.Done()
				port := 6000 + calculatePortOffset(partyID, otherID)
				address := fmt.Sprintf("%s:%d", ips[otherID], port)
				if partyID < otherID {
					ListenTCP(comm, fmt.Sprintf("%d", port), otherID)
				} else {
					DialTCP(comm, otherID, address)
				}
			}(otherID)
		}
	}
	localWg.Wait()
}

func calculatePortOffset(partyID, otherID int) int {
	if partyID < otherID {
		return partyID*100 + otherID
	}
	return otherID*100 + partyID
}
