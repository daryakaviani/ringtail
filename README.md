# Ringtail

This is a pure Golang implementation of Ringtail, a practical two-round threshold signature scheme from LWE.

**WARNING:** This implementation is an academic proof-of-concept prototype, has not received careful code review, and is not ready for production use.

### Codebase Overview
- `networking/`
    - `networking.go`: Includes the networking stack which allows signers to form peer-to-peer network connections with other parties. Each party concurrently communicates with every other party by serializing and sending its messages through outgoing TCP sockets, while simultaneously receiving and processing incoming messages.
- `primitives/`
    - `hash.go`: Hashes, MACs, PRFs involved in the scheme.
    - `shamir.go`: Shamir secret-sharing for secret key vector.
- `sign/`
    - `config.go`: Parameters for concrete instantiation.
    - `local.go`: Locally runs the scheme on a single machine for a given number of parties.
    - `sign.go`: Core functionality of the scheme.
- `utils/`
    - `utils.go`: Helpers related to NTT and Montgomery conversions, multiplying, and initializing matrices and vectors of ring elements.
    - `utils-naive.go`: This is note used in the current version, but can be used for testing. It implements convolution-based naive ring-element multiplication.
- `main.go`: Run the code with `go run main.go id iters parties` where `id` is the party ID of the signer running the code (use `l` if you want to run the scheme locally), `iters` is the number of iterations to average the latencies over if you are benchmarking (if not, just use 1), and `parties` is the total number of parties. This is currently a full-threshold implementation. For testing a smaller threshold, set the `Threshold` config parameter with a different value, and use `ShamirSecretSharingGeneral`.

### License

Ringtail is licensed under the Apache 2.0 License. See [LICENSE](https://github.com/daryakaviani/ringtail/blob/main/LICENSE).
