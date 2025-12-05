# Pollard's Kangaroo Algorithm for Memoryless Discrete Logarithm Computation

## Topic

Implementation of Pollard's Kangaroo algorithm for memoryless discrete logarithm computation in bounded intervals. The algorithm solves the discrete logarithm problem for elements in a cyclic subgroup $\mathbb{G} < \mathbb{F}^{\times}_{2^{115}-85}$ where the exponent is known to lie in the interval $[0, 2^{64}-1]$.

## Project Architecture
>src/
├── README.md # This file
├── kangaroos.c # Main implementation of Pollard's Kangaroo algorithm
├── mul11585.h # Header file for multiplication in GF(2^115-85)
├── param_test.sh # Bash script for running parameter tests
├── parameter_tester.py # Python script for automated parameter testing
└── parameter_tests/ # Directory for parameter testing
>>├── test_baseline.c # Baseline configuration test
├── test_k_small.c # Test with k=16
├── test_k_large.c # Test with k=64
├── test_mu_small.c # Test with μ=2^30
├── test_mu_large.c # Test with μ=2^32
├── test_d_frequent.c # Test with 24-bit distinguished points
├── test_d_infrequent.c # Test with 28-bit distinguished points
├── test_start_early.c # Test with early starting point (W/4)
├── test_start_late.c # Test with late starting point (3W/4)
└── mul11585.h # Copy of multiplication header



## Key Files

1. **kangaroos.c** - Main implementation containing:
   - `gexp()` function for exponentiation in the subgroup
   - `dlog64()` function implementing Pollard's Kangaroo algorithm
   - Hash table for storing distinguished points
   - Test cases and verification functions

2. **mul11585.h** - Implements multiplication in $\mathbb{F}_{2^{115}-85}$:
   - Specialized modular multiplication for the prime $2^{115}-85$
   - Uses 128-bit integer operations for efficiency
   - Implements reduction using the congruence $2^{115} \equiv 85$

3. **parameter_tests/** - Contains test programs for Question 7:
   - Each file tests a specific parameter configuration
   - Used for analyzing parameter impact on performance

## Compilation Instructions

#### Main Program
```bash
# Compile with full optimizations
gcc -O3 -march=native kangaroos.c -o kangaroos

# Or without architecture-specific optimizations
gcc -O3 kangaroos.c -o kangaroos
```
#### Parameter Tests
```bash
# Run the testing script
python3 parameter_tester.py

# Or run manually
gcc -O3 -march=native parameter_tests/test_baseline.c -o parameter_tests/baseline
gcc -O3 -march=native parameter_tests/test_k_small.c -o parameter_tests/k_small
# ... similarly for other test files
```
## Execution Instructions
#### Main Program
```bash
./kangaroos
```
This runs two tests:

1. Verification with known exponent 257

2. Computation of the discrete logarithm for the target element 0x71AC72AF7B138B6263BF2908A7B09

#### Parameter Testing
```bash
# Run all parameter tests (takes ~15 minutes)
python3 parameter_tester.py

# Or run individual tests
./parameter_tests/baseline
./parameter_tests/k_small
./parameter_tests/k_large
# ... etc.
```
## Requirements
- GCC compiler (supporting 128-bit integers with __int128)

- Python 3 for running parameter tests (optional)

- Standard C library and math library

## Expected Output
- Successful execution should output:

- Verification of gexp() with test vectors

- Successful recovery of exponent 257 (Test 1)

- Successful computation of the target's discrete logarithm (Test 2)

- Verification that $g^{\text{result}}$ matches the target element

Typical running time: 70-140 seconds depending on hardware and parameters.
