#!/bin/bash
# param_test.sh - Manual parameter testing script

echo "Compiling all test configurations..."

# Compile all configurations
gcc -O3 -march=native parameter_tests/test_baseline.c -o parameter_tests/baseline
gcc -O3 -march=native parameter_tests/test_k_small.c -o parameter_tests/k_small
gcc -O3 -march=native parameter_tests/test_k_large.c -o parameter_tests/k_large
gcc -O3 -march=native parameter_tests/test_mu_small.c -o parameter_tests/mu_small
gcc -O3 -march=native parameter_tests/test_mu_large.c -o parameter_tests/mu_large
gcc -O3 -march=native parameter_tests/test_d_frequent.c -o parameter_tests/d_frequent
gcc -O3 -march=native parameter_tests/test_d_infrequent.c -o parameter_tests/d_infrequent
gcc -O3 -march=native parameter_tests/test_start_early.c -o parameter_tests/start_early
gcc -O3 -march=native parameter_tests/test_start_late.c -o parameter_tests/start_late

echo ""
echo "================================================"
echo "Run tests individually and take screenshots:"
echo "================================================"
echo ""
echo "1. Baseline (optimal):"
echo "   $ ./parameter_tests/baseline"
echo ""
echo "2. Small k (k=16):"
echo "   $ ./parameter_tests/k_small"
echo ""
echo "3. Large k (k=64):"
echo "   $ ./parameter_tests/k_large"
echo ""
echo "4. Small μ (μ=2^30):"
echo "   $ ./parameter_tests/mu_small"
echo ""
echo "5. Large μ (μ=2^32):"
echo "   $ ./parameter_tests/mu_large"
echo ""
echo "6. Frequent distinguished points (24 bits):"
echo "   $ ./parameter_tests/d_frequent"
echo ""
echo "7. Infrequent distinguished points (28 bits):"
echo "   $ ./parameter_tests/d_infrequent"
echo ""
echo "8. Early start (W/4):"
echo "   $ ./parameter_tests/start_early"
echo ""
echo "9. Late start (3W/4):"
echo "   $ ./parameter_tests/start_late"
echo ""
echo "================================================"
echo "Note: Each test takes 1-2 minutes to run."
echo "Take screenshots of the output for each configuration."