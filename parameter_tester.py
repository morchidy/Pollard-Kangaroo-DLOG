#!/usr/bin/env python3
"""
Parameter Testing Script for Pollard's Kangaroo Algorithm
Generates multiple test configurations and runs them
"""

import subprocess
import time
import csv
import os
from datetime import datetime

# Base template with configurable parameters
BASE_TEMPLATE = '''#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mul11585.h"

// Configurable parameters
#define K_VALUE {K_VALUE}
#define MU_VALUE {MU_VALUE}
#define D_BITS {D_BITS}
#define START_FRAC {START_FRAC}

/* generator g = 4398046511104 = 2^42 in num128 form */
static const num128 G_GEN = {{.t = {{4398046511104ULL, 0ULL}}}};

static inline num128 num128_one(void)
{{
    num128 r = {{.t = {{1ULL, 0ULL}}}};
    return r;
}}

num128 gexp(uint64_t x)
{{
    num128 result = num128_one();
    num128 base   = G_GEN;

    while (x > 0) {{
        if (x & 1ULL) {{
            result = mul11585(result, base);
        }}
        base = mul11585(base, base);
        x >>= 1;
    }}
    return result;
}}

/* Hash table implementation */
typedef struct hash_entry {{
    num128 point;
    uint64_t exponent;
    int is_tame;
    struct hash_entry *next;
}} hash_entry;

#define HASH_SIZE 1024
static hash_entry *hash_table[HASH_SIZE] = {{0}};

static uint32_t hash128(num128 x)
{{
    return (x.t[0] ^ x.t[1]) % HASH_SIZE;
}}

static void hash_add(num128 point, uint64_t exponent, int is_tame)
{{
    uint32_t idx = hash128(point);
    hash_entry *entry = malloc(sizeof(hash_entry));
    entry->point = point;
    entry->exponent = exponent;
    entry->is_tame = is_tame;
    entry->next = hash_table[idx];
    hash_table[idx] = entry;
}}

static int64_t hash_lookup(num128 point, int is_tame, int *found_is_tame)
{{
    uint32_t idx = hash128(point);
    hash_entry *entry = hash_table[idx];
    
    while (entry) {{
        if (entry->point.t[0] == point.t[0] && entry->point.t[1] == point.t[1]) {{
            *found_is_tame = entry->is_tame;
            return entry->exponent;
        }}
        entry = entry->next;
    }}
    return -1;
}}

static void hash_clear(void)
{{
    for (int i = 0; i < HASH_SIZE; i++) {{
        hash_entry *entry = hash_table[i];
        while (entry) {{
            hash_entry *next = entry->next;
            free(entry);
            entry = next;
        }}
        hash_table[i] = NULL;
    }}
}}

/* Configurable distinguished points */
static int is_distinguished(num128 x)
{{
    /* Check last D_BITS bits are 0 */
    return (x.t[0] & ((1UL << D_BITS) - 1)) == 0;
}}

static void jump(num128 *point, uint64_t *exponent_sum, 
                 const num128 *jump_powers, const uint64_t *jump_sizes, int k)
{{
    uint32_t h = (point->t[0] ^ point->t[1]) % k;
    *point = mul11585(*point, jump_powers[h]);
    *exponent_sum += jump_sizes[h];
}}

uint64_t dlog64_configurable(num128 target)
{{
    const uint64_t W = 0xFFFFFFFFFFFFFFFFULL;
    const uint64_t W_half = (uint64_t)(W * START_FRAC);
    
    const int k = K_VALUE;
    const uint64_t mu = MU_VALUE;
    
    srand(time(NULL));
    
    uint64_t jump_sizes[k];
    num128 jump_powers[k];
    
    for (int i = 0; i < k; i++) {{
        jump_sizes[i] = mu + (rand() % (mu/10)) - (mu/20);
        jump_powers[i] = gexp(jump_sizes[i]);
    }}
    
    num128 tame = gexp(W_half);
    uint64_t tame_exp = W_half;
    
    num128 wild = target;
    uint64_t wild_exp = 0;
    
    uint64_t iterations = 0;
    int found_tame, found_wild;
    int64_t other_exp;
    
    while (1) {{
        iterations++;
        
        jump(&tame, &tame_exp, jump_powers, jump_sizes, k);
        
        if (is_distinguished(tame)) {{
            other_exp = hash_lookup(tame, 1, &found_wild);
            if (other_exp != -1 && found_wild == 0) {{
                uint64_t dlog = (tame_exp > other_exp) ? 
                               (tame_exp - other_exp) : 
                               (other_exp - tame_exp);
                hash_clear();
                return dlog;
            }}
            hash_add(tame, tame_exp, 1);
        }}
        
        jump(&wild, &wild_exp, jump_powers, jump_sizes, k);
        
        if (is_distinguished(wild)) {{
            other_exp = hash_lookup(wild, 0, &found_tame);
            if (other_exp != -1 && found_tame == 1) {{
                uint64_t dlog = (other_exp > wild_exp) ?
                               (other_exp - wild_exp) :
                               (wild_exp - other_exp);
                hash_clear();
                return dlog;
            }}
            hash_add(wild, wild_exp, 0);
        }}
        
        if (iterations > (1ULL << 35)) {{
            fprintf(stderr, "Error: Too many iterations\\n");
            hash_clear();
            return 0;
        }}
    }}
}}

num128 hex_to_num128(const char *hex)
{{
    num128 result = {{0}};
    int len = strlen(hex);
    
    for (int i = 0; i < len; i++) {{
        char c = hex[i];
        uint64_t digit;
        
        if (c >= '0' && c <= '9') digit = c - '0';
        else if (c >= 'A' && c <= 'F') digit = c - 'A' + 10;
        else if (c >= 'a' && c <= 'f') digit = c - 'a' + 10;
        else continue;
        
        result.t[1] = (result.t[1] << 4) | (result.t[0] >> 60);
        result.t[0] = (result.t[0] << 4) | digit;
    }}
    
    return result;
}}

int main(void)
{{
    printf("\\n=== Parameter Test Configuration ===\\n");
    printf("K (jump categories): %d\\n", K_VALUE);
    printf("μ (avg jump size): %lu (2^%d)\\n", MU_VALUE, (int)(log2(MU_VALUE)));
    printf("D (distinguished bits): %d (prob: 2^-%d)\\n", D_BITS, D_BITS);
    printf("Start fraction: %.2f\\n", START_FRAC);
    
    num128 target = hex_to_num128("71AC72AF7B138B6263BF2908A7B09");
    printf("Target: 71AC72AF7B138B6263BF2908A7B09\\n");
    
    clock_t start = clock();
    uint64_t result = dlog64_configurable(target);
    clock_t end = clock();
    double time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("\\nResult: %lu (0x%lX)\\n", result, result);
    printf("Time: %.2f seconds\\n", time_elapsed);
    
    num128 check = gexp(result);
    printf("Verification: %s\\n", 
           (check.t[0] == target.t[0] && check.t[1] == target.t[1]) ? "SUCCESS" : "FAILED");
    
    return 0;
}}
'''

# Test configurations
test_configs = [
    # Baseline (optimal from theory)
    {
        'name': 'baseline',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 31,  # 2^31
        'D_BITS': 26,
        'START_FRAC': 0.5
    },
    # Test different k values
    {
        'name': 'k_small',
        'K_VALUE': 16,
        'MU_VALUE': 1 << 31,
        'D_BITS': 26,
        'START_FRAC': 0.5
    },
    {
        'name': 'k_large',
        'K_VALUE': 64,
        'MU_VALUE': 1 << 31,
        'D_BITS': 26,
        'START_FRAC': 0.5
    },
    # Test different μ values
    {
        'name': 'mu_small',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 30,  # 2^30
        'D_BITS': 26,
        'START_FRAC': 0.5
    },
    {
        'name': 'mu_large',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 32,  # 2^32
        'D_BITS': 26,
        'START_FRAC': 0.5
    },
    # Test different distinguished point frequencies
    {
        'name': 'd_frequent',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 31,
        'D_BITS': 24,  # More frequent
        'START_FRAC': 0.5
    },
    {
        'name': 'd_infrequent',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 31,
        'D_BITS': 28,  # Less frequent
        'START_FRAC': 0.5
    },
    # Test different starting points
    {
        'name': 'start_early',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 31,
        'D_BITS': 26,
        'START_FRAC': 0.25
    },
    {
        'name': 'start_late',
        'K_VALUE': 32,
        'MU_VALUE': 1 << 31,
        'D_BITS': 26,
        'START_FRAC': 0.75
    },
]

def generate_test_programs():
    """Generate C programs for each configuration"""
    os.makedirs("parameter_tests", exist_ok=True)
    
    for config in test_configs:
        filename = f"parameter_tests/test_{config['name']}.c"
        
        code = BASE_TEMPLATE.format(
            K_VALUE=config['K_VALUE'],
            MU_VALUE=config['MU_VALUE'],
            D_BITS=config['D_BITS'],
            START_FRAC=config['START_FRAC']
        )
        
        with open(filename, 'w') as f:
            f.write(code)
        
        print(f"Generated {filename}")

def compile_and_run_tests():
    """Compile and run all test programs"""
    results = []
    
    for config in test_configs:
        source_file = f"parameter_tests/test_{config['name']}.c"
        executable = f"parameter_tests/{config['name']}"
        
        print(f"\n{'='*60}")
        print(f"Testing: {config['name']}")
        print(f"  K={config['K_VALUE']}, μ={config['MU_VALUE']}, D={config['D_BITS']} bits, start={config['START_FRAC']}")
        
        # Compile
        compile_cmd = ["gcc", "-O3", "-march=native", source_file, "-o", executable]
        compile_result = subprocess.run(compile_cmd, capture_output=True, text=True)
        
        if compile_result.returncode != 0:
            print(f"  Compilation failed: {compile_result.stderr}")
            continue
        
        # Run
        start_time = time.time()
        run_result = subprocess.run([f"./{executable}"], capture_output=True, text=True)
        end_time = time.time()
        
        if run_result.returncode != 0:
            print(f"  Execution failed: {run_result.stderr}")
            continue
        
        # Parse output
        output = run_result.stdout
        time_match = None
        for line in output.split('\n'):
            if "Time:" in line:
                time_match = float(line.split(':')[1].strip().split()[0])
                break
        
        # Extract result
        result_line = None
        for line in output.split('\n'):
            if "Result:" in line:
                result_line = line
                break
        
        # Store results
        result = {
            'config': config['name'],
            'K': config['K_VALUE'],
            'μ': config['MU_VALUE'],
            'D_bits': config['D_BITS'],
            'start_frac': config['START_FRAC'],
            'time': time_match,
            'output': output
        }
        results.append(result)
        
        print(f"  Time: {time_match:.2f} seconds")
        if result_line:
            print(f"  {result_line}")
    
    return results

def create_summary_table(results):
    """Create a summary table of all results"""
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    print(f"{'Configuration':<20} {'K':<5} {'μ':<12} {'D_bits':<8} {'Start':<8} {'Time (s)':<10}")
    print("-"*80)
    
    for result in results:
        print(f"{result['config']:<20} {result['K']:<5} {result['μ']:<12} "
              f"{result['D_bits']:<8} {result['start_frac']:<8.2f} {result['time']:<10.2f}")

def run_single_quick_test():
    """Run a quick test with just 3 configurations for demonstration"""
    quick_configs = [
        test_configs[0],  # baseline
        test_configs[1],  # k_small
        test_configs[3],  # mu_small
    ]
    
    print("\n=== QUICK PARAMETER TEST ===")
    print("Testing 3 configurations for demonstration\n")
    
    for config in quick_configs:
        source_file = f"parameter_tests/test_{config['name']}.c"
        executable = f"parameter_tests/{config['name']}"
        
        print(f"\n--- Testing {config['name']} ---")
        print(f"Parameters: K={config['K_VALUE']}, μ={config['MU_VALUE']}, "
              f"D={config['D_BITS']} bits, start={config['START_FRAC']}")
        
        # Run
        run_result = subprocess.run([f"./{executable}"], capture_output=True, text=True)
        
        # Display relevant output
        for line in run_result.stdout.split('\n'):
            if any(keyword in line for keyword in ['===', 'K', 'μ', 'D', 'Start', 'Time:', 'Result:']):
                print(f"  {line}")
    
    print("\n" + "="*50)
    print("Quick test completed!")
    print("Take screenshots of each configuration's output.")
    print("="*50)

def main():
    """Main function"""
    print("Pollard's Kangaroo Parameter Testing")
    print("="*60)
    
    # Step 1: Generate all test programs
    print("\n1. Generating test programs...")
    generate_test_programs()
    
    # Step 2: Choose testing mode
    print("\n2. Choose testing mode:")
    print("   a) Full test (all 9 configurations)")
    print("   b) Quick test (3 configurations for screenshots)")
    
    choice = input("\nEnter choice (a/b): ").strip().lower()
    
    if choice == 'a':
        # Full test
        print("\nRunning full test suite (this may take ~10-15 minutes)...")
        results = compile_and_run_tests()
        create_summary_table(results)
        
        # Save results to CSV
        with open('parameter_test_results.csv', 'w', newline='') as csvfile:
            fieldnames = ['config', 'K', 'μ', 'D_bits', 'start_frac', 'time']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for result in results:
                writer.writerow({k: result[k] for k in fieldnames})
        
        print("\nResults saved to 'parameter_test_results.csv'")
        
    else:
        # Quick test for screenshots
        print("\nRunning quick test for screenshots...")
        
        # Make sure the programs are compiled
        for config in [test_configs[0], test_configs[1], test_configs[3]]:
            source_file = f"parameter_tests/test_{config['name']}.c"
            executable = f"parameter_tests/{config['name']}"
            subprocess.run(["gcc", "-O3", "-march=native", source_file, "-o", executable], 
                         capture_output=True)
        
        run_single_quick_test()
        
        # Generate instructions for manual testing
        print("\n" + "="*60)
        print("MANUAL TESTING INSTRUCTIONS:")
        print("="*60)
        print("\nTo get clean screenshots, run each test separately:")
        print("\n1. Baseline configuration:")
        print("   $ ./parameter_tests/baseline")
        print("\n2. Small k (k=16):")
        print("   $ ./parameter_tests/k_small")
        print("\n3. Small μ (μ=2^30):")
        print("   $ ./parameter_tests/mu_small")
        print("\n4. (Optional) Other configurations:")
        print("   $ ./parameter_tests/k_large")
        print("   $ ./parameter_tests/mu_large")
        print("   etc...")
        print("\nTake screenshots of each run for your report.")

if __name__ == "__main__":
    main()