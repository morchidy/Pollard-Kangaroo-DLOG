#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mul11585.h"

// Configurable parameters
#define K_VALUE 32
#define MU_VALUE 2147483648
#define D_BITS 28
#define START_FRAC 0.5

/* generator g = 4398046511104 = 2^42 in num128 form */
static const num128 G_GEN = {.t = {4398046511104ULL, 0ULL}};

static inline num128 num128_one(void)
{
    num128 r = {.t = {1ULL, 0ULL}};
    return r;
}

num128 gexp(uint64_t x)
{
    num128 result = num128_one();
    num128 base   = G_GEN;

    while (x > 0) {
        if (x & 1ULL) {
            result = mul11585(result, base);
        }
        base = mul11585(base, base);
        x >>= 1;
    }
    return result;
}

/* Hash table implementation */
typedef struct hash_entry {
    num128 point;
    uint64_t exponent;
    int is_tame;
    struct hash_entry *next;
} hash_entry;

#define HASH_SIZE 1024
static hash_entry *hash_table[HASH_SIZE] = {0};

static uint32_t hash128(num128 x)
{
    return (x.t[0] ^ x.t[1]) % HASH_SIZE;
}

static void hash_add(num128 point, uint64_t exponent, int is_tame)
{
    uint32_t idx = hash128(point);
    hash_entry *entry = malloc(sizeof(hash_entry));
    entry->point = point;
    entry->exponent = exponent;
    entry->is_tame = is_tame;
    entry->next = hash_table[idx];
    hash_table[idx] = entry;
}

static int64_t hash_lookup(num128 point, int is_tame, int *found_is_tame)
{
    uint32_t idx = hash128(point);
    hash_entry *entry = hash_table[idx];
    
    while (entry) {
        if (entry->point.t[0] == point.t[0] && entry->point.t[1] == point.t[1]) {
            *found_is_tame = entry->is_tame;
            return entry->exponent;
        }
        entry = entry->next;
    }
    return -1;
}

static void hash_clear(void)
{
    for (int i = 0; i < HASH_SIZE; i++) {
        hash_entry *entry = hash_table[i];
        while (entry) {
            hash_entry *next = entry->next;
            free(entry);
            entry = next;
        }
        hash_table[i] = NULL;
    }
}

/* Configurable distinguished points */
static int is_distinguished(num128 x)
{
    /* Check last D_BITS bits are 0 */
    return (x.t[0] & ((1UL << D_BITS) - 1)) == 0;
}

static void jump(num128 *point, uint64_t *exponent_sum, 
                 const num128 *jump_powers, const uint64_t *jump_sizes, int k)
{
    uint32_t h = (point->t[0] ^ point->t[1]) % k;
    *point = mul11585(*point, jump_powers[h]);
    *exponent_sum += jump_sizes[h];
}

uint64_t dlog64_configurable(num128 target)
{
    const uint64_t W = 0xFFFFFFFFFFFFFFFFULL;
    const uint64_t W_half = (uint64_t)(W * START_FRAC);
    
    const int k = K_VALUE;
    const uint64_t mu = MU_VALUE;
    
    srand(time(NULL));
    
    uint64_t jump_sizes[k];
    num128 jump_powers[k];
    
    for (int i = 0; i < k; i++) {
        jump_sizes[i] = mu + (rand() % (mu/10)) - (mu/20);
        jump_powers[i] = gexp(jump_sizes[i]);
    }
    
    num128 tame = gexp(W_half);
    uint64_t tame_exp = W_half;
    
    num128 wild = target;
    uint64_t wild_exp = 0;
    
    uint64_t iterations = 0;
    int found_tame, found_wild;
    int64_t other_exp;
    
    while (1) {
        iterations++;
        
        jump(&tame, &tame_exp, jump_powers, jump_sizes, k);
        
        if (is_distinguished(tame)) {
            other_exp = hash_lookup(tame, 1, &found_wild);
            if (other_exp != -1 && found_wild == 0) {
                uint64_t dlog = (tame_exp > other_exp) ? 
                               (tame_exp - other_exp) : 
                               (other_exp - tame_exp);
                hash_clear();
                return dlog;
            }
            hash_add(tame, tame_exp, 1);
        }
        
        jump(&wild, &wild_exp, jump_powers, jump_sizes, k);
        
        if (is_distinguished(wild)) {
            other_exp = hash_lookup(wild, 0, &found_tame);
            if (other_exp != -1 && found_tame == 1) {
                uint64_t dlog = (other_exp > wild_exp) ?
                               (other_exp - wild_exp) :
                               (wild_exp - other_exp);
                hash_clear();
                return dlog;
            }
            hash_add(wild, wild_exp, 0);
        }
        
        if (iterations > (1ULL << 35)) {
            fprintf(stderr, "Error: Too many iterations\n");
            hash_clear();
            return 0;
        }
    }
}

num128 hex_to_num128(const char *hex)
{
    num128 result = {0};
    int len = strlen(hex);
    
    for (int i = 0; i < len; i++) {
        char c = hex[i];
        uint64_t digit;
        
        if (c >= '0' && c <= '9') digit = c - '0';
        else if (c >= 'A' && c <= 'F') digit = c - 'A' + 10;
        else if (c >= 'a' && c <= 'f') digit = c - 'a' + 10;
        else continue;
        
        result.t[1] = (result.t[1] << 4) | (result.t[0] >> 60);
        result.t[0] = (result.t[0] << 4) | digit;
    }
    
    return result;
}

int main(void)
{
    printf("\n=== Parameter Test Configuration ===\n");
    printf("K (jump categories): %d\n", K_VALUE);
    printf("μ (avg jump size): %lu (2^%d)\n", MU_VALUE, (int)(log2(MU_VALUE)));
    printf("D (distinguished bits): %d (prob: 2^-%d)\n", D_BITS, D_BITS);
    printf("Start fraction: %.2f\n", START_FRAC);
    
    num128 target = hex_to_num128("71AC72AF7B138B6263BF2908A7B09");
    printf("Target: 71AC72AF7B138B6263BF2908A7B09\n");
    
    clock_t start = clock();
    uint64_t result = dlog64_configurable(target);
    clock_t end = clock();
    double time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("\nResult: %lu (0x%lX)\n", result, result);
    printf("Time: %.2f seconds\n", time_elapsed);
    
    num128 check = gexp(result);
    printf("Verification: %s\n", 
           (check.t[0] == target.t[0] && check.t[1] == target.t[1]) ? "SUCCESS" : "FAILED");
    
    return 0;
}
