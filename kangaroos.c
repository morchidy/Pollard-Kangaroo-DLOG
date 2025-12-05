#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mul11585.h"

/* generator g = 4398046511104 = 2^42 in num128 form */
static const num128 G_GEN = {.t = {4398046511104ULL, 0ULL}};

/* return 1 as num128 */
static inline num128 num128_one(void)
{
    num128 r = {.t = {1ULL, 0ULL}};
    return r;
}

/*
 * Compute g^x in the subgroup G modulo (2^115 - 85),
 * using binary (square-and-multiply) exponentiation.
 */
num128 gexp(uint64_t x)
{
    num128 result = num128_one();  // result = 1
    num128 base   = G_GEN;         // base = g

    while (x > 0) {
        if (x & 1ULL) { //if lowest bit of x is 1
            result = mul11585(result, base);
        }
        base = mul11585(base, base);
        x >>= 1;
    }

    return result;
}

/* Hash table for distinguished points */
typedef struct hash_entry {
    num128 point;
    uint64_t exponent;
    int is_tame;
    struct hash_entry *next;
} hash_entry;

#define HASH_SIZE 1024
static hash_entry *hash_table[HASH_SIZE] = {0};

/* Simple hash function for 128-bit numbers */
static uint32_t hash128(num128 x)
{
    return (x.t[0] ^ x.t[1]) % HASH_SIZE;
}

/* Add a distinguished point to hash table */
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

/* Look for a point in hash table, returns exponent if found, -1 otherwise */
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

/* Clear hash table */
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

/* Check if point is distinguished (last 26 bits are 0) */
static int is_distinguished(num128 x)
{
    /* Check last 26 bits are 0 */
    return (x.t[0] & 0x3FFFFFF) == 0;
}

/* Jump function: deterministic mapping based on point */
static void jump(num128 *point, uint64_t *exponent_sum, 
                 const num128 *jump_powers, const uint64_t *jump_sizes, int k)
{
    /* Simple hash to select jump */
    uint32_t h = (point->t[0] ^ point->t[1]) % k;
    
    *point = mul11585(*point, jump_powers[h]);
    *exponent_sum += jump_sizes[h];
}

/*
 * Pollard's Kangaroo algorithm for discrete logarithm in interval [0, 2^64-1]
 */
uint64_t dlog64(num128 target)
{
    const uint64_t W = 0xFFFFFFFFFFFFFFFFULL;  /* 2^64 - 1 */
    const uint64_t W_half = W / 2;
    
    /* Parameters */
    const int k = 32;                    /* log2(W)/2 */
    const uint64_t mu = 1ULL << 31;      /* sqrt(W)/2 = 2^31 */
    const double d_prob = 64.0 / (1ULL << 32); /* log2(W)/sqrt(W) */
    
    /* Initialize random seed */
    srand(time(NULL));
    
    /* Create jump table */
    uint64_t jump_sizes[k];
    num128 jump_powers[k];
    
    for (int i = 0; i < k; i++) {
        /* Jump sizes around mu */
        jump_sizes[i] = mu + (rand() % (mu/10)) - (mu/20);
        jump_powers[i] = gexp(jump_sizes[i]);
    }
    
    /* Initialize kangaroos */
    num128 tame = gexp(W_half);      /* Tame starts at g^(W/2) */
    uint64_t tame_exp = W_half;      /* Known exponent of tame */
    
    num128 wild = target;            /* Wild starts at target */
    uint64_t wild_exp = 0;           /* Unknown exponent, track jumps */
    
    uint64_t iterations = 0;
    int found_tame, found_wild;
    int64_t other_exp;
    
    while (1) {
        iterations++;
        
        /* Tame kangaroo jump */
        jump(&tame, &tame_exp, jump_powers, jump_sizes, k);
        
        /* Check if tame landed on distinguished point */
        if (is_distinguished(tame)) {
            other_exp = hash_lookup(tame, 1, &found_wild);
            if (other_exp != -1 && found_wild == 0) {
                /* Collision with wild kangaroo! */
                uint64_t dlog = (tame_exp > other_exp) ? 
                               (tame_exp - other_exp) : 
                               (other_exp - tame_exp);
                hash_clear();
                return dlog;
            }
            hash_add(tame, tame_exp, 1);
        }
        
        /* Wild kangaroo jump */
        jump(&wild, &wild_exp, jump_powers, jump_sizes, k);
        
        /* Check if wild landed on distinguished point */
        if (is_distinguished(wild)) {
            other_exp = hash_lookup(wild, 0, &found_tame);
            if (other_exp != -1 && found_tame == 1) {
                /* Collision with tame kangaroo! */
                uint64_t dlog = (other_exp > wild_exp) ?
                               (other_exp - wild_exp) :
                               (wild_exp - other_exp);
                hash_clear();
                return dlog;
            }
            hash_add(wild, wild_exp, 0);
        }
        
        /* Safety check - shouldn't happen with good parameters */
        if (iterations > (1ULL << 35)) {  /* ~34 billion iterations */
            fprintf(stderr, "Error: Too many iterations\n");
            hash_clear();
            return 0;
        }
        
        /* Progress indicator */
        if ((iterations & 0xFFFFFFF) == 0) {
            printf("Iterations: %lu\n", iterations);
        }
    }
}

/* Helper function to convert hex string to num128 */
static num128 hex_to_num128(const char *hex)
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
        
        /* Shift and add digit */
        result.t[1] = (result.t[1] << 4) | (result.t[0] >> 60);
        result.t[0] = (result.t[0] << 4) | digit;
    }
    
    return result;
}

int main(void)
{
    printf("Testing gexp function:\n");
    
    num128 r;
    r = gexp(257ULL);
    printf("g^257 = ");
    print_num128(r);
    printf("\n");
    
    r = gexp(112123123412345ULL);
    printf("g^112123123412345 = ");
    print_num128(r);
    printf("\n");
    
    r = gexp(18014398509482143ULL);
    printf("g^18014398509482143 = ");
    print_num128(r);
    printf("\n");
    
    printf("\nTesting dlog64:\n");
    
    /* Test with a known value first */
    printf("\nTest 1: Known exponent (257)\n");
    num128 h1 = gexp(257ULL);
    printf("Target (g^257): ");
    print_num128(h1);
    printf("\n");
    
    clock_t start = clock();
    uint64_t result1 = dlog64(h1);
    clock_t end = clock();
    double time1 = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("Computed exponent: %lu\n", result1);
    printf("Time: %.2f seconds\n", time1);
    printf("Correct: %s\n", result1 == 257ULL ? "Yes" : "No");
    
    /* Test with the given target from TP */
    printf("\nTest 2: Target from TP\n");
    num128 target = hex_to_num128("71AC72AF7B138B6263BF2908A7B09");
    printf("Target: ");
    print_num128(target);
    printf("\n");
    
    start = clock();
    uint64_t result2 = dlog64(target);
    end = clock();
    double time2 = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("Computed exponent: %lu (0x%lX)\n", result2, result2);
    printf("Time: %.2f seconds\n", time2);
    
    /* Verify the result */
    num128 check = gexp(result2);
    printf("Check g^result: ");
    print_num128(check);
    printf("\n");
    if(check.t[0] == target.t[0] && check.t[1] == target.t[1]){
        printf("Matches target: Yes\n");
        printf("Discret log computed succesfuly!\n");
    }
    else {
        printf("Matches target: No\n");
        printf("Failed to compute Discret log!\n");
    }
    // printf("Matches target: %s\n", 
    //        (check.t[0] == target.t[0] && check.t[1] == target.t[1]) ? "Yes" : "No");
    
    return 0;
}