#ifndef FP_H
#define FP_H

// Include statements
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <tutil.h>
#include <fp_constants.h>

#include "gf27500.h"

// Type for elements of GF(p)
#define fp_t gf27500

// Constants (Assumed to be in Montgomery form)
#define ZERO gf27500_ZERO
#define ONE gf27500_ONE

// Operations in fp
#define fp_neg gf27500_neg
#define fp_add gf27500_add
#define fp_sub gf27500_sub
#define fp_mul gf27500_mul
#define fp_sqr gf27500_square
#define fp_half gf27500_half

// Conditional swapping
#define fp_swap gf27500_cswap

// Comparisons for fp elements
#define fp_is_zero gf27500_iszero
#define fp_is_equal gf27500_equals

// Set a uint32 to an Fp value
#define fp_set_small gf27500_set_small

// Copy a convert to and from little endian u64 words
#define fp_to_w64 gf27500_to_w64
#define fp_from_w64 gf27500_from_w64

// Encoding and decoding of bytes
#define fp_encode gf27500_encode
#define fp_decode gf27500_decode
#define fp_decode_reduce gf27500_decode_reduce

// These functions are essentially useless because we can just
// use = for the shallow copies we need, but they're here for
// now until we do a larger refactoring
static inline void
fp_copy(fp_t * out, const fp_t * a)
{
    *out = *a;
}

static inline void
fp_set_zero(fp_t * a)
{
    *a = ZERO;
}

static inline void
fp_set_one(fp_t * a)
{
    *a = ONE;
}

// Functions defined in low level code but with different API
void fp_inv(fp_t * a);
void fp_sqrt(fp_t * a);
bool fp_is_square(const fp_t * a);

// TODO
extern const uint64_t NQR_TABLE[20][2][NWORDS_FIELD];
extern const uint64_t Z_NQR_TABLE[20][2][NWORDS_FIELD];

#endif
