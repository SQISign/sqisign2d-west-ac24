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

#include "gf65376.h"

// Type for elements of GF(p)
#define fp_t gf65376

// Constants (Assumed to be in Montgomery form)
#define ZERO gf65376_ZERO
#define ONE gf65376_ONE

// Operations in fp
#define fp_neg gf65376_neg
#define fp_add gf65376_add
#define fp_sub gf65376_sub
#define fp_mul gf65376_mul
#define fp_sqr gf65376_square
#define fp_half gf65376_half

// Conditional swapping
#define fp_swap gf65376_cswap

// Comparisons for fp elements
#define fp_is_zero gf65376_iszero
#define fp_is_equal gf65376_equals

// Set a uint32 to an Fp value
#define fp_set_small gf65376_set_small

// Copy a convert to and from little endian u64 words
#define fp_to_w64 gf65376_to_w64
#define fp_from_w64 gf65376_from_w64

// Encoding and decoding of bytes
#define fp_encode gf65376_encode
#define fp_decode gf65376_decode
#define fp_decode_reduce gf65376_decode_reduce

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
