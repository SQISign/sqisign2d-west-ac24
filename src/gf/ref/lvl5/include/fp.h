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


/**********************                                    ***********************/
/********************** The below should probably be moved ***********************/
/**********************                                    ***********************/

digit_t mp_shiftr(digit_t* x, const unsigned int shift, const unsigned int nwords);
void multiple_mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords);
void mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords);
void mp_add(digit_t* c, const digit_t* a, const digit_t* b, const unsigned int nwords);
void MUL(digit_t* out, const digit_t a, const digit_t b);

/********************** Constant-time unsigned comparisons ***********************/

// The following functions return 1 (TRUE) if condition is true, 0 (FALSE) otherwise

static inline unsigned int is_digit_nonzero_ct(digit_t x)
{ // Is x != 0?
    return (unsigned int)((x | (0 - x)) >> (RADIX - 1));
}

static inline unsigned int is_digit_zero_ct(digit_t x)
{ // Is x = 0?
    return (unsigned int)(1 ^ is_digit_nonzero_ct(x));
}

static inline unsigned int is_digit_lessthan_ct(digit_t x, digit_t y)
{ // Is x < y?
    return (unsigned int)((x ^ ((x ^ y) | ((x - y) ^ y))) >> (RADIX - 1));
}

/********************** Platform-independent macros for digit-size operations **********************/

// Digit addition with carry
#define ADDC(sumOut, carryOut, addend1, addend2, carryIn)                                         \
    { digit_t tempReg = (addend1) + (digit_t)(carryIn);                                           \
    (sumOut) = (addend2) + tempReg;                                                               \
    (carryOut) = (is_digit_lessthan_ct(tempReg, (digit_t)(carryIn)) | is_digit_lessthan_ct((sumOut), tempReg)); }

// Digit subtraction with borrow
#define SUBC(differenceOut, borrowOut, minuend, subtrahend, borrowIn)                             \
    { digit_t tempReg = (minuend) - (subtrahend);                                                 \
    unsigned int borrowReg = (is_digit_lessthan_ct((minuend), (subtrahend)) | ((borrowIn) & is_digit_zero_ct(tempReg)));  \
    (differenceOut) = tempReg - (digit_t)(borrowIn);                                              \
    (borrowOut) = borrowReg; }

// Shift right with flexible datatype
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << (DigitSize - (shift)));

// Digit shift left
#define SHIFTL(highIn, lowIn, shift, shiftOut, DigitSize)                                         \
    (shiftOut) = ((highIn) << (shift)) ^ ((lowIn) >> (RADIX - (shift)));

#endif
