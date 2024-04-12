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

#include "gf5248.h"

// Type for elements of GF(p)
#define fp_t gf5248

// Constants (Assumed to be in Montgomery form)
#define ZERO gf5248_ZERO
#define ONE gf5248_ONE

// Operations in fp
#define fp_neg gf5248_neg
#define fp_add gf5248_add
#define fp_sub gf5248_sub
#define fp_mul gf5248_mul
#define fp_sqr gf5248_square
#define fp_half gf5248_half

// Conditional swapping
#define fp_swap gf5248_cswap

// Comparisons for fp elements
#define fp_is_zero gf5248_iszero
#define fp_is_equal gf5248_equals

// Set a uint32 to an Fp value
#define fp_set_small gf5248_set_small

// For zero and one, we use predefined constants
void fp_set_zero(fp_t * a);
void fp_set_one(fp_t * a);

// Copy a value
void fp_w64(fp_t *out, const uint64_t data[4]);
void fp_copy(fp_t * out, const fp_t * a);

// Encoding and decoding of bytes
#define fp_encode gf5248_encode
#define fp_decode gf5248_decode
#define fp_decode_reduce gf5248_decode_reduce

// Functions defined in low level code but with different API
void fp_inv(fp_t * a);
void fp_sqrt(fp_t * a);
bool fp_is_square(const fp_t * a);

// TODO: I believe these functions should not be available
// I believe the API should exist without the user thinking
// about or knowing of the Montgomery form. If we need to set
// small integers then we have `gf5248_set_small` which takes
// and integer and internally does the conversion for us.
//
// KILL these functions
// fp_to_mont()
// fp_from_mont()
// fp_set_one_mont()

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
