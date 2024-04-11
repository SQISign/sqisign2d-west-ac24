// Include statements
#include <stdbool.h>
#include <stdint.h>
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

// Comparisons for fp elements
#define fp_is_zero gf5248_iszero
#define fp_equals gf5248_equals

// Set a uint32 to an Fp value
#define fp_set_small gf5248_set_small

// For zero and one, we use predefined constants
void fp_set_zero(fp_t * a);
void fp_set_one(fp_t * a);

// Copy a value
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
