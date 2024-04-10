// Include statements
#include <stdint.h>
#include <stdbool.h>

// Fp definitions
#define NWORDS_FIELD 4

typedef struct fp_t {
  uint64_t w[NWORDS_FIELD];
} fp_t;

// Constants for Montgomery multiplication
// Modulus: p = 5*2^248 - 1
static const fp_t MODULUS = {
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x04FFFFFFFFFFFFFF,
};

// 2*p
static const fp_t MODULUS_X2 = {0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF,
                                   0xFFFFFFFFFFFFFFFF, 0x09FFFFFFFFFFFFFF};

// Constants 0, 1, -1 in Montgomery representation.
static const fp_t ZERO = {0, 0, 0, 0};
static const fp_t ONE = {
    0x0000000000000033,
    0x0000000000000000,
    0x0000000000000000,
    0x0100000000000000,
};

static const fp_t MINUS_ONE = {
    0xFFFFFFFFFFFFFFCC,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x03FFFFFFFFFFFFFF,
};

// 1/2^244 in the field, in Montgomery representation.
static const fp_t INVT244 = {
    0x0000000000001000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
};

// Montgomery representation of 2^256 (i.e. 2^512 mod q).
static const fp_t R2 = {
    0x3333333333333D70,
    0x3333333333333333,
    0x3333333333333333,
    0x0333333333333333,
};

// Most basic prime field operations
void fp_neg(fp_t *out, const fp_t *a);
void fp_add(fp_t *out, const fp_t *a, const fp_t *b);
void fp_sub(fp_t *out, const fp_t *a, const fp_t *b);
void fp_mul(fp_t *out, const fp_t *a, const fp_t *b);
void fp_sqr(fp_t *out, const fp_t *a);
void fp_n_sqr(fp_t * out, const fp_t * x, int n);

void fp_to_mont(fp_t *out, const fp_t *a);
void fp_from_mont(fp_t *out, const fp_t *a);

void fp_copy(fp_t *out, const fp_t *a);

void fp_set_zero(fp_t *out);
void fp_set_one(fp_t *out);

uint32_t fp_is_zero(fp_t *x);
uint32_t fp_equal(fp_t *x, fp_t *y);

static void fp_exp3div4(fp_t* out, const fp_t* a);
void fp_sqrt(fp_t *a);
void fp_inv(fp_t *a);
bool fp_is_square(const fp_t* a);

