#ifndef FP2_H
#define FP2_H

#include "fp.h"
#include <stdio.h>

// Structure for representing elements in GF(p^2)
typedef struct fp2_t {
    fp_t re, im;
} fp2_t;

void fp2_set_small(fp2_t* x, const digit_t val);
void fp2_set_one(fp2_t* x);
void fp2_set_zero(fp2_t* x);
bool fp2_is_zero(const fp2_t* a);
bool fp2_is_equal(const fp2_t* a, const fp2_t* b);
bool fp2_is_one(const fp2_t* a);
void fp2_copy(fp2_t* x, const fp2_t* y);
fp2_t fp2_non_residue();
void fp2_add(fp2_t* x, const fp2_t* y, const fp2_t* z);
void fp2_sub(fp2_t* x, const fp2_t* y, const fp2_t* z);
void fp2_neg(fp2_t* x, const fp2_t* y);
void fp2_mul(fp2_t* x, const fp2_t* y, const fp2_t* z);
void fp2_sqr(fp2_t* x, const fp2_t* y);
void fp2_inv(fp2_t* x);
bool fp2_is_square(const fp2_t* x);
void fp2_frob(fp2_t* x, const fp2_t* y);
void fp2_sqrt(fp2_t* x);
void fp2_tomont(fp2_t* x, const fp2_t* y);
void fp2_frommont(fp2_t* x, const fp2_t* y);
int fp2_cmp(fp2_t* x, fp2_t* y);
void fp2_half(fp2_t* x, const fp2_t* y);
void fp2_batched_inv(fp2_t *x,int len);
void fp2_pow(fp2_t *out,const fp2_t * x,const digit_t *exp,const int size);
void fp2_print(char *name, fp2_t const a);
void digit_print(char *name,digit_t *a);
void fp2_encode(void *dst, const fp2_t *a);
void fp2_cswap(fp2_t *a, fp2_t *b, uint32_t ctl);

extern const uint64_t NQR_TABLE[20][2][NWORDS_FIELD];
extern const uint64_t Z_NQR_TABLE[20][2][NWORDS_FIELD];

#endif
