#ifndef FP2_H
#define FP2_H

#include "fp.h"
#include "mp.h" // TODO this should be its own import rather than hiding in here...
#include <stdio.h>

// Structure for representing elements in GF(p^2)
typedef struct fp2_t {
    fp_t re, im;
} fp2_t;

static inline void 
fp2_set_small(fp2_t* x, const uint32_t val)
{
    fp_set_small(&(x->re), val);
    fp_set_zero(&(x->im));
}

static inline void
fp2_set_zero(fp2_t* x)
{
    fp_set_zero(&(x->re));
    fp_set_zero(&(x->im));
}

static inline void
fp2_set_one(fp2_t* x)
{
    fp_set_one(&(x->re));
    fp_set_zero(&(x->im));
}

static inline bool
fp2_is_equal(const fp2_t* a, const fp2_t* b)
{ // Compare two GF(p^2) elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise

    return fp_is_equal(&(a->re), &(b->re)) & fp_is_equal(&(a->im), &(b->im));
}

static inline bool
fp2_is_zero(const fp2_t* a)
{ // Is a GF(p^2) element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise

    return fp_is_zero(&(a->re)) & fp_is_zero(&(a->im));
}

static inline bool
fp2_is_one(const fp2_t* a)
{ // Is a GF(p^2) element one?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    return fp_is_equal(&(a->re), &ONE) & fp_is_zero(&(a->im));
}

static inline void
fp2_half(fp2_t* x, const fp2_t* y)
{
    fp_half(&(x->re), &(y->re));
    fp_half(&(x->im), &(y->im));
}

static inline void 
fp2_add(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_add(&(x->re), &(y->re), &(z->re));
    fp_add(&(x->im), &(y->im), &(z->im));
}

static inline void 
fp2_sub(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_sub(&(x->re), &(y->re), &(z->re));
    fp_sub(&(x->im), &(y->im), &(z->im));
}

static inline void
fp2_neg(fp2_t* x, const fp2_t* y)
{
    fp_neg(&(x->re), &(y->re));
    fp_neg(&(x->im), &(y->im));
}

static inline void
fp2_mul(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_t t0, t1;

    fp_add(&t0, &(y->re), &(y->im));
    fp_add(&t1, &(z->re), &(z->im));
    fp_mul(&t0, &t0, &t1);
    fp_mul(&t1, &(y->im), &(z->im));
    fp_mul(&(x->re), &(y->re), &(z->re));
    fp_sub(&(x->im), &t0, &t1);
    fp_sub(&(x->im), &(x->im), &(x->re));
    fp_sub(&(x->re), &(x->re), &t1);
}

static inline void
fp2_sqr(fp2_t* x, const fp2_t* y)
{
    fp_t sum, diff;

    fp_add(&sum, &(y->re), &(y->im));
    fp_sub(&diff, &(y->re), &(y->im));
    fp_mul(&(x->im), &(y->re), &(y->im));
    fp_add(&(x->im), &(x->im), &(x->im));
    fp_mul(&(x->re), &sum, &diff);
}

static inline void
fp2_cswap(fp2_t *a, fp2_t *b, uint32_t ctl){
    fp_swap(&(a->re), &(b->re), ctl);
    fp_swap(&(a->im), &(b->im), ctl);
}

static inline void 
fp2_copy(fp2_t* x, const fp2_t* y){
    *x = *y;
}

// New functions
void fp2_from_w64(fp2_t* out, const uint64_t data[2][NWORDS_FIELD]);
void fp2_to_w64(uint64_t data[2][NWORDS_FIELD], const fp2_t* a);
void fp2_encode(void *dst, const fp2_t *a);
uint32_t fp2_decode(fp2_t *d, const void *src);
fp2_t fp2_non_residue();
void fp2_inv(fp2_t* x);
bool fp2_is_square(const fp2_t* x);
void fp2_frob(fp2_t* x, const fp2_t* y);
void fp2_sqrt(fp2_t* x);
int fp2_cmp(fp2_t* x, fp2_t* y);
void fp2_batched_inv(fp2_t *x,int len);
void fp2_pow(fp2_t *out,const fp2_t * x,const uint64_t *exp,const int size);
void fp2_print(char *name, fp2_t const a);

#endif
