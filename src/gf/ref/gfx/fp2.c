#include <fp2.h>

extern const uint64_t NQR_TABLE[20][2][NWORDS_FIELD];
extern const uint64_t Z_NQR_TABLE[20][2][NWORDS_FIELD];

/* Arithmetic modulo X^2 + 1 */

void fp2_set_small(fp2_t* x, const digit_t val)
{
    fp_set_small(&(x->re), val);
    fp_set_zero(&(x->im));
}

void fp2_set_one(fp2_t* x)
{
    fp_set_one(&(x->re));
    fp_set_zero(&(x->im));
}

void fp2_set_zero(fp2_t* x)
{
    fp_set_zero(&(x->re));
    fp_set_zero(&(x->im));
}

bool fp2_is_zero(const fp2_t* a)
{ // Is a GF(p^2) element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise

    return fp_is_zero(&(a->re)) & fp_is_zero(&(a->im));
}

bool fp2_is_equal(const fp2_t* a, const fp2_t* b)
{ // Compare two GF(p^2) elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise

    return fp_is_equal(&(a->re), &(b->re)) & fp_is_equal(&(a->im), &(b->im));
}

bool fp2_is_one(const fp2_t* a)
{ // Is a GF(p^2) element one?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    fp2_t t;
    fp2_set_one(&t);
    return fp2_is_equal(&t, a);
}

void
fp2_cswap(fp2_t *a, fp2_t *b, uint32_t ctl){
    fp_cswap(&(a->re), &(b->re), ctl);
    fp_cswap(&(a->im), &(b->im), ctl);
}


void fp2_copy(fp2_t* x, const fp2_t* y)
{
    fp_copy(&(x->re), &(y->re));
    fp_copy(&(x->im), &(y->im));
}

void fp2_encode(void *dst, const fp2_t *a){
	uint8_t *buf = dst;
	fp_encode(buf, &(a->re));
	fp_encode(buf + 8*NWORDS_FIELD, &(a->im));
}

void fp2_decode(fp2_t *d, const void *src){
    const uint8_t *buf = src;
    uint32_t re, im;
    
    fp_decode(&(d->re), buf);
    fp_decode(&(d->im), buf + 8*NWORDS_FIELD);

}

// TODO: this is toxic waste in a gfx file!!!
fp2_t fp2_non_residue()
{ // 2 + i is a quadratic non-residue for p1913
    fp2_t res;

    fp_set_small(&res.re, 2);
    fp_set_one(&res.im);
    return res;
}

void
fp2_half(fp2_t* x, const fp2_t* y)
{
    fp_half(&(x->re), &(y->re));
    fp_half(&(x->im), &(y->im));
}

void fp2_add(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_add(&(x->re), &(y->re), &(z->re));
    fp_add(&(x->im), &(y->im), &(z->im));
}

void fp2_sub(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_sub(&(x->re), &(y->re), &(z->re));
    fp_sub(&(x->im), &(y->im), &(z->im));
}

void fp2_neg(fp2_t* x, const fp2_t* y)
{
    fp_neg(&(x->re), &(y->re));
    fp_neg(&(x->im), &(y->im));
}

void
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

void
fp2_sqr(fp2_t* x, const fp2_t* y)
{
    fp_t sum, diff;

    fp_add(&sum, &(y->re), &(y->im));
    fp_sub(&diff, &(y->re), &(y->im));
    fp_mul(&(x->im), &(y->re), &(y->im));
    fp_add(&(x->im), &(x->im), &(x->im));
    fp_mul(&(x->re), &sum, &diff);
}

void fp2_inv(fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);
    fp_inv(&t0);
    fp_mul(&(x->re), &(x->re), &t0);
    fp_mul(&(x->im), &(x->im), &t0);
    fp_neg(&(x->im), &(x->im));
}

void fp2_batched_inv(fp2_t *x, int len) {

    fp2_t t1[len],t2[len];
    fp2_t inverse;

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    fp2_copy(&t1[0],&x[0]);
    for (int i=1;i<len;i++) {
        fp2_mul(&t1[i],&t1[i-1],&x[i]);
    }

    // inverse = 1/ (x0 * x1 * ... * xn)
    fp2_copy(&inverse,&t1[len-1]);
    fp2_inv(&inverse);

    fp2_copy(&t2[0],&inverse);
    // t2 = 1/ (x0 * x1 * ... * xn), 1/ (x0 * x1 * ... * x(n-1)) , ... , 1/xO
    for (int i=1;i<len;i++) {
        fp2_mul(&t2[i],&t2[i-1],&x[len-i]);
    }

    fp2_copy(&x[0],&t2[len-1]);
    
    for (int i=1;i<len;i++){
        fp2_mul(&x[i],&t1[i-1],&t2[len-i-1]);
    }

}

bool fp2_is_square(const fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);

    return fp_is_square(&t0);
}

void fp2_frob(fp2_t* x, const fp2_t* y)
{
    fp_copy(&(x->re), &(y->re));
    fp_neg(&(x->im), &(y->im));
}

// These functions are never needed with the new API
void fp2_tomont(fp2_t* x, const fp2_t* y)
{ 
    fp_tomont(&(x->re), &(y->re));
    fp_tomont(&(x->im), &(y->im));
}

// These functions are never needed with the new API
void fp2_frommont(fp2_t* x, const fp2_t* y)
{
    fp_frommont(&(x->re), &(y->re));
    fp_frommont(&(x->im), &(y->im));
}

// NOTE: old, non-constant-time implementation. Could be optimized
void fp2_sqrt(fp2_t* x)
{
    fp_t sdelta, re, tmp1, tmp2, im;

    if (fp_is_zero(&(x->im))) {
        if (fp_is_square(&(x->re))) {
            fp_sqrt(&(x->re));
            return;
        } else {
            fp_neg(&(x->im), &(x->re));
            fp_sqrt(&(x->im));
            fp_set_zero(&(x->re));
            return;
        }
    }

    // sdelta = sqrt(re^2 + im^2)
    fp_sqr(&sdelta, &(x->re));
    fp_sqr(&tmp1, &(x->im));
    fp_add(&sdelta, &sdelta, &tmp1);
    fp_sqrt(&sdelta);

    fp_add(&re, &(x->re), &sdelta);
    fp_half(&re, &re);
    fp_copy(&tmp2, &re);

    if (!fp_is_square(&tmp2)) {
        fp_sub(&re, &(x->re), &sdelta);
        fp_half(&re, &re);
    }

    fp_sqrt(&re);
    fp_copy(&im, &re);

    fp_inv(&im);
    fp_half(&im, &im);
    fp_mul(&(x->im), &im, &(x->im));   
    fp_copy(&(x->re), &re); 
}


// Lexicographic comparison of two field elements. Returns +1 if x > y, -1 if x < y, 0 if x = y
int fp2_cmp(fp2_t* x, fp2_t* y){
    fp2_t a, b;
    fp2_frommont(&a, x);
    fp2_frommont(&b, y);
    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a.re[i] > b.re[i])
            return 1;
        if(a.re[i] < b.re[i])
            return -1;
    }
    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a.im[i] > b.im[i])
            return 1;
        if(a.im[i] < b.im[i])
            return -1;
    }
    return 0;
}


// exponentiation 
// TODO could be improved
void fp2_pow(fp2_t *out,const fp2_t * x,const digit_t *exp,const int size) {

    fp2_t acc;
    digit_t exp_tmp[size];
    digit_t bit;

    memcpy((digit_t*)exp_tmp, (digit_t*)exp, size*RADIX/8);
    memcpy((digit_t*)acc.re, (digit_t*)x->re, NWORDS_FIELD*RADIX/8);
    memcpy((digit_t*)acc.im, (digit_t*)x->im, NWORDS_FIELD*RADIX/8);
    fp2_set_one(out);

    for (int i = 0; i < NWORDS_FIELD*RADIX; i++) {
        bit = exp_tmp[0] & 1;
        mp_shiftr(exp_tmp, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp2_mul(out, out, &acc);
        }
        fp2_sqr(&acc, &acc);
    }

}

void digit_print(char *name,digit_t *a) {
    printf("%s0x = ", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016llx", a[i]);
    printf("\n");
}

void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set_one(&b);
    fp2_mul(&b, &b, &a);
    printf("%s0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016llx", b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016llx", b.im[i]);
    printf("\n");
}
