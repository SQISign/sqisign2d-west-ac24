#include "include/fp.h"

// TODO should we modify fp_sqrt to take an out
// variable?
void 
fp_sqrt(fp_t * x)
{   
    fp_t tmp = *x;
    (void)gf27500_sqrt(x, &tmp);
}

bool
fp_is_square(const fp_t * a)
{
    return -1 != gf27500_legendre(a);
}

// TODO should we modify fp_inv to take an out
// variable?
void
fp_inv(fp_t * x)
{
    fp_t tmp = *x;
    (void)gf27500_invert(x, &tmp);
}


void
fp_copy(fp_t * out, const fp_t * a)
{
    out->v0 = a->v0;
    out->v1 = a->v1;
    out->v2 = a->v2;
    out->v3 = a->v3;
    out->v4 = a->v4;
    out->v5 = a->v5;
    out->v6 = a->v6;
    out->v7 = a->v7;

}

void fp_set_zero(fp_t * a){
    fp_copy(a, &ZERO);
}

void fp_set_one(fp_t * a){
    fp_copy(a, &ONE);
}

/**********************                                    ***********************/
/********************** The below should probably be moved ***********************/
/**********************                                    ***********************/

void MUL(digit_t* out, const digit_t a, const digit_t b)
{ // Digit multiplication, digit*digit -> 2-digit result 
  // Inputs: a, b in [0, 2^w-1], where w is the computer wordsize 
  // Output: 0 < out < 2^(2w)-1    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t)*4);            // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t)*4);

    albl = al * bl;
    albh = al * bh;
    ahbl = ah * bl;
    ahbh = ah * bh;
    out[0] = albl & mask_low;                 // out00

    res1 = albl >> (sizeof(digit_t)*4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t)*4);
    out[0] ^= temp << (sizeof(digit_t)*4);    // out01   

    res1 = ahbl >> (sizeof(digit_t)*4);
    res2 = albh >> (sizeof(digit_t)*4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    out[1] = temp & mask_low;                 // out10 
    carry = temp & mask_high;
    out[1] ^= (ahbh & mask_high) + carry;     // out11
}


void mp_add(digit_t* c, const digit_t* a, const digit_t* b, const unsigned int nwords)
{ // Multiprecision addition
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(c[i], carry, a[i], b[i], carry);
    }
}

digit_t mp_shiftr(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision right shift
    digit_t bit_out = x[0] & 1;

    for (unsigned int i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], shift, x[i], RADIX);
    }
    x[nwords-1] >>= shift;
    return bit_out;
}

void mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision left shift of at most 64

    for (int i = nwords-1; i > 0; i--) {
        SHIFTL(x[i], x[i-1], shift, x[i], RADIX);
    }
    x[0] <<= shift;
}

void multiple_mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords) {
    int t = shift;
    while (t>60) {
        mp_shiftl(x,60,nwords);
        t = t-60;
    }
    mp_shiftl(x,t,nwords);
}
