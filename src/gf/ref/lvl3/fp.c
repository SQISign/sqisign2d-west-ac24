#include "include/fp.h"

// TODO should we modify fp_sqrt to take an out
// variable?
void 
fp_sqrt(fp_t * x)
{
    (void)gf65376_sqrt(x, x);
}

// TODO
//
// bool
// fp_is_square(const fp_t * a)
// {
//     return -1 != gf65376_legendre(a);
// }

// // TODO should we modify fp_inv to take an out
// // variable?
// void
// fp_inv(fp_t * x)
// {
//     (void)gf65376_invert(x, x);
// }

// **************** THIS IS SLOW AND WILL BE REMOVED **************** //
// **************** THIS IS SLOW AND WILL BE REMOVED **************** //
// **************** THIS IS SLOW AND WILL BE REMOVED **************** //

const uint64_t p[NWORDS_FIELD] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x40ffffffffffffff };


static void fp_exp3div4(fp_t* out, const fp_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
    fp_t acc;
    uint64_t p_t[NWORDS_FIELD];
    digit_t bit;

    memcpy(p_t, p, NWORDS_FIELD*RADIX/8);
    // memcpy((digit_t*)acc, (digit_t*)a, NWORDS_FIELD*RADIX/8);
    fp_copy(&acc, a);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    fp_set_one(out);

    for (int i = 0; i < NWORDS_FIELD*RADIX-2; i++) {
        bit = p_t[0] & 1;
        mp_shiftr(p_t, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp_mul(out, out, &acc);
        }
        fp_sqr(&acc, &acc);
    }
}

void fp_inv(fp_t* a)
{ // Modular inversion, out = x^-1*R mod p, where R = 2^(w*nwords), w is the computer wordsize and nwords is the number of words to represent p
  // Input: a=xR in [0, p-1] 
  // Output: out in [0, p-1]. It outputs 0 if the input does not have an inverse  
  // Requirement: Ceiling(Log(p)) < w*nwords
    fp_t t;

    fp_exp3div4(&t, a);
    fp_sqr(&t, &t);
    fp_sqr(&t, &t);
    fp_mul(a, &t, a);    // a^(p-2)
}

bool fp_is_square(const fp_t* a)
{ // Is field element a square?
  // Output: out = 0 (false), 1 (true)
    fp_t t;

    fp_exp3div4(&t, a);
    fp_sqr(&t, &t);
    fp_mul(&t, &t, a);    // a^((p-1)/2)
    return fp_is_equal(&t, &ONE);
}

// ^^^^^^^^^^^^^^^^ THIS IS SLOW AND WILL BE REMOVED ^^^^^^^^^^^^^^^^ //
// **************** THIS IS SLOW AND WILL BE REMOVED **************** //
// **************** THIS IS SLOW AND WILL BE REMOVED **************** //

void
fp_copy(fp_t * out, const fp_t * a)
{
    out->v0 = a->v0;
    out->v1 = a->v1;
    out->v2 = a->v2;
    out->v3 = a->v3;
    out->v4 = a->v4;
    out->v5 = a->v5;

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
