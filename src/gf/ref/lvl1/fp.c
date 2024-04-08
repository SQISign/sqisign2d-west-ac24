#include "include/fp.h"

const uint64_t p[NWORDS_FIELD] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x4ffffffffffffff };

void fp_set(digit_t* x, const digit_t val)
{ // Set field element x = val, where val has wordsize

    x[0] = val;
    for (unsigned int i = 1; i < NWORDS_FIELD; i++) {
        x[i] = 0;
    }
}

bool fp_is_equal(const digit_t* a, const digit_t* b)
{ // Compare two field elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ b[i];

    return (bool)is_digit_zero_ct(r);
}

bool fp_is_zero(const digit_t* a)
{ // Is a field element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void fp_copy(digit_t* out, const digit_t* a)
{
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}

void fp_neg(digit_t* out, const digit_t* a)
{ // Modular negation, out = -a mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(out[i], borrow, ((digit_t*)p)[i], a[i], borrow);
    }
    fp_sub(out, out, (digit_t*)p);
}

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

void fp_n_sqr(digit_t* out, digit_t* in, int n)
{ // Repeated squaring of an element n times
    fp_copy(out, in);
    for (int i = 0; i < n; i++){
        fp_sqr(out, out);
    }
}

static void fp_exp3div4(digit_t* out, const digit_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
  //


  // We optimise this by using the shape of the prime
  // to avoid almost all multiplications:
  //
  // We write:
  //     (p - 3) / 4 = (5*2^248 - 4) / 4
  //                 = 5*2^246 - 1
  //                 = 5*(2^246 - 1) + 4

  // Then we first compute:
  //     a246 = a**(2^246 - 1)

  // Then from this we get the desired result as:
  //     a**((p-3)/4) = a246**5 * a**4

  // We can compute this with 12 multiplications and
  // 247 squares, so we expect this to be about 1/2 the cost
  // of the naive method using square and multiply
  int i;
  fp_t z3, z4, t3, t6, tmp;

  // Compute a**3 and a**4
  fp_sqr(z4, a);
  fp_mul(z3, a, z4);
  fp_sqr(z4, z4);

  // Compute a**(2^3 - 1) = a**7
  fp_mul(t3, z3, z4);

  // Compute a**(2^6 - 1)
  fp_n_sqr(t6, t3, 3);
  fp_mul(t6, t6, t3);

  // Compute a**(2^12 - 1)
  fp_n_sqr(out, t6, 6);
  fp_mul(out, out, t6);

  // Compute a**(2^15 - 1)
  fp_n_sqr(out, out, 3);
  fp_mul(out, out, t3);

  // Compute a**(2^30 - 1)
  fp_n_sqr(tmp, out, 15);
  fp_mul(out, out, tmp);

  // Compute a**(2^60 - 1)
  fp_n_sqr(tmp, out, 30);
  fp_mul(out, out, tmp);

  // Compute a**(2^120 - 1)
  fp_n_sqr(tmp, out, 60);
  fp_mul(out, out, tmp);

  // Compute a**(2^240 - 1)
  fp_n_sqr(tmp, out, 120);
  fp_mul(out, out, tmp);

  // Compute a**(2^246 - 1)
  fp_n_sqr(out, out, 6);
  fp_mul(out, out, t6);

  // Compute a**(5*(2^246 - 1))
  fp_sqr(tmp, out);
  fp_sqr(tmp, tmp);
  fp_mul(out, tmp, out);

  // Compute a**(5*(2^246 - 1) + 4) 
  fp_mul(out, out, z4);
}


void fp_inv(digit_t* a)
{ // Modular inversion, out = x^-1*R mod p, where R = 2^(w*nwords), w is the computer wordsize and nwords is the number of words to represent p
  // Input: a=xR in [0, p-1] 
  // Output: out in [0, p-1]. It outputs 0 if the input does not have an inverse  
  // Requirement: Ceiling(Log(p)) < w*nwords
    fp_t t;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_sqr(t, t);
    fp_mul(a, t, a);    // a^(p-2)
}

bool fp_is_square(const digit_t* a)
{ // Is field element a square?
  // Output: out = 0 (false), 1 (true)
    fp_t t, one;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_mul(t, t, a);    // a^((p-1)/2)
    fp_frommont(t, t);
    fp_set(one, 1);

    return fp_is_equal(t, one);
}

// Square root computation, out = a^((p+1)/4) mod p
// Uses that p+1 = 5*2**248 so we can compute the
// sqrt with only 248 squares and one multiplication
void fp_sqrt(digit_t* a)
{
    int i;
    fp_t t;

    // Compute a^5
    fp_copy(t, a);
    fp_sqr(a, a);
    fp_sqr(a, a);
    fp_mul(a, a, t);

    // Compute (a^(5)) ^ 2^246
    for (i = 0; i < 246; i++){
        fp_sqr(a, a);
    }
}
