#include <fp.h>

extern uint64_t p[NWORDS_FIELD];

/*
 * If ctl == 0x00000000, then *a and *b are unchanged.
 * If ctl == 0xFFFFFFFF, then the contents of *a and *b are swapped.
 * ctl MUST be either 0x00000000 or 0xFFFFFFFF.
 */
void
fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
{
	uint64_t cw = (uint64_t)*(int32_t *)&ctl;
	uint64_t t;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++) {
        t = cw & ((*a)[i] ^ (*b)[i]); (*a)[i] ^= t; (*b)[i] ^= t;
    }

}

void
fp_set_zero(fp_t* a)
{
    for (unsigned int i = 0; i < NWORDS_FIELD; i++) {
        (*a)[i] = 0;
    }
}

void
fp_set_one(fp_t* a)
{
    fp_mont_setone(a);
}

void 
fp_set_small(fp_t* x, const digit_t val)
{ // Set field element x = val, where val has wordsize
  // automatically converted to Montgomery form
  (*x)[0] = val;
  for (unsigned int i = 1; i < NWORDS_FIELD; i++) {
    (*x)[i] = 0;
  }
  fp_tomont(x, x);
}

bool fp_is_equal(const fp_t* a, const fp_t* b)
{ // Compare two field elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= (*a)[i] ^ (*b)[i];

    return (bool)is_digit_zero_ct(r);
}

bool fp_is_zero(const fp_t* a)
{ // Is a field element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= (*a)[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void fp_copy(fp_t* out, const fp_t* a)
{
    memcpy(*out, *a, NWORDS_FIELD*RADIX/8);
}

void fp_neg(fp_t* out, const fp_t* a)
{ // Modular negation, out = -a mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC((*out)[i], borrow, p[i], (*a)[i], borrow);
    }
    fp_sub(out, out, (fp_t*)&p);
}

void fp_half(fp_t* out, const fp_t* a){
    fp_t inv_two;
    fp_set_small(&inv_two, 2);
    fp_inv(&inv_two);
    fp_mul(out, a, &inv_two);
}

// Little-endian encoding of a 64-bit integer.
static inline void
enc64le(void *dst, uint64_t x)
{
	uint8_t *buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
	buf[4] = (uint8_t)(x >> 32);
	buf[5] = (uint8_t)(x >> 40);
	buf[6] = (uint8_t)(x >> 48);
	buf[7] = (uint8_t)(x >> 56);
}

// Little-endian decoding of a 64-bit integer.
static inline uint64_t
dec64le(const void *src)
{
	const uint8_t *buf = src;
	return (uint64_t)buf[0]
		| ((uint64_t)buf[1] << 8)
		| ((uint64_t)buf[2] << 16)
		| ((uint64_t)buf[3] << 24)
		| ((uint64_t)buf[4] << 32)
		| ((uint64_t)buf[5] << 40)
		| ((uint64_t)buf[6] << 48)
		| ((uint64_t)buf[7] << 56);
}

// Encode elements to bytes
void
fp_encode(void *dst, const fp_t *a)
{
	uint8_t *buf = dst;
	fp_t x;

	fp_frommont(&x, a);
    for (int i = 0; i < NWORDS_FIELD; i++) {
        enc64le(buf + i * 8, x[i]);
    }
}

// TODO: handle lengths
void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
	const uint8_t *buf = src;

    for (int i = 0; i < NWORDS_FIELD; i++) {
        (*d)[i] = dec64le(buf + i * 8);
    }

    fp_tomont(d, d);
}

void
fp_decode(fp_t *d, const void *src)
{
	const uint8_t *buf = src;

    for (int i = 0; i < NWORDS_FIELD; i++) {
        (*d)[i] = dec64le(buf + i * 8);
    }

    fp_tomont(d, d);
}


static void fp_exp3div4(fp_t* out, const fp_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
    fp_t p_t, acc;
    digit_t bit;

    memcpy((digit_t*)p_t, (digit_t*)p, NWORDS_FIELD*RADIX/8);
    memcpy((digit_t*)acc, (digit_t*)*a, NWORDS_FIELD*RADIX/8);
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
    fp_t t, one;

    fp_exp3div4(&t, a);
    fp_sqr(&t, &t);
    fp_mul(&t, &t, a);    // a^((p-1)/2)
    fp_set_one(&one);

    return fp_is_equal(&t, &one);
}

void fp_sqrt(fp_t* a)
{ // Square root computation, out = a^((p+1)/4) mod p
    fp_t t;

    fp_exp3div4(&t, a);
    fp_mul(a, &t, a);    // a^((p+1)/4)
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
