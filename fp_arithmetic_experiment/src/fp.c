#include "fp.h"
#include "utils.h"

// Computation of (x + y) mod p
void fp_add(fp_t *out, const fp_t *x, const fp_t *y) {
  uint64_t d0, d1, d2, d3, f;
  uint8_t cc;

  // Raw addition of x and y
  cc = addcarry_u64(&d0, x->w[0], y->w[0], 0);
  cc = addcarry_u64(&d1, x->w[1], y->w[1], cc);
  cc = addcarry_u64(&d2, x->w[2], y->w[2], cc);
  cc = addcarry_u64(&d3, x->w[3], y->w[3], cc); // cc here is unused

  // Inputs were up to 2^251 - 1; sum can be up to 2^252 - 2.
  // We subtract q if the value is not lower than 2^251. Subtraction
  // of q is done by adding -q (modulo 2^256) (because -q has two
  // limbs equal to zero, which saves a few operations).
  f = d3 >> 59;
  cc = addcarry_u64(&d0, d0, f, 0);
  cc = addcarry_u64(&d1, d1, 0, cc);
  cc = addcarry_u64(&d2, d2, 0, cc);
  cc = addcarry_u64(&d3, d3, ((uint64_t)0xFB << 56) & (-f), cc); // cc here is unused

  // If the value is still not lower than 2^251 then we subtract q again.
  f = d3 >> 59;
  cc = addcarry_u64(&d0, d0, f, 0);
  cc = addcarry_u64(&d1, d1, 0, cc);
  cc = addcarry_u64(&d2, d2, 0, cc);
  cc = addcarry_u64(&d3, d3, ((uint64_t)0xFB << 56) & (-f), cc); // cc here is unused

  out->w[0] = d0;
  out->w[1] = d1;
  out->w[2] = d2;
  out->w[3] = d3;
}

void fp_sub(fp_t *out, const fp_t *x, const fp_t *y){
  uint64_t d0, d1, d2, d3, f, m;
  uint8_t cc;

  // Raw subtraction of x and y
  cc = subborrow_u64(&d0, x->w[0], y->w[0], 0);
  cc = subborrow_u64(&d1, x->w[1], y->w[1], cc);
  cc = subborrow_u64(&d2, x->w[2], y->w[2], cc);
  cc = subborrow_u64(&d3, x->w[3], y->w[3], cc);

  // If the result is negative, we add 2*q (by subtracting -2*q).
  cc = subborrow_u64(&m, 0, 0, cc); // cc is unused here
  cc = subborrow_u64(&d0, d0, m & 2, 0);
  cc = subborrow_u64(&d1, d1, 0, cc);
  cc = subborrow_u64(&d2, d2, 0, cc);
  cc = subborrow_u64(&d3, d3, ((uint64_t)0xF6 << 56) & m, cc); // cc here is unused

  // If the value is not lower than 2^251 then we subtract q (i.e.
  // we add -q).
  f = d3 >> 59;
  cc = addcarry_u64(&d0, d0, f, 0);
  cc = addcarry_u64(&d1, d1, 0, cc);
  cc = addcarry_u64(&d2, d2, 0, cc);
  cc = addcarry_u64(&d3, d3, ((uint64_t)0xFB << 56) & (-f), cc); // cc here is unused

  out->w[0] = d0;
  out->w[1] = d1;
  out->w[2] = d2;
  out->w[3] = d3;
}

void fp_neg(fp_t *out, const fp_t *x){
  uint64_t d0, d1, d2, d3, f;
  uint8_t cc;

  // 2*q - x
  cc = subborrow_u64(&d0, MODULUS_X2.w[0], x->w[0], 0);
  cc = subborrow_u64(&d1, MODULUS_X2.w[1], x->w[1], cc);
  cc = subborrow_u64(&d2, MODULUS_X2.w[2], x->w[2], cc);
  cc = subborrow_u64(&d3, MODULUS_X2.w[3], x->w[3], cc); // cc here is unused

  // If the value is not lower than 2^251 then we subtract q (i.e.
  // we add -q).
  f = d3 >> 59;
  cc = addcarry_u64(&d0, d0, f, 0);
  cc = addcarry_u64(&d1, d1, 0, cc);
  cc = addcarry_u64(&d2, d2, 0, cc);
  cc = addcarry_u64(&d3, d3, ((uint64_t)0xFB << 56) & (-f), cc); // cc here is unused

  out->w[0] = d0;
  out->w[1] = d1;
  out->w[2] = d2;
  out->w[3] = d3;
}

// Computation of a * b mod p
// (Assumes a, b in Montgomery representation)
void fp_mul(fp_t *out, const fp_t *a, const fp_t *b) {
  uint64_t a0, a1, a2, a3;
  uint64_t b0, b1, b2, b3;
  uint64_t d0, d1, d2, d3;
  uint64_t e0, e1, e2, e3, e4, e5, e6, e7;
  uint64_t f0, f1, f2, f3;
  uint64_t g0, g1, g2, g3, g4, g5, g6, g7;
  uint64_t lo, hi, lo2, hi2, tt;
  uint8_t cc;

  // Extract out words
  a0 = a->w[0];
  a1 = a->w[1];
  a2 = a->w[2];
  a3 = a->w[3];

  b0 = b->w[0];
  b1 = b->w[1];
  b2 = b->w[2];
  b3 = b->w[3];

  // Compute the product 502 bits
  e0 = umull(&e1, a0, b0);
  e2 = umull(&e3, a1, b1);
  e4 = umull(&e5, a2, b2);
  e6 = umull(&e7, a3, b3);

  lo = umull(&hi, a0, b1);
  cc = addcarry_u64(&e1, e1, lo, 0);
  cc = addcarry_u64(&e2, e2, hi, cc);
  lo = umull(&hi, a0, b3);
  cc = addcarry_u64(&e3, e3, lo, cc);
  cc = addcarry_u64(&e4, e4, hi, cc);
  lo = umull(&hi, a2, b3);
  cc = addcarry_u64(&e5, e5, lo, cc);
  cc = addcarry_u64(&e6, e6, hi, cc);
  cc = addcarry_u64(&e7, e7, 0, cc); // cc here is unused

  lo = umull(&hi, a1, b0);
  cc = addcarry_u64(&e1, e1, lo, 0);
  cc = addcarry_u64(&e2, e2, hi, cc);
  lo = umull(&hi, a3, b0);
  cc = addcarry_u64(&e3, e3, lo, cc);
  cc = addcarry_u64(&e4, e4, hi, cc);
  lo = umull(&hi, a3, b2);
  cc = addcarry_u64(&e5, e5, lo, cc);
  cc = addcarry_u64(&e6, e6, hi, cc);
  cc = addcarry_u64(&e7, e7, 0, cc); // cc here is unused

  lo = umull(&hi, a0, b2);
  cc = addcarry_u64(&e2, e2, lo, 0);
  cc = addcarry_u64(&e3, e3, hi, cc);
  lo = umull(&hi, a1, b3);
  cc = addcarry_u64(&e4, e4, lo, cc);
  cc = addcarry_u64(&e5, e5, hi, cc);
  cc = addcarry_u64(&e6, e6, 0, cc);
  cc = addcarry_u64(&e7, e7, 0, cc); // cc here is unused

  lo = umull(&hi, a2, b0);
  cc = addcarry_u64(&e2, e2, lo, 0);
  cc = addcarry_u64(&e3, e3, hi, cc);
  lo = umull(&hi, a3, b1);
  cc = addcarry_u64(&e4, e4, lo, cc);
  cc = addcarry_u64(&e5, e5, hi, cc);
  cc = addcarry_u64(&e6, e6, 0, cc);
  cc = addcarry_u64(&e7, e7, 0, cc); // cc here is unused

  lo = umull(&hi, a1, b2);
  lo2 = umull(&hi2, a2, b1);
  cc = addcarry_u64(&lo, lo, lo2, 0);
  tt = addcarry_u64(&hi, hi, hi2, cc);
  cc = addcarry_u64(&e3, e3, lo, 0);
  cc = addcarry_u64(&e4, e4, hi, cc);
  cc = addcarry_u64(&e5, e5, (uint64_t)tt, cc);
  cc = addcarry_u64(&e6, e6, 0, cc);
  cc = addcarry_u64(&e7, e7, 0, cc); // cc here is unused

  // Montgomery reduction.
  //
  // The low part is lo(e) = e0..e3 (256 bits).
  // Let m = -1/q mod 2^256; the Montgomery reduction is adding
  // (lo(e)*m mod 2^256)*q to the high part g = e4..e7 (246 bits).
  //
  // We have m = 5*2^248 + 1; the product f = lo(e)*m mod 2^256 is
  // obtained by simply added e0*5 (mod 2^8) to the high byte of e3.
  f0 = e0;
  f1 = e1;
  f2 = e2;
  f3 = e3 + ((e0 * 5) << 56);

  // Compute f * q
  g3 = umull(&hi, f0, (uint64_t)5 << 56);
  g4 = umull_add(&hi, f1, (uint64_t)5 << 56, hi);
  g5 = umull_add(&hi, f2, (uint64_t)5 << 56, hi);
  g6 = umull_add(&g7, f3, (uint64_t)5 << 56, hi);

  cc = subborrow_u64(&g0, 0, f0, 0);
  cc = subborrow_u64(&g1, 0, f1, cc);
  cc = subborrow_u64(&g2, 0, f2, cc);
  cc = subborrow_u64(&g3, g3, f3, cc);
  cc = subborrow_u64(&g4, g4, 0, cc);
  cc = subborrow_u64(&g5, g5, 0, cc);
  cc = subborrow_u64(&g6, g6, 0, cc);
  cc = subborrow_u64(&g7, g7, 0, cc); // cc here is unused

  // We add g = f*q to e0..e7. Since e0..e7 < 2^502, and f < 2^256,
  // we know that the result is not less than
  // 2^502 + 2^256*5*2^248 < 6*2^504; it is also a multiple of
  // 2^256. After dividing by 2^256, we get a value which is
  // less than 6*2^248, i.e. already in our proper range.

  uint64_t tmp; // this is never used...
  cc = addcarry_u64(&tmp, e0, g0, 0);
  cc = addcarry_u64(&tmp, e1, g1, cc);
  cc = addcarry_u64(&tmp, e2, g2, cc);
  cc = addcarry_u64(&tmp, e3, g3, cc);
  cc = addcarry_u64(&d0, e4, g4, cc);
  cc = addcarry_u64(&d1, e5, g5, cc);
  cc = addcarry_u64(&d2, e6, g6, cc);
  cc = addcarry_u64(&d3, e7, g7, cc);

  out->w[0] = d0;
  out->w[1] = d1;
  out->w[2] = d2;
  out->w[3] = d3;
}

void fp_sqr(fp_t *out, const fp_t *a){
  uint64_t a0, a1, a2, a3;
  uint64_t d0, d1, d2, d3;
  uint64_t e0, e1, e2, e3, e4, e5, e6, e7;
  uint64_t f0, f1, f2, f3;
  uint64_t g0, g1, g2, g3, g4, g5, g6, g7;
  uint64_t lo, hi;
  uint8_t cc;

  // Extract out words
  a0 = a->w[0];
  a1 = a->w[1];
  a2 = a->w[2];
  a3 = a->w[3];

  // 1. Non-square products. Max intermediate value:
  //   a0*a1            * 2^64
  //   a0*a2            * 2^128
  //   (a0*a3 + a1*a2)  * 2^192
  //   a1*a3            * 2^256
  //   a2*a3            * 2^320
  // for a total which is stlightly below 2^448, which means that
  // the value fits on e1..e6 (no possible carry into e7).
  e1 = umull(&e2, a0, a1);
  e3 = umull(&e4, a0, a3);
  e5 = umull(&e6, a2, a3);

  lo = umull(&hi, a0, a2);
  cc = addcarry_u64(&e2, e2, lo, 0);
  cc = addcarry_u64(&e3, e3, hi, cc);

  lo = umull(&hi, a1, a3);
  cc = addcarry_u64(&e4, e4, lo, cc);
  cc = addcarry_u64(&e5, e5, hi, cc);
  cc = addcarry_u64(&e6, e6, 0, cc);

  lo = umull(&hi, a1, a2);
  cc = addcarry_u64(&e3, e3, lo, 0);
  cc = addcarry_u64(&e4, e4, hi, cc);
  cc = addcarry_u64(&e5, e5, 0, cc);
  cc = addcarry_u64(&e6, e6, 0, cc); // cc is unused

  // Double the intermediate value then add squares
  e7 = e6 >> 63;
  e6 = (e6 << 1) | (e5 >> 63);
  e5 = (e5 << 1) | (e4 >> 63);
  e4 = (e4 << 1) | (e3 >> 63);
  e3 = (e3 << 1) | (e2 >> 63);
  e2 = (e2 << 1) | (e1 >> 63);
  e1 = e1 << 1;

  e0 = umull(&hi, a0, a0);
  cc = addcarry_u64(&e1, e1, hi, 0);
  lo = umull(&hi, a1, a1);
  cc = addcarry_u64(&e2, e2, lo, cc);
  cc = addcarry_u64(&e3, e3, hi, cc);
  lo = umull(&hi, a2, a2);
  cc = addcarry_u64(&e4, e4, lo, cc);
  cc = addcarry_u64(&e5, e5, hi, cc);
  lo = umull(&hi, a3, a3);
  cc = addcarry_u64(&e6, e6, lo, cc);
  cc = addcarry_u64(&e7, e7, hi, cc); // cc is unused

  // Reduction, see fp_mul for details
  f0 = e0;
  f1 = e1;
  f2 = e2;
  f3 = e3 + ((e0 * 5) << 56);

  // Compute f * q
  g3 = umull(&hi, f0, (uint64_t)5 << 56);
  g4 = umull_add(&hi, f1, (uint64_t)5 << 56, hi);
  g5 = umull_add(&hi, f2, (uint64_t)5 << 56, hi);
  g6 = umull_add(&g7, f3, (uint64_t)5 << 56, hi);

  cc = subborrow_u64(&g0, 0, f0, 0);
  cc = subborrow_u64(&g1, 0, f1, cc);
  cc = subborrow_u64(&g2, 0, f2, cc);
  cc = subborrow_u64(&g3, g3, f3, cc);
  cc = subborrow_u64(&g4, g4, 0, cc);
  cc = subborrow_u64(&g5, g5, 0, cc);
  cc = subborrow_u64(&g6, g6, 0, cc);
  cc = subborrow_u64(&g7, g7, 0, cc); // cc here is unused

  uint64_t tmp; // this is never used...
  cc = addcarry_u64(&tmp, e0, g0, 0);
  cc = addcarry_u64(&tmp, e1, g1, cc);
  cc = addcarry_u64(&tmp, e2, g2, cc);
  cc = addcarry_u64(&tmp, e3, g3, cc);
  cc = addcarry_u64(&d0, e4, g4, cc);
  cc = addcarry_u64(&d1, e5, g5, cc);
  cc = addcarry_u64(&d2, e6, g6, cc);
  cc = addcarry_u64(&d3, e7, g7, cc); // cc here is unused

  out->w[0] = d0;
  out->w[1] = d1;
  out->w[2] = d2;
  out->w[3] = d3;
}

void fp_n_sqr(fp_t * out, const fp_t * x, int n) {
  int i;
  
  fp_sqr(out, x);
  for (i=1; i < n; i++){
    fp_sqr(out, out);
  }
}

void fp_from_mont(fp_t *out, const fp_t *x) {
  // Let m = -1/q mod 2^256 = 5*2^248 + 1.
  // For input x, we compute f = x*m mod 2^256, then
  // h = x + f*q, which is a multiple of 2^256. The output
  // is then h/2^256. Note that if x < 2^256, then:
  //   h <= 2^256 - 1 + (2^256 - 1)*q
  //   h <= q*2^256 + 2^256 - q - 1
  // Since h = 0 mod 2^256 and, this implies that h <= q*2^256.
  // The output h/2^256 is thus in the 0 to q range (inclusive).

  uint64_t x0, x1, x2, x3;
  uint64_t f0, f1, f2, f3;
  uint64_t d0, d1, d2, d3;
  uint64_t g0, g1, g2, g3, g4, g5, g6, g7;
  uint64_t hi, t, w;
  uint8_t cc;

  // Grab the input a
  x0 = x->w[0];
  x1 = x->w[1];
  x2 = x->w[2];
  x3 = x->w[3];

  // f = x*(-1/p) mod 2**256
  f0 = x0;
  f1 = x1;
  f2 = x2;
  f3 = x3 + ((x0 * 5) << 56);

  // g = f*p
  g3 = umull(&hi, f0, (uint64_t)5 << 56);
  g4 = umull_add(&hi, f1, (uint64_t)5 << 56, hi);
  g5 = umull_add(&hi, f2, (uint64_t)5 << 56, hi);
  g6 = umull_add(&g7, f3, (uint64_t)5 << 56, hi);

  cc = subborrow_u64(&g0, 0, f0, 0);
  cc = subborrow_u64(&g1, 0, f1, cc);
  cc = subborrow_u64(&g2, 0, f2, cc);
  cc = subborrow_u64(&g3, g3, f3, cc);
  cc = subborrow_u64(&g4, g4, 0, cc);
  cc = subborrow_u64(&g5, g5, 0, cc);
  cc = subborrow_u64(&g6, g6, 0, cc);
  cc = subborrow_u64(&g7, g7, 0, cc); // cc here is unused

  // h = x + f*q
  // We drop the lower 256 bits
  uint64_t tmp; // this is never used...
  cc = addcarry_u64(&tmp, g0, x0, 0);
  cc = addcarry_u64(&tmp, g1, x1, cc);
  cc = addcarry_u64(&tmp, g2, x2, cc);
  cc = addcarry_u64(&tmp, g3, x3, cc);
  cc = addcarry_u64(&d0, g4, 0, cc);
  cc = addcarry_u64(&d1, g5, 0, cc);
  cc = addcarry_u64(&d2, g6, 0, cc);
  cc = addcarry_u64(&d3, g7, 0, cc); // cc here is unused

  // h is in [0..q]
  // To normalise, we set h to zero if h = p
  t = d0 & d1 & d2 & (d3 ^ ~(uint64_t)(MODULUS.w[3]));
  cc = addcarry_u64(&tmp, t, 1, 0); // TODO we don't use the value of tmp here
  cc = subborrow_u64(&w, 0, 0, cc); //      we don't use the value of cc here
  w = ~(uint64_t)(w);

  out->w[0] = d0 & w;
  out->w[1] = d1 & w;
  out->w[2] = d2 & w;
  out->w[3] = d3 & w;
}

void fp_to_mont(fp_t *out, const fp_t *a) {
  fp_mul(out, a, (fp_t *)&R2);
}

void fp_copy(fp_t *out, const fp_t *a){
  out->w[0] = a->w[0];
  out->w[1] = a->w[1];
  out->w[2] = a->w[2];
  out->w[3] = a->w[3];
}

void fp_set_zero(fp_t *out){
  out->w[0] = ZERO.w[0];
  out->w[1] = ZERO.w[1];
  out->w[2] = ZERO.w[2];
  out->w[3] = ZERO.w[3];
}

void fp_set_one(fp_t *out){
  out->w[0] = ONE.w[0];
  out->w[1] = ONE.w[1];
  out->w[2] = ONE.w[2];
  out->w[3] = ONE.w[3];
}

uint32_t fp_is_zero(fp_t *x) {
  uint64_t x0, x1, x2, x3;
  uint64_t t0, t1, r;

  x0 = x->w[0];
  x1 = x->w[1];
  x2 = x->w[2];
  x3 = x->w[3];

  t0 = x0 | x1 | x2 | x3;
  t1 = (x0 ^ MODULUS.w[0]) | (x1 ^ MODULUS.w[1]) | (x2 ^ MODULUS.w[2]) | (x3 ^ MODULUS.w[3]);
  r = (t0 | (-t0)) & (t1 | (-t1));
  return ((uint32_t)(r >> 63) - 1);
}

uint32_t fp_equal(fp_t *x, fp_t *y) {
  fp_t check;
  fp_sub(&check, x, y);
  return fp_is_zero(&check);
}

static void fp_exp3div4(fp_t* out, const fp_t* a)
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
  fp_t z3, z4, t3, t6, tmp;

  // Compute a**3 and a**4
  fp_sqr(&z4, a);
  fp_mul(&z3, a, &z4);
  fp_sqr(&z4, &z4);

  // Compute a**(2^3 - 1) = a**7
  fp_mul(&t3, &z3, &z4);

  // Compute a**(2^6 - 1)
  fp_n_sqr(&t6, &t3, 3);
  fp_mul(&t6, &t6, &t3);

  // Compute a**(2^12 - 1)
  fp_n_sqr(out, &t6, 6);
  fp_mul(out, out, &t6);

  // Compute a**(2^15 - 1)
  fp_n_sqr(out, out, 3);
  fp_mul(out, out, &t3);

  // Compute a**(2^30 - 1)
  fp_n_sqr(&tmp, out, 15);
  fp_mul(out, out, &tmp);

  // Compute a**(2^60 - 1)
  fp_n_sqr(&tmp, out, 30);
  fp_mul(out, out, &tmp);

  // Compute a**(2^120 - 1)
  fp_n_sqr(&tmp, out, 60);
  fp_mul(out, out, &tmp);

  // Compute a**(2^240 - 1)
  fp_n_sqr(&tmp, out, 120);
  fp_mul(out, out, &tmp);

  // Compute a**(2^246 - 1)
  fp_n_sqr(out, out, 6);
  fp_mul(out, out, &t6);

  // Compute a**(5*(2^246 - 1))
  fp_sqr(&tmp, out);
  fp_sqr(&tmp, &tmp);
  fp_mul(out, &tmp, out);

  // Compute a**(5*(2^246 - 1) + 4) 
  fp_mul(out, out, &z4);
}

/// Assumes the input is a square
void fp_sqrt(fp_t* a)
{ // Square root computation, out = a^((p+1)/4) mod p
  // Uses that p+1 = 5*2**248 so we can compute the
  // sqrt with only 248 squares and one multiplication
    fp_t t;

    // Compute a^5
    fp_copy(&t, a);
    fp_sqr(a, a);
    fp_sqr(a, a);
    fp_mul(a, a, &t);

    // Compute (a^(5)) ^ 2^246
    for (int i = 0; i < 246; i++){
        fp_sqr(a, a);
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
    return fp_equal(&t, (fp_t*)&ONE);
}
