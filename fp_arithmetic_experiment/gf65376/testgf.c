#define _GNU_SOURCE
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * gmplib is used in tests to verify that the functions compute things
 * correctly.
 */
#include <gmp.h>

/*
 * A custom SHA-3 / SHAKE implementation is used for pseudorandom (but
 * reproducible) generation of test values.
 */
#include "sha3.h"

#include "gf65376.h"

static void check_true(int v, const char *banner) {
  if (!v) {
    fprintf(stderr, "ERR: %s\n", banner);
    exit(EXIT_FAILURE);
  }
}

#define ST(x) ST_(x)
#define ST_(x) #x
#define CHECK(op) check_true((op), "line " ST(__LINE__))

static void check_gf_ops(const uint8_t *va, const uint8_t *vb,
                         const uint8_t *vx) {
  mpz_t za, zb, zc, zd, zq;
  gf65376 a, b, c, d;
  uint8_t vc[48];

  mpz_init(za);
  mpz_init(zb);
  mpz_init(zc);
  mpz_init(zd);
  mpz_init(zq);

  mpz_set_ui(zq, 65);
  mpz_mul_2exp(zq, zq, 376);
  mpz_sub_ui(zq, zq, 1);


  gf65376_decode_reduce(&a, va, 48);
	gf65376_decode_reduce(&b, vb, 48);
	mpz_import(za, 48, -1, 1, -1, 0, va);
	mpz_import(zb, 48, -1, 1, -1, 0, vb);

	gf65376_encode(vc, &a);
	mpz_import(zc, 48, -1, 1, -1, 0, vc);
	CHECK(mpz_cmp(zc, zq) < 0);
	mpz_fdiv_r(zd, za, zq);
	CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_add(&c, &a, &b);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_add(zd, za, zb);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_sub(&c, &a, &b);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_sub(zd, za, zb);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_neg(&c, &a);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_neg(zd, za);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_mul(&c, &a, &b);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_mul(zd, za, zb);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_square(&c, &a);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_mul(zd, za, za);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_mul2(&c, &a);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_mul_ui(zd, za, 2);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  // TODO: not working with gmp, but seems fine
	gf65376_mul2(&c, &a);
	gf65376_half(&d, &c);
	CHECK(gf65376_equals(&a, &d) == 0xFFFFFFFF);

  uint32_t f = (uint32_t)vx[0]
  	| ((uint32_t)vx[1] << 8)
  	| ((uint32_t)vx[2] << 16)
  	| ((uint32_t)vx[3] << 24);
  gf65376_mul_small(&c, &a, f);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_mul_ui(zd, za, f);
  mpz_fdiv_r(zd, zd, zq);
  CHECK(mpz_cmp(zc, zd) == 0);

  gf65376_set_small(&c, f);
  gf65376_encode(vc, &c);
  mpz_import(zc, 48, -1, 1, -1, 0, vc);
  CHECK(mpz_cmp(zc, zq) < 0);
  mpz_set_ui(zd, f);
  CHECK(mpz_cmp(zc, zd) == 0);

  uint8_t tmp[144];
  memcpy(tmp, va, 48);
  tmp[31] &= 0x07;
  mpz_import(zd, 48, -1, 1, -1, 0, tmp);
  if (mpz_cmp(zd, zq) >= 0) {
  	CHECK(gf65376_decode(&c, tmp) == 0);
  	CHECK(gf65376_iszero(&c) == 0xFFFFFFFF);
  } else {
  	CHECK(gf65376_decode(&c, tmp) == 0xFFFFFFFF);
  	gf65376_encode(tmp, &c);
  	mpz_import(zc, 48, -1, 1, -1, 0, tmp);
  	CHECK(mpz_cmp(zc, zd) == 0);
  }

  memcpy(tmp, va, 48);
  memcpy(tmp + 48, vb, 48);
  memcpy(tmp + 96, vx, 48);
  for (size_t k = 0; k <= 144; k ++) {
  	gf65376_decode_reduce(&c, tmp, k);
  	gf65376_encode(vc, &c);
  	mpz_import(zc, 48, -1, 1, -1, 0, vc);
  	CHECK(mpz_cmp(zc, zq) < 0);
  	mpz_import(zd, k, -1, 1, -1, 0, tmp);
  	mpz_fdiv_r(zd, zd, zq);
  	CHECK(mpz_cmp(zc, zd) == 0);
  }

  // if (gf65376_iszero(&b)) {
  // 	CHECK(gf65376_div(&c, &a, &b) == 0);
  // 	CHECK(gf65376_iszero(&c) == 0xFFFFFFFF);
  // 	CHECK(gf65376_invert(&c, &b) == 0);
  // 	CHECK(gf65376_iszero(&c) == 0xFFFFFFFF);
  // } else {
  // 	CHECK(gf65376_div(&c, &a, &b) == 0xFFFFFFFF);
  // 	gf65376_mul(&d, &c, &b);
  // 	CHECK(gf65376_equals(&a, &d));
  // 	CHECK(gf65376_invert(&c, &b) == 0xFFFFFFFF);
  // 	gf65376_mul(&c, &c, &b);
  // 	CHECK(gf65376_equals(&c, &gf65376_ONE) == 0xFFFFFFFF);
  // 	gf65376_encode(vc, &c);
  // 	CHECK(vc[0] == 1);
  // 	for (size_t k = 1; k < 48; k ++) {
  // 		CHECK(vc[k] == 0);
  // 	}
  // }
}

static void test_gf65376(void) {
  uint8_t va[48], vb[48], vx[48];
  gf65376 a, b;

  memset(va, 0, sizeof va);
  memset(vb, 0, sizeof vb);
  memset(vx, 0, sizeof vx);
  check_gf_ops(va, vb, vx);
  gf65376_decode_reduce(&a, va, 48);
  CHECK(gf65376_iszero(&a) == 0xFFFFFFFF);
  // CHECK(gf65376_legendre(&a) == 0);
  // CHECK(gf65376_sqrt(&a, &a) == 0xFFFFFFFF);
  CHECK(gf65376_iszero(&a) == 0xFFFFFFFF);
  printf(".");
  fflush(stdout);

  memset(va, 0xFF, sizeof va);
  memset(vb, 0xFF, sizeof vb);
  memset(vx, 0xFF, sizeof vx);
  check_gf_ops(va, vb, vx);
  gf65376_decode_reduce(&a, va, 48);
  gf65376_decode_reduce(&b, vb, 48);
  CHECK(gf65376_iszero(&a) == 0);
  CHECK(gf65376_iszero(&b) == 0);
  CHECK(gf65376_equals(&a, &b) == 0xFFFFFFFF);
  printf(".");
  fflush(stdout);

  va[47] = 0x40;
  check_gf_ops(va, vb, vx);
  gf65376_decode_reduce(&a, va, 48);
  CHECK(gf65376_iszero(&a) == 0xFFFFFFFF);
  CHECK(gf65376_equals(&a, &b) == 0);
  printf(".");
  fflush(stdout);

  for (int i = 0; i < 10000; i++) {
    shake_context rng;
    uint8_t tmp[4];
    shake_init(&rng, 256);
    tmp[0] = (uint8_t)i;
    tmp[1] = (uint8_t)(i >> 8);
    tmp[2] = (uint8_t)(i >> 16);
    tmp[3] = (uint8_t)(i >> 24);
    shake_inject(&rng, tmp, sizeof tmp);
    shake_flip(&rng);
    shake_extract(&rng, va, sizeof va);
    shake_extract(&rng, vb, sizeof vb);
    shake_extract(&rng, vx, sizeof vx);
    check_gf_ops(va, vb, vx);
    gf65376_decode_reduce(&a, va, 48);
    gf65376_decode_reduce(&b, vb, 48);
    CHECK(gf65376_iszero(&a) == 0);
    CHECK(gf65376_equals(&a, &b) == 0);

    gf65376 s, s2, t;
    uint8_t vt[48];
    gf65376_square(&s, &a);
    gf65376_neg(&s2, &s);
    // CHECK(gf65376_legendre(&s) == 1);
    // CHECK(gf65376_legendre(&s2) == -1);
    CHECK(gf65376_sqrt(&t, &s) == 0xFFFFFFFF);
    gf65376_encode(vt, &t);
    CHECK((vt[0] & 0x01) == 0);
    gf65376_square(&t, &t);
    CHECK(gf65376_equals(&t, &s) == 0xFFFFFFFF);
    CHECK(gf65376_sqrt(&t, &s2) == 0);
    gf65376_encode(vt, &t);
    CHECK((vt[0] & 0x01) == 0);
    gf65376_square(&t, &t);
    CHECK(gf65376_equals(&t, &s) == 0xFFFFFFFF);

    if (i % 256 == 0) {
      printf(".");
      fflush(stdout);
    }
  }
  printf("\n");
}

#if (defined _MSC_VER && defined _M_X64) \
        || (defined __x86_64__ && (defined __GNUC__ || defined __clang__))

/*
 * We use __rdtsc() to access the cycle counter.
 *
 * Ideally we should use __rdpmc(0x40000001) to get the actual cycle
 * counter (independently of frequency scaling), but this requires access
 * to have been enabled (on Linux at least):
 *    echo 2 > /sys/bus/event_source/devices/cpu/rdpmc
 * This command requires root privileges, but once it is done, normal
 * user processes can use __rdpmc() to read the cycle counter (the
 * setting is reset at boot time).
 * (If that file contains 1 instead of 2, then we may read the counter
 * as well, but only if we first enable the performance events, which
 * means more code.)
 *
 * For the kind of things we do, though, there should be no practical
 * difference between TSC and the raw cycle counter as long as frequency
 * scaling does not happen after some warmup. In particular, TurboBoost
 * should be disabled, which is usually done in the BIOS screen (if you
 * are testing this on a VM, you might be out of luck).
 *
 * SMT ("HyperThreading") should also be disabled. This can be done
 * (with root privileges) through the following command:
 *    echo off > /sys/devices/system/cpu/smt/control
 */
#include <immintrin.h>

static inline uint64_t
core_cycles(void)
{
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}

// Make n random-ish field elements (for tests only!).
static void
rand_gfs(gf65376 *aa, size_t n)
{
	shake_context sc;
	uint8_t tmp[48];
	uint64_t z;

	z = core_cycles();
	shake_init(&sc, 256);
	for (int i = 0; i < 8; i ++) {
		tmp[i] = (uint8_t)(z >> (8 * i));
	}
	shake_inject(&sc, tmp, 8);
	shake_flip(&sc);
	while (n -- > 0) {
		shake_extract(&sc, tmp, 48);
		gf65376_decode_reduce(aa, tmp, 48);
		aa ++;
	}
}

static int
cmp_u64(const void *v1, const void *v2)
{
	uint64_t x1 = *(const uint64_t *)v1;
	uint64_t x2 = *(const uint64_t *)v2;
	if (x1 < x2) {
		return -1;
	} else if (x1 == x2) {
		return 0;
	} else {
		return 1;
	}
}

static void
speed_add(void)
{
	gf65376 xx[2], x, y;
	uint64_t tt[20];
	uint8_t tmp[48];

	rand_gfs(xx, 2);
	x = xx[0];
	y = xx[1];
	for (int i = 0; i < 20; i ++) {
		uint64_t begin = core_cycles();
		for (int j = 0; j < 1000; j ++) {
			gf65376_add(&x, &x, &y);
			gf65376_add(&y, &x, &y);
			gf65376_add(&x, &x, &y);
			gf65376_add(&y, &x, &y);
			gf65376_add(&x, &x, &y);
			gf65376_add(&y, &x, &y);
		}
		uint64_t end = core_cycles();
		tt[i] = end - begin;
	}
	gf65376_encode(tmp, &y);
	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
	printf("GF(65*2^376 - 1) add:                %11.2f   (%u)\n",
		(double)tt[15] / 6000.0, tmp[0]);
}

static void
speed_sub(void)
{
	gf65376 xx[2], x, y;
	uint64_t tt[20];
	uint8_t tmp[48];

	rand_gfs(xx, 2);
	x = xx[0];
	y = xx[1];
	for (int i = 0; i < 20; i ++) {
		uint64_t begin = core_cycles();
		for (int j = 0; j < 1000; j ++) {
			gf65376_sub(&x, &x, &y);
			gf65376_sub(&y, &x, &y);
			gf65376_sub(&x, &x, &y);
			gf65376_sub(&y, &x, &y);
			gf65376_sub(&x, &x, &y);
			gf65376_sub(&y, &x, &y);
		}
		uint64_t end = core_cycles();
		tt[i] = end - begin;
	}
	gf65376_encode(tmp, &y);
	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
	printf("GF(65*2^376 - 1) sub:                %11.2f   (%u)\n",
		(double)tt[15] / 6000.0, tmp[0]);
}

static void
speed_mul(void)
{
	gf65376 xx[2], x, y;
	uint64_t tt[20];
	uint8_t tmp[48];

	rand_gfs(xx, 2);
	x = xx[0];
	y = xx[1];
	for (int i = 0; i < 20; i ++) {
		uint64_t begin = core_cycles();
		for (int j = 0; j < 1000; j ++) {
			gf65376_mul(&x, &x, &y);
			gf65376_mul(&y, &x, &y);
			gf65376_mul(&x, &x, &y);
			gf65376_mul(&y, &x, &y);
			gf65376_mul(&x, &x, &y);
			gf65376_mul(&y, &x, &y);
		}
		uint64_t end = core_cycles();
		tt[i] = end - begin;
	}
	gf65376_encode(tmp, &y);
	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
	printf("GF(65*2^376 - 1) mul:                %11.2f   (%u)\n",
		(double)tt[15] / 6000.0, tmp[0]);
}

static void
speed_square(void)
{
	gf65376 x;
	uint64_t tt[20];
	uint8_t tmp[48];

	rand_gfs(&x, 1);
	for (int i = 0; i < 20; i ++) {
		uint64_t begin = core_cycles();
		gf65376_xsquare(&x, &x, 6000);
		uint64_t end = core_cycles();
		tt[i] = end - begin;
	}
	gf65376_encode(tmp, &x);
	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
	printf("GF(65*2^376 - 1) square:             %11.2f   (%u)\n",
		(double)tt[15] / 6000.0, tmp[0]);
}

// static void
// speed_div(void)
// {
// 	gf65376 xx[2], x, y;
// 	uint64_t tt[20];
// 	uint8_t tmp[48];

// 	rand_gfs(xx, 2);
// 	x = xx[0];
// 	y = xx[1];
// 	for (int i = 0; i < 20; i ++) {
// 		uint64_t begin = core_cycles();
// 		for (int j = 0; j < 1000; j ++) {
// 			gf65376_div(&x, &x, &y);
// 			gf65376_div(&y, &x, &y);
// 			gf65376_div(&x, &x, &y);
// 			gf65376_div(&y, &x, &y);
// 			gf65376_div(&x, &x, &y);
// 			gf65376_div(&y, &x, &y);
// 		}
// 		uint64_t end = core_cycles();
// 		tt[i] = end - begin;
// 	}
// 	gf65376_encode(tmp, &y);
// 	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
// 	printf("GF(65*2^376 - 1) div:                %11.2f   (%u)\n",
// 		(double)tt[15] / 6000.0, tmp[0]);
// }

static void
speed_sqrt(void)
{
	gf65376 x;
	uint64_t tt[20];
	uint8_t tmp[48];

	rand_gfs(&x, 1);
	for (int i = 0; i < 20; i ++) {
		uint64_t begin = core_cycles();
		for (int j = 0; j < 6000; j ++) {
			gf65376_sqrt(&x, &x);
		}
		uint64_t end = core_cycles();
		tt[i] = end - begin;
	}
	gf65376_encode(tmp, &x);
	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
	printf("GF(65*2^376 - 1) sqrt:               %11.2f   (%u)\n",
		(double)tt[15] / 6000.0, tmp[0]);
}

// static void
// speed_legendre(void)
// {
// 	gf65376 x;
// 	uint64_t tt[20];
// 	uint8_t tmp[48];

// 	rand_gfs(&x, 1);
// 	for (int i = 0; i < 20; i ++) {
// 		uint64_t begin = core_cycles();
// 		for (int j = 0; j < 6000; j ++) {
// 			int32_t ls = gf65376_legendre(&x);
// 			x.v0 += 2 + ls;
// 		}
// 		uint64_t end = core_cycles();
// 		tt[i] = end - begin;
// 	}
// 	gf65376_encode(tmp, &x);
// 	qsort(tt + 10, 10, sizeof tt[0], &cmp_u64);
// 	printf("GF(65*2^376 - 1) Legendre:           %11.2f   (%u)\n",
// 		(double)tt[15] / 6000.0, tmp[0]);
// }
#endif

int main(void) {
  test_gf65376();
  #if (defined _MSC_VER && defined _M_X64) \
        || (defined __x86_64__ && (defined __GNUC__ || defined __clang__))
  	speed_add();
  	speed_sub();
  	speed_mul();
  	speed_square();
  // 	speed_div();
  	speed_sqrt();
  // 	speed_legendre();
  #endif
  return 0;
}
