#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <fp2.h>
#include <inttypes.h>

static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 512;       // Number of iterations per test

// VERY NOT SECURE (testing only)
void fp2_random(fp2_t *a){

	fp_set_small(&a->re, rand());
	fp_set_small(&a->im, rand());
	fp2_neg(a, a);

    // Update seed
	uint8_t tmp[8*NWORDS_FIELD];
	fp_encode(&tmp, &(a->re));
	unsigned seed = (unsigned) tmp[0] | (unsigned)tmp[1] << 8 | (unsigned)tmp[2] << 16 | (unsigned)tmp[3] << 24;
    srand((unsigned) seed);
}

int main(int argc, char* argv[])
{
	if (argc > 1) {
		TEST_LOOPS = atoi(argv[1]);
	}

	fp2_t fp2_0, fp2_1;
	// ------------
	fp2_set_zero(&fp2_0);
	fp2_set_one(&fp2_1);
	// ------------

	int i;
	fp2_t a, b, c, d;
	fp_t e;

	for (i = 0; i < TEST_LOOPS; i++)
	{
		printf("[%3d%%] Testing fp2_t arithmetic", 100 * i / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");
                
		// Random elements of fp
		fp2_random(&a);
		fp2_random(&b);
		fp2_copy(&c, &a);
		fp_add(&c.re, &c.re, &ONE);
		fp2_copy(&d, &b);
		fp_add(&d.re, &d.re, &ONE);

		assert(fp2_is_equal(&a,&b) == 0);		// different values check --> (a != b)
		assert(fp2_is_equal(&c,&c) == 1);		// equal values check --> 1 (c == c)

		// Testing neg
		fp2_set_zero(&b);
		fp2_copy(&c, &a);
		fp2_neg(&a, &a);
		fp2_sub(&c, &b, &c);
		assert(fp2_is_equal(&a,&c) == 1);

		fp2_set_one(&a);	// Now a == 1
		fp2_set_zero(&b);	// Now b == 0

		assert(fp2_is_zero(&a) == 0);
		assert(fp2_is_zero(&b) == 1);

		// testing c - c
		fp2_sub(&d, &c, &c);
		assert(fp2_is_zero(&d) == 1);

		// tetsing c * 0
		fp2_mul(&d, &c, &b);
		assert(fp2_is_zero(&d) == 1);

		// tetsing c * 1 ... recall, in Montgomery domain R mod p plays the role of the 1
		fp2_set_one(&a);
		fp2_mul(&d, &c, &a);
		assert(fp2_is_equal(&d, &c) == 1);

		// fp_set(e, 1);	// Now e == 1
		// fp2_pow(d, e, c);
		// assert(fp2_is_equal(&d,& c) == 1);
		
		// fp_set(e, 0);	// Now e == 0
		// fp2_pow(d, e, c);
		// assert(fp2_is_one(&d) == 1);

		// fp2_set(a, 1);	// Now e == R mod p
		// fp_random(e);
		// fp2_pow(d, e, a);
		// assert(fp2_is_one(&d) == 1);

		// Testing 1/a by computing (1/a) x a
		fp2_random(&a);
		fp2_copy(&b, &a);
		fp2_inv(&a);
		fp2_mul(&c, &a, &b);
		assert(fp2_is_one(&c) == 1);

		fp2_random(&a);
		fp2_sqr(&b, &a);
		assert( fp2_is_square(&b) );

	};

	if(TEST_LOOPS){
		printf("[%2d%%] Tested fp2_t arithmetic:\tNo errors!\n", 100 * i /TEST_LOOPS);
	}
	printf("-- All tests passed.\n");
	return 0;
}
