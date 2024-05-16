/*
 * A custom SHA-3 / SHAKE implementation is used for pseudorandom (but
 * reproducible) generation of test values.
 */
#include "sha3.h"
#include "test_utils.h"

// Make n random-ish field elements (for tests only!).
void
fp_random_test(fp_t *a)
{
	shake_context sc;
	uint8_t tmp[32];
	uint64_t z;

	z = core_cycles();
	shake_init(&sc, 256);
	for (int i = 0; i < 8; i ++) {
		tmp[i] = (uint8_t)(z >> (8 * i));
	}
	shake_inject(&sc, tmp, 8);
	shake_flip(&sc);
	shake_extract(&sc, tmp, 32);
	fp_decode_reduce(a, tmp, 32);
}

void
fp2_random_test(fp2_t *a)
{
	fp_random_test(&(a->re));
	fp_random_test(&(a->im));
}

int
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

#if (defined _MSC_VER && defined _M_X64) \
        || (defined __x86_64__ && (defined __GNUC__ || defined __clang__))

#include <immintrin.h>

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
uint64_t
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
#endif
