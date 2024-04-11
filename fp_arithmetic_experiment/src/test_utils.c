#include "test_utils.h"
#include <stdio.h>

void fp_random_test(fp_t* a)
{ // Generating a pseudo-random field element in [0, p-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i;
    uint64_t a0, a1;

    for (i = 0; i < 4; i++) {
        a0 = rand();
        a1 = rand();
        a->w[i] = a0 * a1;             // Obtain 256-bit number
        // printf("%llu\n", a->w[i]);
    }
    // printf("\n");
   
    // Clear the top few bits
    int top_bits = 5;
    a->w[3] &= (((uint64_t)(-1) << top_bits) >> top_bits);

}

// https://cplusplus.com/reference/cstdlib/qsort/
// For sorting arrays during the benchmarks
int compare (const void * a, const void * b)
{
  return ( *(long*)a - *(long*)b );
}

void fp_random_test_mont(fp_t* a){
  fp_random_test(a);
  fp_to_mont(a, a);
}
