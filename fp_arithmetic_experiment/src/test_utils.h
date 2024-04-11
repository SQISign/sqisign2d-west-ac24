#include "fp.h"
#include <stdlib.h>

#define PASSED    0
#define FAILED    1

// Generating a pseudo-random field element in [0, p-1] 
void fp_random_test(fp_t* a);

// Comparison of u64 for qsort
int cmp_u64(const void *v1, const void *v2);

uint64_t core_cycles(void);
