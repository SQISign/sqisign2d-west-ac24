#include "fp.h"
#include "fp2.h"
#include <stdlib.h>

#define PASSED    0
#define FAILED    1

void fp_random_test(fp_t* a);
void fp2_random_test(fp2_t* a);

// Comparison of u64 for qsort
int cmp_u64(const void *v1, const void *v2);

uint64_t core_cycles(void);
