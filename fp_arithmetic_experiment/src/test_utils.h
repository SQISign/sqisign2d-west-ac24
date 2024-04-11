#ifndef TEST_EXTRAS_H
#define TEST_EXTRAS_H

#include <time.h>
#include <stdlib.h>
#include "fp.h"

#define PASSED    0
#define FAILED    1

int compare (const void * a, const void * b);

// Generating a pseudo-random field element in [0, p-1] 
void fp_random_test(fp_t* a);

// Generating a pseudo-random field element in [0, p-1]
// then transfer into Montgomery representation
void fp_random_test_mont(fp_t * a);

#endif
