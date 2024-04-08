#include "test_extras.h"
#include <stdio.h>
#include <string.h>
#include <bench.h>

// Global constants
extern const digit_t p[NWORDS_FIELD];

// Benchmark and test parameters  
static int BENCH_LOOPS = 100000;       // Number of iterations per bench
static int TEST_LOOPS  = 100000;       // Number of iterations per test


bool fp2_test()
{ // Tests for the GF(p^2) arithmetic
    bool OK = true;
    int n, passed;
    fp2_t a, b, c, d, e, f, ma, mb, mc, md, me, mf;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing arithmetic over GF(p^2): \n\n"); 

    // Addition in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random_test(&a); fp2random_test(&b); fp2random_test(&c); fp2random_test(&d); 

        fp2_add(&d, &a, &b); fp2_add(&e, &d, &c);                 // e = (a+b)+c
        fp2_add(&d, &b, &c); fp2_add(&f, &d, &a);                 // f = a+(b+c)
        if (compare_words((digit_t*)&e, (digit_t*)&f, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_add(&d, &a, &b);                                      // d = a+b 
        fp2_add(&e, &b, &a);                                      // e = b+a
        if (compare_words((digit_t*)&d, (digit_t*)&e, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_set(&b, 0);
        fp2_add(&d, &a, &b);                                      // d = a+0 
        if (compare_words((digit_t*)&a, (digit_t*)&d, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_set(&b, 0);   
        fp2_neg(&d, &a);                      
        fp2_add(&e, &a, &d);                                      // e = a+(-a)
        if (compare_words((digit_t*)&e, (digit_t*)&b, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) addition tests ............................................ PASSED");
    else { printf("  GF(p^2) addition tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Subtraction in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random_test(&a); fp2random_test(&b); fp2random_test(&c); fp2random_test(&d);

        fp2_sub(&d, &a, &b); fp2_sub(&e, &d, &c);                 // e = (a-b)-c
        fp2_add(&d, &b, &c); fp2_sub(&f, &a, &d);                 // f = a-(b+c)
        if (compare_words((digit_t*)&e, (digit_t*)&f, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_sub(&d, &a, &b);                                      // d = a-b 
        fp2_sub(&e, &b, &a);
        fp2_neg(&e, &e);                                          // e = -(b-a)
        if (compare_words((digit_t*)&d, (digit_t*)&e, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_set(&b, 0);
        fp2_sub(&d, &a, &b);                                      // d = a-0 
        if (compare_words((digit_t*)&a, (digit_t*)&d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        fp2_set(&b, 0);              
        fp2_sub(&e, &a, &a);                                      // e = a+(-a)
        if (compare_words((digit_t*)&e, (digit_t*)&b, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) subtraction tests ......................................... PASSED");
    else { printf("  GF(p^2) subtraction tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Multiplication in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fp2random_test(&a); fp2random_test(&b); fp2random_test(&c);
        
        fp2_tomont(&ma, &a);
        fp2_frommont(&c, &ma);
        if (compare_words((digit_t*)&a, (digit_t*)&c, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_tomont(&ma, &a); fp2_tomont(&mb, &b); fp2_tomont(&mc, &c);
        fp2_mul(&md, &ma, &mb); fp2_mul(&me, &md, &mc);                          // e = (a*b)*c
        fp2_mul(&md, &mb, &mc); fp2_mul(&mf, &md, &ma);                          // f = a*(b*c)
        fp2_frommont(&e, &me);
        fp2_frommont(&f, &mf);
        if (compare_words((digit_t*)&e, (digit_t*)&f, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_tomont(&ma, &a); fp2_tomont(&mb, &b); fp2_tomont(&mc, &c); 
        fp2_add(&md, &mb, &mc); fp2_mul(&me, &ma, &md);                          // e = a*(b+c)
        fp2_mul(&md, &ma, &mb); fp2_mul(&mf, &ma, &mc); fp2_add(&mf, &md, &mf);  // f = a*b+a*c
        fp2_frommont(&e, &me);
        fp2_frommont(&f, &mf);
        if (compare_words((digit_t*)&e, (digit_t*)&f, 2*NWORDS_FIELD)!=0) { passed=0; break; }
     
        fp2_tomont(&ma, &a); fp2_tomont(&mb, &b);
        fp2_mul(&md, &ma, &mb);                                                  // d = a*b 
        fp2_mul(&me, &mb, &ma);                                                  // e = b*a 
        fp2_frommont(&d, &md);
        fp2_frommont(&e, &me);
        if (compare_words((digit_t*)&d, (digit_t*)&e, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_tomont(&ma, &a);
        fp2_set(&b, 1); fp2_tomont(&mb, &b);
        fp2_mul(&md, &ma, &mb);                                                  // d = a*1  
        fp2_frommont(&a, &ma);
        fp2_frommont(&d, &md);                
        if (compare_words((digit_t*)&a, (digit_t*)&d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
       
        fp2_set(&b, 0);
        fp2_tomont(&mb, &b);
        fp2_mul(&md, &ma, &mb);                                                  // d = a*0 
        fp2_frommont(&d, &md);                
        if (compare_words((digit_t*)&b, (digit_t*)&d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) multiplication tests ...................................... PASSED");
    else { printf("  GF(p^2) multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Squaring in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random_test(&a);
        
        fp2_tomont(&ma, &a);
        fp2_sqr(&mb, &ma);                                          // b = a^2
        fp2_mul(&mc, &ma, &ma);                                     // c = a*a 
        fp2_frommont(&b, &mb);
        fp2_frommont(&c, &mc);
        if (compare_words((digit_t*)&b, (digit_t*)&c, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2_set(&a, 0); fp2_tomont(&ma, &a);
        fp2_sqr(&md, &ma);                                          // d = 0^2 
        if (compare_words((digit_t*)&ma, (digit_t*)&md, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) squaring tests............................................. PASSED");
    else { printf("  GF(p^2) squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Inversion in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random_test(&a);

        fp2_tomont(&ma, &a);
        fp2_set(&d, 1);
        memcpy(&mb, &ma, RADIX/8 * 2*NWORDS_FIELD);
        fp2_inv(&ma);
        fp2_mul(&mc, &ma, &mb);                                     // c = a*a^-1 
        fp2_frommont(&c, &mc);
        if (compare_words((digit_t*)&c, (digit_t*)&d, 2*NWORDS_FIELD) != 0) { passed = 0; break; }

        fp2_set(&a, 0);
        fp2_set(&d, 0);
        fp2_inv(&a);                                                // c = 0^-1
        if (compare_words((digit_t*)&a, (digit_t*)&d, 2*NWORDS_FIELD) != 0) { passed = 0; break; }
    }
    if (passed == 1) printf("  GF(p^2) inversion tests............................................ PASSED");
    else { printf("  GF(p^2) inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Square root and square detection in GF(p^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random_test(&a);

        fp2_tomont(&ma, &a);
        fp2_sqr(&mc, &ma);
        fp2_frommont(&c, &mc);                                      // c = a^2
        if (fp2_is_square(&mc) != 1) { passed = 0; break; }        

        fp2_sqrt(&mc);                                              // c = a = sqrt(c) 
        fp2_neg(&md, &mc);
        fp2_frommont(&c, &mc);
        fp2_frommont(&d, &md);
        if ((compare_words((digit_t*)&a, (digit_t*)&c, 2*NWORDS_FIELD) != 0) & (compare_words((digit_t*)&a, (digit_t*)&d, 2*NWORDS_FIELD) != 0)) { passed = 0; break; }
    }
    if (passed == 1) printf("  Square root, square tests.......................................... PASSED");
    else { printf("  Square root, square tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    return OK;
}

// https://cplusplus.com/reference/cstdlib/qsort/
// For sorting arrays during the benchmarks
int compare (const void * a, const void * b)
{
  return ( *(long*)a - *(long*)b );
}


bool fp2_run()
{
    bool OK = true;
    int n, i;
    unsigned long long cycles, cycles1, cycles2;
    fp2_t a, b;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking GF(p^2) field arithmetic: \n\n"); 

    // GF(p^2) addition
    fp2random_test(&a);
    fp2random_test(&b);
    long cycle_runs[10];

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p^2) addition runs in .......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / (6 * BENCH_LOOPS), a.re[0]);

    // GF(p^2) subtraction
    fp2random_test(&a);
    fp2random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p^2) subtraction runs in ....................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / (6 * BENCH_LOOPS), a.re[0]);

    // GF(p^2) multiplication
    fp2random_test(&a);
    fp2random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p^2) multiplication runs in .................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / (6 * BENCH_LOOPS), a.re[0]);

    // GF(p^2) inversion
    fp2random_test(&a);
    fp2random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_inv(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p^2) inversion runs in ......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / BENCH_LOOPS, a.re[0]);

    // GF(p^2) sqrt
    fp2random_test(&a);
    fp2random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_sqrt(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p^2) sqrt runs in .............................................. %ld cycles, (%llu ignore me)\n", cycle_runs[4] / BENCH_LOOPS, a.re[0]);

    // GF(p^2) is_square
    fp2random_test(&a);
    fp2random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp2_is_square(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  Square checking runs in ........................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / BENCH_LOOPS, a.re[0]);

    return OK;
}


int main(int argc, char* argv[])
{
    if (argc < 3) {
        printf("Please enter an argument: 'test' or 'bench' and <reps>\n");
        exit(1);
    }
    if (!strcmp(argv[1], "test")) {
        TEST_LOOPS = atoi(argv[2]);
        return !fp2_test();
    } else if (!strcmp(argv[1], "bench")) {
        BENCH_LOOPS = atoi(argv[2]);
        return !fp2_run();
    } else {
        exit(1);
    }
}
