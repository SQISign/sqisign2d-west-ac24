#include "bench.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>

// Benchmark and test parameters  
static int BENCH_LOOPS = 1000;         // Number of iterations per bench
static int TEST_LOOPS  = 100000;       // Number of iterations per test


bool fp_test(void)
{ // Tests for the field arithmetic
    bool OK = true;
    int n, passed;
    fp_t a, b, c, d, e, f;
    // fp_t ma, mb;


    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing field arithmetic over GF(p): \n\n"); 

    // Test equality
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp_random_test(&a);
        fp_add(&b, &a, (fp_t *)&ONE);
        fp_set_zero(&c);

        if (fp_equal(&a, &a) == 0) { passed=0; break; }
        if (fp_equal(&a, &b) != 0) { passed=0; break; }
        if (fp_equal(&c, (fp_t *)&ZERO) == 0) { passed=0; break; }

        if (fp_is_zero((fp_t *)&ZERO) == 0) { passed=0; break; }
        if (fp_is_zero((fp_t *)&ONE) != 0)  { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) equality tests ............................................ PASSED");
    else { printf("  GF(p) equality tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field addition
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp_random_test(&a); fp_random_test(&b); fp_random_test(&c); fp_random_test(&d); 

        fp_add(&d, &a, &b); fp_add(&e, &d, &c);                 // e = (a+b)+c
        fp_add(&d, &b, &c); fp_add(&f, &d, &a);                 // f = a+(b+c)
        if (fp_equal(&e, &f) == 0) { passed=0; break; }

        fp_add(&d, &a, &b);                                  // d = a+b 
        fp_add(&e, &b, &a);                                  // e = b+a
        if (fp_equal(&d, &e) == 0 ) { passed=0; break; }

        fp_set_zero(&b);
        fp_add(&d, &a, &b);                                  // d = a+0 
        if (fp_equal(&a, &d) == 0 ) { passed=0; break; }

        fp_set_zero(&b);   
        fp_neg(&d, &a);                      
        fp_add(&e, &a, &d);                                  // e = a+(-a)
        if (fp_equal(&e, &b) == 0 ) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) addition tests ............................................ PASSED");
    else { printf("  GF(p) addition tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field subtraction
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp_random_test(&a); fp_random_test(&b); fp_random_test(&c); fp_random_test(&d);

        fp_sub(&d, &a, &b); fp_sub(&e, &d, &c);                 // e = (a-b)-c
        fp_add(&d, &b, &c); fp_sub(&f, &a, &d);                 // f = a-(b+c)
        if (fp_equal(&e, &f) == 0) { passed=0; break; }

        fp_sub(&d, &a, &b);                                  // d = a-b 
        fp_sub(&e, &b, &a);
        fp_neg(&e, &e);                                      // e = -(b-a)
        if (fp_equal(&d, &e) == 0) { passed=0; break; }

        fp_set_zero(&b);
        fp_sub(&d, &a, &b);                                  // d = a-0 
        if (fp_equal(&a, &d) == 0) { passed=0; break; }
                
        fp_sub(&e, &a, &a);                                  // e = a+(-a)
        if (fp_is_zero(&e) == 0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) subtraction tests ......................................... PASSED");
    else { printf("  GF(p) subtraction tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Montgomery conversion
    // TODO: this fails currently
    // I must have a stupid bug
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fp_random_test(&a); fp_random_test(&b); fp_random_test(&c);

        // fp_from_mont(&b, &a);
        // fp_to_mont(&c, &b);
        // if (fp_equal(&a, &c) == 0) { passed=0; break; }

        fp_to_mont(&b, &a);
        fp_from_mont(&c, &b);
        if (fp_equal(&a, &c) == 0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) montgomery conversion tests ............................... PASSED");
    else { printf("  GF(p) montgomery conversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    

    // Field multiplication
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fp_random_test_mont(&a); fp_random_test_mont(&b); fp_random_test_mont(&c);
        fp_mul(&d, &a, &b); fp_mul(&e, &d, &c);                          // e = (a*b)*c
        fp_mul(&d, &b, &c); fp_mul(&f, &d, &a);                          // f = a*(b*c)
        if (fp_equal(&e, &f) == 0) { passed=0; break; }

        fp_add(&d, &b, &c); fp_mul(&e, &a, &d);                          // e = a*(b+c)
        fp_mul(&d, &a, &b); fp_mul(&f, &a, &c); fp_add(&f, &d, &f);      // f = a*b+a*c
        if (fp_equal(&e, &f) == 0) { passed=0; break; }
     
        fp_mul(&d, &a, &b);                                              // d = a*b 
        fp_mul(&e, &b, &a);                                              // e = b*a 
        if (fp_equal(&d, &e) == 0) { passed=0; break; }

        fp_set_one(&b);
        fp_mul(&d, &a, &b);                                              // d = a*1               
        if (fp_equal(&a, &d) == 0) { passed=0; break; }
       
        fp_set_zero(&b);
        fp_mul(&d, &a, &b);                                              // d = a*0           
        if (fp_equal(&b, &d) == 0) { passed=0; break; } 
    }
    if (passed==1) printf("  GF(p) multiplication tests ...................................... PASSED");
    else { printf("  GF(p) multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field squaring
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp_random_test_mont(&a);

        fp_sqr(&b, &a);                                    // b = a^2
        fp_mul(&c, &a, &a);                               // c = a*a 
        if (fp_equal(&b, &c) == 0) { passed=0; break; } 

        fp_set_zero(&a); 
        fp_sqr(&d, &a);                                  // d = 0^2 
        if (fp_is_zero(&d) == 0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) squaring tests............................................. PASSED");
    else { printf("  GF(p) squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field inversion
    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {
        fp_random_test(&a);
        
        fp_set(&b, &a);
        fp_inv(&b);
        fp_mul(&c, &a, &b);                               // c = a*a^-1 
        if (fp_equal(&c, (fp_t*)&ONE) == 0) { passed=0; break; } 


        fp_set_zero(&a);
        fp_inv(&a);                                        // c = 0^-1
        if (fp_equal(&a, (fp_t*)&ZERO) == 0) { passed=0; break; } 
    }
    if (passed == 1) printf("  GF(p) inversion tests............................................ PASSED");
    else { printf("  GF(p) inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Square root and square detection
    passed = 1;
    for (n = 0; n < TEST_LOOPS; n++)
    {
        fp_random_test(&a);

        fp_sqr(&c, &a);                       // c = a^2
        if (fp_is_square(&c) != 1) { passed = 0; break; }

        fp_sqrt(&c);                           // c, d = ±sqrt(c) 
        fp_neg(&d, &c);
        if ((fp_equal(&a, &c) == 0) && (fp_equal(&a, &d) == 0)) { passed = 0; break; }
    }
    if (passed == 1) printf("  Square root, square tests........................................ PASSED");
    else { printf("  Square root, square tests... FAILED"); printf("\n"); return false; }
    printf("\n");
 
    return OK;
}


bool fp_run(void)
{
    bool OK = true;
    int n, i;
    unsigned long long cycles, cycles1, cycles2;
    fp_t a, b;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking GF(p) field arithmetic: \n\n"); 

    // GF(p) addition
    fp_random_test(&a);
    fp_random_test(&b);
    long cycle_runs[10];

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) addition runs in .......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / (6 * BENCH_LOOPS), a.w[0]);

    // GF(p) subtraction
    fp_random_test(&a);
    fp_random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) subtraction runs in ....................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] /  (6 * BENCH_LOOPS), a.w[0]);

    // GF(p) multiplication
    fp_random_test(&a);
    fp_random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) multiplication runs in .................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] /  (6 * BENCH_LOOPS), a.w[0]);

    // GF(p) squaring
    fp_random_test(&a);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_sqr(&a, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) squaring runs in .......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] /  (BENCH_LOOPS), a.w[0]);

  // GF(p) inversion
    fp_random_test(&a);
    fp_random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_inv(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) inversion runs in ......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] /  BENCH_LOOPS, a.w[0]);


    // GF(p) sqrt
    fp_random_test(&a);
    fp_random_test(&b);

    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_sqrt(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  GF(p) sqrt runs in .............................................. %ld cycles, (%llu ignore me)\n", cycle_runs[4] / BENCH_LOOPS, a.w[0]);

    // GF(p) is_square
    fp_random_test(&a);
    fp_random_test(&b);
    
    for (i=0; i<10; i++){
        cycles = 0;
        cycles1 = cpucycles(); 
        for (n=0; n<BENCH_LOOPS; n++)
        {
            fp_is_square(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2-cycles1;
    }
    qsort(cycle_runs, 10, sizeof(long), compare);
    printf("  Square checking runs in ......................................... %ld cycles, (%llu ignore me)\n", cycle_runs[4] / BENCH_LOOPS, a.w[0]);

    return OK;

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
        return !fp_test();
    } else if (!strcmp(argv[1], "bench")) {
        BENCH_LOOPS = atoi(argv[2]);
        return !fp_run();
    } else {
        exit(1);
    }
}
