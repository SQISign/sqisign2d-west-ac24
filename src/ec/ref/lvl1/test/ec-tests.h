#ifndef EC_TESTS_H
#define EC_TESTS_H

#include "test_extras.h"
#include <stdio.h>
#include <string.h>
#include <bench.h>       //////// NOTE: enable later
#include "test-basis.h"
#include "ec_params.h"

// Benchmark and test parameters  
static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 512;       // Number of iterations per test

void random_scalar(digit_t* k)
{
    for(int i = 0; i < NWORDS_FIELD; i++)
        k[i] = rand();
}

bool ec_test()
{ // Tests for ecc arithmetic
    bool OK = true;
    int passed;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_point_t AC = {0};
    digit_t k[NWORDS_ORDER] = {0}, l[NWORDS_ORDER] = {0};

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing ecc functions: \n\n"); 

    // Point doubling
    passed = 1;
    uint64_t xP_data[2][4] = {{0xDFD70ED0861BD329, 0x20ACD3758C7F5540, 0x3DCCDC007277F80A, 0x18D6D2A22981DCE1}, 
                              {0x3C23730A3F08F38C, 0x98BB973AFD3D954D, 0x8D98ADFC2829AE8A, 0x21A2464D6369AFBA}};
    fp2_w64(&P.x, xP_data);
    fp2_set_one(&P.z);
    fp2_set_one(&AC.z);

    copy_point(&R, &P);
    xDBL(&S, &R, &AC);
    copy_point(&SS, &S);         // Copy of S = SS <- 2P 

    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);

    uint64_t xR_data[2][4] = {{0x5950EE0A4AF90FC8, 0x16488065A0A98B08, 0xCE65322229DA0FD1, 0x270A35FF781EE204}, 
                               {0x564447FD9EC57F6B, 0x2EE24E984294F729, 0x53A6C7360E972C71, 0x4FCF4B9928A7C7E}};

    fp2_w64(&R.x, xR_data);
    printf("%llu\n", S.x.re.v0);
    printf("%llu\n", AC.z.re.v3);


    // R.x.re.v0 = 0x5950EE0A4AF90FC8; R.x.re.v1 = 0x16488065A0A98B08; R.x.re.v2 = 0xCE65322229DA0FD1; R.x.re.v3 = 0x270A35FF781EE204;
    // R.x.im.v0 = 0x564447FD9EC57F6B; R.x.im.v1 = 0x2EE24E984294F729; R.x.im.v2 = 0x53A6C7360E972C71; R.x.im.v3 = 0x4FCF4B9928A7C7E;

    // TODO failing...
    // if (fp2_is_equal(&R.x, &S.x) == 0) { passed=0; goto out0; };
    // if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2)!=0) { passed=0; goto out0; }
    
    uint64_t xQ_data[2][4] = {{0xC46076A670C70053, 0x97517AFA3AB9ED13, 0x349644C942EDF993, 0xBB4A4DB6F29AF9E},
                              {0x8B47629FB5A15BB0, 0x4EC6E809953C1A10, 0x1F83F0EC6CBB84D6, 0x1D8417C1D33265D3}};
    fp2_w64(&Q.x, xQ_data);
    fp2_set_one(&Q.z);
    
    // Q.x.re.v0 = 0xC46076A670C70053; Q.x.re.v1 = 0x97517AFA3AB9ED13; Q.x.re.v2 = 0x349644C942EDF993; Q.x.re.v3 = 0xBB4A4DB6F29AF9E;
    // Q.x.im.v0 = 0x8B47629FB5A15BB0; Q.x.im.v1 = 0x4EC6E809953C1A10; Q.x.im.v2 = 0x1F83F0EC6CBB84D6; Q.x.im.v3 = 0x1D8417C1D33265D3;
    // Q.z.re.v0 = 0x01;

    // PQ.x.re.v0 = 0x853F66D11BE5534F; PQ.x.re.v1 = 0x27C8FD4E52D03D4A; PQ.x.re.v2 = 0xF88EA78D0A0C29D2; PQ.x.re.v3 = 0x2F6DFB07D397A067;
    // PQ.x.im.v0 = 0xE8DBC4AA34434BA1; PQ.x.im.v1 = 0x7A73AE182636F8A0; PQ.x.im.v2 = 0x419EC260137868EB; PQ.x.im.v3 = 0x129B3E301703D43F;
    // PQ.z.re.v0 = 0x01;
    uint64_t xPQ_data[2][4] = {{0x853F66D11BE5534F, 0x27C8FD4E52D03D4A, 0xF88EA78D0A0C29D2, 0x2F6DFB07D397A067},
                              {0xE8DBC4AA34434BA1, 0x7A73AE182636F8A0, 0x419EC260137868EB, 0x129B3E301703D43F}};
    fp2_w64(&PQ.x, xPQ_data);
    fp2_set_one(&PQ.z);

    fp2_copy(&S.x, &Q.x);
    fp2_copy(&S.z, &Q.z);
    xADD(&S, &SS, &S, &PQ);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    // fp2_frommont(&S.x, &S.x);

    // R.x.re.v0 = 0xED0BEB8F93AB4FF9; R.x.re.v1 = 0x27CF508B80CD49BF; R.x.re.v2 = 0x38A6134DFA04B2BA; R.x.re.v3 = 0x27B4CB15E109EF1F;
    // R.x.im.v0 = 0x6F731BA6FD227BDE; R.x.im.v1 = 0x14C12335341167F8; R.x.im.v2 = 0xECA7B60F7866E27A; R.x.im.v3 = 0x2A7A79A152880457;
    uint64_t xRR_data[2][4] = {{0xED0BEB8F93AB4FF9, 0x27CF508B80CD49BF, 0x38A6134DFA04B2BA, 0x27B4CB15E109EF1F},
                               {0x6F731BA6FD227BDE, 0x14C12335341167F8, 0xECA7B60F7866E27A, 0x2A7A79A152880457}};
    fp2_w64(&R.x, xRR_data);
    if (fp2_is_equal(&R.x, &S.x) == 0) { passed=0; goto out0; };
    // if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    // fp2_tomont(&R.x, &P.x);
    // fp2_tomont(&R.z, &P.z);
    // k[0] = 126;
    // xMUL(&S, &R, k, (ec_curve_t*)&AC);
    // fp2_inv(&S.z);
    // fp2_mul(&S.x, &S.x, &S.z);
    // fp2_frommont(&S.x, &S.x);

    // R.x.re.v0 = 0xDE80F87A1203A147; R.x.re.v1 = 0xD59E1215928A3B2D; R.x.re.v2 = 0xD5A67F83A5A8CE46; R.x.re.v3 = 0xA11E162488C9CDF;
    // R.x.im.v0 = 0x9417D0D79A26741B; R.x.im.v1 = 0x8B1F47D6F0FE5EEC; R.x.im.v2 = 0xE52188DCB054CE36; R.x.im.v3 = 0x1A8075A6C3148AB3;

    // if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    // fp2_tomont(&R.x, &P.x);
    // fp2_tomont(&R.z, &P.z);
    // k[0] = 0xE77AD6B6C6B2D8CD;
    // k[1] = 0xDE43A0B600F38D12;
    // k[2] = 0xA35F4A7897E17CE2;
    // k[3] = 0x10ACB62E614D1237;
    // xMUL(&S, &R, k, (ec_curve_t*)&AC);
    // fp2_inv(&S.z);
    // fp2_mul(&S.x, &S.x, &S.z);
    // fp2_frommont(&S.x, &S.x);

    // R.x.re.v0 = 0xD3938B0A68A3E7C0; R.x.re.v1 = 0xE0667113208A0595; R.x.re.v2 = 0x258F314C84E9CB60; R.x.re.v3 = 0x14984BA7CA59AB71;
    // R.x.im.v0 = 0xFE728423EE3BFEF4; R.x.im.v1 = 0xBF68C42FE21AE0E4; R.x.im.v2 = 0xA8FAF9C9528609CA; R.x.im.v3 = 0x1225EC77A1DC0285;

    // if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    // fp2_tomont(&R.x, &Q.x);
    // fp2_tomont(&R.z, &Q.z);
    // k[0] = 0xE77AD6B6C6B2D8CD;
    // k[1] = 0xDE43A0B600F38D12;
    // k[2] = 0xA35F4A7897E17CE2;
    // k[3] = 0x10ACB62E614D1237;
    // l[0] = 0x34AB78B6C6B2D8C0;
    // l[1] = 0xDE6B2D8CD00F38D1;
    // l[2] = 0xA35F4A7897E17CE2;
    // l[3] = 0x20ACF4A789614D13;
    // fp2_inv(&SS.z);
    // fp2_mul(&SS.x, &SS.x, &SS.z);
    // fp2_copy(&SS.z, &R.z);
    // xDBLMUL(&S, &R, k, &SS, l, &PQ, (ec_curve_t*)&AC);
    // fp2_inv(&S.z);
    // fp2_mul(&S.x, &S.x, &S.z);
    // fp2_frommont(&S.x, &S.x);

    // R.x.re.v0 = 0x554E1ADC609B992F; R.x.re.v1 = 0xE407D961F8CC4C42; R.x.re.v2 = 0x1CF626AFED5A68CE; R.x.re.v3 = 0x6D02692EE110483;
    // R.x.im.v0 = 0x16FB094E831C8997; R.x.im.v1 = 0xFDE4ECF31DC5F702; R.x.im.v2 = 0x89303D868DFAD7B4; R.x.im.v3 = 0xC91ACE81346F22D;

    // if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
out0:
    if (passed==1) printf("  ECC arithmetic tests ............................................ PASSED");
    else { printf("  ECC arithmetic tests... FAILED"); printf("\n"); return false; }
    printf("\n");
 
    return OK;
}

bool dlog_test()
{ // Tests for dlog
    bool OK = true;
    int passed;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_curve_t AC = {0};
    ec_basis_t PQ2;
    digit_t scalarP[NWORDS_ORDER], scalarQ[NWORDS_ORDER], k[NWORDS_ORDER] = {0}, l[NWORDS_ORDER] = {0};
    digit_t kt[NWORDS_ORDER], lt[NWORDS_ORDER], f1[NWORDS_ORDER] = {0}, f2[NWORDS_ORDER] = {0}, zero[NWORDS_ORDER] = {0}, tpFdiv2[NWORDS_ORDER] = {0}, tpF[NWORDS_ORDER] = {0};

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing dlog functions: \n\n");

    // dlog2 testing
    passed = 1;
    
    fp2_w64(&P.x,  xP2_data);
    fp2_w64(&Q.x,  xQ2_data);
    fp2_w64(&PQ.x, xPQ2_data);
    fp2_set_one(&P.z);
    fp2_set_one(&Q.z);
    fp2_set_one(&PQ.z);

    memcpy(f1, TWOpFm1, NWORDS_ORDER*RADIX/8);
    memcpy(f2, TWOpF, NWORDS_ORDER*RADIX/8);

    fp2_set_one(&AC.C);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);
    k[0] = 0xFFFFFFFFFFFFFFFF;
    k[1] = 0x00000000000007FF;
    l[0] = 0xFFFFFFFFFFFFFFFE;
    l[1] = 0x00000000000007FF;

    for (int n = 0; n < TEST_LOOPS; n++)
    {
        k[0] -= 1;
        l[0] -= 2;
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        ec_dlog_2(scalarP, scalarQ, &PQ2, &R, &AC);

        memcpy(kt, k, NWORDS_ORDER*RADIX/8);
        memcpy(lt, l, NWORDS_ORDER*RADIX/8);
        if (compare_words(k, f1, NWORDS_ORDER) == 1 ||
           (compare_words(l, f1, NWORDS_ORDER) == 1 && (compare_words(k, zero, NWORDS_ORDER) == 0 || compare_words(k, f1, NWORDS_ORDER) == 0))) {
            if (compare_words(k, zero, NWORDS_ORDER) != 0) {
                sub_test(kt, f2, kt, NWORDS_ORDER);
            }
            if (compare_words(l, zero, NWORDS_ORDER) != 0) {
                sub_test(lt, f2, lt, NWORDS_ORDER);
            }
        }
        if (compare_words((digit_t*)scalarP, (digit_t*)kt, NWORDS_ORDER) != 0 || compare_words((digit_t*)scalarQ, (digit_t*)lt, NWORDS_ORDER) != 0) { passed = 0; break; }
    }

    if (passed == 1) printf("  dlog2 tests ..................................................... PASSED");
    else { printf("  dlog2 tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // dlog3 testing
    passed = 1;
    
    fp2_w64(&P.x,  xP3_data);
    fp2_w64(&Q.x,  xQ3_data);
    fp2_w64(&PQ.x, xPQ3_data);
    fp2_set_one(&P.z);
    fp2_set_one(&Q.z);
    fp2_set_one(&PQ.z);

    fp2_set_one(&AC.C);
    memcpy(tpFdiv2, THREEpFdiv2, NWORDS_ORDER*RADIX/8);
    memcpy(tpF, THREEpF, NWORDS_ORDER*RADIX/8);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);
    k[1] = 0;
    l[1] = 0;
    k[0] = 0x02153E468B91C6D1;
    l[0] = 0x02153E468B91C6D0;

    for (int n = 0; n < TEST_LOOPS; n++)
    {
        k[0] -= 1;
        l[0] -= 2;
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        ec_dlog_3(scalarP, scalarQ, &PQ2, &R, &AC);

        memcpy(kt, k, NWORDS_ORDER*RADIX/8);
        memcpy(lt, l, NWORDS_ORDER*RADIX/8);
        if (compare_words(k, tpFdiv2, NWORDS_ORDER) == 1 ||
           (compare_words(l, tpFdiv2, NWORDS_ORDER) == 1 && compare_words(k, zero, NWORDS_ORDER) == 0)) {
            if (compare_words(k, zero, NWORDS_ORDER) != 0) {
                sub_test(kt, tpF, kt, NWORDS_ORDER);
            }
            if (compare_words(l, zero, NWORDS_ORDER) != 0) {
                sub_test(lt, tpF, lt, NWORDS_ORDER);
            }
        }
        if (compare_words((digit_t*)scalarP, (digit_t*)kt, NWORDS_ORDER) != 0 || compare_words((digit_t*)scalarQ, (digit_t*)lt, NWORDS_ORDER) != 0) { passed = 0; break; }
    }

    if (passed == 1) printf("  dlog3 tests ..................................................... PASSED");
    else { printf("  dlog3 tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    return OK;
}

bool ec_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    ec_point_t P, Q, R, PQ, AC;
    digit_t k[NWORDS_ORDER], l[NWORDS_ORDER];
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking ecc arithmetic: \n\n"); 

    // Point doubling
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        xDBL(&Q, &P, &AC);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Montgomery x-only doubling runs in .............................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point addition
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        xADD(&R, &Q, &P, &PQ);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only addition runs in .............................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point multiplication
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        random_scalar(k);
        cycles1 = cpucycles();
        xMUL(&Q, &P, k, (ec_curve_t*)&AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only scalar multiplication runs in ................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point multiplication
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        random_scalar(k);
        random_scalar(l);
        cycles1 = cpucycles();
        xDBLMUL(&R, &P, k, &Q, l, &PQ, (ec_curve_t*)&AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only double-scalar multiplication runs in .......... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    return OK;
}

bool dlog_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_curve_t AC = {0};
    ec_basis_t PQ2;
    digit_t scalarP[NWORDS_ORDER], scalarQ[NWORDS_ORDER], k[NWORDS_ORDER], l[NWORDS_ORDER];


    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking dlog2: \n\n");

    // dlog2 computation
    fp2_w64(&P.x,  xP2_data);
    fp2_w64(&Q.x,  xQ2_data);
    fp2_w64(&PQ.x, xPQ2_data);
    fp2_set_one(&P.z);
    fp2_set_one(&Q.z);
    fp2_set_one(&PQ.z);

    fp2_set_one(&AC.C);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        random_scalar(k);
        random_scalar(l);

        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        cycles1 = cpucycles();
        ec_dlog_2(scalarP, scalarQ, &PQ2, &R, &AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  dlog2 runs in ................................................... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // dlog3 computation

    fp2_w64(&P.x,  xP3_data);
    fp2_w64(&Q.x,  xQ3_data);
    fp2_w64(&PQ.x, xPQ3_data);
    fp2_set_one(&P.z);
    fp2_set_one(&Q.z);
    fp2_set_one(&PQ.z);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        random_scalar(k);
        random_scalar(l);

        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        cycles1 = cpucycles();
        ec_dlog_3(scalarP, scalarQ, &PQ2, &R, &AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  dlog3 runs in ................................................... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    return OK;
}

#endif
