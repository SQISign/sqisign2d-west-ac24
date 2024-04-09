#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include "test_sqisigndim2.h"
#include <tools.h>

static inline int64_t cpucycles(void) {
#if (defined(TARGET_AMD64) || defined(TARGET_X86))
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d"(hi));
    return ((int64_t) lo) | (((int64_t) hi) << 32);
#elif (defined(TARGET_S390X))
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q" (tod) : : "cc");
    return (tod * 1000 / 4096);
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}

static __inline__ uint64_t rdtsc(void)
{
    return (uint64_t) cpucycles();
}




bool curve_is_canonical(ec_curve_t const *E)
{
    ec_curve_t EE;
    ec_isom_t isom;
    ec_curve_normalize(&EE, &isom, E);

    fp2_t lhs, rhs;
    fp2_mul(&lhs, &E->A, &EE.C);
    fp2_mul(&rhs, &E->C, &EE.A);
    return fp2_is_equal(&lhs, &rhs);
}

int bench_fp2_operations(int repeat) {

    uint64_t t0, t1;

    fp2_t a = BASIS_EVEN.P.x;
    fp2_t b = BASIS_EVEN.PmQ.x;
    t0 = rdtsc();
    for (int i=0;i<repeat;i++) {
        fp2_add(&a,&b,&a);
    }
    t1=  rdtsc();
    printf("\x1b[34m fp2_add %'" PRIu64 " cycles\x1b[0m\n", (t1-t0)/repeat );

    a = BASIS_EVEN.P.x;
    b = BASIS_EVEN.Q.x;

    t0 = rdtsc();
    for (int i=0;i<repeat;i++) {
        fp2_sqr(&a,&a);
    }
    t1=  rdtsc();
    printf("\x1b[34m fp2_sqr %'" PRIu64 " cycles\x1b[0m\n", (t1-t0)/repeat );

    a = BASIS_EVEN.P.x;
    b = BASIS_EVEN.PmQ.x;
    
    t0 = rdtsc();
    for (int i=0;i<repeat;i++) {
        fp2_mul(&a,&b,&a);
    }
    t1=  rdtsc();
    printf("\x1b[34m fp2_mul %'" PRIu64 " cycles\x1b[0m\n", (t1-t0)/repeat);

    t0 = rdtsc();
    for (int i=0;i<repeat;i++) {
        fp2_inv(&a);
        fp2_mul(&a,&b,&a);
    }
    t1=  rdtsc();
    printf("\x1b[34m fp2_inv %'" PRIu64 " cycles\x1b[0m\n", (t1-t0)/repeat);

}


int test_sqisign(int repeat)
{
    int res = 1;

    int num_sig =10;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    unsigned char msg[32] = { 0 };

    public_key_init(&pk);
    secret_key_init(&sk);
    secret_sig_init(&(sig));
    

    


    // printf("Printing details of first signature\n");
    // int val = protocols_sign(&sig, &pk, &sk, msg, 32, 1);
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1;
    clock_t t;


    printf("\n\nTesting signatures\n");
    for (int i = 0; i < repeat; ++i)
    {   

        printf("#%d \n",i);
        t = tic();
        t0=rdtsc();
        protocols_keygen(&pk, &sk);
        t1 = rdtsc();
        TOC(t,"Keygen");
        printf("\x1b[34mkeygen  took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);

        t = tic();
        t0 = rdtsc();
        int val = protocols_sign(&sig, &pk, &sk, msg, 32, 0);
        t1 = rdtsc();
        TOC(t,"Signing");
        printf("\x1b[34msigning took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);

        t = tic();
        t0 = rdtsc();
        int check = protocols_verif(&sig,&pk,msg,32);
        if (!check) {
            printf("verif failed ! \n");
        } 
        t1 = rdtsc();
        TOC(t,"Verification");
        printf("\x1b[34mverif   took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);    
    
        printf(" \x1b[35mfull\x1b[0m signature was: %s\n\n", check ? "\x1b[32mvalid\x1b[0m" : "\x1b[31minvalid\x1b[0m");
    }


    // TOC(t, "protocols_sign");
    // float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    // printf("average signing time [%.2f ms]\n", (float) (ms/repeat));

    // clock_t t_verif = tic();
    // for (int i = 0; i < repeat; ++i)
    // {   
        
    // }
    // TOC(t_verif,"protocols_verif");  
    // ms = (1000. * (float) (clock() - t_verif) / CLOCKS_PER_SEC);
    // printf("average verif time [%.2f ms]\n", (float) (ms/repeat));

    

    public_key_finalize(&pk);
    secret_key_finalize(&sk);
    secret_sig_finalize(&sig);

    return res;
}

// run all tests in module
int main(){
    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();

    printf("\nRunning sqisigndim2 tests\n \n");

    bench_fp2_operations(1000);

    res &= test_sqisign(3);

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("All tests passed!\n");
    }
    return(!res);
}
