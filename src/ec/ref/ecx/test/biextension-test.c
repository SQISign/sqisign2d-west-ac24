#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "ec.h"
#include "biextension.h"

// VERY NOT SECURE (testing only)
void fp2_random(fp2_t *a){
    for(int i = 0; i < NWORDS_FIELD; i++){
        a->re[i] = rand();
        a->im[i] = rand();
    }
    // Normalize
    fp2_t one;
    fp_mont_setone(one.re);fp_set(one.im,0);
    fp2_mul(&*a, &*a, &one);
    // Update seed
    srand((unsigned) a->re[0]);
}


int biextension_test() {
    E0=CURVE_E0;
    A24=CURVE_E0_A24
    ec_point_t P, Q, PmQ, A24;
    copy_point(&A24, &CURVE_E0_A24);
    copy_point(&P, &BASIS_EVEN.P);
    copy_point(&Q, &BASIS_EVEN.Q);
    copy_point(&PmQ, &BASIS_EVEN.PmQ);
}

int main() {

    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    printf("Running hd module unit tests\n");

    res = res & hd_chain_test();

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("\nAll tests passed!\n");
    }
    return(!res);
}
