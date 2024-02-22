#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "ec.h"
#include "biextension.h"
#include <endomorphism_action.h>

void fp2_exp_2e(int e, fp2_t* r, fp2_t const* x)
{
    fp2_copy(r, x);
    for (int i=0; i<e; i++) {
        fp2_sqr(r, r);
    }
}

int biextension_test() 
{
    ec_curve_t E0 = CURVE_E0; 
    int e = TORSION_PLUS_EVEN_POWER;
    ibz_t two_pow, tmp;
    fp2_t one, r1, rr1, r2, r3;
    ec_point_t P, Q, PmQ, A24;
    ec_point_t PQ, PP, QQ, PPQ, PQQ;

    ibz_init(&two_pow); ibz_init(&tmp);
    ibz_pow(&two_pow,&ibz_const_two,length);

    E0=CURVE_E0;
    A24=CURVE_E0_A24
    copy_point(&A24, &CURVE_E0_A24);
    copy_point(&P, &BASIS_EVEN.P);
    copy_point(&Q, &BASIS_EVEN.Q);
    copy_point(&GPmQ, &BASIS_EVEN.PmQ);

    xADD(&PQ, &P, &Q, &PmQ);
    weil(&r1, &P, &Q, &PQ, &A24);
    fp2_setone(&one);
    fp2_exp_2e(&r2, e, &r1);
    assert(fp2_is_equal(&r2, &one));

    xDBLv2(&PP, &P, &A24);
    xDBLv2(&QQ, &Q, &A24);
    xADD(&PPQ, &PQ, &P, &Q);
    xADD(&PQQ, &PQ, &Q, &P);

    weil(&r2, &PP, &Q, &PPQ, &A24);
    weil(&r3, &P, &QQ, &PQQ, &A24);
    assert(fp2_is_equal(&r2, &r3));
    fp2_sqr(&rr1, &r1);
    assert(fp2_is_equal(&rr1, &r2));

    return 0;
}

int main() 
{
    int res = 1;

    printf("Running biextension unit tests\n");

    res = res & biextension_test();

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("\nAll tests passed!\n");
    }
    return(!res);
}
