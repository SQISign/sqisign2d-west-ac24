#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "biextension.h"
#include "curve_extras.h"
#include <endomorphism_action.h>
#include <torsion_constants.h>

void fp2_exp_2e(fp2_t* r, uint64_t e, fp2_t const* x)
{
    fp2_copy(r, x);
    for (uint64_t i=0; i<e; i++) {
        fp2_sqr(r, r);
    }
}

void dbl_2e(ec_point_t* R, uint64_t e, ec_point_t const* P, ec_point_t const* A24)
{
    copy_point(R, P);
    for (uint64_t i=0; i<e; i++) {
        xDBLv2(R, R, A24);
    }
}

void cubical_2e(ec_point_t* R, uint64_t e, ec_point_t const* P, ec_point_t const* A24)
{
    copy_point(R, P);
    for (uint64_t i=0; i<e; i++) {
        cubicalDBL(R, R, A24);
    }
}

int biextension_test() 
{
    ec_curve_t E0 = CURVE_E0; 
    uint64_t e = TORSION_PLUS_EVEN_POWER;
    // ibz_t two_pow, tmp;
    fp2_t one, r1, rr1, r2, r3;
    ec_point_t P, Q, PmQ, A24, AC;
    ec_point_t PQ, PP, QQ, PPQ, PQQ;
    ec_point_t Pe, Qe, PmQe;

    // ibz_init(&two_pow); ibz_init(&tmp);
    // ibz_pow(&two_pow,&ibz_const_two,length);

    copy_point(&AC, &CURVE_E0_A24); //Warning, this is AC, not A24!
    A24_from_AC(&A24, &AC);
    copy_point(&P, &BASIS_EVEN.P);
    copy_point(&Q, &BASIS_EVEN.Q);
    copy_point(&PmQ, &BASIS_EVEN.PmQ);

    printf("Testing order of points\n");
    dbl_2e(&Pe, e, &P, &A24);
    dbl_2e(&Qe, e, &Q, &A24);
    dbl_2e(&PmQe, e, &PmQ, &A24);
    assert(ec_is_zero(&Pe));
    assert(ec_is_zero(&Qe));
    assert(ec_is_zero(&PmQe));

    for (uint64_t i=0; i<e; i++) {
      dbl_2e(&Pe, i, &P, &A24);
      cubical_2e(&Qe, i, &P, &A24);
      printf("i=%d: r=%d\n", i, is_point_equal(&Pe, &Qe));
    }
    cubical_2e(&Pe, e, &P, &A24);
    assert(ec_is_zero(&Pe));

    printf("Testing order of Weil pairing\n");
    xADD(&PQ, &P, &Q, &PmQ);
    weil(&r1, e, &P, &Q, &PQ, &A24);
    fp2_setone(&one);
    fp2_exp_2e(&r2, e, &r1);
    assert(fp2_is_equal(&r2, &one));

    printf("Bilinearity tests\n");
    xDBLv2(&PP, &P, &A24);
    xDBLv2(&QQ, &Q, &A24);
    xADD(&PPQ, &PQ, &P, &Q);
    xADD(&PQQ, &PQ, &Q, &P);

    weil(&r2, e, &PP, &Q, &PPQ, &A24);
    weil(&r3, e, &P, &QQ, &PQQ, &A24);
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
