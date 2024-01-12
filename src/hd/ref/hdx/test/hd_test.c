#include "hd_test.h"
#include <curve_extras.h>


// p = 13 * 2^126 * 3^78 − 1 
// 2^122 - 3 = 2274117029588772901^2 + 381187249096818150^2;

static int test_point_order_twof(const ec_point_t *P, const ec_curve_t *E) {
    ec_point_t test;
    copy_point(&test, P);
    if (fp2_is_zero(&test.z)) return 0;
    for (int i = 0;i<TORSION_PLUS_EVEN_POWER-1;i++) {
        ec_dbl(&test,E,&test);
    }
    if (fp2_is_zero(&test.z)) return 0;
    ec_dbl(&test,E,&test);
    return (fp2_is_zero(&test.z));
}

static int test_point_order_threef(const ec_point_t *P, const ec_curve_t *E) {
    ec_point_t test;
    copy_point(&test, P);
    digit_t three[NWORDS_ORDER] = {0};
    three[0] = 3;
    if (fp2_is_zero(&test.z)) return 0;
    for (int i = 0;i<TORSION_PLUS_ODD_POWERS[0]-1;i++) {
        ec_mul(&test, E, three, &test);
    }
    if (fp2_is_zero(&test.z)) return 0;
    ec_mul(&test, E, three, &test);
    return (fp2_is_zero(&test.z));
}



int hd_chain_test() {


    // var declaration & init
    ibz_t x,y,tmp,twopow;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&tmp);
    ibz_init(&twopow);
    ibz_pow(&twopow,&ibz_const_two,TORSION_PLUS_EVEN_POWER);

    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);


    digit_t scalars[2][NWORDS_ORDER] = {0};

    ec_curve_t E0,E1;
    ec_basis_t B0_two,B1_two,B0_three;
    ec_isog_odd_t phi;


    // setting the coefficient 
    ibz_set(&x,2274117029588772901);
    ibz_set(&y,381187249096818150);

    // copying the basis
    copy_point(&B0_two.P,&BASIS_EVEN.P);
    copy_point(&B0_two.Q,&BASIS_EVEN.Q);
    copy_point(&B0_two.PmQ,&BASIS_EVEN.PmQ);

    copy_point(&B0_three.P,&BASIS_THREE.P);
    copy_point(&B0_three.Q,&BASIS_THREE.Q);
    copy_point(&B0_three.PmQ,&BASIS_THREE.PmQ);

    E0=CURVE_E0;

    #ifndef NDEBUG
        assert(test_point_order_twof(&B0_two.P,&E0));
        assert(test_point_order_twof(&B0_two.Q,&E0));
        assert(test_point_order_twof(&B0_two.PmQ,&E0));
        assert(test_point_order_threef(&B0_three.P,&E0));        
    #endif

    ibz_set(&mat[0][0],0);ibz_set(&mat[0][1],0);ibz_set(&mat[1][0],0);ibz_set(&mat[1][1],0);

    // constructing the matrix corresponding to the endomorphism x + iy
    for (unsigned i = 0; i < 2; ++i) {
            ibz_add(&mat[i][i], &mat[i][i], &x);
            for (unsigned j = 0; j < 2; ++j) {
                ibz_mul(&tmp, &ACTION_GEN2[i][j], &y);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mod(&mat[i][j],&mat[i][j],&twopow);
        }
    }    

    // applying the matrix on the two torsion

    // first basis element
    ibz_to_digit_array(scalars[0],&mat[0][0]);
    // ibz_set(&mat[0][1],0);
    ibz_to_digit_array(scalars[1],&mat[1][0]);
    
    ec_biscalar_mul(&B1_two.P,&CURVE_E0,scalars[0],scalars[1],&B0_two);
    ibz_to_digit_array(scalars[0],&mat[0][1]);
    ibz_to_digit_array(scalars[1],&mat[1][1]);
    ec_biscalar_mul(&B1_two.Q,&CURVE_E0,scalars[0],scalars[1],&B0_two);

    assert(test_point_order_twof(&B1_two.P,&E0));
    assert(test_point_order_twof(&B1_two.Q,&E0));


    // setting up the isogeny 
    phi.curve = E0;
    uint8_t tab[3]= {1,0,0};
    phi.degree[0]=tab[0];
    phi.degree[1]=tab[1];
    phi.degree[2]=tab[2];
    printf("%d %d %d \n",phi.degree[0],phi.degree[1],phi.degree[2]);
    ec_set_zero(&phi.ker_minus);
    // preparating the kernel of phi as a point of 3 torsion
    digit_t three[NWORDS_ORDER] = {0};
    three[0] = 3;
    phi.ker_plus=B0_three.P;
    for (int i=0;i<TORSION_ODD_POWERS[0]-1;i++) {
        ec_mul(&phi.ker_plus,&E0,three,&phi.ker_plus);
    }
    ec_eval_odd_basis(&E1,&phi,&B0_two,1);

    assert(test_point_order_twof(&B0_two.P,&E1));
    assert(test_point_order_twof(&B0_two.Q,&E1));

    // ready to make the dim 2 computation
    theta_couple_curve_t E01;
    theta_couple_point_t T1;
    theta_couple_point_t T2;
    theta_chain_t dimtwo_chain;

    // setting the couples
    E01.E1=E0;
    E01.E2=E1;
    T1.P1=B1_two.P;
    T1.P2=B0_two.P;
    T2.P1=B1_two.Q;
    T2.P2=B0_two.Q;

    assert(test_point_order_twof(&T1.P1,&E01.E1));
    assert(test_point_order_twof(&T1.P2,&E01.E2));
    assert(test_point_order_twof(&T2.P1,&E01.E1));
    assert(test_point_order_twof(&T2.P2,&E01.E2));

    #ifndef NDEBUG
    theta_couple_point_t C1,C2;
        C1=T1;
        C2=T2;
        double_couple_point_iter(&C1,TORSION_PLUS_EVEN_POWER,&E01,&C1);
        double_couple_point_iter(&C2,TORSION_PLUS_EVEN_POWER,&E01,&C2);
        assert(fp2_is_zero(&C1.P1.z));
        assert(fp2_is_zero(&C1.P2.z));
        assert(fp2_is_zero(&C2.P1.z));
        assert(fp2_is_zero(&C2.P2.z));
    #endif 

    // multiplying by 2
    double_couple_point_iter(&T1,2,&E01,&T1);
    double_couple_point_iter(&T2,2,&E01,&T2);

    int length=122;

    #ifndef NDEBUG

        // checking that the points have the correct order 
        theta_couple_point_t P1,P2;
        P1=T1;
        P2=T2;
        double_couple_point_iter(&P1,length+1,&E01,&P1);
        double_couple_point_iter(&P2,length+1,&E01,&P2);

        ec_point_t test1,test2;
        test1=P1.P1;
        test2=P1.P2;
        assert(!fp2_is_zero(&test1.z));
        assert(!fp2_is_zero(&test2.z));
        ec_dbl(&test1,&E01.E1,&test1);
        ec_dbl(&test2,&E01.E2,&test2);
        assert(fp2_is_zero(&test1.z));
        assert(fp2_is_zero(&test2.z));
        test1=P2.P1;
        test2=P2.P2;
        assert(!fp2_is_zero(&test1.z));
        assert(!fp2_is_zero(&test2.z));
        ec_dbl(&test1,&E01.E1,&test1);
        ec_dbl(&test2,&E01.E2,&test2);
        assert(fp2_is_zero(&test1.z));
        assert(fp2_is_zero(&test2.z));
    #endif

    theta_chain_comput(&dimtwo_chain,length,&E01,&T1,&T2);


    return 1;
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
        printf("All tests passed!\n");
    }
    return(!res);
}