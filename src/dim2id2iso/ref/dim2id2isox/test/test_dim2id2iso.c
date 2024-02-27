#include <inttypes.h>
#include "dim2id2iso_tests.h"
#include <tools.h>


static void fp2_print(char *name, fp2_t const a){
    fp_t b1,b2;
    fp_frommont(b1,a.re);
    fp_frommont(b2,a.im);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b1[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b2[i]);
    printf("\n");
}

int dim2id2iso_test_fixed_degree_isogeny() {

    ibz_t u,two_pow;
    quat_left_ideal_t lideal;
    ibz_init(&u);
    ibz_init(&two_pow);   
    quat_left_ideal_init(&lideal);

    // u is 3 times a random prime 
    generate_random_prime(&u,1,TORSION_PLUS_EVEN_POWER-10);
    ibz_mul(&u,&u,&ibz_const_three);

    theta_chain_t F;
    fixed_degree_isogeny(&F,&lideal,&u,1);

    // now we check that we get a correct codomain in the end
    // by evaluating some point and checking the pairing
    theta_couple_jac_point_t Teval1,Teval2,Teval3;
    jac_point_t temp;
    theta_couple_point_t Tev1,Tev2,Tev3;
    theta_couple_curve_t E00;
    ec_curve_t E0;
    E0 = CURVE_E0;
    E00.E1 = CURVE_E0;
    E00.E2 = CURVE_E0;
    ec_basis_t bas = BASIS_EVEN;
    lift_basis(&Teval1.P1,&Teval2.P1,&bas,&E0);
    jac_neg(&temp,&Teval2.P1);
    ADD(&Teval3.P1,&Teval1.P1,&temp,&E0);
    fp2_set(&Teval1.P2.x,0);
    fp2_set(&Teval1.P2.y,1);
    fp2_set(&Teval1.P2.z,0);
    fp2_set(&Teval2.P2.x,0);
    fp2_set(&Teval2.P2.y,1);
    fp2_set(&Teval2.P2.z,0);
    fp2_set(&Teval3.P2.x,0);
    fp2_set(&Teval3.P2.y,1);
    fp2_set(&Teval3.P2.z,0);
    theta_chain_eval_no_help(&Tev1,&F,&Teval1,&E00);
    theta_chain_eval_no_help(&Tev2,&F,&Teval2,&E00);
    theta_chain_eval_no_help(&Tev3,&F,&Teval3,&E00);



    assert(test_point_order_twof(&Tev1.P1,&F.codomain.E1,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Tev1.P2,&F.codomain.E2,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Tev2.P1,&F.codomain.E1,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Tev2.P2,&F.codomain.E2,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Tev3.P1,&F.codomain.E1,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Tev3.P2,&F.codomain.E2,TORSION_PLUS_EVEN_POWER));

    fp2_t w0,w1,w2;
    ec_point_t AC,A24;
    copy_point(&AC, &CURVE_E0_A24); //Warning, this is AC, not A24!
    A24_from_AC(&A24, &AC);
    weil(&w0,TORSION_PLUS_EVEN_POWER,&bas.P,&bas.Q,&bas.PmQ,&A24);
    // Changing the AC
    fp2_copy(&AC.x,&F.codomain.E1.A);
    fp2_copy(&AC.z,&F.codomain.E1.C);
    A24_from_AC(&A24, &AC);
    weil(&w1,TORSION_PLUS_EVEN_POWER,&Tev1.P1,&Tev2.P1,&Tev3.P1,&A24);
    fp2_copy(&AC.x,&F.codomain.E2.A);
    fp2_copy(&AC.z,&F.codomain.E2.C);
    A24_from_AC(&A24, &AC);
    weil(&w2,TORSION_PLUS_EVEN_POWER,&Tev1.P2,&Tev2.P2,&Tev3.P2,&A24);
    ibz_pow(&two_pow,&ibz_const_two,F.length);
    ibz_sub(&two_pow,&two_pow,&u);
    // now we are checking that one of the two is equal to the correct value 
    digit_t digit_u[NWORDS_ORDER]={0};
    ibz_to_digit_array(digit_u,&u);
    fp2_t test_pow;
    fp2_pow(&test_pow,&w0,digit_u,NWORDS_ORDER);

    // it seems like we always get the second curve 
    assert(fp2_is_equal(&test_pow,&w2));
    ibz_to_digit_array(digit_u,&two_pow);
    fp2_pow(&test_pow,&w0,digit_u,NWORDS_ORDER);
    assert(fp2_is_equal(&test_pow,&w1));

    if (fp2_is_equal(&test_pow,&w1)) {
        // printf("first curve!  \n");
        ibz_to_digit_array(digit_u,&two_pow);
        fp2_pow(&test_pow,&w0,digit_u,NWORDS_ORDER);
        fp2_is_equal(&test_pow,&w2);
    }
    else if (fp2_is_equal(&test_pow,&w2)) {
        // printf(" \n \n \n \n second curve! \n \n \n \n");
        assert(fp2_is_equal(&test_pow,&w2));
        ibz_to_digit_array(digit_u,&two_pow);
        fp2_pow(&test_pow,&w0,digit_u,NWORDS_ORDER);
        assert(fp2_is_equal(&test_pow,&w1));
    }
    else {
        assert(0);
    }
    

    ibz_finalize(&u);
    ibz_finalize(&two_pow);
    quat_left_ideal_finalize(&lideal);
    return 1;

}


int dim2id2iso_test_dimid2iso() {

    // var dec
    int found =1;  
    ibz_t temp,remainder,n1,n2;
    ibq_t ibq_norm;
    quat_alg_elem_t gen;

    quat_left_ideal_t lideal_small;
    quat_order_t right_order;
    ibz_mat_4x4_t reduced,gram;
    ibz_vec_4_t coeffs;
    quat_alg_elem_t beta1,beta2;
    ibz_t u,v,au,bu,d1,d2;

    // theta stuff
    theta_chain_t Phi;
    theta_chain_t phiu,phiv;
 
    // var init
    ibq_init(&ibq_norm);
    ibz_init(&temp);ibz_init(&remainder);ibz_init(&n1);ibz_init(&n2);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal_small);
    quat_order_init(&right_order);
    ibz_mat_4x4_init(&reduced);ibz_mat_4x4_init(&gram);
    ibz_vec_4_init(&coeffs);

    quat_alg_elem_init(&beta1);
    quat_alg_elem_init(&beta2);

    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&d1);ibz_init(&d2);

    // computation of lideal_small
    generate_random_prime(&n1,1,128);
    generate_random_prime(&n2,1,256);
    ibz_mul(&temp,&n1,&n2);
    found = found && represent_integer(&gen,&temp,&QUATALG_PINFTY);
    assert(found);
    quat_lideal_create_from_primitive(&lideal_small,&gen,&n1,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);

    clock_t tt = tic();
    found = dim2id2iso_ideal_to_isogeny_clapotis(&Phi,&beta1,&beta2,&u,&v,&coeffs,&phiu,&phiv,&d1,&d2,&lideal_small,&QUATALG_PINFTY);
    TOC(tt,"total time for ideal to isogeny clapotis \n");

    // // now we verify that Phi correctly translates the isogeny we are looking for by evaluating and checking weil pairing
    // theta_couple_jac_point_t Teval1,Teval2,Teval3;
    // jac_point_t temp_neg;
    // theta_couple_point_t Tev1,Tev2,Tev3;
    // theta_couple_curve_t E00;
    // ec_curve_t E0;
    // E0 = CURVE_E0;
    // E00.E1 = CURVE_E0;
    // E00.E2 = CURVE_E0;
    // ec_basis_t bas = BASIS_EVEN;
    // lift_basis(&Teval1.P1,&Teval2.P1,&bas,&E0);
    // jac_neg(&temp_neg,&Teval2.P1);
    // ADD(&Teval3.P1,&Teval1.P1,&temp_neg,&E0);
    // fp2_set(&Teval1.P2.x,0);
    // fp2_set(&Teval1.P2.y,1);
    // fp2_set(&Teval1.P2.z,0);
    // fp2_set(&Teval2.P2.x,0);
    // fp2_set(&Teval2.P2.y,1);
    // fp2_set(&Teval2.P2.z,0);
    // fp2_set(&Teval3.P2.x,0);
    // fp2_set(&Teval3.P2.y,1);
    // fp2_set(&Teval3.P2.z,0);
    // theta_chain_eval_no_help(&Tev1,&Phi,&Teval1,&E00);
    // theta_chain_eval_no_help(&Tev2,&Phi,&Teval2,&E00);
    // theta_chain_eval_no_help(&Tev3,&Phi,&Teval3,&E00);

    // assert(test_point_order_twof(&Tev1.P1,&Phi.codomain.E1,TORSION_PLUS_EVEN_POWER));
    // assert(test_point_order_twof(&Tev1.P2,&Phi.codomain.E2,TORSION_PLUS_EVEN_POWER));
    // assert(test_point_order_twof(&Tev2.P1,&Phi.codomain.E1,TORSION_PLUS_EVEN_POWER));
    // assert(test_point_order_twof(&Tev2.P2,&Phi.codomain.E2,TORSION_PLUS_EVEN_POWER));
    // assert(test_point_order_twof(&Tev3.P1,&Phi.codomain.E1,TORSION_PLUS_EVEN_POWER));
    // assert(test_point_order_twof(&Tev3.P2,&Phi.codomain.E2,TORSION_PLUS_EVEN_POWER));

    // fp2_t w0,w1,w2;
    // ec_point_t AC,A24;
    // copy_point(&AC, &CURVE_E0_A24); //Warning, this is AC, not A24!
    // A24_from_AC(&A24, &AC);
    // weil(&w0,TORSION_PLUS_EVEN_POWER,&bas.P,&bas.Q,&bas.PmQ,&A24);
    // // Changing the AC
    // fp2_copy(&AC.x,&Phi.codomain.E1.A);
    // fp2_copy(&AC.z,&Phi.codomain.E1.C);
    // A24_from_AC(&A24, &AC);
    // weil(&w1,TORSION_PLUS_EVEN_POWER,&Tev1.P1,&Tev2.P1,&Tev3.P1,&A24);
    // fp2_copy(&AC.x,&Phi.codomain.E2.A);
    // fp2_copy(&AC.z,&Phi.codomain.E2.C);
    // A24_from_AC(&A24, &AC);
    // weil(&w2,TORSION_PLUS_EVEN_POWER,&Tev1.P2,&Tev2.P2,&Tev3.P2,&A24);
    // // now we are checking that one of the two is equal to the correct value 
    // digit_t digit_d[NWORDS_ORDER]={0};
    // ibz_to_digit_array(digit_d,&d1);
    // fp2_t test_pow;
    // fp2_pow(&test_pow,&w0,digit_d,NWORDS_ORDER);

    // // it seems like we always get the second curve 
    // // assert(fp2_is_equal(&test_pow,&w2));
    // // ibz_to_digit_array(digit_u,&two_pow);
    // // fp2_pow(&test_pow,&w0,digit_u,NWORDS_ORDER);
    // // assert(fp2_is_equal(&test_pow,&w1));

    // if (fp2_is_equal(&test_pow,&w1)) {
    //     printf("first curve encodes degree d1! \n");
    //     ibz_to_digit_array(digit_d,&d2);
    //     fp2_pow(&test_pow,&w0,digit_d,NWORDS_ORDER);
    //     assert(fp2_is_equal(&test_pow,&w2));
    // }
    // else if (fp2_is_equal(&test_pow,&w2)) {
    //     printf(" second curve encodes degree d1! \n");
    //     assert(fp2_is_equal(&test_pow,&w2));
    //     ibz_to_digit_array(digit_d,&d2);
    //     fp2_pow(&test_pow,&w0,digit_d,NWORDS_ORDER);
    //     assert(fp2_is_equal(&test_pow,&w1));
    // }
    // else {
    //     assert(0);
    // }



    ibq_finalize(&ibq_norm);
    ibz_finalize(&temp);ibz_finalize(&remainder);ibz_finalize(&n1);ibz_finalize(&n2);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal_small);
    quat_order_finalize(&right_order);
    ibz_mat_4x4_finalize(&reduced);ibz_mat_4x4_finalize(&gram);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&d1);ibz_finalize(&d2);
    ibz_vec_4_finalize(&coeffs);

    quat_alg_elem_finalize(&beta1);
    quat_alg_elem_finalize(&beta2);

    return found;
}

int main() {

    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    printf("\nRunning dim2id2iso module unit tests\n");

    printf("\nRunning fixed degree tests\n");

    res = res & dim2id2iso_test_fixed_degree_isogeny();

    printf("\nRunning id2iso_clapotis tests\n");

    res = res & dim2id2iso_test_dimid2iso();

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("\nAll tests passed!\n");
    }
    return(!res);
}