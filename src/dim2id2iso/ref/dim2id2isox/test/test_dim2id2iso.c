// #include <inttypes.h>

#include "dim2id2iso_tests.h"

int dim2id2iso_test_fixed_degree_isogeny() {

    ibz_t u;
    quat_left_ideal_t lideal;
    ibz_init(&u);
    quat_left_ideal_init(&lideal);

    // u is 3 times a random prime 
    generate_random_prime(&u,1,TORSION_PLUS_EVEN_POWER-10);
    ibz_mul(&u,&u,&ibz_const_three);

    theta_chain_t F;
    fixed_degree_isogeny(&F,&lideal,&u,1);

    ibz_finalize(&u);
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
    theta_chain_t phiv;
 
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

    found = dim2id2iso_ideal_to_isogeny_clapotis(&Phi,&beta1,&beta2,&u,&v,&coeffs,&phiv,&d1,&d2,&lideal_small,&QUATALG_PINFTY);

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

    printf("Running dim2id2iso module unit tests\n");

    res = res & dim2id2iso_test_fixed_degree_isogeny();

    res = res & dim2id2iso_test_dimid2iso();

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("\nAll tests passed!\n");
    }
    return(!res);
}