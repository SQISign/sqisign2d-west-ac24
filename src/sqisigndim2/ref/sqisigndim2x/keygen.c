#include <sqisigndim2.h>
#include <curve_extras.h>


#include <inttypes.h>





void secret_key_init(secret_key_t *sk) {
    quat_left_ideal_init(&(sk->secret_ideal));
    ibz_mat_2x2_init(&(sk->mat_BAcan_to_BA0_two));
}

void secret_key_finalize(secret_key_t *sk) {
    quat_left_ideal_finalize(&(sk->secret_ideal));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_two));
}


void protocols_keygen(public_key_t *pk, secret_key_t *sk) {
    
    int found=1;

    ibz_t n,n_temp;
    quat_alg_elem_t gen;
    quat_left_ideal_t lideal;
    ec_point_t list_points[3];
    ec_basis_t B_can_two,B_0_two;
    ibz_vec_4_t coeffs;

    quat_alg_elem_t beta1,beta2;
    ibz_t u,v,au,bu,d1,d2;

    // theta stuff
    theta_chain_t Phi;
    theta_chain_t phiu,phiv;


    ibz_init(&n);
    ibz_init(&n_temp);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal); 
    quat_alg_elem_init(&beta1);
    quat_alg_elem_init(&beta2);

    ibz_vec_4_init(&coeffs);

    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&d1);ibz_init(&d2);

    // generate a random ideal of random norm for the secret ideal
    // TODO make a clean function for all of that and  
    generate_random_prime(&n,1,128);
    generate_random_prime(&n_temp,1,256);
    ibz_mul(&n_temp,&n,&n_temp);
    found = found && represent_integer(&gen,&n_temp,&QUATALG_PINFTY);
    assert(found);
    quat_lideal_create_from_primitive(&lideal,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);


    // ideal to isogeny clapotis
    found = dim2id2iso_ideal_to_isogeny_clapotis(&Phi,&beta1,&beta2,&u,&v,&coeffs,&phiu,&phiv,&d1,&d2,&sk->curve,&B_0_two,&lideal,&QUATALG_PINFTY);
    assert(found);

    // B_O_two is equal to u * phi1 (BASIS_EVEN) where phi1 is an isogeny of degree d1
    // computation of sk->secret_ideal = lideal* conj(beta1) /norm(lideal)
    quat_alg_conj(&beta1,&beta1);
    ibz_mul(&beta1.denom,&beta1.denom,&lideal.norm); 
    found = quat_lideal_mul(&sk->secret_ideal,&lideal,&beta1,&QUATALG_PINFTY,0); 
    assert(found);


    assert(test_point_order_twof(&(B_0_two.P), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.Q), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.PmQ), &(sk->curve),TORSION_PLUS_EVEN_POWER));


    // sk->two_to_three_transporter = conj(gamma)/power_of_2
    // quat_alg_conj(&(sk->two_to_three_transporter), &gamma);
    // ibz_mul(&(sk->two_to_three_transporter.denom), &(sk->two_to_three_transporter.denom) , &(sk->secret_ideal_two.norm));


    copy_curve(&(pk->curve), &(sk->curve));


    ec_curve_to_basis_2(&B_can_two, &(sk->curve)); // canonical 
    assert(test_point_order_twof(&(B_can_two.P), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.Q), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.PmQ), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    
    // point_print("pk->curve->basis_two->P = ", B_can_two.P);
    // point_print("pk->curve->basis_two->Q = ", B_can_two.Q);
    // point_print("pk->curve->basis_two->PmQ = ", B_can_two.PmQ);
    printf("in of change matrix \n");
    change_of_basis_matrix_two(&(sk->mat_BAcan_to_BA0_two), &B_can_two, &B_0_two, &(sk->curve));
    
    printf("out of change matrix \n");

    // ec_curve_to_basis_3(&B_can_three, &(pk->curve)); // canonical 
    // assert(test_point_order_threef(&(B_can_three.P), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    // assert(test_point_order_threef(&(B_can_three.Q), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    // assert(test_point_order_threef(&(B_can_three.PmQ), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));

    // point_print("pk->curve->basis_three->P = ", B_can_three.P);
    // point_print("pk->curve->basis_three->Q = ", B_can_three.Q);
    // point_print("pk->curve->basis_three->PmQ = ", B_can_three.PmQ);
    // change_of_basis_matrix_three(&(sk->mat_BAcan_to_BA0_three), &B_can_three, &B_0_three, &(sk->curve));


    ibz_finalize(&n);
    ibz_finalize(&n_temp);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal); 

    ibz_vec_4_finalize(&coeffs);

    quat_alg_elem_finalize(&beta1);
    quat_alg_elem_finalize(&beta2);

    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&d1);ibz_finalize(&d2);
    return;
}

