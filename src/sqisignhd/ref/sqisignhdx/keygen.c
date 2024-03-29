#include <sqisignhd.h>
#include <curve_extras.h>


#include <inttypes.h>


void secret_key_init(secret_key_t *sk) {
    quat_left_ideal_init(&(sk->secret_ideal_two));
    quat_alg_elem_init(&(sk->two_to_three_transporter));
    ibz_mat_2x2_init(&(sk->mat_BAcan_to_BA0_two));
    ibz_mat_2x2_init(&(sk->mat_BAcan_to_BA0_three));
}

void secret_key_finalize(secret_key_t *sk) {
    quat_left_ideal_finalize(&(sk->secret_ideal_two));
    quat_alg_elem_finalize(&(sk->two_to_three_transporter));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_two));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_three));
}


void protocols_keygen(public_key_t *pk, secret_key_t *sk) {
    quat_left_ideal_t lideal_odd;
    quat_alg_elem_t gamma;
    ec_isog_even_t two_isogeny_first_half, two_isogeny_second_half;
    ec_isog_odd_t phi_first_half, phi_second_half;
    ec_point_t list_points[3];
    ec_basis_t B_can_two, B_can_three, B_0_two, B_0_three;

    quat_left_ideal_init(&lideal_odd); 
    quat_alg_elem_init(&gamma);

    doublepath(&gamma, &(sk->secret_ideal_two), &lideal_odd, 
    &B_0_three, 
    &B_0_two, &(sk->curve),
    0);



    assert(test_point_order_twof(&(B_0_two.P), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.Q), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.PmQ), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    
    assert(test_point_order_threef(&(B_0_three.P), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    assert(test_point_order_threef(&(B_0_three.Q), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    assert(test_point_order_threef(&(B_0_three.PmQ), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));


    // sk->two_to_three_transporter = conj(gamma)/power_of_2
    quat_alg_conj(&(sk->two_to_three_transporter), &gamma);
    ibz_mul(&(sk->two_to_three_transporter.denom), &(sk->two_to_three_transporter.denom) , &(sk->secret_ideal_two.norm));



    #ifndef NDEBUG
    {
        quat_left_ideal_t lideal_odd_again;
        quat_left_ideal_init(&lideal_odd_again); 

        quat_lideal_mul(&lideal_odd_again, &(sk->secret_ideal_two), &(sk->two_to_three_transporter), &QUATALG_PINFTY, 0); 

        assert(quat_lideal_equals(&lideal_odd_again, &lideal_odd, &QUATALG_PINFTY));

        quat_left_ideal_finalize(&lideal_odd_again); 
    }
    #endif


    copy_curve(&(pk->curve), &(sk->curve));


    ec_curve_to_basis_2(&B_can_two, &(pk->curve),TORSION_PLUS_EVEN_POWER); // canonical 
    assert(test_point_order_twof(&(B_can_two.P), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.Q), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.PmQ), &(sk->curve),TORSION_PLUS_EVEN_POWER));
    
    // point_print("pk->curve->basis_two->P = ", B_can_two.P);
    // point_print("pk->curve->basis_two->Q = ", B_can_two.Q);
    // point_print("pk->curve->basis_two->PmQ = ", B_can_two.PmQ);
    change_of_basis_matrix_two(&(sk->mat_BAcan_to_BA0_two), &B_can_two, &B_0_two, &(sk->curve),TORSION_PLUS_EVEN_POWER);
    

    ec_curve_to_basis_3(&B_can_three, &(pk->curve)); // canonical 
    assert(test_point_order_threef(&(B_can_three.P), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    assert(test_point_order_threef(&(B_can_three.Q), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));
    assert(test_point_order_threef(&(B_can_three.PmQ), &(sk->curve),TORSION_PLUS_ODD_POWERS[0]));

    // point_print("pk->curve->basis_three->P = ", B_can_three.P);
    // point_print("pk->curve->basis_three->Q = ", B_can_three.Q);
    // point_print("pk->curve->basis_three->PmQ = ", B_can_three.PmQ);
    change_of_basis_matrix_three(&(sk->mat_BAcan_to_BA0_three), &B_can_three, &B_0_three, &(sk->curve));


    quat_left_ideal_finalize(&lideal_odd); 
    quat_alg_elem_finalize(&gamma);
    return;
}

