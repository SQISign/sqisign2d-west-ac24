#include "theta_isogenies.h"




// TODO check that we computed the matrix the right way (and not the transpose)
void get_matrix(fp2_t *a00,fp2_t *a01,fp2_t *a10,fp2_t *a11,const ec_point_t *P, const ec_curve_t *E) {
    ec_point_t P2;
    fp2_t det,x,u,invz,temp;
    
    // temp = [2]P
    ec_dbl(&P2,E,P);

    // invz = 1/P.z  
    invz = P->z;
    fp2_inv(&invz);
    fp2_inv(&x);

    // P = (x,z) P2 = (u,w)
    // det = x w - u z
    fp2_mul(&det,&P->x,&P2.z);
    fp2_mul(&temp,&P->z,&P2.x);
    fp2_sub(&det,&det,&temp);

    // det = 1/det 
    fp2_inv(&det);

    // a10 = ux /det - x/z
    fp2_mul(&temp,&P->x,&invz);
    fp2_mul( a10,&P->x,&P2.x);
    fp2_mul( a10,a10,&det);
    fp2_sub(a10,a10,&temp); 

    // a11 = u z * det
    fp2_mul(a11,&P2.x,&det);
    fp2_mul(a11,a11,&P->z);

    // a00 = -a11 
    fp2_neg(a00,a11);

    // a01 = - w z det 
    fp2_mul(a01,&P2.z,&det);
    fp2_mul(a01,a01,&P->z);
    fp2_neg(a01,a01);
}


/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing 
 * @param K1_8 a couple point
 * @param L1_8 a point in E1[8]
 * @param E2 an elliptic curve
 * @param K2_8 a point in E2[8]
 * @param L2_8 a point in E2[8]  
 * 
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8) 
 *  
   */
void gluing_comput(theta_gluing_t *out,const theta_couple_curve_t *E12,const theta_couple_point_t *K1_8,const theta_couple_point_t *K2_8)
{

    // var init
    fp2_t M1100;fp2_t M1101;fp2_t M1110;fp2_t M1111;
    fp2_t M1200;fp2_t M1201;fp2_t M1210;fp2_t M1211;
    fp2_t M2100;fp2_t M2101;fp2_t M2110;fp2_t M2111;
    fp2_t M2200;fp2_t M2201;fp2_t M2210;fp2_t M2211;
    fp2_t t001,t101,t002,t102,temp;
    

    // K1_4 = [2] K1_8  and K2_4 = [2] K2_8
    double_couple_point(&out->K1_4,E12,K1_8);
    double_couple_point(&out->K2_4,E12,K2_8);

    // computing the base change matrix
    // we start by computing the matrices for each points in the kernel
    get_matrix(&M1100,&M1101,&M1110,&M1111,&out->K1_4.P1,&E12->E1);
    get_matrix(&M1200,&M1201,&M1210,&M1211,&out->K1_4.P2,&E12->E2);
    get_matrix(&M2100,&M2101,&M2110,&M2111,&out->K2_4.P1,&E12->E1);
    get_matrix(&M2200,&M2201,&M2210,&M2211,&out->K2_4.P2,&E12->E2);

    // multiplication of the matrices
    // t001,t101 (resp t002,t102) first column of M11 * M21 (resp M12 * M22)
    fp2_mul(&t001,&M1100,&M2100);
    fp2_mul(&temp,&M1101,&M2110); 
    fp2_mul(&t001,&t001,&temp);  

    fp2_mul(&t101,&M1110,&M2100);
    fp2_mul(&temp,&M1111,&M2110); 
    fp2_mul(&t101,&t101,&temp);

    fp2_mul(&t002,&M1200,&M2200);
    fp2_mul(&temp,&M1201,&M2210); 
    fp2_mul(&t002,&t002,&temp);  

    fp2_mul(&t102,&M1210,&M2200);
    fp2_mul(&temp,&M1211,&M2210); 
    fp2_mul(&t102,&t102,&temp); 

    // trace for the first row 
    fp2_set(&out->M00,1);
    fp2_mul(&temp,&t001,&t002);
    fp2_add(&out->M00,&out->M00,&temp);
    fp2_mul(&temp,&M2100,&M2200);
    fp2_add(&out->M00,&out->M00,&temp);
    fp2_mul(&temp,&M1100,&M1200);
    fp2_add(&out->M00,&out->M00,&temp);

    fp2_mul(&out->M01,&t001,&t102);
    fp2_mul(&temp,&M2100,&M2210);
    fp2_add(&out->M01,&out->M01,&temp);
    fp2_mul(&temp,&M1100,&M1210);
    fp2_add(&out->M01,&out->M01,&temp);


    fp2_mul(&out->M02,&t101,&t002);
    fp2_mul(&temp,&M2110,&M2200);
    fp2_add(&out->M02,&out->M02,&temp);
    fp2_mul(&temp,&M1110,&M1200);
    fp2_add(&out->M02,&out->M02,&temp);

    fp2_mul(&out->M03,&t101,&t102);
    fp2_mul(&temp,&M2110,&M2210);
    fp2_add(&out->M03,&out->M03,&temp);
    fp2_mul(&temp,&M1110,&M1210);
    fp2_add(&out->M03,&out->M03,&temp);

    // Compute the action of (0,out.K2_4.P2) for the second row 
    fp2_mul(&temp,&M2201,&out->M01);
    fp2_mul(&out->M10,&M2200,&out->M00);
    fp2_add(&out->M10,&out->M10,&temp);

    fp2_mul(&temp,&M2211,&out->M01);
    fp2_mul(&out->M11,&M2210,&out->M00);
    fp2_add(&out->M11,&out->M11,&temp);

    fp2_mul(&temp,&M2201,&out->M03);
    fp2_mul(&out->M12,&M2200,&out->M02);
    fp2_add(&out->M12,&out->M12,&temp);

    fp2_mul(&temp,&M2211,&out->M03);
    fp2_mul(&out->M13,&M2210,&out->M02);
    fp2_add(&out->M13,&out->M13,&temp);

    // compute the action of (K1_4.P1,0) for the third row
    fp2_mul(&temp,&M1101,&out->M02);
    fp2_mul(&out->M20,&M1100,&out->M00);
    fp2_add(&out->M20,&out->M20,&temp);

    fp2_mul(&temp,&M1101,&out->M03);
    fp2_mul(&out->M21,&M1100,&out->M01);
    fp2_add(&out->M11,&out->M21,&temp);

    fp2_mul(&temp,&M1111,&out->M02);
    fp2_mul(&out->M22,&M1110,&out->M00);
    fp2_add(&out->M22,&out->M22,&temp);

    fp2_mul(&temp,&M1111,&out->M03);
    fp2_mul(&out->M23,&M1110,&out->M01);
    fp2_add(&out->M23,&out->M23,&temp);

    // compute the action of (K1_4.P1,K2_4.P2) for the final row
    fp2_mul(&temp,&M1101,&out->M12);
    fp2_mul(&out->M30,&M1100,&out->M10);
    fp2_add(&out->M30,&out->M30,&temp);

    fp2_mul(&temp,&M1101,&out->M13);
    fp2_mul(&out->M31,&M1100,&out->M11);
    fp2_add(&out->M31,&out->M31,&temp);

    fp2_mul(&temp,&M1111,&out->M12);
    fp2_mul(&out->M32,&M1110,&out->M10);
    fp2_add(&out->M32,&out->M32,&temp);

    fp2_mul(&temp,&M1111,&out->M13);
    fp2_mul(&out->M33,&M1110,&out->M11);
    fp2_add(&out->M33,&out->M33,&temp);

    

}