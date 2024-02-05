#include "theta_structure.h"
#include <assert.h>





/**
 * @brief Perform the hadamard transform on a theta point
 *
 * @param out Output: the theta_point 
 * @param in a theta point*  
 * in = (x,y,z,t)
 * out = (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)
 *  
   */
void hadamard(theta_point_t *out, const theta_point_t *in) {
    fp2_t t1,t2,t3,t4;
    
    // t1 = x + y
    fp2_add(&t1,&in->x,&in->y);
    // t2 = x-y
    fp2_sub(&t2,&in->x,&in->y);
    // t3 = z+t
    fp2_add(&t3,&in->z,&in->t);
    // t4 = z-t
    fp2_sub(&t4,&in->z,&in->t);

    fp2_add(&out->x,&t1,&t3);
    fp2_add(&out->y,&t2,&t4);
    fp2_sub(&out->z,&t1,&t3);
    fp2_sub(&out->t,&t2,&t4);
}


/**
 * @brief Square the coordinates and then perform the hadamard transform
 *
 * @param out Output: the theta_point 
 * @param in a theta point*  
 * in = (x,y,z,t)
 * out = (x^2+y^2+z^2+t^2, x^2-y^2+z^2-t^2, x^2+y^2-z^2-t^2, x^2-y^2-z^2+t^2)
 *  
   */
void to_squared_theta(theta_point_t *out,const theta_point_t *in) {
    fp2_sqr(&out->x,&in->x);
    fp2_sqr(&out->y,&in->y);
    fp2_sqr(&out->z,&in->z);
    fp2_sqr(&out->t,&in->t);
    hadamard(out,out);
}


/**
 * @brief Perform the theta structure precomputation 
 *
 * @param A Output: the theta_structure  
 * 
 * if A.null_point = (x,y,z,t)
 * if (xx,yy,zz,tt) = to_squared_theta(A.null_point)
 * Computes y0,z0,t0,Y0,Z0,T0 = x/y,x/z,x/t,XX/YY,XX/ZZ,XX/TT 
 *  
   */
void theta_precomputation(theta_structure_t *A){

    if (!A->precomputation) {
         // temp = (xx,yy,zz,tt)
    theta_point_t temp;
    to_squared_theta(&temp,&A->null_point);
    // Computes t1,t2,t3,t4,t5,t6 = 1/y,1/z,1/t,1/YY,1/ZZ,1/TT
    // TODO : batch inversion? 
    fp2_t t1,t2,t3,t4,t5,t6;
    t1 = A->null_point.y;
    t2 = A->null_point.z;
    t3 = A->null_point.t;
    t4 = temp.y;
    t5 = temp.z;
    t6 = temp.t;
    fp2_inv(&t1);
    fp2_inv(&t2);
    fp2_inv(&t3);
    fp2_inv(&t4);
    fp2_inv(&t5);
    fp2_inv(&t6);

    // y0,z0,t0,Y0,Z0,T0 = x/y,x/z,x/t,XX/YY,XX/ZZ,XX/TT 
    fp2_mul(&A->y0,&t1,&A->null_point.x);
    fp2_mul(&A->z0,&t2,&A->null_point.x);
    fp2_mul(&A->t0,&t3,&A->null_point.x);
    fp2_mul(&A->Y0,&t4,&temp.x);
    fp2_mul(&A->Z0,&t5,&temp.x);
    fp2_mul(&A->T0,&t6,&temp.x);

    A->precomputation=1;
    }
   
}

/**
 * @brief Compute the double of the theta point in on the theta struc A
 *
 * @param out Output: the theta_point 
 * @param A a theta structure
 * @param in a theta point in the theta structure A   
 * in = (x,y,z,t)
 * out = [2] (x,y,z,t)
 * /!\ assumes that no coordinates is zero and that A has been precomputed
 *  
   */
void double_point(theta_point_t *out,theta_structure_t *A,const theta_point_t *in) {

    to_squared_theta(out,in);
    fp2_sqr(&out->x,&out->x);
    fp2_sqr(&out->y,&out->y);
    fp2_sqr(&out->z,&out->z);
    fp2_sqr(&out->t,&out->t);
    if (!A->precomputation) {
        theta_precomputation(A);
    }
    fp2_mul(&out->y,&out->y,&A->Y0);
    fp2_mul(&out->z,&out->z,&A->Z0);
    fp2_mul(&out->t,&out->t,&A->T0);


    hadamard(out,out);
    fp2_mul(&out->y,&out->y,&A->y0);
    fp2_mul(&out->z,&out->z,&A->z0);
    fp2_mul(&out->t,&out->t,&A->t0);

}

/**
 * @brief Compute the iterated double of the theta point in on the theta struc A
 *
 * @param out Output: the theta_point 
 * @param A a theta structure
 * @param in a theta point in the theta structure A  
 * @param exp the exponent
 * in = (x,y,z,t)
 * out = [2^2] (x,y,z,t)
 * /!\ assumes that no coordinates is zero and that the precomputation of A has been done
 *  
   */
void double_iter(theta_point_t *out,theta_structure_t *A,const theta_point_t *in, int exp){
    if (exp==0) {
        fp2_copy(&out->x,&in->x);
        fp2_copy(&out->y,&in->y);
        fp2_copy(&out->z,&in->z);
        fp2_copy(&out->t,&in->t);
    }
    else {
        double_point(out,A,in);
        for (int i=1;i<exp;i++){
            double_point(out,A,out);
        }
    }
    
}

/**
 * @brief Compute the differential addition of two theta points in the theta struc A
 *
 * @param R Output: the theta_point 
 * @param A a theta structure
 * @param P a theta point in the theta structure A
 * @param Q a theta point in the theta structure A
 * @param PQ the theta point P-Q in the theta structure A      
 * R = P+Q
 * /!\ assumes that no coordinates is zero and that the precomputation of A has been done
 *  
   */
void diff_add(theta_point_t *R,const theta_structure_t *A,const theta_point_t *P,const theta_point_t *Q,const theta_point_t *PQ) {

    theta_point_t tP,tQ;
    // tP = H(P), tQ = H(Q)
    hadamard(&tP,P);
    hadamard(&tQ,Q);

    fp2_mul(&R->x,&tP.x,&tQ.x);
    fp2_mul(&R->y,&tP.y,&tQ.y);
    fp2_mul(&R->y,&R->y,&A->Y0);
    fp2_mul(&R->z,&tP.z,&tQ.z);
    fp2_mul(&R->z,&R->z,&A->Z0);
    fp2_mul(&R->t,&tP.t,&tQ.t);
    fp2_mul(&R->t,&R->t,&A->T0);

    fp2_t t1,t2;
    fp2_mul(&t1,&PQ->x,&PQ->y);
    fp2_mul(&t2,&PQ->z,&PQ->t);

    hadamard(R,R);
    fp2_mul(&R->x,&R->x,&t2);
    fp2_mul(&R->x,&R->x,&PQ->y);
    fp2_mul(&R->y,&R->y,&t2);
    fp2_mul(&R->y,&R->y,&PQ->x);
    fp2_mul(&R->z,&R->z,&t1);
    fp2_mul(&R->z,&R->z,&PQ->t);
    fp2_mul(&R->t,&R->t,&t1);
    fp2_mul(&R->t,&R->t,&PQ->z);

}




/// TODO : see if this needed, and if yes, finish to implement correctly
// /**
//  * @brief Compute the scalar multiplication of a point in the theta struc A
//  *
//  * @param R Output: the theta_point 
//  * @param A a theta structure
//  * @param P a theta point in the theta structure A
//  * @param k a scalar 
//  * R = [k] P
//  * /!\ assumes that no coordinates is zero and that the precomputation of A has been done
//  *  
//    */
// void mul(theta_point_t* R, theta_structure_t const* A, digit_t const* k, theta_point_t const* P)
// {
//     theta_point_t R0, R1, R2;
//     digit_t mask;
//     unsigned int bit = 0, prevbit = 0, swap;

//     // R0 <- P, R1 <- P, R2 <- [2]P
//     fp2_copy(&R1.x, &P->x);
//     fp2_copy(&R1.z, &P->z);
//     fp2_copy(&R1.y, &P->y);
//     fp2_copy(&R1.t, &P->t);
//     fp2_copy(&R0.x, &P->x);
//     fp2_copy(&R0.z, &P->z);
//     fp2_copy(&R0.y, &P->y);
//     fp2_copy(&R0.t, &P->t);

//     double_point(&R2,&A,P);

//     // Main loop
//     for (int i = BITS-1; i >= 0; i--) {                                          
//         bit = (k[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;                         
//         swap = bit ^ prevbit;
//         prevbit = bit;
//         mask = 0 - (digit_t)swap;

//         swap_points(&R0, &R1, mask);
//         xDBLADD(&R0, &R1, &R0, &R1, P, &A24);
//     }
//     swap = 0 ^ prevbit;
//     mask = 0 - (digit_t)swap;
//     swap_points(&R0, &R1, mask);

//     fp2_copy(&Q->x, &R0.x);
//     fp2_copy(&Q->z, &R0.z);
// }