#include "theta_isogenies.h"
#include <assert.h>


void choose_index_theta_point(fp2_t *res,int ind, const theta_point_t *T) {
    int t = ind%4;
    if (t == 0) {
        fp2_copy(res,&T->x);
    }
    else if (t==1) {
        fp2_copy(res,&T->y);
    }
    else if (t==2) {
        fp2_copy(res,&T->z);
    }
    else {
        fp2_copy(res,&T->t);
    }
}

void set_index_theta_point(theta_point_t *res,int ind, const fp2_t *val) {
    int t = ind%4;
    if (t == 0) {
        fp2_copy(&res->x,val);
    }
    else if (t==1) {
        fp2_copy(&res->y,val);
    }
    else if (t==2) {
        fp2_copy(&res->z,val);
    }
    else {
        fp2_copy(&res->t,val);
    }
}

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

// compute the theta_point corresponding to the couple of point T on an elliptic product
void base_change( theta_point_t *out, const theta_gluing_t *phi, const theta_couple_point_t *T) {
    fp2_t a,b,c,d,x1,x2;
    if (fp2_is_zero(&T->P1.z) && fp2_is_zero(&T->P1.x)) {
        fp2_set(&x1,1);       
    }
    else {
        x1 = T->P1.x;
    }
    if (fp2_is_zero(&T->P2.z) && fp2_is_zero(&T->P2.x)) {
        fp2_set(&x2,1);       
    }
    else {
        x2 = T->P2.x;
    }
    // a = P1.x P2.x, b = P1.x P2.z, c =P1.z P2.x, d = P1.z P2.z  
    fp2_mul(&a,&x1,&x2);
    fp2_mul(&b,&x1,&T->P2.z);
    fp2_mul(&c,&x2,&T->P1.z);
    fp2_mul(&d,&T->P1.z,&T->P2.z);

    // Apply the matrix
    fp2_mul(&out->x,&a,&phi->M00);
    fp2_mul(&x1,&b,&phi->M01);
    fp2_add(&out->x,&out->x,&x1);
    fp2_mul(&x1,&c,&phi->M02);
    fp2_add(&out->x,&out->x,&x1);
    fp2_mul(&x1,&d,&phi->M03);
    fp2_add(&out->x,&out->x,&x1);

    fp2_mul(&out->y,&a,&phi->M10);
    fp2_mul(&x1,&b,&phi->M11);
    fp2_add(&out->y,&out->y,&x1);
    fp2_mul(&x1,&c,&phi->M12);
    fp2_add(&out->y,&out->y,&x1);
    fp2_mul(&x1,&d,&phi->M13);
    fp2_add(&out->y,&out->y,&x1);

    fp2_mul(&out->z,&a,&phi->M20);
    fp2_mul(&x1,&b,&phi->M21);
    fp2_add(&out->z,&out->z,&x1);
    fp2_mul(&x1,&c,&phi->M22);
    fp2_add(&out->z,&out->z,&x1);
    fp2_mul(&x1,&d,&phi->M23);
    fp2_add(&out->z,&out->z,&x1);

    fp2_mul(&out->t,&a,&phi->M30);
    fp2_mul(&x1,&b,&phi->M31);
    fp2_add(&out->t,&out->t,&x1);
    fp2_mul(&x1,&c,&phi->M32);
    fp2_add(&out->t,&out->t,&x1);
    fp2_mul(&x1,&d,&phi->M33);
    fp2_add(&out->t,&out->t,&x1);



}

/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing 
 * @param K1_8 a couple point
 * @param E12 an elliptic curve product
 * @param K2_8 a point in E2[8]
 * 
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8) 
 *  
   */
void gluing_comput(theta_gluing_t *out,const theta_couple_curve_t *E12,const theta_couple_point_t *K1_8,const theta_couple_point_t *K2_8) {

    // var init
    fp2_t M1100;fp2_t M1101;fp2_t M1110;fp2_t M1111;
    fp2_t M1200;fp2_t M1201;fp2_t M1210;fp2_t M1211;
    fp2_t M2100;fp2_t M2101;fp2_t M2110;fp2_t M2111;
    fp2_t M2200;fp2_t M2201;fp2_t M2210;fp2_t M2211;
    fp2_t t001,t101,t002,t102,temp;
    
    theta_point_t TT1,TT2;

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

    // apply the base change 
    base_change(&out->T1_8,out,K1_8);
    base_change(&out->T2_8,out,K2_8);

    // computing the codomain 

    // computation of the zero index
    to_squared_theta(&TT1,&out->T1_8);
    to_squared_theta(&TT2,&out->T2_8);

    if (fp2_is_zero(&TT1.x)) { out->zero_idx = 0;}
    else if (fp2_is_zero(&TT1.y)) {out->zero_idx = 1;}
    else if (fp2_is_zero(&TT1.z)) {out->zero_idx = 2;}
    else {out->zero_idx=3;}

    #ifndef NDEBUG
        fp2_t a1,a2;
        choose_index_theta_point(&a1,out->zero_idx,&TT1);
        choose_index_theta_point(&a2,out->zero_idx,&TT2);
        assert(fp2_is_zero(&a1) && fp2_is_zero(&a2));
    #endif 

    choose_index_theta_point(&t001,1+out->zero_idx,&TT2);
    choose_index_theta_point(&t002,2+out->zero_idx,&TT1);
    choose_index_theta_point(&t101,3+out->zero_idx,&TT2);
    choose_index_theta_point(&t102,3+out->zero_idx,&TT1);

    fp2_t t1,t2,t3,t4;
    t1 = t001;t2=t002;t3=t101;t4=t102;

    // TODO batch inversion?
    fp2_inv(&t001);fp2_inv(&t002);fp2_inv(&t102);fp2_inv(&t101);

    // Compute A,B,C,D
    fp2_set(&temp,0);
    set_index_theta_point(&out->codomain,out->zero_idx,&temp);
    fp2_mul(&temp,&t101,&t1);
    set_index_theta_point(&out->codomain,1+out->zero_idx,&temp);
    fp2_mul(&temp,&t2,&t102);
    set_index_theta_point(&out->codomain,2+out->zero_idx,&temp);
    fp2_set(&temp,1);
    set_index_theta_point(&out->codomain,3+out->zero_idx,&temp);

    // compute precomp 
    fp2_set(&temp,0);
    set_index_theta_point(&out->precomputation,out->zero_idx,&temp);
    fp2_mul(&temp,&t001,&t3);
    set_index_theta_point(&out->precomputation,1+out->zero_idx,&temp);
    fp2_mul(&temp,&t4,&t002);
    set_index_theta_point(&out->precomputation,2+out->zero_idx,&temp);
    fp2_set(&temp,1);
    set_index_theta_point(&out->precomputation,3+out->zero_idx,&temp);

    // compute the final codomain 
    hadamard(&out->codomain,&out->codomain);
}

void gluing_eval(theta_point_t *image,const theta_couple_point_t *P,const theta_couple_curve_t *E12, const theta_gluing_t *phi){
    
    theta_couple_point_t Pt;
    theta_point_t T,Tt;

    fp2_t x,y,z,t,temp;

    // First we translate the point 
    add_couple_point(&Pt,E12,P,&phi->K1_4);

    // apply the basis change
    base_change(&T,phi,P);
    base_change(&Tt,phi,&Pt);

    // apply the to_squared_theta transform
    to_squared_theta(&T,&T);
    to_squared_theta(&Tt,&Tt);

    // compute y,z,t
    choose_index_theta_point(&y,1+phi->zero_idx,&T);
    choose_index_theta_point(&temp,1+phi->zero_idx,&phi->precomputation);
    fp2_mul(&y,&y,&temp);
    choose_index_theta_point(&z,2+phi->zero_idx,&T);
    choose_index_theta_point(&temp,2+phi->zero_idx,&phi->precomputation);
    fp2_mul(&z,&z,&temp);
    choose_index_theta_point(&z,3+phi->zero_idx,&T);

    //  normalize
    if (!fp2_is_zero(&z)) {
        fp2_copy(&x,&z);
        choose_index_theta_point(&temp,3+phi->zero_idx,&Tt);
        fp2_inv(&temp);
        fp2_mul(&x,&x,&temp);
    }
    else {
        choose_index_theta_point(&temp,2+phi->zero_idx,&Tt);
        choose_index_theta_point(&x,2+phi->zero_idx,&phi->precomputation);
        fp2_mul(&temp,&temp,&x);
        fp2_inv(&temp);
        fp2_copy(&x,&t);
        fp2_mul(&x,&x,&temp);
    }

    // recover x
    choose_index_theta_point(&temp,1+phi->zero_idx,&Tt);
    fp2_mul(&x,&x,&temp);
    choose_index_theta_point(&temp,1+phi->zero_idx,&phi->precomputation);
    fp2_mul(&x,&x,&temp);

    // fill the image coordinates 
    set_index_theta_point(image,phi->zero_idx,&x);
    set_index_theta_point(image,1+phi->zero_idx,&y);
    set_index_theta_point(image,2+phi->zero_idx,&z);
    set_index_theta_point(image,3+phi->zero_idx,&t);

    // hadamard
    hadamard(image,image);
}


/**
 * @brief Compute  a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_gluing 
 * @param A a theta null point for the domain
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param bool1 a boolean
 * @param boo2 a boolean
 *   
 * out : A -> B of kernel [4](T1_8,T2_8)
 * bool1 controls if the domain is in standard or dual coordinates
 * bool2 controls if the codomain is in standard or dual coordinates 
 *  
   */
void theta_isogeny_comput(theta_isogeny_t *out,const theta_structure_t *A,const theta_point_t *T1_8,const theta_point_t *T2_8,int bool1, int bool2) {
    out->bool1=bool1;
    out->bool2=bool2;
    out->domain=*A;
    out->T1_8=*T1_8;
    out->T2_8=*T2_8;
    out->codomain.precomputation=0;

    theta_point_t TT1,TT2;

    fp2_t xA_inv,zA_inv,tB_inv;

    if (bool1) {
        // TODO we do not need all the values of the results, so we may save some operations? 
        hadamard(&TT1,T1_8);
        to_squared_theta(&TT1,&TT1);
        hadamard(&TT2,T2_8);
        to_squared_theta(&TT2,&TT2);
    } 
    else {
        to_squared_theta(&TT1,T1_8);
        to_squared_theta(&TT2,T2_8);
    }

    if (!bool1 && A->precomputation) {
        xA_inv = TT1.x;
        zA_inv = TT2.x;
        tB_inv = TT2.y;
        //TODO bach_inversion?
        fp2_inv(&xA_inv);
        fp2_inv(&zA_inv);
        fp2_inv(&tB_inv);

        fp2_set(&out->precomputation.x,1);
        fp2_set(&out->codomain.null_point.x,1);

        // computation of B,C,D for the codomain
        fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT1.t);
        fp2_mul(&out->codomain.null_point.t,&out->codomain.null_point.t,&out->codomain.null_point.y);

        // computation B_inv,C_inv,D_inv for the precomputation
        fp2_mul(&out->precomputation.y,&out->codomain.null_point.y,&out->domain.Y0);
        fp2_mul(&out->precomputation.z,&out->codomain.null_point.z,&out->domain.Z0);
        fp2_mul(&out->precomputation.t,&out->codomain.null_point.t,&out->domain.T0);
    }
    else {
        fp2_t xB_inv,zC_inv,tD_inv;
        xA_inv = TT1.x;
        zA_inv = TT2.x;
        tB_inv = TT2.y;
        xB_inv=TT1.y;
        zC_inv=TT2.z;
        tD_inv=TT2.t;
        //TODO bach_inversion?
        fp2_inv(&xA_inv);
        fp2_inv(&zA_inv);
        fp2_inv(&tB_inv);
        fp2_inv(&xB_inv);
        fp2_inv(&zC_inv);
        fp2_inv(&tD_inv);


        fp2_set(&out->precomputation.x,1);
        fp2_set(&out->codomain.null_point.x,1);

        // computation of B,C,D for the codomain
        fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT1.t);
        fp2_mul(&out->codomain.null_point.t,&out->codomain.null_point.t,&out->codomain.null_point.y);

        // computation of B_inv,C_inv,D_inv for the precomputation
        fp2_mul(&out->precomputation.y,&xB_inv,&TT1.x);
        fp2_mul(&out->precomputation.z,&zC_inv,&TT2.x);
        fp2_mul(&out->precomputation.t,&tD_inv,&TT2.y);
        fp2_mul(&out->precomputation.t,&out->precomputation.t,&out->precomputation.y);

        if (bool2) {
            hadamard(&out->codomain.null_point,&out->codomain.null_point);
        }

    }

}

/**
 * @brief Evaluate a theta isogeny
 *
 * @param out Output: the evaluating point 
 * @param phi a theta isogeny
 * @param P a point in the domain of phi
 *   
 * out = phi(P) 
 *  
   */
void theta_isogeny_eval(theta_point_t *out,const theta_isogeny_t *phi,const theta_point_t *P) {

    if (phi->bool1) {
        hadamard(out,P);
        to_squared_theta(out,out);
    }
    else {
        to_squared_theta(out,P);
    }
    fp2_mul(&out->y,&out->y,&phi->precomputation.y);
    fp2_mul(&out->z,&out->z,&phi->precomputation.z);
    fp2_mul(&out->t,&out->t,&phi->precomputation.t);

    if (phi->bool2) {
        hadamard(out,out);
    }
    
}