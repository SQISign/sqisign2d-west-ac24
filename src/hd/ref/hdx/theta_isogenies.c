#include "theta_isogenies.h"
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>


//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.im[i]);
    printf("\n");
}

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
    fp2_t det,u,invz,temp;
    
    // temp = [2]P
    ec_dbl(&P2,E,P);

    // invz = 1/P.z  
    invz = P->z;
    fp2_inv(&invz);

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
    fp2_add(&t001,&t001,&temp);  

    fp2_mul(&t101,&M1110,&M2100);
    fp2_mul(&temp,&M1111,&M2110); 
    fp2_add(&t101,&t101,&temp);

    fp2_mul(&t002,&M1200,&M2200);
    fp2_mul(&temp,&M1201,&M2210); 
    fp2_add(&t002,&t002,&temp);  

    fp2_mul(&t102,&M1210,&M2200);
    fp2_mul(&temp,&M1211,&M2210); 
    fp2_add(&t102,&t102,&temp); 

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
        char s1[2000],s2[2000];
        fp2_print(s1,a1);
        fp2_print(s2,a2);
        printf("%s %s %d \n",s1,s2,out->zero_idx);
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

void apply_isomorphism(theta_point_t *res,const theta_splitting_t *out,const theta_point_t *P) {
        fp2_t x1;
        fp2_mul(&res->x,&P->x,&out->M00);
        fp2_mul(&x1,&P->y,&out->M01);
        fp2_add(&res->x,&res->x,&x1);
        fp2_mul(&x1,&P->z,&out->M02);
        fp2_add(&res->x,&res->x,&x1);
        fp2_mul(&x1,&P->t,&out->M03);
        fp2_add(&res->x,&res->x,&x1);

        fp2_mul(&res->y,&P->x,&out->M10);
        fp2_mul(&x1,&P->y,&out->M11);
        fp2_add(&res->y,&res->y,&x1);
        fp2_mul(&x1,&P->z,&out->M12);
        fp2_add(&res->y,&res->y,&x1);
        fp2_mul(&x1,&P->t,&out->M13);
        fp2_add(&res->y,&res->y,&x1);

        fp2_mul(&res->z,&P->x,&out->M20);
        fp2_mul(&x1,&P->y,&out->M21);
        fp2_add(&res->z,&res->z,&x1);
        fp2_mul(&x1,&P->z,&out->M22);
        fp2_add(&res->z,&res->z,&x1);
        fp2_mul(&x1,&P->t,&out->M23);
        fp2_add(&res->z,&res->z,&x1);

        fp2_mul(&res->t,&P->x,&out->M30);
        fp2_mul(&x1,&P->y,&out->M31);
        fp2_add(&res->t,&res->t,&x1);
        fp2_mul(&x1,&P->z,&out->M32);
        fp2_add(&res->t,&res->t,&x1);
        fp2_mul(&x1,&P->t,&out->M33);
        fp2_add(&res->t,&res->t,&x1);
}


/** 
 * @brief Compute the splitting isomorphism from a theta structure to the product theta structure, returns false if the given theta structure is not isomorphic to an elliptic product
 * 
 * @return a boolean indicating if A is isomorphic to an elliptic product
 * @param out: the splitting isomorphism
 * @param A : the theta_structure in consideration
 *
 * out : A -> B where B is theta product associated to ExF an elliptic product 
*/
int splitting_comput(theta_splitting_t *out, const theta_structure_t *A) {

    // init 
    int good=0;
    int even_index[10][2] = { {0,0},{0,1},{0,2},{0,3},{1,0},{1,2},{2,0},{2,1},{3,0},{3,3}};
    int chi_eval[4][4] = {{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}};

    fp2_t U_cst,temp,temp2;

    // enumerate through all indices
    for (int i=0;i<10;i++) {
       fp2_set(&U_cst,0);
       for (int t=0;t<4;t++) {
        choose_index_theta_point(&temp2,t,&A->null_point);
        choose_index_theta_point(&temp,t+even_index[i][1],&A->null_point);
        fp2_mul(&temp,&temp,&temp2);
        fp2_set(&temp2,chi_eval[even_index[i][0]][t]);
        fp2_mul(&temp,&temp,&temp2);
        fp2_add(&U_cst,&U_cst,&temp);    
       }
       if (fp2_is_zero(&U_cst)) {
        good =1+i;
        break;
       }
    }

    // temp = sqrt{-1}
    // TODO precompute this?
    fp2_set(&temp,-1);
    fp2_sqrt(&temp);

    // compute the matrix
    if (good==1) {
        // (0, 0): [1, i, 1, i,
        //          1, -i, -1, i, 
        //          1, i, -1, -i, 
        //          -1, i, -1, i],
        fp2_set(&out->M00,1);
        fp2_copy(&out->M01,&temp);
        fp2_set(&out->M02,1);
        fp2_copy(&out->M03,&temp);
        fp2_set(&out->M10,1);
        fp2_neg(&out->M11,&temp);
        fp2_set(&out->M12,-1);
        fp2_copy(&out->M13,&temp);
        fp2_set(&out->M20,1);
        fp2_copy(&out->M21,&temp);
        fp2_neg(&out->M23,&temp);
        fp2_set(&out->M22,-1);
        fp2_set(&out->M30,-1);
        fp2_copy(&out->M31,&temp);
        fp2_set(&out->M32,-1);
        fp2_copy(&out->M33,&temp);
        
    }
    else if (good==2) {
        // (0, 1): [1, 0, 0, 0,
        //          0, 0, 0, 1, 
        //          0, 0, 1, 0, 
        //          0, -1, 0, 0],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,0);
        fp2_set(&out->M02,0);
        fp2_set(&out->M03,0);
        fp2_set(&out->M10,0);
        fp2_set(&out->M11,0);
        fp2_set(&out->M12,0);
        fp2_set(&out->M13,1);
        fp2_set(&out->M20,0);
        fp2_set(&out->M21,0);
        fp2_set(&out->M22,1);
        fp2_set(&out->M23,0);
        fp2_set(&out->M30,0);
        fp2_set(&out->M31,-1);
        fp2_set(&out->M32,0);
        fp2_set(&out->M33,0);
    }
    else if (good==3) {
        // (0, 2): [1, 0, 0, 0, 
        //          0, 1, 0, 0,
        //          0, 0, 0, 1,
            //      0, 0, -1, 0],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,0);
        fp2_set(&out->M02,0);
        fp2_set(&out->M03,0);
        fp2_set(&out->M10,0);
        fp2_set(&out->M11,1);
        fp2_set(&out->M12,0);
        fp2_set(&out->M13,0);
        fp2_set(&out->M20,0);
        fp2_set(&out->M21,0);
        fp2_set(&out->M22,0);
        fp2_set(&out->M23,1);
        fp2_set(&out->M30,0);
        fp2_set(&out->M31,0);
        fp2_set(&out->M32,-1);
        fp2_set(&out->M33,0);
    }
    else if (good==4) {
        // (0, 3): [1, 0, 0, 0,
        //          0, 1, 0, 0, 
        //          0, 0, 1, 0, 
        //          0, 0, 0, -1],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,0);
        fp2_set(&out->M02,0);
        fp2_set(&out->M03,0);
        fp2_set(&out->M10,0);
        fp2_set(&out->M11,1);
        fp2_set(&out->M12,0);
        fp2_set(&out->M13,0);
        fp2_set(&out->M20,0);
        fp2_set(&out->M21,0);
        fp2_set(&out->M22,1);
        fp2_set(&out->M23,0);
        fp2_set(&out->M30,0);
        fp2_set(&out->M31,0);
        fp2_set(&out->M32,0);
        fp2_set(&out->M33,-1);

    }
    else if (good==7) {
        // (2, 0): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, -1, -1, 1, 
        //          -1, -1, 1, 1],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,1);
        fp2_set(&out->M02,1);
        fp2_set(&out->M03,1);
        fp2_set(&out->M10,1);
        fp2_set(&out->M11,-1);
        fp2_set(&out->M12,1);
        fp2_set(&out->M13,-1);
        fp2_set(&out->M20,1);
        fp2_set(&out->M21,-1);
        fp2_set(&out->M22,-1);
        fp2_set(&out->M23,1);
        fp2_set(&out->M30,-1);
        fp2_set(&out->M31,-1);
        fp2_set(&out->M32,1);
        fp2_set(&out->M33,1);
    }
    else if (good==8) {
        //(2, 1): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, -1, -1, 1, 
        //          1, 1, -1, -1],

        fp2_set(&out->M00,1);
        fp2_set(&out->M01,1);
        fp2_set(&out->M02,1);
        fp2_set(&out->M03,1);
        fp2_set(&out->M10,1);
        fp2_set(&out->M11,-1);
        fp2_set(&out->M12,1);
        fp2_set(&out->M13,-1);
        fp2_set(&out->M20,1);
        fp2_set(&out->M21,-1);
        fp2_set(&out->M22,-1);
        fp2_set(&out->M23,1);
        fp2_set(&out->M30,1);
        fp2_set(&out->M31,1);
        fp2_set(&out->M32,-1);
        fp2_set(&out->M33,-1);
    }
    else if (good==5){
        // (1, 0): [1, 1, 1, 1, 
        //          1, -1, -1, 1, 
        //          1, 1, -1, -1, 
        //          -1, 1, -1, 1],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,1);
        fp2_set(&out->M02,1);
        fp2_set(&out->M03,1);
        fp2_set(&out->M10,1);
        fp2_set(&out->M11,-1);
        fp2_set(&out->M12,-1);
        fp2_set(&out->M13,1);
        fp2_set(&out->M20,1);
        fp2_set(&out->M21,1);
        fp2_set(&out->M22,-1);
        fp2_set(&out->M23,-1);
        fp2_set(&out->M30,-1);
        fp2_set(&out->M31,1);
        fp2_set(&out->M32,-1);
        fp2_set(&out->M33,1);

    }
    else if (good==6) {
        //(1, 2): [1, 0, 0, 0, 
        //          0, 1, 0, 0, 
        //          0, 0, 0, 1, 
        //          0, 0, 1, 0],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,0);
        fp2_set(&out->M02,0);
        fp2_set(&out->M03,0);
        fp2_set(&out->M10,0);
        fp2_set(&out->M11,1);
        fp2_set(&out->M12,0);
        fp2_set(&out->M13,0);
        fp2_set(&out->M20,0);
        fp2_set(&out->M21,0);
        fp2_set(&out->M22,0);
        fp2_set(&out->M23,1);
        fp2_set(&out->M30,0);
        fp2_set(&out->M31,0);
        fp2_set(&out->M32,1);
        fp2_set(&out->M33,0);

    }
    else if (good==9) {
        // (3, 0): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, 1, -1, -1, 
        //          -1, 1, 1, -1],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,1);
        fp2_set(&out->M02,1);
        fp2_set(&out->M03,1);
        fp2_set(&out->M10,1);
        fp2_set(&out->M11,-1);
        fp2_set(&out->M12,1);
        fp2_set(&out->M13,-1);
        fp2_set(&out->M20,1);
        fp2_set(&out->M21,1);
        fp2_set(&out->M22,-1);
        fp2_set(&out->M23,-1);
        fp2_set(&out->M30,-1);
        fp2_set(&out->M31,1);
        fp2_set(&out->M32,1);
        fp2_set(&out->M33,-1);
    }
    else if (good==10) {
        // (3, 3): [1, 0, 0, 0, 
        //          0, 1, 0, 0, 
        //          0, 0, 1, 0,
        //          0, 0, 0, 1],
        fp2_set(&out->M00,1);
        fp2_set(&out->M01,0);
        fp2_set(&out->M02,0);
        fp2_set(&out->M03,0);
        fp2_set(&out->M10,0);
        fp2_set(&out->M11,1);
        fp2_set(&out->M12,0);
        fp2_set(&out->M13,0);
        fp2_set(&out->M20,0);
        fp2_set(&out->M21,0);
        fp2_set(&out->M22,1);
        fp2_set(&out->M23,0);
        fp2_set(&out->M30,0);
        fp2_set(&out->M31,0);
        fp2_set(&out->M32,0);
        fp2_set(&out->M33,1);

    }

    // now we apply the isomorphism if it was computed
    if (good) {
        apply_isomorphism(&out->B.null_point,out,&A->null_point);
    }
    return good;
}


void theta_product_structure_to_elliptic_product(theta_couple_curve_t *E12, theta_structure_t *A) {
    fp2_t xx,yy,temp1,temp2;

    // xx = x², yy = y² 
    fp2_sqr(&xx,&A->null_point.x);
    fp2_sqr(&yy,&A->null_point.y);

    // A1 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1,&xx,&yy);
    fp2_sub(&temp2,&xx,&yy);
    fp2_mul(&E12->E1.A,&temp1,&temp2);
    fp2_sqr(&temp1,&temp1);
    fp2_sqr(&temp2,&temp2);
    fp2_add(&temp1,&temp1,&temp2);
    fp2_inv(&E12->E1.A);
    fp2_mul(&E12->E1.A,&E12->E1.A,&temp1);
    fp2_neg(&E12->E1.A,&E12->E1.A);

    // same with t,z
    fp2_sqr(&xx,&A->null_point.z);
    fp2_sqr(&yy,&A->null_point.t);

    // A2 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1,&xx,&yy);
    fp2_sub(&temp2,&xx,&yy);
    fp2_mul(&E12->E2.A,&temp1,&temp2);
    fp2_sqr(&temp1,&temp1);
    fp2_sqr(&temp2,&temp2);
    fp2_add(&temp1,&temp1,&temp2);
    fp2_inv(&E12->E2.A);
    fp2_mul(&E12->E2.A,&E12->E2.A,&temp1);
    fp2_neg(&E12->E2.A,&E12->E2.A);
}

/**
 * @brief Compute  a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the theta_chain
 * @param n : the length of the isogeny chain
 * @param E12 an elliptic curve product 
 * @param T1 a couple point on E12[2^(n+2)]
 * @param T2 a couple point on E12[2^(n+2)]
 *   
 * out : E1xE2 -> E3xE4 of kernel [4](T1,T2) 
 *  
   */
void theta_chain_comput(theta_chain_t *out,int n,const theta_couple_curve_t *E12,const theta_couple_point_t *T1,const theta_couple_point_t *T2) {

    theta_couple_point_t P1,P2;
    theta_point_t Q1,Q2,R1,R2;
    theta_isogeny_t steps[n-1];
    theta_structure_t codomain;

    // TODO use a better strategy

    // init of the isogeny chain
    out->domain=*E12;
    out->length=n;
    out->T1=*T1;
    out->T2=*T2;
    out->steps=malloc((n-1)*sizeof(theta_isogeny_t));

    // First, we compute the first step  
    // multiply by 2^n-1
    double_couple_point_iter(&P1,n-1,E12,T1);
    double_couple_point_iter(&P2,n-1,E12,T2);

    #ifndef NDEBUG
        // checking that the points have order 8 
        ec_point_t test1,test2;
        test1=P1.P1;
        test2=P1.P2;
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test2,&E12->E2,&test2);
        ec_dbl(&test2,&E12->E2,&test2);
        assert(!fp2_is_zero(&test1.z));
        assert(!fp2_is_zero(&test2.z));
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test2,&E12->E2,&test2);
        assert(fp2_is_zero(&test1.z));
        assert(fp2_is_zero(&test2.z));
        test1=P2.P1;
        test2=P2.P2;
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test2,&E12->E2,&test2);
        ec_dbl(&test2,&E12->E2,&test2);
        assert(!fp2_is_zero(&test1.z));
        assert(!fp2_is_zero(&test2.z));
        ec_dbl(&test1,&E12->E1,&test1);
        ec_dbl(&test2,&E12->E2,&test2);
        assert(fp2_is_zero(&test1.z));
        assert(fp2_is_zero(&test2.z));
    #endif
    

    // compute the gluing isogeny 
    gluing_comput(&out->first_step,E12,&P1,&P2);

    // push the kernel through the gluing isogeny
    gluing_eval(&Q1,T1,E12,&out->first_step);
    gluing_eval(&Q2,T2,E12,&out->first_step);

    // set-up the theta_structure for the first codomain 
    codomain.null_point=out->first_step.codomain;
    codomain.precomputation=0;
    theta_precomputation(&codomain);

    for (int i=0;i<n-1;i++) {

        // computing the kernel of the next step
        double_iter(&R1,&codomain,&Q1,n-i-1);
        double_iter(&R2,&codomain,&Q2,n-i-1);
    
        // computing the next step
        if (i==n-3) {
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,0,0);
        }
        else if (i==n-2) {
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,1,0);
        }
        else {
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,0,1);
        }

        // updating the codomain
        codomain=steps[i].codomain;

        // pushing the kernel
        if (i< n-2) {
            theta_isogeny_eval(&Q1,&steps[i],&Q1);
            theta_isogeny_eval(&Q2,&steps[i],&Q2);
        }
    
    }

    // copying the steps
    out->steps=steps;

    //final splitting step
    int is_split = splitting_comput(&out->last_step,&steps[n-2].codomain);
    assert(is_split);

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain,&out->last_step.B);

}


void theta_point_to_montgomery_point(theta_couple_point_t *P12, const theta_point_t *P, const theta_structure_t *A){

    fp2_t temp;

    // P1.X = A.null_point.x * P.x + A.null_point.y * P.y
    // P1.Z = A.null_point.x * P.x - A.null_point.y * P.y
    fp2_mul(&P12->P1.x,&A->null_point.x,&P->x);
    fp2_mul(&temp,&A->null_point.y,&P->y);
    fp2_sub(&P12->P1.z,&P12->P1.x,&temp);
    fp2_add(&P12->P1.x,&P12->P1.x,&temp);

    // P2.X = A.null_point.z * P.z + A.null_point.t * P.t
    // P2.Z = A.null_point.z * P.z - A.null_point.t * P.t
    fp2_mul(&P12->P2.x,&A->null_point.z,&P->z);
    fp2_mul(&temp,&A->null_point.t,&P->t);
    fp2_sub(&P12->P2.z,&P12->P2.x,&temp);
    fp2_add(&P12->P2.x,&P12->P2.x,&temp);
    

}

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple point on E12, 
 *   
 * phi : E1xE2 -> E3xE4 of kernel 
 * P12 in E1xE2
 * out = phi(P12) in E3xE4 
 *  
   */
void theta_chain_eval(theta_couple_point_t *out,theta_chain_t *phi,theta_couple_point_t *P12) {
    
    theta_point_t temp;

    // first, we apply the gluing 
    gluing_eval(&temp,P12,&phi->domain,&phi->first_step);

    // then, we apply the successive isogenies
    for (int i=0;i<phi->length-1;i++) {
        theta_isogeny_eval(&temp,&phi->steps[i],&temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp,&phi->last_step,&temp);

    // finaly the send the result to the elliptic product
    theta_point_to_montgomery_point(out,&temp,&phi->last_step.B);

}