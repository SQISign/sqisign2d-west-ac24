#include "theta_isogenies.h"
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>


//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
void fp2_print(char *name, fp2_t const a){
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

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}

static void theta_print(char *name, theta_point_t P) {
    fp2_t a;
    assert(!fp2_is_zero(&P.x));
    fp2_copy(&a,&P.x);
    fp2_inv(&a);
    fp2_mul(&a,&a,&P.t);
    fp2_print(name,a);
}


// set the finite field element in fp2 to 1
void fp2_setone(fp2_t *a) {
    fp_set(a->im,0);
    fp_mont_setone(a->re);
}

void fp2_sett(fp2_t *a,int t) {
    if (t==1) {
        fp2_setone(a);
    }
    else if (t==0) {
        fp2_set(a,0);
    }
    else if (t==-1) {
        fp2_setone(a);
        fp2_neg(a,a);
    }
    else {
        fp2_set(a,t);
        fp_tomont(a->re,a->re);
        fp_tomont(a->im,a->im);
    }
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
    else if (t==3) {
        fp2_copy(res,&T->t);
    }
    else {
        assert(0);
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
    else if (t==3) {
        fp2_copy(&res->t,val);
    }
    else {
        assert(0);
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
        fp2_sett(&x1,1);       
    }
    else {
        x1 = T->P1.x;
    }
    if (fp2_is_zero(&T->P2.z) && fp2_is_zero(&T->P2.x)) {
        fp2_sett(&x2,1);       
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
    fp2_sett(&out->M00,1);
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
    fp2_add(&out->M21,&out->M21,&temp);

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

    theta_print("T1_8",out->T1_8);
    theta_print("T2_8",out->T2_8);
    printf("\n");
    // computing the codomain 

    // computation of the zero index
    to_squared_theta(&TT1,&out->T1_8);
    to_squared_theta(&TT2,&out->T2_8);

    theta_print("HS(T1)",TT1);
    theta_print("HS(T2)",TT2);
    printf("\n");

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

    choose_index_theta_point(&t001,1^out->zero_idx,&TT2);
    choose_index_theta_point(&t002,2^out->zero_idx,&TT1);
    choose_index_theta_point(&t101,3^out->zero_idx,&TT2);
    choose_index_theta_point(&t102,3^out->zero_idx,&TT1);

    fp2_t t1,t2,t3,t4;
    // t1 = t001;t2=t002;t3=t101;t4=t102;
    fp2_copy(&t1,&t001);
    fp2_copy(&t2,&t002);
    fp2_copy(&t3,&t101);
    fp2_copy(&t4,&t102);

    // TODO batch inversion?
    fp2_inv(&t001);fp2_inv(&t002);fp2_inv(&t102);fp2_inv(&t101);

    // Compute A,B,C,D
    fp2_sett(&temp,0);
    set_index_theta_point(&out->codomain,0^out->zero_idx,&temp);
    // 
    fp2_mul(&temp,&t101,&t1);
    set_index_theta_point(&out->codomain,1^out->zero_idx,&temp);
    fp2_mul(&temp,&t2,&t102);
    // fp2_mul(&temp,&t4,&t002);
    set_index_theta_point(&out->codomain,2^out->zero_idx,&temp);
    fp2_sett(&temp,1);
    set_index_theta_point(&out->codomain,3^out->zero_idx,&temp);
    
    choose_index_theta_point(&temp,2^out->zero_idx,&out->codomain);
    fp2_print("ABCD",temp);
    // compute precomp 
    fp2_sett(&temp,0);
    set_index_theta_point(&out->precomputation,out->zero_idx,&temp);
    fp2_mul(&temp,&t001,&t3);
    set_index_theta_point(&out->precomputation,1^out->zero_idx,&temp);
    // fp2_mul(&temp,&t2,&t102);
    fp2_mul(&temp,&t4,&t002);
    set_index_theta_point(&out->precomputation,2^out->zero_idx,&temp);
    fp2_sett(&temp,1);
    set_index_theta_point(&out->precomputation,3^out->zero_idx,&temp);
    choose_index_theta_point(&temp,2^out->zero_idx,&out->precomputation);
    fp2_print("precomp",temp);
    printf("\n");

    // compute the final codomain 
    hadamard(&out->codomain,&out->codomain);
}


// sub routine of the gluing eval
void gluing_eval_point(theta_point_t *image1, const theta_couple_point_t *P, const theta_couple_point_t *Pt,const theta_gluing_t *phi) {

    theta_point_t T,Tt;
    fp2_t x,y,z,t,temp;

    // apply the basis change
    base_change(&T,phi,P);
    base_change(&Tt,phi,Pt);

    // apply the to_squared_theta transform
    to_squared_theta(&T,&T);
    to_squared_theta(&Tt,&Tt);

    // compute y,z,t
    choose_index_theta_point(&y,1^phi->zero_idx,&T);
    choose_index_theta_point(&temp,1^phi->zero_idx,&phi->precomputation);
    fp2_mul(&y,&y,&temp);
    choose_index_theta_point(&z,2^phi->zero_idx,&T);
    choose_index_theta_point(&temp,2^phi->zero_idx,&phi->precomputation);
    fp2_mul(&z,&z,&temp);
    choose_index_theta_point(&t,3^phi->zero_idx,&T);

    //  normalize
    if (!fp2_is_zero(&z)) {
        fp2_copy(&x,&z);
        choose_index_theta_point(&temp,3^phi->zero_idx,&Tt);
        fp2_inv(&temp);
        fp2_mul(&x,&x,&temp);
    }
    else {
        choose_index_theta_point(&temp,2^phi->zero_idx,&Tt);
        choose_index_theta_point(&x,2^phi->zero_idx,&phi->precomputation);
        fp2_mul(&temp,&temp,&x);
        fp2_inv(&temp);
        fp2_copy(&x,&t);
        fp2_mul(&x,&x,&temp);
    }

    // recover x
    choose_index_theta_point(&temp,1^phi->zero_idx,&Tt);
    fp2_mul(&x,&x,&temp);
    choose_index_theta_point(&temp,1^phi->zero_idx,&phi->precomputation);
    fp2_mul(&x,&x,&temp);

    // fill the image coordinates 
    set_index_theta_point(image1,0^phi->zero_idx,&x);
    set_index_theta_point(image1,1^phi->zero_idx,&y);
    set_index_theta_point(image1,2^phi->zero_idx,&z);
    set_index_theta_point(image1,3^phi->zero_idx,&t);

    // hadamard
    hadamard(image1,image1);
}


/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a basis
 *
 * @param image1 Output: the theta_point of the image of the first couple of points
 * @param image2 Output : the theta point of the image of the second couple of points
 * @param P : a couple point in E12
 * @param Q : a couple point in E12
 * @param PmQ : a couple point in E12 corresponding to P-Q
 * @param a : ibz
 * @param b : ibz 
 * @param E12 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A 
 * 
 * The integers a,b are such that the phi.K1_4 = a P + b Q
 * out : phi( P ), Phi (Q)  
 *  
   */
void gluing_eval_basis(theta_point_t *image1, theta_point_t *image2, const theta_couple_point_t *P,const theta_couple_point_t *Q, const theta_couple_point_t *PmQ,const ibz_t *a, const ibz_t *b,const theta_couple_curve_t *E12, const theta_gluing_t *phi) {
    
    theta_couple_point_t Pt;
    
    ec_basis_t bas1,bas2;
    
    digit_t scalars[2][NWORDS_ORDER] = {0};    
    ibz_t ibz_temp;
    ibz_init(&ibz_temp);


    // set-up the basis
    bas1.P = P->P1;
    bas1.Q = Q->P1;
    bas1.PmQ = PmQ->P1;
    bas2.P = P->P2;
    bas2.Q = Q->P2;
    bas2.PmQ = PmQ->P2;

    // scalars = a+1,b
    ibz_add(&ibz_temp,a,&ibz_const_one);
    ibz_to_digit_array(scalars[0],&ibz_temp);
    ibz_to_digit_array(scalars[1],b);

    // First we translate the point by phi->K1_4
    // add_couple_point(&Pt,E12,P,&phi->K1_4);
    ec_biscalar_mul(&Pt.P1,&E12->E1,scalars[0],scalars[1],&bas1);
    ec_biscalar_mul(&Pt.P2,&E12->E2,scalars[0],scalars[1],&bas2);

    // then we evaluate the gluing 
    gluing_eval_point(image1,P,&Pt,phi);    

    // we do the same on the second point 

    //scalars = a,b+1
    ibz_add(&ibz_temp,b,&ibz_const_one);
    ibz_to_digit_array(scalars[0],a);
    ibz_to_digit_array(scalars[1],&ibz_temp);

    // First we translate the point by phi->K1_4
    // add_couple_point(&Pt,E12,P,&phi->K1_4);
    ec_biscalar_mul(&Pt.P1,&E12->E1,scalars[0],scalars[1],&bas1);
    ec_biscalar_mul(&Pt.P2,&E12->E2,scalars[0],scalars[1],&bas2);

    // then we evaluate the gluing 
    gluing_eval_point(image2,Q,&Pt,phi); 

    ibz_finalize(&ibz_temp);
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

        fp2_sett(&out->precomputation.x,1);
        fp2_sett(&out->codomain.null_point.x,1);

        // computation of B,C,D for the codomain
        fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT2.t);
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


        fp2_sett(&out->precomputation.x,1);
        fp2_sett(&out->codomain.null_point.x,1);

        // computation of B,C,D for the codomain
        fp2_mul(&out->codomain.null_point.y,&xA_inv,&TT1.y);
        fp2_mul(&out->codomain.null_point.z,&zA_inv,&TT2.z);
        fp2_mul(&out->codomain.null_point.t,&tB_inv,&TT2.t);
        fp2_mul(&out->codomain.null_point.t,&out->codomain.null_point.t,&out->codomain.null_point.y);

        // computation of B_inv,C_inv,D_inv for the precomputation
        fp2_mul(&out->precomputation.y,&xB_inv,&TT1.x);
        fp2_mul(&out->precomputation.z,&zC_inv,&TT2.x);
        fp2_mul(&out->precomputation.t,&tD_inv,&TT2.y);
        fp2_mul(&out->precomputation.t,&out->precomputation.t,&out->precomputation.y);

        
    }
    if (bool2) {
        hadamard(&out->codomain.null_point,&out->codomain.null_point);
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

        theta_point_t temp;

        fp2_t x1;
        fp2_mul(&temp.x,&P->x,&out->M00);
        fp2_mul(&x1,&P->y,&out->M01);
        fp2_add(&temp.x,&temp.x,&x1);
        fp2_mul(&x1,&P->z,&out->M02);
        fp2_add(&temp.x,&temp.x,&x1);
        fp2_mul(&x1,&P->t,&out->M03);
        fp2_add(&temp.x,&temp.x,&x1);

        fp2_mul(&temp.y,&P->x,&out->M10);
        fp2_mul(&x1,&P->y,&out->M11);
        fp2_add(&temp.y,&temp.y,&x1);
        fp2_mul(&x1,&P->z,&out->M12);
        fp2_add(&temp.y,&temp.y,&x1);
        fp2_mul(&x1,&P->t,&out->M13);
        fp2_add(&temp.y,&temp.y,&x1);

        fp2_mul(&temp.z,&P->x,&out->M20);
        fp2_mul(&x1,&P->y,&out->M21);
        fp2_add(&temp.z,&temp.z,&x1);
        fp2_mul(&x1,&P->z,&out->M22);
        fp2_add(&temp.z,&temp.z,&x1);
        fp2_mul(&x1,&P->t,&out->M23);
        fp2_add(&temp.z,&temp.z,&x1);

        fp2_mul(&temp.t,&P->x,&out->M30);
        fp2_mul(&x1,&P->y,&out->M31);
        fp2_add(&temp.t,&temp.t,&x1);
        fp2_mul(&x1,&P->z,&out->M32);
        fp2_add(&temp.t,&temp.t,&x1);
        fp2_mul(&x1,&P->t,&out->M33);
        fp2_add(&temp.t,&temp.t,&x1);

        fp2_copy(&res->x,&temp.x);
        fp2_copy(&res->y,&temp.y);
        fp2_copy(&res->z,&temp.z);
        fp2_copy(&res->t,&temp.t);
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
       fp2_sett(&U_cst,0);
       for (int t=0;t<4;t++) {
        choose_index_theta_point(&temp2,t,&A->null_point);
        choose_index_theta_point(&temp,t^even_index[i][1],&A->null_point);
        fp2_mul(&temp,&temp,&temp2);
        fp2_sett(&temp2,chi_eval[even_index[i][0]][t]);
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
    fp2_sett(&temp,-1);
    fp2_sqrt(&temp);

    // compute the matrix
    if (good==1) {
        // (0, 0): [1, i, 1, i,
        //          1, -i, -1, i, 
        //          1, i, -1, -i, 
        //          -1, i, -1, i],
        fp2_sett(&out->M00,1);
        fp2_copy(&out->M01,&temp);
        fp2_sett(&out->M02,1);
        fp2_copy(&out->M03,&temp);
        fp2_sett(&out->M10,1);
        fp2_neg(&out->M11,&temp);
        fp2_sett(&out->M12,-1);
        fp2_copy(&out->M13,&temp);
        fp2_sett(&out->M20,1);
        fp2_copy(&out->M21,&temp);
        fp2_neg(&out->M23,&temp);
        fp2_sett(&out->M22,-1);
        fp2_sett(&out->M30,-1);
        fp2_copy(&out->M31,&temp);
        fp2_sett(&out->M32,-1);
        fp2_copy(&out->M33,&temp);
        
    }
    else if (good==2) {
        // (0, 1): [1, 0, 0, 0,
        //          0, 0, 0, 1, 
        //          0, 0, 1, 0, 
        //          0, -1, 0, 0],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,0);
        fp2_sett(&out->M02,0);
        fp2_sett(&out->M03,0);
        fp2_sett(&out->M10,0);
        fp2_sett(&out->M11,0);
        fp2_sett(&out->M12,0);
        fp2_sett(&out->M13,1);
        fp2_sett(&out->M20,0);
        fp2_sett(&out->M21,0);
        fp2_sett(&out->M22,1);
        fp2_sett(&out->M23,0);
        fp2_sett(&out->M30,0);
        fp2_sett(&out->M31,-1);
        fp2_sett(&out->M32,0);
        fp2_sett(&out->M33,0);
    }
    else if (good==3) {
        // (0, 2): [1, 0, 0, 0, 
        //          0, 1, 0, 0,
        //          0, 0, 0, 1,
            //      0, 0, -1, 0],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,0);
        fp2_sett(&out->M02,0);
        fp2_sett(&out->M03,0);
        fp2_sett(&out->M10,0);
        fp2_sett(&out->M11,1);
        fp2_sett(&out->M12,0);
        fp2_sett(&out->M13,0);
        fp2_sett(&out->M20,0);
        fp2_sett(&out->M21,0);
        fp2_sett(&out->M22,0);
        fp2_sett(&out->M23,1);
        fp2_sett(&out->M30,0);
        fp2_sett(&out->M31,0);
        fp2_sett(&out->M32,-1);
        fp2_sett(&out->M33,0);
    }
    else if (good==4) {
        // (0, 3): [1, 0, 0, 0,
        //          0, 1, 0, 0, 
        //          0, 0, 1, 0, 
        //          0, 0, 0, -1],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,0);
        fp2_sett(&out->M02,0);
        fp2_sett(&out->M03,0);
        fp2_sett(&out->M10,0);
        fp2_sett(&out->M11,1);
        fp2_sett(&out->M12,0);
        fp2_sett(&out->M13,0);
        fp2_sett(&out->M20,0);
        fp2_sett(&out->M21,0);
        fp2_sett(&out->M22,1);
        fp2_sett(&out->M23,0);
        fp2_sett(&out->M30,0);
        fp2_sett(&out->M31,0);
        fp2_sett(&out->M32,0);
        fp2_sett(&out->M33,-1);

    }
    else if (good==7) {
        // (2, 0): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, -1, -1, 1, 
        //          -1, -1, 1, 1],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,1);
        fp2_sett(&out->M02,1);
        fp2_sett(&out->M03,1);
        fp2_sett(&out->M10,1);
        fp2_sett(&out->M11,-1);
        fp2_sett(&out->M12,1);
        fp2_sett(&out->M13,-1);
        fp2_sett(&out->M20,1);
        fp2_sett(&out->M21,-1);
        fp2_sett(&out->M22,-1);
        fp2_sett(&out->M23,1);
        fp2_sett(&out->M30,-1);
        fp2_sett(&out->M31,-1);
        fp2_sett(&out->M32,1);
        fp2_sett(&out->M33,1);
    }
    else if (good==8) {
        //(2, 1): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, -1, -1, 1, 
        //          1, 1, -1, -1],

        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,1);
        fp2_sett(&out->M02,1);
        fp2_sett(&out->M03,1);
        fp2_sett(&out->M10,1);
        fp2_sett(&out->M11,-1);
        fp2_sett(&out->M12,1);
        fp2_sett(&out->M13,-1);
        fp2_sett(&out->M20,1);
        fp2_sett(&out->M21,-1);
        fp2_sett(&out->M22,-1);
        fp2_sett(&out->M23,1);
        fp2_sett(&out->M30,1);
        fp2_sett(&out->M31,1);
        fp2_sett(&out->M32,-1);
        fp2_sett(&out->M33,-1);
    }
    else if (good==5){
        // (1, 0): [1, 1, 1, 1, 
        //          1, -1, -1, 1, 
        //          1, 1, -1, -1, 
        //          -1, 1, -1, 1],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,1);
        fp2_sett(&out->M02,1);
        fp2_sett(&out->M03,1);
        fp2_sett(&out->M10,1);
        fp2_sett(&out->M11,-1);
        fp2_sett(&out->M12,-1);
        fp2_sett(&out->M13,1);
        fp2_sett(&out->M20,1);
        fp2_sett(&out->M21,1);
        fp2_sett(&out->M22,-1);
        fp2_sett(&out->M23,-1);
        fp2_sett(&out->M30,-1);
        fp2_sett(&out->M31,1);
        fp2_sett(&out->M32,-1);
        fp2_sett(&out->M33,1);

    }
    else if (good==6) {
        //(1, 2): [1, 0, 0, 0, 
        //          0, 1, 0, 0, 
        //          0, 0, 0, 1, 
        //          0, 0, 1, 0],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,0);
        fp2_sett(&out->M02,0);
        fp2_sett(&out->M03,0);
        fp2_sett(&out->M10,0);
        fp2_sett(&out->M11,1);
        fp2_sett(&out->M12,0);
        fp2_sett(&out->M13,0);
        fp2_sett(&out->M20,0);
        fp2_sett(&out->M21,0);
        fp2_sett(&out->M22,0);
        fp2_sett(&out->M23,1);
        fp2_sett(&out->M30,0);
        fp2_sett(&out->M31,0);
        fp2_sett(&out->M32,1);
        fp2_sett(&out->M33,0);

    }
    else if (good==9) {
        // (3, 0): [1, 1, 1, 1, 
        //          1, -1, 1, -1, 
        //          1, 1, -1, -1, 
        //          -1, 1, 1, -1],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,1);
        fp2_sett(&out->M02,1);
        fp2_sett(&out->M03,1);
        fp2_sett(&out->M10,1);
        fp2_sett(&out->M11,-1);
        fp2_sett(&out->M12,1);
        fp2_sett(&out->M13,-1);
        fp2_sett(&out->M20,1);
        fp2_sett(&out->M21,1);
        fp2_sett(&out->M22,-1);
        fp2_sett(&out->M23,-1);
        fp2_sett(&out->M30,-1);
        fp2_sett(&out->M31,1);
        fp2_sett(&out->M32,1);
        fp2_sett(&out->M33,-1);
    }
    else if (good==10) {
        // (3, 3): [1, 0, 0, 0, 
        //          0, 1, 0, 0, 
        //          0, 0, 1, 0,
        //          0, 0, 0, 1],
        fp2_sett(&out->M00,1);
        fp2_sett(&out->M01,0);
        fp2_sett(&out->M02,0);
        fp2_sett(&out->M03,0);
        fp2_sett(&out->M10,0);
        fp2_sett(&out->M11,1);
        fp2_sett(&out->M12,0);
        fp2_sett(&out->M13,0);
        fp2_sett(&out->M20,0);
        fp2_sett(&out->M21,0);
        fp2_sett(&out->M22,1);
        fp2_sett(&out->M23,0);
        fp2_sett(&out->M30,0);
        fp2_sett(&out->M31,0);
        fp2_sett(&out->M32,0);
        fp2_sett(&out->M33,1);

    }

    // now we apply the isomorphism if it was computed
    if (good) {
        apply_isomorphism(&out->B.null_point,out,&A->null_point);
    }
    return good;
}


void theta_product_structure_to_elliptic_product(theta_couple_curve_t *E12, theta_structure_t *A) {
    fp2_t xx,yy,temp1,temp2;

    printf("\n");
    fp2_print("a",A->null_point.x);
    fp2_print("b",A->null_point.y);
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
    fp2_sett(&E12->E1.C,1);
    fp2_print("A",E12->E1.A);
    printf("\n");
    fp2_print("a",A->null_point.y);
    fp2_print("b",A->null_point.t);

    // same with y,t
    fp2_sqr(&xx,&A->null_point.y);
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
    fp2_sett(&E12->E2.C,1);
    fp2_print("A",E12->E2.A);
}

/**
 * @brief Compute  a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the theta_chain
 * @param n : the length of the isogeny chain
 * @param E12 an elliptic curve product 
 * @param T1 a couple point on E12[2^(n+2)]
 * @param T2 a couple point on E12[2^(n+2)]
 * @param T1m2 a couple point on E12[2^(n+2)] equal to T1-T2
 *   
 * out : E1xE2 -> E3xE4 of kernel [4](T1,T2) 
 *  
   */
void theta_chain_comput(theta_chain_t *out,int n,const theta_couple_curve_t *E12,const theta_couple_point_t *T1,const theta_couple_point_t *T2, const theta_couple_point_t *T1m2) {

    theta_couple_point_t P1,P2,P1m2;
    theta_point_t Q1,Q2,R1,R2;
    theta_isogeny_t steps[n-1];
    theta_structure_t codomain;

    ibz_t a,b;
    ibz_init(&a);
    ibz_init(&b);

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
    double_couple_point_iter(&P1m2,n-1,E12,T1m2);

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

    printf("\n");
    theta_print("T0",out->first_step.codomain);


    // set-up the theta_structure for the first codomain 
    codomain.null_point=out->first_step.codomain;
    codomain.precomputation=0;
    theta_precomputation(&codomain);


    // push the kernel through the gluing isogeny
    // need to setup the input before
    ibz_pow(&a,&ibz_const_two,n);
    ibz_set(&b,0);
    gluing_eval_basis(&Q1,&Q2,T1,T2,T1m2,&a,&b,E12,&out->first_step);
    
    printf("\n");
    theta_print("phi0(K1)",Q1);
    theta_print("phi0(K2)",Q2);

    for (int i=0;i<n-1;i++) {

        // computing the kernel of the next step
        double_iter(&R1,&codomain,&Q1,n-i-2);
        double_iter(&R2,&codomain,&Q2,n-i-2);
        if ((i < 4)|(i > 114)) {
            printf("\nK1_%d",i+1);
            theta_print("",R1);
            printf("K2_%d",i+1);
            theta_print("",R2);
        }
    
        // computing the next step
        if (i==n-3) {
            printf("avant-dernier %d \n",i+1);
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,0,0);
        }
        else if (i==n-2) {
            printf("dernier %d \n",i+1);
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,1,0);
        }
        else {
            theta_isogeny_comput(&steps[i],&codomain,&R1,&R2,0,1);
        }

        // updating the codomain
        codomain=steps[i].codomain;
        if ((i < 4)|(i > 114)) {
            printf("\n");
            printf("T%d",i+1);
            theta_print("",codomain.null_point);
            fp2_print("precomp",steps[i].precomputation.z);
            
        }
        

        // pushing the kernel
        if (i< n-2) {
            theta_isogeny_eval(&Q1,&steps[i],&Q1);
            theta_isogeny_eval(&Q2,&steps[i],&Q2);
        }
        if ((i < 4)|(i > 114 && i<n-2)) {
            printf("\nphi%d",i+1);
            theta_print("(K1)",Q1);
            printf("phi%d",i+1);
            theta_print("(K2)",Q2);
        }
    
    }

    // copying the steps
    out->steps=steps;

    printf("\n");
    fp2_print("x",steps[n-2].codomain.null_point.x);
    fp2_print("y",steps[n-2].codomain.null_point.y);
    fp2_print("z",steps[n-2].codomain.null_point.z);
    fp2_print("t",steps[n-2].codomain.null_point.t);


    printf("\n");
    //final splitting step
    int is_split = splitting_comput(&out->last_step,&steps[n-2].codomain);
    printf("\n");
    // theta_print("T_fin",out->last_step.B.null_point);
    // fp2_t inv_x;
    // fp2_copy(&inv_x,&out->last_step.B.null_point.x);
    // fp2_inv(&inv_x);
    // fp2_mul(&out->last_step.B.null_point.x,&out->last_step.B.null_point.x,&inv_x);
    // fp2_mul(&out->last_step.B.null_point.y,&out->last_step.B.null_point.y,&inv_x);
    // fp2_mul(&out->last_step.B.null_point.z,&out->last_step.B.null_point.z,&inv_x);
    // fp2_mul(&out->last_step.B.null_point.t,&out->last_step.B.null_point.t,&inv_x);
    fp2_print("x",out->last_step.B.null_point.x);
    fp2_print("y",out->last_step.B.null_point.y);
    fp2_print("z",out->last_step.B.null_point.z);
    fp2_print("t",out->last_step.B.null_point.t);
    assert(is_split);

    // computing the curves of the codomain
    theta_product_structure_to_elliptic_product(&out->codomain,&out->last_step.B);

    fp2_t j2,j3;
    ec_j_inv(&j2,&out->codomain.E1);
    ec_j_inv(&j3,&out->codomain.E2);

    printf("\n");

    fp2_print("j2",j2);
    fp2_print("j3",j3);

    ibz_finalize(&a);
    ibz_finalize(&b);

}


void theta_point_to_montgomery_point(theta_couple_point_t *P12, const theta_point_t *P, const theta_structure_t *A){

    fp2_t temp;

    // P1.X = A.null_point.y * P.x + A.null_point.x * P.y
    // P1.Z = - A.null_point.y * P.x + A.null_point.x * P.y
    fp2_mul(&P12->P1.x,&A->null_point.y,&P->x);
    fp2_mul(&temp,&A->null_point.x,&P->y);
    fp2_sub(&P12->P1.z,&temp,&P12->P1.x);
    fp2_add(&P12->P1.x,&P12->P1.x,&temp);

    // P2.X = A.null_point.t * P.y + A.null_point.y * P.t
    // P2.Z = A.null_point.t * P.y - A.null_point.y * P.t
    fp2_mul(&P12->P2.x,&A->null_point.t,&P->y);
    fp2_mul(&temp,&A->null_point.y,&P->t);
    fp2_sub(&P12->P2.z,&temp,&P12->P2.x);
    fp2_add(&P12->P2.x,&P12->P2.x,&temp);
    

}

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple point on E12, 
 * @param Help a couple point on E12 
 *   
 * phi : E1xE2 -> E3xE4 of kernel 
 * P12 in E1xE2
 * out = phi(P12) in E3xE4 
 * Help is equal to phi.first_step.K1_4 + P12
 *  
   */
void theta_chain_eval(theta_couple_point_t *out,theta_chain_t *phi,const theta_couple_point_t *P12,const theta_couple_point_t *Help) {
    
    theta_point_t temp;

    // first, we apply the gluing
    gluing_eval_point(&temp,P12,Help,&phi->first_step);
    printf("\n");
    theta_print("P0",temp);

    // then, we apply the successive isogenies
    for (int i=0;i<phi->length-1;i++) {
        theta_isogeny_eval(&temp,&phi->steps[i],&temp);
        if ((i<4)||(i>118) ) {
            printf("\nP%d",i+1);
            theta_print("",temp);
        }
        
    }


    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp,&phi->last_step,&temp);
    printf("\nP%d",phi->length);
    theta_print("",temp);

    // finaly the send the result to the elliptic product
    theta_point_to_montgomery_point(out,&temp,&phi->last_step.B);

    printf("\n");
    point_print("F1",out->P1);
    point_print("F2",out->P2);

}