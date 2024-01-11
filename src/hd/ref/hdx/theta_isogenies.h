#include <ec.h> 
#include <fp2.h>
#include "theta_structure.h"

/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief the theta isogeny header 
 */

#ifndef THETA_ISOGENY_H
#define THETA_ISOGENY_H




/** @brief Type for gluing (2,2) theta isogeny * 
 * @typedef theta_gluing_t
 *
 * @struct theta_gluing
 * 
 * the  theta_gluing structure  
*/
typedef struct theta_gluing {  

    theta_couple_point_t K1_4;
    theta_couple_point_t K2_4;
    theta_point_t T1_8;
    theta_point_t T2_8;

    fp2_t M00;
    fp2_t M01;
    fp2_t M02;
    fp2_t M03;

    fp2_t M10;
    fp2_t M11;
    fp2_t M12;
    fp2_t M13;

    fp2_t M20;
    fp2_t M21;
    fp2_t M22;
    fp2_t M23;

    fp2_t M30;
    fp2_t M31;
    fp2_t M32;
    fp2_t M33;
    int32_t zero_idx;
    theta_point_t precomputation;
    theta_point_t codomain;
} theta_gluing_t;


/** @brief Type for standard (2,2) theta isogeny * 
 * @typedef theta_isogeny_t
 *
 * @struct theta_isogeny
 * 
 * the  theta_isogeny structure  
*/
typedef struct theta_isogeny {  
    theta_point_t T1_8;
    theta_point_t T2_8;
    int bool1;
    int bool2;
    theta_structure_t domain;
    theta_point_t precomputation;
    theta_structure_t codomain;
} theta_isogeny_t;


/** @brief Type for chain of (2,2) theta isogenies with preimage by 4 multiplication above * 
 * @typedef theta_chain_t
 *
 * @struct theta_chain
 * 
 * the  theta_chain structure  
*/
typedef struct theta_chain {  
    theta_couple_point_t T1;
    theta_couple_point_t T2;
    int length;
    theta_couple_curve_t domain;
    theta_couple_curve_t codomain;
    theta_gluing_t first_step;
    theta_isogeny_t *steps;
} theta_chain_t;


/*************************** Functions *****************************/

/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing 
 * @param E12 an elliptic curve couple E1 x E2
 * @param K1_8 a point in E1xE2[8]
 * @param K2_8 a point in E1xE2[8]  
 * 
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8) 
 *  
   */
void gluing_comput(theta_gluing_t *out,const theta_couple_curve_t *E12,const theta_couple_point_t *K1_8,const theta_couple_point_t *K2_8);


/**
 * @brief Evaluate a gluing isogeny from an elliptic product
 *
 * @param image Output: the theta_point of the image 
 * @param P : a couple point in E12 
 * @param E12 : an elliptic product
 * @param phi : an isogeny E1 x E2 -> A 
 * 
 * out : phi( P1,P2 ) 
 *  
   */
void gluing_eval(theta_point_t *image,const theta_couple_point_t *P,const theta_couple_curve_t *E12, const theta_gluing_t *phi);

/**
 * @brief Compute  a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny 
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
void theta_isogeny_comput(theta_isogeny_t *out,const theta_structure_t *A,const theta_point_t *T1_8,const theta_point_t *T2_8,int bool1, int bool2);

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
void theta_isogeny_eval(theta_point_t *out,const theta_isogeny_t *phi,const theta_point_t *P);

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
void theta_chain_comput(theta_chain_t *out,int n,const theta_couple_curve_t *E12,const theta_couple_point_t *T1,const theta_couple_point_t *T2);

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
void theta_chain_eval(theta_couple_point_t *out,theta_chain_t *phi,theta_couple_point_t *P12);

#endif