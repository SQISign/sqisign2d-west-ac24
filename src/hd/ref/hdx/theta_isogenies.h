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

    theta_structure_t codomain;
} theta_gluing_t;

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
 * @param P1 : a point in E1 
 * @param E1 : an elliptic curve
 * @param P2 : a point in E2
 * @param E2 : an elliptic curve
 * @param phi : an isogeny E1 x E2 -> A 
 * 
 * out : phi( P1,P2 ) 
 *  
   */
void gluing_eval(theta_point_t *image,const ec_point_t *P1,const ec_curve_t *E1,const ec_point_t *P2,const ec_curve_t *E2,theta_gluing_t *phi);


#endif