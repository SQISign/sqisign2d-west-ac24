#include <ec.h> 
#include <fp2.h>

/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief the theta structure header 
 */

#ifndef THETA_STRUCTURE_H
#define THETA_STRUCTURE_H

/** @brief Type for theta point * 
 * @typedef theta_point_t
 *
 * @struct theta_point
 * 
 * the  theta_point structure used  
*/
typedef struct theta_point {  
    fp2_t x;
    fp2_t y;
    fp2_t z;
    fp2_t t;
} theta_point_t;


/** @brief Type for theta structure * 
 * @typedef theta_structure_t
 *
 * @struct theta_structure
 * 
 * the  theta_structure structure used  
*/
typedef struct theta_structure {
    theta_point_t null_point;
    fp2_t y0;
    fp2_t z0;
    fp2_t t0;
    fp2_t Y0;
    fp2_t Z0;
    fp2_t T0;
    
} theta_structure_t;


/*************************** Functions *****************************/




/**
 * @brief Perform the hadamard transform on a theta point
 *
 * @param out Output: the theta_point 
 * @param in a theta point*  
 * in = (x,y,z,t)
 * out = (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)
 *  
   */
void hadamard(theta_point_t *out,const theta_point_t *in);

/**
 * @brief Square the coordinates and then perform the hadamard transform
 *
 * @param out Output: the theta_point 
 * @param in a theta point*  
 * in = (x,y,z,t)
 * out = (x^2+y^2+z^2+t^2, x^2-y^2+z^2-t^2, x^2+y^2-z^2-t^2, x^2-y^2-z^2+t^2)
 *  
   */
void to_squared_theta(theta_point_t *out,const theta_point_t *in);

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
void theta_precomputation(theta_structure_t *A);

/**
 * @brief Compute the double of the theta point in on the theta struc A
 *
 * @param out Output: the theta_point 
 * @param A a theta structure
 * @param in a theta point in the theta structure A   
 * in = (x,y,z,t)
 * out = [2] (x,y,z,t)
 * /!\ assumes that no coordinates is zero and that the precomputation of A has been done
 *  
   */
void double_point(theta_point_t *out,const theta_structure_t *A,const theta_point_t *in);

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
void double_iter(theta_point_t *out,const theta_structure_t *A,const theta_point_t *in, int exp);

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
void diff_add(theta_point_t *R,const theta_structure_t *A,const theta_point_t *P,const theta_point_t *Q,const theta_point_t *PQ);


//// TODO: see if this is needed (maybe not if we only handle powers of 2)
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
// void mul(theta_point_t* R, theta_structure_t const* A, digit_t const* k, theta_point_t const* P);

#endif