/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief The HD-isogenies algorithm required by the signature
 *   
 */

#ifndef HD_H
#define HD_H 

#include <ec.h>
#include <torsion_constants.h> 
#include <stdio.h>

/** @defgroup hd_iso_types Types for dim 2 isogenies   
 * @{
*/

/** @brief Type for product of an elliptic product
 * 
 * @typedef hd_product_t
 * 
 * Represented as two ec_curve_t
*/
typedef struct hd_product {
    ec_curve_t E1; ///< the domain of the isogeny
    ec_curve_t E2; ///< the codomain of the isogeny
} hd_product_t;

/** @brief Type for kernel of an isogeny between elliptic product
 * 
 * @typedef hd_kernel_t
 * 
 * Represented as two ec_point_t
*/
typedef struct hd_kernel {
    ec_point_t E1; ///< the first point generating the kernel of the isogeny
    ec_point_t E2; ///< the second point generating the kernel of the isogeny
} hd_kernel_t;

/** @brief Type for chain of two isogenies
 * 
 * @typedef hd_two_isog
 * 
 * Represented as a vector ec_isog_even_t
*/
typedef struct hd_two_isog {
    unsigned short length; ///< the length of the 2-isogeny chain
    hd_kernel_t kernel; ///< the kernel of the isogeny
    hd_product_t domain; ///< the domain of the isogeny (which must be an elliptic product)
} hd_two_isog_t;


/** @}
*/


/*************************** Functions *****************************/


#endif