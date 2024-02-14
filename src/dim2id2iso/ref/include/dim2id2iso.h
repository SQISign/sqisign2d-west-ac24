    /** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief The dim2_id2iso algorithms 
 */

#ifndef DIM2ID2ISO_H
#define DIM2ID2ISO_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt.h>
#include <ec.h>
#include <hd.h>
#include <id2iso.h>



/*************************** Functions *****************************/

/**
 * @brief Computes an arbitrary isogeny of fixed degree starting from E0
 * 
 * @param isog Output : a dim2 isogeny encoding an isogeny of degree u
 * @param lideal Output : an ideal of norm u
 * @param u : integer
 * @returns a bit indicating if the computation succeeded  
 * 
 * F is an isogeny encoding an isogeny phi : E0 -> Eu of degree u
*/
int fixed_degree_isogeny(theta_chain_t *isog, quat_left_ideal_t *lideal, ibz_t *u);



/** @defgroup dim2id2iso_others Other functions needed for id2iso
 * @{
*/

/**
 * @brief Translating an ideal into a representation of the corresponding isogeny
 *
 * @param isog Output : dim 2 isogeny  
 * @param beta1 Output : quaternion element 
 * @param beta2 Output : quaternion element
 * @param u Output : integer
 * @param v Output : integer
 * @param au Output : integer
 * @param bu Output : integer
 * @param phiv Output : dim2 isogeny representation of an isogeny of degree v
 * @param d1 Output : integer 
 * @param d2 Output : integer
 * @param lideal : ideal
 * @param Bpoo : the quaternion algebra
 * @returns a bit indicating if the computation succeeded
 *  
 * beta1 and beta2 are elements in lideal of norm n(lideal)d1 and n(lideal) d2 respectively
 * u,v are integers such that 2^e = d1 u + d2 v and u = au^2 + bu^2 
 * phiv is a dim 2 isogeny representing an isogeny of degree v : E0 -> Ev
 * F is a dim2 2^e - isogeny between E0 x Ev -> E_I x E 
 * that encodes an isogeny E0 -> E_I corresponding to the ideal lideal
 */
int dim2id2iso_ideal_to_isogeny_clapotis(theta_chain_t *isog, quat_alg_elem_t *beta1, quat_alg_elem_t *beta2, ibz_t *u, ibz_t *v, ibz_t *au, ibz_t *bu, theta_chain_t *phiv,ibz_t *d1,ibz_t *d2, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo);



#endif
