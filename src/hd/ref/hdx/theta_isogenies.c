#include "theta_isogenies.h"


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

    // K1_4 = [2] K1_8  and K2_4 = [2] K2_8
    double_couple_point(&out->K1_4,E12,K1_8);
    double_couple_point(&out->K2_4,E12,K2_8);

}