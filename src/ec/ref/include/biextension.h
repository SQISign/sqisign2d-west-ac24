#ifndef _BIEXT_H_
#define _BIEXT_H_

#include <fp2.h>
#include <ec.h>

void non_reduced_tate(fp2_t* r, int e, ec_point_t const* P_, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24_)
void weil_n(fp2_t* r, int e, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24)
void weil(fp2_t* r, int e, ec_point_t const* P_, ec_point_t const* Q_, ec_point_t const* PQ, ec_point_t const* A24_)
#endif

