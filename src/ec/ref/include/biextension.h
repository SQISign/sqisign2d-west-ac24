#ifndef _BIEXT_H_
#define _BIEXT_H_

#include <fp2.h>
#include "ec.h"

void A24_from_AC(ec_point_t* A24, ec_point_t const* AC);

void weil(fp2_t* r, int e, ec_point_t* P, ec_point_t* Q, ec_point_t* PQ, ec_point_t* A24);
#endif

