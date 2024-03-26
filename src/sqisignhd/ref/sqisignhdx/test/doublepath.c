
#include <ec.h>
#include <quaternion.h>
#include <sqisignhd.h>
#include <inttypes.h>
#include "test_sqisignhd.h"



bool curve_is_canonical(ec_curve_t const *E);   // test_sqisign.c

int test_doublepath()
{
    int res = 1;

    quat_left_ideal_t ideal;
    quat_left_ideal_init(&ideal);

    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;

    res &= 1;

    quat_left_ideal_finalize(&ideal);

    return res;
}

