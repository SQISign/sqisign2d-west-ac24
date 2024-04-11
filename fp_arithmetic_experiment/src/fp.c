#include "fp.h"

// TODO should we modify fp_sqrt to take an out
// variable?
void 
fp_sqrt(fp_t * x)
{
    (void)gf5248_sqrt(x, x);
}

bool
fp_is_square(const fp_t * a)
{
    return -1 != gf5248_legendre(a);
}

void
fp_inv(fp_t * x)
{
    (void)gf5248_invert(x, x);
}

void
fp_copy(fp_t * out, const fp_t * a)
{
    out->v0 = a->v0;
    out->v1 = a->v1;
    out->v2 = a->v2;
    out->v3 = a->v3;
}

void fp_set_zero(fp_t * a){
    fp_copy(a, &ZERO);
}

void fp_set_one(fp_t * a){
    fp_copy(a, &ONE);
}
