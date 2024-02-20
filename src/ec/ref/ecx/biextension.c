#include "curve_extras.h"
#include "tedwards.h"
#include <ec_params.h>
#include <assert.h>

/*
 * We implement the biextension arithmetic by using the cubical torsor representation. For now only implement the 2^e-ladder.
 *
 * Warning: both cubicalDBL and cubicalADD are off by a factor x4 with respect
 * to the cubical arithmetic. 
 * Since this factor is the same, this means that the biextension
 * arithmetic is correct, so the pairings are ok (they only rely on the
 * biextension arithmetic).
 * In the special case where Q=P (self pairings), we use the cubical ladder
 * rather than the biextension ladder because this is faster. In that case,
 * when we do a ladder we are off by a factor 4^m, m the number of bits.
 * This factor thus disappear in the Weil pairing since we take a quotient,
 * and also in the Tate pairing due to the final exponentiation; so
 * everything is ok too.
 * (Note that when the curves are supersingular as in our case, the Tate
 * self pairing is always trivial anyway because the Galois structure of the
 * isogeneous curves are all the same, so the étale torsor representing the
 * Tate pairing has to be trivial).
 */

// this is exactly like xDBLv2
// Warning: for now we need to assume that A24 is normalised, ie C24=1.
// (maybe add an assert?)
// Anyway, we won't use this function but directly xDBLv2 and xADD
void cubicalDBL(ec_point_t* Q, ec_point_t const* P, ec_point_t const* A24)
{
    // A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

// this would be exactly like xADD if PQ was 'antinormalised' as (1,z)
// Not used yet, we use xADD directly
void cubicalADD(ec_point_t* R, ec_point_t const* P, ec_point_t const* Q, fp2_t const* ixPQ)
{
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&t2, &t2);
    fp2_sqr(&R->z, &t3);
    fp2_mul(&R->x, ixPQ, &t2);
}

void ec_anti_normalize(ec_point_t* P){
    fp2_inv(&P->x);
    fp2_mul(&P->z, &P->z, &P->x);
    fp_mont_setone(P->x.re);
    fp_set(P->x.im, 0);
}

// given cubical reps of P+Q, Q, P, return P+2Q, 2Q
void biextDBL(ec_point_t* R, ec_point_t* S, ec_point_t const* PQ, ec_point_t const* Q, ec_point_t const* P, ec_point_t const* A24)
{
    //cubicalADD(R, PQ, Q, P);
    //cubicalDBL(S, Q, A24);
    xADD(R, PQ, Q, P);
    xDBLv2(S, Q, A24);
}

void biext_ladder_2e(int e, ec_point_t* R, ec_point_t* S, ec_point_t const* PQ, ec_point_t const* Q, ec_point_t const* P, ec_point_t const* A24)
{
    copy_point(R, PQ);
    copy_point(S, Q);
    for (int i=0; i<e; i++) {
        biextDBL(R, S, R, S, P, A24);
    }
}

void ratio(fp2_t* r, ec_point_t const* PnQ, ec_point_t const* nQ, ec_point_t const* P) 
{
    fp2_t t0;
    fp2_mul(r, &nQ->x, &P->x);
    fp2_copy(&t0, &PnQ->x);
    fp2_inv(&t0);
    fp2_mul(r, r, &t0);
}

void translate(ec_point_t* P, ec_point_t const* T)
{
    fp2_t t0, t1;
    if (fp2_is_zero(&T->z)) {
        // do nothing
    }
    else if (fp2_is_zero(&T->x)) {
        fp2_copy(&t0, &P->x);
        fp2_copy(&P->x, &P->z);
        fp2_copy(&P->z, &t0);
    }
    else {
        fp2_mul(&t0, &T->x, &P->x);
        fp2_mul(&t1, &T->z, &P->z);
        fp2_sub(&P->x, &t0, &t1);
        fp2_mul(&t0, &T->z, &P->x);
        fp2_mul(&t1, &T->x, &P->z);
        fp2_sub(&P->z, &t0, &t1);
    }
}

void monodromy(fp2_t* r, int e, ec_point_t const* PQ, ec_point_t const* Q, ec_point_t const* P, ec_point_t const* A24)
{
    ec_point_t R0, R1;
    biext_ladder_2e(e-1, &R0, &R1, PQ, Q, P, A24);
    translate(&R0, &R1);
    translate(&R1, &R1);
    ratio(r, &R0, &R1, P);
}

// TODO: use only one inversion
void inline_to_cubical(ec_point_t* P, ec_point_t* A24) {
    ec_normalize(A24);
    ec_anti_normalize(P);
}

void to_cubical(ec_point_t* P, ec_point_t* A24, ec_point_t const* P_, ec_point_t const* A24_) {
    copy_point(P, P_);
    copy_point(A24, A24_);
    inline_to_cubical(P, A24);
}

// TODO: use only one inversion
void inline_to_cubical_3(ec_point_t* P, ec_point_t* Q, ec_point_t* A24) {
    ec_normalize(A24);
    ec_anti_normalize(P);
    ec_anti_normalize(Q);
}

void to_cubical_3(ec_point_t* P, ec_point_t* Q, ec_point_t* A24, ec_point_t const* P_, ec_point_t const* Q_, ec_point_t const* A24_) {
    copy_point(P, P_);
    copy_point(Q, Q_);
    copy_point(A24, A24_);
    inline_to_cubical_3(P, Q, A24);
}

// non reduced Tate pairing, PQ should be P+Q in (X:Z) coordinates
// If the points are already normalized correctly, use 'monodromy'
void non_reduced_tate(fp2_t* r, int e, ec_point_t const* P_, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24_) {
    ec_point_t P, A24;
    to_cubical(&P, &A24, P_, A24_);
    monodromy(r, e, PQ, Q, &P, &A24);
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
// We assume the points are normalised correctly
void weil_n(fp2_t* r, int e, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24) {
    fp2_t t0;
    monodromy(&t0, e, PQ, Q, P, A24);
    fp2_inv(&t0);
    monodromy(r, e, PQ, P, Q, A24);
    // TODO: check if that's the Weil pairing or its inverse
    fp2_mul(r, r, &t0);
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
void weil(fp2_t* r, int e, ec_point_t const* P_, ec_point_t const* Q_, ec_point_t const* PQ, ec_point_t const* A24_) {
    ec_point_t P, Q, A24;
    to_cubical_3(&P, &Q, &A24, P_, Q_, A24_);
    weil_n(r, e, &P, &Q, PQ, &A24);
}
