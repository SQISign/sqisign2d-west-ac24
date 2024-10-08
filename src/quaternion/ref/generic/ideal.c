#include <quaternion.h>
#include <stdlib.h>
#include "internal.h"

void
quat_lideal_create_principal(quat_left_ideal_t *lideal,
                             const quat_alg_elem_t *x,
                             const quat_order_t *order,
                             const quat_alg_t *alg)
{
    ibz_mat_4x4_t mulmat;
    ibq_t norm;
    ibz_mat_4x4_init(&mulmat);
    ibq_init(&norm);

    // Multiply order basis on the right by x
    quat_alg_rightmul_mat(&mulmat, x, alg);
    ibz_mat_4x4_mul(&lideal->lattice.basis, &mulmat, &order->basis);

    // Adjust denominators
    ibz_mul(&lideal->lattice.denom, &x->denom, &order->denom);
    quat_lattice_reduce_denom(&lideal->lattice, &lideal->lattice);

    // Compute HNF
    quat_lattice_hnf(&lideal->lattice);

    // Compute norm and check it's integral
    quat_alg_norm(&norm, x, alg);
    assert(ibq_is_ibz(&norm));
    ibq_to_ibz(&lideal->norm, &norm);

    // Set order
    lideal->parent_order = order;

    ibq_finalize(&norm);
    ibz_mat_4x4_finalize(&mulmat);
}

void
quat_lideal_create_from_primitive(quat_left_ideal_t *lideal,
                                  const quat_alg_elem_t *x,
                                  const ibz_t *N,
                                  const quat_order_t *order,
                                  const quat_alg_t *alg)
{
    quat_lattice_t ON;
    quat_lattice_init(&ON);

    // Compute ideal generated by x
    quat_lideal_create_principal(lideal, x, order, alg);

    // Compute norm
    ibz_gcd(&lideal->norm, &lideal->norm, N);

    // Compute ideal generated by N (without reducing denominator)
    ibz_mat_4x4_scalar_mul(&ON.basis, N, &order->basis);
    ibz_copy(&ON.denom, &order->denom);

    // Add lattices (reduces denominators)
    quat_lattice_add(&lideal->lattice, &lideal->lattice, &ON);

#ifndef NDEBUG
    // verify norm
    ibz_t tmp;
    ibz_init(&tmp);
    quat_lattice_index(&tmp, &(lideal->lattice), (lideal->parent_order));
    int ok = ibz_sqrt(&tmp, &tmp);
    assert(ok);
    ok = (0 == ibz_cmp(&tmp, &(lideal->norm)));
    assert(ok);
    ibz_finalize(&tmp);
#endif

    // Set order
    lideal->parent_order = order;

    quat_lattice_finalize(&ON);
}

void
quat_lideal_make_primitive_then_create(quat_left_ideal_t *lideal,
                                       const quat_alg_elem_t *x,
                                       const ibz_t *N,
                                       const quat_order_t *order,
                                       const quat_alg_t *alg)
{
    quat_alg_elem_t prim;
    ibz_t imprim, n1;
    quat_alg_elem_init(&prim);
    ibz_init(&imprim);
    ibz_init(&n1);

    // Store the primitive part of x in elem
    quat_alg_make_primitive(&prim.coord, &imprim, x, order, alg);
    ibz_mat_4x4_eval(&prim.coord, &order->basis, &prim.coord);
    ibz_copy(&prim.denom, &order->denom);

    // imprim = gcd(imprimitive part of x, N)
    // n1 = N / imprim
    ibz_gcd(&imprim, &imprim, N);
    ibz_div(&n1, &imprim, N, &imprim);

    // Generate the ideal (elem, n1)
    quat_lideal_create_from_primitive(lideal, &prim, &n1, order, alg);

    quat_alg_elem_finalize(&prim);
    ibz_finalize(&imprim);
    ibz_finalize(&n1);
}

int
quat_lideal_random_2e(quat_left_ideal_t *lideal,
                      const quat_order_t *order,
                      const quat_alg_t *alg,
                      int64_t e,
                      unsigned char n)
{
    ibq_t norm;
    ibz_t norm_num;
    ibq_init(&norm);
    ibz_init(&norm_num);

    quat_alg_elem_t gen;
    quat_alg_elem_init(&gen);
    ibz_set(&gen.coord[0], 1);

    // Start with the trivial left ideal of O
    quat_lideal_create_principal(lideal, &gen, order, alg);

    // Create the lattice 2·O
    quat_lattice_t O2;
    quat_lattice_init(&O2);
    if (ibz_get(&order->denom) % 2 == 0) {
        ibz_mat_4x4_copy(&O2.basis, &order->basis);
        ibz_div_2exp(&O2.denom, &order->denom, 1);
    } else {
        ibz_mat_4x4_scalar_mul(&O2.basis, &ibz_const_two, &order->basis);
        ibz_copy(&O2.denom, &order->denom);
    }

    for (int i = 0; i < e; ++i) {
        // Take random gen ∈ lideal until one is found such that
        // val₂(N(gen)) > i and gen ∉ 2·O
        do {
            if (!quat_lattice_random_elem(&gen, &lideal->lattice, n))
                return 0;

            quat_alg_norm(&norm, &gen, alg);
            int ok = ibq_to_ibz(&norm_num, &norm);
            assert(ok);
            // N(gen)/N(I)
            ibz_div_2exp(&norm_num, &norm_num, i);
        }
        // Check that 2-norm has increased, do not backtrack
        while ((ibz_get(&norm_num) % 2 != 0) || quat_lattice_contains(NULL, &O2, &gen, alg));

        // Redefine lideal as (gen, 2^(i+1))
        ibz_mul(&norm_num, &lideal->norm, &ibz_const_two);
        quat_lideal_create_from_primitive(lideal, &gen, &norm_num, order, alg);
        assert(ibz_cmp(&lideal->norm, &norm_num) == 0);
    }

    quat_lattice_finalize(&O2);
    ibq_finalize(&norm);
    ibz_finalize(&norm_num);
    quat_alg_elem_finalize(&gen);
    return 1;
}

int
quat_lideal_generator(quat_alg_elem_t *gen,
                      const quat_left_ideal_t *lideal,
                      const quat_alg_t *alg,
                      int bound)
{
    return (quat_lideal_generator_coprime(gen, lideal, &ibz_const_one, alg, bound));
}

int
quat_lideal_generator_coprime(quat_alg_elem_t *gen,
                              const quat_left_ideal_t *lideal,
                              const ibz_t *n,
                              const quat_alg_t *alg,
                              int bound)
{
    ibq_t norm;
    ibz_t norm_int, norm_n, gcd, r, q, n2, gcd_ref;
    ibz_vec_4_t vec;
    ibz_vec_4_init(&vec);
    ibq_init(&norm);
    ibz_init(&norm_int);
    ibz_init(&norm_n);
    ibz_init(&gcd_ref);
    ibz_init(&n2);
    ibz_init(&r);
    ibz_init(&q);
    ibz_init(&gcd);
    ibz_mul(&norm_n, &lideal->norm, n);
    ibz_mul(&n2, n, n);
    ibz_gcd(&gcd_ref, n, &(lideal->norm));
    int a, b, c, d, int_norm;
    int used_bound = QUATERNION_lideal_generator_search_bound;
    if (bound != 0)
        used_bound = bound;
    assert(used_bound > 0);
    int found = 0;
    for (int_norm = 1; int_norm < used_bound; int_norm++) {
        for (a = -int_norm; a <= int_norm; a++) {
            for (b = -int_norm + abs(a); b <= int_norm - abs(a); b++) {
                for (c = -int_norm + abs(a) + abs(b); c <= int_norm - abs(a) - abs(b); c++) {
                    d = int_norm - abs(a) - abs(b) - abs(c);
                    ibz_vec_4_set(&vec, a, b, c, d);
                    ibz_content(&gcd, &vec);
                    if (ibz_is_one(&gcd)) {
                        ibz_mat_4x4_eval(&(gen->coord), &(lideal->lattice.basis), &vec);
                        ibz_copy(&(gen->denom), &(lideal->lattice.denom));
                        quat_alg_norm(&norm, gen, alg);
                        found = ibq_to_ibz(&norm_int, &norm);
                        if (found) {
                            ibz_div(&q, &r, &norm_int, &(lideal->norm));
                            found = ibz_is_zero(&r);
                            if (found) {
                                ibz_gcd(&gcd, &norm_n, &q);
                                found = (0 == ibz_cmp(&gcd, &ibz_const_one));
                                if (found) {
                                    ibz_gcd(&gcd, &n2, &norm_int);
                                    found = (0 == ibz_cmp(&gcd, &gcd_ref));
                                }
                            }
                        }
                        if (found)
                            goto fin;
                    }
                }
            }
        }
    }
fin:;
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibq_finalize(&norm);
    ibz_finalize(&norm_int);
    ibz_finalize(&norm_n);
    ibz_finalize(&gcd_ref);
    ibz_finalize(&n2);
    ibz_vec_4_finalize(&vec);
    ibz_finalize(&gcd);
    return (found);
}

int
quat_lideal_mul(quat_left_ideal_t *product,
                const quat_left_ideal_t *lideal,
                const quat_alg_elem_t *alpha,
                const quat_alg_t *alg,
                int bound)
{
    ibq_t norm, norm_lideal;
    ibz_t norm_int, num, denom;
    quat_alg_elem_t gen;
    ibq_init(&norm);
    ibq_init(&norm_lideal);
    ibz_init(&norm_int);
    ibz_init(&num);
    ibz_init(&denom);
    quat_alg_elem_init(&gen);

    quat_alg_norm(&norm, alpha, alg);
    ibq_num(&num, &norm);
    ibq_denom(&denom, &norm);
    ibz_gcd(&denom, &denom, &num);
    ibz_div(&norm_int, &denom, &num, &denom);
    // Find a random generator gen of norm coprime to N(alpha)
    // Warning: hardcoded the default constant, and this call can fail!
    int found = quat_lideal_generator_coprime(&gen, lideal, &norm_int, alg, 0);
    if (found) {
        // Define the ideal (gen·α, N(lideal)·N(α))
        quat_alg_mul(&gen, &gen, alpha, alg);
        ibz_copy(&num, &(lideal->norm));
        ibz_set(&denom, 1);
        ibq_set(&norm_lideal, &num, &denom);
        ibq_mul(&norm, &norm, &norm_lideal);
        int ok = ibq_to_ibz(&norm_int, &norm);
        assert(ok);
        quat_lideal_create_from_primitive(product, &gen, &norm_int, lideal->parent_order, alg);
    }

    ibq_finalize(&norm);
    ibz_finalize(&norm_int);
    ibq_finalize(&norm_lideal);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    quat_alg_elem_finalize(&gen);
    return (found);
}

void
quat_lideal_add(quat_left_ideal_t *sum,
                const quat_left_ideal_t *I1,
                const quat_left_ideal_t *I2,
                const quat_alg_t *alg)
{
    assert(I1->parent_order == I2->parent_order);
    quat_lattice_add(&sum->lattice, &I1->lattice, &I2->lattice);
    quat_lattice_index(&sum->norm, &sum->lattice, I1->parent_order);
    int ok = ibz_sqrt(&sum->norm, &sum->norm);
    assert(ok);
    sum->parent_order = I1->parent_order;
}

void
quat_lideal_inter(quat_left_ideal_t *inter,
                  const quat_left_ideal_t *I1,
                  const quat_left_ideal_t *I2,
                  const quat_alg_t *alg)
{
    assert(I1->parent_order == I2->parent_order);
    quat_lattice_intersect(&inter->lattice, &I1->lattice, &I2->lattice);
    quat_lattice_index(&inter->norm, &inter->lattice, I1->parent_order);
    int ok = ibz_sqrt(&inter->norm, &inter->norm);
    assert(ok);
    inter->parent_order = I1->parent_order;
}

int
quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg)
{
    return I1->parent_order == I2->parent_order && ibz_cmp(&I1->norm, &I2->norm) == 0 &&
           quat_lattice_equal(&I1->lattice, &I2->lattice);
}

int
quat_lideal_isom(quat_alg_elem_t *iso,
                 const quat_left_ideal_t *I1,
                 const quat_left_ideal_t *I2,
                 const quat_alg_t *alg)
{
    // Only accept strict equality of parent orders
    if (I1->parent_order != I2->parent_order)
        return 0;

    quat_lattice_t trans;
    ibz_mat_4x4_t lll;
    ibq_t norm, norm_bound;
    quat_lattice_init(&trans);
    ibz_mat_4x4_init(&lll);
    ibq_init(&norm_bound);
    ibq_init(&norm);

    quat_lattice_right_transporter(&trans, &I1->lattice, &I2->lattice, alg);

    int err = quat_lattice_lll(&lll, &trans, &alg->p);
    assert(!err);
    // The shortest vector found by lll is a candidate
    for (int i = 0; i < 4; i++)
        ibz_copy(&iso->coord[i], &lll[i][0]);
    ibz_copy(&iso->denom, &trans.denom);

    // bound on the smallest vector in transporter: N₂ / N₁
    ibq_set(&norm_bound, &I2->norm, &I1->norm);
    quat_alg_norm(&norm, iso, alg);

    int isom = ibq_cmp(&norm, &norm_bound) <= 0;

    quat_lattice_finalize(&trans);
    ibz_mat_4x4_finalize(&lll);
    ibq_finalize(&norm_bound);
    ibq_finalize(&norm);

    return isom;
}

void
quat_lideal_right_order(quat_order_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    quat_lattice_right_transporter(order, &lideal->lattice, &lideal->lattice, alg);
}

void
quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced,
                         ibz_mat_4x4_t *gram,
                         const quat_left_ideal_t *lideal,
                         const quat_alg_t *alg)
{
    ibz_mat_4x4_t prod;
    ibz_mat_4x4_init(&prod);
    quat_lattice_lll(reduced, &(lideal->lattice), &(alg->p));
    ibz_mat_4x4_transpose(&prod, reduced);
    ibz_mat_4x4_mul(&prod, &prod, &(alg->gram));
    ibz_mat_4x4_mul(gram, &prod, reduced);
    ibz_mat_4x4_finalize(&prod);
}

/***************************** Function from quaternion_tools.c
 * ***************************************/

void
quat_connecting_ideal(quat_left_ideal_t *connecting_ideal,
                      const quat_order_t *O1,
                      const quat_order_t *O2,
                      const quat_alg_t *alg)
{
    quat_lattice_t inter;
    ibz_t norm, one;
    ibz_mat_4x4_t gens;
    quat_alg_elem_t gen;
    quat_left_ideal_t I[4];
    ibz_init(&norm);
    ibz_init(&one);
    quat_lattice_init(&inter);
    ibz_mat_4x4_init(&gens);
    quat_alg_elem_init(&gen);
    for (int i = 0; i < 4; i++) {
        quat_left_ideal_init(I + i);
    }
    ibz_set(&one, 1);

    quat_lattice_intersect(&inter, O1, O2);
    quat_lattice_index(&norm, &inter, O1);
    ibz_mat_4x4_scalar_mul(&gens, &norm, &(O2->basis));
    quat_alg_scalar(&gen, &norm, &one);
    quat_lideal_create_principal(connecting_ideal, &gen, O1, alg);
    for (int i = 0; i < 4; i++) {
        ;
        quat_alg_elem_copy_ibz(
            &gen, &(O2->denom), &(gens[0][i]), &(gens[1][i]), &(gens[2][i]), &(gens[3][i]));
        quat_lideal_create_principal(I + i, &gen, O1, alg);
    }

    quat_lideal_add(connecting_ideal, connecting_ideal, I + 0, alg);
    quat_lideal_add(connecting_ideal, connecting_ideal, I + 1, alg);
    quat_lideal_add(connecting_ideal, connecting_ideal, I + 2, alg);
    quat_lideal_add(connecting_ideal, connecting_ideal, I + 3, alg);

    ibz_finalize(&norm);
    ibz_finalize(&one);
    quat_lattice_finalize(&inter);
    ibz_mat_4x4_finalize(&gens);
    quat_alg_elem_finalize(&gen);
    for (int i = 0; i < 4; i++) {
        quat_left_ideal_finalize(I + i);
    }
}
