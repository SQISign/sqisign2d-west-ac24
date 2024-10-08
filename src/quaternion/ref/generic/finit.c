#include <quaternion.h>

void
quat_alg_init_set(quat_alg_t *alg, const ibz_t *p)
{
    ibz_init(&(*alg).p);
    ibz_mat_4x4_init(&(*alg).gram);
    ibz_copy(&(*alg).p, p);
    ibz_set(&(*alg).gram[0][0], 1);
    ibz_set(&(*alg).gram[1][1], 1);
    ibz_copy(&(*alg).gram[2][2], p);
    ibz_copy(&(*alg).gram[3][3], p);
}
void
quat_alg_finalize(quat_alg_t *alg)
{
    ibz_finalize(&(*alg).p);
    ibz_mat_4x4_finalize(&(*alg).gram);
}

void
quat_alg_elem_init(quat_alg_elem_t *elem)
{
    quat_alg_coord_init(&(*elem).coord);
    ibz_init(&(*elem).denom);
    ibz_set(&(*elem).denom, 1);
}
void
quat_alg_elem_finalize(quat_alg_elem_t *elem)
{
    quat_alg_coord_finalize(&(*elem).coord);
    ibz_finalize(&(*elem).denom);
}

void
quat_alg_coord_init(quat_alg_coord_t *coord)
{
    for (int i = 0; i < 4; i++) {
        ibz_init(&(*coord)[i]);
    }
}
void
quat_alg_coord_finalize(quat_alg_coord_t *coord)
{
    for (int i = 0; i < 4; i++) {
        ibz_finalize(&(*coord)[i]);
    }
}

void
ibz_vec_4_init(ibz_vec_4_t *vec)
{
    for (int i = 0; i < 4; i++) {
        ibz_init(&(*vec)[i]);
    }
}
void
ibz_vec_4_finalize(ibz_vec_4_t *vec)
{
    for (int i = 0; i < 4; i++) {
        ibz_finalize(&(*vec)[i]);
    }
}

void
ibz_vec_5_init(ibz_vec_5_t *vec)
{
    for (int i = 0; i < 5; i++) {
        ibz_init(&(*vec)[i]);
    }
}
void
ibz_vec_5_finalize(ibz_vec_5_t *vec)
{
    for (int i = 0; i < 5; i++) {
        ibz_finalize(&(*vec)[i]);
    }
}

void
ibz_mat_2x2_init(ibz_mat_2x2_t *mat)
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_init(&(*mat)[i][j]);
        }
    }
}
void
ibz_mat_2x2_finalize(ibz_mat_2x2_t *mat)
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_finalize(&(*mat)[i][j]);
        }
    }
}

void
ibz_mat_4x4_init(ibz_mat_4x4_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_init(&(*mat)[i][j]);
        }
    }
}
void
ibz_mat_4x4_finalize(ibz_mat_4x4_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_finalize(&(*mat)[i][j]);
        }
    }
}

void
ibz_mat_4x5_init(ibz_mat_4x5_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            ibz_init(&(*mat)[i][j]);
        }
    }
}
void
ibz_mat_4x5_finalize(ibz_mat_4x5_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            ibz_finalize(&(*mat)[i][j]);
        }
    }
}

void
ibz_mat_4x8_init(ibz_mat_4x8_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_init(&(*mat)[i][j]);
        }
    }
}
void
ibz_mat_4x8_finalize(ibz_mat_4x8_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_finalize(&(*mat)[i][j]);
        }
    }
}

void
ibz_mat_init(int rows, int cols, ibz_t mat[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            ibz_init(&mat[i][j]);
}

void
ibz_mat_finalize(int rows, int cols, ibz_t mat[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            ibz_finalize(&mat[i][j]);
}

void
quat_lattice_init(quat_lattice_t *lat)
{
    ibz_mat_4x4_init(&(*lat).basis);
    ibz_init(&(*lat).denom);
    ibz_set(&(*lat).denom, 1);
}
void
quat_lattice_finalize(quat_lattice_t *lat)
{
    ibz_finalize(&(*lat).denom);
    ibz_mat_4x4_finalize(&(*lat).basis);
}

void
quat_order_init(quat_order_t *order)
{
    quat_lattice_init(order);
}
void
quat_order_finalize(quat_order_t *order)
{
    quat_lattice_finalize(order);
}

void
quat_left_ideal_init(quat_left_ideal_t *lideal)
{
    quat_lattice_init(&(*lideal).lattice);
    ibz_init(&(*lideal).norm);
    (*lideal).parent_order = NULL;
}
void
quat_left_ideal_finalize(quat_left_ideal_t *lideal)
{
    ibz_finalize(&(*lideal).norm);
    quat_lattice_finalize(&(*lideal).lattice);
}
