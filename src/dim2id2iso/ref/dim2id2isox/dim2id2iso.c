#include <quaternion.h>
#include <ec.h>
#include <endomorphism_action.h>
#include <dim2id2iso.h>
#include <inttypes.h>
#include <locale.h> 
#include <bench.h>
#include <curve_extras.h>


void swap(ibz_t *a, ibz_t *b,ibz_vec_4_t *va,ibz_vec_4_t *vb) {
    ibz_t temp;
    ibz_vec_4_t vtemp;
    ibz_init(&temp);
    ibz_vec_4_init(&vtemp);
    ibz_copy(&temp,a);
    ibz_copy(&vtemp[0],&(*va)[0]);
    ibz_copy(&vtemp[1],&(*va)[1]);
    ibz_copy(&vtemp[2],&(*va)[2]);
    ibz_copy(&vtemp[3],&(*va)[3]);
    ibz_copy(a,b);
    ibz_copy(&(*va)[0],&(*vb)[0]);
    ibz_copy(&(*va)[1],&(*vb)[1]);
    ibz_copy(&(*va)[2],&(*vb)[2]);
    ibz_copy(&(*va)[3],&(*vb)[3]);
    ibz_copy(b,&temp);
    ibz_copy(&(*vb)[0],&vtemp[0]);
    ibz_copy(&(*vb)[1],&vtemp[1]);
    ibz_copy(&(*vb)[2],&vtemp[2]);
    ibz_copy(&(*vb)[3],&vtemp[3]);
    ibz_finalize(&temp);
    ibz_vec_4_finalize(&vtemp);
}


int partition(ibz_t arr[],ibz_vec_4_t varr[], int low, int high) {
    ibz_t pivot;
    ibz_init(&pivot);
    ibz_copy(&pivot,&arr[high]);

    int i = low - 1;

    for (int j = low; j <= high - 1; j++) {
        if (ibz_cmp(&arr[j], &pivot) < 0) {
            i++;
            swap(&arr[i], &arr[j],&varr[i],&varr[j]);
        }
    }

    swap(&arr[i + 1], &arr[high],&varr[i+1],&varr[high]);

    ibz_finalize(&pivot);
    return i + 1;
}


void quicksort(ibz_t arr[],ibz_vec_4_t varr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr,varr, low, high);

        quicksort(arr,varr, low, pi - 1);
        quicksort(arr,varr, pi + 1, high);
    }
}

/** 
 * @brief Find good equivalent ideals
 * 
 * @param u Output : integer
 * @param v Output : integer 
 * @param beta1 Output : quaternion element 
 * @param beta2 Output : quaternion element
 * @param d1 Output : integer 
 * @param d2 Output : integer
 * @param number_sum_square : int 
 * @param lideal : left quaternion ideal
*/
int find_uv(ibz_t *u,ibz_t *v,quat_alg_elem_t *beta1,quat_alg_elem_t *beta2,ibz_t *d1,ibz_t *d2,int number_sum_square,const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    

    int m = 4;
    int m4 = m*m*m*m-1;

    // variable declaration
    ibz_vec_4_t small_vecs[m4];
    ibz_t small_norms[m4];
    ibz_t remain,adjusted_norm;

    ibz_init(&remain);
    ibz_init(&adjusted_norm);

    // multiplying the norm by the denominator squared 
    ibz_mul(&adjusted_norm,&lideal->norm,&lideal->lattice.denom);
    ibz_mul(&adjusted_norm,&adjusted_norm,&lideal->lattice.denom);

    // computing a reduced basis
    ibz_mat_4x4_t gram, reduced;
    quat_lideal_reduce_basis(&reduced,&gram,lideal,Bpoo);

    // we start by enumerating a set of small vectors
    for (int i1=0;i1<m;i1++) {
        for (int i2=0;i2<m;i1++) {
            for (int i3=0;i3<m;i3++) {
                for (int i4=0;i4<m;i4++) {
                    
                    if (!(i1==0 && i2==0 && i3==0 && i4==0)) {
                        int index = i1+m*i2+m*m*i3+m*m*m*i4 -1;
                        ibz_init(&small_norms[index]);
                        ibz_vec_4_init(&small_vecs[index]);
                        ibz_set(&small_vecs[index][0],i1);
                        ibz_set(&small_vecs[index][1],i2);
                        ibz_set(&small_vecs[index][2],i3);
                        ibz_set(&small_vecs[index][3],i4);

                        quat_qf_eval(&small_norms[index],&gram,&small_vecs[index]);
                        // TODO : do this on the gram matrix to avoid the need to divide inside the loop? 
                        ibz_div(&small_norms[index],&remain,&small_norms[index],&adjusted_norm);

                        assert(ibz_is_zero(&remain));

                    }
                    
                }
            }

        }
    }

    // sorting the list 
    quicksort(small_norms,small_vecs,0,m4-1);
    
    
    return 0;
} 


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
int dim2id2iso_ideal_to_isogeny_clapotis(theta_chain_t *isog, quat_alg_elem_t *beta1, quat_alg_elem_t *beta2, ibz_t *u, ibz_t *v, ibz_t *au, ibz_t bu, theta_chain_t *phiv,ibz_t *d1,ibz_t *d2, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    
    int number_sum_quare=1;

    // first, we find u,v,d1,d2,beta1,beta2
    find_uv(u,v,beta1,beta2,d1,d2,number_sum_quare,lideal,Bpoo);

    return 0;

}
