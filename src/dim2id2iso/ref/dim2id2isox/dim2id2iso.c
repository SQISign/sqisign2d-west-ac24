#include <quaternion.h>
#include <ec.h>
#include <endomorphism_action.h>
#include <dim2id2iso.h>
#include <inttypes.h>
#include <locale.h> 
#include <bench.h>
#include <curve_extras.h>
#include <tools.h>


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
int find_uv(ibz_t *u,ibz_t *v,ibz_t *au, ibz_t *bu,quat_alg_elem_t *beta1,quat_alg_elem_t *beta2,ibz_t *d1,ibz_t *d2, const ibz_t *target,int number_sum_square,const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {


    // variable declaration
    ibz_vec_4_t vec;
    ibz_t n;
    
    ibz_t remain,adjusted_norm;
    ibz_mat_4x4_t gram, reduced;

    ibz_init(&n);

    int prime_list_length; ibz_t prod_bad_primes; // TODO make this a precomputation
    short prime_list[12] = {2,5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97};
    prime_list_length = 12;  
    ibz_init(&prod_bad_primes);ibz_copy(&prod_bad_primes,&ibz_const_one); 
    ibz_set(&n,140227657289781369);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);
    ibz_set(&n,8695006970070847579);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);
    ibz_set(&n,4359375434796427649);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);
    ibz_set(&n,221191130330393351);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);
    ibz_set(&n,1516192381681334191);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);
    ibz_set(&n,5474546011261709671);
    ibz_mul(&prod_bad_primes,&prod_bad_primes,&n);

    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&reduced);
    
    ibz_vec_4_init(&vec);

    ibz_init(&remain);
    ibz_init(&adjusted_norm);

    // multiplying the norm by the denominator squared 
    ibz_set(&adjusted_norm,1);
    ibz_mul(&adjusted_norm,&adjusted_norm,&lideal->lattice.denom);
    ibz_mul(&adjusted_norm,&adjusted_norm,&lideal->lattice.denom);

    // computing a reduced basis
    clock_t t = tic();
    quat_lideal_reduce_basis(&reduced,&gram,lideal,Bpoo);
    TOC(t,"\nbasis reduction");

    // adjusting the sign if needed 
    // TODO : sometimes this fails
    if (ibz_cmp(&reduced[0][0],&reduced[1][1])!=0) {
        for (int i=0;i<4;i++) {
            ibz_neg(&reduced[1][i],&reduced[1][i]);
            ibz_neg(&gram[i][1],&gram[i][1]);
            ibz_neg(&gram[1][i],&gram[1][i]);
        }
        assert(ibz_cmp(&reduced[0][0],&reduced[1][1])==0);
    }
    if (ibz_cmp(&reduced[2][0],&reduced[3][1])!=0) {
        for (int i=0;i<4;i++) {
            ibz_neg(&reduced[3][i],&reduced[3][i]);
            ibz_neg(&gram[i][3],&gram[i][3]);
            ibz_neg(&gram[3][i],&gram[3][i]);
        }
        assert(ibz_cmp(&reduced[2][0],&reduced[3][1])==0);
    }

    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            ibz_div(&gram[i][j],&remain,&gram[i][j],&lideal->norm);
            assert(ibz_is_zero(&remain));
        }
    }

     // printing the size of the elements
    printf("p : %d n : %d adj_n : %d \n",ibz_bitsize(&QUATALG_PINFTY.p),ibz_bitsize(&lideal->norm),ibz_bitsize(&adjusted_norm));
    printf("n1 : %d n2 : %d n3 :  %d n4 : %d \n",ibz_bitsize(&gram[0][0])-ibz_bitsize(&adjusted_norm),ibz_bitsize(&gram[1][1])-ibz_bitsize(&adjusted_norm),ibz_bitsize(&gram[2][2])-ibz_bitsize(&adjusted_norm),ibz_bitsize(&gram[3][3])-ibz_bitsize(&adjusted_norm));





    // we start by enumerating a set of small vectors 
    // global parameters for the enumerate
    int m = 2+number_sum_square;
    int m4 = (2*m+1)*(2*m+1)*(2*m+1)*(2*m+1)-1;
    int m3 = (2*m+1)*(2*m+1)*(2*m+1);
    int m2 = (2*m+1)*(2*m+1);
    int m1 = 2*m+1;
    m4 = m4/4;
    ibz_vec_4_t small_vecs[m4];
    ibz_t small_norms[m4];

    t = tic();
    int index = 0; 
    for (int i1=-m;i1<m;i1++) {
        for (int i2=-m;i2<m+1;i2++) {
            for (int i3=-m;i3<m+1;i3++) {
                for (int i4=-m;i4<m+1;i4++) {
                    if (!(i1>0)&&!(i1==0 && i2>0)&&!(i1==0 && i2==0 && i3>0)&&!(i1==0 && i2==0 & i3==0 && i4>=0)
                    // &&!(i1>i2 && i2>-i1 && i3>i4 && i4>-i3)&&!(i1>-i2 && i2>i1 && i3>-i4 && i4>i3)
                    &&!((m+i4)+m1*(m+i3)+m2*(m+i2)+m3*(m+i1)>(m-i3)+m1*(m+i4)+m2*(m-i1)+m3*(m+i2))
                    &&!((m+i4)+m1*(m+i3)+m2*(m+i2)+m3*(m+i1)>(m+i3)+m1*(m-i4)+m2*(m+i1)+m3*(m-i2)) 
                    &&!(i1%2==0 && i2%2==0 && i3%2==0 && i4%2==0)
                    &&!(i1%3==0 && i2%3==0 && i3%3==0 && i4%3==0)
                    ) { 

                        ibz_set(&vec[0],i1);
                        ibz_set(&vec[1],i2);
                        ibz_set(&vec[2],i3);
                        ibz_set(&vec[3],i4);
                        quat_qf_eval(&n,&gram,&vec);
                        ibz_div(&n,&remain,&n,&adjusted_norm);
                        // if (ibz_mod_ui(&n,2)==1) {
                        if (1) {
                            
                            ibz_init(&small_norms[index]);
                            ibz_vec_4_init(&small_vecs[index]);
                            ibz_set(&small_vecs[index][0],i1);
                            ibz_set(&small_vecs[index][1],i2);
                            ibz_set(&small_vecs[index][2],i3);
                            ibz_set(&small_vecs[index][3],i4);
                            ibz_copy(&small_norms[index],&n);
                            index++;
    
                        } 
                    }
                    
                }
            }

        }
    }
    index--;
    printf("number of elements = %d / %d \n",index,m4);
    TOC(t,"\nenum");
    // printf("size = %d/%d \n",index,m4);
    
    
    // sorting the list
    quicksort(small_norms,small_vecs,0,index);
    TOC(t,"sorting + enum");

   
    int found=0;
    int count=0;
    // starting to try solutions

    // int i1=0;
    // int i2=0;
    // TODO try to go through d1,d2 by increasing size of the products d1*d2
    // TODO allow for even d1, d2 
    printf("\n");
    t= tic();
    clock_t tot_spec =0;
    clock_t tot_cor = 0;
    clock_t tot_ops = 0;
    clock_t tot_inv = 0;
    clock_t t_loc,t_loc2;

    // precomputing a list of small_norms[i]/target
    ibz_t quotients[index];
    for (int i=0;i<index;i++) {
        ibz_init(&quotients[i]);
    } 

    ibz_copy(&n,target);

    t_loc = tic();    
    for (int i=0;i<index;i++) {
        ibz_div(&quotients[i],&remain,&n,&small_norms[i]);
    }
    TOC(t,"list of quotients"); 
    tot_spec= tot_spec + dclock(t_loc);

    int cmp;

    index = 100*(2*number_sum_square+1);

    for (int i1=0;i1<index;i1++) {
        t_loc = tic();
        ibz_mod(&adjusted_norm,&n,&small_norms[i1]);
        tot_ops = tot_ops + dclock(t_loc);
        for (int i2=i1;i2<index;i2++) {        

            

            t_loc = tic();
            if (ibz_is_even(&small_norms[i1]) && ibz_is_even(&small_norms[i2])) {
                break;
            }
            // u = target / d1 mod d2
            t_loc2=tic();
            if (!ibz_invmod(&remain,&small_norms[i2],&small_norms[i1])) {
                continue;
            }
            tot_inv = tot_inv + dclock(t_loc2);
            ibz_mul(v,&remain,&n); 
            ibz_mod(v,v,&small_norms[i1]);
            t_loc2 =tic();
            // ibz_mul(&remain,v,&small_norms[i2]);
            // cmp=ibz_cmp(&remain,target);
            cmp = ibz_cmp(v,&quotients[i2]);            
            tot_spec = tot_spec + dclock(t_loc2); 
            tot_ops = tot_ops + dclock(t_loc);
            while (!found && cmp<0) {

                count++;
                t_loc = tic();
                if (number_sum_square>0) {
                    found = ibz_cornacchia_extended(au,bu,v,prime_list,prime_list_length,KLPT_primality_num_iter,&prod_bad_primes);
                    tot_cor = tot_cor + dclock(t_loc);
                }
                else {found = 1;}
                
                if (found) {
                    t_loc=tic();
                    ibz_mul(&remain,v,&small_norms[i2]);
                    ibz_sub(u,&n,&remain);
                    ibz_div(u,&remain,u,&small_norms[i1]);
                    assert(ibz_is_zero(&remain));
                    tot_ops = tot_ops+dclock(t_loc);
                    if (number_sum_square==2) {
                        found = ibz_cornacchia_extended(&remain,&n,v,prime_list,prime_list_length,KLPT_primality_num_iter,&prod_bad_primes);
                    }
                    
                }
                if (!found) {
                    t_loc = tic();
                    ibz_add(v,v,&small_norms[i1]);
                    t_loc2= tic();
                    cmp = ibz_cmp(v,&quotients[i2]);  
                    tot_ops = tot_ops +dclock(t_loc);
                    tot_spec = tot_spec + dclock(t_loc2);
                }
            }
            if (found) {
                printf("Found after %d attempts and %d deg one hit out of %d \n",i1*index + i2,count,index*(index-1)/2);
                break;
            }

        }
        if (found) {
            break;
        }
    }
    clock_to_time(tot_cor,"cornacchias :");
    clock_to_time(tot_ops,"all other ops :");
    clock_to_time(tot_inv,"among which modular inversions :");
    clock_to_time(tot_spec,"among which cmp");
    TOC(t,"total time searching for u,v");
    clock_to_time(dclock(t)-(tot_ops+tot_cor),"difference between total and recorded ");
    // printf(" difference between total clock cycles and recorded one %ld \n",dclock(t)-(tot_ops+tot_cor));



    // var finalize
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);

    ibz_finalize(&n);
    ibz_vec_4_finalize(&vec);

    ibz_finalize(&remain);
    ibz_finalize(&adjusted_norm);
    for (int i=0;i<index;i++) {
        ibz_finalize(&quotients[i]);
        ibz_finalize(&small_norms[i]);
        ibz_vec_4_finalize(&small_vecs[i]);
    }

    
    return found;
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
int dim2id2iso_ideal_to_isogeny_clapotis(theta_chain_t *isog, quat_alg_elem_t *beta1, quat_alg_elem_t *beta2, ibz_t *u, ibz_t *v, ibz_t *au, ibz_t *bu, theta_chain_t *phiv,ibz_t *d1,ibz_t *d2, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    
    ibz_t target;
    ibz_init(&target);
    int exp =2*TORSION_PLUS_EVEN_POWER;
    ibz_pow(&target,&ibz_const_two,exp);

    printf("looking for 2^%d \n",exp);

    int number_sum_quare=0;
    clock_t t = tic();
    // first, we find u,v,d1,d2,beta1,beta2
    int found = find_uv(u,v,au,bu,beta1,beta2,d1,d2,&target,number_sum_quare,lideal,Bpoo);

    TOC(t,"\n \ntotal time to find u,v");
    ibz_finalize(&target);

    return found;

}
