#include <quaternion.h>
#include <ec.h>
#include <endomorphism_action.h>
#include <dim2id2iso.h>
#include <inttypes.h>
#include <locale.h> 
#include <bench.h>
#include <curve_extras.h>
#include <id2iso.h>

//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
static void fp2_print(char *name, fp2_t const a){
    fp_t b1,b2;
    fp_frommont(b1,a.re);
    fp_frommont(b2,a.im);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b1[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b2[i]);
    printf("\n");
}

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}

static void theta_print(char *name, theta_point_t P) {
    fp2_t a;
    assert(!fp2_is_zero(&P.x));
    fp2_copy(&a,&P.x);
    fp2_inv(&a);
    fp2_mul(&a,&a,&P.y);
    fp2_print(name,a);
}



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
 * @brief Computes an arbitrary isogeny of fixed degree starting from E0
 * 
 * @param isog Output : a dim2 isogeny encoding an isogeny of degree u
 * @param lideal Output : an ideal of norm u
 * @param u : integer
 * @param extra_info : bit indicating if we want to use some the torsion with extra information
 * @returns a bit indicating if the computation succeeded  
 * 
 * F is an isogeny encoding an isogeny phi : E0 -> Eu of degree u
*/
int fixed_degree_isogeny(theta_chain_t *isog, quat_left_ideal_t *lideal, ibz_t *u, int extra_info) {

    // var declaration
    int found;
    ibz_t two_pow,tmp;
    quat_alg_elem_t theta;
    ec_curve_t E0 = CURVE_E0; 

    int length = TORSION_PLUS_EVEN_POWER-2;


    // var init 
    ibz_init(&two_pow);ibz_init(&tmp);
    quat_alg_elem_init(&theta);
    

    ibz_pow(&two_pow,&ibz_const_two,length);

    // check that u is not big 
    if (ibz_cmp(&two_pow,u)<0) {
        printf("too small \n");
        assert(0);
    }
    // check that u is odd 
    if (ibz_is_even(u)) {
        printf("even \n");
        assert(0);
    }

    // computing the endomorphism theta of norm u (2^(TORSION_PLUS_EVEN_POWER-2) -u)
    // TODO this is by default for now, we may want to change that later
    ibz_sub(&two_pow,&two_pow,u);
    ibz_mul(&two_pow,&two_pow,u);

    found = represent_integer_non_diag(&theta,&two_pow,&QUATALG_PINFTY);


    if (!found) {
        printf("not found \n");
        return 0;
    }
    quat_lideal_create_from_primitive(lideal,&theta,u,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);


    ec_basis_t B0_two;
    // copying the basis
    copy_point(&B0_two.P,&BASIS_EVEN.P);
    copy_point(&B0_two.Q,&BASIS_EVEN.Q);
    copy_point(&B0_two.PmQ,&BASIS_EVEN.PmQ);

    assert(test_point_order_twof(&B0_two.P,&E0,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&B0_two.Q,&E0,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&B0_two.PmQ,&E0,TORSION_PLUS_EVEN_POWER));
    

    // applying theta
    endomorphism_application_even_basis(&B0_two,&theta,TORSION_PLUS_EVEN_POWER);    

    assert(test_point_order_twof(&B0_two.P,&E0,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&B0_two.Q,&E0,TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&B0_two.PmQ,&E0,TORSION_PLUS_EVEN_POWER));

    // now we set-up the kernel 
    theta_couple_curve_t E01;
    theta_couple_point_t T1;
    theta_couple_point_t T2,T1m2;
    E01.E1=E0;
    E01.E2=E0;
    T1.P2=B0_two.P;
    T2.P2=B0_two.Q;
    T1m2.P2=B0_two.PmQ; 
      
    // multiplication by u 
    ec_biscalar_mul_ibz(&T1.P1,&E0,u,&ibz_const_zero,&BASIS_EVEN);
    ec_biscalar_mul_ibz(&T2.P1,&E0,&ibz_const_zero,u,&BASIS_EVEN);
    ibz_neg(&tmp,u);
    ibz_pow(&two_pow,&ibz_const_two,TORSION_PLUS_EVEN_POWER);
    ibz_mod(&tmp,&tmp,&two_pow);
    ec_biscalar_mul_ibz(&T1m2.P1,&E0,u,&tmp,&BASIS_EVEN);

   

    if (extra_info) {
        // computing the isogeny
        theta_chain_comput_strategy(isog,length,&E01,&T1,&T2,&T1m2,strategies[TORSION_PLUS_EVEN_POWER-length],extra_info);
    }
    else {
        
        // TODO change 
        int strategy[243] = {89, 55, 34, 31, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 10, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 2, 1, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1};
        double_couple_point_iter(&T1,2,&E01,&T1);
        double_couple_point_iter(&T2,2,&E01,&T2);
        double_couple_point_iter(&T1m2,2,&E01,&T1m2);
        // computing the isogeny
        theta_chain_comput_strategy(isog,length,&E01,&T1,&T2,&T1m2,strategies[TORSION_PLUS_EVEN_POWER-length+2],extra_info);
    }


    #ifndef NDEBUG
            if (extra_info) {
            ec_isom_t isom1,isom2;
            theta_chain_t isog2;
                int strategy2[243] = {89, 55, 34, 31, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 10, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 2, 1, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1};
            double_couple_point_iter(&T1,2,&E01,&T1);
            double_couple_point_iter(&T2,2,&E01,&T2);
            double_couple_point_iter(&T1m2,2,&E01,&T1m2);
            // computing the isogeny
            theta_chain_comput_strategy(&isog2,length,&E01,&T1,&T2,&T1m2,strategy2,0);

        

            fp2_t j,j2;
            ec_j_inv(&j,&isog->codomain.E1);
            ec_j_inv(&j2,&isog2.codomain.E1);
            ec_isom_t isom,isom_t;
            theta_couple_point_t im,im2,Help;
            ibz_t scal;
            ibz_init(&scal);
            digit_t scal_dig[NWORDS_ORDER] = {0};
            ibz_pow(&scal,&ibz_const_two,244);
            ibz_add(&scal,&ibz_const_one,&scal);
            ibz_to_digit_array(scal_dig,&scal);
            ec_mul(&Help.P1,&E0,scal_dig,&T1.P1);
            copy_point(&Help.P2,&isog->first_step.K1_4.P2);
            ibz_finalize(&scal);
            im.P1 = T1.P1;
            im2.P1 = T1.P1;
            // multiplying by the correct factor
            ec_set_zero(&im.P2);
            ec_set_zero(&im2.P2);
            theta_chain_eval(&im,isog,&im,&Help);
            theta_chain_eval(&im2,&isog2,&im2,&Help);

            assert(test_point_order_twof(&im.P1,&isog->codomain.E1,246));
            assert(test_point_order_twof(&im.P2,&isog->codomain.E2,246));
            assert(test_point_order_twof(&im2.P1,&isog2.codomain.E1,246));
            assert(test_point_order_twof(&im2.P2,&isog2.codomain.E2,246));
            if (fp2_is_equal(&j,&j2)) {
                ec_isomorphism(&isom,&isog->codomain.E1,&isog2.codomain.E1);
                ec_isomorphism(&isom2,&isog->codomain.E2,&isog2.codomain.E2);
                ec_iso_eval(&im.P1,&isom);
                ec_iso_eval(&im.P2,&isom2);
                assert(test_point_order_twof(&im.P1,&isog2.codomain.E1,246));
                assert(test_point_order_twof(&im.P2,&isog2.codomain.E2,246));
                assert(ec_is_equal(&im.P1,&im2.P1));
                assert(ec_is_equal(&im.P2,&im2.P2));
            }
            else {
                ec_isomorphism(&isom,&isog->codomain.E1,&isog2.codomain.E2);
                ec_isomorphism(&isom2,&isog->codomain.E2,&isog2.codomain.E1);
                ec_iso_eval(&im.P1,&isom);
                ec_iso_eval(&im.P2,&isom2);
                assert(test_point_order_twof(&im.P1,&isog2.codomain.E2,246));
                assert(test_point_order_twof(&im.P2,&isog2.codomain.E1,246));
                assert(ec_is_equal(&im.P1,&im2.P2));
                assert(ec_is_equal(&im.P2,&im2.P1));
            }


        }
    #endif
    

    // var finalize
    ibz_finalize(&two_pow);ibz_finalize(&tmp);
    quat_alg_elem_finalize(&theta);

    return 1;
    
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
int find_uv(ibz_t *u,ibz_t *v,ibz_vec_4_t *coeffs,quat_alg_elem_t *beta1,quat_alg_elem_t *beta2,ibz_t *d1,ibz_t *d2, const ibz_t *target,int number_sum_square,const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {


    // variable declaration
    ibz_vec_4_t vec;
    ibz_t n;
    ibz_t au,bu,av,bv;
    ibz_init(&au);
    ibz_init(&bu);
    ibz_init(&av);
    ibz_init(&bv);
    
    ibz_t remain,adjusted_norm;
    ibz_mat_4x4_t gram, reduced;

    ibz_init(&n);

    // TODO this could be much simpler (use ibz_set_str for a start)
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

    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            ibz_div(&gram[i][j],&remain,&gram[i][j],&lideal->norm);
            assert(ibz_is_zero(&remain));
        }
    }

    // reordering the basis if needed
    if (ibz_cmp(&gram[0][0],&gram[2][2])==0) {
        for (int i=0;i<4;i++) {
            ibz_swap(&reduced[i][1],&reduced[i][2]);
        }
        ibz_swap(&gram[0][2],&gram[0][1]);
        ibz_swap(&gram[2][0],&gram[1][0]);
        ibz_swap(&gram[3][2],&gram[3][1]);
        ibz_swap(&gram[2][3],&gram[1][3]);
        ibz_swap(&gram[2][2],&gram[1][1]);
    }
    else if (ibz_cmp(&gram[0][0],&gram[3][3])==0) {
        for (int i=0;i<4;i++) {
            ibz_swap(&reduced[i][1],&reduced[i][3]);
        }
        ibz_swap(&gram[0][3],&gram[0][1]);
        ibz_swap(&gram[3][0],&gram[1][0]);
        ibz_swap(&gram[2][3],&gram[2][1]);
        ibz_swap(&gram[3][2],&gram[1][2]);
        ibz_swap(&gram[3][3],&gram[1][1]);
    }
    assert(ibz_cmp(&gram[0][0],&gram[1][1])==0);
    // adjusting the sign if needed 
    if (ibz_cmp(&reduced[0][0],&reduced[1][1])!=0) {
        for (int i=0;i<4;i++) {
            ibz_neg(&reduced[i][1],&reduced[i][1]);
            ibz_neg(&gram[i][1],&gram[i][1]);
            ibz_neg(&gram[1][i],&gram[1][i]);
        }
        assert(ibz_cmp(&reduced[0][0],&reduced[1][1])==0);
    }
    if (ibz_cmp(&reduced[0][2],&reduced[1][3])!=0) {
        for (int i=0;i<4;i++) {
            ibz_neg(&reduced[i][3],&reduced[i][3]);
            ibz_neg(&gram[i][3],&gram[i][3]);
            ibz_neg(&gram[3][i],&gram[3][i]);
        }
        assert(ibz_cmp(&reduced[0][2],&reduced[1][3])==0);
    }

     // printing the size of the elements
    printf("p : %d n : %d adj_n : %d \n",ibz_bitsize(&QUATALG_PINFTY.p),ibz_bitsize(&lideal->norm),ibz_bitsize(&adjusted_norm));
    printf("n1 : %d n2 : %d n3 :  %d n4 : %d \n",ibz_bitsize(&gram[0][0]),ibz_bitsize(&gram[1][1])-ibz_bitsize(&adjusted_norm),ibz_bitsize(&gram[2][2])-ibz_bitsize(&adjusted_norm),ibz_bitsize(&gram[3][3])-ibz_bitsize(&adjusted_norm));


    // we start by enumerating a set of small vectors 
    // global parameters for the enumerate
    // TODO may be overshot
    int m = 2+2*number_sum_square;
    int m4 = (2*m+1)*(2*m+1)*(2*m+1)*(2*m+1)-1;
    int m3 = (2*m+1)*(2*m+1)*(2*m+1);
    int m2 = (2*m+1)*(2*m+1);
    int m1 = 2*m+1;
    m4 = m4/4;
    ibz_vec_4_t small_vecs[m4];
    ibz_t small_norms[m4];
    for (int i=0;i<m4;i++) {
        ibz_init(&small_norms[i]);
        ibz_vec_4_init(&small_vecs[i]);
    }

    t = tic();
    int index = 0; 
    for (int i1=-m;i1<m;i1++) {
        for (int i2=-m;i2<m+1;i2++) {
            for (int i3=-m;i3<m+1;i3++) {
                for (int i4=-m;i4<m+1;i4++) {
                    // we ensure that we don't record the same norm in the list
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
                        // TODECIDE : do we use only odd norms ?
                        if (ibz_mod_ui(&n,2)==1) {
                        // if (1) {
                            
                            
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
    printf("number of elements = %d / %d (value of m=%d) \n",index,m4,m);
    TOC(t,"\nenum");
    
    
    // sorting the list
    quicksort(small_norms,small_vecs,0,index);
    TOC(t,"sorting + enum");

   
    int found=0;
    int count=0;
    // starting to try solutions

    // TODO try to go through d1,d2 by increasing size of the products d1*d2
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
    int cnt_missed=0;

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
            cmp = ibz_cmp(v,&quotients[i2]);            
            tot_spec = tot_spec + dclock(t_loc2); 
            tot_ops = tot_ops + dclock(t_loc);
            while (!found && cmp<0) {

                // testing if we have
                ibz_t tt;
                ibz_init(&tt);
                ibz_mul(&tt,v,&ibz_const_two);
                if (ibz_cmp(&tt,&quotients[i2])<0) {
                    cnt_missed++;
                }
                ibz_finalize(&tt);

                count++;
                t_loc = tic();
                if (number_sum_square>0) {
                    found = ibz_cornacchia_extended(&av,&bv,v,prime_list,prime_list_length,KLPT_primality_num_iter,&prod_bad_primes);
                    tot_cor = tot_cor + dclock(t_loc);
                }
                else if (number_sum_square==0) {
                    found = 1;
                }
                
                if (found) {
                    t_loc=tic();
                    ibz_mul(&remain,v,&small_norms[i2]);
                    ibz_sub(u,&n,&remain);
                    ibz_div(u,&remain,u,&small_norms[i1]);
                    assert(ibz_is_zero(&remain));
                    tot_ops = tot_ops+dclock(t_loc);
                    if (number_sum_square==2) {
                        found = ibz_cornacchia_extended(&au,&bu,u,prime_list,prime_list_length,KLPT_primality_num_iter,&prod_bad_primes);
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

                // recording the solution that we found
                ibz_copy(&beta1->denom,&lideal->lattice.denom);
                ibz_copy(&beta2->denom,&lideal->lattice.denom);
                ibz_copy(d1,&small_norms[i1]);
                ibz_copy(d2,&small_norms[i2]);
                ibz_mat_4x4_eval(&beta1->coord,&reduced,&small_vecs[i1]);
                ibz_mat_4x4_eval(&beta2->coord,&reduced,&small_vecs[i2]);
                #ifndef NDEBUG
                    ibq_t norm;
                    ibq_init(&norm);
                    quat_alg_norm(&norm,beta1,&QUATALG_PINFTY);
                    ibq_to_ibz(&remain,&norm);
                    ibz_mul(&n,d1,&lideal->norm);
                    assert(ibz_cmp(&n,&remain)==0);
                    quat_alg_norm(&norm,beta2,&QUATALG_PINFTY);
                    ibq_to_ibz(&remain,&norm);
                    ibz_mul(&n,d2,&lideal->norm);
                    assert(ibz_cmp(&n,&remain)==0);
                    ibq_finalize(&norm);

                    // testing the values of coeffs
                    if (number_sum_square==2) {
                        ibz_mul(&n,&au,&au);
                        ibz_mul(&remain,&bu,&bu);
                        ibz_add(&n,&n,&remain);
                        assert(ibz_cmp(&n,u)==0);
                        ibz_mul(&n,&av,&av);
                        ibz_mul(&remain,&bv,&bv);
                        ibz_add(&n,&n,&remain);
                        assert(ibz_cmp(&n,v)==0);
                        
                    }

                #endif
                ibz_copy(&((*coeffs)[0]),&au);
                ibz_copy(&((*coeffs)[1]),&bu);
                ibz_copy(&((*coeffs)[2]),&av);
                ibz_copy(&((*coeffs)[3]),&bv);
                ibz_printf("%Zd %Zd %Zd %Zd \n",au,bu,av,bv);
                ibz_vec_4_print(coeffs);

                printf("Found after %d attempts and %d deg one hit out of %d \n",i1*index + i2,count,index*(index-1)/2);
                break;
            }

        }
        if (found) {
            break;
        }
    }
    if (!found) {
        printf("NOT Found after %d attempts and %d deg one hit out of %d \n",index*(index-1)/2,count,index*(index-1)/2);
    }
    printf("Number of missed deg 1 : hit %d \n ",cnt_missed);
    clock_to_time(tot_cor,"cornacchias :");
    clock_to_time(tot_ops,"all other ops :");
    clock_to_time(tot_inv,"among which modular inversions :");
    clock_to_time(tot_spec,"among which cmp");
    TOC(t,"total time searching for u,v");
    clock_to_time(dclock(t)-(tot_ops+tot_cor),"difference between total and recorded ");

    // var finalize
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);

    ibz_finalize(&n);
    ibz_vec_4_finalize(&vec);
    ibz_finalize(&au);
    ibz_finalize(&bu);
    ibz_finalize(&av);
    ibz_finalize(&bv);
    ibz_finalize(&remain);
    ibz_finalize(&adjusted_norm);
    for (int i=0;i<index;i++) {
        ibz_finalize(&quotients[i]);
    }
    for (int i=0;i<m4;i++) {
        
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
int dim2id2iso_ideal_to_isogeny_clapotis(theta_chain_t *isog, quat_alg_elem_t *beta1, quat_alg_elem_t *beta2, ibz_t *u, ibz_t *v, ibz_vec_4_t *coeffs, theta_chain_t *phiv,ibz_t *d1,ibz_t *d2, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    
    ibz_t target,tmp,two_pow;
    quat_alg_elem_t theta;
    quat_alg_elem_t quat_tmp;
    quat_alg_elem_t quat_gcd_remove;

    ibq_t norm;
    ibq_init(&norm);
    ibz_t test1,test2;
    ibz_init(&test1);
    ibz_init(&test2);


    ibz_init(&target);ibz_init(&tmp);ibz_init(&two_pow);
    // TODECIDE allow for a smaller exponent ?
    int exp = TORSION_PLUS_EVEN_POWER;
    ibz_pow(&target,&ibz_const_two,exp);
    quat_alg_elem_init(&theta);
    quat_alg_elem_init(&quat_tmp);
    quat_alg_elem_init(&quat_gcd_remove);

    int number_sum_square=2;
    clock_t t = tic();
    // first, we find u,v,d1,d2,beta1,beta2
    int found = find_uv(u,v,coeffs,beta1,beta2,d1,d2,&target,number_sum_square,lideal,Bpoo);
    if (!found) {
        return 0;
    }
    TOC(t,"\n \ntotal time to find u,v");
    ibz_printf("u = %Zd  v= %Zd d1=%Zd d2=%Zd \n",*u,*v,*d1,*d2);

    // TODO the following works only when d1,d2 are odd
    assert(ibz_get(d1)%2==1 && ibz_get(d2)%2==1); 
    // compute the valuation of the GCD of u,v 
    ibz_gcd(&tmp,u,v);
    assert(ibz_get(&tmp)!=0);
    int exp_gcd = two_adic_valuation(ibz_get(&tmp));
    exp = TORSION_PLUS_EVEN_POWER - exp_gcd;
    // removing the power of 2 from u and v
    ibz_div(u,&test1,u,&tmp);
    assert(ibz_cmp(&test1,&ibz_const_zero)==0);
    ibz_div(v,&test1,v,&tmp);
    assert(ibz_cmp(&test1,&ibz_const_zero)==0);

    // setting-up the element to remove the power of two from phiu, phiv
    // quat_gcd_remove = (i+1)^(exp_gcd)
    ibz_set(&quat_gcd_remove.denom,1);
    if (exp_gcd >0) {
        ibz_set(&quat_gcd_remove.coord[1],1);
        printf("exp sup %d \n",exp_gcd);
    }
    else {
        printf("exp 0 \n");
        ibz_set(&quat_gcd_remove.coord[1],0);
    }
    ibz_set(&quat_gcd_remove.coord[0],1);
    ibz_set(&quat_gcd_remove.coord[2],0);
    ibz_set(&quat_gcd_remove.coord[3],0);
    for (int i=0;i<exp_gcd-1;i++) {
        quat_alg_mul(&quat_gcd_remove,&quat_gcd_remove,&quat_gcd_remove,&QUATALG_PINFTY);
    }
    ibz_pow(&quat_gcd_remove.denom,&ibz_const_two,exp_gcd);

    #ifndef NDEBUG 
    
            quat_alg_norm(&norm,&quat_gcd_remove,&QUATALG_PINFTY);
            ibq_denom(&test1,&norm);
            ibz_printf("%Zd %Zd \n",test1,quat_gcd_remove.denom);
            assert(ibz_cmp(&test1,&quat_gcd_remove.denom)==0);

        #endif

    // now we compute the dimension 2 isogeny 
    // F : Eu x Ev -> E x E' 
    // where we have phi_u : Eu -> E0 and phi_v : Ev -> E0
    // if we have phi1 : E0 -> E of degree d1 
    // and phi2 : E0 -> E of degree d2 
    // we can define theta = phi2 o hat{phi1}  
    // and the kernel of F is given by 
    // ( [ud1](P), phiv o theta o hat{phiu} (P)),( [ud1](Q), phiv o theta o hat{phiu} (Q)) where P,Q is a basis of E0[2e]


    // now we set-up the kernel 
    ec_curve_t E0 = CURVE_E0;
    ec_basis_t bas;
    theta_couple_curve_t E01;
    theta_couple_point_t T1;
    theta_couple_point_t T2,T1m2;
    E01.E1=E0;
    E01.E2=E0;


    copy_point(&bas.P,&BASIS_EVEN.P);
    copy_point(&bas.Q,&BASIS_EVEN.Q);
    copy_point(&bas.PmQ,&BASIS_EVEN.PmQ);

    for (int i=0;i<TORSION_PLUS_EVEN_POWER-exp;i++) {
        ec_dbl(&bas.P,&E0,&bas.P);
        ec_dbl(&bas.Q,&E0,&bas.Q);
        ec_dbl(&bas.PmQ,&E0,&bas.PmQ);
    }


    // compute the first point of the basis by scalar multiplication by u * d1 
    ibz_mul(&target,u,d1);
    ec_biscalar_mul_ibz(&T1.P1,&E0,&target,&ibz_const_zero,&bas);
    ec_biscalar_mul_ibz(&T2.P1,&E0,&ibz_const_zero,&target,&bas);
    ibz_neg(&tmp,&target);
    ibz_pow(&two_pow,&ibz_const_two,exp);
    ibz_mod(&tmp,&tmp,&two_pow);
    ec_biscalar_mul_ibz(&T1m2.P1,&E0,&target,&tmp,&bas);

    if (number_sum_square==0) {

    }
    else if (number_sum_square==1) {

    }
    else {
        assert(number_sum_square==2);
        // in this case we have phiu = coeffs[0] + i coeffs[1] and phiv = coeffs[2] + i*coeffs[3]
        // and theta = beta2 \hat{beta1}/n 
        // we start by computing theta
        ibz_set(&theta.denom,1);
        quat_alg_conj(&theta,beta1);
        quat_alg_mul(&theta,beta2,&theta,&QUATALG_PINFTY);
        ibz_mul(&theta.denom,&theta.denom,&lideal->norm);

        #ifndef NDEBUG 
    
            quat_alg_norm(&norm,&theta,&QUATALG_PINFTY);
            ibq_to_ibz(&test1,&norm);
            ibz_mul(&test2,d1,d2);
            assert(ibz_cmp(&test1,&test2)==0);

        #endif

        // phiu 
        ibz_set(&quat_tmp.denom,1);
        ibz_copy(&quat_tmp.coord[0],&((*coeffs)[0]));
        ibz_copy(&quat_tmp.coord[1],&((*coeffs)[1]));
        ibz_set(&quat_tmp.coord[2],0);
        ibz_set(&quat_tmp.coord[3],0);
        quat_alg_mul(&quat_tmp,&quat_tmp,&quat_gcd_remove,&QUATALG_PINFTY);

        #ifndef NDEBUG 
    
            quat_alg_norm(&norm,&quat_tmp,&QUATALG_PINFTY);
            assert(ibq_to_ibz(&test1,&norm));
            assert(ibz_cmp(&test1,u)==0);

        #endif

        // theta <- phiu * theta
        quat_alg_mul(&theta,&quat_tmp,&theta,&QUATALG_PINFTY);

        // phiv    
        ibz_set(&quat_tmp.denom,1);
        ibz_copy(&quat_tmp.coord[0],&((*coeffs)[2]));
        ibz_copy(&quat_tmp.coord[1],&((*coeffs)[3]));
        ibz_set(&quat_tmp.coord[2],0);
        ibz_set(&quat_tmp.coord[3],0); 
        quat_alg_mul(&quat_tmp,&quat_tmp,&quat_gcd_remove,&QUATALG_PINFTY);

        // theta <- theta * phiv
        quat_alg_mul(&theta,&theta,&quat_tmp,&QUATALG_PINFTY);
        quat_alg_normalize(&theta);
        quat_alg_elem_print(&theta);
        assert(test_point_order_twof(&bas.P,&E0,exp));
        assert(test_point_order_twof(&bas.Q,&E0,exp));
        assert(test_point_order_twof(&bas.PmQ,&E0,exp));


        #ifndef NDEBUG 
            quat_alg_norm(&norm,&theta,&QUATALG_PINFTY);
            ibq_to_ibz(&test1,&norm);
            ibz_mul(&test2,&test2,u);
            ibz_mul(&test2,&test2,v);
            assert(ibz_cmp(&test1,&test2)==0);
            assert(ibz_get(&test1)%2==1);
            
        #endif

        // applying theta
        endomorphism_application_even_basis(&bas,&theta,exp); 

        // copying  
        copy_point(&T1.P2,&bas.P);
        copy_point(&T2.P2,&bas.Q);
        copy_point(&T1m2.P2,&bas.PmQ); 

        assert(test_point_order_twof(&bas.P,&E0,exp));
        assert(test_point_order_twof(&bas.Q,&E0,exp));
        assert(test_point_order_twof(&bas.PmQ,&E0,exp));

        
    }


    // TODO change 
    int strategy[245] = {89, 55, 34, 33, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 12, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1};


    assert(test_point_order_twof(&T1.P1,&E0,exp));
    assert(test_point_order_twof(&T1.P2,&E0,exp));
    assert(test_point_order_twof(&T2.P1,&E0,exp));
    assert(test_point_order_twof(&T2.P2,&E0,exp));
    assert(test_point_order_twof(&T1m2.P1,&E0,exp));
    assert(test_point_order_twof(&T1m2.P2,&E0,exp));

    theta_chain_comput_strategy(isog,exp,&E01,&T1,&T2,&T1m2,strategies[TORSION_PLUS_EVEN_POWER-exp+2],0);


    ibq_finalize(&norm);
    ibz_finalize(&test1);
    ibz_finalize(&test2);

    ibz_finalize(&target);
    ibz_finalize(&tmp);
    ibz_finalize(&two_pow);
    quat_alg_elem_finalize(&theta);
    quat_alg_elem_finalize(&quat_tmp);
    quat_alg_elem_init(&quat_gcd_remove);
    return found;

}
