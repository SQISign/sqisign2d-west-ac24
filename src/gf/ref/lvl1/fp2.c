#include "fp2.h"
#include <encoded_sizes.h>

/* Arithmetic modulo X^2 + 1 */

void fp2_from_w64(fp2_t * out, const uint64_t data[2][NWORDS_FIELD]){
    fp_from_w64(&(out->re), data[0]);
    fp_from_w64(&(out->im), data[1]);
}

void fp2_to_w64(uint64_t data[2][NWORDS_FIELD], const fp2_t * a){
    fp_to_w64(data[0], &(a->re));
    fp_to_w64(data[0], &(a->im));
}

// TODO test these!
void fp2_encode(void *dst, const fp2_t *a){
	uint8_t *buf = dst;
	fp_encode(buf, &(a->re));
	fp_encode(buf + FP_ENCODED_BYTES, &(a->im));
}

// TODO test these!
uint32_t fp2_decode(fp2_t *d, const void *src){
    const uint8_t *buf = src;
    uint32_t re, im;
    
    re = fp_decode(&(d->re), buf);
    im = fp_decode(&(d->im), buf + FP_ENCODED_BYTES);
    return re & im;
}

// TODO: we should be more careful hardcoding things 
// like this??
fp2_t fp2_non_residue()
{ // 2 + i is a quadratic non-residue for p1913
    fp2_t res;
    fp_set_small(&res.re, 2);
    fp_set_one(&res.im);
    return res;
}

void fp2_inv(fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);
    fp_inv(&t0);
    fp_mul(&(x->re), &(x->re), &t0);
    fp_mul(&(x->im), &(x->im), &t0);
    fp_neg(&(x->im), &(x->im));
}

void fp2_batched_inv(fp2_t *x, int len) {
    fp2_t t1[len],t2[len];
    fp2_t inverse;

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    t1[0] = x[0];
    for (int i=1;i<len;i++) {
        fp2_mul(&t1[i], &t1[i-1], &x[i]);
    }

    // inverse = 1/ (x0 * x1 * ... * xn)
    inverse = t1[len-1];
    fp2_inv(&inverse);
    t2[0] = inverse;

    // t2 = 1/ (x0 * x1 * ... * xn), 1/ (x0 * x1 * ... * x(n-1)) , ... , 1/xO
    for (int i=1;i<len;i++) {
        fp2_mul(&t2[i], &t2[i-1], &x[len-i]);
    }

    x[0] = t2[len-1];
    for (int i=1;i<len;i++){
        fp2_mul(&x[i], &t1[i-1], &t2[len-i-1]);
    }

}

bool fp2_is_square(const fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);

    return fp_is_square(&t0);
}

void fp2_frob(fp2_t* x, const fp2_t* y)
{
    x->re = y->re;
    fp_neg(&(x->im), &(y->im));
}

// NOTE: old, non-constant-time implementation. Could be optimized
void fp2_sqrt(fp2_t* x)
{
    fp_t sdelta, re, tmp1, tmp2, im;

    if (fp_is_zero(&(x->im))) {
        if (fp_is_square(&(x->re))) {
            fp_sqrt(&(x->re));
            return;
        } else {
            fp_neg(&(x->im), &(x->re));
            fp_sqrt(&(x->im));
            fp_set_zero(&(x->re));
            return;
        }
    }

    // sdelta = sqrt(re^2 + im^2)
    fp_sqr(&sdelta, &(x->re));
    fp_sqr(&tmp1, &(x->im));
    fp_add(&sdelta, &sdelta, &tmp1);
    fp_sqrt(&sdelta);

    fp_add(&re, &(x->re), &sdelta);
    fp_half(&re, &re);
    tmp2 = re;

    if (!fp_is_square(&tmp2)) {
        fp_sub(&re, &(x->re), &sdelta);
        fp_half(&re, &re);
    }

    fp_sqrt(&re);
    im = re;

    fp_inv(&im);
    fp_half(&im, &im);
    fp_mul(&(x->im), &im, &(x->im));   
    x->re = re; 
}


// Lexicographic comparison of two field elements. Returns +1 if x > y, -1 if x < y, 0 if x = y
int fp2_cmp(fp2_t* a, fp2_t* b){
    digit_t a_arr[NWORDS_FIELD];
    digit_t b_arr[NWORDS_FIELD];
    fp_to_w64(a_arr, &(a->re));
    fp_to_w64(b_arr, &(b->re));

    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a_arr[i] > b_arr[i])
            return 1;
        if(a_arr[i] < b_arr[i])
            return -1;
    }

    fp_to_w64(a_arr, &(a->im));
    fp_to_w64(b_arr, &(b->im));

    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a_arr[i] > b_arr[i])
            return 1;
        if(a_arr[i] < b_arr[i])
            return -1;
    }
    return 0;
}


// exponentiation 
// TODO could be improved
void fp2_pow(fp2_t *out,const fp2_t * x,const digit_t *exp,const int size) {

    fp2_t acc;
    digit_t exp_tmp[size];
    digit_t bit;

    memcpy((digit_t*)exp_tmp, (digit_t*)exp, size*RADIX/8);
    acc = *x;
    fp2_set_one(out);

    for (int i = 0; i < NWORDS_FIELD*RADIX; i++) {
        bit = exp_tmp[0] & 1;
        mp_shiftr(exp_tmp, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp2_mul(out, out, &acc);
        }
        fp2_sqr(&acc, &acc);
    }

}

void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set_one(&b);
    fp2_mul(&b, &b, &a);
    printf("%s0x", name);

    printf("%016llx", b.re.v0);
    printf("%016llx", b.re.v1);
    printf("%016llx", b.re.v2);
    printf("%016llx", b.re.v3);

    printf(" + i*0x");
    
    printf("%016llx", b.im.v0);
    printf("%016llx", b.im.v1);
    printf("%016llx", b.im.v2);
    printf("%016llx", b.im.v3);
    printf("\n");
}

/*
CAUTION: HAZARDOUS MATERIALS

The following tables are constants used for torsion basis generation. All values have already been converted into Montgomery
form, and so if the internal API for the field changes, these values may also need to change

NQR_TABLE:
    is a table of 20 values all of which are not squares in GF(p^2) with modulus x^2 + 1 and prime p = 5*2**248 - 1
Z_NQR_TABLE:
    is a table of 20 values for which the value is a square and (value - 1) is not a square in the field GF(p^2) 
    with modulus x^2 + 1 and prime p = 5*2**248 - 1 
*/

const uint64_t NQR_TABLE[20][2][NWORDS_FIELD] = {
    {{0xD79D529B2DD58744, 0xD3411D68DF181489, 0xC9E1ED24451A63AD, 0x025F87E62FE3C539}, {0xA5975731B068466B, 0x52046D01C46E4C61, 0xDCE9E28E58A80160, 0x0261AAF78311C5AE}},
    {{0x0E6A024EA4503262, 0x5190912777311F71, 0x8EB85BC3AC5F9702, 0x03DCA14E130C1AAB}, {0xD5CEEBB8EE1FAF65, 0x1C3427038C9F767C, 0x0924BDD76AFE7979, 0x03D48A9212A29F39}},
    {{0x5E811F355E7DBB4A, 0x3D22B2502135811A, 0x7D37635280CA238C, 0x00F0BD396DA777F4}, {0x3811BF126F111957, 0x65B6AFBAE6EF0EEF, 0x79B06A08C5779191, 0x02A68C0E60659080}},
    {{0x60874A127BCB51DD, 0xD8E89FA1633B9A8E, 0x9CB4C4C37991A9AF, 0x02689A703E4A0277}, {0xD1575D9EED835634, 0x316BD254FBEADF9F, 0x700D12636A701239, 0x01BB44C9C53CC72E}},
    {{0xDE52250706A3DE1D, 0x19832EFE0CE1B018, 0xD8F8EAF624BDBB3D, 0x0150FCCCD15EBEB2}, {0x582EF2F55BB7C536, 0xB6613515D8233D01, 0x0E4144CC346AB3D6, 0x03D0614DD22CF21E}},
    {{0x5ED5111A7C58D23E, 0x93FB2107C1B31496, 0x1A0F5E36622BA614, 0x02A769A468E38518}, {0x81AD3CA10011BFDC, 0xF661CF9A369612A4, 0x3F3620178903660D, 0x016027C14F273E5C}},
    {{0xAD1DFA971D40E696, 0x6A1A475CDB7A8245, 0x2F9A6E90C5B6D88B, 0x0086B68756C4588E}, {0x4296599EB1C5F9AF, 0x22EC1AA4613AF693, 0x9F2D2D048047F975, 0x02FDEA4728AEB163}},
    {{0x213469CECCF23114, 0x1060A84F970954B5, 0x2C78F4737E09AB18, 0x048A4F9446B25BA7}, {0x75C90430872CC6D0, 0x1EE78C48B0F07A94, 0x72FDF7203015AC24, 0x03E6EF1DFCEAF978}},
    {{0x857D645C406F4D32, 0xAD11D93ABB4B9FA6, 0xB501525C6069D8E7, 0x01838E7E5DCE8161}, {0x8E3EA655567C3D69, 0x0A6C22577CD93373, 0xBAD286A80E122524, 0x035A62C2E4457A42}},
    {{0x5705BE913508D784, 0xABACF1EC1ECF7350, 0xBA206371790005FB, 0x01C42951F68EB6AF}, {0x415E2613E59AECC1, 0x9BA371492C116484, 0xF435C3539765A167, 0x00D303BC18F4F95D}},
    {{0x40A6CAA7ED2CA710, 0x6F6FBD87203C1AE2, 0x96A57897B3D5D716, 0x00CAB4E28FA872A6}, {0x259C20FD3F755E95, 0xEFE89B19BE7D933C, 0xCBAC857CE54F2877, 0x03F9FFA018D683DE}},
    {{0xE300278584FA889D, 0x071D74F413388E21, 0xF6024900509632E3, 0x02542EA82ED08B4A}, {0xB222756871492878, 0xCA2ACC9D0EFF06AF, 0x088DDA4A4FD53661, 0x0059E27A9B18C6A0}},
    {{0x825F83E5978D47B8, 0x3AAE6F358C56B70A, 0x16A013FB97EA4D46, 0x033D96D383A01920}, {0x697E5CFA6D794BEB, 0xE30AB1D4206FAED2, 0x4856786BBF175880, 0x0044887541DF9430}},
    {{0x68A6BA2011B5087B, 0xFF877D82FE3548E9, 0xF5FC5E016D42FA12, 0x0467F3DFBFC50B5D}, {0x7CCA15CA1607AECA, 0xAB71852734A29B9C, 0x9919FF91A32BE585, 0x0408D6979F7C9812}},
    {{0xBA3F4CF1DDC8DD7B, 0x6B59B8B512BD7A8E, 0x64EDE6BDD2FEBF98, 0x0444B585F9D88ACB}, {0x5B50BED3D167CCF2, 0x0D383A4217837CF7, 0x84EC734848AFD580, 0x00EDB45A99B3DCD1}},
    {{0x6CEF4D6616B6DE5D, 0xBF794773A154C753, 0xF2CE06963929BEF2, 0x02EBE095969AF014}, {0xD587AD6CBB5AB4AA, 0x3417F99E28BA67D0, 0x74E6F5A14EB8E7BD, 0x03042CE5AB494A3E}},
    {{0xE071760697F2B7DB, 0x2B29D69C1669C627, 0xBED20542FC68BC98, 0x02B8DC22A3F74A1D}, {0x0C54F2B915024155, 0x9A28E4B30CF43B83, 0xFA26E4E7F43F39E4, 0x017CD5A53A16B363}},
    {{0xF66B090E03C74978, 0xE3379A190D081E81, 0xFF03C161FF538E7F, 0x034D5097027078A2}, {0x36CB99358808FD22, 0xF3E17C3EDEE70D1B, 0xB76A2643E0AE9ECF, 0x02FD556CA0BB489C}},
    {{0xC4E781B564FC3372, 0x98D83790DB4A10E7, 0xFB54DA62093A5120, 0x0253B6028B18E114}, {0xCA17FD29DD0A0CBB, 0x7D7538566E282F35, 0xF3D61A28034991F6, 0x04AA30374A64AE7B}},
    {{0xD573B3168F378A04, 0x44DF508E999D2249, 0x0BB9F32AA012ECD2, 0x00E9622238488D5D}, {0xD8AD7CBFD04EC933, 0xC4C5DECF64A8BFE0, 0x9BE611EF2DE8421A, 0x046B0407DADCAACA}},
};

const uint64_t Z_NQR_TABLE[20][2][NWORDS_FIELD] = {
    {{0xE163104B8C955FE5, 0x6E74A3F16600C0E2, 0xEB36AD6E3BF23D0A, 0x02C6C2ABD743B695}, {0x18BBEA0ED19572EF, 0x05956DDA26C2B1D7, 0xE96DEA307617B70E, 0x010FDA4043D24C07}},
    {{0xD01C09994BE5FD04, 0x28F5B080E7EC165F, 0x238AAE66B4F62FDA, 0x03762A6DAAC5002E}, {0xC0A2406503787419, 0x36BEDB1502531F59, 0x45E49A9F88DC2E1A, 0x0351699BA159EF17}},
    {{0x067CCB7203FFE696, 0x84381E8913B39671, 0xA94117CD96B1BC4B, 0x02B4731A7B31C7AA}, {0xA397B2917DF0112D, 0x8A0655B4731C6401, 0x7FFD6EF665B21913, 0x00C071D8C6F51BCD}},
    {{0x9D011911A0B9B553, 0x6289F5F0B9206DB4, 0x9527ACAB82590D8B, 0x01A3AC9A3F7CB985}, {0x0005938EA1855B31, 0x9A6DB02B018E9108, 0x4A8F5595C1B11D16, 0x01383DB6C01FCA79}},
    {{0x157DE916212CC8E8, 0x45979967D442DFA7, 0x32262D1D488F4368, 0x03595745F4E84B93}, {0x0B170146BD07303D, 0xE1C89F94CBC10791, 0x1DD359F9915B1734, 0x032CFECF470A8A72}},
    {{0x5B9A8CF60874D474, 0xC4DC6B630995A4D6, 0x90726993E4011012, 0x02EFC608252127C8}, {0xDC1ED0C2BFE5DD8C, 0xCBEEDED58382D793, 0x7DE995C91504EF69, 0x0215D25C489C4F25}},
    {{0xB1E2933AD4544166, 0xD5E37302990A35DB, 0x1BF14C96D3F8ED18, 0x0133F37BD9DC01C6}, {0x957110FC3AC92431, 0xB4EC21066371DC7C, 0xEC8C6F5DBA4A7B83, 0x04BF69D5F1C22D60}},
    {{0x4023C4DF3F3AD785, 0x5C278E339D6EEB5D, 0x3D702554EEB1A891, 0x0273EC769B63114A}, {0xF1D11E4F0DF3BAF6, 0xE42F7D76D3576070, 0x767903EE8DB8B828, 0x02B32C48ACF03878}},
    {{0x4512D062F9EEE551, 0x1157BFA169D37AE3, 0x0C5E15D2355770EA, 0x00DD2C83CA7B2DB0}, {0x4CB87A91ECE7E9E7, 0x3EBE403E5957A284, 0xEAE2308D11BC0E09, 0x02FD9B95AA57BA61}},
    {{0x56234818F4DB7F1B, 0x0E640687DC3E82DE, 0xD60AEAF1FF989E22, 0x00612701102ED552}, {0x575A4F89C81C1EF5, 0xF197A92A6AF0E61D, 0x82CF32EB6C4CC434, 0x007CCF08DCD0E5CF}},
    {{0xB4C361570C9F60E9, 0xB7E07EDFD33AD98F, 0x0498041C7B4B4994, 0x01E2CE69231D301F}, {0x1CC6B479823ACA2E, 0x86E4C727B16DBA61, 0xA636F605001FC8BA, 0x03AEFA11B6B21E78}},
    {{0x65166510D615FBF9, 0x7D27F65B20A53877, 0x659BC903F0F451F9, 0x044BA48CB284F333}, {0xC0D6D8D0F26A00AD, 0x9184894568F3B85C, 0x3AB432F56CE37413, 0x02599A3C52D80613}},
    {{0x1BDF776EAAF8441E, 0x20B2868293C3F2EF, 0xE5C3F1660C9B7C82, 0x00E7A4A4D5FD400B}, {0x52D730AB4B0DDD7F, 0x5F4504E48A8CF66A, 0x38710AF77D521E2A, 0x036BBE147FB618CA}},
    {{0x7AB105A891E776A1, 0x5936AB29CCD08111, 0xFC4E8E5DD517FF69, 0x02A1542DD3B5FC42}, {0x0D13FC3AA6C62882, 0xA8C3DB4819B39853, 0x171707D850E6A9BD, 0x01640D826B0DD1FB}},
    {{0x0DE09B2E611B1F45, 0xA1AD5338584F766D, 0xE80499DEAC81F6D0, 0x01D197A560EC57ED}, {0x3483CD479FA76962, 0x75ADA53E2A8F1698, 0x46C492BD364EA23B, 0x01B9E35195BDA6FD}},
    {{0xE1D7F65B515D4D39, 0x6F81DEA0D7AA9D0A, 0x1565CB4D168202BD, 0x0474DAA1E688C451}, {0x11AF3004AF79A7FA, 0x8B8F14E00EB39028, 0x6C2743E45573C62F, 0x04F7060F42999B8F}},
    {{0xFBD393ACD171415B, 0xD594EFAC07682F8F, 0x0701A5844DC7970B, 0x0247CCBF02E957F4}, {0xE75603B97C6CEFC5, 0xFC0A9EC0B1F04734, 0xD69F1E40BCCD1CBF, 0x048ECD95551C53E5}},
    {{0xFA29206905641163, 0x5B9D43693BB14BDA, 0x5517F1A307DCCB46, 0x027E5D2B05C2DC09}, {0x5013DE983A801BCC, 0x70CACC800221FB39, 0x6747684236EC7485, 0x0138C52BE83FF435}},
    {{0x21CB6AD42458DCE5, 0xDE864E48385426F3, 0xD9C9C15C3E383444, 0x02D37CFA988E99D0}, {0x723FF85F5BA7C940, 0xAF2E06B0F1DC2A07, 0x47BCBF212EC4671E, 0x03CA3F3FD0F69329}},
    {{0xD51CCAF091138758, 0xCE5017334A4D0D41, 0xB8D21552EBE65AE7, 0x04BCBA0A13C5CFD2}, {0x4D935E8ABC5BBFA7, 0x7A910B4AED1A48B5, 0x240499CF828389A7, 0x031F10ACAFC9FEEA}},
};
