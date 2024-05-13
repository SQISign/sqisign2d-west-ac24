#include <stdint.h>

void fpadd5248_asm(uint64_t*, const uint64_t*, const uint64_t*);
void fpsub5248_asm(uint64_t*, const uint64_t*, const uint64_t*);
void mul5248_asm(uint64_t*, const uint64_t*, const uint64_t*);
void rdc5248_asm(uint64_t*, const uint64_t*);

const uint64_t p[4] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x4ffffffffffffff };
const uint64_t R[4] = { 0x0000000000000033, 0, 0, 0x0100000000000000 };
const uint64_t R2[4] = { 0x3333333333333d70, 0x3333333333333333, 0x3333333333333333, 0x0333333333333333 };

void fp_add(uint64_t* out, const uint64_t* a, const uint64_t* b) {
    fpadd5248_asm(out, a, b);
}

void fp_sub(uint64_t* out, const uint64_t* a, const uint64_t* b) {
    fpsub5248_asm(out, a, b);
}

void fp_mul(uint64_t* out, const uint64_t* a, const uint64_t* b) {
    uint64_t prod[8];
    mul5248_asm(prod, a, b);
    rdc5248_asm(out, prod);
}

void fp_sqr(uint64_t* out, const uint64_t* a) {
    fp_mul(out, a, a);
}

void fp_tomont(uint64_t* out, const uint64_t* a) {
    fp_mul(out, a, R2);
}

void fp_frommont(uint64_t* out, const uint64_t* a) {
    uint64_t one[4] = { 1 };

    fp_mul(out, one, a);
}

void fp_mont_setone(uint64_t* out) {
    out[0] = R[0];
    out[1] = R[1];
    out[2] = R[2];
    out[3] = R[3];
}
