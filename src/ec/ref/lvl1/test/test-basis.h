#ifndef TEST_BASIS_H
#define TEST_BASIS_H

#include "fp2.h"
// Full-torsion basis for A=0 (excluding 2^f and huge prime factors)
const uint64_t xPA_data[2][4] = {{0x344d1323cfceb9ff, 0xc7e92573b75b5fe4, 0x8b8a89fe49e2b24e, 0x45dc1b23143ce7f}, {0x770c0337ce832b58, 0x652d84ed2575bd1d, 0x6e96808c6fa2427c, 0x3b2d5015963565b}};
const uint64_t xQA_data[2][4] = {{0x7a3e24ac8aba9b98, 0xa68fc229bc4bdefb, 0x1dc9b16bf44ef377, 0x39d4144914f837}, {0x42e93a79ab10fe61, 0x3c61e50259425854, 0x9857e70cad607794, 0x34ada95a9693a5}};
const uint64_t xPQA_data[2][4] = {{0x65c2f0ce93da7630, 0xb5ba908cdf3d4e5b, 0x52db7c9e59caa1ec, 0x4edd26ed929ed4}, {0x9dc593aa424d6a8a, 0x9a462d8d7d2ed069, 0xaa95bdc036c57c48, 0x1ab8b1cc6801d80}};

const uint64_t xPB_data[2][4] = {{0x20968ebaa4b5d92, 0x2304a67673a00c59, 0x8f8f0f1968b59de4, 0x216843883c4f302}, {0xb52e25acb80d569d, 0xf02f8bdfa1488a61, 0x98b90130587b7a2f, 0x20ca99496add0e5}};
const uint64_t xQB_data[2][4] = {{0x7fa0440b5bbca752, 0x733fc6be1fd0b39d, 0x44b3fc5c987d415a, 0x3de3aeb3ad647a2}, {0x81c73b731a77165a, 0x15f26e4cb4664b1c, 0x637f10bcd2bb6ae9, 0x3d97517c12bc38a}};
const uint64_t xPQB_data[2][4] = {{0x5a166bb0102b8cc3, 0x32bce1f33a42a563, 0xc077e696f0338828, 0x30114343fb6641f}, {0xa9d8caecde9cf0f3, 0x33073a0782cbd891, 0x1c28fbc7055a67db, 0x5b614522205725}};

// 2^f-torsion basis for A=0
const uint64_t xP2_data[2][4] = {{0x50bef4c64b32a4f4, 0x5d7f9ce4d6d125ff, 0xc3c8fa96a800b8bb, 0x580cacb59b98a8}, {0xd03b7cf11fd096c, 0xb24f1f2f1126e474, 0xcb7f5ea7a9315aa, 0x456d2450ba2cc47}};
const uint64_t xQ2_data[2][4] = {{0x4081eacd4a3a0ac9, 0xb691c8b60a7200e2, 0x9ba3d980202f7b53, 0xe83693b1566335}, {0xf546abc4a7d8163a, 0xe316fd34b7634357, 0xa6c0cc7067250ffb, 0xac261d12dae56e}};
const uint64_t xPQ2_data[2][4] = {{0xecb666794189b874, 0x2eff89910481a488, 0xb555b3d755c55cdb, 0x1caa20697adf0f1}, {0x2216989fdfeeba06, 0x57f9236f4e63870f, 0xa53a1f716880e766, 0xa229b9dbefcd92}};

// 3^g-torsion basis for A=0
const uint64_t xP3_data[2][4] = {{0xca848d4a0970d02d, 0xa3b8dedce55d85b5, 0xc18eed9410b2d43f, 0x290d8feb2c74fb9}, {0x584284be6a861930, 0x834dd986a7cdc187, 0x824dd9e3f38a1b90, 0x492921816f3a637}};
const uint64_t xQ3_data[2][4] = {{0xe0b365dc1e0b76e0, 0x9ed3ba0df43113cf, 0x6c3fa356d2130d69, 0x2913172347ea846}, {0x4736d8e19acde7f1, 0x30765bb93466e02a, 0xdd5fdf324d183383, 0x26746db524d830b}};
const uint64_t xPQ3_data[2][4] = {{0xf557d152e0ee767, 0x931172eea1df3860, 0xc257bcb60a6b8344, 0x415b9a337f5367}, {0x2270de5a4594a23f, 0x785639d03fd46509, 0xf8dcaa9626e37909, 0x18febb33096ffec}};

#endif
