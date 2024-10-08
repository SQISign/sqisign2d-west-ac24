#ifndef TEST_BASIS_H
#define TEST_BASIS_H

#include "fp2.h"
// Full-torsion basis for A=0 (excluding 2^f and huge prime factors)
const fp2_t xPA = {
    { 0x3c780e636a5869dc,
     0xb8a1d106332efe8e, 0x7dd946e490e6578e,
     0x71d1fadbea881f88, 0xb94912baba3999f0,
     0x85343be0a74ca9e1, 0x22ae01775a9f7fa4,
     0x001032ffab70a66e },
    { 0x15908a4b85221a67,
     0x342f82e6a1db4e1d, 0x3d7c806a0d47b041,
     0x693830fad798c598, 0xcfa244134a61827a,
     0x7f723d6f5d9628cf, 0x10da657833d4d027,
     0x000c48499df01216 }
};
const fp2_t xQA = {
    { 0x79a766df9c10c642,
     0x7677cb85097be8be, 0x2a21c7f9b84b9deb,
     0xb263e837f57210ce, 0x551d6636b7c7e061,
     0x78d332581bee10b2, 0xce30a9926772e06c,
     0x00150b5009b1d6ed },
    { 0xbb2f097dae470eb9,
     0x53940c6df1eb93a9, 0x7786a4bab87320c1,
     0x89d32acc1c91db18, 0x733ef7f139fb7f9b,
     0x7bc336ee25a3901b, 0xf7dfe8f5559eeeb1,
     0x00210555ab63e7f3 }
};
const fp2_t xPQA = {
    { 0x315ead6fadc8b0d6,
     0x7da37e8b7e94de95, 0xcc6a9e206f513651,
     0x84fa9fab584acf3d, 0x293b25689ac50519,
     0xe3222bd1c8154964, 0x8ad7f39d04a8274f,
     0x000898edca69c223 },
    { 0x3e6c3e1864851e7e,
     0x01807c724f75ad5e, 0xe9cd50eff4e66fb7,
     0x6c7c19a88fed9707, 0x3ab57d0499386a40,
     0x6b5fd53c6efdc0b5, 0x092fe030da27bc43,
     0x00076f2f409c5f8e }
};

const fp2_t xPB = {
    { 0x229e388475511856,
     0x2f6b17e9ec9258c0, 0x0cb28c568697f9f4,
     0xca039e28512c9f9b, 0xd52d823761b0daa2,
     0xa09c3800e22c5e3b, 0x2971022668c3b76a,
     0x0006e91c4415afd1 },
    { 0xbd5059b7406e1dcd,
     0x9da456ed8c11f1a3, 0x1fb30e9cf66f928e,
     0x867c348b2f488d26, 0x9d4b03d8aa4229bc,
     0x1c01ca1088d145a8, 0xc9d6a201d77644a1,
     0x000a0d45131bf5b0 }
};
const fp2_t xQB = {
    { 0x712f0e5d0e3b4dfa,
     0x52260082dda1a07e, 0x5a7513dcfd273829,
     0xc686f0976cbb5dcf, 0xf5fc3df004cc7efc,
     0x615d0c2da4f2fb9f, 0x796efbb3f65aede8,
     0x00028176c42e1d9f },
    { 0xb8779b5a7bd2436b,
     0x4067b7e09d0ca56c, 0xfdbaee6ff27ebe38,
     0x69310e98174025de, 0x71960a10fa15706e,
     0x08ffb4b3f6efafbf, 0xb7116ca162211ea3,
     0x00253c0f60765f1f }
};
const fp2_t xPQB = {
    { 0x0e90506c89b46e0c,
     0x24ec65d5deb4e5b9, 0x8477f7e141db8725,
     0xf76957ec1940dbd3, 0xc2857af32534e715,
     0x06820654c6bae5f4, 0x5ac928ef3c90c1f8,
     0x0024f724366faeed },
    { 0xf6d7d2fdb06b91c4,
     0xe603cf05ce3f7555, 0x8a0876277637415c,
     0xa1ef891f00155f8f, 0x159db3ac93d39d57,
     0x5a05683aeaa453ff, 0x180c38da2402f6fc,
     0x000b69d01dcb9107 }
};

// 2^f-torsion basis for A=0
const fp2_t xP2 = {
    { 0x5d453ee3e6de9bf6,
     0xb5e51a5e88d8bbf3, 0xc91ce6ef41eda957,
     0x4e0ba74e86fd3385, 0xeff87c1def35e01f,
     0xedcd6c20496988a5, 0x91a2c14abdb955fe,
     0x000be92a3f4de175 },
    { 0xa8a13d8e0022a825,
     0xb26bb70885d42bef, 0x2533c31e799596b4,
     0xc41d58b247fb5ac9, 0x8d45fa188fd5cb65,
     0x1b0593f6e4af948d, 0x0ede22e4fcbe17ca,
     0x0014f54c5d5e1308 }
};
const fp2_t xQ2 = {
    { 0x90414b2365f868cd,
     0x68af18688f73fe25, 0x46ca4c4b4ca19114,
     0xadae5e2564f79c98, 0xfe3e09af9d00eb08,
     0x6856810a298a57bf, 0x170d41ba9327205d,
     0x001d588b6744b4ea },
    { 0xfb94e978bcf29be5,
     0x136700c07b264bd6, 0x62a3c89d8466b8f9,
     0x9f990ca7d3084bd8, 0xaab6fb1040e242d0,
     0x9e9325c5a5c20740, 0xa9a6ee97f376e198,
     0x0003c8eee3581511 }
};
const fp2_t xPQ2 = {
    { 0x873d426c501eafe6,
     0xdeb1e87769484669, 0x57c38f42bd1fef4d,
     0x53ca12d14b2ded18, 0xb72ef4a808fc9d70,
     0x59d9a54b1844cca1, 0x6ca7ccb15b6a9e49,
     0x00132a12929654f7 },
    { 0xffc6b824b6603270,
     0xb4152cbd3b607298, 0xbe97764acdcb16ce,
     0x5205b1ec222c3be9, 0x0cf5ac18d1eb4984,
     0xf5233664fd72c328, 0x492e775887a3367c,
     0x001ce6bdfc847b45 }
};

// 3^g-torsion basis for A=0
const fp2_t xP3 = {
    { 0x807a6abcb56d1915,
     0x3ab8ff7df809ea8f, 0x2bd4f1eba48b23ac,
     0xeb32542370dde5ff, 0xe6c50551eaaf2329,
     0x545dceaf98f07f09, 0x90bfb0e10f3e5b48,
     0x000cc0084da1b367 },
    { 0xbd6f9c82cd4acc13,
     0x9b39d0711267d8a2, 0x0ff31ab9fd38bb36,
     0xccc169cd75c1a58b, 0xd943ad3571e304b4,
     0xfc3cda0859595d00, 0xabda66362732b019,
     0x00070c5abcf1f329 }
};
const fp2_t xQ3 = {
    { 0x2b46bbfa6e57a9db,
     0xa7a5881479d3aaff, 0x5c8106d57698b7cb,
     0xde0ccd3c436cd1ad, 0xed351e8fbc28fd8f,
     0xe18a9a18e4f5bf03, 0x9a98961a81073911,
     0x001ed93f47abe8f2 },
    { 0x5dc96ddee6e9a9eb,
     0x5e8905d15b918006, 0xe89cecdc3f9b48f1,
     0x9d1a98543001e35e, 0x0795c7b134dadeba,
     0x8050c48376f36d87, 0xe9f364f7c6fbee1f,
     0x00061cb05b384f81 }
};
const fp2_t xPQ3 = {
    { 0xd44970f662987227,
     0x4c8eda7256920e8d, 0x857f42e972e25a0e,
     0xc66a5b62daa3644d, 0x6ab4ded74a464c38,
     0x4157cc1048b85a3a, 0x9916ab1ee4e2305a,
     0x000c6943137ffba1 },
    { 0x0c5118f818e5279d,
     0xacb0c4a011613c7a, 0xb87b4a9cb16a7565,
     0xc997ccbe0159f318, 0x6fc50720bce6f45f,
     0xbd1916a5ca7789d7, 0x3f48f437fdeccc64,
     0x000674d925340bc4 }
};

#endif