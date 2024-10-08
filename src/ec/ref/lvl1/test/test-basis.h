#ifndef TEST_BASIS_H
#define TEST_BASIS_H

#include "fp2.h"
// Full-torsion basis for A=0 (excluding 2^f and huge prime factors)
const uint64_t xPA_data[2][4] = {
    { 0x7505815fb30f099e, 0x89e78dbb4294c8df, 0x7db9b4b1f7716d7b, 0x13fcd4c87af65308 },
    { 0x93533c1017088fd4, 0x6df9e398a1bb4cb1, 0xc928f082be2e2b4c, 0x17aa7e2906bef0af }
};
const uint64_t xQA_data[2][4] = {
    { 0xe96336b75eb5a505, 0x5640cecad0ad7b5a, 0x1394f0771bc58ac1, 0x18d92124656d68d9 },
    { 0xa54e8e24605754f0, 0xe52de9790bbe4bb9, 0x3bf9b7833f62e255, 0x277a07644ec4f0e2 }
};
const uint64_t xPQA_data[2][4] = {
    { 0xc8fcceb408e3444c, 0x9f8ca4d2c05c3287, 0x259e496f17c0f529, 0x0eb18a51c2a3dd1a },
    { 0x1014dbe2534b8310, 0x6b035ee3c371ea12, 0x8354ecb4c111db6d, 0x178259b78fe08093 }
};

const uint64_t xPB_data[2][4] = {
    { 0xbd0a2f0c9a5378ca, 0x74af17405042203d, 0x0ccdcb4b7f0b8c15, 0x314c70951a92d8bf },
    { 0xe889e6bc5f9842af, 0xefb0edbb5e266ab3, 0x7bfb9d05f1ba6962, 0x0a5f3f4fe6f16514 }
};
const uint64_t xQB_data[2][4] = {
    { 0x137e215438caaf3b, 0xc4403ee1b69f1382, 0x2b5783edcefa7246, 0x3015572698262f66 },
    { 0x8e88e4293f84536e, 0x8d6dbc277f85ff77, 0xb3f17b53b01da916, 0x08dd3f4976c5dad1 }
};
const uint64_t xPQB_data[2][4] = {
    { 0xf0c2701a7050d9b9, 0xc8fdb069c0234d3a, 0x9ec25780f2b101a8, 0x221a0565053e8ff4 },
    { 0xd8513bf6a05910ae, 0x47ff2422258dfb3a, 0xb98ccceae31ac407, 0x21bcc8e659aaa1b3 }
};

// 2^f-torsion basis for A=0
const uint64_t xP2_data[2][4] = {
    { 0xfc93bac7df77fd30, 0xa8d37e10783215bd, 0x4bd2ece4f148039b, 0x2bd5b83f5f8c09fb },
    { 0x444112970b59f12f, 0x557b8b9beb55c276, 0x633f97cd9464df6c, 0x00a1b21b593a2dfd }
};
const uint64_t xQ2_data[2][4] = {
    { 0x6b4289960273222c, 0xa290d8eb8e343a04, 0x0c0a333f80a0ed68, 0x31a58910e276aff0 },
    { 0xb7ca615ad7473865, 0xeb6f72f20794f050, 0x2941c3fe3203b94f, 0x32ad5cbe915e467b }
};
const uint64_t xPQ2_data[2][4] = {
    { 0xac9f90005e47b095, 0x47eafdafd5168836, 0xb88aac8334acdad0, 0x1a5cf52a20f665b4 },
    { 0x4baa70fb1f5fa99c, 0xffb7ddb12c87f1a3, 0xdd3a229d370a8484, 0x1e992ad0a14baf03 }
};

// 3^g-torsion basis for A=0
const uint64_t xP3_data[2][4] = {
    { 0x8cf496c2722f340d, 0x3e329c5a507ad39c, 0xa0c7caa3e4537e25, 0x1371d43cf97de48e },
    { 0xa4b94c97b8149e7d, 0xd290853fa14704c7, 0x158b854173c1b289, 0x04c6dcda7872c23f }
};
const uint64_t xQ3_data[2][4] = {
    { 0x0f6380fd4c963950, 0x101a22a245c4f563, 0x601d3e30b21a5f43, 0x0becd5f73b067949 },
    { 0xd364123c6806057e, 0x8ff24fca9e060260, 0x3b52df5bfb817901, 0x30950462489b838f }
};
const uint64_t xPQ3_data[2][4] = {
    { 0xe04cab7169e64a82, 0x56df573ea9295c19, 0x06cbb6af8e341990, 0x0f1046ca03017ca1 },
    { 0x2dac3457c35be728, 0x2f59af21113f25f9, 0xa0dc4f54eec2715d, 0x102ecf9a7ff2f2ff }
};

#endif
