#ifndef TEST_BASIS_H
#define TEST_BASIS_H

#include "fp2.h"
// Full-torsion basis for A=0 (excluding 2^f and huge prime factors)
const fp2_t xPA = {
    { 0x35b53c72e7494775,
     0x5791b499bc29710d, 0x2060f3aca68fa4ff,
     0x81150c19a14f523a, 0x08af6c81a906d44a,
     0x00cca2a93efb536e },
    { 0x14eaac356375af76,
     0x5655011e771be3b4, 0x6273ccee274d7754,
     0x440d6b5b4496c183, 0xa3d7f80e9f9111ba,
     0x0302e153bee01a18 }
};
const fp2_t xQA = {
    { 0x80c0767d1b7b5fd8,
     0x24e9039d430ca3b5, 0x26485254625dc85a,
     0x612eaebc345b64d1, 0x59669fbd946a4409,
     0x004c3a8564e16101 },
    { 0x0e1eac4e38449c54,
     0x752c042b4c6675cb, 0x88ec0e75c8e9ea0e,
     0xbf7c4cdbfc4483f0, 0xd594cb5474bbc264,
     0x02f5e2345a9b4654 }
};
const fp2_t xPQA = {
    { 0x1f5accaff9a7da90,
     0x91884964774d4cb2, 0x0e938e13dd088e63,
     0x453c9af09879a724, 0xb2bd09ec3740312b,
     0x0007a5837e23aaa1 },
    { 0x8e1ac4b319787bd4,
     0x7cb9fba402f67bfe, 0x370b2951f9ec29cf,
     0x7a020172566f9d17, 0x063e31753d703130,
     0x01551136265bade6 }
};

const fp2_t xPB = {
    { 0xb702a70a8ae132ad,
     0x56d8804c83a8e696, 0x5ac3e12f4df1792e,
     0x0a89da435664746e, 0xd8758765206844bd,
     0x01a92f6e9e0e9296 },
    { 0x8aaab711b76b0959,
     0x210e6695ca5e5fdd, 0x593be0d75909ca12,
     0xfbc074d8ebdeb927, 0xb61fcc328d3756bc,
     0x0198a5942855c8bf }
};
const fp2_t xQB = {
    { 0x2b6b82b950b61fda,
     0x0ef2dd717daed334, 0x99dee4db0b268ac9,
     0x3534eb384e1fcaf0, 0xbaf112845a4f2d81,
     0x037f1492d8d815a1 },
    { 0x97e80590f9a0556b,
     0x7d9b4b87a22a7792, 0xda4534fe75595b4b,
     0xbe1092a2733c03e1, 0xbf5b1bd147b0d630,
     0x0125721476e5267f }
};
const fp2_t xPQB = {
    { 0xb7d459a56d4aebec,
     0x6ac7f10ba20e1e71, 0x9a95a8928507f7ef,
     0xc4c5aff6b97f3dfe, 0x644beb3e86806b77,
     0x022319eb6eaf072a },
    { 0x8ad0f6b18934790e,
     0xdad82b7b38e166bf, 0xcb08f5a3ab53d9a9,
     0xd2ff39b401ba8aba, 0xbff9b5e40ed9e5ce,
     0x03c1773791f554c0 }
};

// 2^f-torsion basis for A=0
const fp2_t xP2 = {
    { 0x7a26fdb0e5844206,
     0x0752b2ba140f7dfd, 0x1728013f8f5fe257,
     0xd05f129975ed6bba, 0xe736dbce707ad5a8,
     0x01f861715896d0be },
    { 0xdac046927a0c5352,
     0x5a42474ac156ff18, 0xe887982ff4c5a9ea,
     0x3875be6432251f1c, 0xdfae47315af877ee,
     0x005627f085582ecc }
};
const fp2_t xQ2 = {
    { 0xc4f03ab3db57331b,
     0xf04261fc3b713778, 0xa99b82430c7e40d1,
     0x5fe52b1324c2a091, 0xfcaa2a7049d0f657,
     0x021f2caa09302141 },
    { 0x4a92a1d5ff9f6730,
     0x6dcd5f600f33783e, 0xdb8b4e2e5149b45e,
     0x993458635c01d0c0, 0x5f9bc7d3bb307f91,
     0x01fcc7eae4712b6a }
};
const fp2_t xPQ2 = {
    { 0x7f4ee9c86c4341a2,
     0x0c867f482063bdfc, 0xe46fb7b0fbd479c7,
     0xddaa716e091be9ad, 0x29239eadddf5dc59,
     0x0231c09c660f0a89 },
    { 0xde64fa344dd64237,
     0xa89aaaed3dd84555, 0xbb70924d8fb73f27,
     0x0869ec018b3366dc, 0x47a0356ce742bcbc,
     0x00547dbda6dc094d }
};

// 3^g-torsion basis for A==0
const fp2_t xP3 = {
    { 0x7c878d0ceaa821f0,
     0xf94db4cab7186625, 0x7cff6d5fb0ca7867,
     0x4e3f5bd19cbca9d6, 0x05ec8273d0042548,
     0x0233a79cf87040b3 },
    { 0x060e9f3dcab8192c,
     0xa94e86d063a46398, 0x0e5cc403bfb60867,
     0x3ea1277f98087283, 0xaff1fd95bb094917,
     0x025041b12719d3b8 }
};
const fp2_t xQ3 = {
    { 0xb25aaa192bd351b7,
     0xc5db1962aed7e543, 0x1f722ab174319947,
     0xd1c9bb4a0a5d8aa3, 0x351415ec64f88921,
     0x0288ae044d62c930 },
    { 0xb41ede1724f8e06a,
     0xfb10ce5a83c66629, 0x9846173e31a9d448,
     0x35c94966192f08db, 0x72f7252946af3f9c,
     0x02ea05c971e7b34c }
};
const fp2_t xPQ3 = {
    { 0x674703cc3134d90b,
     0x507e338e496b8f75, 0x0c8cb1f138346e4c,
     0x54cb7ad5ba580da7, 0x65750f0bcd0a9857,
     0x038b435f51669e87 },
    { 0xdcdc0116c67589a0,
     0x45ce94f4d345c827, 0x0f2cbfb3c53b73ea,
     0x03e7951bc98efbb8, 0x3335ad0991864858,
     0x01e151a64210f74f }
};

#endif