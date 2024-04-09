
#include <stdint.h>

// unsigned 128-bit type integer
typedef unsigned __int128 uint128_t;

#include <stdint.h>

/// Addition of (x + y + c) where x, y < 2**64, c < 2**8
/// Output r, cc where (x + y + c) = 2**64 * cc + r
inline uint8_t addcarry_u64(uint64_t *r, uint64_t x, uint64_t y, uint8_t c) {
  uint128_t z;
  z = (uint128_t)x + (uint128_t)y + (uint128_t)c;
  *r = (uint64_t)z;
  return (uint8_t)(z >> 64);
}

/// Computation of (x - y - c) where x, y < 2**64, c < 2**8
/// Output r, cc where (x - y - c) = 2**64 * cc + r
inline uint8_t subborrow_u64(uint64_t *r, uint64_t x, uint64_t y, uint8_t c) {
  uint128_t z;
  z = ((uint128_t)x - (uint128_t)y) - (uint128_t)c;
  *r = (uint64_t)z;
  return (uint8_t)(z >> 127);
}

/// Computation of x * y where x, y < 2**64
/// Output hi, low where x * y = 2**64 * hi + low
inline uint64_t umull(uint64_t *hi, uint64_t x, uint64_t y) {
  uint128_t z;
  z = (uint128_t)x * (uint128_t)y;
  *hi = (uint64_t)(z >> 64);
  return (uint64_t)z;
}

/// Computation of x * y + z where x, y, z < 2**64
/// Output hi, low where (x * y + z) = 2**64 * hi + low
inline uint64_t umull_add(uint64_t *hi, uint64_t x, uint64_t y, uint64_t z) {
  uint128_t t;
  t = (uint128_t)x * (uint128_t)y + (uint128_t)z;
  *hi = (uint64_t)(t >> 64);
  return (uint64_t)t;
}

/// Computation of x * y + z1 + z2 where x, y, zi < 2**64
/// Output hi, low where (x * y + z1 + z2) = 2**64 * hi + low
inline uint64_t umull_add2(uint64_t *hi, uint64_t x, uint64_t y, uint64_t z1,
                           uint64_t z2) {
  uint128_t t;
  t = (uint128_t)x * (uint128_t)y + (uint128_t)z1 + (uint128_t)z2;
  *hi = (uint64_t)(t >> 64);
  return (uint64_t)t;
}
