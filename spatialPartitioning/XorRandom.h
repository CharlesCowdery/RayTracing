
#pragma once

#include <math.h>
#include <random>
#include <stdint.h>
#include <immintrin.h>

/*Xorshiro256+ pseudorandom start. https://prng.di.unimi.it/*/

class XorGen {
public:
    void seed(uint64_t n1, uint64_t n2, uint64_t n3, uint64_t n4);
    uint64_t xe_next(void);
    void xe_jump(void);
    void regen_pool(void);
    void xe_long_jump(void);
    float xe_frand();
    float rand_frand();
    float fRand(float fMin, float fMax);
    float pool_rand(int index);
    __m256 pool;
private:
    uint64_t xe_seed[4] = { 0,0,0,0 };
    __m256i xe_seed_AVX[4];
    __m256i float_mask = _mm256_set1_epi32(0x007fffff); //clear sign and exponent bits
    __m256i float_exponent = _mm256_set1_epi32(0x3f800000); //sets exponent to 0111110 -> makes a random mantissa range from 0 to 1
    inline uint64_t rotl(const uint64_t x, int k);
};

