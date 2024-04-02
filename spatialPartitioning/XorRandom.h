
#pragma once

#include <math.h>
#include <random>
#include <stdint.h>

/*Xorshiro256+ pseudorandom start. https://prng.di.unimi.it/*/

class XorGen {
public:
    void seed(uint64_t n1, uint64_t n2, uint64_t n3, uint64_t n4);
    uint64_t xe_next(void);
    void xe_jump(void);
    void xe_long_jump(void);
    float xe_frand();
    float rand_frand();
    float fRand(float fMin, float fMax);
private:
    uint64_t xe_seed[4] = { 0,0,0,0 };
    inline uint64_t rotl(const uint64_t x, int k);
};

