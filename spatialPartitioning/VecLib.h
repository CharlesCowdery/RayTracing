#pragma once

#include "XYZ.h"
#include "XorRandom.h"
#include <vector>
#include <immintrin.h>

#define PI 3.14159265358979323846264

#define LOOKUP_SIZE_TRIG 1024
#define LOOKUP_TRIG_MASK (LOOKUP_SIZE_TRIG-1)
#define LOOKUP_SIZE_BIASED_HEMI 1024*32
#define LOOKUP_BIASED_HEMI_MASK (LOOKUP_SIZE_BIASED_HEMI-1)

namespace VecLib {
    __m256 sgn_fast(__m256& x);
    void cross_avx(const m256_vec3& v1, const m256_vec3& v2, m256_vec3& output);
    void dot_avx(const m256_vec3& v1, const m256_vec3& v2, __m256& output);
    double poly_acos(float x);
    __m256 inline_dot_avx(const m256_vec3& v1, const m256_vec3& v2);
    void dot_mul_avx(const m256_vec3& v1, const m256_vec3& v2, const __m256& v3, __m256& output);
    void prep(XorGen& G);

    float Lsin(float input);

    float Lcos(float input);

    float Lacos(float input);

    XYZ gen_random_point_on_sphere(XorGen& G);

    XYZ generate_biased_random_hemi_v2(XorGen& G, float r1_min = 0, float r1_max = 1, float r2_min = 0, float r2_max = 1);

    XYZ lookup_biased_random_hemi(XorGen& G);

    XYZ generate_biased_random_hemi(XorGen& G);

    XYZ generate_unbiased_random_hemi(XorGen& G);

    XYZ biased_random_hemi(XorGen& G, float r1_min = 0, float r1_max = 1, float r2_min = 0, float r2_max = 1);

    XYZ unbiased_random_hemi(XorGen& G);

    XYZ lookup_random_cone(XorGen& G, float spread);

    XYZ y_random_cone(XorGen& G, float spread);

    XYZ aligned_random(XorGen& G, float spread, const Quat& r);

    float surface_area(const XYZ& max, const XYZ& min);

    float volume(const XYZ& max, const XYZ& min);

    bool between(const XYZ& p1, const XYZ& p2, const XYZ test);

    bool betweenX(const XYZ& p1, const XYZ& p2, const XYZ test);

    bool betweenY(const XYZ& p1, const XYZ& p2, const XYZ test);

    bool betweenZ(const XYZ& p1, const XYZ& p2, const XYZ test);

    bool volumeClip(const XYZ& max_1, const XYZ& min_1, const XYZ& max_2, const XYZ& min_2);

    bool volumeContains(const XYZ& max, const XYZ& min, const XYZ& point);
}