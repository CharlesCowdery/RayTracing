#include "VecLib.h"

float arccos_lookup[LOOKUP_SIZE_TRIG];
float sin_lookup[LOOKUP_SIZE_TRIG];
float cos_lookup[LOOKUP_SIZE_TRIG];
XYZ biased_hemi_lookup[LOOKUP_SIZE_BIASED_HEMI];

const float coefs[5] = { 
0.0000231539316590538762175742441588523467,
-0.00138537043082318983893723662479142648,
0.0416635846931078386653947196040757567,
-0.49999905347076729097546897993796764,
0.999999953466670136306412430924463351 }; // gotten via genetic polynomial approximation

const __m128 coefs_ss[5] = {
    _mm_set_ss(coefs[0]),
    _mm_set_ss(coefs[1]),
    _mm_set_ss(coefs[2]),
    _mm_set_ss(coefs[3]),
    _mm_set_ss(coefs[4]) 
};

const float coefs_wc[7] = {
    1.7094536e-09f,
    -2.7113737e-07f,
    2.478466e-05f,
    -0.001388931f,
    0.04166709f,
    -0.5000007f,
    1.0000001f
};
const float coefs_ws[6] = {
    -1.926853e-08f,
    2.666575e-06f,
    -0.0001977116f,
    0.008330478f,
    -0.1666613f,
    0.9999967f
};

const __m256 coefs_AVX[5] = {
    _mm256_set1_ps(coefs[0]),
    _mm256_set1_ps(coefs[1]),
    _mm256_set1_ps(coefs[2]),
    _mm256_set1_ps(coefs[3]),
    _mm256_set1_ps(coefs[4])
};

const __m128 coefs_wc_ss_[7] = {
    _mm_set_ss(1.0f/ 6227020800.0f),
    _mm_set_ss(-1.0f / 39916800.0f),
    _mm_set_ss(1.0f / 362880.0f),
    _mm_set_ss(-1.0f / 5040.0f),
    _mm_set_ss(1.0f / 120.0f),
    _mm_set_ss(-1.0f / 6.0f),
    _mm_set_ss(1.0f / 1.0f)
};

const __m128 coefs_wc_ss[7] = {
    _mm_set_ss(-coefs_wc[0]),
    _mm_set_ss(-coefs_wc[1]),
    _mm_set_ss(-coefs_wc[2]),
    _mm_set_ss(-coefs_wc[3]),
    _mm_set_ss(-coefs_wc[4]),
    _mm_set_ss(-coefs_wc[5]),
    _mm_set_ss(-coefs_wc[6])
};
const __m128 coefs_ws_ss[6] = {
    _mm_set_ss(coefs_ws[0]),
    _mm_set_ss(coefs_ws[1]),
    _mm_set_ss(coefs_ws[2]),
    _mm_set_ss(coefs_ws[3]),
    _mm_set_ss(coefs_ws[4]),
    _mm_set_ss(coefs_ws[5])
};
const __m256 coefs_ws_AVX[6] = {
    _mm256_set1_ps(coefs_ws[0]),
    _mm256_set1_ps(coefs_ws[1]),
    _mm256_set1_ps(coefs_ws[2]),
    _mm256_set1_ps(coefs_ws[3]),
    _mm256_set1_ps(coefs_ws[4]),
    _mm256_set1_ps(coefs_ws[5])
};

__m256 VecLib::full_cos_AVX(const __m256& x) { 
    //performs horrible math tragedies to stich together a full 0-2pi sin out of a 0-pi/2 cosine
    //its mainly because cosine goes negative in the middle, which is just a pain.
    //as such its a bit slower, theoretically <42 clocks.

    //creates a triangular drive signal, and flips the middle part of the output.
    // the drive is       /\/\
    // the flip signal is _--_
    // the result is      \  /
    //                     \/

    __m256 reverse_mask = _mm256_cmp_ps(x, AVX_PI, _CMP_GE_OQ); //if x > pi, reset the drive signal and flip the output
    __m256 k = _mm256_and_ps( AVX_ABS_MASK,
        _mm256_sub_ps(x,
        _mm256_and_ps(AVX_2PI, reverse_mask)
    ));
    __m256 flip_mask = _mm256_cmp_ps(k, AVX_HALF_PI, _CMP_GE_OQ); //if x > pi, reset the drive signal and flip the output
    __m256 driver = _mm256_and_ps(AVX_ABS_MASK,
        _mm256_sub_ps(
            k,
            _mm256_and_ps(flip_mask, AVX_PI)
        )
    );
    return _mm256_xor_ps(
        _mm256_and_ps(flip_mask, AVX_SIGN_MASK),
        quarter_cosine_AVX(driver)
    );
}

__m256 VecLib::full_sin_AVX(__m256& x) { 
    const __m256 fpiavx = _mm256_set1_ps(fPI);
    x = _mm256_sub_ps(fpiavx,x);
    __m256 x2 = _mm256_mul_ps(x, x);
    __m256 v = _mm256_fmadd_ps(x2, coefs_ws_AVX[0], coefs_ws_AVX[1]);
    v = _mm256_fmadd_ps(v, x2, coefs_ws_AVX[2]);
    v = _mm256_fmadd_ps(v, x2, coefs_ws_AVX[3]);
    v = _mm256_fmadd_ps(v, x2, coefs_ws_AVX[4]);
    v = _mm256_fmadd_ps(v, x2, coefs_ws_AVX[5]);
    x = _mm256_mul_ps(x, v);
    return x;
}

__m256 VecLib::full_sin_AVX_(const __m256& x) {
    //theoretically this is <35 clocks! thats about one regular sine call.
    //there might be room to shave about 3 more off with better logic

    //this overall creates a triangular drive signal, and flips the second half of the output.
    // drive is       \/\/
    // flip signal is __--
    // result is      /\
    //                  \/

    __m256 flip_mask = _mm256_cmp_ps(x, AVX_PI, _CMP_GE_OQ); //if x > pi, reset the drive signal and flip the output
    __m256 inv_flip_mask = _mm256_xor_ps(flip_mask, AVX_PACKED_ONES);
    __m256 driver = _mm256_and_ps(AVX_ABS_MASK, //take the abs of the result
        _mm256_sub_ps(x,
            _mm256_or_ps(//select which value to subtract from x.
                _mm256_and_ps(flip_mask, AVX_3O2_PI),
                _mm256_and_ps(inv_flip_mask, AVX_HALF_PI)
            ))
    );
    return _mm256_xor_ps( //calculation and flipping of output
        _mm256_and_ps(flip_mask, AVX_SIGN_MASK),
        quarter_cosine_AVX(driver)
    );
}

__m256 VecLib::quarter_cosine_AVX(const __m256& x) {
    //this very closely approximates sine, and importantly 0 at 0,
    //and within like 5 decimals at 1 (its like 0.9999997 or smth). 
    //a lookup table _might_ be faster, but I sincerely doubt it.
    //theoretically less than 20 clocks! fmadd is insanely goated.
    __m256 x2 = _mm256_mul_ps(x,x);
    return
        _mm256_fmadd_ps(
            x2,
            _mm256_fmadd_ps(
                x2,
                _mm256_fmadd_ps(
                    x2,
                    _mm256_fmadd_ps(
                        x2,
                        coefs_AVX[0],
                        coefs_AVX[1]
                    ),
                    coefs_AVX[2]
                ),
                coefs_AVX[3]
            ),
            coefs_AVX[4]
        );
}
static const uint32_t PACKED_ONES = 0xffffffff;
static const uint32_t PACKED_ZEROS = 0x00000000;
static const uint32_t ABS_MASK = ~0x80000000;
static const uint32_t SIGN_MASK = 0x80000000;
#define BITWISE_SELECT(c,a,b) std::bit_cast<float>(((c) & std::bit_cast<uint32_t>(a)) | ((~(c)) & std::bit_cast<uint32_t>(b)))
#define BITWISE_SELECT0(c,a) std::bit_cast<float>((c)&std::bit_cast<uint32_t>(a))
#define BITWISE_ABS(a) std::bit_cast<float>(ABS_MASK&std::bit_cast<uint32_t>(a))
#define BITWISE_NEGATE_SELECT(c,a) std::bit_cast<float>(((c)&SIGN_MASK)^std::bit_cast<uint32_t>(a))
#define CREATE_BITMASK(c) (PACKED_ZEROS-(c))
float VecLib::full_sin(const float& x) {
    //return Lsin(x);
    __m128 xss = _mm_set_ss(fPI-x); 
    whole_sine_ss(xss);
    return _mm_cvtss_f32(xss);
}
float VecLib::full_sin_(const float& x) {
    __m128 xss = _mm_set_ss(x);
    __m128 piss = _mm_set_ss(fPI);
    __m128 sel_ps = _mm_set_ps(3.0f / 2.0f * fPI,0.0f,0.0f, fPI / 2.0f);
    __m128 abs_maskss = _mm_set_ss(std::bit_cast<float>(ABS_MASK));
    __m128 flip_mask = _mm_cmpgt_ss(xss, piss);
    __m128 driver_sub = _mm_permutevar_ps(sel_ps, std::bit_cast<__m128i>(flip_mask));
    xss = _mm_sub_ss(xss, driver_sub);
    xss = _mm_and_ps(xss, abs_maskss);
    quarter_cosine_ss(xss);
    __m128 xor_mask = _mm_andnot_ps(abs_maskss, flip_mask);
    xss = _mm_xor_ps(xss, xor_mask);
    return _mm_cvtss_f32(xss);
}
float VecLib::full_sin__(const float& x) { //this is just a readable variant
    //if (x > 2 * fPI || x < 0) throw std::exception();
    bool flip_mask = x > fPI;//CREATE_BITMASK(x > PI); //if x > pi, reset the drive signal and flip the output
    float driver_sub = 0.5f*fPI;
    if (flip_mask) driver_sub = 1.5f*fPI;
    float driver = abs(x - driver_sub); //take the abs of the result
    return std::bit_cast<float>((flip_mask << 31) ^ std::bit_cast<unsigned int>(quarter_cosine(driver)));
}

float VecLib::full_cos(const float& x) {
    //return Lcos(x);
    __m128 xss = _mm_set_ss(x-fPI);
    whole_cosine_ss(xss);
    return _mm_cvtss_f32(xss);
}
float VecLib::tan(const float& x) {
    return full_sin(x + fPI) / full_cos(x + fPI);
    //__m128 sine = _mm_set_ss(full_sin(x + fPI));
    //__m128 cosine = _mm_mul_ss(sine, sine);
    //__m128 one_ss = _mm_set_ss(1);
    //cosine = _mm_sub_ss(one_ss, cosine);
    //cosine = _mm_rsqrt_ss(cosine);
    //__m128 tan = _mm_mul_ss(sine, cosine);
    //return -_mm_cvtss_f32(tan);
}
float VecLib::full_cos_(const float& x) {
    //if (x > 2 * fPI || x < 0) throw std::exception();
    //return cos(x);
    //I am so sincerely sorry about this function. Manual bit manipulation was the only way to make this efficient
    unsigned int reverse_mask = CREATE_BITMASK(x>fPI); //if x > pi, reset the drive signal and flip the output
    float k_sub = BITWISE_SELECT0(reverse_mask,2.0f * fPI);
    float k = abs(x - k_sub);
    unsigned int flip_mask = CREATE_BITMASK(k >= fPI/2.0f); //if x > pi, reset the drive signal and flip the output
    float driver_sub = BITWISE_SELECT0(flip_mask,fPI);
    float driver = abs(k - driver_sub);
    return BITWISE_NEGATE_SELECT(flip_mask,quarter_cosine(driver));
}
void VecLib::quarter_cosine_ss(__m128& x) {
    x = _mm_mul_ss(x,x);
    __m128 v = _mm_fmadd_ss(x, coefs_ss[0], coefs_ss[1]);
    v = _mm_fmadd_ss(v, x, coefs_ss[2]);
    v = _mm_fmadd_ss(v, x, coefs_ss[3]);
    x = _mm_fmadd_ss(x, v, coefs_ss[4]);
}
void VecLib::whole_sine_ss(__m128& x) {
    __m128 x2 = _mm_mul_ss(x, x);
    __m128 v = _mm_fmadd_ss(x2, coefs_ws_ss[0], coefs_ws_ss[1]);
    v = _mm_fmadd_ss(v, x2, coefs_ws_ss[2]);
    v = _mm_fmadd_ss(v, x2, coefs_ws_ss[3]);
    v = _mm_fmadd_ss(v, x2, coefs_ws_ss[4]);
    x2 = _mm_fmadd_ss(v, x2, coefs_ws_ss[5]);
    x = _mm_mul_ss(x, x2);
}
void VecLib::whole_cosine_ss(__m128& x) {
    x = _mm_mul_ss(x, x);
    __m128 v = _mm_fmadd_ss(x, coefs_wc_ss[0], coefs_wc_ss[1]);
    v = _mm_fmadd_ss(v, x, coefs_wc_ss[2]);
    v = _mm_fmadd_ss(v, x, coefs_wc_ss[3]);
    v = _mm_fmadd_ss(v, x, coefs_wc_ss[4]);
    v = _mm_fmadd_ss(v, x, coefs_wc_ss[5]);
    x = _mm_fmadd_ss(x, v, coefs_wc_ss[6]);
}
float VecLib::quarter_cosine(const float& x) {
    float x2 = x*x;
    return coefs[4] + x2 * (coefs[3] + x2 * (coefs[2] + x2 * (coefs[1] + x2 * coefs[0])));
}
__m256 VecLib::sgn_fast(__m256& x) //credit to Peter Cordes https://stackoverflow.com/a/41353450
{
    __m256 negzero = _mm256_set1_ps(-0.0f);

    // using _mm_setzero_ps() here might actually be better without AVX, since xor-zeroing is as cheap as a copy but starts a new dependency chain
    //__m128 nonzero = _mm_cmpneq_ps(x, negzero);  // -0.0 == 0.0 in IEEE floating point

    __m256 x_signbit = _mm256_and_ps(x, negzero);
    return _mm256_or_ps(_mm256_set1_ps(1.0f), x_signbit);
}
void VecLib::cross_avx(const m256_vec3& v1, const m256_vec3& v2, m256_vec3& output) {
    __m256 c1 = _mm256_mul_ps(v1.Z, v2.Y);
    __m256 c2 = _mm256_mul_ps(v1.X, v2.Z);
    __m256 c3 = _mm256_mul_ps(v1.Y, v2.X);
    output.X = _mm256_fmsub_ps(v1.Y, v2.Z, c1);
    output.Y = _mm256_fmsub_ps(v1.Z, v2.X, c2);
    output.Z = _mm256_fmsub_ps(v1.X, v2.Y, c3);
}
void VecLib::dot_avx(const m256_vec3& v1, const m256_vec3& v2, __m256& output) {
    output = _mm256_fmadd_ps(
        v1.Z,
        v2.Z,
        _mm256_fmadd_ps(
            v1.Y,
            v2.Y,
            _mm256_mul_ps(
                v1.X,
                v2.X
            )
        )
    );
}
__m256 VecLib::inline_dot_avx(const m256_vec3& v1, const m256_vec3& v2) {
    __m256 output;
    dot_avx(v1, v2, output);
    return output;
}
void VecLib::dot_mul_avx(const m256_vec3& v1, const m256_vec3& v2, const __m256& v3, __m256& output) {
    output = _mm256_mul_ps(
        v3,
        _mm256_fmadd_ps(
            v1.Z,
            v2.Z,
            _mm256_fmadd_ps(
                v1.Y,
                v2.Y,
                _mm256_mul_ps(
                    v1.X,
                    v2.X
                )
            )
        )
    );
}
double VecLib::poly_acos(float x) {
    return (-0.69813170079773212 * x * x - 0.87266462599716477) * x + 1.5707963267948966;
}
float interpolatedLsin(float input) {
    if (input < 0) input = fPI - input;
    float v = (input * (LOOKUP_SIZE_TRIG / (fPI * 2.0f)));
    float l_v = v - floor(v);
    int i1 = ((int)v) & LOOKUP_TRIG_MASK;
    int i2 = (i1 + 1) & LOOKUP_TRIG_MASK;
    return ((1.0f - l_v) * sin_lookup[i1] + l_v * sin_lookup[i2]);
}
float VecLib::Lsin(float input) {
    return interpolatedLsin(input);
    if (input < 0) input = fPI - input;
    float v = (input * (LOOKUP_SIZE_TRIG / (fPI * 2.0f)));
    int i1 = ((int)round(v)) & LOOKUP_TRIG_MASK;
    return sin_lookup[i1];
}

float interpolatedLcos(float input) {
    if (input < 0) input = fPI - input;
    float v = (input * (LOOKUP_SIZE_TRIG / (fPI * 2.0f)));
    float l_v = v - floor(v);
    int i1 = ((int)v) & LOOKUP_TRIG_MASK;
    int i2 = (i1 + 1) & LOOKUP_TRIG_MASK;
    return ((1.0f - l_v) * cos_lookup[i1] + l_v * cos_lookup[i2]);
}

float VecLib::Lcos(float input) {
    return interpolatedLcos(input);
    if (input < 0) input = 2 * PI - input;
    int index = ((int)(input * (LOOKUP_SIZE_TRIG / (PI * 2)))) & LOOKUP_TRIG_MASK;
    return cos_lookup[index];
}
float VecLib::Lacos(float input) {
    int index = ((int)((input + 1) / 2.0 * LOOKUP_SIZE_TRIG)) & LOOKUP_TRIG_MASK;
    return arccos_lookup[index];
}
XYZ VecLib::gen_random_point_on_sphere(XorGen& G) {
    double r1 = G.fRand(0, 1);
    double r2 = G.fRand(0, 1);
    double f1 = Lacos(2 * r1 - 1) + PI / 2.0;
    double f2 = r2 * 2 * PI;
    return XYZ(-full_cos(f1) * full_cos(f2),- full_sin(f1), -full_cos(f1) * full_sin(f2));
}
//r1 is the vertical component, r2 is the rotational component
XYZ VecLib::generate_biased_random_hemi_v2(XorGen& G, float r1_min, float r1_max, float r2_min, float r2_max) {//https://www.rorydriscoll.com/2009/01/07/better-sampling/
    const float r1 = G.fRand(r1_min, r1_max);
    const float r2 = G.fRand(r2_min, r2_max);
    const float r = sqrt(r1);
    const float theta = 2 * PI * r2;

    const float x = r * VecLib::full_cos(theta);
    const float y = r * VecLib::full_sin(theta);

    return XYZ(x, sqrt(1 - r1), y);
}
XYZ VecLib::lookup_biased_random_hemi(XorGen& G) {
    int index = G.fRand(0, 1) * LOOKUP_BIASED_HEMI_MASK;
    return biased_hemi_lookup[index];
}
XYZ VecLib::generate_biased_random_hemi(XorGen& G) {
    XYZ pos = XYZ(0, 1, 0) + gen_random_point_on_sphere(G);
    return XYZ::slope(XYZ(0), pos);
}
void VecLib::prep(XorGen& G) {
    for (int i = 0; i < LOOKUP_SIZE_TRIG; i++) {
        double pi2_scalar = ((double)LOOKUP_SIZE_TRIG) / 2.0 / PI;
        double scalar2 = (double)LOOKUP_SIZE_TRIG / 2;
        sin_lookup[i] = sin((double)i / pi2_scalar);
        cos_lookup[i] = cos((double)i / pi2_scalar);
        arccos_lookup[i] = acos((double)i / scalar2 - 1);
    }
    for (int i = 0; i < LOOKUP_SIZE_BIASED_HEMI; i++) {
        biased_hemi_lookup[i] = generate_biased_random_hemi(G);
    }
}
XYZ VecLib::generate_unbiased_random_hemi(XorGen& G) {
    XYZ pos = gen_random_point_on_sphere(G);
    if (pos.Y < 0) pos.Y *= -1;
    return XYZ::slope(0, pos);
}
XYZ VecLib::biased_random_hemi(XorGen& G, float r1_min, float r1_max, float r2_min, float r2_max) {
    return generate_biased_random_hemi_v2(G, r1_min, r1_max, r2_min, r2_max);
}
XYZ VecLib::unbiased_random_hemi(XorGen& G) {
    return generate_unbiased_random_hemi(G);
}
XYZ VecLib::lookup_random_cone(XorGen& G, float spread) {
    float y_f = G.fRand(0, spread);
    float z_f = G.fRand(0, 2 * PI);
    return XYZ(full_sin(y_f) * full_cos(z_f), full_cos(y_f), full_sin(y_f) * full_sin(z_f));
}
XYZ VecLib::y_random_cone(XorGen& G, float spread) {
    float y = G.fRand(1, spread);
    float rot = G.fRand(0, PI);
    float f_y = sqrt(1 - y * y);
    return XYZ(f_y * cos(rot), y, f_y * sin(rot));
}
XYZ VecLib::aligned_random(XorGen& G, float spread, const Quat& r) {
    return Quat::applyRotation(lookup_random_cone(G, spread), r);
}
float VecLib::surface_area(const XYZ& max, const XYZ& min) {
    float l_x = max.X - min.X;
    float l_y = max.Y - min.Y;
    float l_z = max.Z - min.Z;
    return (l_x * l_y + l_x * l_z + l_y * l_z) * 2;
}
float VecLib::volume(const XYZ& max, const XYZ& min) {
    float l_x = max.X - min.X;
    float l_y = max.Y - min.Y;
    float l_z = max.Z - min.Z;
    return (l_x * l_y * l_z);
}
bool VecLib::between(const XYZ& p1, const XYZ& p2, const XYZ test) {
    return (
        (p1.X < test.X && test.X < p1.X) &&
        (p1.Y < test.Y && test.Y < p1.Y) &&
        (p1.Z < test.Z && test.Z < p1.Z)
        );
}
bool VecLib::betweenX(const XYZ& p1, const XYZ& p2, const XYZ test) {
    return (
        (p1.X < test.X && test.X < p1.X)
        );
}
bool VecLib::betweenY(const XYZ& p1, const XYZ& p2, const XYZ test) {
    return (
        (p1.Y < test.Y && test.Y < p1.Y)
        );
}
bool VecLib::betweenZ(const XYZ& p1, const XYZ& p2, const XYZ test) {
    return (
        (p1.Z < test.Z && test.Z < p1.Z)
        );
}
bool VecLib::volumeClip(const XYZ& max_1, const XYZ& min_1, const XYZ& max_2, const XYZ& min_2) {
    bool term1 = (max_1.X >= min_2.X && max_2.X >= min_1.X);
    bool term2 = (max_1.Y >= min_2.Y && max_2.Y >= min_1.Y);
    bool term3 = (max_1.Z >= min_2.Z && max_2.Z >= min_1.Z);
    return term1 && term2 && term3;
    /*bool term1 = (
        (((max_2.X < max_1.X) && (max_2.X > min_1.X)) || ((min_2.X < max_1.X) && (min_2.X > min_1.X))) &&
        (((max_2.Y < max_1.Y) && (max_2.Y > min_1.Y)) || ((min_2.Y < max_1.Y) && (min_2.Y > min_1.Y))) &&
        (((max_2.Z < max_1.Z) && (max_2.Z > min_1.Z)) || ((min_2.Z < max_1.Z) && (min_2.Z > min_1.Z)))
        );
    return term1;*/
}
bool VecLib::volumeContains(const XYZ& max, const XYZ& min, const XYZ& point) {
    bool t1 = (point.X<max.X && point.X>min.X);
    bool t2 = (point.Y<max.Y && point.Y>min.Y);
    bool t3 = (point.Z<max.Z && point.Z>min.Z);
    return t1 && t2 && t3;
}