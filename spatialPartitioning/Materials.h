#pragma once

#include <string>
#include <vector>;
#include <math.h>;
#include <assert.h>
#include "Matrix.h";
#include "Texture.h"
#include "commons.h"



class MaterialSample {
public:
    XYZ color;
    float roughness;
    float metallic;
    float specularFactor;
    XYZ emissive;
    XYZ normalPertubation;
    float transmission;
    float IOR;

    XYZ normal = normal;
    XYZ out;
    XYZ in;
    XYZ half;

    float n1;
    float n2;
    float R0;
    XYZ R0C;

    float i_metal;
    float i_roughness;

    float sigma; // oren nayar standard deviation
    XYZ albedo;

    float k = 0;
    float a = 0;
    float a_2 = 0;
    float f_spec = 0;
    float f_diff = 0;
    float f_tran = 0;
    float diff_spread = 0;
    XYZ diff_c = XYZ();
    XYZ diff_t = XYZ();
    XYZ spec_color = XYZ();

    float dot_NO;
    float dot_NI;
    float dot_NH;
    float dot_HI;
    float dot_HO;

    float Ti;
    float To;

    float dot_AD; //azimuth delta dot

    bool use_fresnel = true;

    MaterialSample(XYZ _color, float _roughness, float _metallic, float _specular, float _transmission, float _IOR, XYZ _emissive, XYZ _normal) {
        color = _color;
        if (isnan(color.X)) {
            throw std::exception();
        }
        roughness = _roughness;
        metallic = _metallic;
        specularFactor = _specular;
        emissive = _emissive;
        normalPertubation = _normal;
        transmission = _transmission;
        IOR = _IOR;
    }
    void config(XYZ _out, XYZ _normal, float _n1, float _n2) {
        out = _out;
        normal = _normal;

        dot_NO = std::min(XYZ::dot(out,normal),1.0f); //might not be a bad idea to use avx here

        n1 = _n1;
        n2 = _n2;

        i_metal = 1*metallic;             //interpolate metallic for fresnel curve
        i_roughness = roughness*0.999 + 0.001;  //interpolate roughness to conform to "normal" values

        sigma = roughness * 0.3;

        To = acos(dot_NO);

        R0 = get_R0();   //fresnel 0
        R0C = get_R0C(); //fresnel 0 with color

        //k = pow(i_roughness + 1, 2) / 8;
        
        k = i_roughness * i_roughness / 2;
        a = i_roughness;
        a_2 = a*a;

        f_spec = get_specular_factor();
        f_diff = get_diffuse_factor();
        f_tran = get_transmissive_factor();
        
        diff_c = get_diffuse_color();
        diff_t = diff_c * f_diff;
        float r_f = roughness - 0.2;
        diff_spread = (r_f * r_f) * PI / 2;
    }
    float acos(float x) {
        float negate = float(x < 0);
        x = abs(x);
        float ret = -0.0187293;
        ret = ret * x;
        ret = ret + 0.0742610;
        ret = ret * x;
        ret = ret - 0.2121144;
        ret = ret * x;
        ret = ret + 1.5707288;
        ret = ret * sqrt(1.0 - x);
        ret = ret - 2 * negate * ret;
        return negate * 3.14159265358979 + ret;
    }
    void set_input(XYZ _in) {
        in = _in;
        half = XYZ::normalize((out + in) / 2);


        dot_NI = std::min(1.0f,XYZ::dot(in, normal));
        dot_NH = std::min(1.0f,XYZ::dot(half, normal));
        dot_HI = std::min(1.0f,XYZ::dot(half, in));
        dot_HO = std::min(1.0f,XYZ::dot(half, out));



        Ti = VecLib::Lacos(dot_NI);

        XYZ projected_in = XYZ::normalize(in - normal * XYZ::dot(in, normal));
        XYZ projected_out = XYZ::normalize(out - normal * XYZ::dot(out, normal));

        dot_AD = XYZ::dot(projected_in, projected_out);
    }
    float get_R0() {
        float spec = pow((n1 - n2) / (n1 + n2), 2);
        return ((1 - i_metal) * spec + i_metal);
    }
    XYZ get_R0C() {
        return XYZ::linear_mix(metallic, XYZ(R0), i_metal*color);
    }
    float get_specular_factor() const {
        float g = 1 - dot_NO;
        return R0 + (1 - R0) * (g * g * g * g * g);
    }
    float get_diffuse_factor() {
        return (1 - f_spec) * (1 - transmission);
    }
    float get_transmissive_factor() {
        return (1 - f_diff - f_spec);
    }
    XYZ get_diffuse_reflectance() const {
        return color * (1 - metallic);
    }
    XYZ get_specular_factor_v2(float dot_NI) const {
        return fast_fresnel(dot_NI);
    }
    XYZ get_diffuse_factor_v2(float dot_NI) const {
        return XYZ(1) - get_specular_factor_v2(dot_NI);
    }
    XYZ get_diffuse_color() const {
        return get_diffuse_reflectance() / PI;
    }
    
    XYZ fast_fresnel(float dot_NI) const {
        float g = 1 - dot_NI;
        return R0C + (1 - R0C) * (g * g * g * g * g);
    }
    float fast_fresnel(XYZ normal, XYZ vector) const{
        float g = 1 - XYZ::dot(normal,vector);
        return R0 + (1 - R0) * (g * g * g * g * g);
    }
    float get_normal_distribution_GGXTR(const XYZ& normal, const XYZ& half_vector) const {
        float dot = XYZ::dot(normal, half_vector);
        float final = (a * a)
            /
            (PI * pow((dot * dot) * (a * a - 1) + 1, 2));

        return final;   
    }
    float normal_dist_ggx() const {
        float c = dot_NH;
        float s = sqrt(1 - dot_NH * dot_NH);
        return pow(1 + pow(s / c, 2) / a_2, -2);
    }
    float fast_normal_dist(const float dot_NH) const {
        float g = dot_NH * dot_NH * (a_2 - 1) + 1;
        return dot_NH * (a_2) / (PI * g * g);
    }
    float geoSchlickGGX(const XYZ& normal, const XYZ& vector, float k) const {

        float dot = XYZ::dot(normal, vector);

        return dot / (dot * (1 - k) + k);
        //return 1;

    }
    float unkownG1(float dot) const {
        float c = dot;
        float s = sqrt(1 - dot * dot);
        return (-1 + sqrt(1 + a_2 * a_2 * pow(s / c, 2))) / 2;
    }
    float G1_poly(const XYZ& vec) const { //https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf pg 7
        float c = std::min(XYZ::dot(vec, normal),0.99f);
        float s = sqrt(1 - c * c);
        float _a = 1/(a * s / c);
        float v = (3.535 * _a + 2.181 * _a * _a) / (1 + 2.276 * _a + 2.577 * _a * _a);
        if (_a < 1.6) return v;
        else return 1;
    }
    float G1(const XYZ& vec) const {
        float c = std::min(XYZ::dot(vec, normal), 0.99f);
        float s_2 = 1 - c * c;
        float t_2 = s_2 / (c * c);
        return 1 / (1 + (-1 + sqrt(1 + a * a * t_2)) / 2);
    }
    float G_GGX() const {
        return G1(out) * G1(in);
        //return 1 / (1 + G1(dot_HO) + G1(dot_HI));
    }
    float fastGeo_both(const float dot_NO, const float dot_NI) const {
        return dot_NO / (dot_NO * (1 - k) + k) * dot_NI / (dot_NI * (1 - k) + k);
    }
    XYZ lambertian_diffuse_BRDF() const {
        return color / PI;
    }
    XYZ oren_nayar_diffuse_BRDF() const {
        float s = sigma;
        float s_2 = sigma * sigma;

        float a = std::max(Ti, To);
        float b = std::min(Ti, To);

        float A = 1 - 0.5 * s_2 / (s_2 + 0.33);
        float B = 0.45 * s_2 / (s_2 + 0.09) * (VecLib::full_sin(a) - ((dot_AD >= 0) * pow(2 * b / PI, 3)));
        float C = 0.125 * s_2 / (s_2 + 0.09) * pow(4 * a * b / PI, 2);

        XYZ modifier = color / PI;
        float L1 = (A + B * std::max(0.0f, dot_AD) * VecLib::tan(b) + C * (1 - abs(dot_AD)) * VecLib::tan((a + b) / 2));
        float L2 = 0.17 * s_2 / (s_2 + 0.13) * (1 - dot_AD) * pow(2 * b / PI, 2);

        XYZ final = modifier * L1 + color * modifier * L2;
        return final;
    }
    XYZ diffuse_BRDF() const {
        return oren_nayar_diffuse_BRDF();
        //return lambertian_diffuse_BRDF();
    }
    XYZ specular_BRDF() const {
        if (dot_NH < 0 || dot_NO < 0 || dot_NI < 0 || dot_HI < 0 || dot_HO < 0) return XYZ(0, 0, 0);
        XYZ fresnel = (R0C*(!use_fresnel)) + use_fresnel * fast_fresnel(dot_HO);
        float geo = G_GGX();
        float normal_dist = normal_dist_ggx();//fast_normal_dist(dot_NH);
        float divisor = 4 * dot_NO * dot_NI;

        return geo * normal_dist * fresnel / divisor;// / divisor;
    }
    XYZ PDF_specular_BRDF() const {
        if (dot_NH < 0 || dot_NO < 0 || dot_NI < 0 || dot_HI < 0 || dot_HO < 0) return XYZ(0, 0, 0);
        XYZ tint = XYZ::linear_mix(metallic, XYZ(1, 1, 1), color);
        XYZ weight = G_GGX()/G1(in);
        return tint * weight;
    }
    XYZ fast_BRDF_co(bool apply_scalar) const {
        XYZ specular_return = specular_BRDF();
        XYZ light_remainder = XYZ(1) - specular_return;
        XYZ diffuse_return = light_remainder*diffuse_BRDF();

        return (specular_return+diffuse_return) * ((apply_scalar)?dot_NI:1);

    }
    XYZ PDF_BRDF() const{
        XYZ specular_return = PDF_specular_BRDF();
        XYZ light_remainder = XYZ(1) - specular_return;
        XYZ diffuse_return = light_remainder * diffuse_BRDF();

        return (specular_return + diffuse_return);
    }
    XYZ random_bounce(XorGen& G, const Matrix3x3& diffuse_rotation, const Matrix3x3& reflection_rotation) const {
        float prob = G.fRand(0, 1);
        if (prob < f_diff) {
            return Matrix3x3::aligned_random(G, PI / 2, diffuse_rotation);
        }
        else {
            return Matrix3x3::aligned_random(G, diff_spread, reflection_rotation);
        }
    }
    XYZ biased_diffuse_bounce(XorGen& G, const Matrix3x3& diffuse_rotation, float r1_min = 0, float r1_max = 1, float r2_min = 0, float r2_max = 1) const {
        return Matrix3x3::aligned_biased_hemi(G, diffuse_rotation, r1_min, r1_max, r2_min, r2_max);
    }
    XYZ unbiased_diffuse_bounce(XorGen& G, const Matrix3x3& diffuse_rotation) const {
        return Matrix3x3::aligned_unbiased_hemi(G, diffuse_rotation);
    }
    XYZ reflective_bounce(XorGen& G, const Matrix3x3& reflection_rotation) const {
        XYZ dir = VecLib::generate_unbiased_random_hemi(G);
        double r_v = roughness;
        float dist = 1 / (pow(r_v+1,5) - 1)-0.05;
        dir.Y += dist;
        XYZ slope = XYZ::normalize(dir);
        return Matrix3x3::applyRotationMatrix(slope, reflection_rotation);
    }

    XYZ calculate_emissions() const {
        return emissive;
    }
};

class Material {
public:
    Parameter<XYZ>   color        = Parameter<XYZ>(0);
    Parameter<float> roughness    = Parameter<float>(0.1);
    Parameter<float> metallic     = Parameter<float>(0.);
    Parameter<float> specular     = Parameter<float>(0.);
    Parameter<XYZ>   RMS          = Parameter<XYZ>(XYZ(0.1,0,0));
    Parameter<XYZ>   emissive     = Parameter<XYZ>(0);
    Parameter<XYZ>   normal       = Parameter<XYZ>(0);
    Parameter<float> transmission = Parameter<float>(0.);
    Parameter<float> IOR          = Parameter<float>(1);
    Parameter<float> alpha        = Parameter<float>(1);

    std::string name;

    int UV_map_index = 0;
    bool use_normals = false;

    Material() {}

    void prep() {
        color.prep();
        roughness.prep();
        metallic.prep();
        transmission.prep();
        specular.prep();
        emissive.prep();
        normal.prep();
    }
    MaterialSample sample_UV(float U, float V) {
        return MaterialSample(
            color.getValue(U, V),
            roughness.getValue(U, V),
            metallic.getValue(U, V),
            specular.getValue(U, V),
            transmission.getValue(U,V),
            IOR.getValue(U,V),
            emissive.getValue(U, V),
            normal.getValue(U, V)
        );
    }
};

class AVX_MaterialSample {
public:
    m256_vec3 color;
    __m256 roughness;
    __m256 metallic;
    __m256 specularFactor;
    m256_vec3 emissive;
    m256_vec3 normalPertubation;
    __m256 transmission;
    __m256 IOR;
    __m256 use_normal;

    __m256 fresnel_0;
    __m256 max_fresnel;
    void R0() {
        __m256 spec = _mm256_div_ps( // n1-n2 / n1+n2
            _mm256_sub_ps(AVX_ONES, IOR),
            _mm256_add_ps(AVX_ONES, IOR)
        );
        fresnel_0 = _mm256_mul_ps(spec, spec); //square it
    }
    __m256 G1(const __m256& dot) const {
        __m256 c_2 = _mm256_mul_ps(dot, dot);
        __m256 a_2 = _mm256_mul_ps(roughness, roughness);
        __m256 s_2 = _mm256_sub_ps(AVX_ONES, c_2);
        __m256 t_2 = _mm256_div_ps(c_2, s_2);
        __m256 g1rt = _mm256_sqrt_ps(_mm256_fmadd_ps(a_2, t_2, AVX_ONES));
        __m256 denom = _mm256_add_ps(AVX_ONES, _mm256_fmsubadd_ps(g1rt, AVX_HALFS, AVX_HALFS));
        __m256 g1final = _mm256_div_ps(AVX_ONES, denom);
        return g1final;
    }
    void PDF_specular_BRDF(m256_vec3& outputs, const __m256& dot_NO, const __m256& dot_NI) {
        
        R0();
        __m256 f_a = _mm256_sub_ps(AVX_ONES, metallic);
        f_a = _mm256_mul_ps(fresnel_0, f_a);
        outputs.X = _mm256_add_ps(f_a, _mm256_mul_ps(metallic, color.X));
        outputs.Y = _mm256_add_ps(f_a, _mm256_mul_ps(metallic, color.Y));
        outputs.Z = _mm256_add_ps(f_a, _mm256_mul_ps(metallic, color.Z)); //calculating metallic fresnel tint
        max_fresnel = _mm256_mul_ps( //gets the max color channel, and multiplies it by the fresnel value
            _mm256_max_ps(
                _mm256_max_ps(
                    outputs.X,
                    outputs.Y
                ),
                outputs.Z
            ),
            fresnel_0
        );
        outputs.X = _mm256_div_ps(outputs.X, max_fresnel);//normalize fresnel return such that the max is 1
        outputs.Y = _mm256_div_ps(outputs.X, max_fresnel);
        outputs.Z = _mm256_div_ps(outputs.X, max_fresnel);
        __m256 g2 = _mm256_mul_ps(//geometric term
            G1(dot_NO),
            G1(dot_NI)
        );
        outputs.X = _mm256_mul_ps(outputs.X, g2);
        outputs.Y = _mm256_mul_ps(outputs.Y, g2);
        outputs.Z = _mm256_mul_ps(outputs.Z, g2);
    }
    void lambertian_BRDF(m256_vec3& outputs) {
        outputs.X = _mm256_div_ps(color.X, AVX_PI);
        outputs.Y = _mm256_div_ps(color.Y, AVX_PI);
        outputs.Z = _mm256_div_ps(color.Z, AVX_PI);
    }
    AVX_MaterialSample(Material** materials, float* UV_U, float* UV_V) {
        for (int i = 0; i < 8; i++) color.set(i, materials[i]->color.getValue(UV_U[i], UV_V[i]));
        for (int i = 0; i < 8; i++) emissive.set(i, materials[i]->emissive.getValue(UV_U[i], UV_V[i]));
        for (int i = 0; i < 8; i++) normalPertubation.set(i, materials[i]->normal.getValue(UV_U[i], UV_V[i]));
        for (int i = 0; i < 8; i++) ((float*)&roughness)[i] =  materials[i]->roughness.getValue(UV_U[i], UV_V[i]);
        for (int i = 0; i < 8; i++) ((float*)&metallic)[i]  = materials[i]->metallic.getValue(UV_U[i], UV_V[i]);
        for (int i = 0; i < 8; i++) ((float*)&specularFactor)[i] = materials[i]->specular.getValue(UV_U[i], UV_V[i]);
        for (int i = 0; i < 8; i++) ((float*)&IOR)[i] = materials[i]->IOR.getValue(UV_U[i], UV_V[i]);
        for (int i = 0; i < 8; i++) ((float*)&use_normal)[i] = std::bit_cast<float>(materials[i]->use_normals * 0xffffffff);
        roughness = _mm256_fmadd_ps(
            _mm256_set1_ps(0.999),
            roughness,
            _mm256_set1_ps(0.001)
        );
    }
    void transfer(int origin, int target) {
        color.transfer(origin, target);
        emissive.transfer(origin, target);
        normalPertubation.transfer(origin, target);
        ((float*)&roughness)[target] = ((float*)&roughness)[origin];
        ((float*)&metallic)[target] = ((float*)&metallic)[origin];
        ((float*)&specularFactor)[target] = ((float*)&specularFactor)[origin];
        ((float*)&IOR)[target] = ((float*)&IOR)[origin];
        ((float*)&use_normal)[target] = ((float*)&use_normal)[origin];

    }
};