
#include <string>
#include <vector>;
#include <math.h>;
#include <assert.h>
#include <malloc.h>
#include <immintrin.h>

import XYZ;
import Matrix;
import VecLib;
import XorRandom;
export module Materials;

#define PI 3.14159265358979323846264
#define USE_MORTON 0
using namespace std;

export template <typename T>
class Texture {
public:
    int resolution_x;
    int resolution_y;
    double U_multiplier;
    double V_multiplier;
    XY scale = XY(1,1);
    XY offset = XY(0,0);
    T* data;//was originally a vector, but the red tape surrounding vector made texture loading very slow in debug mode


    Texture() {
        data = nullptr;
    }
    Texture(T* data_ptr) {
        data = data_ptr;
    }
    Texture(unsigned char* img, int width, int height, int components, T (*converter)(float,float,float)) {
        allocate(width, height);
        unsigned int position = 0;
        unsigned int data_pos = 0;
        for (unsigned int y = 0; y < resolution_y; y++) {
            for (unsigned int x = 0; x < resolution_x; x++) {
                position+=components;
                float r = ((float)img[position + 0])/255.0;
                float g = ((float)img[position + 1]) / 255.0;
                float b = ((float)img[position + 2]) / 255.0;
                unsigned int address = to_memory(x, y);
                data[address] = (converter(r, g, b));
                data_pos++;
            }
        }
    }
    ~Texture() {
        if (data != nullptr) {
            free(data);
        }
    }
    Texture<T>* clone_shallow() {
        Texture<T>* tex = new Texture<T>();
        tex->resolution_x = resolution_x;
        tex->resolution_y = resolution_y;
        tex->set_transform(scale, offset);
        tex->data = data;
        return tex;
    }


    unsigned int direct_size(unsigned int rX, unsigned int rY) {
        return rX * rY;
    }
    unsigned int morton_size(unsigned int rX, unsigned rY) {
        int largest_axis = max(rX, rY);
        int power_2_axes = ceil(log2(largest_axis));
        return pow(pow(2, power_2_axes),2);
    }
    unsigned int get_memory_size(unsigned int rX, unsigned int rY) {
#if USE_MORTON 
        return morton_size(rX, rY);
#else          
        return direct_size(rX, rY);
#endif
    }
    void allocate(int rx, int ry) {
        resolution_x = rx;
        resolution_y = ry;
        unsigned int memory_size = get_memory_size(rx, ry);
        data = (T*)_aligned_malloc(sizeof(T) * memory_size, 64);
    }

    void set_transform(XY _scale, XY _offset) {
        scale = _scale;
        offset = _offset;
    }
    int size() {
        return resolution_x * resolution_y;
    }

    Texture<float>* export_channel(int channel) {
        Texture<float>* out = new Texture<float>();
        out->allocate(resolution_x, resolution_y);
        out->set_transform(scale, offset);
        int pos = 0;
        for (int y = 0; y < resolution_y; y++) {
            for (int x = 0; x < resolution_x; x++) {
                out->data[pos] = lookup(x, y)[channel];
                pos++;
            }
        }
        return out;
    }
    void prep() {
        U_multiplier = resolution_x*scale.X;
        V_multiplier = resolution_y*scale.Y;
    }
    unsigned int to_morton(unsigned int x, unsigned int y) {
        return _pdep_u64(x, 0x5555555555555555) | _pdep_u64(y, 0xaaaaaaaaaaaaaaaa);
    }
    unsigned int to_direct_address(unsigned int x, unsigned int y) {
        return y * resolution_x + x;
    }
    unsigned int to_memory(unsigned int x, unsigned int y) {
#if USE_MORTON
        return to_morton(x, y);
#else
        return to_direct_address(x, y);
#endif
    }

    T lookup(int x, int y) {
        x = x % resolution_x;
        y = y % resolution_y;
        if (x < 0) x = resolution_x + x;
        if (y < 0) y = resolution_y + y;
        unsigned int address = to_memory(x, y);
        return data[address];
    }
    T getUV(double U, double V) {
        V = 1 - V;
        int pos_x = (int)floor(U_multiplier * U)%resolution_x;
        int pos_y = (int)floor(V_multiplier * V)%resolution_y;
        return lookup(pos_x, pos_y);
    }
    T getPixelLinear(double U, double V) {
        double x_1 = U_multiplier * U - 0.5 * U_multiplier;
        double x_2 = U_multiplier * U + 0.5 * U_multiplier;
        double y_1 = V_multiplier * V - 0.5 * V_multiplier;
        double y_2 = V_multiplier * V + 0.5 * V_multiplier;
        int pos_x_1 = std::max(0, std::min(resolution_x - 1, (int)floor(x_1)));
        int pos_x_2 = std::max(0, std::min(resolution_x - 1, (int)floor(x_2)));
        int pos_y_1 = std::max(0, std::min(resolution_y - 1, (int)floor(y_1)));
        int pos_y_2 = std::max(0, std::min(resolution_y - 1, (int)floor(y_2)));
        T c1 = lookup(pos_x_1, pos_y_1);
        T c2 = lookup(pos_x_2, pos_y_1);
        T c3 = lookup(pos_x_1, pos_y_2);
        T c4 = lookup(pos_x_2, pos_y_2);
        double f1 = pos_x_2 - x_1;
        double f2 = x_2 - pos_x_2;
        double f3 = pos_y_2 - y_1;
        double f4 = y_2 - pos_y_2;
        T c_m1 = f1 * c1 + f2 * c2;
        T c_m2 = f1 * c3 + f2 * c4;
        return f3 * c_m1 + f4 * c_m2;
    }
};

export template <typename T>
class Parameter {
public:
    T static_value = 0;
    Texture<T>* texture = nullptr;
    Parameter(T value) : static_value(value) {}
    Parameter(Texture<T>* texture_ptr) : texture(texture_ptr) {}
    //Parameter(string file_name) : texture(new XYZTexture(file_name)) {}
    void set_static(T value) {
        static_value = value;
    }
    void prep() {
        if (texture != nullptr) {
            texture->prep();
        }
    }
    void set_texture(Texture<T>* texture_ptr, bool free = false) {
        if (texture != nullptr) {
            if (free) {
                delete texture;
            }
        }
        texture = texture_ptr;
    }
    void set_texture(string file_name, bool free = false) {
        //set_texture(new XYZTexture(file_name), free);
    }
    T getValue(double U = 0, double V = 0) {
        if (texture != nullptr) {
            return texture->getUV(U, V);
        }
        else {
            return static_value;
        }
    }
};

export class MaterialSample {
public:
    XYZ color;
    float roughness;
    float metallic;
    float specular;
    XYZ emissive;
    XYZ normal;

    float k = 0;
    float a_2 = 0;
    float spec_f = 0;
    float diff_f = 0;
    float diff_spread = 0;
    XYZ diff_c = XYZ();
    XYZ diff_t = XYZ();
    XYZ spec_color = XYZ();
    XYZ I_spec = XYZ();

    MaterialSample(XYZ _color, float _roughness, float _metallic, float _specular, XYZ _emissive, XYZ _normal) {
        color = _color;
        roughness = _roughness * 0.8 + 0.2;
        metallic = _metallic;
        specular = _specular;
        emissive = _emissive;
        normal = _normal;

        k = pow(roughness + 1, 2) / 8;
        a_2 = roughness * roughness;
        spec_color = get_specular_color();
        I_spec = 1 - get_fresnel_0();
        spec_f = get_specular_factor();
        diff_f = get_diffuse_factor();
        diff_c = get_diffuse_color();
        diff_t = diff_c * diff_f;
        float r_f = roughness - 0.2;
        diff_spread = (r_f * r_f) * PI / 2;
    }
    float get_diffuse_factor() const {
        return 1 - get_specular_factor();
    }
    XYZ get_fresnel_0() const {
        return XYZ::linear_mix(metallic, XYZ(0.04), color);
    }
    XYZ get_diffuse_reflectance() const {
        return color * (1 - metallic);
    }
    float get_specular_factor() const {
        return std::min(std::max(metallic, specular), (float)1.0);
    }
    XYZ get_specular_factor_v2(float dot_NI) const {
        return fast_fresnel(dot_NI);
    }
    XYZ get_diffuse_factor_v2(float dot_NI) const {
        return XYZ(1) - get_specular_factor_v2(dot_NI);
    }
    XYZ get_specular_color() const {
        return XYZ::linear_mix(metallic, specular * XYZ(1, 1, 1), color);
    }
    XYZ get_diffuse_color() const {
        return get_diffuse_reflectance() / PI;
    }
    XYZ get_fresnel(XYZ light_slope, XYZ normal) const {
        XYZ specular_color = get_fresnel_0();
        auto second_term = (1 - specular_color) * pow(1 - XYZ::dot(light_slope, normal), 5);
        return specular_color + second_term;
    }
    XYZ fast_fresnel(float dot_NI) const {
        float g = 1 - dot_NI;
        return spec_color+(1-spec_color) * (g * g * g * g * g);
    }
    float get_normal_distribution_beta(const XYZ& normal, const XYZ& half_vector) const {
        float a = roughness;
        float dot = XYZ::dot(normal, half_vector);
        float exponent = -(1 - pow(dot, 2)) / (pow(a * dot, 2));
        float base = 1 / (PI * pow(a, 2) * pow(dot, 4));
        float final = base * exp(exponent);

        return final;
    }
    float get_normal_distribution_GGXTR(const XYZ& normal, const XYZ& half_vector) const {
        float a = roughness;
        float dot = XYZ::dot(normal, half_vector);
        float final = (a * a)
            /
            (PI * pow((dot * dot) * (a * a - 1) + 1, 2));

        return final;
    }
    float fast_normal_dist(const float dot_NH) const {
        float a_2 = roughness * roughness;
        float g = dot_NH * dot_NH * (a_2 - 1) + 1;
        return dot_NH * (a_2) / (PI * g * g);
    }
    float geoSchlickGGX(const XYZ& normal, const XYZ& vector, float k) const {

        float dot = XYZ::dot(normal, vector);

        return dot / (dot * (1 - k) + k);
        //return 1;

    }
    float fastGeo_both(const float dot_NO, const float dot_NI) const {
        return dot_NO / (dot_NO * (1 - k) + k) * dot_NI / (dot_NI * (1 - k) + k);
    }
    XYZ diffuse_BRDF() const {

        //return get_diffuse_factor_v2(dot_NI) * get_diffuse_reflectance();
        return get_diffuse_reflectance();
    }
    XYZ diffuse_BRDF(float dot_NI) const {

        //return get_diffuse_factor_v2(dot_NI) * get_diffuse_reflectance();
        return diffuse_BRDF();
    }
    XYZ diffuse_BRDF_unweighted(float dot_NI) const {
        return get_diffuse_factor_v2(dot_NI) * get_diffuse_reflectance();
    }
    XYZ specular_BRDF(const XYZ& normal, const XYZ& input_slope, XYZ& output_slope) const {
        float dot_NI = XYZ::dot(normal, input_slope);
        float dot_NO = XYZ::dot(normal, output_slope);
        if (dot_NI <= 0 || dot_NO <= 0) {
            return XYZ(0, 0, 0);
        }
        XYZ half_vector = XYZ::normalize(XYZ::add(input_slope, output_slope));
        float dot_HO = XYZ::dot(half_vector, output_slope);
        float dot_NH = XYZ::dot(normal, half_vector);
        XYZ fresnel = fast_fresnel(dot_HO);
        float geo = fastGeo_both(dot_NO, dot_NI);
        float normal_dist = fast_normal_dist(dot_NH);
        float divisor = 4 * dot_NO;

        return geo * normal_dist * fresnel;// / divisor;
    }
    XYZ fast_BRDF_co(const XYZ& normal, const XYZ& input_slope, XYZ& output_slope) const {
        float dot_NI = XYZ::dot(normal, input_slope);
        float dot_NO = XYZ::dot(normal, output_slope);
        if (dot_NI <= 0 || dot_NO <= 0) {
            return XYZ(0, 0, 0);
        }
        XYZ half_vector = XYZ::normalize(XYZ::add(input_slope, output_slope));
        float dot_HO = XYZ::dot(half_vector, output_slope);
        float dot_NH = XYZ::dot(normal, half_vector);
        XYZ fresnel = fast_fresnel(dot_HO);
        float geo = fastGeo_both(dot_NO, dot_NI);
        float normal_dist = fast_normal_dist(dot_NH);
        float divisor = 4 * dot_NI * dot_NO;

        XYZ specular_return = (geo * normal_dist * fresnel) / divisor;

        XYZ light_remainer = XYZ(1) - specular_return;

        return (specular_return + light_remainer*diffuse_BRDF()) * dot_NI;

    }
    XYZ random_bounce(XorGen& G, const Matrix3x3& diffuse_rotation, const Matrix3x3& reflection_rotation) const {
        float prob = G.fRand(0, 1);
        if (prob < diff_f) {
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
        double r_v = (roughness - 0.2) / 0.8;
        float dist = 1 / (pow(r_v+1,4) - 1)-0.05;
        dir.Y += dist;
        XYZ slope = XYZ::normalize(dir);
        return Matrix3x3::applyRotationMatrix(slope, reflection_rotation);
    }

    XYZ calculate_emissions() const {
        return emissive;
    }
};

export class Material {
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

    int UV_map_index = 0;
    bool use_normals = false;

    Material() {}

    void prep() {
        color.prep();
        roughness.prep();
        metallic.prep();
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
            emissive.getValue(U, V),
            normal.getValue(U, V)
        );
    }
};