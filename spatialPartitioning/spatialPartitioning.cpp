﻿// spatialPartitioning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _AMD64_ true
#define _WIN64 true
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>

#include <tuple>
#include <math.h>
#include <format>
#include <string>
#include <ctime>
#include <chrono>
#include <functional>
#include <fstream>

#include <SFML/Graphics.hpp>
#include <OpenColorIO/OpenColorIO.h>
namespace OCIO = OCIO_NAMESPACE;

#include <thread>
#include <vector>
#include <array>
#include <unordered_set>
#include <set>
#include <queue>
#include <processthreadsapi.h>

#include <cassert>
#include <immintrin.h>

#define TINYGLTF_IMPLEMENTATION
#include "tiny_gltf.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//#include "imgui-master/imgui.h"
//#include "imgui-master/backends/imgui_impl_win32.h"
//#include "imgui-master/backends/imgui_impl_dx12.h"

import XYZ;
import XorRandom;
import VecLib;
import Matrix;
import Materials;

#define assertm(exp, msg) assert(((void)msg, exp))
#define SMALL 0.001
#define NEAR_THRESHOLD 0.000001
#define SCENE_BOUNDS 10000
#define RED_MASK 255<<16
#define GREEN_MASK 255<<8
#define BLUE_MASK 255

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define PI 3.14159265358979323846264

#define kEpsilon 0.000000001

#define DRAW_REFLECTIVE_COEFFICIENTS 0
#define DRAW_DIFFUSE_COEFFICIENTS 0
#define DRAW_LIGHT_EXPOSURE 0
#define DRAW_COLOR 0
#define DEPTH_VIEW 0
#define DRAW_UV 0
#define DRAW_BOUNCE_DIRECTION 0
#define DRAW_NORMAL 0
#define DRAW_NORMAL_DELTAS 0
#define DRAW_NORMAL_FACING 0
#define DRAW_EMISSIVE 0


#define DRAW_LIGHTS   (1 && !DRAW_REFLECTIVE_COEFFICIENTS && !DRAW_DIFFUSE_COEFFICIENTS)
#define DRAW_DIFFUSE  (1 && !DRAW_REFLECTIVE_COEFFICIENTS)
#define DRAW_SPECULAR (0 && !DRAW_DIFFUSE_COEFFICIENTS)

#define USE_ADVANCED_BVH false

//defines program precision

#define decimal float

#define USE_AVX_BVH 1
#define USE_AVX_TRI (1 && USE_AVX_BVH)
#define ENABLE_MIS  0
#define USE_PACKET_TRAVERSAL 1

#if USE_AVX_TRI
#define LEAF_SIZE 8
#define PENALIZE_UNFILLED_LEAFS 1
#else
#define LEAF_SIZE 2
#define PENALIZE_UNFILLED_LEAFS 0
#endif


using namespace std;
using time_point = chrono::steady_clock::time_point;


static thread_local const __m256 AVX_ZEROS = _mm256_set1_ps(0);
int one = 1;
float float_one = *((float*)&one);
static thread_local const __m256 AVX_FINT_ONES = _mm256_set1_ps(float_one);


static thread_local XorGen gen;


/*Xorshiro256+ pseudorandom end*/


decimal clamp(decimal value, decimal low, decimal high) {
    return max(min(value, high), low);
}

decimal lerp(decimal a, decimal b, decimal factor) {
    return a + factor * (b - a);
}

vector<string> eng_levels = vector<string>{
    "k",
    "M",
    "B",
    "T",
    "qu",
    "Qu"
};

string padString(string s, string padder, int size) {
    string out = s;
    while (out.size() < size) {
        out += padder;
    }
    return out;
}

string iToFixedLength(int input, int length, string fill = " ") {
    string in = to_string(input);
    int spaces = length - in.size();
    for (int i = 0; i < spaces;i++) {
        in = fill + in;
    }
    return in;
}
string intToEng(int inum, bool force_small = false) {
    double num = inum;
    string prefix = "";
    if (num < 1000) {
        return iToFixedLength(inum, 4) + " ";
    }
    for (auto level : eng_levels) {
        if (num < 1000) {
            string number = to_string(num);
            if (num >= 100) {
                number = number.substr(0, 3);
                if (!force_small) {
                    number = " " + number;
                }
            }
            else {
                if (!force_small) {
                    number = number.substr(0, 4);
                }
                else {
                    number = number.substr(0, 3);
                }
            }
            return number + prefix;
        }
        num /= 1000;
        prefix = level;
    }
}

vector<string> split_string(string input, char splitter) {
    vector<string> output;
    string s;
    int size = input.size();
    for (int i = 0; i < size; i++) {
        char c = input[i];
        if (c == splitter) {
            output.push_back(s);
            s = "";
        }
        else {
            s += c;
        }
    }
    output.push_back(s);
    return output;
}

struct CastResults {
    static const decimal default_distance;
    XYZ normal;
    Material* material;
    decimal distance;
    XYZ TUV = XYZ(-1);
    XY texUV = XY(-1);
    int parent_index = 0;
    int tri_index = 0;
    CastResults() : normal(XYZ()), material(nullptr), distance(default_distance) {};
    CastResults(XYZ _normal, Material* _mat) : normal(_normal), material(_mat), distance(default_distance) {}
};

const decimal CastResults::default_distance = 99999999999999999;

class Primitive {
public:
    Material* material;
    //pair<XYZ, XYZ> bounds;
    XYZ origin;
    int obj_type = 0;
    Primitive(Material* _material) :material(_material) {}
    virtual pair<XYZ, XYZ> get_bounds() { return pair<XYZ, XYZ>(); }

    virtual bool check_backface(XYZ& position) { return false; }

    virtual decimal distance(XYZ& position) {
        cout << "[ERROR] Default distance function called! Something is misconfigured!" << endl;
        return 0;
    } //how to fix if this gets called: alcohol and crying
    virtual Material* material_properties_at(XYZ) {
        return material;
    }
};

class Light {
public:
    XYZ position;
    XYZ rotation;
    XYZ emission;
    virtual bool canVector(XYZ& pos2) {
        return true;
    }
    virtual XYZ vectorTo(XYZ& pos2) {
        return XYZ::slope(pos2, position);
    }
    virtual bool intersects(XYZ& pos, const XYZ& slope, CastResults& res) {
        throw exception("This function is not implemented");
        return 0;
    }
};

class SunLight :public Light {
public:
    XYZ vectorTo(XYZ& pos2) override {
        return rotation;
    }
    bool intersects(XYZ& pos, const XYZ& slope, CastResults& res) override {
        if (res.distance > 9999) {
            res.distance = 9999;
            return true;
        }
        return false;
    }
};

class PointLight : public Light {
public:
    XYZ radius;
    XYZ vectorTo(XYZ& pos2) override {
        return XYZ::slope(pos2, position);
    }
};


struct PackagedTri { //reduced memory footprint to improve cache perf
    XYZ p1;
    XYZ normal;
    XYZ n1;
    XYZ n2;
    XYZ n3;
    XYZ t1;
    XYZ t2;
    XYZ t3;
    XYZ b1;
    XYZ b2;
    XYZ b3;
    XYZ p1p3;
    XYZ p1p2;
    XY UV_basis;
    XY U_delta;
    XY V_delta;
    Material* material;
    XYZ origin_offset;
    PackagedTri() {}
    PackagedTri(const XYZ& _p1, const XYZ& _p2, const XYZ& _p3, XY UV_base, XY UV_U, XY UV_V, Material* _material) {
        material = _material;
        p1 = _p1;
        p1p3 = _p3 - p1;
        p1p2 = _p2 - p1;
        U_delta = UV_U - UV_base;
        V_delta = UV_V - UV_base;
        UV_basis = UV_base;
        normal = XYZ::normalize(XYZ::cross(p1p3, p1p2));
        origin_offset = XYZ::dot((_p1 + _p2 + _p3) / 3, normal) * normal;
    }
    XYZ get_normal(float u, float v) {
        float k = 1 - u - v;
        return  XYZ::normalize(n1 * k + n2 * u + n3 * v);
    }
    XYZ get_tangent(float u, float v) {
        float k = 1 - u - v;
        return  XYZ::normalize(t1 * k + t2 * u + t3 * v);
    }
    XYZ get_bitangent(float u, float v) {
        float k = 1 - u - v;
        return  XYZ::normalize(b1 * k + b2 * u + b3 * v);
    }
    void calculate_tangents() {
        t1 = XYZ::normalize(p1p3 * U_delta.X - p1p2 * V_delta.X);;
        t2 = t1;
        t3 = t1;
    }
    void calculate_bitangents() {
        b1 = XYZ::normalize(p1p2 * V_delta.Y - p1p3 * U_delta.Y);
        b2 = b1;
        b3 = b1;
    }
    void reorthogonalize() {
        t1 = XYZ::normalize(t1-n1*XYZ::dot(t1, n1));
        t2 = XYZ::normalize(t2-n2*XYZ::dot(t2, n1));
        t3 = XYZ::normalize(t3-n3*XYZ::dot(t3, n1));
        b1 = XYZ::cross(n1, t1);
        b2 = XYZ::cross(n2, t2);
        b3 = XYZ::cross(n3, t3);
    }
};

struct PTri_AVX { //literally 0.7kB of memory per object. What have I done.
    m256_vec3 p1;
    m256_vec3 p1p2;
    m256_vec3 p1p3;
    XYZ normal[8];
    XYZ origin_offset[8];
    Material* materials[8];
    int index = 0;
    char size = 0;
    PTri_AVX(vector<PackagedTri> tris) {
        vector<XYZ> _p1;
        vector<XYZ> _p1p2;
        vector<XYZ> _p1p3;
        vector<XYZ> _normal;
        vector<XYZ> _origin_offset;
        size = tris.size();
        for (int i = 0; i < tris.size() && i < 8; i++) {
            PackagedTri& tri = tris[i];
            _p1.push_back(tri.p1);
            _p1p2.push_back(tri.p1p2);
            _p1p3.push_back(tri.p1p3);
            _normal.push_back(tri.normal);
            _origin_offset.push_back(tri.origin_offset);
            materials[i] = tri.material;
        }
        p1 = m256_vec3(_p1);
        p1p2 = m256_vec3(_p1p2);
        p1p3 = m256_vec3(_p1p3);
        for (int i = 0; i < tris.size(); i++) {
            origin_offset[i] = _origin_offset[i];
            normal[i] = _normal[i];
        }
    }
    int get_index(int i) const {
        return index + i;
    }
    static void intersection_check(const PTri_AVX& T, const m256_vec3& position, const m256_vec3& slope, m256_vec3& output) { //my FUCKING cache 
        m256_vec3 pvec; VecLib::cross_avx(slope, T.p1p3, pvec);
        __m256 det;     VecLib::dot_avx(T.p1p2, pvec, det);

        __m256 invDet = _mm256_div_ps(_mm256_set1_ps(1), det);

        m256_vec3 tvec; m256_vec3::sub(position, T.p1, tvec);
        m256_vec3 qvec; VecLib::cross_avx(tvec, T.p1p2, qvec);
        

        m256_vec2 uv;
        VecLib::dot_mul_avx(tvec, pvec, invDet, uv.X);
        VecLib::dot_mul_avx(slope, qvec, invDet, uv.Y);
        VecLib::dot_mul_avx(T.p1p3, qvec, invDet, output.X); //t

        __m256 u_pass = _mm256_cmp_ps(uv.X, AVX_ZEROS, _CMP_GE_OQ);
        __m256 v_pass = _mm256_cmp_ps(uv.Y, AVX_ZEROS, _CMP_GE_OQ);
        __m256 uv_pass = _mm256_cmp_ps(
            _mm256_add_ps(
                uv.X,
                uv.Y
            ),
            _mm256_set1_ps(1),
            _CMP_LE_OQ
        );

        output.X = _mm256_and_ps(
            uv_pass,
            output.X
        );

        output.X = _mm256_and_ps(
            u_pass,
            output.X
        );

        output.X = _mm256_and_ps(
            v_pass,
            output.X
        );

        output.Y = uv.X;
        output.Z = uv.Y;
    }
    //static void get_normal(const PTri_AVX& T, const m256_vec3& position, m256_vec3& output) {
    //    __m256 dot;
    //    m256_vec3 relative_position;
    //    m256_vec3::sub(position, T.origin_offset, relative_position);
    //    VecLib::dot_avx(T.normal, relative_position, dot); //this can probably be done faster via xor of sign bit rather than multiplication, but whatever
    //    __m256 sgn = VecLib::sgn_fast(dot);
    //    output.X = _mm256_mul_ps(
    //        T.normal.X,
    //        sgn
    //    );
    //    output.Y = _mm256_mul_ps(
    //        T.normal.Y,
    //        sgn
    //    );
    //    output.Z = _mm256_mul_ps(
    //        T.normal.Z,
    //        sgn
    //    );
    //}
};

class Tri :public Primitive {
public:
    XYZ p1;
    XYZ p2;
    XYZ p3;
    XY UV_1;
    XY UV_2;
    XY UV_3;
    XYZ n1;
    XYZ n2;
    XYZ n3;
    XYZ midpoint;
    XYZ AABB_max;
    XYZ AABB_min;
    Tri(XYZ _p1, XYZ _p2, XYZ _p3, XY _UV_1, XY _UV_2, XY _UV_3, Material* _material) : Primitive(_material),
        p1(_p1), p2(_p2), p3(_p3) {
        midpoint = (p1 + p2 + p3) / 3;
        AABB_max = XYZ::max(p1, XYZ::max(p2, p3));
        AABB_min = XYZ::min(p1, XYZ::min(p2, p3));
        UV_1 = _UV_1;
        UV_2 = _UV_2;
        UV_3 = _UV_3;
        XYZ p1p3 = _p3 - p1;
        XYZ p1p2 = _p2 - p1;
        XYZ normal = XYZ::normalize(XYZ::cross(p1p3, p1p2));
        n1 = normal;
        n2 = normal;
        n3 = normal;
    }
    Tri(XYZ _p1, XYZ _p2, XYZ _p3, Material* _material) : Primitive(_material),
        p1(_p1), p2(_p2), p3(_p3) {
        midpoint = (p1 + p2 + p3) / 3;
        AABB_max = XYZ::max(p1, XYZ::max(p2, p3));
        AABB_min = XYZ::min(p1, XYZ::min(p2, p3));
    }
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
    static XYZ intersection_check_MoTru(const PackagedTri& T, const XYZ& position, const XYZ& slope) {
        XYZ pvec = XYZ::cross(slope, T.p1p3);
        decimal det = XYZ::dot(T.p1p2, pvec);
#if BACKFACE_CULLING
        // if the determinant is negative, the triangle is 'back facing'
        // if the determinant is close to 0, the ray misses the triangle
        (det < kEpsilon) return XYZ(-1);
#else
        // ray and triangle are parallel if det is close to 0
        if (fabs(det) < kEpsilon) return -1;
#endif
        decimal invDet = 1 / det;

        XYZ tvec = position - T.p1;
        decimal u = XYZ::dot(tvec, pvec) * invDet;
        if (u < 0 || u > 1) return XYZ(-1);

        XYZ qvec = XYZ::cross(tvec, T.p1p2);
        decimal v = XYZ::dot(slope, qvec) * invDet;
        if (v < 0 || u + v > 1) return XYZ(-1);

        decimal t = XYZ::dot(T.p1p3, qvec) * invDet;
        decimal U = u * T.U_delta.X + v * T.V_delta.X + T.UV_basis.X;
        decimal V = u * T.U_delta.Y + v * T.V_delta.Y + T.UV_basis.Y;
        return XYZ(t, U, V);
    }
    static XYZ intersection_check(const PackagedTri& T, const XYZ& position, const XYZ& slope) {
        return intersection_check_MoTru(T, position, slope);
    }
    static XYZ random(const PackagedTri& T) {
        float r1 = gen.fRand(0, 1);
        float r2 = gen.fRand(0, 1);
        if (r1 + r2 > 1) {
            r1 = 1 - r1;
            r2 = 1 - r2;
        }
        return T.p1 + r1 * T.p1p2 + r2 * T.p1p3;
    }
    static XYZ get_normal(XYZ normal, XYZ relative_position) {
        return (XYZ::dot(normal, relative_position) > 0) ? normal : -normal;
    }
    PackagedTri pack() {
        PackagedTri PT = PackagedTri(
            p1, p2, p3, UV_1, UV_2, UV_3, material
        );
        PT.n1 = n1;
        PT.n2 = n2;
        PT.n3 = n3;
        PT.calculate_tangents();
        PT.calculate_bitangents();
        PT.reorthogonalize();
        return PT;
    }
};



class Sphere : public Primitive {
public:
    decimal radius = 0;
    Sphere(decimal _radius, XYZ _origin, Material* _material) : Primitive(_material) {
        radius = _radius;
        origin = _origin;
        //bounds = get_bounds();
        obj_type = 1;
    }
    pair<XYZ, XYZ> get_bounds() {
        XYZ pos = origin + radius;
        XYZ neg = origin - radius;
        return pair<XYZ, XYZ>(pos, neg);
    }

    bool check_backface(XYZ& position) {
        return true;
    }

    decimal distance(XYZ& position) {
        return XYZ::distance(origin, position) - radius;
    }

    static decimal distance(const XYZ& origin, const decimal radius, const XYZ& position) {
        return XYZ::distance(origin, position) - radius;
    }

    static XYZ normal(const XYZ& origin, XYZ point) {
        return XYZ::slope(origin, point);
    }

    XYZ normal(XYZ point) {
        return XYZ::normalize(XYZ::slope(origin, point));
    }
    static decimal intersection_check(const XYZ& origin, const decimal radius, const XYZ& position, const XYZ& slope) {
        XYZ j = origin - position;
        decimal j_2 = XYZ::dot(j, j);
        decimal slope_2 = XYZ::dot(slope, slope);
        decimal k = XYZ::dot(j, slope);
        decimal sqrt_comp = k * k - slope_2 * (j_2 - radius * radius);
        if (sqrt_comp < 0) {
            return -1;
        }
        decimal g = sqrt(sqrt_comp);
        decimal d_1 = (k - g) / slope_2;
        decimal d_2 = (k + g) / slope_2;
        if (d_1 < 0) {
            return d_2;
        }
        if (d_2 < 0) {
            return d_1;
        }
        if (d_1 < d_2) {
            return d_1;
        }
        else {
            return d_2;
        }

    }
};

struct PackagedSphere { //reduced memory footprint to improve cache perf
    decimal radius;
    XYZ origin;
    Material* material;
    PackagedSphere(Sphere& sphere) {
        radius = sphere.radius;
        origin = sphere.origin;
        material = sphere.material;
    }
};

class Plane : public Primitive {
public:
    XYZ normal;
    XYZ origin_offset;
    Plane(XYZ _normal, XYZ _origin, Material* _material) : Primitive(_material) {
        origin = _origin;
        normal = XYZ::normalize(_normal);
        origin_offset = XYZ::dot(origin, normal) * normal;
    }
    void set_origin(XYZ _origin) {
        origin = _origin;
        origin_offset = XYZ::dot(origin, normal) * normal;
    }
    void set_normal(XYZ _normal) {
        normal = XYZ::normalize(_normal);
        origin_offset = XYZ::dot(origin, normal) * normal;
    }
    void recalc() {
        normal = XYZ::normalize(normal);
        origin_offset = XYZ::dot(origin, normal) * normal;
    }
    bool check_backface(XYZ& position) {
        return true;
    }
    static decimal distance(const XYZ& normal, const XYZ& origin_offset, const XYZ& position) {
        return XYZ::dot(position - origin_offset, normal);
    }
    static XYZ project_point(XYZ& normal, XYZ& position) {
        return (position - XYZ::dot(position, normal) * normal);
    }
    decimal distance(XYZ& position) {
        //XYZ projected_point = (position - XYZ::dot(position, normal) * normal)+origin_offset;
        //return XYZ::distance(position, projected_point);
        return XYZ::dot(position - origin_offset, normal);//this is like wildly simplified from position-(position after projection onto plane)
        //just ends in a dot product. Aint that weird
    }
    static decimal intersection_check(const XYZ& normal, const XYZ& origin_offset, const XYZ& position, const XYZ& slope) {
        decimal denom = XYZ::dot(normal, slope);
        if (denom == 0) {
            return -1;
        }
        return (XYZ::dot(origin_offset, normal) - XYZ::dot(position, normal)) / denom;
    }
};

struct PackagedPlane {
    XYZ normal;
    XYZ origin_offset;
    Material* material;
    PackagedPlane(Plane& plane) {
        normal = plane.normal;
        origin_offset = plane.origin_offset;
        material = plane.material;
    }
};

class Transformation {
private:
    Quat* rot_transform = nullptr;
    XYZ* XYZ_transform = nullptr;
public:

    Transformation(XYZ transformation, Quat rotation) {
        rot_transform = new Quat(rotation);
        XYZ_transform = new XYZ(transformation);
    }
    Transformation(Quat rotation) {
        rot_transform = new Quat(rotation);
    }
    Transformation(XYZ transformation) {
        XYZ_transform = new XYZ(transformation);
    }
    Transformation rotation() {
        return Transformation(*rot_transform);
    }
    void stack(XYZ& XYZ_transf, Quat& rotation) {
        if (XYZ_transform != nullptr) {
            XYZ_transf += *XYZ_transform;
        }
        if (rot_transform != nullptr) {
            rotation = Quat::multiply(rotation, *rot_transform);
            XYZ_transf = Quat::applyRotation(XYZ_transf, *rot_transform);
        }
    }
    //takes an XYZ position reference, and applies itself to it
    void apply(XYZ& position) {
        if (rot_transform != nullptr) {
            position = Quat::applyRotation(position, *rot_transform);
        }
        if (XYZ_transform != nullptr) {
            position += *XYZ_transform;
        }
    }
    static Transformation collpase(vector<Transformation>& transforms) {
        XYZ position = XYZ(0, 0, 0);
        Quat rotation = Quat(0, 0, 0, 1);
        for (Transformation& t : transforms) {
            t.stack(position, rotation);
        }
        return Transformation(position, rotation);
    }
};

namespace Packers { //I wanted to put these methods in the primitives objects, but due to me not using header files atm, its not feasible
    PackagedTri transformedPack(Tri& t, Transformation T, XYZ origin, XYZ scale) {
        XYZ np1 = t.p1 - origin;
        XYZ np2 = t.p2 - origin;
        XYZ np3 = t.p3 - origin;
        np1 *= scale;
        np2 *= scale;
        np3 *= scale;
        T.apply(np1);
        T.apply(np2);
        T.apply(np3);
        XYZ _n1 = t.n1;
        XYZ _n2 = t.n2;
        XYZ _n3 = t.n3;
        T.rotation().apply(_n1);
        T.rotation().apply(_n2);
        T.rotation().apply(_n3);
        np1 += origin;
        np2 += origin;
        np3 += origin;
        PackagedTri PT = PackagedTri(
            np1, np2, np3, t.UV_1, t.UV_2, t.UV_3, t.material
        );
        PT.n1 = _n1;
        PT.n2 = _n2;
        PT.n3 = _n3;
        PT.calculate_tangents();
        PT.calculate_bitangents();
        PT.reorthogonalize();
        return PT;
    }
    Tri transformT(Tri& t, Transformation T, XYZ origin, XYZ scale) {
        XYZ np1 = t.p1 - origin;
        XYZ np2 = t.p2 - origin;
        XYZ np3 = t.p3 - origin;
        np1 *= scale;
        np2 *= scale;
        np3 *= scale;
        T.apply(np1);
        T.apply(np2);
        T.apply(np3);
        XYZ _n1 = t.n1;
        XYZ _n2 = t.n2;
        XYZ _n3 = t.n3;
        T.rotation().apply(_n1);
        T.rotation().apply(_n2);
        T.rotation().apply(_n3);
        np1 += origin;
        np2 += origin;
        np3 += origin;
        Tri Tr = Tri(
            np1, np2, np3, t.UV_1, t.UV_2, t.UV_3, t.material
        );
        Tr.n1 = _n1;
        Tr.n2 = _n2;
        Tr.n3 = _n3;
        return Tr;
    }
}

class Mesh {
public:
    vector<Tri> tris;
    string name = "";
    int primitive_count;
    Mesh() {}
    void prep() {
        for (Tri& tri : tris) {
            tri.material->prep();
        }
    }
    void addTri(Tri& T) {
        tris.push_back(T);
        primitive_count++;
    }
};

class Object {
public:
    Object* parent = nullptr;

    XYZ origin;
    XYZ scale;
    vector<Transformation> transformations;

    vector<Object*> children;

    vector<Mesh*> meshes;
    vector<Sphere*> spheres;
    vector<Plane*> planes;

    string name = "";

    Object(XYZ _origin, XYZ _scale) :
        origin(_origin), scale(_scale) {

    }
    void addMesh(Mesh* M) {
        meshes.push_back(M);
    }
    void addSphere(Sphere* S) {
        spheres.push_back(S);
    }
    void addPlane(Plane* P) {
        planes.push_back(P);
    }
    void addChild(Object* object) {
        children.push_back(object);
        object->parent = this;
    }
    void prep() {
        for (Mesh* m : meshes) {
            m->prep();
        }
        for (Plane* p : planes) {
            p->material->prep();
        }
        for (Sphere* s : spheres) {
            s->material->prep();
        }
    }
    void fetchData(vector<Sphere>& spheres_vec, vector<Plane>& planes_vec, vector<Tri>& tris_vec) {
        XYZ scale_final = compile_scale();
        Transformation T = final_transform();
        for (Sphere* S : spheres) {
            Sphere new_sphere = Sphere(*S);
            T.apply(new_sphere.origin);
            new_sphere.origin *= scale_final;
            spheres_vec.push_back(new_sphere);
        }
        for (Plane* P : planes) {
            Plane new_plane = Plane(*P);
            T.apply(new_plane.origin);
            new_plane.origin *= scale_final;
            new_plane.normal *= scale_final;
            new_plane.recalc();
            planes_vec.push_back(new_plane);
        }
        for (Mesh* M : meshes) {
            for (Tri& tri : M->tris) {
                tris_vec.push_back(Packers::transformT(tri, T, origin, scale_final));
            }
        }
        for (Object* child : children) {
            child->fetchData(spheres_vec, planes_vec, tris_vec);
        }
    }
    void _registerRotation(Quat rot) {
        transformations.push_back(Transformation(rot));
    }
    void _registerMove(XYZ move) {
        transformations.push_back(Transformation(move));
    }
    void applyTransformXYZ(decimal x, decimal y, decimal z) {
        _registerMove(XYZ(x, y, z));
    }
    void applyTransformXYZ(XYZ xyz) {
        applyTransformXYZ(xyz.X, xyz.Y, xyz.Z);
    }
    void rotateX(decimal rotation) {
        XYZ orig = XYZ(0, 1, 0);
        XYZ pointing = XYZ(0, cos(rotation), sin(rotation));
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    void rotateY(decimal rotation) {
        XYZ orig = XYZ(0, 0, 1);
        XYZ pointing = XYZ(sin(rotation), 0, cos(rotation));
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    void rotateZ(decimal rotation) {
        XYZ orig = XYZ(1, 0, 0);
        XYZ pointing = XYZ(cos(rotation), sin(rotation), 0);
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    XYZ compile_scale(XYZ _scale = XYZ(1,1,1)) {
        _scale *= scale;
        if (parent != nullptr) {
           _scale*= parent->compile_scale(_scale);
        }
        return _scale;
    }
    void compile_transforms(vector<Transformation>& in) {
        
        for (auto& T : transformations) {
            in.push_back(T);
        }
        if (parent != nullptr) {
            parent->compile_transforms(in);
        }
        
    }
    Transformation final_transform() {
        vector<Transformation> total_transform;
        compile_transforms(total_transform);
        return Transformation::collpase(total_transform);
    }
};





//http://bannalia.blogspot.com/2015/06/cache-friendly-binary-search.html
int tree_total = 0;
int leaf_total = 0;
decimal shortest = 99;

class BVH {
public:
    XYZ max = XYZ(0, 0, 0);
    XYZ min = XYZ(0, 0, 0);
    vector<Tri*> elements;
    BVH* c1 = nullptr;
    BVH* c2 = nullptr;
    BVH* parent;
    int _count = -1;
    BVH() {    }
    BVH(vector<Tri*>* T_vec) {
        elements = vector<Tri*>(*T_vec);
        leaf_total += T_vec->size();
        for (Tri* T_ptr : *T_vec) {
            max = XYZ::max(max, T_ptr->AABB_max);
            min = XYZ::min(min, T_ptr->AABB_min);
        }

    }
    ~BVH() {
        delete c1;
        delete c2;
    }
    void printDesmos() {
        //cout << "B(" + min.to_string() + "," + max.to_string() + ")" << endl;
        if (c1 != nullptr) c1->printDesmos();
        if (c2 != nullptr) c2->printDesmos();
        for (auto& T : elements) {
            cout << "triangle(" + T->p1.to_string() + "," + T->p2.to_string() + "," + T->p3.to_string() + ")" << endl;
        }
    }
    static decimal intersection(const XYZ& max, const XYZ& min, const XYZ& origin, const XYZ& inv_slope) {
        decimal tx1 = (min.X - origin.X) * inv_slope.X;
        decimal tx2 = (max.X - origin.X) * inv_slope.X;

        decimal tmin = std::min(tx1, tx2);
        decimal tmax = std::max(tx1, tx2);

        decimal ty1 = (min.Y - origin.Y) * inv_slope.Y;
        decimal ty2 = (max.Y - origin.Y) * inv_slope.Y;

        tmin = std::max(tmin, std::min(ty1, ty2));
        tmax = std::min(tmax, std::max(ty1, ty2));

        decimal tz1 = (min.Z - origin.Z) * inv_slope.Z;
        decimal tz2 = (max.Z - origin.Z) * inv_slope.Z;

        tmin = std::max(tmin, std::min(tz1, tz2));
        tmax = std::min(tmax, std::max(tz1, tz2));

        if (tmax >= tmin) { //introduces a branch, but itll probably be fine. Lets me order search checks by nearest intersection
            if (tmin >= 0 && tmin < tmax) {
                return tmin;
            }
            else {
                return tmax;
            }
        }
        else {
            return -1;
        }
    }
    decimal intersection(const XYZ& origin, const XYZ& inv_slope) { //https://tavianator.com/2011/ray_box.html
        return BVH::intersection(max, min, origin, inv_slope);
    }

    struct WorkPacket {
        BVH* target;
        vector<Tri*>* contents;
        WorkPacket(BVH* t, vector<Tri*>* c) : target(t), contents(c) {}
    };
    //very dangerous. only ever use if BVH is being used as an intermediate
    //will turn every tri into a hanging pointer if the passed vector moves out of scope
    int construct_dangerous(vector<Tri>& initial_geo) {
        vector<Tri*>* working_data = new vector<Tri*>();
        working_data->reserve(initial_geo.size());
        for (Tri& T : initial_geo) {
            working_data->push_back(&T);
        }
        int return_value = construct(working_data);
        return return_value;
    }
    int construct(vector<Tri*>* initial_geo) {
        queue<WorkPacket> packets;
        assertm((initial_geo->size() > 0), "BVH was creation was attempted with zero geometry. Did everything load right?");
        auto initial_packet = WorkPacket(this, initial_geo);
        packets.push(initial_packet);
        int i = 0;
        while (packets.size() > 0) {
            i++;
            WorkPacket packet = packets.front();
            vector<Tri*>* geo = packet.contents;
            BVH* target = packet.target;

            assertm(geo->size() > 0, "zero size bin");

            AABB my_bounds = getBounds(geo);

            target->max = my_bounds.max;
            target->min = my_bounds.min;

            XYZ extent = target->max - target->min;
            int axis = 0;
            if (extent.Y > extent.X) axis = 1;
            if (extent.Z > extent[axis]) axis = 2;
            float splitPos = target->min[axis] + extent[axis] * 0.5f;

            Split* split = new Split(splitPos, axis);
            get_stats(geo, *split);
            evaluate_split(*split, 1);
            Split* probe_split = probe(geo);
            get_stats(geo, *probe_split);
            evaluate_split(*probe_split, 1);


            if (probe_split->score < split->score) {
                delete split;
                split = probe_split;
            }
            auto bins = bin(split, geo);
            //if (geo->size() > 200) {
            //    auto binned_split = binned_split_probe(geo);
            //    if (binned_split.first->score < split->score) {
            //        bins = binned_split.second;
            //        split = binned_split.first;
            //    }
            //}


            if (((size_t)abs((int)(bins.p_geo->size() - bins.n_geo->size()))) == geo->size()) {
                target->elements = vector<Tri*>(*geo);
                leaf_total += geo->size();
                packets.pop();
                continue;
            }

            if (bins.p_geo->size() <= LEAF_SIZE) {
                assertm(bins.p_geo->size() > 0, "attempted bin creation of zero elements");
                target->c1 = new BVH(bins.p_geo);
                target->c1->parent = target;
            }
            else {
                BVH* leaf = new BVH();
                target->c1 = leaf;
                target->c1->parent = target;
                WorkPacket k = WorkPacket(leaf, bins.p_geo);
                packets.push(k);
            }
            if (bins.n_geo->size() <= LEAF_SIZE) {
                assertm(bins.n_geo->size() > 0, "attempted bin creation of zero elements");
                target->c2 = new BVH(bins.n_geo);
                target->c2->parent = target;
            }
            else {
                BVH* leaf = new BVH();
                target->c2 = leaf;
                target->c2->parent = target;
                WorkPacket k = WorkPacket(leaf, bins.n_geo);
                packets.push(k);
            }

            delete geo;
            packets.pop();
        }
        return i;
    }
    pair<vector<BVH*>*, vector<Tri*>*> flatten(int depth = 0) {
        vector<BVH*>* l = new vector<BVH*>();
        vector<Tri*>* t = new vector<Tri*>();

        flatten(l, t, depth);
        return pair<vector<BVH*>*, vector<Tri*>*>(l, t);
    }
    void flatten(vector<BVH*>* l, vector<Tri*>* t, int depth) {
        depth = depth - 1;
        l->push_back(this);
        if (depth != 0) {
            if (c1 != nullptr) {
                c1->flatten(l, t, depth);
            }
            if (c2 != nullptr) {
                c2->flatten(l, t, depth);
            }
            if (elements.size() > 0) {
                for (auto T : elements) {
                    t->push_back(T);
                }
            }
        }
    }
    int size() {
        int s = 0;
        size(s);
        return s;
    }
    void size(int& s) {
        s++;
        if (c1 != nullptr) {
            c1->size(s);
            c2->size(s);
        }
    }
    int count() {
        if (_count == -1) {
            count(_count);
        }
        return _count;
    }
    void count(int& s) const {
        s += this->elements.size();
        if (c1 != nullptr) {
            c1->count(s);
            c2->count(s);
        }
    }
    float avg_path() const {
        if (c1 == nullptr) return 0;
        return (c1->avg_path() + c2->avg_path()) / 2.0;
    }
    vector<PackagedTri> get_emissive_tris() {
        vector<PackagedTri> out;
        get_emissive_tris(out);
        return out;
    }
    void get_emissive_tris(vector<PackagedTri>& vec) {
        if (c1 != nullptr) {
            c1->get_emissive_tris(vec);
            c2->get_emissive_tris(vec);
        }
        else {
            for (int i = 0; i < elements.size(); i++) {
                Tri* t = elements[i];
                if (t->material->emissive.getValue().magnitude() > 0) {
                    vec.push_back(t->pack());
                }
            }
        }
    }

private:
    struct BinResults {
        vector<Tri*>* p_geo;
        vector<Tri*>* n_geo;
    };
    struct Stats {
        XYZ min = XYZ(99999999999);
        XYZ max = XYZ();
        int count = 0;
        decimal SA() {
            return VecLib::surface_area(max, min);
        }
        decimal volume() {
            return VecLib::volume(max, min);
        }
    };
    struct Split {
        Stats p;
        Stats n;
        decimal score;
        XYZ placement;
        int facing;
        Split(XYZ place, int _facing) :placement(place), facing(_facing) {};
    };
    struct AABB {
        XYZ max;
        XYZ min;
        AABB() : max(-99999999999999), min(99999999999999) {}
        AABB(XYZ _max, XYZ _min) :max(_max), min(_min) {}
        void expand(XYZ point) {
            max = XYZ::max(max, point);
            min = XYZ::min(min, point);
        }
    };
    struct reduced_tri {
        XYZ aabb_max;
        XYZ aabb_min;
        XYZ midpoint;
        reduced_tri(Tri T) {
            aabb_max = T.AABB_max;
            aabb_min = T.AABB_min;
            midpoint = T.midpoint;
        }
        struct less_than_axis_operator {
            int axis;
            less_than_axis_operator(int _axis) {
                axis = _axis;
            }
            inline bool operator() (const reduced_tri t1, const reduced_tri t2)
            {
                return (t1.midpoint[axis] < t2.midpoint[axis]);
            }
        };
    };
    AABB getBounds(vector<Tri*>* geo) {
        AABB results;
        for (Tri* T_ptr : *geo) {
            results.expand(T_ptr->AABB_max);
            results.expand(T_ptr->AABB_min);
        }
        return results;
    }
    BinResults bin(Split* split, vector<Tri*>* geo) {
        BinResults results;
        results.p_geo = new vector<Tri*>();
        results.n_geo = new vector<Tri*>();
        int axis = split->facing;
        for (Tri* T_ptr : *geo) {
            if (T_ptr->midpoint[axis] > split->placement[axis]) {
                results.p_geo->push_back(T_ptr);
            }
            else {
                results.n_geo->push_back(T_ptr);
            }
        }
        return results;
    }
    struct relevant_value {
        XYZ midpoint;
        XYZ max;
        XYZ min;
        XYZ right_max;
        XYZ right_min;
        Tri* parent;
        int value_selector = 0;
        relevant_value(Tri* source, int v_selector = 0) {
            parent = source;
            midpoint = parent->midpoint;
            max = parent->AABB_max;
            min = parent->AABB_min;
            value_selector = v_selector;
        }
        XYZ get_value() const {
            switch (value_selector) {
            case 0:
                return min;
            case 1:
                return midpoint;
            case 2:
                return max;
            default:
                throw exception("out of bounds value selector");
            }
        }
        struct comparer {
        public:
            comparer(char compare_index, bool _use_selector = false) {
                internal_index = compare_index;
                use_selector = _use_selector;
            }
            inline bool operator()(const relevant_value& v1, const relevant_value& v2) {
                if (use_selector) {
                    XYZ v1_value = v1.get_value();
                    XYZ v2_value = v2.get_value();
                    return v1_value[internal_index] < v2_value[internal_index];
                }
                return v1.midpoint[internal_index] < v2.midpoint[internal_index];
            }
        private:
            char internal_index;
            bool use_selector;
        };
    };
    pair<Split*, BinResults> binned_split_probe(vector<Tri*>* geo) {
        int min_bins = 10;
        int max_bins = 1000;
        int bin_count = std::max(min_bins, std::min(max_bins, (int)geo->size()));
        int adaptive_bin_count = 3;
        int adpative_sweep_recusion_count = 1;

        int geo_count = geo->size();
        auto sorted = vector<relevant_value>();

        Split operator_split = Split(XYZ(), 0);
        Split best = operator_split;
        best.score = 999999999999999999999999.0;

        for (Tri* T_ptr : *geo) {
            sorted.push_back(relevant_value(T_ptr, 0));
            //sorted.push_back(relevant_value(T_ptr, 1));
            sorted.push_back(relevant_value(T_ptr, 2));
        }
        for (int i = 0; i < 3; i++) {
            operator_split.facing = i;
            sort(sorted.begin(), sorted.end(), relevant_value::comparer(i, true));

            XYZ left_min = sorted.front().min;
            XYZ left_max = sorted.front().min;
            XYZ right_min = sorted.back().min;
            XYZ right_max = sorted.back().max;
            int left_count = 0;
            int right_count = geo_count;



            for (int j = sorted.size() - 1; j >= 0; j--) {
                relevant_value& v = sorted[j];
                v.right_min = right_min;
                v.right_max = right_max;
                if (v.value_selector == 0) {
                    right_min = XYZ::min(right_min, v.min);
                    right_max = XYZ::max(right_max, v.max);
                }
            }

            double span = right_max[i] - right_min[i];
            double bin_interval = span / (bin_count);
            vector<double> bin_scores;

            operator_split.n.min = XYZ();
            operator_split.n.max = XYZ();
            operator_split.n.count = 0;
            operator_split.p.min = right_min;
            operator_split.p.max = right_max;
            operator_split.p.count = geo_count;
            evaluate_split(operator_split, 1);
            if (operator_split.score < best.score) best = operator_split;

            double pos = left_min[i];
            int index = 0;
            set<Tri*> partials;
            for (int j = 0; j < bin_count; j++) {
                while (index < sorted.size() && sorted[index].get_value()[i] < pos) {
                    relevant_value& v = sorted[index];
                    assert(v.parent != nullptr);
                    if (v.value_selector == 0) {
                        partials.insert(v.parent);
                        left_count++;
                    }
                    if (v.value_selector == 2) {
                        partials.erase(v.parent);
                        right_min = v.right_min;
                        right_max = v.right_max;
                        left_min = XYZ::min(v.min, left_min);
                        left_max = XYZ::max(v.max, left_max);
                        right_count--;
                    }
                    index++;
                }

                XYZ temp_right_min = right_min;
                XYZ temp_right_max = right_max;
                XYZ temp_left_min = left_min;
                XYZ temp_left_max = left_max;

                vector<XYZ> planar;
                for (Tri* T : partials) {
                    vector<XYZ> points;
                    vector<XYZ> left; //note, this is probably pretty inefficient, but its a simple solution
                    vector<XYZ> right;
                    points.push_back(T->p1);
                    points.push_back(T->p2);
                    points.push_back(T->p3);
                    for (int ti = 0; ti < 3; ti++) {
                        XYZ p = points[ti];
                        if (p[i] < pos)left.push_back(p);
                        if (p[i] == pos)planar.push_back(p);
                        if (p[i] > pos)right.push_back(p);
                    }

                    for (int li = 0; li < left.size(); li++) {
                        for (int ri = 0; ri < right.size(); ri++) {
                            XYZ v1 = left[li];
                            XYZ v2 = right[ri];
                            XYZ slope = XYZ::slope(v1, v2);
                            double t_span = v2[i] - v1[i];
                            double t = (pos - v1[i]) / t_span;
                            XYZ planar_point = slope * t + v1;
                            planar.push_back(planar_point);
                        }
                    }
                }
                for (XYZ& p : planar) {
                    temp_right_min = XYZ::min(p, temp_right_min);
                    temp_right_max = XYZ::max(p, temp_right_max);
                    temp_left_min = XYZ::min(p, temp_left_min);
                    temp_left_max = XYZ::max(p, temp_left_max);
                }
                operator_split.n.min = temp_left_min;
                operator_split.n.max = temp_left_max;
                operator_split.n.count = left_count;
                operator_split.p.min = temp_right_min;
                operator_split.p.max = temp_right_max;
                operator_split.p.count = right_count;
                operator_split.placement = temp_left_max;
                evaluate_split(operator_split, 1);
                if (operator_split.score < 0) {
                    assert(0);
                }
                if (operator_split.score < best.score) {
                    best = operator_split;
                }
                bin_scores.push_back(operator_split.score);
                pos += bin_interval;

            }
            for (int i = 0; i < bin_scores.size(); i++) {
                cout << bin_scores[i] << endl;
            }
        }
        exit(0);

        BinResults out_bin = BinResults();
        out_bin.n_geo = new vector<Tri*>();
        out_bin.p_geo = new vector<Tri*>();
        int facing = best.facing;
        for (Tri* T : *geo) {
            if (T->AABB_max[facing] < best.placement[facing]) {
                out_bin.n_geo->push_back(T);
                continue;
            }
            if (T->AABB_min[facing] > best.placement[facing]) {
                out_bin.p_geo->push_back(T);
                continue;
            }

            vector<XYZ> points;
            vector<XYZ> left; //note, this is probably pretty inefficient, but its a simple solution
            vector<XYZ> right;
            points.push_back(T->p1);
            points.push_back(T->p2);
            points.push_back(T->p3);
            Tri* left_T = new Tri(*T);
            Tri* right_T = new Tri(*T);
            left_T->AABB_min = T->AABB_max;
            left_T->AABB_max = T->AABB_min;
            right_T->AABB_min = T->AABB_max;
            right_T->AABB_max = T->AABB_min;
            double pos = best.placement[facing];
            for (int ti = 0; ti < 3; ti++) {
                XYZ p = points[ti];
                if (p[facing] < pos) {
                    left.push_back(p);
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                }
                if (p[facing] == pos) {
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
                if (p[facing] > pos) {
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
            }

            for (int li = 0; li < left.size(); li++) {
                for (int ri = 0; ri < right.size(); ri++) {
                    XYZ v1 = left[li];
                    XYZ v2 = right[ri];
                    XYZ slope = XYZ::slope(v1, v2);
                    double t_span = v2[facing] - v1[facing];
                    double t = (pos - v1[facing]) / t_span;
                    XYZ p = slope * t + v1;
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
            }
            out_bin.n_geo->push_back(left_T);
            out_bin.p_geo->push_back(right_T);
        }
        return pair<Split*, BinResults>(new Split(best), out_bin);
    }
    Split* probe(vector<Tri*>* geo) {

        auto sorted = vector<relevant_value>();
        for (Tri* T_ptr : *geo) {
            sorted.push_back(relevant_value(T_ptr));
        }
        Split operator_split = Split(XYZ(), 0);
        Split best = operator_split;
        best.score = 999999999999999999999999.0;
        int geo_count = geo->size();
        int iteration = 0;
        double last_read = 0;
        for (int i = 0; i < 3; i++) {
            operator_split.facing = i;
            sort(sorted.begin(), sorted.end(), relevant_value::comparer(i));

            XYZ left_min = sorted.front().min;
            XYZ left_max = sorted.front().max;
            XYZ right_min = sorted.back().min;
            XYZ right_max = sorted.back().max;
            int count = 0;
            for (int j = sorted.size() - 1; j >= 0; j--) {
                relevant_value& v = sorted[j];
                right_min = XYZ::min(right_min, v.min);
                right_max = XYZ::max(right_max, v.max);
                v.right_min = right_min;
                v.right_max = right_max;
            }

            operator_split.n.min = XYZ();
            operator_split.n.max = XYZ();
            operator_split.n.count = 0;
            operator_split.p.min = right_min;
            operator_split.p.max = right_max;
            operator_split.p.count = geo_count;
            evaluate_split(operator_split, 1);
            if (operator_split.score < best.score) best = operator_split;

            double prev_pos = -99999999;
            int segments = 100;
            double span = right_max[i] - right_min[i];
            double segment = segments / span;


            for (int j = 0; j < sorted.size(); j++) {
                auto v = sorted[j];
                left_min = XYZ::min(left_min, v.min);
                left_max = XYZ::max(left_max, v.max);
                count++;
                prev_pos = v.midpoint[i];
                operator_split.n.min = left_min;
                operator_split.n.max = left_max;
                operator_split.n.count = count;
                operator_split.p.min = v.right_min;
                operator_split.p.max = v.right_max;
                operator_split.p.count = geo_count - count;
                operator_split.placement = v.midpoint;
                evaluate_split(operator_split, 1);
                if (operator_split.score < best.score) {
                    best = operator_split;
                }
                if (abs(operator_split.score - last_read) / operator_split.score > 0.05) {
                    //cout << iteration << endl;
                    //cout << operator_split.score << endl;
                    last_read = operator_split.score;
                }
                iteration++;
            }
        }
        //exit(0);

        return new Split(best);
    }
    static void operate_stats(Stats& stat, const XYZ& point_max, const XYZ& point_min) {
        stat.count++;
        stat.max = XYZ::max(point_max, stat.max);
        stat.min = XYZ::min(point_min, stat.min);
    }
    static void get_stats(vector<Tri*>* geo, Split& split) {
        split.p.count = 0;
        split.n.count = 0;
        split.p.max = XYZ(-99999999);
        split.n.max = XYZ(-99999999);
        split.p.min = XYZ(9999999);
        split.n.min = XYZ(9999999);
        int axis = split.facing;
        for (const Tri* T_ptr : *geo) {
            if (T_ptr->midpoint[axis] > split.placement[axis]) {
                operate_stats(split.p, T_ptr->AABB_max, T_ptr->AABB_min);
            }
            if (T_ptr->midpoint[axis] <= split.placement[axis]) {
                operate_stats(split.n, T_ptr->AABB_max, T_ptr->AABB_min);
            }
        }
    }
    static decimal evaluate_split(Split& split, decimal SA_parent) {
        decimal traversal_cost = 1;
        decimal intersect_cost = 1;
        decimal Sa = split.p.SA() / SA_parent;
        decimal Sb = split.n.SA() / SA_parent;
        int count_a = split.p.count;
        int count_b = split.n.count;
#if PENALIZE_UNFILLED_LEAFS
        count_a = ceil(count_a / (float)LEAF_SIZE);
        count_b = ceil(count_b / (float)LEAF_SIZE);
        intersect_cost = 2;
#endif
        decimal Ha = Sa * count_a * intersect_cost;
        decimal Hb = Sb * count_b * intersect_cost;
        decimal final_h = traversal_cost + Ha + Hb;
        split.score = final_h;
        return final_h;
        //return (split.p.count+split.n.count)/geo->size()*(split.p.count+split.n.count);
    }

};



class OctBVH {
public:
    OctBVH* children[8];
    XYZ max;
    XYZ min;
    vector<Tri*> elements;
    OctBVH* parent;

    struct bounding_packet {
        Tri** tris;
        XYZ max;
        XYZ min;
        void reserve(unsigned int count) {
            tris = (Tri**) _aligned_malloc(sizeof(Tri*) * count, 64);
        }
        ~bounding_packet() {
            free(tris);
        }
    };
    void construct(vector<Tri>* tris) {
        bounding_packet* entry_packet = new bounding_packet();
        entry_packet->reserve(tris->size());
        for (int i = 0; i < tris->size(); i++) {
            entry_packet->tris[i] = &(*tris)[i];
        }
        recurse_construct(entry_packet, this);
    }
private:
    void recurse_construct(bounding_packet* BP, OctBVH* parent) {

    }
};

class PackagedBVH { //memory optimized. produced once the BVH tree is finalized
public:
    int32_t index;
    XYZ sMax;
    XYZ sMin;
    int32_t leaf_size;
    PackagedBVH(BVH* target) {
        sMax = target->max;
        sMin = target->min;
        index = -1;
    }
    static pair<PackagedBVH*, vector<PackagedTri>*> collapse(BVH* top) {
        vector<PackagedTri>* tri_vec = new vector<PackagedTri>();
        int tri_count = top->count();
        int tree_size = top->size();
        PackagedBVH* out = (PackagedBVH*)_aligned_malloc((tree_size + 1) * sizeof(PackagedBVH), 32);
        tri_vec->reserve(tri_count + 1);
        PackagedBVH::collapse(out, tri_vec, top);
        return pair<PackagedBVH*, vector<PackagedTri>*>(out, tri_vec);
    }
private:
    static void collapse(PackagedBVH* vec, vector<PackagedTri>* tri_vec, BVH* top) {
        vector<pair<BVH*, int>> current;
        vector<pair<BVH*, int>> next;

        int accumulated_index = 0;
        current.push_back(pair<BVH*, int>(top, -1));
        while (current.size() > 0) {
            for (int i = 0; i < current.size(); i++) {
                auto selection = current.at(i);
                BVH* target = selection.first;
                int parent_index = selection.second;
                auto PBVH = PackagedBVH(target);
                XYZ a = 0.5 * (PBVH.sMax + PBVH.sMin);
                //PBVH.sMax += 0.5*(PBVH.sMax-a);
                //PBVH.sMin += 0.5*(PBVH.sMin-a);
                PBVH.leaf_size = target->elements.size();
                if (PBVH.leaf_size > 0) {
                    PBVH.index = tri_vec->size();
                    for (int i = 0; i < PBVH.leaf_size; i++) {
                        tri_vec->push_back(target->elements[i]->pack());
                    }
                }
                if (target->c1 != nullptr) {
                    next.push_back(pair<BVH*, int>(target->c1, accumulated_index));
                    next.push_back(pair<BVH*, int>(target->c2, accumulated_index));
                }
                vec[accumulated_index] = PBVH;
                if (parent_index != -1) {
                    if (vec[parent_index].index == -1) {
                        vec[parent_index].index = accumulated_index;
                    }

                }
                else accumulated_index++;
                accumulated_index++;
            }
            current = next;
            next.clear();
        }
    }
};



int total_tris = 0;
int nodes_traversed = 0;

#define IOR_STACK_SIZE 5

struct PackagedRay {
    XYZ position;
    XYZ slope;
    char generation = 0;
    XYZ* output;
    XYZ coefficient;
    Material* IOR_stack[IOR_STACK_SIZE];
    PackagedRay() {};
    PackagedRay(XYZ _position, XYZ _slope, char gen) :
        position(_position),
        slope(_slope),
        generation(gen)
    {}
    void move(decimal distance) {
        position += slope * distance;
    }
};

struct MinimumRay {
    XYZ position;
    XYZ inv_slope;
    uint64_t index;
};

class BVH_AVX {
public:
    m256_vec3 max;
    m256_vec3 min;
    unsigned int leaf_size[8];
    unsigned int indexes[8];
    int parent_index;
    int self_index;
    vector<MinimumRay>* pool;
    BVH_AVX(vector<BVH*>& batch) {
        vector<XYZ> max_vec;
        vector<XYZ> min_vec;
        for (int i = 0; i < batch.size(); i++) {
            max_vec.push_back(batch[i]->max);
            min_vec.push_back(batch[i]->min);
            leaf_size[i] = batch[i]->elements.size();
            total_tris += leaf_size[i];
        }
        for (int i = batch.size(); i < 8;i++) {
            max_vec.push_back(XYZ());
            min_vec.push_back(XYZ());
            leaf_size[i] = 0;
            indexes[i] = 0;
        }
        max = m256_vec3(max_vec);
        min = m256_vec3(min_vec);

        pool = new vector<MinimumRay>();

    }
    __m256 intersectionv1(const m256_vec3& fusedorigin, const m256_vec3& inv_slope) const { //holy mother of moving data. I pray for you, my cpu, I pray
        __m256 tx1 = _mm256_fmadd_ps(min.X, inv_slope.X, fusedorigin.X);
        __m256 tx2 = _mm256_fmadd_ps(max.X, inv_slope.X, fusedorigin.X);
        __m256 ty1 = _mm256_fmadd_ps(min.Y, inv_slope.Y, fusedorigin.Y);
        __m256 ty2 = _mm256_fmadd_ps(max.Y, inv_slope.Y, fusedorigin.Y);
        __m256 tz1 = _mm256_fmadd_ps(min.Z, inv_slope.Z, fusedorigin.Z);
        __m256 tz2 = _mm256_fmadd_ps(max.Z, inv_slope.Z, fusedorigin.Z);

        __m256 tmin = _mm256_max_ps(
            _mm256_min_ps(tz1, tz2),
            _mm256_max_ps(
                _mm256_min_ps(ty1, ty2),
                _mm256_min_ps(tx1, tx2)
            ));
        __m256 tmax = _mm256_min_ps(
            _mm256_max_ps(tz1, tz2),
            _mm256_min_ps(
                _mm256_max_ps(ty1, ty2),
                _mm256_max_ps(tx1, tx2)
            ));
        __m256 diff = _mm256_sub_ps(tmax, tmin);
        __m256 masked_diff = _mm256_max_ps(
            diff,
            AVX_ZEROS
        );
        __m256 return_mask = _mm256_cmp_ps(
            masked_diff,
            diff,
            _CMP_EQ_OQ
        );
        __m256 final = _mm256_and_ps(
            _mm256_add_ps(
                tmin,
                masked_diff
            ),
            return_mask
        );
        //__m256 diff_is_pos = _mm256_cmp_ps(diff, AVX_ZEROS, _CMP_GT_OQ);

        return final;
    }
    __m256 intersection(const m256_vec3& fusedorigin, const m256_vec3& inv_slope) const { //holy mother of moving data. I pray for you, my cpu, I pray
        __m256 tx1 = _mm256_fmadd_ps(min.X, inv_slope.X, fusedorigin.X);
        __m256 tx2 = _mm256_fmadd_ps(max.X, inv_slope.X, fusedorigin.X);
        __m256 ty1 = _mm256_fmadd_ps(min.Y, inv_slope.Y, fusedorigin.Y);
        __m256 ty2 = _mm256_fmadd_ps(max.Y, inv_slope.Y, fusedorigin.Y);
        __m256 tz1 = _mm256_fmadd_ps(min.Z, inv_slope.Z, fusedorigin.Z);
        __m256 tz2 = _mm256_fmadd_ps(max.Z, inv_slope.Z, fusedorigin.Z);

        __m256 tminx = _mm256_min_ps(tx1, tx2);
        __m256 tminy = _mm256_min_ps(ty1, ty2);
        __m256 tminz = _mm256_min_ps(tz1, tz2);
        __m256 tmaxx = _mm256_max_ps(tx1, tx2);
        __m256 tmaxy = _mm256_max_ps(ty1, ty2);
        __m256 tmaxz = _mm256_max_ps(tz1, tz2);

        __m256 tmin = _mm256_max_ps(
            tminz,
            _mm256_max_ps(
                tminy,
                tminx
            ));
        __m256 tmax = _mm256_min_ps(
            tmaxz,
            _mm256_min_ps(
                tmaxy,
                tmaxx
            ));
        __m256 return_mask = _mm256_and_ps(
            _mm256_cmp_ps(
                tmax,
                AVX_ZEROS,
                _CMP_GE_OQ
            ),
            _mm256_cmp_ps(
                tmin,
                tmax,
                _CMP_LT_OQ
            )
        );


        __m256 final = _mm256_and_ps(
            tmin,
            return_mask
        );
        //__m256 diff_is_pos = _mm256_cmp_ps(diff, AVX_ZEROS, _CMP_GT_OQ);

        return final;
    }
    static pair<vector<BVH_AVX>*, vector<PackagedTri>*> collapse(BVH* top) {
        vector<PackagedTri>* tri_vec = new vector<PackagedTri>();
        int tri_count = top->count();
        int tree_size = top->size();
        vector<BVH_AVX> v;
        v.reserve(top->size());
        tri_vec->reserve(tri_count + 1);
        BVH_AVX::collapse(v, tri_vec, top);
        vector<BVH_AVX>* out = new vector<BVH_AVX>();
        out->swap(v);//trimming overallocation of vector
        return pair<vector<BVH_AVX>*, vector<PackagedTri>*>(out, tri_vec);
    }
private:
    struct comparer {
        bool operator()(BVH* b1, BVH* b2) {
            //return b1->avg_path() < b2->avg_path();
            return VecLib::surface_area(b1->max, b1->min) * b1->count() > VecLib::surface_area(b2->max, b2->min) * b2->count();
        }
    };
    static void collapse(vector<BVH_AVX>& vec, vector<PackagedTri>* tri_vec, BVH* top) {
        vector<pair<BVH*, pair<int, int>>> current;
        vector<pair<BVH*, pair<int, int>>> next;

        int accumulated_index = 0;
        current.push_back(pair<BVH*, pair<int, int>>(top, pair<int, int>(-1, -1)));
        while (current.size() > 0) {
            for (int i = 0; i < current.size(); i++) {
                auto selection = current.at(i);
                vector<BVH*> batch;
                vector<BVH*> has_tris;
                batch.push_back(selection.first);
                while (true) {
                    while (batch.size() + has_tris.size() < 8) {
                        BVH* at_index = batch[0];
                        if (at_index->c1 != nullptr) {
                            batch.push_back(at_index->c1);
                            batch.push_back(at_index->c2);
                            nodes_traversed += 2;
                            batch.erase(batch.begin());
                        }
                        else {
                            has_tris.push_back(at_index);
                            batch.erase(batch.begin());
                            if (batch.size() == 0) {
                                break;
                            }
                        }
                        sort(batch.begin(), batch.end(), comparer());
                    }
                    if (has_tris.size() + batch.size() >= 8 || batch.size() == 0) {
                        if (has_tris.size() + batch.size() > 8) {
                            cout << "error occured in AVX BVH construction: over allocated" << endl;
                            throw exception();
                        }
                        for (int i = 0; i < has_tris.size(); i++) {
                            batch.push_back(has_tris[i]);
                        }
                        break;
                    }
                }
                BVH* target = selection.first;
                int parent_index = selection.second.first;
                int point_back_index = selection.second.second;
                BVH_AVX ABVH(batch);
                ABVH.parent_index = parent_index;
                ABVH.self_index = point_back_index;
                for (int i = 0; i < batch.size(); i++) {
                    if (ABVH.leaf_size[i] > 0) {
                        ABVH.indexes[i] = tri_vec->size();
                        for (int j = 0; j < ABVH.leaf_size[i]; j++) {
                            tri_vec->push_back(batch[i]->elements[j]->pack());
                        }
                    }
                    else {
                        next.push_back(pair<BVH*, pair<int, int>>(batch[i], pair<int, int>(accumulated_index, i)));
                    }
                }
                vec.push_back(ABVH);
                if (parent_index != -1) {
                    BVH_AVX& parent = vec[parent_index];
                    for (int k = 0; k < 8; k++) {
                        parent.indexes[point_back_index] = accumulated_index;
                    }
                }
                accumulated_index++;
            }
            current = next;
            next.clear();
        }
    }
};



class Lens {
public:
    int resolution_x;
    int resolution_y;
    int subdivision_size;
    void (*_prep)(Lens* self, int resolution_x, int resolution_y, int subdivision_size);
    Lens* (*_clone)(Lens* self);
    XYZ(*_at)(Lens* self, int p_x, int p_y, int sample_index);
    XYZ(*_random_at)(Lens* self, int p_x, int p_y);
    void prep(int resolution_x, int resolution_y, int subdivision_size) {
        _prep(this, resolution_x, resolution_y, subdivision_size);
    }
    XYZ at(int p_x, int p_y, int sample_index) {
        return _at(this, p_x, p_y, sample_index);
    }
    XYZ random_at(int p_x, int p_y) {
        return _random_at(this, p_x, p_y);
    }
    Lens* clone() {
        Lens* obj = _clone(this);
        obj->resolution_x = resolution_x;
        obj->resolution_y = resolution_y;
        obj->subdivision_size = subdivision_size;
        return obj;
    }

};

class RectLens : public Lens {
public:
    decimal width;
    decimal height;
    vector<XY> subdiv_offsets;
    vector<vector<XY>> outputs;
    float pixel_width;
    float pixel_height;
    RectLens(decimal _width, decimal _height) :
        width(_width), height(_height) {
        _random_at = _random_at_function;
        _at = _at_function;
        _prep = _prep_function;
        _clone = _clone_function;
    }
private:
    static void _prep_function(Lens* self, int resolution_x, int resolution_y, int subdivision_count) {
        RectLens* self_rect = (RectLens*)self;
        decimal width = self_rect->width;
        decimal height = self_rect->height;
        decimal partial_width = width / resolution_x;
        decimal partial_height = height / resolution_y;
        self_rect->pixel_width = partial_width;
        self_rect->pixel_height = partial_height;
        decimal offset_x = -width / 2;
        decimal offset_y = -height / 2;
        decimal subdiv_width_x = partial_width / subdivision_count;
        decimal subdiv_width_y = partial_height / subdivision_count;
        for (int y = 0; y < subdivision_count; y++) {
            for (int x = 0; x < subdivision_count; x++) {
                self_rect->subdiv_offsets.push_back(XY(
                    subdiv_width_x * (x + 0.5),
                    subdiv_width_y * (y + 0.5)
                ));
            }
        }
        for (int y = 0; y < resolution_y; y++) {
            vector<XY> row;
            for (int x = 0; x < resolution_x; x++) {
                decimal final_x = offset_x + partial_width * x;
                decimal final_y = offset_y + partial_height * y;
                row.push_back(XY(
                    final_x,
                    final_y
                ));
            }
            self_rect->outputs.push_back(row);
        }

    }
    static XYZ _at_function(Lens* self, int p_x, int p_y, int sample_index) {
        RectLens* self_rect = (RectLens*)self;
        XY pos = self_rect->outputs[p_y][p_x] + self_rect->subdiv_offsets[sample_index];
        return XYZ(pos.X, 0, pos.Y);
    }
    static XYZ _random_at_function(Lens* self, int p_x, int p_y) {
        RectLens* self_rect = (RectLens*)self;
        XY pos = self_rect->outputs[p_y][p_x] + XY(gen.fRand(0, self_rect->pixel_width), gen.fRand(0, self_rect->pixel_height));
        return XYZ(pos.X, 0, pos.Y);
    }
    static Lens* _clone_function(Lens* self) {
        RectLens* rect_self = (RectLens*)self;
        Lens* obj = new RectLens(rect_self->width, rect_self->height);
        RectLens* rect_obj = (RectLens*)obj;
        rect_obj->outputs = rect_self->outputs;
        rect_obj->subdiv_offsets = rect_self->subdiv_offsets;
        return obj;
    }
};

class Camera {
    //Ill note this is a fairly non-direct class. the original was more streamlined, but I decided to update it to remove gunk
    //and to allow progressive output updates during monte carlo rendering
public:
    XYZ position;
    Quat rotation = Quat(0, 0, 0, 1);

    Lens* lens;
    XYZ focal_position;

    int current_resolution_x = 0;
    int current_resolution_y = 0;

    Camera(XYZ _position, Lens* _lens, decimal focal_distance) :
        lens(_lens), position(_position), focal_position(XYZ(0, 0, -focal_distance)) {}

    Camera(XYZ _position, Lens* _lens, XYZ _focal_position) :
        lens(_lens), position(_position), focal_position(_focal_position) {}

    void prep(int resolution_x, int resolution_y, int samples) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;
        lens->prep(resolution_x, resolution_y, samples);
    }

    XYZ slope_at(int p_x, int p_y, int sample_index) {
        return Quat::applyRotation(XYZ::slope(XYZ(0, 0, 0), lens->at(p_x, p_y, sample_index) + focal_position), rotation);
    }

    XYZ random_slope_at(int p_x, int p_y) {
        return Quat::applyRotation(XYZ::slope(XYZ(0, 0, 0), lens->random_at(p_x, p_y) + focal_position), rotation);
    }

    Camera* clone() {
        Lens* new_lens = lens->clone();
        Camera* obj = new Camera(position, new_lens, focal_position);
        obj->current_resolution_x = current_resolution_x;
        obj->current_resolution_y = current_resolution_y;
        return obj;
    }

    ~Camera() {
        delete lens;
    }
};

struct Casting_Diagnostics {
    long reflections_cast = 0;
    long shadows_cast = 0;
    long diffuses_cast = 0;
    long long rays_processed = 0;
    chrono::nanoseconds duration;
};

class PackagedScene {
public:
    vector<PackagedSphere> sphere_data; //stores primitives in a more packed data format for better cache optimizations.
    vector<PackagedPlane> plane_data; //Ive made quite a few mistakes throughout this project but I think Im finally on track with this
    vector<PackagedTri> tri_data;
    vector<PTri_AVX> avx_tri_data;
    vector<Light*> lights;
    vector<PackagedTri> emissive_tris;
    PackagedBVH* flat_bvh = nullptr;
    vector<BVH_AVX> avx_bvh;
    short monte_carlo_generations = 2;
    short max_generations = 5;
    short monte_carlo_max = 64;
    float monte_carlo_modifier = 1.0 / 16;
};

class Scene {
public:

    int object_count = 0;
    int primitive_count = 0;

    vector<Object*> objects;
    vector<Camera*> cameras;
    Camera* camera;

    vector<Light*> lights;

    int current_resolution_x;
    int current_resolution_y;
    Scene() {}

    void register_object(Object* obj) {
        object_count++;
        for (Mesh* M : obj->meshes) {
            primitive_count += M->primitive_count;
        }
        objects.push_back(obj);
    }
    void register_camera(Camera* cam) {
        cameras.push_back(cam);
    }
    void register_light(Light* li) {
        lights.push_back(li);
    }
    void merge(Scene& mergee) {
        for (Object* O : mergee.objects) {
            objects.push_back(O);
        }
        for (Camera* cam : mergee.cameras) {
            cameras.push_back(cam);
        }
        for (Light* light : mergee.lights) {
            lights.push_back(light);
        }
    }

    void prep(int res_x, int res_y, int subdiv_count) {
        auto start = chrono::high_resolution_clock::now();
        cout << padString("[Scene] Prepping", ".", 100) << flush;
        for (Object* O : objects) {
            O->prep();
        }
        camera = cameras[0];
        for (Camera* camera : cameras) {
            camera->prep(res_x, res_y, subdiv_count);
        }
        auto end = chrono::high_resolution_clock::now();
        cout << "[Done][" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms]" << endl;
    }

    PackagedScene* package() {
        vector<Sphere> spheres;
        vector<Plane> planes;
        vector<Tri> tris;
        BVH* bvh = new BVH();


        PackagedScene* PS = new PackagedScene();

        auto data_group_start = chrono::high_resolution_clock::now();
        cout << padString("[Scene] Regrouping data", ".", 100) << flush;

        for (auto O : objects) {
            O->fetchData(spheres, planes, tris);
        }
        for (Sphere& S : spheres) {
            PS->sphere_data.push_back(PackagedSphere(S));
        }
        for (Plane& P : planes) {
            PS->plane_data.push_back(PackagedPlane(P));
        }
        for (Light* L : lights) {
            PS->lights.push_back(L);
        }

        auto data_group_end = chrono::high_resolution_clock::now();
        cout << "[Done][" << chrono::duration_cast<chrono::milliseconds>(data_group_end - data_group_start).count() << "ms]" << endl;

        //for (int i = 0; i < tris.size();i++) {
        //    Tri T = tris[i];
        //    cout << "triangle(" << T.p1.to_string() << "," << T.p2.to_string() << "," << T.p3.to_string() << ")" << endl;
        //}

        if (tris.size() > 0) {
            auto BVH_con_start = chrono::high_resolution_clock::now();
            cout << padString("[BVH] Constructing", ".", 100) << flush;

            bvh->construct_dangerous(tris); //note, this will cause the BVH to break when tris moves out of scope

            auto BVH_con_end = chrono::high_resolution_clock::now();
            cout << "[Done][" << chrono::duration_cast<chrono::milliseconds>(BVH_con_end - BVH_con_start).count() << "ms]" << endl;
            auto BVH_collapse_start = chrono::high_resolution_clock::now();
            cout << padString("[BVH] Flattening", ".", 100) << flush;

            PS->emissive_tris = bvh->get_emissive_tris();
#if USE_AVX_BVH
            auto collapse_out = BVH_AVX::collapse(bvh);
            PS->avx_bvh = *collapse_out.first;
            PS->tri_data = *collapse_out.second;

#else
            auto collapse_out = PackagedBVH::collapse(bvh);
            PS->flat_bvh = collapse_out.first;
            PS->tri_data = *collapse_out.second;
#endif
            //bvh->printDesmos();
            delete bvh;

#if USE_AVX_TRI
            map<int, pair<BVH_AVX*, int>> order;
            for (int i = 0; i < PS->avx_bvh.size(); i++) {
                BVH_AVX& current = PS->avx_bvh[i];
                for (int i = 0; i < 8; i++) {
                    int leaf_size = current.leaf_size[i];
                    int index = current.indexes[i];
                    if (leaf_size > 0) {
                        order[index] = pair<BVH_AVX*, int>(&current, i);
                    }
                }
            }
            for (auto iter = order.begin(); iter != order.end(); ++iter)
            {
                int index = iter->first;
                auto datapoint = iter->second;
                BVH_AVX& current = *datapoint.first;
                int selection_index = datapoint.second;
                int leaf_size = current.leaf_size[selection_index];
                current.indexes[selection_index] = PS->avx_tri_data.size();
                current.leaf_size[selection_index] = ceil(((float)leaf_size) / 8.0);
                int i = 0;
                while (i < ceil(((float)leaf_size) / 8.0)) {
                    vector<PackagedTri> PTris;
                    for (int j = 0; j < 8 && i * 8 + j < leaf_size; j++) {
                        PTris.push_back(PS->tri_data[index + 8 * i + j]);
                    }
                    PS->avx_tri_data.push_back(PTri_AVX(PTris));
                    PS->avx_tri_data.back().index = index + 8 * i;
                    i++;
                }
            }
#endif
            auto BVH_collapse_end = chrono::high_resolution_clock::now();
            cout << "[Done][" << chrono::duration_cast<chrono::milliseconds>(BVH_collapse_end - BVH_collapse_start).count() << "ms]" << endl;
        }
        return PS;
    }
};
#define AVX_STACK_SIZE 256



struct AVX_AABB_ray {
    m256_vec3 fusedposition;
    m256_vec3 inv_slope;
    unsigned int index;
};

class RayEngine {
public:
    PackagedScene* data;
    vector<vector<PackagedRay>> ray_packet_queue;
    const int max_queue = pow(2, 20); //bit over a million
    const int sub_packet_size = 1024;
    RayEngine() {}
    void load_scene(PackagedScene* PS) {
        data = PS;
    }
    void traverse_packets() {
        vector<PackagedRay>& packet = ray_packet_queue.front();
        vector<MinimumRay> minimal_sub_packet;
        vector<AVX_AABB_ray> sub_packet;
        sub_packet.resize(sub_packet_size);
        minimal_sub_packet.resize(sub_packet_size);
        unsigned int packet_ray_index = 0;
        while (packet_ray_index < packet.size()) {
            for (unsigned int i = 0; i < sub_packet_size; i++) {
                auto& ray_data = packet[packet_ray_index];
                minimal_sub_packet[i] = MinimumRay{ ray_data.position,1 / ray_data.slope,packet_ray_index };
                sub_packet[i] = AVX_AABB_ray{ m256_vec3(-1 * ray_data.position / ray_data.slope),m256_vec3(1 / ray_data.slope),i };
                packet_ray_index++;
            }
        }


    }
    void iterativeBVH_AVX(CastResults& res, const XYZ& position, const XYZ& slope) {
        //global_counter++;
        m256_vec3 position_avx(position);
        m256_vec3 slope_avx(slope);
        m256_vec3 fusedposition(-1 * position / slope);
        m256_vec3 inv_slope(1 / slope);
        int stack[AVX_STACK_SIZE];
        float distances[AVX_STACK_SIZE];
        int tri_count[AVX_STACK_SIZE];
        int stack_index = 0;
        stack[0] = 0;
        distances[0] = 0;
        tri_count[0] = 0;
        while (stack_index >= 0) {
            const float dist = distances[stack_index];
            const int self_triangles = tri_count[stack_index];
            if (dist >= res.distance) {
                stack_index--;
                continue;
            };
            if (self_triangles == 0) {
                const BVH_AVX& current = data->avx_bvh[stack[stack_index]];
                stack_index--;
                auto m256_results = current.intersection(fusedposition, inv_slope);
                
                
                //float results[8];
                float* results = (float*)&m256_results;
                int sorted_index[8];
                int sorted_index_size = 0;
                for (int i = 0; i < 8; i++) {
                    float& num = results[i];
                    if (num != 0) {
                        int j = sorted_index_size - 1;
                        for (; j >= 0; j--) {
                            if (results[sorted_index[j]] > num) {
                                sorted_index[j + 1] = sorted_index[j];
                            }
                            else break;
                        }
                        sorted_index[j + 1] = i;
                        sorted_index_size++;
                    }
                }
                for (int i = sorted_index_size - 1; i >= 0; i--) {
                    int selection_index = sorted_index[i];
                    float dist2 = results[selection_index];
                    stack_index++;
                    stack[stack_index] = current.indexes[selection_index];
                    distances[stack_index] = dist2;
                    tri_count[stack_index] = current.leaf_size[selection_index];
                
                }

                //stack_index++; places the bvhs in without order. Slower overall, but the simpler method makes it still viable.
                //for (int i = 0; i < 8; i++) {
                //    bool sgnbit = results[i] != 0;
                //    stack[stack_index]      = current.indexes[i];
                //    tri_count[stack_index]  = current.leaf_size[i];
                //    distances[stack_index] = results[i];
                //    stack_index+=sgnbit;
                //}
                //stack_index--;

                //__m256 temp_m256 = simd_sort_1V(_mm256_max_ps(m256_results,AVX_ZEROS));
                //float* temp = (float*) & temp_m256;
                //
                //float values[8];
                //int pointer[8];
                //int size = 0;
                //for (int i = 0; i < 8; i++) {
                //    if (results[i] > 0) {
                //        values[size] = results[i];
                //        pointer[size] = i;
                //        size++;
                //    }
                //}
                //
                //for (int i = 7; i > 7-size; i--) {
                //    float dist2 = temp[i];
                //    int index = 0;
                //    for (int k = size-1; k >= 0; k--) {
                //        if (values[k] == dist2) {
                //            index = pointer[k];
                //            values[k] = -1;
                //            break;
                //        }
                //    }
                //    stack_index++;
                //    stack[stack_index] = current.indexes[index];
                //    distances[stack_index] = dist2;
                //    tri_count[stack_index] = current.leaf_size[index];
                //
                //}
            }
            else {
                int start_index = stack[stack_index];
                stack_index--;
                for (int i = 0; i < self_triangles; i++) {
#if !USE_AVX_TRI

                    const PackagedTri& t = (data->tri_data)[i + start_index];
                    XYZ distance = Tri::intersection_check_MoTru(t, position, slope);
                    if (distance.X >= 0) {
                        if (distance.X < res.distance) {
                            res.distance = distance.X;
                            res.normal = Tri::get_normal(t.normal, position  - t.origin_offset);
                            res.material = t.material;
                            res.UV = distance;
                        }
                    }
#else

                    const PTri_AVX& PTri_AVX_pack = (data->avx_tri_data)[i + start_index];
                    m256_vec3 intersection_stats;
                    //m256_vec3 normal_stats;
                    //vector<PackagedTri> pTris;
                    //vector<XYZ> real_res;
                    //vector<XYZ> real_normal;
                    //for (int l = 0; l < 8; l++) {
                    //    PackagedTri pTri = PackagedTri();
                    //    pTri.p1 = PTri_AVX_pack.p1.at(l);
                    //    pTri.p1p2 = PTri_AVX_pack.p1p2.at(l);
                    //    pTri.p1p3 = PTri_AVX_pack.p1p3.at(l);
                    //    pTri.normal = PTri_AVX_pack.normal.at(l);
                    //    pTri.origin_offset = PTri_AVX_pack.origin_offset.at(l);
                    //    pTri.UV_basis = PTri_AVX_pack.UV_basis.at(l);
                    //    pTri.U_delta = PTri_AVX_pack.U_delta.at(l);
                    //    pTri.V_delta = PTri_AVX_pack.V_delta.at(l);
                    //    pTri.material = PTri_AVX_pack.materials[l];
                    //    real_res.push_back(Tri::intersection_check(pTri, position, slope));
                    //    real_normal.push_back(Tri::get_normal(pTri.normal, position - pTri.origin_offset));
                    //    pTris.push_back(pTri);
                    //}
                    //
                    //PTri_AVX::intersection_check(PTri_AVX_pack, position_avx, slope_avx, intersection_stats);
                    //PTri_AVX::get_normal(PTri_AVX_pack, position_avx, normal_stats);
                    //
                    //float* dists = (float*)&intersection_stats.X;
                    //for (int i = 0; i < 8; i++) {
                    //    if (dists[i] > 0 && dists[i] < res.distance) {
                    //        res.distance = dists[i];
                    //        res.normal = normal_stats.at(i);
                    //        res.material = PTri_AVX_pack.materials[i];
                    //        res.UV = intersection_stats.at(i);
                    //        res.tri_index = PTri_AVX_pack.get_index(i);
                    //    }
                    //}
                    PTri_AVX::intersection_check(PTri_AVX_pack, position_avx, slope_avx, intersection_stats);

                    float* dists = (float*)&intersection_stats.X;
                    for (int i = 0; i < 8; i++) {
                        if (dists[i] > 0 && dists[i] < res.distance) {
                            res.distance = dists[i];
                            res.normal = Tri::get_normal(PTri_AVX_pack.normal[i], position - PTri_AVX_pack.origin_offset[i]);
                            res.material = PTri_AVX_pack.materials[i];
                            res.TUV = intersection_stats.at(i);
                            res.tri_index = PTri_AVX_pack.get_index(i);
                        }
                    }
#endif
                }

            }
        }
    }
    void iterativeBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        int stack[64];
        float dist[64];
        int stack_index = 0;
        stack[0] = 0;
        while (stack_index >= 0) {
            const PackagedBVH& current = data->flat_bvh[stack[stack_index]];
            stack_index--;
            if (dist[stack_index + 1] > res.distance) continue;
            //if (distance <= 0 || distance >= res.distance) {
            //    continue;
            //}
            if (current.leaf_size == 0) {
                int c1_index = current.index;
                int c2_index = c1_index + 1;
                const PackagedBVH& c1 = data->flat_bvh[c1_index];
                const PackagedBVH& c2 = data->flat_bvh[c2_index];
                decimal t1 = BVH::intersection(c1.sMax, c1.sMin, position, inv_slope);
                decimal t2 = BVH::intersection(c2.sMax, c2.sMin, position, inv_slope);

                bool f1 = t1 < 0;//t1 < 0 || t1>res.distance;
                bool f2 = t2 < 0;//(t2 < 0 || t2>res.distance);

                if (f1 && f2) {
                    continue;
                }
                if (f1) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                    continue;
                }
                if (f2) {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                    continue;
                }
                if (t1 > t2 && abs(t1-t2)>0.1) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                }
                else {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                }
            }
            else {
                for (int i = 0; i < current.leaf_size; i++) {
                    const PackagedTri& t = (data->tri_data)[i + current.index];
                    XYZ distance = Tri::intersection_check(t, position, slope);
                    if (distance.X >= 0) {
                        if (distance.X < res.distance) {
                            res.distance = distance.X;
                            res.normal = Tri::get_normal(t.normal, position  - t.origin_offset);
                            res.material = t.material;
                            res.TUV = distance;
                        }
                    }
                }
            }
        }
    }
    void navBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
#if USE_AVX_BVH
        //iterativeBVH_AVX_unordered(starting_index, res, position, slope);
        iterativeBVH_AVX(res, position, slope);
#else
        iterativeBVH(res, position, slope, inv_slope);
#endif
    }
    CastResults execute_bvh_cast(XYZ& position, const XYZ& slope) {
        CastResults returner;
        for (const PackagedSphere& s : data->sphere_data) {
            decimal distance = Sphere::intersection_check(s.origin, s.radius, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = Sphere::normal(s.origin, position + distance * slope);
                    returner.material = s.material;
                }
            }
        }
        for (const PackagedPlane& p : data->plane_data) {
            decimal distance = Plane::intersection_check(p.normal, p.origin_offset, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = p.normal;
                    returner.material = p.material;
                }
            }
        }
        XYZ inv_slope = XYZ(1) / slope;
        navBVH(returner, position, slope, inv_slope);
        position += returner.distance * slope * 0.999;
        return returner;
    }
    CastResults execute_naive_cast(XYZ& position, const XYZ& slope) {
#define default_smallest_distance 9999999;
        decimal smallest_distance = default_smallest_distance;
        CastResults returner = CastResults(XYZ(-1, -1, -1), nullptr);
        for (const PackagedSphere& s : data->sphere_data) {
            decimal distance = Sphere::intersection_check(s.origin, s.radius, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = Sphere::normal(s.origin, position + distance * slope);
                    returner.material = s.material;
                }
            }
        }
        for (const PackagedPlane& p : data->plane_data) {
            decimal distance = Plane::intersection_check(p.normal, p.origin_offset, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = p.normal;
                    returner.material = p.material;
                }
            }
        }
        for (const PackagedTri& t : data->tri_data) {
            XYZ distance = Tri::intersection_check(t, position, slope);
            if (distance.X >= 0) {
                if (distance.X < returner.distance) {
                    returner.distance = distance.X;
                    returner.normal = Tri::get_normal(t.normal, position - t.origin_offset);
                    returner.material = t.material;
                }
            }
        }
        position += returner.distance * slope;
        return returner;
    }
    CastResults execute_ray_cast(XYZ& position, const XYZ& slope) {
        return execute_bvh_cast(position, slope);
        //return execute_naive_cast(position,slope);

    }
    bool execute_lighting_cast(XYZ position, Light* L) {
        XYZ slope = L->vectorTo(position);
        auto res = execute_bvh_cast(position, slope);
        if (L->intersects(position, slope, res)) {
            return true;
        }
        return false;
    }
    XYZ process_ray(Casting_Diagnostics& stats, PackagedRay& ray_data) {
        stats.rays_processed++;

        XYZ orig = ray_data.position;
        CastResults results = execute_ray_cast(ray_data.position, ray_data.slope);
        PackagedTri& Tri = data->tri_data[results.tri_index];
        results.texUV.X = results.TUV.Y * Tri.U_delta.X + results.TUV.Z * Tri.V_delta.X + Tri.UV_basis.X;
        results.texUV.Y = results.TUV.Y * Tri.U_delta.Y + results.TUV.Z * Tri.V_delta.Y + Tri.UV_basis.Y;

        if (results.material == nullptr) {
            return XYZ(0);
        }
        MaterialSample material = results.material->sample_UV(results.texUV.X, results.texUV.Y);
        XYZ aggregate = XYZ();
        XYZ normal = Tri::get_normal(Tri.get_normal(results.TUV.Y,results.TUV.Z),orig-Tri.origin_offset);
        if (!XYZ::equals(XYZ(0, 0, 0), material.normalPertubation)) {
            XYZ T  = Tri.get_tangent(results.TUV.Y,results.TUV.Z);
            XYZ BT = Tri.get_bitangent(results.TUV.Y, results.TUV.Z);
            normal = material.normalPertubation.X * T + material.normalPertubation.Y * BT + material.normalPertubation.Z * normal;

        }
        normal = XYZ::normalize(normal);
        {
            if (DRAW_UV) {
                if (results.TUV.X != -1) {
                    XYZ color;
                    color.X = results.texUV.X;
                    color.Y = results.texUV.Y;
                    color.Z = 1 - color.X - color.Y;
                    return color * 1;
                }
            }
            if (DRAW_NORMAL) {
                XYZ res = (normal);
                return (normal / 2 + XYZ(0.5, 0.5, 0.5)) * 1;

            }
            if (DRAW_NORMAL_DELTAS) {
                XYZ res = (results.normal - normal);
                if (res.magnitude() < 0.1) {
                    return XYZ(0, 0, 0);
                }
                return (XYZ::normalize(res) / 2 + XYZ(0.5, 0.5, 0.5)) * 1;

            }
            if (DRAW_NORMAL_FACING) return XYZ::dot(normal, orig) / abs(XYZ::dot(normal, orig));
            if (DRAW_COLOR) return material.color;
            if (false) {
                if (ray_data.generation == 1)
                    aggregate += ray_data.position / 2 + XYZ(0.5, 0.5, 0.5);
            }
            //if (DRAW_EMISSIVE) return material.calculate_emissions();
        }

        aggregate += material.calculate_emissions();

        int bounce_count = data->monte_carlo_max * pow(data->monte_carlo_modifier, ray_data.generation);
        XYZ flipped_output = XYZ::flip(ray_data.slope);
        XYZ reflection_slope = XYZ::normalize(XYZ::reflect(XYZ::flip(ray_data.slope), normal));

        Quat normal_rot = Quat::makeRotationFromY(normal);
        Quat reflection_rot = Quat::makeRotationFromY(reflection_slope);
        Matrix3x3 normal_rot_m = Matrix3x3::quatToMatrix(normal_rot);
        Matrix3x3 reflection_rot_m = Matrix3x3::quatToMatrix(reflection_rot);
        //ray_data.PreLL.value = XYZ(0, 0, 100);

        float dot_NO = XYZ::dot(normal, flipped_output);
        material.config(flipped_output, normal, 1.0003, material.IOR);

#if DRAW_LIGHTS
        for (Light* L : data->lights) {
            XYZ pos = ray_data.position + 0.001 * results.normal;
            XYZ slope = L->vectorTo(pos);
            material.set_input(slope);
            XYZ return_coefficient = material.fast_BRDF_co(true);
            if (execute_lighting_cast(pos, L)) {
                aggregate += return_coefficient * L->emission;
#if DRAW_LIGHT_EXPOSURE
                return aggregate;
            }
            return XYZ();
#else
        }
#endif
        }
#endif
        if (ray_data.generation >= data->max_generations) return aggregate;
        if (ray_data.generation < data->monte_carlo_generations) {


            //int diffuse_bounces = bounce_count;
            //int specular_bounces = 0;
            int diffuse_bounces  = bounce_count * material.f_diff;
            int specular_bounces = bounce_count * material.f_spec;
#if !DRAW_DIFFUSE
            diffuse_bounces = 0;
            specular_bounces = bounce_count;
#endif
#if !DRAW_SPECULAR
            diffuse_bounces = bounce_count;
            specular_bounces = 0;
#endif
#if !DRAW_SPECULAR && !DRAW_DIFFUSE
            diffuse_bounces = 0;
            specular_bounces = 0;
#endif
            if (diffuse_bounces < 4) diffuse_bounces = 4;
            int diffuse_Vslices = floor(log2(diffuse_bounces) - 1);
            int diffuse_Rslices = floor(diffuse_bounces / diffuse_Vslices);
            diffuse_bounces = diffuse_Vslices * diffuse_Rslices;

            double diffuse_multiplier = 1.0 / bounce_count;
            double specular_multiplier = 1.0 / bounce_count;



            float diffuse_V_increment = 1.0 / diffuse_Vslices;
            float diffuse_R_increment = 1.0 / diffuse_Rslices;
            float diffuse_V_position = 0;
            float diffuse_R_position = 0;


            for (int i = 0; i < diffuse_bounces; i++) {
                if (i % diffuse_Rslices == 0 && i != 0) {
                    diffuse_V_position += diffuse_V_increment;
                }
                diffuse_R_position += diffuse_R_increment;
                if (diffuse_R_position > 1) diffuse_R_position -= 1;
                XYZ diffuse_slope = material.biased_diffuse_bounce(
                    gen,
                    normal_rot_m
                    , diffuse_V_position
                    , diffuse_V_position + diffuse_V_increment
                    , diffuse_R_position
                    , diffuse_R_position + diffuse_R_increment
                );
                //XYZ diffuse_slope = material.biased_diffuse_bounce(
                //    normal_rot_m
                //);
                if (DRAW_BOUNCE_DIRECTION) {
                    if (ray_data.generation == 0) {
                        aggregate += diffuse_multiplier * (diffuse_slope / 2 + XYZ(0.5, 0.5, 0.5)) * 1;
                        continue;
                    }
                }

                //XYZ return_coefficient = ray_data.coefficient * material.fast_BRDF_co(normal, diffuse_slope, flipped_output);
                material.set_input(diffuse_slope);
                XYZ return_coefficient = material.fast_BRDF_co(false);
#if DRAW_DIFFUSE_COEFFICIENTS
                aggregate += return_coefficient * diffuse_multiplier;
                continue;
#endif
                if (XYZ::magnitude(return_coefficient) < 0.0001) continue;

                return_coefficient = return_coefficient * diffuse_multiplier;

                auto ray = PackagedRay(
                    ray_data.position + 0.001 * results.normal,
                    diffuse_slope,
                    ray_data.generation + 1
                );
                ray.output = ray_data.output;
                stats.diffuses_cast++;
                XYZ returned_light = process_ray(stats, ray);
                aggregate += return_coefficient * returned_light;
            }
            for (int i = 0; i < specular_bounces;i++) {
                XYZ specular_slope;
                int tries = 0;
                do {
                    if (tries > 10) throw exception();
                    if (material.roughness > 0) {
                        specular_slope = material.reflective_bounce(gen, reflection_rot_m);
                    }
                    else {
                        specular_slope = reflection_slope;
                    }
                } while (false);
                material.set_input(specular_slope);
                XYZ return_coefficient = material.fast_BRDF_co(true);
                //XYZ return_coefficient = XYZ(1, 1, 1);
#if DRAW_REFLECTIVE_COEFFICIENTS
                aggregate += return_coefficient*specular_multiplier;
                continue;
#endif
#if DRAW_DOT_RETURNS
    
#endif
                //stringstream s;
                //float f = XYZ::dot(flipped_output, specular_slope);
                //s << left << setw(10) << ((f>0) ? "+" + to_string(f) : to_string(f)) << flipped_output.to_string() << " " << specular_slope.to_string() << " " << normal.to_string() << " " << endl;
                //cout << s.str();

                if (DRAW_BOUNCE_DIRECTION) {
                    if (ray_data.generation == 0) {
                        aggregate += specular_multiplier * (specular_slope / 2 + XYZ(0.5, 0.5, 0.5)) * 1;
                        continue;
                    }

                }

                if (XYZ::magnitude(return_coefficient) < 0.0001) continue;
                return_coefficient = return_coefficient * specular_multiplier;


                auto ray = PackagedRay(
                    ray_data.position + 0.001 * results.normal,
                    specular_slope,
                    ray_data.generation + 1
                );

                ray.output = ray_data.output;
                stats.diffuses_cast++;
                XYZ returned_light = process_ray(stats, ray);
                aggregate += return_coefficient * returned_light;

            }
        }
        else {
            material.set_input(reflection_slope);
            XYZ return_coefficient = material.fast_BRDF_co(true);
            
            if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                return aggregate;
            }
            auto ray = PackagedRay(
                ray_data.position + NEAR_THRESHOLD * results.normal * 10,
                reflection_slope,
                ray_data.generation + 1
            );
            ray.output = ray_data.output;
            aggregate += return_coefficient * process_ray(stats, ray);

            stats.reflections_cast++;
        }
        return aggregate;
    }
};



class GUIHandler {
public:
    int current_resolution_x = 0;
    int current_resolution_y = 0;

    sf::RenderWindow* window;
    sf::Image canvas;
    sf::Texture render_texture;
    sf::Sprite render_sprite;

    map<string, sf::Texture> textures;
    map<string, sf::Sprite> sprites;
    map<string, sf::Image> images;
    map<string, sf::RectangleShape> rects;

    double scalar_exponent = 0;

    XY focus_size;
    vector<XY> focuses;

    GUIHandler() {}

    GUIHandler(int resX, int resY, int block_size) {
        current_resolution_x = resX;
        current_resolution_y = resY;
        canvas.create(current_resolution_x, current_resolution_y, sf::Color::Black);
        render_sprite.setScale(make_scale(), -make_scale());
        render_sprite.setPosition(0, make_scale() * current_resolution_y);
        focus_size = XY(block_size, block_size);
    }

    void hold_window() {
        auto frame_time = chrono::high_resolution_clock::now();
        while (window->isOpen())
        {
            sf::Event event;
            while (window->pollEvent(event))
            {
                handle_events(event);
            }
            partial_render();
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }
    }
    XY get_relative_position(int x, int y) {
        double relative_x = ((double)x) / current_resolution_x;
        double relative_y = ((double)y) / current_resolution_y;
        return XY(make_scale() * relative_x, make_scale() * relative_y);

    }
    void create_window() {
        cout << "Opening Window...." << flush;
        window = new sf::RenderWindow(sf::VideoMode(current_resolution_x * make_scale(), current_resolution_y * make_scale()), "Render Window!");
        cout << "Done" << endl;
    }
    void commit_pixel(XYZ color, int x, int y) {
        auto sfcolor = sf::Color(color.X, color.Y, color.Z);
        canvas.setPixel(x, y, sfcolor);
    }
    void partial_render() {
        window->clear();
        draw_render_sprite();
        window->display();
    }
    void render_pass() {
        window->clear();
        update_focuses();
        draw_render_sprite();
        draw_secondary_sprites();
        draw_shapes();
        window->display();
    }
    void update_focuses() {
        XY scaled_size = focus_size * make_scale();
        for (int i = 0; i < focuses.size(); i++) {
            XY focus_position = focuses.at(i);
            focus_position.Y = current_resolution_y - focus_position.Y - focus_size.Y;
            XY scaled_position = focus_position * make_scale();
            string focus_id = "render_focus_" + to_string(i);
            rects[focus_id].setSize(sf::Vector2f(scaled_size.X, scaled_size.Y));
            rects[focus_id].setPosition(scaled_position.X, scaled_position.Y);
        }
    }
    void draw_render_sprite() {
        window->draw(render_sprite);
    }
    void draw_secondary_sprites() {
        for (auto i = sprites.begin(); i != sprites.end(); i++) {
            sf::Sprite& sprite = i->second;
            window->draw(sprite);

        }
    }
    void draw_shapes() {
        for (auto i = rects.begin(); i != rects.end(); i++) {
            sf::RectangleShape& rect = i->second;
            window->draw(rect);
        }
    }
    void update_image_textures() {
        for (auto i = images.begin(); i != images.end(); i++) {
            string key = i->first;
            sf::Image& image = i->second;
            getTexture(key).loadFromImage(image);
        }
    }
    void commit_canvas() {
        render_texture.loadFromImage(canvas);
        render_sprite.setTexture(render_texture, false);
        render_pass();
    }
    void add_texture_sprite(string name) {
        textures[name] = sf::Texture();
        sprites[name] = sf::Sprite();
        getSprite(name).setTexture(getTexture(name));
    }
    void add_texture_sprite_image(string name) {
        textures[name] = sf::Texture();
        sprites[name] = sf::Sprite();
        images[name] = sf::Image();
        getSprite(name).setTexture(getTexture(name));
    }
    sf::RectangleShape& add_rect(string name) {
        rects[name] = sf::RectangleShape();
        return rects[name];
    }
    void add_focus() {
        string focus_id = "render_focus_" + to_string(focuses.size());
        add_rect(focus_id);
        focuses.push_back(XY(0, 0));

        sf::RectangleShape& focus = rects[focus_id];
        focus.setFillColor(sf::Color::Transparent);
        focus.setOutlineThickness(1);
        focus.setOutlineColor(sf::Color(255, 211, 110));
    }
    void set_focus_position(int index, int x, int y) {
        focuses[index] = XY(x, y);
    }
    sf::Sprite& getSprite(string name) {
        return sprites[name];
    }
    sf::Texture& getTexture(string name) {
        return textures[name];
    }
    sf::Image& getImage(string name) {
        return images[name];
    }
    void handle_events() {
        sf::Event event;
        while (window->pollEvent(event)) {
            handle_events(event);
        }
    }
    double make_scale() {
        return pow(2, scalar_exponent);
    }
    void handle_events(sf::Event event) {
        if (event.type == sf::Event::Closed) window->close();
        if (event.type == sf::Event::Resized)
        {
            // update the view to the new size of the window
            sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
            window->setView(sf::View(visibleArea));
        }
        if (event.type == sf::Event::MouseWheelMoved)
        {
            // display number of ticks mouse wheel has moved
            double amount = event.mouseWheel.delta;
            scalar_exponent += amount / 5;
            double scale = make_scale();
            render_sprite.setScale(scale, -scale);
            render_sprite.setPosition(0, make_scale() * current_resolution_y);
        }
    }
};
namespace ImageHandler {
    OCIO::ConstConfigRcPtr config = OCIO::Config::CreateFromFile("C:\\Users\\Charlie\\Libraries\\OCIOConfigs\\AgX-main\\config.ocio");
    OCIO::ConstProcessorRcPtr processor = config->getProcessor("Linear Rec.709", "AGX Base sRGB");
    auto compute = processor->getDefaultCPUProcessor();
    static decimal luminance(XYZ v)
    {
        return XYZ::dot(v, XYZ(0.2126, 0.7152, 0.0722));
    }

    static XYZ change_luminance(XYZ c_in, decimal l_out)
    {
        decimal l_in = luminance(c_in);
        return c_in * (l_out / l_in);
    }

    //decimal realistic_response(decimal f, decimal iso) // https://graphics-programming.org/resources/tonemapping/index.html

    static XYZ apply_gamma(XYZ in, float gamma) {
        float lums = luminance(in);
        float post_lums = pow(lums, 1.0/gamma);
        return in * post_lums / lums;
    }

    XYZ post_process_pixel(XYZ lums) {  //lums is raw light data via RGB
        lums = apply_gamma(lums, 1);
        float pixel[3] = { lums[0],lums[1],lums[2] };
        compute->applyRGB(pixel);
        lums = XYZ(pixel[0], pixel[1], pixel[2]);
        lums = apply_gamma(lums, 2.2);

        return XYZ::clamp(lums, 0, 1) * 255;//scale to 24 bit rgb
    }
    static void postProcessRaw(vector<vector<XYZ*>>* data) {
        for (vector<XYZ*>& data_row : *data) {
            for (XYZ*& pixel_ptr : data_row) {
                *pixel_ptr = ImageHandler::post_process_pixel(*pixel_ptr);
            }
        }
    }
    static GUIHandler* openImageVector(vector<vector<XYZ*>>* data) {
        int resX = data->at(0).size();
        int resY = data->size();
        cout << resX << " " << resY << endl;
        GUIHandler* GUI = new GUIHandler(resX, resY, 0);
        for (int y = 0; y < resY;y++) {
            for (int x = 0; x < resX; x++) {
                GUI->commit_pixel(*(data->at(y)[x]), x, y);
            }
        }
        GUI->create_window();
        GUI->commit_canvas();
        return GUI;
    }
};

class RenderThread {
public:
    struct block {
        int size;
        int block_index_x;
        int block_index_y;
        iXY offset;
        vector<XYZ*> raws;
        atomic<int> pixels_done;
        atomic<int> iterations_done{ 0 };
        int pixels_last;
        atomic<bool> done = false;
        block(int _block_size, int _block_x, int _block_y) {
            size = _block_size;
            block_index_x = _block_x;
            block_index_y = _block_y;
            int x_offset = block_index_x * size;
            int y_offset = block_index_y * size;
            offset = iXY(x_offset, y_offset);
            pixels_done = 0;
            pixels_last = 0;
            for (int i = 0; i < size * size; i++) {
                raws.push_back(new XYZ());
            }
        }
    };
    struct processing_info {
        processing_info() {}

        Casting_Diagnostics CD;
        atomic<long long> parent_rays_cast{ 0 };
        atomic<long long> child_rays_cast{ 0 };

        int block_size = 0;
        int res_y = 0;
        int res_x = 0;
        int y_increment = 0;
        int x_increment = 0;
        int pixels_per_block = 0;
        int samples_per_pixel = 0;

        atomic<int> pixels_done{ 0 };
        atomic<int> iterations_done{ 0 };

        Camera* camera = nullptr;
        RayEngine* RE = nullptr;

        XYZ emit_coord;
        XYZ starting_coefficient;


        void add_casts() {
            parent_rays_cast++;
            child_rays_cast += CD.rays_processed - 1;
            CD.rays_processed = 0;
        }
        long long rays_cast() {
            return parent_rays_cast + child_rays_cast;
        }
        double percent() {
            return ((double)pixels_done) / (res_x * res_y);
        }

        ~processing_info() {
            //delete camera;
        }
    };
    atomic<bool> idle = true;
    atomic<bool> stop_thread = false;
    atomic<bool> halt_ops = false;
    atomic<int> mode = 0;

    thread worker;
    block* current;
    processing_info* process_info;
    vector<block*> blocks;

    RenderThread(processing_info* config) :process_info(config), worker(&RenderThread::main, this) {

    }

    void main() {
        process_info->parent_rays_cast = 0;
        process_info->child_rays_cast = 0;
        process_info->pixels_done = 0;
        srand(chrono::high_resolution_clock::now().time_since_epoch().count());
        gen.seed(rand(), rand(), rand(), rand());
        PackagedScene* PS = process_info->RE->data;
        while (true) {
            if (stop_thread.load(std::memory_order_relaxed)) {
                return;
            }
            if (!idle.load(std::memory_order_relaxed)) {
                switch (mode.load(std::memory_order_relaxed)) {
                case 0:
                    current = (blocks.back());
                    process_linearly();
                    current->done = true;
                    idle = true;
                    break;
                case 1:
                    int block_id = 0;
                    while (!halt_ops.load(std::memory_order_relaxed)) {
                        current = blocks[block_id];
                        process_iteratively();
                        block_id++;
                        block_id %= blocks.size();
                    }
                }
            }
            this_thread::sleep_for(chrono::milliseconds(1));
        }
    }
    void process_linearly() {
        for (int pixel_index = 0; pixel_index < process_info->pixels_per_block; pixel_index++) {
            iXY pixel_coordinate = get_coord(pixel_index);
            XYZ* output_link = get_output_link(pixel_index);
            for (int sub_index = 0; sub_index < process_info->samples_per_pixel; sub_index++) {
                if (halt_ops.load(std::memory_order_relaxed)) {
                    return;
                };
                XYZ ray_slope = slope_at(pixel_coordinate, sub_index);
                auto ray = PackagedRay(
                    process_info->emit_coord,
                    ray_slope,
                    0
                );
                (*output_link) += process_info->starting_coefficient * process_info->RE->process_ray(process_info->CD, ray);
                process_info->add_casts();
            }
            process_info->pixels_done++;
            current->pixels_done.store(pixel_index + 1, std::memory_order_relaxed);
        }
        current->done.store(true, std::memory_order_relaxed);
    }
    void process_iteratively() {
        int iterations_done = current->iterations_done.load(std::memory_order_relaxed);
        for (int pixel_index = 0; pixel_index < process_info->pixels_per_block; pixel_index++) {
            if (halt_ops.load(std::memory_order_relaxed)) {
                return;
            };
            iXY pixel_coordinate = get_coord(pixel_index);
            XYZ* output_link = get_output_link(pixel_index);
            XYZ ray_slope = random_slope_at(pixel_coordinate);
            auto ray = PackagedRay(
                process_info->emit_coord,
                ray_slope,
                0
            );
            (*output_link) *= (double)iterations_done / (iterations_done + 1);
            (*output_link) += 1.0 / (iterations_done + 1) * process_info->RE->process_ray(process_info->CD, ray);
            process_info->add_casts();
        }
        process_info->pixels_done++;
        current->iterations_done.store(iterations_done + 1, std::memory_order_relaxed);
    }
    iXY get_coord(int pixel_index) {
        int block_size = process_info->block_size;
        iXY inside_offset = iXY(pixel_index % block_size, block_size - pixel_index / block_size - 1);
        return current->offset + inside_offset;
    }
    XYZ* get_output_link(int pixel_index) {
        return current->raws[pixel_index];
    }
    XYZ slope_at(iXY coord, int sub_index) {
        return process_info->camera->slope_at(coord.X, coord.Y, sub_index);
    }
    XYZ random_slope_at(iXY coord) {
        return process_info->camera->random_slope_at(coord.X, coord.Y);
    }



};

class SceneManager {
public:
    Scene* scene;
    RayEngine RE;
    GUIHandler GUI;
    vector<vector<XYZ*>> raw_output;

    int thread_count = 12;
    vector<RenderThread*> threads;

    int current_resolution_x = 0;
    int current_resolution_y = 0;
    int subdivision_count = 0;
    int current_samples_per_pixel = 0;
    const int block_size = 8;

    int y_increment;
    int x_increment;

    chrono::steady_clock::time_point render_start;

    SceneManager(Scene* _scene) : scene(_scene) {}

    void render(int resolution_x, int resolution_y, int _subdivision_count = 1, int mode = 0) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;

        y_increment = ceil(current_resolution_y / block_size);
        x_increment = ceil(current_resolution_x / block_size);

        subdivision_count = _subdivision_count;
        current_samples_per_pixel = subdivision_count * subdivision_count;

        GUI = GUIHandler(current_resolution_x, current_resolution_y, block_size);

        render_start = chrono::high_resolution_clock::now();
        prep();
        //enqueue_rays();
        GUI.create_window();
        if (mode == 0) prepGUI();
        cout << endl << "+RAYCASTING+" << endl;
        spawn_threads();
        run_threaded_engine(mode);
        //_run_engine(true);

    }
    void hold_window() {
        GUI.hold_window();
    }

private:
    void prep() {
        cout << endl << "+PREPROCESSING+" << endl;

        auto prep_start = chrono::high_resolution_clock::now();

        for (int y = 0; y < current_resolution_y; y++) {
            vector<XYZ*> row;
            for (int x = 0; x < current_resolution_x; x++) {
                row.push_back(new XYZ(0, 0, 0));
            }
            raw_output.push_back(row);
        }

        scene->prep(current_resolution_x, current_resolution_y, subdivision_count);
        RE.load_scene(scene->package());

        auto prep_end = chrono::high_resolution_clock::now();

        cout << "Preprocessing done [" << chrono::duration_cast<chrono::milliseconds>(prep_end - prep_start).count() << "ms]" << endl;

    }
    void spawn_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread::processing_info* config = create_thread_config();
            RenderThread* thread = new RenderThread(config);
            string name = "RenderingThread-" + to_string(i);
            wstring wname = wstring(name.begin(), name.end());
            const wchar_t* wcname = wname.c_str();
            threads.push_back(thread);
            SetThreadPriority(thread->worker.native_handle(), 0);
            SetThreadDescription(thread->worker.native_handle(), wcname);

        }
    }
    struct render_stats {
        time_point start_time;
        time_point now;
        time_point last_refresh = chrono::high_resolution_clock::now();
        long long parent_rays = 0;
        long long child_rays = 0;
        long long parent_rays_last = 0;
        long long child_rays_last = 0;
        double percent_last = 0;
        long long pixels_done = 0;
        double percent_done = 0;
        vector<long long> parent_rays_rate_history;
        vector<long long> child_rays_rate_history;
        vector<double>    percent_rate_history;
        int remaining_threads;
        render_stats(int history_size_rays, int history_size_rate) {
            parent_rays_rate_history = vector<long long>(history_size_rays, 0);
            child_rays_rate_history = vector<long long>(history_size_rays, 0);
            percent_rate_history = vector<double>(history_size_rate, 0);
        }
    };
    struct pixel_strip {
        int strip_start;
        iXY block;
        vector<XYZ> outputs;
        pixel_strip(int start, iXY _block) {
            strip_start = start;
            block = _block;
        }
    };
    void assign_iterative_groups(vector<RenderThread::block*>& job_stack) {
        double jobs_per = (double)job_stack.size() / threads.size();
        double amount_given = 0;
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->mode = 1;
            double num_now = amount_given + jobs_per;
            double to_do = num_now - ceil(amount_given);
            amount_given = num_now;
            if (i == thread_count - 1) num_now = ceil(num_now);
            for (int j = 0; j < to_do; j++) {
                int job_index = rand() % job_stack.size();
                RenderThread::block* job = job_stack[job_index];
                render_thread->blocks.push_back(job);
                job_stack.erase(job_stack.begin() + job_index);
            }
        }
    }
    void start_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->idle = false;
        }
    }
    void halt_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->halt_ops = true;
        }
    }
    void unhalt_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->halt_ops = false;
        }
    }
    void reset_blocks(vector<RenderThread::block*>& blocks) {
        for (int i = 0;i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            block->iterations_done = 0;
            for (int j = 0; j < block->raws.size(); j++) {
                *block->raws[j] = XYZ();
            }
        }
    }
    void gather_thread_stats(render_stats& rStat) {
        rStat.parent_rays = 0;
        rStat.child_rays = 0;
        rStat.pixels_done = 0;
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            rStat.parent_rays += render_thread->process_info->parent_rays_cast;
            rStat.child_rays += render_thread->process_info->child_rays_cast;
            rStat.pixels_done += render_thread->process_info->pixels_done;
        }
    }
    void manage_threads(render_stats& rStat, vector<RenderThread::block*>& job_stack) {
        time_point now = chrono::high_resolution_clock::now();
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            bool is_idle = render_thread->idle.load(std::memory_order_relaxed);
            bool is_stopped = render_thread->stop_thread.load(std::memory_order_relaxed);
            if (is_idle && !is_stopped) {
                if (job_stack.size() > 0) {
                    RenderThread::block* job_block = job_stack.back();
                    GUI.set_focus_position(i, job_block->offset.X, job_block->offset.Y);
                    render_thread->blocks.push_back(job_block);
                    job_stack.pop_back();
                    render_thread->idle = false;
                }
                else {
                    render_thread->stop_thread = true;
                    rStat.remaining_threads--;
                    GUI.set_focus_position(i, -10000, -10000);
                }
            }
        }
    }
    iXY get_pixel_offset(int index) {
        return iXY(index % block_size, block_size - index / block_size - 1);
    }
    void get_block_updates(vector<RenderThread::block*>& blocks, vector<pixel_strip>& outputs) {
        for (int i = 0; i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            int pixels_done = block->pixels_done;
            int pixels_last = block->pixels_last;
            if (pixels_last < pixels_done) {
                iXY block_offset = block->offset;
                pixel_strip p_strip(pixels_last, block_offset);
                for (int pixel_index = pixels_last; pixel_index < pixels_done; pixel_index++) {
                    iXY internal_offset = get_pixel_offset(pixel_index);
                    iXY final_position = block_offset + internal_offset;
                    (*raw_output[final_position.Y][final_position.X]) = *block->raws[pixel_index];
                    XYZ color = ImageHandler::post_process_pixel(*raw_output[final_position.Y][final_position.X]);
                    p_strip.outputs.push_back(color);
                    //GUI.commit_pixel(color, final_position.X, final_position.Y);

                }
                outputs.push_back(p_strip);
            }
            block->pixels_last = pixels_done;
        }
    }
    void get_screen_update(vector<RenderThread::block*>& blocks) {
        for (int i = 0;i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            iXY block_offset = block->offset;
            for (int pixel_index = 0; pixel_index < block_size * block_size; pixel_index++) {
                iXY internal_offset = get_pixel_offset(pixel_index);
                iXY final_coord = block_offset + internal_offset;
                (*raw_output[final_coord.Y][final_coord.X]) = *block->raws[pixel_index];
                XYZ color = ImageHandler::post_process_pixel(*raw_output[final_coord.Y][final_coord.X]);
                GUI.commit_pixel(color, final_coord.X, final_coord.Y);
            }
        }
    }
    void pre_process_stats(render_stats& rStat) {
        int seconds = chrono::duration_cast<chrono::seconds>(rStat.now - rStat.start_time).count();
        int minutes = seconds / 60;
        int hours = minutes / 60;
        seconds %= 60;
        minutes %= 60;
        int micro_delta = chrono::duration_cast<chrono::microseconds>(rStat.now - rStat.last_refresh).count();
        int milli_delta = micro_delta / 1000;
        double second_ratio = ((double)1000000) / micro_delta;

        long long total_pixels = current_resolution_x * current_resolution_y;
        double percent = ((double)rStat.pixels_done) / total_pixels;
        rStat.percent_done = percent * 100;
        double percent_rate = (rStat.percent_done - rStat.percent_last) * second_ratio;


        long long delta_parent = rStat.parent_rays - rStat.parent_rays_last;
        long long delta_child = rStat.child_rays - rStat.child_rays_last;

        long long parent_rate = delta_parent * second_ratio;
        long long child_rate = delta_child * second_ratio;

        rStat.percent_rate_history.push_back(percent_rate);
        rStat.percent_rate_history.erase(rStat.percent_rate_history.begin());
        rStat.parent_rays_rate_history.push_back(parent_rate);
        rStat.parent_rays_rate_history.erase(rStat.parent_rays_rate_history.begin());
        rStat.child_rays_rate_history.push_back(child_rate);
        rStat.child_rays_rate_history.erase(rStat.child_rays_rate_history.begin());
    }
    void post_process_stats(render_stats& rStat) {
        rStat.parent_rays_last = rStat.parent_rays;
        rStat.child_rays_last = rStat.child_rays;
        rStat.percent_last = rStat.percent_done;
        rStat.last_refresh = chrono::high_resolution_clock::now();
    }
    void display_stats(render_stats& rStat) {
        int seconds = chrono::duration_cast<chrono::seconds>(rStat.now - rStat.start_time).count();
        int minutes = seconds / 60;
        int hours = minutes / 60;
        seconds %= 60;
        minutes %= 60;
        int micro_delta = chrono::duration_cast<chrono::microseconds>(rStat.now - rStat.last_refresh).count();
        int milli_delta = micro_delta / 1000;
        double second_ratio = ((double)1000000) / micro_delta;

        double percent_rate_average = 0;
        long long parent_rate_average = 0;
        long long child_rate_average = 0;

        for (int i = 0; i < rStat.parent_rays_rate_history.size(); i++) {
            parent_rate_average += rStat.parent_rays_rate_history[i];
        }
        parent_rate_average /= rStat.parent_rays_rate_history.size();
        for (int i = 0; i < rStat.child_rays_rate_history.size(); i++) {
            child_rate_average += rStat.child_rays_rate_history[i];
        }
        child_rate_average /= rStat.child_rays_rate_history.size();
        for (int i = 0; i < rStat.percent_rate_history.size(); i++) {
            percent_rate_average += rStat.percent_rate_history[i];
        }
        percent_rate_average /= rStat.percent_rate_history.size();

        string primary_per_sec = intToEng(parent_rate_average);
        string secondary_per_sec = intToEng(child_rate_average);
        string total_per_sec = intToEng(parent_rate_average + child_rate_average);



        string format_string = "[%02i:%02i:%02i][%7ims] - ([%.5srays/s Primary][%.5srays/s Secondary])[%.5srays/s] - [%4.1f%% (+%4.3f/s)]\r";
        printf(format_string.c_str(), hours, minutes, seconds, milli_delta, primary_per_sec.c_str(), secondary_per_sec.c_str(), total_per_sec.c_str(), rStat.percent_done, percent_rate_average);

    }
    void process_strips(vector<pixel_strip> p_strips) {
        for (int o_index = 0; o_index < p_strips.size(); o_index++) {
            pixel_strip& p_strip = p_strips[o_index];
            iXY block_offset = p_strip.block;
            int strip_start = p_strip.strip_start;
            vector<XYZ>& outputs = p_strip.outputs;
            for (int i = 0; i < outputs.size(); i++) {
                XYZ color = outputs[i];
                iXY internal_offset = get_pixel_offset(i + strip_start);
                iXY final_coord = internal_offset + block_offset;
                GUI.commit_pixel(color, final_coord.X, final_coord.Y);
            }

        }
        p_strips.clear();
    }
    void refresh_display() {
        GUI.commit_canvas();
        GUI.handle_events();
    }
    void run_threaded_engine(int mode = 0) {
        render_stats rStat(120, 300);
        rStat.start_time = chrono::high_resolution_clock::now();
        render_start = rStat.start_time;
        rStat.remaining_threads = thread_count;
        vector<RenderThread::block*> blocks;
        vector<RenderThread::block*> job_stack;
        vector<pixel_strip> update_strips;
        for (int block_index_y = 0; block_index_y < y_increment; block_index_y++) {
            for (int block_index_x = 0; block_index_x < x_increment; block_index_x++) {
                job_stack.push_back(new RenderThread::block(block_size, block_index_x, block_index_y));
                blocks.push_back(job_stack.back());
            }
        }
        switch (mode) {
        case 0:
            while (true) {
                rStat.now = chrono::high_resolution_clock::now();
                manage_threads(rStat, job_stack);
                if (chrono::duration_cast<chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 16) {
                    gather_thread_stats(rStat);
                    get_block_updates(blocks, update_strips);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    process_strips(update_strips);
                    refresh_display();
                }
                if (rStat.remaining_threads <= 0) {
                    break;
                }
            }
            break;
        case 1:
            assign_iterative_groups(job_stack);
            start_threads();
            while (true) {
                rStat.now = chrono::high_resolution_clock::now();
                if (chrono::duration_cast<chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 16) {
                    gather_thread_stats(rStat);
                    get_screen_update(blocks);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    refresh_display();
                }

            }
            break;
        case 2:
            assign_iterative_groups(job_stack);
            start_threads();
            while (true) {
                rStat.now = chrono::high_resolution_clock::now();
                if (chrono::duration_cast<chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 64) {
                    halt_threads();
                    gather_thread_stats(rStat);
                    get_screen_update(blocks);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    refresh_display();
                    reset_blocks(blocks);
                    unhalt_threads();
                }
            }
        }
        for (int i = 0; i < thread_count; i++) {
            threads[i]->stop_thread = true;
            threads[i]->worker.join();
        }
        for (int i = 0; i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            iXY block_offset = block->offset;
            for (int pixel_index = 0; pixel_index < block->raws.size(); pixel_index++) {
                iXY internal_offset = iXY(pixel_index % block_size, block_size - pixel_index / block_size - 1);
                iXY final_position = block_offset + internal_offset;
                (*raw_output[final_position.Y][final_position.X]) = *block->raws[pixel_index];
                XYZ color = ImageHandler::post_process_pixel(*raw_output[final_position.Y][final_position.X]);
                GUI.commit_pixel(color, final_position.X, final_position.Y);
            }
        }
        refresh_canvas();
        cout << endl;

        auto end_time = chrono::high_resolution_clock::now();
        double millis = chrono::duration_cast<chrono::milliseconds>(end_time - rStat.start_time).count();
        cout << "Ray Casting Done [" << to_string(millis) << "ms]" << endl;
        //cout << "Approximate speed: " << intToEng((rays_cast / (millis / 1000.0))) << " rays/s";
    }
    void prepGUI() {
        for (int i = 0; i < thread_count; i++) {
            GUI.add_focus();
        }
    }
    RenderThread::processing_info* create_thread_config() {
        RenderThread::processing_info* process_info = new RenderThread::processing_info();

        int max_bounces = 2;

        process_info->emit_coord = scene->camera->position;
        process_info->starting_coefficient = XYZ(1, 1, 1) / current_samples_per_pixel;

        process_info->res_y = current_resolution_y;
        process_info->res_x = current_resolution_x;
        process_info->y_increment = y_increment;
        process_info->x_increment = x_increment;
        process_info->block_size = block_size;
        process_info->samples_per_pixel = current_samples_per_pixel;
        process_info->pixels_per_block = block_size * block_size;

        process_info->camera = scene->camera;//camera->clone();
        process_info->RE = new RayEngine(RE);

        return process_info;
    }


    void refresh_canvas() {
        for (int y = 0; y < current_resolution_y; y++) {
            for (int x = 0; x < current_resolution_x; x++) {
                XYZ color = ImageHandler::post_process_pixel(*raw_output[y][x]);
                GUI.commit_pixel(color, x, y);
            }
        }
        GUI.commit_canvas();
    }
};

class FileManager {
public:
    FileManager() {}
    static void writeRawFile(vector<vector<XYZ*>>* data, string fName) {
        ofstream file(fName, ios::out | ios::binary | ios::trunc);
        FileManager::writeRaw(data, file);
        file.close();
    }
    static void writeRaw(vector<vector<XYZ*>>* data, ofstream& file) {
        unsigned int resX = data->at(0).size();
        unsigned int resY = data->size();
        cout << resX << " " << resY << endl;
        file.write((char*)&resX, sizeof(unsigned int));
        file.write((char*)&resY, sizeof(unsigned int));
        for (vector<XYZ*>& data_row : *data) {
            for (XYZ* pixel_ptr : data_row) {
                file.write((char*)pixel_ptr, sizeof(XYZ));
            }
        }
    }
    static GUIHandler* openRawFile(string fName) {
        ifstream file(fName, ios::out | ios::binary);
        return openRaw(file);
    }
    static GUIHandler* openRaw(ifstream& file) {
        vector<vector<XYZ*>>* data = FileManager::readRaw(file);
        file.close();
        ImageHandler::postProcessRaw(data);
        GUIHandler* GUI = ImageHandler::openImageVector(data);
        return GUI;
    }
    static vector<vector<XYZ*>>* readRaw(ifstream& file) {
        vector<vector<XYZ*>>* data = new vector<vector<XYZ*>>();
        int resolutionX = 0;
        int resolutionY = 0;
        file.read((char*)&resolutionX, sizeof(int));
        file.read((char*)&resolutionY, sizeof(int));
        cout << resolutionX << " " << resolutionY << endl;
        for (int y = 0; y < resolutionY; y++) {
            vector<XYZ*> output_row;
            for (int x = 0; x < resolutionX; x++) {
                output_row.push_back(new XYZ(0, 0, 0));
                file.read((char*)output_row[x], sizeof(XYZ));
            }
            data->push_back(output_row);
        }
        return data;
    }
    static Mesh* loadObjFile(string fName, Material* mat) {
        cout << padString("[File] Loading resource: " + fName, ".", 100);
        ifstream file(fName, ios::in);
        Mesh* return_value = loadObj(file, mat);
        cout << "[Done]" << endl;
        return return_value;
    }
    static Mesh* loadObj(ifstream& file, Material* mat) {
        Mesh* mesh = new Mesh();
        vector<XYZ> verts;
        vector<XY> UV_coords;
        for (std::string line; getline(file, line);) {
            string line_type = line.substr(0, line.find(" "));
            vector<string> words = split_string(line, ' ');
            if (line_type == "#") {
                continue;
            }
            if (line_type == "mtllib") {
                continue;
            }
            if (line_type == "o") {
                continue;
            }
            if (line_type == "v") {
                verts.push_back(XYZ(
                    stod(words[1]),
                    stod(words[2]),
                    stod(words[3])
                ));
                //cout << verts.back() << endl;
            }
            if (line_type == "vt") {
                UV_coords.push_back(XY(
                    stod(words[1]), stod(words[2])
                ));
            }
            if (line_type == "f") {
                vector<string> sub_words_1 = split_string(words[1], '/');
                vector<string> sub_words_2 = split_string(words[2], '/');
                vector<string> sub_words_3 = split_string(words[3], '/');
                int v_i0 = stoi(sub_words_1[0]) - 1;
                int v_i1 = stoi(sub_words_2[0]) - 1;
                int v_i2 = stoi(sub_words_3[0]) - 1;
                XY  vt1 = XY(-1, -1);
                XY  vt2 = XY(-1, -1);
                XY  vt3 = XY(-1, -1);
                if (sub_words_1.size() > 1) {
                    if (sub_words_1[1] != "") {
                        int vt_i0 = stoi(sub_words_1[1]) - 1;
                        int vt_i1 = stoi(sub_words_2[1]) - 1;
                        int vt_i2 = stoi(sub_words_3[1]) - 1;
                        vt1 = UV_coords[vt_i0];
                        vt2 = UV_coords[vt_i1];
                        vt3 = UV_coords[vt_i2];
                    }
                }
                Tri f = Tri(
                    verts[v_i0], verts[v_i1], verts[v_i2],
                    vt1, vt2, vt3,
                    mat);
                mesh->addTri(f);
            }
        }
        return mesh;
    }
    static Mesh* loadTriFile(string fName, Material* mat) {
        ifstream file(fName, ios::in);
        return loadTri(file, mat);
    }
    static Mesh* loadTri(ifstream& file, Material* mat) {
        Mesh* mesh = new Mesh();
        vector<XYZ> verts;
        for (std::string line; getline(file, line);) {
            vector<string> words = split_string(line, ' ');
            if (words.size() > 1) {
                double v1_x = stod(words[0]);
                double v1_y = stod(words[1]);
                double v1_z = stod(words[2]);
                double v2_x = stod(words[3]);
                double v2_y = stod(words[4]);
                double v2_z = stod(words[5]);
                double v3_x = stod(words[6]);
                double v3_y = stod(words[7]);
                double v3_z = stod(words[8]);
                Tri f1 = Tri(
                    XYZ(v1_x, v1_y, v1_z),
                    XYZ(v2_x, v2_y, v2_z),
                    XYZ(v3_x, v3_y, v3_z),
                    mat);
                mesh->addTri(f1);
                //cout << verts.back() << endl;
            }
        }
        return mesh;
    }
    static SceneManager* loadGLTFFile(string fName) {
        tinygltf::Model model;
        tinygltf::TinyGLTF loader;
        std::string err;
        std::string warn;

        cout << padString("[File] Loading scene " + fName + "...\n", ".", 0);
        //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, argv[1]);
        bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, fName); // for binary glTF(.glb)

        if (!warn.empty()) {
            printf("Warn: %s\n", warn.c_str());
        } if (!err.empty()) {
            printf("Err: %s\n", err.c_str());
        } if (!ret) {
            throw exception("Failed to parse glTF\n");
        }

        SceneManager* out = parseGLTFData(model);

        cout << "[Load] loading done";

        return out;

    }
    template<typename T> struct iterator_pair {
        _Vector_iterator<_Vector_val<_Simple_types<T>>> begin;
        _Vector_iterator<_Vector_val<_Simple_types<T>>> end;
        T* begin_ptr;
        iterator_pair() {}
    };
    struct buffer_accessor {
        tinygltf::Model& data;
        vector<char> swizzle;
        buffer_accessor(tinygltf::Model& _data) : data(_data) {}
        iterator_pair<unsigned char> get_buffer(int accessor_index) {
            iterator_pair<unsigned char> returner;
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            int buffer_index = view.buffer;
            int start_index = view.byteOffset;
            int end_index = start_index + view.byteLength;

            auto& buffer = data.buffers[buffer_index];
            returner.begin = buffer.data.begin() + start_index;
            returner.end = buffer.data.begin() + end_index;
            returner.begin_ptr = &(*returner.begin);

            return returner;

        }
        void get_float_data(int accessor_index, vector<float>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5126);
            assert(view.byteLength % 4 == 0);

            target.clear();
            target.resize(view.byteLength / 4, 0);
            std::copy_n(reinterpret_cast<float*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_ushort_data(int accessor_index, vector<unsigned short>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5123);
            assert(view.byteLength % 2 == 0);

            target.clear();
            target.resize(view.byteLength / 2, 0);
            std::copy_n(reinterpret_cast<unsigned short*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_uint_data(int accessor_index, vector<unsigned int>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5125);
            assert(view.byteLength % 4 == 0);

            target.clear();
            target.resize(view.byteLength / 4, 0);
            std::copy_n(reinterpret_cast<unsigned int*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_data_any_int(int accessor_index, vector<unsigned int>& target) {
            auto& accessor = data.accessors[accessor_index];
            vector<unsigned short> intermediate;
            switch (accessor.componentType) {
            case 5125:
                get_uint_data(accessor_index, target);
                break;
            case 5123:
                get_ushort_data(accessor_index, intermediate);
                target.clear();
                target.reserve(intermediate.size());
                for (int i = 0; i < intermediate.size(); i++) {
                    target.push_back((int)intermediate[i]);
                }
                break;
            default:
                assert(0);
                break;
            }
        }
        void get_vec4_data(int accessor_index, vector<Quat>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 4);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(Quat(float_data[i * 4], float_data[i * 4 + 1], float_data[i * 4 + 2], float_data[i * 4 + 3]));
            }
        }
        void get_vec4_data_strip(int accessor_index, vector<XYZ>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 4);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XYZ(float_data[i * 4], float_data[i * 4 + 1], float_data[i * 4 + 2]));
            }
        }
        void get_vec3_data(int accessor_index, vector<XYZ>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 3);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XYZ(float_data[i * 3], float_data[i * 3 + 1], float_data[i * 3 + 2]));
            }
        }
        void get_vec2_data(int accessor_index, vector<XY>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 2);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XY(float_data[i * 2], float_data[i * 2 + 1]));
            }
        }
    };
    struct TextureWrapper {
        Texture<XYZ>* xyz;
        Texture<float>* f;
    };
    template <int t> static float converter_float(float x, float y, float z) {
        if (t == 0) return x;
        if (t == 1) return y;
        if (t == 2) return z;
    }
    static XYZ converter_xyz(float x, float y, float z) {
        return XYZ(x, y, z);
    }
    static XYZ converter_xyz_sRGB(float x, float y, float z) {
        return ImageHandler::apply_gamma(XYZ(x, y, z),1/2.2);
    }
    static XYZ converter_xzy(float x, float y, float z) {
        return XYZ(x, z, y);
    }
    template<float scale>
    class converter_normal_scaler {
    public:
        static XYZ func(float x, float y, float z) {
            return (XYZ::normalize(XYZ(x, y, z)) * 2 - 1) * XYZ(scale, scale, 1/2.2);
        }
    
    };
    template <typename T> static Texture<T>* loadTexture(tinygltf::Model& data, int index, T(*converter)(float, float, float)) {
        auto& texture_data = data.textures[index];
        int image_index = texture_data.source;
        auto& image_data = data.images[image_index];
        int len = image_data.image.size();
        int res_x = image_data.width;
        int res_y = image_data.height;
        int component = image_data.component;

        cout << endl;
        cout << "[Load] Loading texture \"" << image_data.name << "\"\r" << flush;

        if (component != 3 && component != 4) throw exception("illegal texture!");

        int is_addressable = -1;
        if (is_same<T, XYZ>()) {
            is_addressable = 1;
        }
        if (is_same<T, float>()) {
            is_addressable = 0;
        }
        if (is_addressable == -1) throw exception("illegal texture typing!");

        Texture<T>* t = new Texture<T>(image_data.image.data(), res_x, res_y, component,converter);

        

        image_data.image.clear(); //freeing memory

        cout << padString("", " ", 100) << "\x1b[A\r";
        return t;
    }

    static pair<XY, XY> fetch_transform(tinygltf::TextureInfo& Ti) {
        XY scale = XY(1, 1);
        XY offset = XY(0, 0);
        if (Ti.extensions.count("KHR_texture_transform")) {
            auto& ext = Ti.extensions["KHR_texture_transform"];
            auto& offset_object = ext.Get("offset");
            auto& scale_object = ext.Get("scale");
            offset.X = offset_object.Get(0).GetNumberAsDouble();
            offset.Y = offset_object.Get(1).GetNumberAsDouble();
            scale.X = scale_object.Get(0).GetNumberAsDouble();
            scale.Y = scale_object.Get(1).GetNumberAsDouble();
        }
        return { scale,offset };
    }
    static pair<XY, XY> fetch_transform(tinygltf::NormalTextureInfo& Ti) {
        XY scale = XY(1, 1);
        XY offset = XY(0, 0);
        if (Ti.extensions.count("KHR_texture_transform")) {
            auto& ext = Ti.extensions["KHR_texture_transform"];
            auto& offset_object = ext.Get("offset");
            auto& scale_object = ext.Get("scale");
            offset.X = offset_object.Get(0).GetNumberAsDouble();
            offset.Y = offset_object.Get(1).GetNumberAsDouble();
            scale.X = scale_object.Get(0).GetNumberAsDouble();
            scale.Y = scale_object.Get(1).GetNumberAsDouble();
        }
        return { scale,offset };
    }
    static SceneManager* parseGLTFData(tinygltf::Model data) {

        const bool print_in_place = true;

        Scene* scene = new Scene();
        SceneManager* SM = new SceneManager(scene);
        vector<Material*> materials;
        vector<Mesh*> meshes;
        vector<Camera*> cameras;
        vector<Object*> objects;
        vector<Light*> lights;
        map<int, pair<int, int>> node_lookup;
        vector<char> swizzle = { 0,-2,1 };
        vector<char> rot_swizzle = { 0,-2, 1, 3 };
        vector<char> scale_swizzle = { 0,2,1 };
        buffer_accessor buf_accessor = buffer_accessor(data);
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Textures loaded" << endl;
        for (auto& material_data : data.materials) {
            string name = material_data.name;
            cout << padString("", " ", 100) << "\r";
            cout << "[Load] Loading material \"" << name << "\"\r" << flush;
            Material* mat = new Material();
            auto pbr = material_data.pbrMetallicRoughness;
            auto emissive = material_data.emissiveFactor;

            mat->color.set_static(XYZ(pbr.baseColorFactor));
            mat->metallic.set_static(pbr.metallicFactor);
            mat->roughness.set_static(pbr.roughnessFactor);
            mat->emissive.set_static(XYZ(emissive) * 25);

            if (pbr.baseColorTexture.index != -1) {
                int t_index = pbr.baseColorTexture.index;
                mat->UV_map_index = pbr.baseColorTexture.texCoord;
                Texture<XYZ>* tex = loadTexture<XYZ>(data, t_index, converter_xyz_sRGB);
                auto transform_data = fetch_transform(pbr.baseColorTexture);
                tex->set_transform(transform_data.first, transform_data.second);
                mat->color.set_texture(tex);
            }
            if (pbr.metallicRoughnessTexture.index != -1) {
                int t_index = pbr.metallicRoughnessTexture.index;
                mat->UV_map_index = pbr.baseColorTexture.texCoord;
                Texture<XYZ>* tex = loadTexture<XYZ>(data, t_index, converter_xyz);
                auto& texture_data = data.textures[t_index];
                auto& image_data = data.images[texture_data.source];
                cout << endl << "[Load] Postprocessing texture \"" << image_data.name << "\"\r" << flush;

                auto transform_data = fetch_transform(pbr.metallicRoughnessTexture);
                tex->set_transform(transform_data.first, transform_data.second);

                Texture<float>* metallic = tex->export_channel(2);
                Texture<float>* roughness = tex->export_channel(1);
                

                if (pbr.metallicFactor != 0) {
                    if (pbr.metallicFactor != 1) {
                        for (int i = 0; i < metallic->size(); i++) {
                            metallic->data[i] *= pbr.metallicFactor;
                        }
                    }
                    mat->metallic.set_texture(metallic);
                }
                if (pbr.roughnessFactor != 0) {
                    if (pbr.roughnessFactor != 1) throw exception("Illegal roughness factor");
                    mat->roughness.set_texture(roughness);
                }
                cout << padString("", " ", 100) << "\x1b[A\r";
            }
            if (material_data.normalTexture.index != -1) {
                int t_index = material_data.normalTexture.index;
                mat->UV_map_index = material_data.normalTexture.texCoord;
                const float scale = material_data.normalTexture.scale;
                Texture<XYZ>* tex = loadTexture<XYZ>(data, t_index, converter_xyz);
                auto& texture_data = data.textures[t_index];
                auto& image_data = data.images[texture_data.source];
                cout << endl << "[Load] Postprocessing normal map \"" << image_data.name << "\"\r" << flush;
                int count = tex->size();
                for (int i = 0; i < count; i++) {
                    XYZ in = tex->data[i];
                    XYZ out = XYZ::normalize((in * 2 - 1) * XYZ(scale, scale, 1.0));
                    tex->data[i] = out;
                }

                auto transform_data = fetch_transform(material_data.normalTexture);
                tex->set_transform(transform_data.first, transform_data.second);

                mat->normal.set_texture(tex);
                mat->use_normals = 1;
                cout << padString("", " ", 100) << "\x1b[A\r";
            }
            if (material_data.extensions.count("KHR_materials_specular")) {
                auto specular_iterator = material_data.extensions.find("KHR_materials_specular");
                if (specular_iterator != material_data.extensions.end()) {
                    auto specular_value = specular_iterator->second.Get("specularColorFactor");
                    if (specular_value.Type() == 0) {
                        mat->specular.set_static(0.0f);
                    }
                    else {
                        auto specular_0 = specular_value.Get(0).GetNumberAsDouble();
                        auto specular_1 = specular_value.Get(1).GetNumberAsDouble();
                        auto specular_2 = specular_value.Get(2).GetNumberAsDouble();
                        mat->specular.set_static(specular_0);
                    }
                }
            }
            else {
                mat->specular = 1;
            }
            if (material_data.extensions.count("KHR_materials_ior")) {
                auto ior_iterator = material_data.extensions.find("KHR_materials_ior");
                if (ior_iterator != material_data.extensions.end()) {
                    auto ior_value = ior_iterator->second.Get("ior");
                    auto ior = ior_value.GetNumberAsDouble();
                    mat->IOR.set_static(ior);
                }
            }
            else {
                mat->IOR = 1.5;
            }
            materials.push_back(mat);
        }
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Materials loaded" << endl;
        for (auto& model_data : data.meshes) {
            auto& primitives = model_data.primitives;
            Mesh* mesh = new Mesh();
            mesh->name = model_data.name;
            cout << padString("", " ", 100) << "\r";
            cout << "[Load] Loading mesh \"" << mesh->name << "\"\r" << flush;
            for (auto& primitive : primitives) {
                auto& attributes = primitive.attributes;
                Material* mat;
                if (primitive.material > -1) {
                    mat = materials[primitive.material];
                }
                else {
                    mat = new Material();
                }

                assert(attributes.count("POSITION") > 0);
                auto  position_accessor_index = attributes["POSITION"];
                auto texcoord_accessor_index = attributes["TEXCOORD_"+to_string(mat->UV_map_index)];
                auto normal_accessor_index = attributes["NORMAL"];
                auto tangent_accessor_index = attributes["TANGENT"];
                auto tangent_accessor = data.accessors[tangent_accessor_index];

                vector<XYZ> position_data;
                buf_accessor.get_vec3_data(position_accessor_index, position_data);

                vector<XY> texcoord_data;
                buf_accessor.get_vec2_data(texcoord_accessor_index, texcoord_data);

                bool smooth_shading = normal_accessor_index != 0;
                smooth_shading &= 1;
                vector<XYZ> normal_data;
                if (smooth_shading) {
                    buf_accessor.get_vec3_data(normal_accessor_index, normal_data);
                }


                //vector<XYZ> tangent_data;
                //if (tangent_accessor.type == 4) {
                //    buf_accessor.get_vec4_data_strip(tangent_accessor_index, tangent_data);
                //}
                //else {
                //    buf_accessor.get_vec3_data(tangent_accessor_index, tangent_data);
                //}

                assert(primitive.indices >= 0);
                auto  indices_accessor_index = primitive.indices;
                auto& indices_accessor = data.accessors[position_accessor_index];

                vector<unsigned int> indices;
                buf_accessor.get_data_any_int(indices_accessor_index, indices);

                assert(indices.size() % 3 == 0);
                int tri_count = indices.size() / 3;

                mesh->tris.reserve(mesh->tris.size() + tri_count);
                for (int i = 0; i < tri_count; i++) {
                    int look_1 = indices[i * 3 + 0];
                    int look_2 = indices[i * 3 + 1];
                    int look_3 = indices[i * 3 + 2];
                    XYZ v1 = position_data[look_1].swizzle(swizzle);
                    XYZ v2 = position_data[look_2].swizzle(swizzle);
                    XYZ v3 = position_data[look_3].swizzle(swizzle);
                    XY uv1 = XY(0,1)-XY(-1,1)*texcoord_data[look_1];
                    XY uv2 = XY(0,1)-XY(-1,1)*texcoord_data[look_2];
                    XY uv3 = XY(0,1)-XY(-1,1)*texcoord_data[look_3];
                    
                    //cout << v1 << " " << v2 << " " << v3 << " " << uv1 << " " << uv2 << " " << uv3 << endl;
                    Tri T = Tri(v1, v2, v3, uv1,uv2,uv3, mat);
                    if (smooth_shading) {
                        T.n1 = normal_data[look_1].swizzle(swizzle);
                        T.n2 = normal_data[look_2].swizzle(swizzle);
                        T.n3 = normal_data[look_3].swizzle(swizzle);
                    }
                    mesh->addTri(T);
                }
            }
            meshes.push_back(mesh);
        }
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Meshes loaded" << endl;
        for (auto& camera_data : data.cameras) {
            assert(camera_data.type == "perspective");
            auto& perspective = camera_data.perspective;
            auto aspect_ratio = perspective.aspectRatio;
            Lens* lens = new RectLens(aspect_ratio, 1);
            float distance = 0.5 / tan(0.5 * perspective.yfov);
            Camera* camera = new Camera(XYZ(), lens, XYZ(0, distance, 0));
            cameras.push_back(camera);
        }
        int node_index = -1;
        for (auto& node_data : data.nodes) {
            node_index++;
            if (node_data.camera > -1) {
                XYZ scale = XYZ(1);
                XYZ translation = XYZ();
                Quat rotation = Quat();
                if (node_data.scale.size() > 0) scale = XYZ(node_data.scale, swizzle);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                cameras[node_data.camera]->position += translation;
                cameras[node_data.camera]->rotation = rotation;
                node_lookup[node_index] = pair<int, int>(0, node_data.camera);
                continue;
            }
            if (node_data.mesh > -1) {
                XYZ scale = XYZ(1);
                XYZ translation = XYZ(0,0,0);
                Quat rotation = Quat(0,0,0,1);
                if (node_data.scale.size() > 0) scale = XYZ(node_data.scale, scale_swizzle);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                Object* obj = new Object(XYZ(), scale);
                obj->name = node_data.name;
                obj->addMesh(meshes[node_data.mesh]);
                obj->_registerRotation(rotation);
                obj->_registerMove(translation);
                objects.push_back(obj);
                node_lookup[node_index] = pair<int, int>(1, objects.size() - 1);
                continue;
            }
            if (node_data.light > -1) {
                XYZ translation = XYZ();
                Quat rotation = Quat(0, 0, 0, 1);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                Light* li= nullptr;
                auto& light_data = data.lights[node_data.light];
                if (light_data.type == "directional") {
                    li = new SunLight();
                    XYZ dir = XYZ(0, -1, 0);
                    dir = Quat::applyRotation(dir, rotation);
                    cout << dir << endl;
                    li->position = translation;
                    li->rotation = dir;
                    li->emission = XYZ(light_data.color) * light_data.intensity / 50000;
                }
                if (li != nullptr) {
                    lights.push_back(li);
                }
                else {
                    cout << "Unrecognized camera type: " << light_data.type;
                }
                node_lookup[node_index] = pair<int, int>(2, lights.size() - 1);
                continue;
            }
        }
        for (auto& object_lookup : node_lookup) {
            int node_index = object_lookup.first;
            int object_index = object_lookup.second.second;
            Object* obj = objects[object_index];
            auto& node_data = data.nodes[node_index];
            for (int child_node_index : node_data.children) {
                auto& lookup_entry = node_lookup[child_node_index];
                int child_object_index = lookup_entry.second;
                Object* child_obj = objects[child_object_index];
                obj->addChild(child_obj);
            }

        }
        int default_scene = data.defaultScene;
        auto scene_data = data.scenes[default_scene];
        //for (auto& node_index : scene_data.nodes) {
        for(int i = 0; i < node_lookup.size(); i++){ //this is 100% incorrect, but I was getting weird import issues
            auto& node_register = node_lookup[i];
            int node_type = node_register.first;
            int lookup_index = node_register.second;
            if (node_type == 0) {
                Camera* camera = cameras[lookup_index];
                scene->register_camera(camera);
            }
            if (node_type == 1) {
                Object* obj = objects[lookup_index];
                if (obj->parent == nullptr) {
                    scene->register_object(obj);
                }
            }
            if (node_type == 2) {
                Light* li = lights[lookup_index];
                scene->register_light(li);
            }
        }

        scene->camera = scene->cameras[0];
        return SM;
    }

};

int main(int argc, char* argv[])
{
    srand(0);

    gen.seed(rand(), rand(), rand(), rand());

    TCHAR buffer[MAX_PATH] = { 0 };
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::wstring::size_type pos = std::wstring(buffer).find_last_of(L"\\/");
    auto str = std::wstring(buffer).substr(0, pos);
    cout << string(str.begin(), str.end()) << endl;

    VecLib::prep(gen);
    map<string, string> arg_bindings;
    if (argc > 1) {
        string prev = argv[1];
        for (int i = 2; i < argc; i++) {
            arg_bindings[prev] = argv[i];
            prev = argv[i];
        }
        arg_bindings[prev] = "";
    }
    int mode = 0;
    string fpath = "C:\\Users\\Charlie\\Documents\\models\\Scenes\\dining.glb";
    //string fpath = "C:\\Users\\Charlie\\Documents\\models\\Scenes\\title.glb";
    string host_ip = "127.0.0.1";
    string host_port = "15413";
    string matcher_ip = "";
    string matcher_port = "15413";
    int scalar = 2;
    int resolution_x = 1080/scalar;
    int resolution_y = 1080/scalar;
    int subdivisions = 1;
    if (arg_bindings.count("-f")) fpath = arg_bindings["-f"];
    if (arg_bindings.count("-H")) mode = 1;
    if (arg_bindings.count("-S")) mode = 2;
    try {
        if (arg_bindings.count("-rX")) resolution_x = stoi(arg_bindings["-rX"]);
        if (arg_bindings.count("-rY")) resolution_y = stoi(arg_bindings["-rY"]);
        if (arg_bindings.count("-sd")) subdivisions = stoi(arg_bindings["-sd"]);
    }
    catch (exception e) {
        cout << "invalid numerical parameter" << endl;
    }
    //GUIHandler* GUI = FileManager::openRawFile("outputs.raw");

    //GUI->hold_window();
    //SceneManager* scene_manager = load_cornell_box();

    SceneManager* scene_manager;


    switch (mode) {
    case 0:
        scene_manager = FileManager::loadGLTFFile(fpath);
        scene_manager->render(resolution_x, resolution_y, subdivisions, 1);
        scene_manager->hold_window();
        break;
    }


    //scene_manager->render(1080/4, 1080/4, 2);

    //cout << endl << "writing out raw......." << flush;
    //FileManager::writeRawFile(&scene_manager->raw_output, "outputs.raw");
    //cout << "done" << endl;

}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
