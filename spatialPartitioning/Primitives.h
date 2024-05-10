#pragma once
#include "commons.h"
#include "Materials.h"
#include "Transformation.h"
#include "VecLib.h"
#include <math.h>
#include "Coredefs.h"

using std::pair;
using std::exception;
using std::string;

class Primitive {
public:
    uint32_t material;
    //pair<XYZ, XYZ> bounds;
    XYZ origin;
    Primitive(short _material) :material(_material) {}
    virtual pair<XYZ, XYZ> get_bounds() { return pair<XYZ, XYZ>(); }
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
    XYZ p1p3;
    XYZ p1p2;
    XY UV_basis;
    XY U_delta;
    XY V_delta;
    uint32_t material;
    XYZ normal;
    XYZ origin_offset;
    XYZ n1;
    XYZ n2;
    XYZ n3;
    XYZ t1;
    XYZ t2;
    XYZ t3;
    XYZ b1;
    XYZ b2;
    XYZ b3;
    PackagedTri() {}
    PackagedTri(const XYZ& _p1, const XYZ& _p2, const XYZ& _p3, XY UV_base, XY UV_U, XY UV_V, short _material) {
        material = _material;
        p1 = _p1;
        p1p3 = _p3 - p1;
        p1p2 = _p2 - p1;
        U_delta = UV_U - UV_base;
        V_delta = UV_V - UV_base;
        if (U_delta == 0) U_delta = XY(1,0);
        if (V_delta == 0) V_delta = XY(0, 1);
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
        t1 = XYZ::normalize(t1 - n1 * XYZ::dot(t1, n1));
        t2 = XYZ::normalize(t2 - n2 * XYZ::dot(t2, n1));
        t3 = XYZ::normalize(t3 - n3 * XYZ::dot(t3, n1));
        b1 = XYZ::cross(n1, t1);
        b2 = XYZ::cross(n2, t2);
        b3 = XYZ::cross(n3, t3);
    }
};

struct PTri_AVX { //literally 0.7kB of memory per object. What have I done.
    m256_vec3 p1;
    m256_vec3 p1p2;
    m256_vec3 p1p3;
    XYZ origin_offset[8];
    uint32_t materials[8];
    int index = 0;
    char size = 0;
    PTri_AVX(vector<PackagedTri> tris) {
        vector<XYZ> _p1;
        vector<XYZ> _p1p2;
        vector<XYZ> _p1p3;
        vector<XYZ> _origin_offset;
        size = tris.size();
        for (int i = 0; i < tris.size() && i < 8; i++) {
            PackagedTri& tri = tris[i];
            _p1.push_back(tri.p1);
            _p1p2.push_back(tri.p1p2);
            _p1p3.push_back(tri.p1p3);
            _origin_offset.push_back(tri.origin_offset);
            materials[i] = tri.material;
        }
        p1 = m256_vec3(_p1);
        p1p2 = m256_vec3(_p1p2);
        p1p3 = m256_vec3(_p1p3);
        for (int i = 0; i < tris.size(); i++) {
            origin_offset[i] = _origin_offset[i];
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
    Tri(XYZ _p1, XYZ _p2, XYZ _p3, XY _UV_1, XY _UV_2, XY _UV_3, short _material) : Primitive(_material),
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
    Tri(XYZ _p1, XYZ _p2, XYZ _p3, short _material) : Primitive(_material),
        p1(_p1), p2(_p2), p3(_p3) {
        midpoint = (p1 + p2 + p3) / 3;
        AABB_max = XYZ::max(p1, XYZ::max(p2, p3));
        AABB_min = XYZ::min(p1, XYZ::min(p2, p3));
    }
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
    static XYZ intersection_check_MoTru(const PackagedTri& T, const XYZ& position, const XYZ& slope) {
        XYZ pvec = XYZ::cross(slope, T.p1p3);
         float det = XYZ::dot(T.p1p2, pvec);
#if BACKFACE_CULLING
        // if the determinant is negative, the triangle is 'back facing'
        // if the determinant is close to 0, the ray misses the triangle
        (det < kEpsilon) return XYZ(-1);
#else
        // ray and triangle are parallel if det is close to 0
        if (fabs(det) < kEpsilon) return -1;
#endif
         float invDet = 1 / det;

        XYZ tvec = position - T.p1;
         float u = XYZ::dot(tvec, pvec) * invDet;
        if (u < 0 || u > 1) return XYZ(-1);

        XYZ qvec = XYZ::cross(tvec, T.p1p2);
         float v = XYZ::dot(slope, qvec) * invDet;
        if (v < 0 || u + v > 1) return XYZ(-1);

         float t = XYZ::dot(T.p1p3, qvec) * invDet;
         float U = u * T.U_delta.X + v * T.V_delta.X + T.UV_basis.X;
         float V = u * T.U_delta.Y + v * T.V_delta.Y + T.UV_basis.Y;
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
    void prep(vector<Material>& materials) {
        for (Tri& tri : tris) {
            materials[tri.material].prep();
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

    string name = "";

    Object(XYZ _origin, XYZ _scale) :
        origin(_origin), scale(_scale) {

    }
    void addMesh(Mesh* M) {
        meshes.push_back(M);
    }
    void addChild(Object* object) {
        children.push_back(object);
        object->parent = this;
    }
    void prep(vector<Material>& materials) {
        for (Mesh* m : meshes) {
            m->prep(materials);
        }
    }
    void fetchData(vector<Tri>& tris_vec) {
        XYZ scale_final = compile_scale();
        Transformation T = final_transform();
        for (Mesh* M : meshes) {
            for (Tri& tri : M->tris) {
                tris_vec.push_back(Packers::transformT(tri, T, origin, scale_final));
            }
        }
        for (Object* child : children) {
            child->fetchData(tris_vec);
        }
    }
    void _registerRotation(Quat rot) {
        transformations.push_back(Transformation(rot));
    }
    void _registerMove(XYZ move) {
        transformations.push_back(Transformation(move));
    }
    void applyTransformXYZ(float x, float y, float z) {
        _registerMove(XYZ(x, y, z));
    }
    void applyTransformXYZ(XYZ xyz) {
        applyTransformXYZ(xyz.X, xyz.Y, xyz.Z);
    }
    void rotateX(float rotation) {
        XYZ orig = XYZ(0, 1, 0);
        XYZ pointing = XYZ(0, cos(rotation), sin(rotation));
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    void rotateY(float rotation) {
        XYZ orig = XYZ(0, 0, 1);
        XYZ pointing = XYZ(sin(rotation), 0, cos(rotation));
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    void rotateZ(float rotation) {
        XYZ orig = XYZ(1, 0, 0);
        XYZ pointing = XYZ(cos(rotation), sin(rotation), 0);
        Quat rot = Quat::makeRotation(orig, pointing);
        _registerRotation(rot);
    }
    XYZ compile_scale(XYZ _scale = XYZ(1, 1, 1)) {
        _scale *= scale;
        if (parent != nullptr) {
            _scale *= parent->compile_scale(_scale);
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