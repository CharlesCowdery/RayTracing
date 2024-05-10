#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <immintrin.h>
#include <math.h>
#include <unordered_set>

#define NOMINMAX 1

#undef max
#undef min

template <typename T> T sign(T& input);
struct XY;
std::ostream& operator<<(std::ostream& stream, const XY& vec2);
XY operator*(const float& self, const XY& coord);
struct m256_vec2 {
    __m256 X;
    __m256 Y;
    m256_vec2();
    m256_vec2(XY fill);
    m256_vec2(std::vector<XY> input);
    XY at(int i) const;
    void set(int i, XY v);
};

struct XY {
    float X;
    float Y;
    XY();
    XY(float v);
    XY(float _x, float _y);
    XY operator+(const XY& other) const;
    XY operator-(const XY& other) const;
    XY operator*(const float& scalar) const;
    XY operator*(const XY& other) const;
    bool operator==(const XY& other) const;
    struct HashFunction;
    typedef std::unordered_set<XY, XY::HashFunction> set;
    std::string to_string() const;
};

struct XYZ;
std::ostream& operator<<(std::ostream& os, XYZ& m);


struct XYZ {
public:
    float X;
    float Y;
    float Z;
    XYZ();
    XYZ(float s);
    XYZ(float _X, float _Y, float _Z);
    XYZ(std::vector<double> attr);
    XYZ(std::vector<float> attr);
    XYZ(std::vector<double> attr, std::vector<char> swizzle);
    XYZ(std::vector<float> attr, std::vector<char> swizzle);
    XYZ(XY xy, float _Z);
    XYZ(XY xy);
    XYZ clone();
    float operator[](const int n) const;
    float& operator[](const int n);
    XYZ swizzle(std::vector<char> swiz);
    void add(float addend);
    void add(const XYZ& other);
    void divide(float divisor);
    void clip_negative(float clip_to = 0);
    float smallest();
    void make_safe();
    float magnitude();
    float magnitude_noRT();
    void normalize();
    bool wellDefined() const;
    std::string to_string() const;
    static float magnitude(const XYZ& v);
    static XYZ reflect(XYZ vector, XYZ pole);
    static XYZ floor(XYZ point);
    static float maxComponent(const XYZ& point);
    static float minComponent(const XYZ& point);
    static XYZ min(const XYZ& point, float value);
    static XYZ min(float value, XYZ point);
    static XYZ min(const XYZ& point, const XYZ& other);
    static XYZ max(const XYZ& point, const XYZ& other);
    static XYZ max(const XYZ& point, const float& value);
    static float length(const XYZ& vector);
    static XYZ normalize(const XYZ& vector);
    static XYZ divide(const XYZ& point, const XYZ& other);
    static XYZ divide(const XYZ& point, float divisor);
    static XYZ add(const XYZ& point, float addend);
    static XYZ add(const XYZ& point, const XYZ& other);
    static float distance_noRt(XYZ& point, XYZ& other);
    static float distance(const XYZ& point, const XYZ& other);
    static XYZ _rtslope(const XYZ& point, const XYZ& other);
    static XYZ _dotslope(const XYZ& point, const XYZ& other);
    static XYZ slope(const XYZ& point, const XYZ& other);
    static XYZ flip(XYZ point);
    static float dot(XYZ point, XYZ other);
    static float cdot(XYZ point, XYZ other);
    static float cosine(XYZ point, XYZ other);
    static XYZ pow(XYZ point, float power);
    static XYZ cross(XYZ point, XYZ other);
    static XYZ log(XYZ point);
    static XYZ clamp(XYZ value, XYZ low, XYZ high);
    static XYZ clamp(XYZ value, float low, float high);
    static XYZ negative(const XYZ& v);
    XYZ operator/(const XYZ& other) const;
    XYZ operator/(const float divisor) const;
    XYZ operator+(const XYZ& other) const;
    XYZ operator+(const float addend) const;
    XYZ operator-(const XYZ& other) const;
    XYZ operator-(const float addend) const;
    XYZ operator*(const float multiplier) const;
    XYZ operator*(const XYZ& other) const;
    XYZ operator+=(const XYZ& other);
    XYZ operator *= (const XYZ& other);
    XYZ operator *= (const float& scalar);
    XYZ operator-();
    bool operator !=(const XYZ& other);
    bool operator !=(XYZ& other);
    static bool equals(const XYZ& point, const XYZ& other);
    bool operator ==(const XYZ& other) const;
    static XYZ linear_mix(float c, const XYZ& first, const XYZ& second);
    struct less_than_x_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2);
    };
    struct less_than_y_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2);
    };
    struct less_than_z_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2);
    };
    struct less_than_x_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2);
    };
    struct less_than_y_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2);
    };
    struct less_than_z_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2);
    };

};

struct X_t : public XYZ {
    float v();
};

struct Y_t : public XYZ {
    float v();
};

struct Z_t : public XYZ {
    float v();
};

XYZ operator*(const float& self, const XYZ& point);
XYZ operator-(const float& self, const XYZ& point);
XYZ operator/(const float& self, const XYZ& point);

struct m256_vec3;
struct Quat;

struct m256_vec3 {
    __m256 X;
    __m256 Y;
    __m256 Z;
    m256_vec3();
    m256_vec3(XYZ fill);
    m256_vec3(std::vector<XYZ> input);
    XYZ at(int i) const;
    void set(int i, XYZ v);
    void transfer(int origin, int target);
    static void sub(const m256_vec3& v1, const m256_vec3& v2, m256_vec3& output);
    static m256_vec3 sub(const m256_vec3& v1, const m256_vec3& v2);
    static m256_vec3 sub_inline(const m256_vec3& v1, const m256_vec3& v2);
    static m256_vec3 mul(const m256_vec3& v1, const m256_vec3& v2);
    static m256_vec3 mul(const m256_vec3& v1, const __m256& v2);
    static m256_vec3 muladd(const m256_vec3& v1, const __m256& v2, const  m256_vec3& v3);
    static m256_vec3 bitwise_xor(const m256_vec3& v1, const __m256& v2);
    static m256_vec3 bitwise_or(const m256_vec3& v1, const  m256_vec3& v2);
    static m256_vec3 bitwise_and(const m256_vec3& v1, const __m256& v2);
    static __m256 dot(const m256_vec3& v1, const m256_vec3& v2);
    static m256_vec3 normalize(const m256_vec3& v);
    static m256_vec3 cross(const m256_vec3& v1, const m256_vec3& v2);
    static m256_vec3 apply_matrix(const m256_vec3 v, const m256_vec3 m1, const m256_vec3 m2, const m256_vec3 m3);
};
struct Quat : public XYZ {
    float W;
    Quat();
    Quat(XYZ _XYZ, float _W);
    Quat(float _X, float _Y, float _Z);
    Quat(float _X, float _Y, float _Z, float _W);
    Quat(std::vector<float> attr);
    Quat(std::vector<double> attr);
    Quat(std::vector<float> attr, std::vector<char> swizzle);
    Quat(std::vector<double> attr, std::vector<char> swizzle);
    Quat clone();
    float magnitude();
    static float dot(const Quat& q1, const Quat& q2);
    Quat operator*(const float& scalar) const;
    static Quat multiply(const Quat& q2, const Quat& q1);
    static Quat normalize(const Quat& in);
    static Quat makeRotation_(const XYZ& start, const XYZ& end);
    static Quat makeRotation(const XYZ& up, const XYZ& direction);
    static Quat makeRotationFromY(const XYZ& direction);
    static XYZ applyRotation(const XYZ& p, const Quat& r);
    Quat operator/(float div) const;
};

struct iXY {
    int X;
    int Y;
    iXY();
    iXY(int _x, int _y);
    iXY operator+(const iXY& other) const;
    iXY operator-(const iXY& other) const;
    bool operator==(const iXY& other) const;
    struct HashFunction;
    typedef std::unordered_set<iXY, iXY::HashFunction> set;
};
