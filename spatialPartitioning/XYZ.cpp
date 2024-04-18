
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <immintrin.h>
#include <math.h>
#include <unordered_set>

#include "XYZ.h"



template <typename T> T sign(T& input) {
    return (T)((input < 0) ? -1 : 1);
}
XY::XY() :X(0), Y(0) {}
XY::XY(float v) : X(v), Y(v) {}
XY::XY(float _x, float _y) : X(_x), Y(_y) {}
XY XY::operator+(const XY& other) const {
    return XY(X + other.X, Y + other.Y);
}
XY XY::operator-(const XY& other) const {
    return XY(X - other.X, Y - other.Y);
}
XY XY::operator*(const float& scalar) const {
    return XY(X * scalar, Y * scalar);
}
XY XY::operator*(const XY& other) const {
    return XY(X * other.X, Y * other.Y);
}
bool XY::operator==(const XY& other) const
{
    if (X == other.X && Y == other.Y) return true;
    else return false;
}

struct HashFunction
{
    size_t operator()(const XY& coord) const
    {
        size_t xHash = std::hash<int>()(coord.X);
        size_t yHash = std::hash<int>()(coord.Y) << 1;
        return xHash ^ yHash;
    }
};
typedef std::unordered_set<XY, XY::HashFunction> set;
std::string XY::to_string() const {
    return "(" + std::to_string(X) + "," + std::to_string(Y) + ")";
}

std::ostream& operator<<(std::ostream& stream, const XY& vec2) {
stream << vec2.to_string();
return stream;
}

XY operator*(const float& self, const XY& coord) {
return coord * self;
}

m256_vec2::m256_vec2() {};
m256_vec2::m256_vec2(XY fill) {
    X = _mm256_set1_ps(fill.X);
    Y = _mm256_set1_ps(fill.Y);
}
m256_vec2::m256_vec2(std::vector<XY> input) {
    std::vector<float> x_vec;
    std::vector<float> y_vec;
    for (int i = 0; i < 8 && i < input.size(); i++) {
        x_vec.push_back(input[i].X);
        y_vec.push_back(input[i].Y);
    }
    while (x_vec.size() < 8) {
        x_vec.push_back(0);
        y_vec.push_back(0);
    }
    X = _mm256_load_ps(x_vec.data());
    Y = _mm256_load_ps(y_vec.data());;
}
XY m256_vec2::at(int i) const {
    double x = ((float*)&X)[i];
    double y = ((float*)&Y)[i];
    return XY(x, y);
}


XYZ::XYZ() : X(0), Y(0), Z(0) {}
XYZ::XYZ(float s) : X(s), Y(s), Z(s) {}
XYZ::XYZ(float _X, float _Y, float _Z) : X(_X), Y(_Y), Z(_Z) {}
XYZ::XYZ(std::vector<double> attr) : X(attr[0]), Y(attr[1]), Z(attr[2]) {}
XYZ::XYZ(std::vector<float> attr) : X(attr[0]), Y(attr[1]), Z(attr[2]) {}
XYZ::XYZ(std::vector<double> attr, std::vector<char> swizzle) : X(attr[sign(swizzle[0]) * abs(swizzle[0])]), Y(sign(swizzle[1])* attr[abs(swizzle[1])]), Z(sign(swizzle[2])* attr[abs(swizzle[2])]) {}
XYZ::XYZ(std::vector<float> attr, std::vector<char> swizzle) : X(attr[sign(swizzle[0]) * abs(swizzle[0])]), Y(sign(swizzle[1])* attr[abs(swizzle[1])]), Z(sign(swizzle[2])* attr[abs(swizzle[2])]) {}
XYZ::XYZ(XY xy, float _Z) : X(xy.X), Y(xy.Y), Z(_Z) {}
XYZ::XYZ(XY xy) : XYZ(xy, 0) {}
XYZ XYZ::clone() {
    return XYZ(X, Y, Z);
}
float XYZ::operator[](const int n) const {
    if (n == 0) return X;
    if (n == 1) return Y;
    if (n == 2) return Z;
    //return *((float*)this + n);
}
float& XYZ::operator[](const int n) {
    if (n == 0) return X;
    if (n == 1) return Y;
    if (n == 2) return Z;
}
XYZ XYZ::swizzle(std::vector<char> swiz) {
    return XYZ(
        sign(swiz[0]) * (*this)[abs(swiz[0])],
        sign(swiz[1]) * (*this)[abs(swiz[1])],
        sign(swiz[2]) * (*this)[abs(swiz[2])]
    );
}
void XYZ::add(float addend) {
    X += addend;
    Y += addend;
    Z += addend;
}
void XYZ::add(const XYZ& other) {
    X += other.X;
    Y += other.Y;
    Z += other.Z;
}
void XYZ::divide(float divisor) {
    X /= divisor;
    Y /= divisor;
    Z /= divisor;
}
void XYZ::clip_negative(float clip_to) {
    X = (X < 0) ? clip_to : X;
    Y = (Y < 0) ? clip_to : Y;
    Z = (Z < 0) ? clip_to : Z;
}
float XYZ::smallest() {
    if (X < Y) {
        if (X < Z) {
            return X;
        }
        else {
            return Z;
        }
    }
    else {
        if (Y < Z) {
            return Y;
        }
        else {
            return Z;
        }
    }
}
void XYZ::make_safe() { //ensures point can always be used for division
    X = (X == 0) ? (float)0.000001 : X;
    Y = (Y == 0) ? (float)0.000001 : Y;
    Z = (Z == 0) ? (float)0.000001 : Z;
}
float XYZ::magnitude() {
    return sqrt(X * X + Y * Y + Z * Z); // 
}
float XYZ::magnitude_noRT() {
    return X * X + Y * Y + Z * Z; // 
}
void XYZ::normalize() {
    float len = magnitude();
    divide(len);
}
float XYZ::magnitude(const XYZ& v) {
    return sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
}
XYZ XYZ::reflect(XYZ vector, XYZ pole) {
    float vector_magnitude = vector.magnitude();
    float pole_magnitude = pole.magnitude();
    float mag_ratio = vector_magnitude / pole_magnitude;
    XYZ reflection_point = pole * (XYZ::cosine(vector, pole) * mag_ratio);
    XYZ pointing = reflection_point - vector;
    //float component = 2*XYZ::cosine(vector, pole)-1; //Weird math shit.
    return vector + (pointing * 2);
}
XYZ XYZ::floor(XYZ point) {
    return XYZ(
        (int)point.X, //using int type coercion since floor was putting an error for some reason
        (int)point.Y,
        (int)point.Z
    );
}
float XYZ::maxComponent(const XYZ& point) {
    return std::max(point.X, std::max(point.Y, point.Z));
}
float XYZ::minComponent(const XYZ& point) {
    return std::min(point.X, std::min(point.Y, point.Z));
}
XYZ XYZ::min(const XYZ& point, float value) {
    return XYZ(
        (point.X < value) ? point.X : value,
        (point.Y < value) ? point.Y : value,
        (point.Z < value) ? point.Z : value
    );

}
XYZ XYZ::min(float value, XYZ point) {
    return XYZ::min(point, value);
}
XYZ XYZ::min(const XYZ& point, const XYZ& other) {
    return XYZ(
        (point.X < other.X) ? point.X : other.X,
        (point.Y < other.Y) ? point.Y : other.Y,
        (point.Z < other.Z) ? point.Z : other.Z
    );
}
XYZ XYZ::max(const XYZ& point, const XYZ& other) {
    return XYZ(
        (point.X > other.X) ? point.X : other.X,
        (point.Y > other.Y) ? point.Y : other.Y,
        (point.Z > other.Z) ? point.Z : other.Z
    );
}
XYZ XYZ::max(const XYZ& point, const float& value) {
    return XYZ(
        (point.X > value) ? point.X : value,
        (point.Y > value) ? point.Y : value,
        (point.Z > value) ? point.Z : value
    );

}
float XYZ::length(const XYZ& vector) {
    return sqrt(vector.X * vector.X + vector.Y * vector.Y + vector.Z * vector.Z);
}

    XYZ XYZ::normalize(const XYZ& vector) {
    float len = XYZ::length(vector);
    return XYZ(vector.X / len, vector.Y / len, vector.Z / len);
}
    XYZ XYZ::divide(const XYZ& point, const XYZ& other) {
    return XYZ(point.X / other.X, point.Y / other.Y, point.Z / other.Z);
}
    XYZ XYZ::divide(const XYZ& point, float divisor) {
    return XYZ(point.X / divisor, point.Y / divisor, point.Z / divisor);
}
    XYZ XYZ::add(const XYZ& point, float addend) {
    return XYZ(point.X + addend, point.Y + addend, point.Z + addend);
}
    XYZ XYZ::add(const XYZ& point, const XYZ& other) {
    return XYZ(point.X + other.X, point.Y + other.Y, point.Z + other.Z);
}
    float XYZ::distance_noRt(XYZ& point, XYZ& other) {
    float f1 = point.X - other.X;
    float f2 = point.Y - other.Y;
    float f3 = point.Z - other.Z;
    return f1 * f1 + f2 * f2 + f3 * f3;
}
    float XYZ::distance(const XYZ& point, const XYZ& other) {
    float f1 = point.X - other.X;
    float f2 = point.Y - other.Y;
    float f3 = point.Z - other.Z;
    return sqrt(f1 * f1 + f2 * f2 + f3 * f3);
}
    XYZ XYZ::_rtslope(const XYZ& point, const XYZ& other) {
    float distance = XYZ::distance(point, other);
    return (other - point) / distance;
}
    XYZ XYZ::_dotslope(const XYZ& point, const XYZ& other) {
    XYZ delta = other - point;
    return delta / XYZ::dot(delta, delta);
}
    XYZ XYZ::slope(const XYZ& point, const XYZ& other) {
    return _rtslope(point, other);
}
    XYZ XYZ::flip(XYZ point) {
    return XYZ(-point.X, -point.Y, -point.Z);
}
    float XYZ::dot(XYZ point, XYZ other) {
    return point.X * other.X + point.Y * other.Y + point.Z * other.Z;
}
    float XYZ::cdot(XYZ point, XYZ other) {
        return std::min(1.0f,std::max(0.0f,point.X * other.X + point.Y * other.Y + point.Z * other.Z));
    }
    float XYZ::cosine(XYZ point, XYZ other) {
    float result = XYZ::dot(point, other) / (point.magnitude() * other.magnitude());
    return result;
}
    XYZ XYZ::pow(XYZ point, float power) {
    return XYZ(std::pow(point.X, power), std::pow(point.Y, power), std::pow(point.Z, power));
}
    XYZ XYZ::cross(XYZ point, XYZ other) {
    return XYZ(
        point.Y * other.Z - point.Z * other.Y,
        point.Z * other.X - point.X * other.Z,
        point.X * other.Y - point.Y * other.X
    );
}
    XYZ XYZ::log(XYZ point) {
    return XYZ(
        log10(point.X),
        log10(point.Y),
        log10(point.Z)
    );
}
    XYZ XYZ::clamp(XYZ value, XYZ low, XYZ high) {
    return XYZ(
        std::min(std::max(value.X, low.X), high.X),
        std::min(std::max(value.Y, low.Y), high.Y),
        std::min(std::max(value.Z, low.Z), high.Z)
    );
}
    XYZ XYZ::clamp(XYZ value, float low, float high) {
    return clamp(value, XYZ(low, low, low), XYZ(high, high, high));
}
    XYZ XYZ::negative(const XYZ& v) {
    return XYZ(-v.X, -v.Y, -v.Z);
}

std::string XYZ::to_string() const {
    std::string xstr = std::to_string(X);
    xstr = ((X > 0) ? " " : "") + xstr;
    std::string ystr = std::to_string(Y);
    ystr = ((Y > 0) ? " " : "") + ystr;
    std::string zstr = std::to_string(Z);
    zstr = ((Z > 0) ? " " : "") + zstr;
    return "(" + xstr + "," + ystr + "," + zstr + ")";
}
XYZ XYZ::operator/(const XYZ& other) const {
    return XYZ(X / other.X, Y / other.Y, Z / other.Z);
}
XYZ XYZ::operator/(const float divisor) const {
    return XYZ(X / divisor, Y / divisor, Z / divisor);
}
XYZ XYZ::operator+(const XYZ& other) const {
    return XYZ(X + other.X, Y + other.Y, Z + other.Z);
}
XYZ XYZ::operator+(const float addend) const {
    return XYZ(X + addend, Y + addend, Z + addend);
}
XYZ XYZ::operator-(const XYZ& other) const {
    return XYZ(X - other.X, Y - other.Y, Z - other.Z);
}
XYZ XYZ::operator-(const float addend) const {
    return XYZ(X - addend, Y - addend, Z - addend);
}
XYZ XYZ::operator*(const float multiplier) const {
    return XYZ(X * multiplier, Y * multiplier, Z * multiplier);
}
XYZ XYZ::operator*(const XYZ& other) const {
    return XYZ(X * other.X, Y * other.Y, Z * other.Z);
}
XYZ XYZ::operator+=(const XYZ& other) {
    X += other.X;
    Y += other.Y;
    Z += other.Z;
    return *this;
}
XYZ XYZ::operator *= (const XYZ& other) {
    X *= other.X;
    Y *= other.Y;
    Z *= other.Z;
    return *this;
}
XYZ XYZ::operator *= (const float& scalar) {
    X *= scalar;
    Y *= scalar;
    Z *= scalar;
    return *this;
}

XYZ XYZ::operator-() {
    return XYZ(-X, -Y, -Z);
}
bool XYZ::operator !=(const XYZ& other) {
    return !((X == other.X) && (Y == other.Y) && (Z == other.Z));
}
bool XYZ::operator !=(XYZ& other) {
    return !((X == other.X) && (Y == other.Y) && (Z == other.Z));
}
    bool XYZ::equals(const XYZ& point, const XYZ& other) {
    return (point.X == other.X) && (point.Y == other.Y) && (point.Z == other.Z);
}
bool XYZ::operator ==(const XYZ& other) const {
    return (X == other.X) && (Y == other.Y) && (Z == other.Z);
}
XYZ XYZ::linear_mix(float c, const XYZ& first, const XYZ& second) {
    float i = 1 - c;
    return XYZ(first.X * i + second.X * c, first.Y * i + second.Y * c, first.Z * i + second.Z * c);
}
inline bool XYZ::less_than_x_operator::operator() (const XYZ& point1, const XYZ& point2)
{
    return (point1.X < point2.X);
}
inline bool XYZ::less_than_y_operator::operator() (const XYZ& point1, const XYZ& point2)
{
    return (point1.Y < point2.Y);
}
inline bool XYZ::less_than_z_operator::operator() (const XYZ& point1, const XYZ& point2)
{
    return (point1.Z < point2.Z);
}
inline bool XYZ::less_than_x_operator_p::operator() (const XYZ* point1, const XYZ* point2)
{
    return (point1->X < point2->X);
}
inline bool XYZ::less_than_y_operator_p::operator() (const XYZ* point1, const XYZ* point2)
{
    return (point1->Y < point2->Y);
}

inline bool XYZ::less_than_z_operator_p::operator() (const XYZ* point1, const XYZ* point2)
{
    return (point1->Z < point2->Z);
}


XYZ operator*(const float& self, const XYZ& point) {
    return point * self;
}
XYZ operator-(const float& self, const XYZ& point) {
    return point + self;
}
XYZ operator/(const float& self, const XYZ& point) {
    return XYZ::divide(self, point);
}

std::ostream& operator<<(std::ostream& os, XYZ& m) {
    return os << m.to_string();
}

float X_t::v() {
    return X;
}

float Y_t::v() {
    return Y;
}

float Z_t::v() {
    return Z;
}


m256_vec3::m256_vec3() {};
m256_vec3::m256_vec3(XYZ fill) {
    X = _mm256_set1_ps(fill.X);
    Y = _mm256_set1_ps(fill.Y);
    Z = _mm256_set1_ps(fill.Z);
}
m256_vec3::m256_vec3(std::vector<XYZ> input) {
    std::vector<float> x_vec;
    std::vector<float> y_vec;
    std::vector<float> z_vec;
    for (int i = 0; i < 8 && i < input.size(); i++) {
        x_vec.push_back(input[i].X);
        y_vec.push_back(input[i].Y);
        z_vec.push_back(input[i].Z);
    }
    while (x_vec.size() < 8) {
        x_vec.push_back(0);
        y_vec.push_back(0);
        z_vec.push_back(0);
    }
    X = _mm256_load_ps(x_vec.data());
    Y = _mm256_load_ps(y_vec.data());;
    Z = _mm256_load_ps(z_vec.data());;
}
XYZ m256_vec3::at(int i) const {
    double x = ((float*)&X)[i];
    double y = ((float*)&Y)[i];
    double z = ((float*)&Z)[i];
    return XYZ(x, y, z);
}
void m256_vec3::sub(const m256_vec3& v1, const m256_vec3& v2, m256_vec3& output) {
    output.X = _mm256_sub_ps(v1.X, v2.X);
    output.Y = _mm256_sub_ps(v1.Y, v2.Y);
    output.Z = _mm256_sub_ps(v1.Z, v2.Z);
}
m256_vec3 m256_vec3::sub_inline(const m256_vec3& v1, const m256_vec3& v2) {
    m256_vec3 output;
    sub(v1, v2, output);
    return output;
}

    Quat::Quat() : XYZ(), W(0) {}
    Quat::Quat(XYZ _XYZ, float _W) : XYZ(_XYZ), W(_W) {}
    Quat::Quat(float _X, float _Y, float _Z) : XYZ(_X, _Y, _Z), W(0) {}
    Quat::Quat(float _X, float _Y, float _Z, float _W) : XYZ(_X, _Y, _Z), W(_W) {}
    Quat::Quat(std::vector<float> attr) : XYZ(attr), W(attr[3]) {}
    Quat::Quat(std::vector<double> attr) : XYZ(attr), W(attr[3]) {}
    Quat::Quat(std::vector<float> attr, std::vector<char> swizzle) : XYZ(attr, swizzle), W(sign(swizzle[3])* attr[abs(swizzle[3])]) {}
    Quat::Quat(std::vector<double> attr, std::vector<char> swizzle) : XYZ(attr, swizzle), W(sign(swizzle[3])* attr[abs(swizzle[3])]) {}
    Quat Quat::clone() {
        return Quat(X, Y, Z, W);
    }
    float Quat::magnitude() {
        return Quat::dot(*this, *this);//inefficient but whatever
    }
    float Quat::dot(const Quat& q1, const Quat& q2) {
        return q1.X * q2.X + q1.Y * q2.Y + q1.Z * q2.Z + q1.W * q2.W;
    }
    Quat Quat::operator*(const float& scalar) const {
        return Quat(X * scalar, Y * scalar, Z * scalar, W * scalar);
    }
    Quat Quat::multiply(const Quat& q2, const Quat& q1) {
        return Quat(
            q1.W * q2.X + q1.X * q2.W + q1.Y * q2.Z - q1.Z * q2.Y,
            q1.W * q2.Y - q1.X * q2.Z + q1.Y * q2.W + q1.Z * q2.X,
            q1.W * q2.Z + q1.X * q2.Y - q1.Y * q2.X + q1.Z * q2.W,
            q1.W * q2.W - q1.X * q2.X - q1.Y * q2.Y - q1.Z * q2.Z
        );
    }
    Quat Quat::normalize(const Quat& in) {
        float d = Quat::dot(in, in);
        float l = sqrt(d);
        return in / l;
    }

    Quat Quat::makeRotation(const XYZ& up, const XYZ& direction) {
        if (XYZ::equals(direction, up)) {
            return Quat(0, 0, 0, 1);
        }
        XYZ a = XYZ::cross(up, direction);
        if (XYZ::equals(direction, XYZ::negative(up))) {
            return Quat(a, 0);
        }

        float m1 = XYZ::magnitude(up);
        float m2 = XYZ::magnitude(direction);

        Quat out = Quat(
            a,
            sqrt(m1 * m1 * m2 * m2) + XYZ::dot(up, direction)
        );

        return Quat::normalize(out);
    }
    Quat Quat::makeRotationFromY(const XYZ& direction) {
        if (XYZ::equals(direction, XYZ(0, 1, 0))) {
            return Quat(0, 0, 0, 1);
        }
        if (XYZ::equals(direction, XYZ(0, -1, 0))) {
            return Quat(0, 0, 1, 0);
        }
        float m1 = 1;
        float m2 = XYZ::magnitude(direction);
        Quat out = Quat(
            XYZ(
                direction.Z,
                0,
                -direction.X
            ),
            m2 + direction.Y
        );
        return Quat::normalize(out);

    }
    XYZ Quat::applyRotation(const XYZ& p, const Quat& r) { //https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
        //I dont understand ANY of this, but by god apparently its fast
        XYZ t = 2 * XYZ::cross(r, p);
        return p + r.W * t + cross(r, t);
    }
    Quat Quat::operator/(float div) const {
        return Quat(
            X / div,
            Y / div,
            Z / div,
            W / div
        );
    }


    iXY::iXY() :X(0), Y(0) {}
    iXY::iXY(int _x, int _y) : X(_x), Y(_y) {}
    iXY iXY::operator+(const iXY & other) const {
        return iXY(X + other.X, Y + other.Y);
    }
    iXY iXY::operator-(const iXY & other) const {
        return iXY(X - other.X, Y - other.Y);
    }
    bool iXY::operator==(const iXY & other) const
    {
        if (X == other.X && Y == other.Y) return true;
        else return false;
    }

    struct iXY::HashFunction
    {
        size_t operator()(const iXY& coord) const
        {
            size_t xHash = std::hash<int>()(coord.X);
            size_t yHash = std::hash<int>()(coord.Y) << 1;
            return xHash ^ yHash;
        }
    };
