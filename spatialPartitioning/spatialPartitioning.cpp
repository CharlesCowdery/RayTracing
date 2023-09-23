﻿// spatialPartitioning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>
#include <format>
#include <string>
#include <ctime>
#include <chrono>
#include <SFML/Graphics.hpp>
#include <functional>
#include <fstream>
#include <array>
#include <unordered_set>
#include <thread>

#include <cassert>
#include <intrin.h>

#define assertm(exp, msg) assert(((void)msg, exp))
#define SMALL 0.001
#define NEAR_THRESHOLD 0.001
#define SCENE_BOUNDS 100
#define RED_MASK 255<<16
#define GREEN_MASK 255<<8
#define BLUE_MASK 255

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define PI 3.1415

#define kEpsilon 0.0001

#define FULL_BRIGHT false

#define OBJECT_CULLING true
#define CULL_RECEDING_OBJECTS true
#define LAZY_OBJECT_CHECKING true
//defines program precision
#define decimal float

using namespace std;

/*Xorshiro256+ pseudorandom start. Not my code: https://prng.di.unimi.it/*/

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}


static uint64_t xe_seed[4];

uint64_t xe_next(void) {
    const uint64_t result = xe_seed[0] + xe_seed[3];

    const uint64_t t = xe_seed[1] << 17;

    xe_seed[2] ^= xe_seed[0];
    xe_seed[3] ^= xe_seed[1];
    xe_seed[1] ^= xe_seed[2];
    xe_seed[0] ^= xe_seed[3];

    xe_seed[2] ^= t;

    xe_seed[3] = rotl(xe_seed[3], 45);

    return result;
}

void xe_jump(void) {
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for (int i = 0; i < sizeof JUMP / sizeof * JUMP; i++)
        for (int b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= xe_seed[0];
                s1 ^= xe_seed[1];
                s2 ^= xe_seed[2];
                s3 ^= xe_seed[3];
            }
            xe_next();
        }

    xe_seed[0] = s0;
    xe_seed[1] = s1;
    xe_seed[2] = s2;
    xe_seed[3] = s3;
}

void xe_long_jump(void) {
    static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for (int i = 0; i < sizeof LONG_JUMP / sizeof * LONG_JUMP; i++)
        for (int b = 0; b < 64; b++) {
            if (LONG_JUMP[i] & UINT64_C(1) << b) {
                s0 ^= xe_seed[0];
                s1 ^= xe_seed[1];
                s2 ^= xe_seed[2];
                s3 ^= xe_seed[3];
            }
            xe_next();
        }

    xe_seed[0] = s0;
    xe_seed[1] = s1;
    xe_seed[2] = s2;
    xe_seed[3] = s3;
}
/*Xorshiro256+ pseudorandom end*/

decimal xe_frand() {
    return (xe_next() >> 11) * 0x1.0p-53;
}

decimal rand_frand() //https://stackoverflow.com/a/2704552
{
    return (decimal)rand() / RAND_MAX;
}

decimal fRand(decimal fMin, decimal fMax) //https://stackoverflow.com/a/2704552
{
    decimal f = xe_frand();
    return fMin + f * (fMax - fMin);
}

decimal clamp(decimal value, decimal low, decimal high) {
    return max(min(value, high), low);
}

decimal lerp(decimal a, decimal b, decimal factor) {
    return a + factor * (b - a);
}

std::vector<std::string> readLinesFromFile(const std::string& fileName) { //chatgpt lol
    std::vector<std::string> lines;
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return lines;
    }

    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    file.close();
    return lines;
}

vector<string> eng_levels = vector<string>{
    "k",
    "M",
    "B",
    "T",
    "qu",
    "Qu"
};


string iToFixedLength(int input,int length, string fill = " ") {
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
        return iToFixedLength(inum,4)+" ";
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

struct iXY {
    int X;
    int Y;
    iXY() :X(0), Y(0) {}
    iXY(int _x, int _y) : X(_x), Y(_y) {}
    iXY operator+(const iXY& other) const {
        return iXY(X + other.X, Y + other.Y);
    }
    iXY operator-(const iXY& other) const {
        return iXY(X - other.X, Y - other.Y);
    }
    bool operator==(const iXY& other) const
    {
        if (X == other.X && Y == other.Y) return true;
        else return false;
    }

    struct HashFunction
    {
        size_t operator()(const iXY& coord) const
        {
            size_t xHash = std::hash<int>()(coord.X);
            size_t yHash = std::hash<int>()(coord.Y) << 1;
            return xHash ^ yHash;
        }
    };
    typedef unordered_set<iXY, iXY::HashFunction> set;
};

struct XY {
    decimal X;
    decimal Y;
    XY() :X(0), Y(0) {}
    XY(decimal _x, decimal _y) : X(_x), Y(_y) {}
    XY operator+(const XY& other) const {
        return XY(X + other.X, Y + other.Y);
    }
    XY operator-(const XY& other) const {
        return XY(X - other.X, Y - other.Y);
    }
    bool operator==(const XY& other) const
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
    typedef unordered_set<XY, XY::HashFunction> set;
};

struct XYZ {
public:
    decimal X;
    decimal Y;
    decimal Z;
    XYZ() : X(0), Y(0), Z(0) {}
    XYZ(decimal _X, decimal _Y, decimal _Z) : X(_X), Y(_Y), Z(_Z) {}
    XYZ(XY xy, decimal _Z) : X(xy.X), Y(xy.Y), Z(_Z) {}
    XYZ(XY xy) : XYZ(xy, 0) {}
    XYZ clone() {
        return XYZ(X, Y, Z);
    }
    void add(decimal addend) {
        X += addend;
        Y += addend;
        Z += addend;
    }
    void add(const XYZ& other) {
        X += other.X;
        Y += other.Y;
        Z += other.Z;
    }
    void divide(decimal divisor) {
        X /= divisor;
        Y /= divisor;
        Z /= divisor;
    }
    void clip_negative(decimal clip_to = 0) {
        X = (X < 0) ? clip_to : X;
        Y = (Y < 0) ? clip_to : Y;
        Z = (Z < 0) ? clip_to : Z;
    }
    decimal smallest() {
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
    void make_safe() { //ensures point can always be used for division
        X = (X == 0) ? (decimal)0.000001 : X;
        Y = (Y == 0) ? (decimal)0.000001 : Y;
        Z = (Z == 0) ? (decimal)0.000001 : Z;
    }
    decimal magnitude() {
        return sqrt(X * X + Y * Y + Z * Z); // 
    }
    decimal magnitude_noRT() {
        return X * X + Y * Y + Z * Z; // 
    }
    void normalize() {
        decimal len = magnitude();
        divide(len);
    }
    static decimal magnitude(const XYZ& v) {
        return sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
    }
    static XYZ reflect(XYZ vector, XYZ pole) {
        decimal vector_magnitude = vector.magnitude();
        decimal pole_magnitude = pole.magnitude();
        decimal mag_ratio = vector_magnitude / pole_magnitude;
        XYZ reflection_point = pole * (XYZ::cosine(vector, pole) * mag_ratio);
        XYZ pointing = reflection_point - vector;
        //decimal component = 2*XYZ::cosine(vector, pole)-1; //Weird math shit.
        return vector + (pointing * 2);
    }
    static XYZ floor(XYZ point) {
        return XYZ(
            std::floor(point.X),
            std::floor(point.Y),
            std::floor(point.Z)
        );
    }
    static XYZ min(XYZ point, decimal value) {
        return XYZ(
            (point.X < value) ? point.X : value,
            (point.Y < value) ? point.Y : value,
            (point.Z < value) ? point.Z : value
        );

    }
    static XYZ min(decimal value, XYZ point) {
        return XYZ::min(point, value);
    }
    static XYZ min(XYZ point, XYZ other) {
        return XYZ(
            (point.X < other.X) ? point.X : other.X,
            (point.Y < other.Y) ? point.Y : other.Y,
            (point.Z < other.Z) ? point.Z : other.Z
        );
    }
    static XYZ max(XYZ point, XYZ other) {
        return XYZ(
            (point.X > other.X) ? point.X : other.X,
            (point.Y > other.Y) ? point.Y : other.Y,
            (point.Z > other.Z) ? point.Z : other.Z
        );
    }
    static XYZ max(XYZ point, decimal value) {
        return XYZ(
            (point.X > value) ? point.X : value,
            (point.Y > value) ? point.Y : value,
            (point.Z > value) ? point.Z : value
        );

    }
    static decimal length(const XYZ& vector) {
        return sqrt(vector.X * vector.X + vector.Y * vector.Y + vector.Z * vector.Z);
    }

    static XYZ normalize(const XYZ& vector) {
        decimal len = XYZ::length(vector);
        return XYZ(vector.X / len, vector.Y / len, vector.Z / len);
    }
    static XYZ divide(const XYZ& point, const XYZ& other) {
        return XYZ(point.X / other.X, point.Y / other.Y, point.Z / other.Z);
    }
    static XYZ divide(const XYZ& point, decimal divisor) {
        return XYZ(point.X / divisor, point.Y / divisor, point.Z / divisor);
    }
    static XYZ add(const XYZ& point, decimal addend) {
        return XYZ(point.X + addend, point.Y + addend, point.Z + addend);
    }
    static XYZ add(const XYZ& point, const XYZ& other) {
        return XYZ(point.X + other.X, point.Y + other.Y, point.Z + other.Z);
    }
    static decimal distance_noRt(XYZ& point, XYZ& other) {
        decimal f1 = point.X - other.X;
        decimal f2 = point.Y - other.Y;
        decimal f3 = point.Z - other.Z;
        return f1 * f1 + f2 * f2 + f3 * f3;
    }
    static decimal distance(const XYZ& point,const XYZ& other){
        decimal f1 = point.X - other.X;
        decimal f2 = point.Y - other.Y;
        decimal f3 = point.Z - other.Z;
        return sqrt(f1 * f1 + f2 * f2 + f3 * f3);
    }
    static XYZ slope(XYZ point, XYZ other) {
        decimal distance = XYZ::distance(point, other);
        return (other - point) / distance;
    }
    static XYZ flip(XYZ point) {
        return XYZ(-point.X, -point.Y, -point.Z);
    }
    static decimal dot(XYZ point, XYZ other) {
        return point.X * other.X + point.Y * other.Y + point.Z * other.Z;
    }
    static decimal cosine(XYZ point, XYZ other) {
        decimal result = XYZ::dot(point, other) / (point.magnitude() * other.magnitude());
        return result;
    }
    static XYZ pow(XYZ point, decimal power) {
        return XYZ(std::pow(point.X, power), std::pow(point.Y, power), std::pow(point.Z, power));
    }
    static XYZ cross(XYZ point, XYZ other) {
        return XYZ(
            point.Y * other.Z - point.Z * other.Y,
            point.Z * other.X - point.X * other.Z,
            point.X * other.Y - point.Y * other.X
        );
    }
    static XYZ log(XYZ point) {
        return XYZ(
            log10(point.X),
            log10(point.Y),
            log10(point.Z)
        );
    }   
    static XYZ clamp(XYZ value, XYZ low, XYZ high) {
        return XYZ(
            std::min(std::max(value.X, low.X), high.X),
            std::min(std::max(value.Y, low.Y), high.Y),
            std::min(std::max(value.Z, low.Z), high.Z)
        );
    }
    static XYZ clamp(XYZ value, decimal low, decimal high) {
        return clamp(value, XYZ(low, low, low), XYZ(high, high, high));
    }
    static XYZ negative(const XYZ& v) {
        return XYZ(-v.X, -v.Y, -v.Z);
    }
    string to_string() {
        return std::to_string(X) + ", " + std::to_string(Y) + ", " + std::to_string(Z);
    }
    XYZ operator/(const XYZ& other) const {
        return XYZ(X / other.X, Y / other.Y, Z / other.Z);
    }
    XYZ operator/(const decimal divisor) const {
        return XYZ(X / divisor, Y / divisor, Z / divisor);
    }
    XYZ operator+(const XYZ& other) const {
        return XYZ(X + other.X, Y + other.Y, Z + other.Z);
    }
    XYZ operator+(const decimal addend) const {
        return XYZ(X + addend, Y + addend, Z + addend);
    }
    XYZ operator-(const XYZ& other) const {
        return XYZ(X - other.X, Y - other.Y, Z - other.Z);
    }
    XYZ operator-(const decimal addend) const {
        return XYZ(X - addend, Y - addend, Z - addend);
    }
    XYZ operator*(const decimal multiplier) const {
        return XYZ(X * multiplier, Y * multiplier, Z * multiplier);
    }
    XYZ operator*(const XYZ& other) const {
        return XYZ(X * other.X, Y * other.Y, Z * other.Z);
    }
    XYZ operator+=(const XYZ& other) {
        X += other.X;
        Y += other.Y;
        Z += other.Z;
        return *this;
    }
    XYZ operator-() {
        return XYZ(-X, -Y, -Z);
    }
    bool operator !=(const XYZ& other) {
        return !((X == other.X) && (Y == other.Y) && (Z == other.Z));
    }
    static bool equals(const XYZ& point, const XYZ& other) {
        return (point.X == other.X) && (point.Y == other.Y) && (point.Z == other.Z);
    }
    bool operator ==(const XYZ& other) {
        return (X == other.X) && (Y == other.Y) && (Z == other.Z);
    }
    static XYZ linear_mix(decimal c, const XYZ& first, const XYZ& second) {
        decimal i = 1 - c;
        return XYZ(first.X * i + second.X * c, first.Y * i + second.Y * c, first.Z * i + second.Z * c);
    }
};
XYZ operator*(decimal self, const XYZ& point) {
    return point * self;
}
XYZ operator-(decimal self, const XYZ& point) {
    return point + self;
}

std::ostream& operator<<(std::ostream& os, XYZ& m) {
    return os << m.to_string();
}

struct Quat : public XYZ {
    decimal W;
    Quat() : XYZ(), W(0) {}
    Quat(XYZ _XYZ, decimal _W) : XYZ(_XYZ), W(_W) {}
    Quat(decimal _X, decimal _Y, decimal _Z) : XYZ(_X,_Y,_Z), W(0) {}
    Quat(decimal _X, decimal _Y, decimal _Z,decimal _W) : XYZ(_X, _Y, _Z), W(_W) {}
    Quat clone() {
        return Quat(X, Y, Z, W);
    }
    static decimal dot(const Quat& q1, const Quat& q2) {
        return q1.X * q2.X + q1.Y * q2.Y + q1.Z * q2.Z + q1.W * q2.W;
    }
    static Quat normalize(const Quat& in) {
        decimal d = Quat::dot(in, in);
        decimal l = sqrt(d);
        return in / l;
    }
    static Quat makeRotation(const XYZ& up, const XYZ& direction) {
        if (XYZ::equals(direction, up)) {
            return Quat(0, 0, 0, 1);
        }
        if (XYZ::equals(direction, XYZ::negative(up))) {
            return Quat(direction, 0);
        }
        
        decimal m1 = XYZ::magnitude(up);
        decimal m2 = XYZ::magnitude(direction);

        Quat out = Quat(
            XYZ::cross(up, direction),
            sqrt(m1 * m1 * m2 * m2) + XYZ::dot(up, direction)
        );

        return Quat::normalize(out);
    }
    static Quat makeRotationFromY(const XYZ& direction) {
        if (XYZ::equals(direction, XYZ(0, 1, 0))) {
            return Quat(0, 0, 0, 1);
        } 
        if (XYZ::equals(direction, XYZ(0, -1, 0))) {
            return Quat(direction, 0);
        }
        decimal m1 = 1;
        decimal m2 = XYZ::magnitude(direction);
        Quat out = Quat(
            XYZ(
                direction.Z,
                0,
                -direction.X
            ),
            m2 * direction.Y
        );
        return Quat::normalize(out);

    }
    static XYZ applyRotation(const XYZ& p, const Quat& r) { //https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
        //I dont understand ANY of this, but by god apparently its fast
        XYZ t = 2 * XYZ::cross(r, p);
        return p + r.W * t + cross(r, t);
    }
    Quat operator/(decimal div) const {
        return Quat(
            X / div,
            Y / div,
            Z / div,
            W / div
        );
    }
};




decimal Ssin(decimal theta) {
    return theta - (theta * theta * theta) / 6 + (theta * theta * theta * theta * theta) / 120;
}

decimal Scos(decimal theta) {
    return 1 - (theta * theta) / 2 + (theta * theta * theta * theta) / 24;
}

#define LOOKUP_SIZE 1024
#define LOOKUP_MASK LOOKUP_SIZE-1

XYZ lookup[LOOKUP_SIZE];



namespace VecLib{
    static void prep() {
        for (int i = 0; i < LOOKUP_SIZE; i++) {
            decimal x;
            decimal y;
            decimal z;
            do {
                x = fRand(-1, 1);
                y = fRand(0, 1);
                z = fRand(-1, 1);
            } while (x * x + y * y + z * z > 1);
           lookup[i]= XYZ(x, y, z) / sqrt(x * x + y * y + z * z);
        }
    }
    static XYZ lookup_hemi() {
        int index = xe_next() & (LOOKUP_MASK);
        cout << index << endl;
        return lookup[index];
    }
    static XYZ random_hemi() {
        decimal x;
        decimal y;
        decimal z;
        //do {
            x = fRand(-1, 1);
            y = fRand(0, 1);
            z = fRand(-1, 1);
        //} while (x * x + y * y + z * z > 1);
        return XYZ(x, y, z)/sqrt(x * x + y * y + z * z);
    }
    static XYZ fast_random_cone(decimal spread) {
        decimal y_f = fRand(0, spread);
        decimal z_f = fRand(-PI, PI);
        decimal mult = 1;
        if (z_f > PI / 2) {
            z_f = PI - z_f;
            mult = -1;
        }
        if (z_f < -PI / 2) {
            z_f = -PI - z_f;
            mult = -1;
        }
        return XYZ(Ssin(y_f) * Scos(z_f) * mult, Scos(z_f), Ssin(y_f) * Ssin(z_f));
    }
    static XYZ y_random_cone(decimal spread) {
        decimal y = fRand(1, spread);
        decimal rot = fRand(0, PI);
        decimal f_y = sqrt(1 - y * y);
        return XYZ(f_y * cos(rot), y, f_y * sin(rot));
    }
    static XYZ random_cone(decimal spread) {
        decimal y_f = fRand(0, spread);
        decimal z_f = fRand(0, 2 * PI);
        return XYZ(sin(y_f)*cos(z_f), cos(y_f), sin(y_f)*sin(z_f));
    }
    static XYZ aligned_random(decimal spread, const Quat& r) {
        return Quat::applyRotation(random_cone(spread), r);
    }
}
struct Matrix3x3 {
public:
    decimal data[3][3];
    Matrix3x3() {
        data[0][0] = 0;
        data[0][1] = 0;
        data[0][2] = 0;
        data[1][0] = 0;
        data[1][1] = 0;
        data[1][2] = 0;
        data[2][0] = 0;
        data[2][1] = 0;
        data[2][2] = 0;
    }
    Matrix3x3(decimal a, decimal b, decimal c, decimal d, decimal e, decimal f, decimal g, decimal h, decimal i) {
        data[0][0] = a;
        data[0][1] = b;
        data[0][2] = c;
        data[1][0] = d;
        data[1][1] = e;
        data[1][2] = f;
        data[2][0] = g;
        data[2][1] = h;
        data[2][2] = i;
    }
    static XYZ multiply_vertically(const Matrix3x3& mat, const XYZ& other) {
        return XYZ(
            other.X * mat.data[0][0] + other.Y * mat.data[1][0] + other.Z * mat.data[2][0],
            other.X * mat.data[0][1] + other.Y * mat.data[1][1] + other.Z * mat.data[2][1],
            other.X * mat.data[0][2] + other.Y * mat.data[1][2] + other.Z * mat.data[2][2]
        );
    }
    static XYZ multiply_horizontally(const Matrix3x3& mat, const XYZ& other) {
        return XYZ(
            other.X * mat.data[0][0] + other.Y * mat.data[0][1] + other.Z * mat.data[0][2],
            other.X * mat.data[1][0] + other.Y * mat.data[1][1] + other.Z * mat.data[1][2],
            other.X * mat.data[2][0] + other.Y * mat.data[2][1] + other.Z * mat.data[2][2]
        );
    }
};

typedef tuple<XYZ, XYZ, XYZ> Face;
typedef tuple<bool, bool, bool> tri_bool;

class Material {
public:
    int material_type = 0;
    XYZ color = XYZ(0,0,0);
    decimal roughness = 0.1;
    decimal metallic  = 0;
    decimal specular  = 0;
    decimal emission  = 0;

    XYZ emissive_color = XYZ(0, 0, 0);

    decimal k = 0;
    decimal a_2 = 0;
    decimal spec_f;
    decimal diff_f;
    decimal diff_spread;
    XYZ diff_c;
    XYZ diff_t;
    XYZ spec_color;
    XYZ I_spec;

    Material() {}
    Material(XYZ _color) : color(_color) {}
    Material* clone() {
        Material* new_self = new Material(color);
        new_self->material_type = 0;
        new_self->roughness = roughness;
        new_self->metallic = metallic;
        new_self->specular = specular;
        new_self->emission = emission;
        new_self->emissive_color = emissive_color.clone();
        return new_self;
    }
    void prep() {
        k = pow(roughness + 1, 2) / 8;
        a_2 = roughness * roughness;
        spec_color = get_specular_color();
        I_spec = 1-get_specular_color();
        spec_f = get_specular_factor();
        diff_f = get_diffuse_factor();
        diff_c = get_diffuse_color();
        diff_t = diff_c * diff_f;
        diff_spread = roughness * roughness*PI/2;
    }
    decimal get_diffuse_factor() {
        return 1-min(max(metallic, specular), (decimal)1.0);
    }
    decimal get_specular_factor() {
        return min(max(metallic,specular),(decimal)1.0);
    }
    XYZ get_specular_color() {
        return XYZ::linear_mix(metallic, specular * XYZ(1, 1, 1), color);
    }
    XYZ get_diffuse_color() {
        return color / PI;
    }
    XYZ get_fresnel(XYZ light_slope,XYZ normal) {
        XYZ specular_color = get_specular_color();
        auto second_term = (1 - specular_color)*pow(1 - XYZ::dot(light_slope, normal), 5);
        return specular_color + second_term;
    }
    XYZ fast_fresnel(decimal dot_NI) {
        decimal g = 1 - dot_NI;
        return I_spec * (g * g * g * g * g);
    }
    decimal get_normal_distribution_beta(const XYZ& normal, const XYZ& half_vector){
        decimal a = roughness;
        decimal dot = XYZ::dot(normal, half_vector);
        decimal exponent = -(1 - pow(dot, 2)) / (pow(a*dot,2));
        decimal base = 1 / (PI * pow(a,2) * pow(dot, 4));
        decimal final = base * exp(exponent);

        return final;
    }
    decimal get_normal_distribution_GGXTR(const XYZ& normal, const XYZ& half_vector) {
        decimal a = roughness;
        decimal dot = XYZ::dot(normal, half_vector);
        decimal final = (a * a)
                          /
            (PI * pow((dot*dot)*(a*a-1)+1,2));

        return final;
    }
    decimal fast_normal_dist(const decimal dot_NH) {
        decimal a_2 = roughness * roughness;
        decimal g = dot_NH*dot_NH*(a_2-1)+1;
        return dot_NH*(a_2) / (PI * g*g);
    }
    decimal geoSchlickGGX(const XYZ& normal, const XYZ& vector,  decimal k) {

        decimal dot = XYZ::dot(normal,vector);

        return dot / (dot * (1 - k) + k);
        //return 1;

    }
    decimal fastGeo_both(const decimal dot_NO, const decimal dot_NI) {
        return dot_NO / (dot_NO * (1 - k) + k) * dot_NI / (dot_NI * (1 - k) + k);
    }
    XYZ fast_BRDF_co(const XYZ& normal, const XYZ& input_slope, XYZ& output_slope) {
        decimal dot_NI = XYZ::dot(normal, input_slope);
        decimal dot_NO = XYZ::dot(normal, output_slope);
        if (dot_NI <= 0 || dot_NO <= 0) {
            return XYZ(0, 0, 0);
        }
        XYZ half_vector = XYZ::normalize(XYZ::add(input_slope, output_slope));
        decimal dot_HO = XYZ::dot(half_vector, output_slope);
        decimal dot_NH = XYZ::dot(normal, half_vector);
        XYZ fresnel = fast_fresnel(dot_HO);
        decimal geo = fastGeo_both(dot_NO, dot_NI);
        decimal normal_dist = fast_normal_dist(dot_NH);
        decimal divisor = 4 * dot_NI * dot_NO;

        return (spec_f * (geo * normal_dist * fresnel) / divisor+diff_t)*dot_NI;

    }
    XYZ calculate_BRDF_coefficient(XYZ normal, XYZ input_slope, XYZ output_slope) {
        
        if (XYZ::dot(normal, input_slope) <= 0 || XYZ::dot(normal, output_slope) <= 0) {
            return XYZ(0,0,0);
        }

        XYZ half_vector = XYZ::add(input_slope, output_slope);
        half_vector.normalize();
        XYZ fresnel = get_fresnel(output_slope, half_vector);
        auto geo = geoSchlickGGX(normal, output_slope, k) * geoSchlickGGX(normal, input_slope, k);//Smiths method
        auto normal_dist = get_normal_distribution_GGXTR(normal, half_vector);
        decimal divisor = (4 * XYZ::dot(normal, input_slope) * XYZ::dot(normal, output_slope));
        //fresnel = XYZ(1,1,1);
        //geo = 1;
        //normal_dist = 1;
        //divisor = 4;
        XYZ specular_output = get_specular_factor()*(geo * normal_dist * fresnel) / divisor; //* ;
        XYZ diffuse_output = get_diffuse_factor() * get_diffuse_color();
        XYZ BRDF_color = specular_output + diffuse_output;
        //return fresnel;//specular_output;
        return BRDF_color* XYZ::dot(normal, input_slope);//
    }
    XYZ calculate_BRDF(XYZ normal, XYZ input_light, XYZ input_slope, XYZ output_slope) {
        return calculate_BRDF_coefficient(normal, input_slope, output_slope);
    }
    XYZ random_bounce(const Quat& diffuse_rotation, const Quat& reflection_rotation) {
        decimal prob = fRand(0, 1);
        if(prob<diff_f){
            return VecLib::aligned_random(PI/2, diffuse_rotation);
        }
        else {
            return VecLib::aligned_random(diff_spread, reflection_rotation);
        }

    }
    XYZ fast_bounce(const XYZ& normal, const XYZ& input_slope) {
        decimal prob = fRand(0, 1);
        if (prob < diff_f) {
            XYZ raw = XYZ::normalize(XYZ(fRand(-1, 1), fRand(-1, 1), fRand(-1, 1)));
            decimal result = XYZ::dot(raw, normal);
            if (result < 0) {
                return raw - 2 * result * normal;
            }
            else {
                return raw;
            }
        }
        else {
            XYZ raw = XYZ::reflect(input_slope, normal);
            decimal a = max(roughness - 0.1, 0.0);
            raw += XYZ(fRand(-a, a), fRand(-a, a), fRand(-a, a));
            raw.normalize();
            return raw;
        }
    }
    Material atUV(decimal x, decimal y) {

    }
    XYZ calculate_emissions() {
        return emission * emissive_color;
    }
};

class Primitive {
public:
    Material material;
    //pair<XYZ, XYZ> bounds;
    XYZ origin;
    int obj_type = 0;
    Primitive(Material _material) {
        material = _material;
    }
    virtual pair<XYZ, XYZ> get_bounds() { return pair<XYZ, XYZ>(); }

    virtual bool check_backface(XYZ& position) { return false; }

    virtual decimal distance(XYZ& position) { 
        cout << "[ERROR] Default distance function called! Something is misconfigured!" << endl;
        return 0;
    } //how to fix if this gets called: alcohol and crying
    virtual Material material_properties_at(XYZ) {
        return material;
    }
};

class PointLikeLight {
public:
    XYZ origin;
    decimal radius;
    Primitive* host;
    PointLikeLight(XYZ _origin, decimal _radius, Primitive* _host) : origin(_origin), radius(_radius), host(_host) {}
};

struct PackagedTri { //reduced memory footprint to improve cache perf
    XYZ p1;
    XYZ p2;
    XYZ p3;
    XYZ normal;
    XYZ p1p3;
    XYZ p1p2;
    Material* material;
    PackagedTri(const XYZ& p1o, const XYZ& p2o, const XYZ& p3o, const XYZ origin, Material& _material) {
        material = &_material;
        p1 = p1o + origin;
        p2 = p2o + origin;
        p3 = p3o + origin;
        p1p3 = p3 - p1;
        p1p2 = p2 - p1;
        normal = XYZ::normalize(XYZ::cross(p1p3, p1p2));
    }
};

class Tri :public Primitive {
public:
    XYZ p1o;
    XYZ p2o;
    XYZ p3o;
    XYZ origin;
    Tri(XYZ _origin, XYZ _p1o, XYZ _p2o, XYZ _p3o, Material _material) : Primitive(_material),
        origin(_origin), p1o(_p1o), p2o(_p2o), p3o(_p3o) {}
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
    static decimal intersection_check_MoTru(const PackagedTri& T, const XYZ& position, const XYZ& slope){
        XYZ pvec = XYZ::cross(slope,T.p1p3);
        decimal det = XYZ::dot(T.p1p2,pvec);
        // if the determinant is negative, the triangle is 'back facing'
        // if the determinant is close to 0, the ray misses the triangle
        //if (det < kEpsilon) return false;
        // ray and triangle are parallel if det is close to 0
        if (fabs(det) < kEpsilon) return -1;
        decimal invDet = 1 / det;

        XYZ tvec = position - T.p1;
        decimal u = XYZ::dot(tvec,pvec) * invDet;
        if (u < 0 || u > 1) return -1;

        XYZ qvec = XYZ::cross(tvec,T.p1p2);
        decimal v = XYZ::dot(slope,qvec) * invDet;
        if (v < 0 || u + v > 1) return -1;

        decimal t = XYZ::dot(T.p1p3,qvec) * invDet;
        return t;
    }
    static decimal intersection_check(const PackagedTri& T, const XYZ& position, const XYZ& slope) {
        return intersection_check_MoTru(T, position, slope);
    }
    static XYZ get_normal(XYZ normal, XYZ position) {
        return (XYZ::dot(normal, position) > 0) ? normal : -normal;
    }
};



class Sphere : public Primitive {
public:
    decimal radius = 0;
    Sphere(decimal _radius, XYZ _origin, Material _material): Primitive(_material){
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
    
    bool check_backface(XYZ &position) {
        return true;
    }

    decimal distance(XYZ& position) {
        return XYZ::distance(origin,position) - radius;
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
        decimal sqrt_comp = k * k - slope_2 * (j_2 - radius*radius);
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
        material = &sphere.material;
    }
};

class Plane : public Primitive {
public:
    XYZ normal;
    XYZ origin_offset;
    Plane( XYZ _normal, XYZ _origin, Material _material) : Primitive(_material) {
        origin = _origin;
        normal = XYZ::normalize(_normal);
        origin_offset = XYZ::dot(origin, normal) * normal;
        obj_type = 2;
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
        return XYZ::dot(position-origin_offset, normal);//this is like wildly simplified from position-(position after projection onto plane)
        //just ends in a dot product. Aint that weird
    }
    static decimal intersection_check(const XYZ& normal, const XYZ& origin_offset, const XYZ& position, const XYZ& slope) {
        decimal denom = XYZ::dot(normal, slope);
        if (denom == 0) {
            return -1;
        }
        return (XYZ::dot(origin_offset,normal) - XYZ::dot(position, normal)) / denom;
    }
};

struct PackagedPlane {
    XYZ normal;
    XYZ origin_offset;
    Material* material;
    PackagedPlane(Plane& plane) {
        normal = plane.normal;
        origin_offset = plane.origin_offset;
        material = &plane.material;
    }
};

class Mesh{
    vector<Face> faces;
    vector<XYZ> planes;
    Mesh() {}
    int import_file(string file_name) {
        return 0;
    }
    vector<Primitive> get_primitives() {

    }
};

class Lens {
public:
    int resolution_x;
    int resolution_y;
    int subdivision_size;
    void (*_prep)(Lens* self, int resolution_x, int resolution_y, int subdivision_size);
    XYZ (*_at)(Lens* self, int p_x, int p_y, int sample_index);
    void prep(int resolution_x, int resolution_y, int subdivision_size) {
        _prep(this, resolution_x, resolution_y, subdivision_size);
    }
    XYZ at(int p_x, int p_y, int sample_index) {
        return _at(this, p_x, p_y, sample_index);
    }

};

class RectLens : public Lens {
public:
    decimal width;
    decimal height;
    vector<XY> subdiv_offsets;
    vector<vector<XY>> outputs;
    RectLens(decimal _width, decimal _height) :
        width(_width), height(_height) {
        _at = _at_function;
        _prep = _prep_function;
    }
private:
    static void _prep_function(Lens* self, int resolution_x, int resolution_y, int subdivision_count) {
        RectLens* self_rect = (RectLens*)self;
        decimal width = self_rect->width;
        decimal height = self_rect->height;
        decimal partial_width =  width / resolution_x;
        decimal partial_height = height / resolution_y;
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
        return XYZ(self_rect->outputs[p_y][p_x] + self_rect->subdiv_offsets[sample_index]);   
    }
};

class Scene {
public:
    
    //Im trying something here. Im going to organize objects into groups, and then explicitly process each one when casting. Pain the ass for me, but in theory it makes for cleaner execution and code.
    //my justification here is that polymorphism can make for some really just nasty fucking code. Doing this will force me to make more modular code and more explicit optimizations.
    //C++ just wasnt made to be flexible :(

    int object_count = 0;
    vector<Sphere*> spheres;
    vector<Plane*> planes;
    vector<Tri*> tris;

    vector<PointLikeLight*> pointlike_lights;

    int current_resolution_x;
    int current_resolution_y;
    Scene() {}
    
    void register_sphere(Sphere* sphere) {
        object_count++;
        sphere->material.prep();
        spheres.push_back(sphere);
        if (sphere->material.emission > 0) {
            pointlike_lights.push_back(new PointLikeLight(sphere->origin, sphere->radius,(Primitive*)sphere));
        }
    }

    void register_plane(Plane* plane) {
        object_count++;
        plane->material.prep();
        planes.push_back(plane);
    }

    void register_tri(Tri* tri) {
        object_count++;
        tri->material.prep();
        tris.push_back(tri);
    }

    
private:
    
};

struct PackagedRay {
    XYZ position;
    XYZ slope;
    XYZ coefficient;
    XYZ* output;
    char remaining_bounces;
    char remaining_monte_carlo;
    int monte_diffuse;
    short monte_shadow;
    bool check_lighting;
    PackagedRay() {};
    PackagedRay(XYZ _position, XYZ _slope, XYZ co, XYZ* out, char bounces, char monte_carlo_bounces, int monte_dif, short monte_shad, bool _check_lighting) :
        position(_position),
        slope(_slope),
        coefficient(co),
        output(out),
        monte_diffuse(monte_dif),
        monte_shadow(monte_shad),
        remaining_bounces(bounces),
        remaining_monte_carlo(monte_carlo_bounces),
        check_lighting(_check_lighting) {}
    void move(decimal distance) {
        position += slope * distance;
    }
};

struct CameraConfig {
    int reflection_bounces = 8;
    int diffuse_bounces = 1;
};

class Camera {
    //Ill note this is a fairly non-direct class. the original was more streamlined, but I decided to update it to remove gunk
    //and to allow progressive output updates during monte carlo rendering
public:
    XYZ position;
    Scene* scene;

    Lens* lens;
    XYZ focal_position;

    decimal pre_scaled_gain = 1;
    decimal post_scaled_gain = 1;

    vector<vector<vector<XYZ*>>> ray_outputs;
    vector<vector<XYZ*>> image_output;

    int current_resolution_x;
    int current_resolution_y;

    Camera(Scene* _scene, XYZ _position, Lens* _lens, decimal focal_distance) : scene(_scene),
        lens(_lens), position(_position), focal_position(XYZ(0,0,-focal_distance)) {}

    Camera(Scene* _scene, XYZ _position, Lens* _lens, XYZ _focal_position) : scene(_scene),
        lens(_lens), position(_position), focal_position(_focal_position) {}

    void prep(int resolution_x, int resolution_y, int samples) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;
        lens->prep(resolution_x, resolution_y, samples);
    }

    XYZ slope_at(int p_x, int p_y, int sample_index) {
        return XYZ::slope(focal_position,lens->at(p_x, p_y, sample_index));
    }
    
    decimal luminance(XYZ v)
    {
        return XYZ::dot(v, XYZ(0.2126, 0.7152, 0.0722));
    }

    XYZ change_luminance(XYZ c_in, decimal l_out)
    {
        decimal l_in = luminance(c_in);
        return c_in * (l_out / l_in);
    }
    /*
    decimal realistic_response(decimal f, decimal iso) // https://graphics-programming.org/resources/tonemapping/index.html
    {
        f = clamp(f, 0.0f, iso); // Clamp to [0, iso]
        f /= iso; // Convert to [0, 1]

        // Returns 1.0 if the index is out-of-bounds
        auto get_or_one = [](const auto& arr, size_t index)
        {
            return index < arr.size() ? arr[index] : 1.0;
        };

        // std::upper_bound uses a binary search to find the position of f in camera_irradiance
        auto upper = 100;//std::upper_bound(camera_irradiance.begin(), camera_irradiance.end(), f);
        size_t idx = std::distance(camera_irradiance.begin(), upper);

        double low_irradiance = camera_irradiance[idx];
        double high_irradiance = get_or_one(camera_irradiance, idx + 1);
        double lerp_param = (f - low_irradiance) / (high_irradiance - low_irradiance);

        double low_val = 1;//camera_intensity[idx];
        double high_val = get_or_one(camera_intensity, idx + 1);

        // LERPing isn't really necessary for RGB8 (as the curve is sampled with 1024 points)
        return clamp(lerp((float)low_val, (float)high_val, (float)lerp_param), 0.0f, 1.0f);
    }
    */
    decimal Filmic_curve(decimal t) {
        return 0.371 * (sqrt(t) + 0.28257 * log(t) + 1.69542);
    }
    XYZ reinhard_extended_luminance(XYZ v, decimal max_white_l)
    {
        decimal l_old = luminance(v);
        decimal numerator = l_old * (1.0 + (l_old / (max_white_l * max_white_l)));
        decimal l_new = numerator / (1.0 + l_old);
        return change_luminance(v, l_new);
    }
    const Matrix3x3 ACESInputMat = Matrix3x3( //https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
            0.59719, 0.35458, 0.04823,
            0.07600, 0.90834, 0.01566,
            0.02840, 0.13383, 0.83777
    );
    const Matrix3x3 ACESOutputMat = Matrix3x3(
         1.60475, -0.53108, -0.07367,
        -0.10208,  1.10813, -0.00605,
        -0.00327, -0.07276,  1.07602
    );
    XYZ RRTAndODTFit(XYZ v)
    {
        XYZ a = v * (v + 0.0245786f) - 0.000090537f;
        XYZ b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
        return a / b;
    }
    XYZ ACESFitted(XYZ color)
    {
        color = Matrix3x3::multiply_horizontally(ACESInputMat,color);

        // Apply RRT and ODT
        color = RRTAndODTFit(color);

        color = Matrix3x3::multiply_horizontally(ACESOutputMat, color);

        return color;
    }
    XYZ FastACES(XYZ x) {
        decimal a = 2.51;
        decimal b = 0.03;
        decimal c = 2.43;
        decimal d = 0.59;
        decimal e = 0.14;
        return (x * (a * x + b)) / (x * (c * x + d) + e);
    }
    XYZ post_process(XYZ luminence) {
        XYZ scaled_return = luminence * pre_scaled_gain;
        decimal scalar_value = scaled_return.magnitude();
        //scaled_return = scaled_return * Filmic_curve(scalar_value);
        //scaled_return = scaled_return*log10(luminance(scaled_return) + 1) / luminance(scaled_return) * post_scaled_gain;//maybe more correct? dunno.
        scaled_return = XYZ::log(scaled_return + 1) * post_scaled_gain;
        //scaled_return = ACESFitted(scaled_return);
        //scaled_return = FastACES(scaled_return);
        //scaled_return = reinhard_extended_luminance(scaled_return, 200);
        scaled_return = XYZ::clamp(scaled_return, 0, 1) * 255;
        return scaled_return;
    }
};

#define BLOCK_BOUNDS_CHECKING false
template <typename T, int max_size> struct DataBlock {
    T* data = new T[max_size];
    int size = 0;
    int tag = 0;

    void enqueue(T& in) {
        data[size] = in;
        size++;
    }

    bool isFull() {
        return size >= max_size;
    }
    T get_val(int index) const {
        return data[index];
    }
    T& get_ref(int index) {
        return data[index];
    }
    T& back() {
        return data[size - 1];
    }
    ~DataBlock() {
        delete [] data;
    }
};

struct CastResults {
    XYZ normal;
    Material* material;
    CastResults(XYZ _normal, Material* _mat) : normal(_normal), material(_mat) {}
};

struct Casting_Diagnostics {
    long reflections_cast = 0;
    long shadows_cast = 0;
    long diffuses_cast = 0;
    long long rays_processed = 0;
    chrono::nanoseconds duration;
};



//#define RAY_BLOCK_SIZE 65536
#define RAY_BLOCK_SIZE 16*4
typedef DataBlock<PackagedRay,RAY_BLOCK_SIZE> Ray_Block;

struct Render_Pair {
    Ray_Block* ray_queue;
    iXY::set* coords;
    void enqueue_ray(PackagedRay& ray, iXY coord) {
        ray_queue->enqueue(ray);
        coords->insert(coord);
    }
    ~Render_Pair() {
        delete ray_queue;
        delete coords;
    }
};

class Block_Manager {
public:
    vector<Render_Pair*> blocks;
    Render_Pair* current_block = nullptr;
    Render_Pair* enqueuement_block = nullptr;
    int max_size = 0;
    void enqueue_ray(PackagedRay& ray, iXY coord) {
        ensure_space();
        enqueuement_block->enqueue_ray(ray, coord);
    }
    void enqueuement_done() {
        commit_enqueuement();
        max_size = size();
    }
    Ray_Block* next() {
        delete_current();
        blocks.pop_back();
        current_block = blocks.back();
        return current_block->ray_queue;
    }
    iXY::set* changed_pixels() {
        return current_block->coords;
    }
    int size() {
        return blocks.size();
    }
private:
    void delete_current() {
        if (current_block == nullptr) {
            return;
        }
        delete current_block;
    }
    void replace_enqueuement() {
        enqueuement_block = new Render_Pair();
        enqueuement_block->ray_queue = new Ray_Block();
        enqueuement_block->coords = new iXY::set();
    }
    void commit_enqueuement() {
        blocks.push_back(enqueuement_block);
        replace_enqueuement();
    }
    void ensure_space() {
        if (enqueuement_block == nullptr) {
            replace_enqueuement();
        }
        if (enqueuement_block->ray_queue->isFull()) {
            commit_enqueuement();
        }
    }
};

class RayEngine {
public:
    vector<PackagedSphere> sphere_data; //stores primitives in a more packed data format for better cache optimizations.
    vector<PackagedPlane> plane_data; //Ive made quite a few mistakes throughout this project but I think Im finally on track with this
    vector<PackagedTri> tri_data;
    vector<PointLikeLight> lights;

    RayEngine() {}
    void prep() {
    }
    void load_scene_objects(Scene* scene) {
        for (auto s : scene->spheres) {
            sphere_data.push_back(PackagedSphere(*s));
        }
        for (auto p : scene->planes) {
            plane_data.push_back(PackagedPlane(*p));
        }
        for (auto t : scene->tris) {
            tri_data.push_back(PackagedTri(t->p1o,t->p2o,t->p3o,t->origin,t->material));
        }
        for (auto PLL : scene->pointlike_lights) {
            lights.push_back(*PLL);
        }
    }
    CastResults execute_ray_cast(XYZ& position,const XYZ& slope) {
        #define default_smallest_distance 9999999;
        decimal smallest_distance = default_smallest_distance;
        CastResults returner = CastResults(XYZ(-1,-1,-1),nullptr);
        for (PackagedSphere& s : sphere_data) {
            decimal distance = Sphere::intersection_check(s.origin, s.radius, position, slope);
            if (distance >= 0) {
                if (distance < smallest_distance) {
                    smallest_distance = distance;
                    returner.normal = Sphere::normal(s.origin, position+distance*slope);
                    returner.material = s.material;
                }
            }
        }
        for (PackagedPlane& p : plane_data) {
            decimal distance = Plane::intersection_check(p.normal, p.origin_offset, position, slope);
            if (distance >= 0) {
                if (distance < smallest_distance) {
                    smallest_distance = distance;
                    returner.normal = p.normal;
                    returner.material = p.material;
                }
            }
        }
        for (PackagedTri& t : tri_data) {
            decimal distance = Tri::intersection_check(t, position, slope);
            if (distance >= 0) {
                if (distance < smallest_distance) {
                    smallest_distance = distance;
                    returner.normal = Tri::get_normal(t.normal,position);
                    returner.material = t.material;
                }
            }
        }
        position += smallest_distance * slope;
        return returner;
    }
    void process_ray(Casting_Diagnostics& stats, PackagedRay& ray_data) {
        stats.rays_processed++;
        XYZ start = ray_data.position;
        CastResults results = execute_ray_cast(ray_data.position, ray_data.slope);
        if (results.material == nullptr) {
            (*ray_data.output) += ray_data.coefficient * XYZ(0, 0, 0);
        }
        else {
            if (!ray_data.check_lighting) {
                (*ray_data.output) += ray_data.coefficient * results.material->calculate_emissions();
            }
            int bounce_count = ray_data.monte_diffuse;
            int shadow_count = ray_data.monte_shadow;
            XYZ flipped_output = XYZ::flip(ray_data.slope);
            XYZ reflection_slope = XYZ::reflect(XYZ::flip(ray_data.slope), results.normal);
            Quat normal_rot = Quat::makeRotation(XYZ(0, 1, 0), results.normal);
            Quat reflection_rot = Quat::makeRotation(XYZ(0,1,0),reflection_slope);
            //ray_data.PreLL.value = XYZ(0, 0, 100);
            
            if (FULL_BRIGHT) {
                (*ray_data.output)+= results.material->color;
            }
            if (ray_data.remaining_bounces > 0) {
                if (ray_data.remaining_monte_carlo > 0) {
                    for (int i = 0; i < bounce_count; i++) {
                        //XYZ monte_slope = results.material->fast_bounce(results.normal, flipped_output);
                        XYZ monte_slope = results.material->random_bounce(normal_rot, reflection_rot);
                        XYZ return_coefficient = ray_data.coefficient * results.material->fast_BRDF_co(results.normal, monte_slope, flipped_output);
                        if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                            return;
                        }
                        return_coefficient = return_coefficient / (bounce_count) * 4;
                        auto ray = PackagedRay(
                            ray_data.position + NEAR_THRESHOLD * results.normal * 1.1,
                            monte_slope,
                            return_coefficient,
                            ray_data.output,
                            ray_data.remaining_bounces - 1,
                            ray_data.remaining_monte_carlo - 1,
                            ray_data.monte_diffuse / 4,
                            max(1, ray_data.monte_shadow / 4),
                            true
                        );
                        //enqueue_ray(ray);
                        process_ray(stats, ray);
                        stats.diffuses_cast++;
                    }
                }
                else {
                    XYZ return_coefficient = ray_data.coefficient * results.material->fast_BRDF_co(results.normal, reflection_slope, flipped_output);
                    if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                        return;
                    }
                    auto ray = PackagedRay(
                        ray_data.position + NEAR_THRESHOLD * results.normal * 1.1,
                        reflection_slope,
                        return_coefficient,
                        ray_data.output,
                        ray_data.remaining_bounces - 1,
                        0,
                        0,
                        1,
                        true
                    );
                    //enqueue_ray(ray);
                    process_ray(stats, ray);

                    stats.reflections_cast++;
                }
            }
            if (ray_data.check_lighting) {
                for (const PointLikeLight& PLL : lights) {
                    XYZ light_slope = XYZ::slope(ray_data.position, PLL.origin);
                    Quat light_rot = Quat::makeRotation(XYZ(0, 1, 0), light_slope);
                    decimal distance = XYZ::distance(ray_data.position, PLL.origin);
                    decimal half_arc = atan(PLL.radius / distance);
                    XYZ light_falloff_coefficient = ray_data.coefficient / (4 * PI * distance * distance);
                    XYZ output = XYZ(0, 0, 0);
                    for (int i = 0; i < shadow_count; i++) {
                        XYZ cast_slope;
                        if (shadow_count > 1) {
                            
                            //cast_slope = light_slope+ XYZ(fRand(-0.1, 0.1), fRand(-0.1, 0.1), fRand(-0.1, 0.1));
                            //cast_slope.normalize();
                            
                            cast_slope = VecLib::aligned_random(half_arc, light_rot);
                        }
                        else {
                            cast_slope = light_slope;
                        }
                        XYZ return_coefficient = light_falloff_coefficient * results.material->fast_BRDF_co(results.normal, cast_slope, flipped_output);
                        if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                            continue;
                        }
                        if (shadow_count > 1) {
                            return_coefficient = return_coefficient / (shadow_count);
                        }
                        /*
                        auto ray = PackagedRay(
                            ray_data.position + NEAR_THRESHOLD * results.normal * 4,
                            cast_slope,
                            return_coefficient,
                            ray_data.output,
                            0,//ray_data.remaining_bounces-1,
                            0,
                            0,
                            200,
                            false
                        );
                        //enqueue_ray(ray);
                        process_ray(stats, ray);
                        */
                        XYZ start_position = ray_data.position + NEAR_THRESHOLD * results.normal * 4;
                        CastResults res = execute_ray_cast(start_position, cast_slope);
                        if (res.material != nullptr) {
                            *ray_data.output += return_coefficient * res.material->calculate_emissions();
                            //output += return_coefficient * res.material->calculate_emissions();
                        }
                      
                        
                        stats.shadows_cast++;
                    }
                    //*ray_data.output += output;
                }
            }

        }
    }
    Casting_Diagnostics process_block(Ray_Block* block) {

        auto start_time = chrono::high_resolution_clock::now();

        Casting_Diagnostics stats;
        
        for (int ray_index = 0; ray_index < block->size; ray_index++) {
            PackagedRay& ray_data = block->get_ref(ray_index);
            process_ray(stats,ray_data);

        }
        auto end_time = chrono::high_resolution_clock::now();
        stats.duration = end_time - start_time;
        return stats;
    }
};


class SceneManager {
public:
    Camera* camera;
    Scene* scene;
    sf::RenderWindow* window;
    sf::Image canvas;
    sf::Texture texture;
    sf::Sprite sprite;
    RayEngine RE;
    Block_Manager BM;
    vector<vector<XYZ*>> raw_output;

    int current_resolution_x = 0;
    int current_resolution_y = 0;
    int current_samples_per_pixel = 0;

    double scalar_exponent = 0;

    chrono::steady_clock::time_point render_start;
    SceneManager(Camera* _camera, Scene* _scene):
    camera(_camera), scene(_scene){}
    void render(int resolution_x, int resolution_y, int samples_per_pixel = 1) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;
        current_samples_per_pixel = samples_per_pixel;

        BM = Block_Manager();

        render_start = chrono::high_resolution_clock::now();
        prep();
        create_window();
        enqueue_rays();
        cout <<endl << "+RAYCASTING+" << endl;
        run_engine(true);
        hold_window();

    }
private:
    void prep() {
        cout << "Prepping.........." << flush;

        auto prep_start = chrono::high_resolution_clock::now();

        for (int y = 0; y < current_resolution_y; y++) {
            vector<XYZ*> row;
            for (int x = 0; x < current_resolution_x; x++) {
                row.push_back(new XYZ(0, 0, 0));
            }
            raw_output.push_back(row);
        }

        camera->prep(current_resolution_x, current_resolution_y, current_samples_per_pixel);
        current_samples_per_pixel *= current_samples_per_pixel;
        RE.load_scene_objects(scene);
        RE.prep();

        canvas.create(current_resolution_x, current_resolution_y, sf::Color::Black);
        sprite.setScale(make_scale(), -make_scale());
        sprite.setPosition(0, make_scale()*current_resolution_y);
        
        auto prep_end = chrono::high_resolution_clock::now();

        cout << "Done [" << chrono::duration_cast<chrono::milliseconds>(prep_end - prep_start).count() << "ms]" << endl;

    }
    void enqueue_rays() {
        cout << "Queueing Rays............" << flush;
        auto queue_start = chrono::high_resolution_clock::now();

        XYZ emit_coord = camera->focal_position + camera->position;

        XYZ start_position = emit_coord;
        XYZ starting_coefficient = XYZ(1, 1, 1)/current_samples_per_pixel;
        int max_bounces = 8;
        int monte_bounce_count = 0;
        int diffuse_emit_count = 4;
        int lighting_emit_count = 16;

        for (int y = 0; y < current_resolution_y; y++) {
            for (int x = 0; x < current_resolution_x; x++) {
                for (int i = 0; i < current_samples_per_pixel; i++) {
                    XYZ* output_link = (raw_output[y][x]);
                    XYZ slope = camera->slope_at(x, y, i);
                    auto ray = PackagedRay(
                        start_position,
                        slope,
                        starting_coefficient,
                        output_link,
                        max_bounces,
                        monte_bounce_count,
                        diffuse_emit_count,
                        lighting_emit_count,
                        true
                    );
                    BM.enqueue_ray(ray, iXY(x, y));
                }
            }
        }
        BM.enqueuement_done();
        auto queue_end = chrono::high_resolution_clock::now();
        cout << "Done [" << chrono::duration_cast<chrono::milliseconds>(queue_end-queue_start).count() << "ms]" << endl;
    }
    void run_engine(bool verbose) {
        auto start_time = chrono::high_resolution_clock::now();
        int blocks_processed = 0;
        long long rays_cast = 0;
        auto frame_time = start_time;
        while (BM.size()>1) { //this returns true if theres still rays to be processed
            process_block(rays_cast,verbose);
            blocks_processed++;
            auto time = chrono::high_resolution_clock::now();
            if (chrono::duration_cast<chrono::milliseconds>(time - frame_time).count() > 16) {
                commit_canvas();
                handle_events();
                frame_time = time;
            }
        }
        cout << endl;

        auto end_time = chrono::high_resolution_clock::now();
        double millis = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        long long approx_rays = blocks_processed * RAY_BLOCK_SIZE;
        cout << "Ray Casting Done [" << to_string(millis) << "ms][" << blocks_processed << " Blocks Processed]" << endl;
        cout << "Approximate speed: " << intToEng((rays_cast / (millis / 1000.0))) << " rays/s";
    }
    void process_block(long long& rays_cast, bool verbose) {
        Ray_Block* block = BM.next();
        fire_casting_event(block, rays_cast, verbose);
        iXY::set* changes = BM.changed_pixels();
        for (iXY change : *changes) {
            commit_pixel(change.X, change.Y);
        }
    }
    void fire_casting_event(Ray_Block* block, long long& rays_cast, bool verbose) {

        auto diagnostics = RE.process_block(block);
        rays_cast += diagnostics.rays_processed;
        if (verbose) {
            auto now = chrono::high_resolution_clock::now();
            int seconds = chrono::duration_cast<chrono::seconds>(now - render_start).count();
            int minutes = seconds / 60;
            int hours = minutes / 60;
            seconds %= 60;
            minutes %= 60;

            string debug_output = "[" + iToFixedLength(hours, 2, "0") + ":" + iToFixedLength(minutes, 2, "0") + ":" + iToFixedLength(seconds, 2, "0") + "]";
            debug_output += "[" + iToFixedLength(chrono::duration_cast<chrono::milliseconds>(diagnostics.duration).count(),6) + "ms::";
            debug_output += iToFixedLength(diagnostics.reflections_cast,9) + "R::";
            debug_output+= iToFixedLength(diagnostics.shadows_cast,9) + "S::";
            debug_output += iToFixedLength(diagnostics.diffuses_cast,9) + "D]";
            int period_count =60 - debug_output.size();
            for (int i = 0; i < period_count; i++) {
                debug_output+=".";
            }
            debug_output += "[" + iToFixedLength(BM.size(), 5) + "/" + iToFixedLength(BM.max_size,5) + "]";
            debug_output += "[" + iToFixedLength(BM.size(), 5) + "]\r";
            cout << debug_output << flush;
        }



    }
    void create_window() {
        cout << "Opening Window...." << flush;
        window = new sf::RenderWindow(sf::VideoMode(current_resolution_x * make_scale(), current_resolution_y * make_scale()), "Render Window!");
        cout << "Done" << endl;
    }
    void commit_pixel(int x, int y) {
        XYZ color = camera->post_process(*raw_output[y][x]);
        auto sfcolor = sf::Color(color.X, color.Y, color.Z);
        canvas.setPixel(x, y, sfcolor);
    }
    void draw_sprite() {
        window->clear();
        window->draw(sprite);
        window->display();
    }
    void commit_canvas() {
        texture.loadFromImage(canvas);
        sprite.setTexture(texture, false);
        draw_sprite();
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
            draw_sprite();
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }
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
            scalar_exponent += amount/5;
            cout << make_scale() << " " << scalar_exponent << endl;
            double scale = make_scale();
            sprite.setScale(scale,-scale);

        }
    }
};


/*
Scene* load_cornell_box() {
    Spatial_Cube* scene = new Spatial_Cube(100, XYZ(0, 0, 0));

    

    decimal size = 1;

    decimal scalar = size / 555;

    decimal box_width = 555*scalar;
    decimal wall_offset = box_width / 2;

    int res_scalar = 10;
    int resolution_x = 108 * res_scalar;
    int resolution_y = 108 * res_scalar;
    Lens lens = RectLens(resolution_x, resolution_y, 0.025*scalar, resolution_y / resolution_x);
    Camera* camera = new Camera(scene, XYZ(0, 0, -800*scalar), lens);

    Material wall_left_mat = Material(XYZ(0, 0.5, 0.1));
    wall_left_mat.roughness = 0.5;
    wall_left_mat.specular = 0;
    Plane* wall_left = new Plane(XYZ(1, 0, 0), XYZ(-wall_offset, 0, 0), wall_left_mat);

    Material wall_right_mat = Material(XYZ(0.7,0,0));
    wall_right_mat.roughness = 0.5;
    wall_right_mat.specular = 0;
    Plane* wall_right = new Plane(XYZ(-1, 0, 0), XYZ(wall_offset, 0, 0), wall_right_mat);

    Material blank_wall_mat = Material(XYZ(1, 1, 1));
    blank_wall_mat.roughness = 0.2;
    blank_wall_mat.specular = 0.9;
    
    Plane* top_wall = new Plane(XYZ(0, -1, 0), XYZ(0, wall_offset, 0), blank_wall_mat);
    Plane* bottom_wall = new Plane(XYZ(0, 1, 0), XYZ(0, -wall_offset, 0), blank_wall_mat);
    Plane* back_wall = new Plane(XYZ(0, 0, -1), XYZ(0, 0, wall_offset), blank_wall_mat);

    Material light_mat;
    light_mat.color = XYZ(1, 1, 1);
    light_mat.emission = 10;
    light_mat.emissive_color = XYZ(1, 0.662, 0.341);
    Sphere* glow_sphere = new Sphere(0.05, XYZ(0, 0, 0), light_mat);

    scene->register_object(wall_right);
    scene->register_object(wall_left);
    scene->register_object(back_wall);
    scene->register_object(top_wall);
    scene->register_object(bottom_wall);
    scene->register_object(glow_sphere);

    return pair<Camera*, Spatial_Cube*>(camera, scene);
}
*/
SceneManager* load_default_scene() {

    Scene* scene = new Scene();

    Lens* lens = new RectLens(1, 0.5625);
    Camera* camera = new Camera(scene, XYZ(0, 0.2, -10), lens , XYZ(0,0,-1));


    vector<Sphere*> spheres;
    int n_sphere = 0;
    decimal max_pos = 7;
    decimal min_radius = 0.3;
    decimal max_radius = 1;
    for (int i = 0; i < n_sphere; i++) {
        decimal radius = fRand(min_radius, max_radius);
        decimal x = fRand(-max_pos, max_pos);
        decimal y = fRand(-max_pos, max_pos);
        decimal z = fRand(-max_pos, max_pos);
        decimal r = fRand(0, 1);
        decimal g = fRand(0, 1);
        decimal b = fRand(0, 1);
        XYZ origin = XYZ(x, y, z);
        Material mat = Material(XYZ(r, g, b));
        mat.roughness = 1;
        Sphere* sphere =new Sphere(radius, origin, mat);
        spheres.push_back(sphere);

    }

    Material sphere_mat;
    sphere_mat.color = XYZ(0.2, 0.2, 1);//XYZ(0.24725, 0.1995, 0.0745);
    sphere_mat.specular = 0;
    sphere_mat.metallic = 0;
    sphere_mat.roughness = 1;
    Sphere* my_sphere =new Sphere(1, XYZ(0, 0, 0), sphere_mat);

    Material sphere_2_mat;
    sphere_2_mat.color = XYZ(1, 0.5, 0.5);//XYZ(0.24725, 0.1995, 0.0745);
    sphere_2_mat.specular = 0;
    sphere_2_mat.metallic = 0;
    sphere_2_mat.roughness = 1;
    Sphere* my_sphere_2 = new Sphere(0.4, XYZ(-0.9, -0.6, -0.7), sphere_2_mat);
    //Sphere* my_sphere_2 = new Sphere(0.4, XYZ(0, -0.6, 0.5), sphere_2_mat);

    Material tri_mat;
    tri_mat.color = XYZ(1, 0.5, 1);//XYZ(0.24725, 0.1995, 0.0745);
    tri_mat.specular = 0.8;
    tri_mat.metallic = 0;
    tri_mat.roughness = 0.05;
    Tri* my_tri = new Tri(XYZ(2, 0, 0), XYZ(3, -0.7, -3), XYZ(0, 3, 0), XYZ(-3, -0.7, 3), tri_mat);

    Material light_mat;
    light_mat.color = XYZ(1, 1, 1);
    light_mat.emission = 5000;
    light_mat.emissive_color = XYZ(1, 1, 1);
    //light_mat.emissive_color = XYZ(1, 0.662, 0.341);
    Sphere* glow_sphere =new Sphere(1, XYZ(-3, 2, 0), light_mat);


    Material floor_mat;
    floor_mat.color = XYZ(1, 1, 1);
    floor_mat.roughness = 1;
    floor_mat.metallic = 0;
    floor_mat.specular = 0;
    Plane* floor =new Plane(XYZ(0, 1, 0), XYZ(0, -1, 0), floor_mat);


    scene->register_sphere(my_sphere);
    scene->register_sphere(my_sphere_2);

    scene->register_tri(my_tri);

    scene->register_sphere(glow_sphere);
    scene->register_plane(floor);
    for (int i = 0; i < spheres.size(); i++) {
        scene->register_sphere(spheres.at(i));
    }

    SceneManager* SM = new SceneManager(camera, scene);

    return SM;
}



int main()
{
    srand(0);

    xe_seed[0] = rand();
    xe_seed[1] = rand();
    xe_seed[2] = rand();
    xe_seed[3] = rand();

    XYZ direction = XYZ(1, 1, 0);
    direction.normalize();

    Quat rot = Quat::makeRotation(XYZ(0,1,0),direction);

    cout << "A::-1::-1::-1::0::2::A::1::0::0::0::0;" << endl;
    cout << "A::1::1::1::0::2::A::1::0::0::0::0;" << endl;
    for (float i = -PI/2; i <= PI/2; i += PI / 100) {
        //cout << "(" << i << "," << Scos(i) << ")" << endl;
    }
    //exit(0);
    /*
    int count = 1000;
    XYZ out = XYZ();
    VecLib::prep();
    auto start_1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count;i++) {
        volatile XYZ k = VecLib::random_cone(0.1);
    }
    auto end_1 = chrono::high_resolution_clock::now();
    auto start_2 = chrono::high_resolution_clock::now();
    decimal f = cos(0.1);
    for (int i = 0; i < count;i++) {
        volatile XYZ k = (VecLib::y_random_cone(f));
    }
    auto end_2 = chrono::high_resolution_clock::now();
    auto start_3 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count;i++) {
        volatile XYZ k = VecLib::random_hemi();
    }
    auto end_3 = chrono::high_resolution_clock::now();
    auto start_4 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count;i++) {
        volatile XYZ k = VecLib::lookup_hemi();
    }
    auto end_4 = chrono::high_resolution_clock::now();
    cout << out << endl;
    cout << intToEng((double)count/chrono::duration_cast<chrono::milliseconds>(end_1 - start_1).count()*1000) << "/s" << endl;
    cout << intToEng((double)count/chrono::duration_cast<chrono::milliseconds>(end_2 - start_2).count()*1000) << "/s" << endl;
    cout << intToEng((double)count / chrono::duration_cast<chrono::milliseconds>(end_3 - start_3).count() * 1000) << "/s" << endl;
    cout << intToEng((double)count / chrono::duration_cast<chrono::milliseconds>(end_4 - start_4).count() * 1000) << "/s" << endl;
    exit(0);*/
    
    for (int i = 0; i < 0; i++) {
        //XYZ p = VecLib::random_cone(0.1);
        XYZ p = VecLib::aligned_random(0.1, rot);
        cout << "A::" << p.X << "::" << p.Y << "::" << p.Z << "::0::2::A::1::0::0::0::0;" << endl;
    }
    //exit(0);

    SceneManager* scene_manager = load_default_scene();
   
    scene_manager->render(1920/2,1080/2,5);

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
