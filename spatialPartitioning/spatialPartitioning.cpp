// spatialPartitioning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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
#include <thread>

#include <vector>
#include <array>
#include <unordered_set>
#include <queue>

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

#define PI 3.14159265358979323846264

#define kEpsilon 0.000000001

#define FULL_BRIGHT false
#define DEPTH_VIEW false

//defines program precision
#define LEAF_SIZE 2

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
    XY operator*(const decimal& scalar) const {
        return XY(X * scalar, Y * scalar);
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
    XYZ(decimal s) : X(s), Y(s), Z(s) {}
    XYZ(decimal _X, decimal _Y, decimal _Z) : X(_X), Y(_Y), Z(_Z) {}
    XYZ(XY xy, decimal _Z) : X(xy.X), Y(xy.Y), Z(_Z) {}
    XYZ(XY xy) : XYZ(xy, 0) {}
    XYZ clone() {
        return XYZ(X, Y, Z);
    }
    decimal operator[](const int n) const {
        if (n == 0) return X;
        if (n == 1) return Y;
        if (n == 2) return Z;
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
    static decimal maxComponent(const XYZ& point) {
        return std::max(point.X, std::max(point.Y, point.Z));
    }
    static decimal minComponent(const XYZ& point) {
        return std::min(point.X, std::min(point.Y, point.Z));
    }
    static XYZ min(const XYZ& point, decimal value) {
        return XYZ(
            (point.X < value) ? point.X : value,
            (point.Y < value) ? point.Y : value,
            (point.Z < value) ? point.Z : value
        );

    }
    static XYZ min(decimal value, XYZ point) {
        return XYZ::min(point, value);
    }
    static XYZ min(const XYZ& point, const XYZ& other) {
        return XYZ(
            (point.X < other.X) ? point.X : other.X,
            (point.Y < other.Y) ? point.Y : other.Y,
            (point.Z < other.Z) ? point.Z : other.Z
        );
    }
    static XYZ max(const XYZ& point, const XYZ& other) {
        return XYZ(
            (point.X > other.X) ? point.X : other.X,
            (point.Y > other.Y) ? point.Y : other.Y,
            (point.Z > other.Z) ? point.Z : other.Z
        );
    }
    static XYZ max(const XYZ& point,const decimal& value) {
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
    static XYZ _rtslope(const XYZ& point, const XYZ& other) {
        decimal distance = XYZ::distance(point, other);
        return (other - point) / distance;
    }
    static XYZ _dotslope(const XYZ& point, const XYZ& other) {
        XYZ delta = other - point;
        return delta / XYZ::dot(delta, delta);
    }
    static XYZ slope(const XYZ& point, const XYZ& other) {
        return _rtslope(point, other);
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
    static XYZ random(decimal nRx, decimal Rx, decimal nRy, decimal Ry, decimal nRz, decimal Rz) {
        return XYZ(
            fRand(nRx, Rx),
            fRand(nRy, Ry),
            fRand(nRz, Rz)
        );
    }
    static XYZ random(decimal Rx, decimal Ry, decimal Rz) {
        return random(
            -Rx, Rx, -Ry, Ry, -Rz, Rz
        );
    }
    static XYZ random(decimal range = 1) {
        return random(range, range, range);
    }
    string to_string() {
        return "("+std::to_string(X) + ", " + std::to_string(Y) + ", " + std::to_string(Z)+")";
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
    XYZ operator *= (const XYZ & other) {
        X *= other.X;
        Y *= other.Y;
        Z *= other.Z;
        return *this;
    }
    XYZ operator *= (const decimal& scalar) {
        X *= scalar;
        Y *= scalar;
        Z *= scalar;
        return *this;
    }
    XYZ operator-() {
        return XYZ(-X, -Y, -Z);
    }
    bool operator !=(const XYZ& other){
        return !((X == other.X) && (Y == other.Y) && (Z == other.Z));
    }
    bool operator !=(XYZ& other){
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
    struct less_than_x_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2)
        {
            return (point1.X < point2.X);
        }
    };
    struct less_than_y_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2)
        {
            return (point1.Y < point2.Y);
        }
    };
    struct less_than_z_operator {
        inline bool operator() (const XYZ& point1, const XYZ& point2)
        {
            return (point1.Z < point2.Z);
        }
    };
    struct less_than_x_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2)
        {
            return (point1->X < point2->X);
        }
    };
    struct less_than_y_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2)
        {
            return (point1->Y < point2->Y);
        }
    };
    struct less_than_z_operator_p {
        inline bool operator() (const XYZ* point1, const XYZ* point2)
        {
            return (point1->Z < point2->Z);
        }
    };
};
XYZ operator*(const decimal& self, const XYZ& point){
    return point * self;
}
XYZ operator-(const decimal& self, const XYZ& point) {
    return point + self;
}
XYZ operator/(const decimal& self, const XYZ& point) {
    return XYZ::divide(self,point);
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
    decimal magnitude() {
        return Quat::dot(*this, *this);//inefficient but whatever
    }
    static decimal dot(const Quat& q1, const Quat& q2) {
        return q1.X * q2.X + q1.Y * q2.Y + q1.Z * q2.Z + q1.W * q2.W;
    }
    static Quat multiply(const Quat& q1, const Quat& q2) {
        return Quat(
            q1.W * q2.X + q1.X * q2.W + q1.Y * q2.Z - q1.Z * q2.Y,
            q1.W * q2.Y - q1.X * q2.Z + q1.Y * q2.W + q1.Z * q2.X,
            q1.W * q2.Z + q1.X * q2.Y - q1.Y * q2.X + q1.Z * q2.W,
            q1.W * q2.W - q1.X * q2.X - q1.Y * q2.Y - q1.Z * q2.Z
        );
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

#define LOOKUP_SIZE_TRIG 1024
#define LOOKUP_TRIG_MASK LOOKUP_SIZE_TRIG-1
#define LOOKUP_SIZE_BIASED_HEMI 4096
#define LOOKUP_BIASED_HEMI_MASK LOOKUP_SIZE_BIASED_HEMI-1

decimal sin_lookup[LOOKUP_SIZE_TRIG];
decimal cos_lookup[LOOKUP_SIZE_TRIG];
XYZ biased_hemi_lookup[LOOKUP_SIZE_BIASED_HEMI];

namespace VecLib{
    static void prep() {
        for (int i = 0; i < LOOKUP_SIZE_TRIG; i++) {
            sin_lookup[i] = sin((double)i / LOOKUP_SIZE_TRIG * 2 * PI);
            cos_lookup[i] = cos((double)i / LOOKUP_SIZE_TRIG * 2 * PI);
        }
        for (int i = 0; i < LOOKUP_SIZE_BIASED_HEMI; i++) {
            double r1 = fRand(0,1);
            double r2 = fRand(0, 1);
            double f1 = acos(2 * r1 - 1) - PI / 2;
            double f2 = 2 * PI * r2;
            XYZ pos = XYZ(0,1,0)+XYZ(sin(f1) * cos(f2), sin(f1), cos(f1)* sin(f2));
            biased_hemi_lookup[i] = XYZ::slope(0, pos);
        }
    }
    static decimal Lsin(decimal input) {
        int index = ((int)(input * LOOKUP_SIZE_TRIG / (PI * 2))) & LOOKUP_TRIG_MASK;
        return sin_lookup[index];
    }
    static decimal Lcos(decimal input) {
        int index = ((int)(input * LOOKUP_SIZE_TRIG / (PI * 2))) & LOOKUP_TRIG_MASK;
        return cos_lookup[index];
    }
    static XYZ biased_random_hemi() {
        int index = xe_next()&LOOKUP_BIASED_HEMI_MASK;
        return biased_hemi_lookup[index];
    }
    static XYZ lookup_random_cone(decimal spread) {
        decimal y_f = fRand(0, spread);
        decimal z_f = fRand(0, 2 * PI);
        return XYZ(Lsin(y_f) * Lcos(z_f), Lcos(y_f), Lsin(y_f) * Lsin(z_f));
    }
    static XYZ y_random_cone(decimal spread) {
        decimal y = fRand(1, spread);
        decimal rot = fRand(0, PI);
        decimal f_y = sqrt(1 - y * y);
        return XYZ(f_y * cos(rot), y, f_y * sin(rot));
    }
    static XYZ aligned_random(decimal spread, const Quat& r) {
        return Quat::applyRotation(lookup_random_cone(spread), r);
    }
    static decimal surface_area(const XYZ& max, const XYZ& min) {
        decimal l_x = max.X - min.X;
        decimal l_y = max.Y - min.Y;
        decimal l_z = max.Z - min.Z;
        return (l_x * l_y + l_x * l_z + l_y * l_z)*2;
    }
    static decimal volume(const XYZ& max, const XYZ& min) {
        decimal l_x = max.X - min.X;
        decimal l_y = max.Y - min.Y;
        decimal l_z = max.Z - min.Z;
        return (l_x * l_y * l_z);
    }
    static bool between(const XYZ& p1, const XYZ& p2, const XYZ test) {
        return (
            (p1.X < test.X && test.X < p1.X) &&
            (p1.Y < test.Y && test.Y < p1.Y) &&
            (p1.Z < test.Z && test.Z < p1.Z)
            );
    }
    static bool betweenX(const XYZ& p1, const XYZ& p2, const XYZ test) {
        return (
            (p1.X < test.X && test.X < p1.X)
            );
    }
    static bool betweenY(const XYZ& p1, const XYZ& p2, const XYZ test) {
        return (
            (p1.Y < test.Y && test.Y < p1.Y)
            );
    }
    static bool betweenZ(const XYZ& p1, const XYZ& p2, const XYZ test) {
        return (
            (p1.Z < test.Z && test.Z < p1.Z)
            );
    }
    static bool volumeClip(const XYZ& max_1, const XYZ& min_1, const XYZ& max_2, const XYZ& min_2) {
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
    static bool volumeContains(const XYZ& max, const XYZ& min, const XYZ& point) {
        bool t1 = (point.X<max.X && point.X>min.X);
        bool t2 = (point.Y<max.Y && point.Y>min.Y);
        bool t3 = (point.Z<max.Z && point.Z>min.Z);
        return t1&&t2&&t3;
    }
}
struct Matrix3x3 {
public:
    decimal data[9];
    Matrix3x3() {
        data[0] = 0;
        data[1] = 0;
        data[2] = 0;
        data[3] = 0;
        data[4] = 0;
        data[5] = 0;
        data[6] = 0;
        data[7] = 0;
        data[8] = 0;
    }
    Matrix3x3(decimal a, decimal b, decimal c, decimal d, decimal e, decimal f, decimal g, decimal h, decimal i) {
        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
        data[4] = e;
        data[5] = f;
        data[6] = g;
        data[7] = h;
        data[8] = i;
    }
    static Matrix3x3 quatToMatrix(Quat q) {
        Matrix3x3 m;
        decimal s = q.magnitude();
        m.data[0] = 1 - 2 * s * (q.Y * q.Y + q.Z * q.Z);
        m.data[1] = 2 * s * (q.X * q.Y - q.Z * q.W);
        m.data[2] = 2 * s * (q.X * q.Z + q.Y * q.W);
        m.data[3] = 2 * s * (q.X * q.Y + q.Z * q.W);
        m.data[4] = 1 - 2 * s * (q.X * q.X + q.Z * q.Z);
        m.data[5] = 2 * s * (q.Y * q.Z - q.X * q.W);
        m.data[6] = 2 * s * (q.X * q.Z - q.Y * q.W);
        m.data[7] = 2 * s * (q.Y * q.Z + q.X * q.W);
        m.data[8] = 1 - 2 * s * (q.X * q.X + q.Y * q.Y);
        return m;
    }
    static XYZ applyRotationMatrix(const XYZ& point, const Matrix3x3& m) {
        return XYZ(
            point.X * m.data[0] + point.Y * m.data[1] + point.Z * m.data[2],
            point.X * m.data[3] + point.Y * m.data[4] + point.Z * m.data[5],
            point.X * m.data[6] + point.Y * m.data[7] + point.Z * m.data[8]
        );
    }
    static XYZ aligned_random(decimal spread, const Matrix3x3& r) {
        return applyRotationMatrix(VecLib::lookup_random_cone(spread),r);
    }
    static XYZ aligned_biased_hemi(const Matrix3x3& r) {
        return applyRotationMatrix(VecLib::biased_random_hemi(), r);
    }
    static XYZ multiply_vertically(const Matrix3x3& mat, const XYZ& other) {
        return XYZ(
            other.X * mat.data[0] + other.Y * mat.data[3] + other.Z * mat.data[6],
            other.X * mat.data[1] + other.Y * mat.data[4] + other.Z * mat.data[7],
            other.X * mat.data[2] + other.Y * mat.data[5] + other.Z * mat.data[8]
        );
    }
    static XYZ multiply_horizontally(const Matrix3x3& mat, const XYZ& other) {
        return XYZ(
            other.X * mat.data[0] + other.Y * mat.data[1] + other.Z * mat.data[2],
            other.X * mat.data[3] + other.Y * mat.data[4] + other.Z * mat.data[5],
            other.X * mat.data[6] + other.Y * mat.data[7] + other.Z * mat.data[8]
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
    decimal spec_f = 0;
    decimal diff_f = 0;
    decimal diff_spread = 0;
    XYZ diff_c = XYZ();
    XYZ diff_t = XYZ();
    XYZ spec_color = XYZ();
    XYZ I_spec = XYZ();

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
        I_spec = 1-get_fresnel_0();
        spec_f = get_specular_factor();
        diff_f = get_diffuse_factor();
        diff_c = get_diffuse_color();
        diff_t = diff_c * diff_f;
        decimal r_f = roughness - 0.2;
        diff_spread = (r_f*r_f)*PI/2;
    }
    decimal get_diffuse_factor() const {
        return 1-get_specular_factor();
    }
    XYZ get_fresnel_0() const {
        return XYZ::linear_mix(metallic, XYZ(0.04), color);
    }
    XYZ get_diffuse_reflectance() const {
        return color * (1 - metallic);
    }
    decimal get_specular_factor() const {
        return min(max(metallic,specular),(decimal)1.0);
    }
    XYZ get_specular_factor_v2(decimal dot_NI) const {
        return fast_fresnel(dot_NI);
    }
    XYZ get_diffuse_factor_v2(decimal dot_NI) const {
        return XYZ(1)-fast_fresnel(dot_NI);
    }
    XYZ get_specular_color() const {
        return XYZ::linear_mix(metallic, specular * XYZ(1, 1, 1), color);
    }
    XYZ get_diffuse_color() const {
        return get_diffuse_reflectance() / PI;
    }
    XYZ get_fresnel(XYZ light_slope,XYZ normal) const {
        XYZ specular_color = get_fresnel_0();
        auto second_term = (1 - specular_color)*pow(1 - XYZ::dot(light_slope, normal), 5);
        return specular_color + second_term;
    }
    XYZ fast_fresnel(decimal dot_NI) const {
        decimal g = 1 - dot_NI;
        return I_spec * (g * g * g * g * g);
    }
    decimal get_normal_distribution_beta(const XYZ& normal, const XYZ& half_vector) const {
        decimal a = roughness;
        decimal dot = XYZ::dot(normal, half_vector);
        decimal exponent = -(1 - pow(dot, 2)) / (pow(a*dot,2));
        decimal base = 1 / (PI * pow(a,2) * pow(dot, 4));
        decimal final = base * exp(exponent);

        return final;
    }
    decimal get_normal_distribution_GGXTR(const XYZ& normal, const XYZ& half_vector) const {
        decimal a = roughness;
        decimal dot = XYZ::dot(normal, half_vector);
        decimal final = (a * a)
                          /
            (PI * pow((dot*dot)*(a*a-1)+1,2));

        return final;
    }
    decimal fast_normal_dist(const decimal dot_NH) const {
        decimal a_2 = roughness * roughness;
        decimal g = dot_NH*dot_NH*(a_2-1)+1;
        return dot_NH*(a_2) / (PI * g*g);
    }
    decimal geoSchlickGGX(const XYZ& normal, const XYZ& vector,  decimal k) const {

        decimal dot = XYZ::dot(normal,vector);

        return dot / (dot * (1 - k) + k);
        //return 1;

    }
    decimal fastGeo_both(const decimal dot_NO, const decimal dot_NI) const {
        return dot_NO / (dot_NO * (1 - k) + k) * dot_NI / (dot_NI * (1 - k) + k);
    }
    XYZ diffuse_BRDF(decimal dot_NI) const {
        return get_diffuse_factor_v2(dot_NI)*get_diffuse_reflectance();
    }
    XYZ specular_BRDF(const XYZ& normal, const XYZ& input_slope, XYZ& output_slope) const {
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

        return geo * normal_dist * fresnel/divisor;
    }
    XYZ fast_BRDF_co(const XYZ& normal, const XYZ& input_slope, XYZ& output_slope) const {
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
    XYZ calculate_BRDF_coefficient(XYZ normal, XYZ input_slope, XYZ output_slope) const {
        
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
    XYZ random_bounce(const Matrix3x3& diffuse_rotation, const Matrix3x3& reflection_rotation) const {
        decimal prob = fRand(0, 1);
        if(prob<diff_f){
            return Matrix3x3::aligned_random(PI/2, diffuse_rotation);
        }
        else {
            return Matrix3x3::aligned_random(diff_spread, reflection_rotation);
        }
    }
    XYZ diffuse_bounce(const Matrix3x3& diffuse_rotation) const {
        return Matrix3x3::aligned_biased_hemi(diffuse_rotation);
    }
    XYZ reflective_bounce(const Matrix3x3& reflection_rotation) const {
        return Matrix3x3::aligned_random(diff_spread, reflection_rotation);
    }
    Material atUV(decimal x, decimal y) const {

    }
    XYZ calculate_emissions() const {
        return emission * emissive_color;
    }
};

class Primitive {
public:
    Material* material;
    //pair<XYZ, XYZ> bounds;
    XYZ origin;
    int obj_type = 0;
    Primitive(Material* _material):material(_material) {}
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

class PointLikeLight {
public:
    XYZ origin;
    decimal radius;
    Primitive* host;
    PointLikeLight(XYZ _origin, decimal _radius, Primitive* _host) : origin(_origin), radius(_radius), host(_host) {}
};

class SpotLight {
    XYZ origin;
    XYZ pointing;
    Primitive* host;
    SpotLight(XYZ _origin, XYZ _pointing, Primitive* _host) :origin(_origin), pointing(_pointing), host(_host) {}
    static SpotLight fromNormal(XYZ origin, XYZ normal, Primitive* host) {
        return SpotLight(origin, XYZ::flip(normal), host);
    }
};

struct PackagedTri { //reduced memory footprint to improve cache perf
    XYZ p1;
    XYZ normal;
    XYZ p1p3;
    XYZ p1p2;
    Material* material;
    PackagedTri() {}
    PackagedTri(const XYZ& p1o, const XYZ& p2o, const XYZ& p3o, const XYZ origin, Material* _material) {
        XYZ _p2 = p2o + origin;
        XYZ _p3 = p3o + origin;
        PackagedTri(p1, _p2, _p3, _material);
    }
    PackagedTri(const XYZ _p1, const XYZ _p2, const XYZ _p3, Material* _material){
        material = _material;
        p1 = _p1;
        p1p3 = _p3 - p1;
        p1p2 = _p2 - p1;
        normal = XYZ::normalize(XYZ::cross(p1p3, p1p2));
    }
    static bool equals(const PackagedTri& t1, const PackagedTri& t2) {
        bool test1 = XYZ::equals(t1.p1,t2.p1);
        bool test2 = XYZ::equals(t1.normal, t2.normal);
        bool test3 = XYZ::equals(t1.p1p3, t2.p1p3);
        bool test4 = XYZ::equals(t1.p1p2, t2.p1p2);
        bool test5 = t1.material == t2.material;
        return test1 && test2 && test3 && test4 && test5;
    }
};

class Tri :public Primitive {
public:
    XYZ p1;
    XYZ p2;
    XYZ p3;
    XYZ midpoint;
    XYZ AABB_max;
    XYZ AABB_min;
    Tri(XYZ _p1, XYZ _p2, XYZ _p3, Material* _material) : Primitive(_material),
        p1(_p1), p2(_p2), p3(_p3) {
        midpoint = (p1 + p2 + p3) / 3;
        AABB_max = XYZ::max(p1, XYZ::max(p2, p3));
        AABB_min = XYZ::min(p1, XYZ::min(p2, p3));
    }
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
    static decimal intersection_check_MoTru(const PackagedTri& T, const XYZ& position, const XYZ& slope){
        XYZ pvec = XYZ::cross(slope,T.p1p3);
        decimal det = XYZ::dot(T.p1p2,pvec);
        // if the determinant is negative, the triangle is 'back facing'
        // if the determinant is close to 0, the ray misses the triangle
        if (det < kEpsilon) return -1;
        // ray and triangle are parallel if det is close to 0
        //if (fabs(det) < kEpsilon) return -1;
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
    PackagedTri pack() {
        return PackagedTri(
            p1, p2, p3, material
        );
    }
};



class Sphere : public Primitive {
public:
    decimal radius = 0;
    Sphere(decimal _radius, XYZ _origin, Material* _material): Primitive(_material){
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
        material = sphere.material;
    }
};

class Plane : public Primitive {
public:
    XYZ normal;
    XYZ origin_offset;
    Plane( XYZ _normal, XYZ _origin, Material* _material) : Primitive(_material) {
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
        material = plane.material;
    }
};

class Transformation {
private:
    Quat* rot_transform = nullptr;
    XYZ*  XYZ_transform = nullptr;
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
    void stack(XYZ& XYZ_transf, Quat& rotation) {
        if (XYZ_transform != nullptr) {
            XYZ_transf += *XYZ_transform;
        }
        if (rot_transform != nullptr) {
            rotation = Quat::multiply(rotation, *rot_transform);
        }
    }
    void apply(XYZ& position) {
        if (rot_transform != nullptr) {
            position = Quat::applyRotation(position, *rot_transform);
        }
        if (XYZ_transform != nullptr) {
            position += *XYZ_transform;
        }
    }
    static Transformation collpase(vector<Transformation>& transforms) {
        XYZ position = XYZ(0,0,0);
        Quat rotation = Quat(0,0,0,1);
        for (Transformation& t : transforms) {
            t.stack(position, rotation);
        }
        return Transformation(position,rotation);
    }
};

class Mesh {
public:
    vector<Tri> tris;
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
    XYZ origin;
    XYZ scale;
    vector<Transformation> transformations;
    
    vector<Mesh*> meshes;
    Object(XYZ _origin, XYZ _scale):
        origin(_origin), scale(_scale){
        
    }
    void addMesh(Mesh* M) {
        meshes.push_back(M);
    }
    void prep() {
        for (Mesh* m : meshes) {
            m->prep();
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
    Transformation final_transform() {
        return Transformation::collpase(transformations);
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
        np1 += origin;
        np2 += origin;
        np3 += origin;
        return PackagedTri(
            np1, np2, np3, t.material
        );
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
        np1 += origin;
        np2 += origin;
        np3 += origin;
        return Tri(
            np1, np2, np3, t.material
        );
    }
}


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
    int primitive_count = 0;
    vector<Sphere*> spheres;
    vector<Plane*> planes;
    vector<Tri*> tris;
    vector<Object*> objects;

    vector<PointLikeLight*> pointlike_lights;

    int current_resolution_x;
    int current_resolution_y;
    Scene() {}
    
    void register_sphere(Sphere* sphere) {
        object_count++;
        primitive_count++;
        sphere->material->prep();
        spheres.push_back(sphere);
        if (sphere->material->emission > 0) {
            pointlike_lights.push_back(new PointLikeLight(sphere->origin, sphere->radius,(Primitive*)sphere));
        }
    }

    void register_plane(Plane* plane) {
        object_count++;
        primitive_count++;
        plane->material->prep();
        planes.push_back(plane);
    }

    void register_tri(Tri* tri) {
        object_count++;
        primitive_count++;
        tri->material->prep();
        tris.push_back(tri);
    }

    void register_object(Object* obj) {
        object_count++;
        for (Mesh* M: obj->meshes) {
            M->prep();
            primitive_count+=M->primitive_count;
        }
        objects.push_back(obj);
    }

    
private:
    
};

//http://bannalia.blogspot.com/2015/06/cache-friendly-binary-search.html
int tree_total = 0;
int leaf_total = 0;
decimal shortest = 99;

class BVH {
public:
    XYZ max = XYZ(0,0,0);
    XYZ min = XYZ(0,0,0);
    vector<Tri*>* elements = nullptr;
    vector<PackagedTri> packedEl;
    BVH* c1 = nullptr;
    BVH* c2 = nullptr;
    BVH* parent;
    BVH() {    }
    BVH(vector<Tri*>* T_vec) {
        elements = new vector<Tri*>(*T_vec);
        leaf_total += T_vec->size();
        for (Tri* T_ptr : *T_vec) {
            max = XYZ::max(max, T_ptr->AABB_max);
            min = XYZ::min(min, T_ptr->AABB_min);
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
            if (tmin >= 0) {
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
    static decimal intersection_alt(const XYZ& max, const XYZ& min, const XYZ& origin, const XYZ& inv_slope) {
        XYZ t0 = (min - origin) * inv_slope;
        XYZ t1 = (max - origin) * inv_slope;
        XYZ tmin = XYZ::min(t0, t1);
        XYZ tmax = XYZ::max(t0, t1);
        decimal max_t = XYZ::minComponent(tmax);
        decimal min_t = XYZ::maxComponent(tmin);
        
        if (max_t >= min_t) { //introduces a branch, but itll probably be fine. Lets me order search checks by nearest intersection
            if (min_t >= 0) {
                return min_t;
            }
            else {
                return max_t;
            }
        }
        else {
            return -1;
        }
    }
    decimal intersection(const XYZ& origin, const XYZ& inv_slope) { //https://tavianator.com/2011/ray_box.html
        return BVH::intersection(max, min, origin, inv_slope);
    }
    void prep() {
        //XYZ avg = (max + min) / 2;
        //XYZ max_ad = max - avg;
        //XYZ min_ad = min - avg;
        //max += max_ad*5;
        //min += min_ad*5;
        if (c1 != nullptr) {
            c1->prep();
            c2->prep();
        }
        else {
            for (Tri* t : *elements) {
                packedEl.push_back(t->pack());
            }
            delete elements;
        }
    }
    struct WorkPacket {
        BVH* target;
        vector<Tri*>* contents;
        WorkPacket(BVH* t, vector<Tri*>* c) : target(t), contents(c) {}
    };
    void construct(vector<Tri*>* initial_geo) {
        queue<WorkPacket> packets;
        assertm((initial_geo->size() > 0), "BVH was creation was attempted with zero geometry. Did everything load right?");
        auto initial_packet = WorkPacket(this,initial_geo);
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
            evaluate_split(geo, *split, 1);
            //Split* split = probe(geo);
            auto bins = bin(split, geo);

            if (((size_t)abs((int)(bins.p_geo->size() - bins.n_geo->size()))) == geo->size()) {
                target->elements = new vector<Tri*>(*geo);
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
        cout << "BVH leaf count: " << i << endl;
    }
    pair<vector<BVH*>*,vector<Tri*>*> flatten(int depth = 0) {
        vector<BVH*>* l = new vector<BVH*>();
        vector<Tri*>* t = new vector<Tri*>();

        flatten(l,t,depth);
        return pair<vector<BVH*>*, vector<Tri*>*>(l,t);
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
            if (elements != nullptr) {
                for (auto T : *elements) {
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
        int s = 0;
        count(s);
        return s;
    }
    void count(int& s) {
        s+=this->packedEl.size();
        if (c1 != nullptr) {
            c1->size(s);
            c2->size(s);
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
            return VecLib::surface_area(max,min);
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
        Split(XYZ place, int _facing):placement(place), facing(_facing) {};
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
    AABB getBounds( vector<Tri*>* geo) {
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
            } else {
                results.n_geo->push_back(T_ptr);
            }
        }
        return results;
    }
    Split* probe(vector<Tri*>* geo) {
       

        //Split* best_split = new Split(0,0);
        //vector<reduced_tri> tris;
        //for (const Tri* T_ptr : *geo) {
        //    tris.push_back(reduced_tri(*T_ptr));
        //}
        //std::sort(tris.begin(), tris.end(), reduced_tri::less_than_axis_operator(0));
        //Split* test_split = new Split(tris.begin()->midpoint,0);
        //get_stats(geo, *test_split);
        //for (int i = 0; i < tris.size(); i++) {
        //    
        //}
        //vector<Split*> tested;
        //while (to_test.size() > 0) {
        //    evaluate_split(geo, *(to_test.back()), 1);
        //    tested.push_back(to_test.back());
        //    to_test.pop_back();
        //}
        //while (tested.size() > 1) {
        //    delete tested.back();
        //    tested.pop_back();
        //}
        //return tested[0];
        //return final;
    }
    static void operate_stats(Stats& stat, const XYZ& point_max, const XYZ& point_min) {
        stat.count++;
        stat.max = XYZ::max(point_max, stat.max);
        stat.min = XYZ::min(point_min, stat.min);   
    }
    static void get_stats(vector<Tri*>* geo, Split& split) {
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
    static decimal evaluate_split(vector<Tri*>* geo, Split& split, decimal SA_parent) {
        const decimal traversal_cost = 1;
        const decimal intersect_cost = 1;
        decimal Sa = split.p.SA()/SA_parent;
        decimal Sb = split.n.SA()/SA_parent;
        decimal Ha = Sa * split.p.count * intersect_cost;
        decimal Hb = Sb * split.n.count * intersect_cost;
        decimal final_h = traversal_cost + Ha + Hb;
        split.score = final_h;
        return final_h;
        //return (split.p.count+split.n.count)/geo->size()*(split.p.count+split.n.count);
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
    static pair<PackagedBVH*,vector<PackagedTri>*> collapse(BVH* top) {
        vector<PackagedTri>* tri_vec = new vector<PackagedTri>();
        int tri_count = top->count();
        int tree_size = top->size();
        PackagedBVH* out = (PackagedBVH*) _aligned_malloc((tree_size+1)*sizeof(PackagedBVH), 32);
        tri_vec->reserve(tri_count+1);
        PackagedBVH::collapse(out, tri_vec, top);
        return pair<PackagedBVH*, vector<PackagedTri>*>(out,tri_vec);
    }
private:
    static void collapse(PackagedBVH* vec, vector<PackagedTri>* tri_vec, BVH* top) {
        vector<pair<BVH*,int>> current;
        vector<pair<BVH*,int>> next;
        
        int accumulated_index = 0;
        current.push_back(pair<BVH*, int>(top, -1));
        while (current.size() > 0) {
            for (int i = 0; i < current.size(); i++) {
                auto selection = current.at(i);
                BVH* target = selection.first;
                int parent_index = selection.second;
                auto PBVH = PackagedBVH(target);
                PBVH.leaf_size = target->packedEl.size();
                if (PBVH.leaf_size > 0) {
                    PBVH.index = tri_vec->size();
                    for (int i = 0; i < PBVH.leaf_size; i++) {
                        tri_vec->push_back(target->packedEl[i]);
                    }
                }
                if (target->c1 != nullptr) {
                    current.push_back(pair<BVH*, int>(target->c1, accumulated_index));
                    current.push_back(pair<BVH*, int>(target->c2, accumulated_index));
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

    

    vector<vector<vector<XYZ*>>> ray_outputs;
    vector<vector<XYZ*>> image_output;

    int current_resolution_x = 0;
    int current_resolution_y = 0;

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
};

const Matrix3x3 ACESOutputMat = Matrix3x3(
    1.60475, -0.53108, -0.07367,
    -0.10208, 1.10813, -0.00605,
    -0.00327, -0.07276, 1.07602
);

const Matrix3x3 ACESInputMat = Matrix3x3( //https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    0.59719, 0.35458, 0.04823,
    0.07600, 0.90834, 0.01566,
    0.02840, 0.13383, 0.83777
);



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
    const Material* material;
    decimal distance;
    CastResults() : normal(XYZ()), material(nullptr), distance(999999999999999) {};
    CastResults(XYZ _normal, Material* _mat) : normal(_normal), material(_mat), distance(9999999999999) {}
};

struct Casting_Diagnostics {
    long reflections_cast = 0;
    long shadows_cast = 0;
    long diffuses_cast = 0;
    long long rays_processed = 0;
    chrono::nanoseconds duration;
};



//#define RAY_BLOCK_SIZE 65536
#define RAY_BLOCK_SIZE 16
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

    BVH* test_bvh = new BVH();
    PackagedBVH* flat_bvh;
    vector<PackagedTri>* flat_bvh_tris;

    RayEngine() {}
    void prep() {
    }
    void load_scene_objects(Scene* scene) {
        vector<Tri*> tris;
        for (auto s : scene->spheres) {
            sphere_data.push_back(PackagedSphere(*s));
        }
        for (auto p : scene->planes) {
            plane_data.push_back(PackagedPlane(*p));
        }
        for (auto t : scene->tris) {
            tri_data.push_back(t->pack());
            tris.push_back(t);
        }
        for (auto O : scene->objects) {
            Transformation final_transform = O->final_transform();
            for (Mesh* M : O->meshes) {
                for (Tri& T : M->tris) {
                    tri_data.push_back(Packers::transformedPack(T, final_transform, O->origin, O->scale));
                    tris.push_back(new Tri(Packers::transformT(T,final_transform,O->origin,O->scale)));
                }
            }
        }
        for (auto PLL : scene->pointlike_lights) {
            lights.push_back(*PLL);
        }
        cout << "constructing BVH....." << flush;
        test_bvh->construct(new vector<Tri*>(tris));
        test_bvh->prep();
        auto collapse_out = PackagedBVH::collapse(test_bvh);
        flat_bvh = collapse_out.first;
        flat_bvh_tris = collapse_out.second;
        int lcf = 0;
        int lcr = 0;
        cout << lcf << " " << lcr << endl;
        cout << "done" << endl;
    }
    void navRawBVH(const BVH* bvh, CastResults& res,const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        if (bvh->c1 != nullptr) {
            decimal t1 = BVH::intersection(bvh->c1->max,bvh->c1->min,position, inv_slope);
            decimal t2 = BVH::intersection(bvh->c2->max,bvh->c2->min,position, inv_slope);
            if (t1 <= 0 && t2 <= 0) {
                return;
            }
            if (t1 < 0) {
                navRawBVH(bvh->c2, res, position, slope, inv_slope);
                return;
            }
            if (t2 < 0) {
                navRawBVH(bvh->c1, res, position, slope, inv_slope);
                return;
            }
            //navRawBVH(bvh->c1, res, position, slope, inv_slope);
            //navRawBVH(bvh->c2, res, position, slope, inv_slope);
            //return;
            if (t1 < t2) {
                navRawBVH(bvh->c1, res, position, slope, inv_slope);
                if (t2<res.distance) {
                    navRawBVH(bvh->c2, res, position, slope, inv_slope);
                }
            }
            else {
                navRawBVH(bvh->c2, res, position, slope, inv_slope);
                if (t1<res.distance) {
                    navRawBVH(bvh->c1, res, position, slope, inv_slope);
                }
            } 
        }
        else {
            for (const PackagedTri& t : bvh->packedEl) {
                decimal distance = Tri::intersection_check(t, position, slope);
                if (distance >= 0) {
                    if (distance < res.distance) {
                        res.distance = distance;
                        res.normal = Tri::get_normal(t.normal, position);
                        res.material = t.material;
                    }
                }
            }
        }
    }
    void iterativeBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        int stack[64];
        int stack_index = 0;
        stack[0] = 0;
        while (stack_index >= 0) {
            const PackagedBVH& current = flat_bvh[stack[stack_index]];
            stack_index--;
            //if (distance <= 0 || distance >= res.distance) {
            //    continue;
            //}
            if (current.leaf_size == 0) {
                int c1_index = current.index;
                int c2_index = c1_index + 1;
                const PackagedBVH& c1 = flat_bvh[c1_index];
                const PackagedBVH& c2 = flat_bvh[c2_index];
                decimal t1 = BVH::intersection(c1.sMax, c1.sMin, position, inv_slope);
                decimal t2 = BVH::intersection(c2.sMax, c2.sMin, position, inv_slope);

                if ((t1 < 0 && t2 < 0) || (t1>res.distance&&t2>res.distance)) {
                    continue;
                }
                if (t1 < 0 || t1>res.distance) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    continue;
                }
                if (t2 < 0 || t2>res.distance) {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    continue;
                }
                if (t1 < t2) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    stack_index++;
                    stack[stack_index] = c1_index;
                }
                else {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    stack_index++;
                    stack[stack_index] = c2_index;
                }
                
            }
            else {
                for (int i = 0; i < current.leaf_size; i++) {
                    const PackagedTri& t = (*flat_bvh_tris)[i+current.index];
                    decimal distance = Tri::intersection_check(t, position, slope);
                    if (distance >= 0) {
                        if (distance < res.distance) {
                            res.distance = distance;
                            res.normal = Tri::get_normal(t.normal, position);
                            res.material = t.material;
                        }
                    }
                }
            }
        }
    }
    /*void navFlatBVH(const PackagedBVH& current, CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        if (current.c1_index!=-1) {
            const PackagedBVH& c1 = flat_bvh->at(current.c1_index);
            const PackagedBVH& c2 = flat_bvh->at(current.c1_index+1);
            decimal t1 = BVH::intersection(c1.sMax,c1.sMin,position, inv_slope);
            decimal t2 = BVH::intersection(c2.sMax,c2.sMin,position, inv_slope);
            if (t1 <= 0 && t2 <= 0) {
                return;
            }
            if (t1 <= 0) {
                navFlatBVH(c2, res, position, slope, inv_slope);
                return;
            }
            if (t2 <= 0) {
                navFlatBVH(c1, res, position, slope, inv_slope);
                return;
            }
            navFlatBVH(c1, res, position, slope, inv_slope);
            navFlatBVH(c2, res, position, slope, inv_slope);
            return;
            if (t1<t2) {
                navFlatBVH(c1, res, position, slope, inv_slope);
                if (t2 < res.distance) {
                    navFlatBVH(c2, res, position, slope, inv_slope);
                }
            }
            else {
                navFlatBVH(c2, res, position, slope, inv_slope);
                if (t1 < res.distance) {
                    navFlatBVH(c1, res, position, slope, inv_slope);
                }
            }
        }
        else {
            for(int i = 0; i < current.leaf_size; i++){
                const PackagedTri& t = current.elements[i];
                decimal distance = Tri::intersection_check(t, position, slope);
                if (distance >= 0) {
                    if (distance < res.distance) {
                        res.distance = distance;
                        res.normal = Tri::get_normal(t.normal, position);
                        res.material = t.material;
                    }
                }
            }
        }
    }*/
    void navBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        const PackagedBVH& top = flat_bvh[0];
        //navFlatBVHI(0, res, position, slope, inv_slope);
        iterativeBVH(res, position, slope, inv_slope);
        //navRawBVH(test_bvh, res, position, slope, inv_slope);

    }
    CastResults execute_bvh_cast(XYZ& position, const XYZ& slope) {
        CastResults returner;
        for (const PackagedSphere& s : sphere_data) {
            decimal distance = Sphere::intersection_check(s.origin, s.radius, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = Sphere::normal(s.origin, position + distance * slope);
                    returner.material = s.material;
                }
            }
        }
        for (const PackagedPlane& p : plane_data) {
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
        position += returner.distance * slope;
        return returner;
    }
    CastResults execute_naive_cast(XYZ& position,const XYZ& slope) {
        #define default_smallest_distance 9999999;
        decimal smallest_distance = default_smallest_distance;
        CastResults returner = CastResults(XYZ(-1,-1,-1),nullptr);
        for (const PackagedSphere& s : sphere_data) {
            decimal distance = Sphere::intersection_check(s.origin, s.radius, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = Sphere::normal(s.origin, position+distance*slope);
                    returner.material = s.material;
                }
            }
        }
        for (const PackagedPlane& p : plane_data) {
            decimal distance = Plane::intersection_check(p.normal, p.origin_offset, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = p.normal;
                    returner.material = p.material;
                }
            }
        }
        for (const PackagedTri& t : tri_data) {
            decimal distance = Tri::intersection_check(t, position, slope);
            if (distance >= 0) {
                if (distance < returner.distance) {
                    returner.distance = distance;
                    returner.normal = Tri::get_normal(t.normal,position);
                    returner.material = t.material;
                }
            }
        }
        position += returner.distance * slope;
        return returner;
    }
    CastResults execute_ray_cast(XYZ& position, const XYZ& slope) {
        return execute_bvh_cast(position, slope);
    }
    void process_ray(Casting_Diagnostics& stats, PackagedRay& ray_data) {
        stats.rays_processed++;
        XYZ start = ray_data.position;
        CastResults results = execute_ray_cast(ray_data.position, ray_data.slope);
        if (DEPTH_VIEW) {
            (*ray_data.output) += XYZ(1, 1, 1) * (4 - ray_data.position.Z) * 5;
            return;
        }
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
            Matrix3x3 normal_rot_m = Matrix3x3::quatToMatrix(normal_rot);
            Matrix3x3 reflection_rot_m = Matrix3x3::quatToMatrix(reflection_rot);
            //ray_data.PreLL.value = XYZ(0, 0, 100);
            
            if (FULL_BRIGHT) {
                (*ray_data.output)+= results.material->color;
            }
            if (ray_data.remaining_bounces > 0) {
                if (ray_data.remaining_monte_carlo > 0) {
                    for (int i = 0; i < bounce_count/2; i++) {
                        XYZ diffuse_slope = results.material->diffuse_bounce(normal_rot_m);
                        XYZ return_coefficient = ray_data.coefficient * results.material->diffuse_BRDF(XYZ::dot(results.normal,diffuse_slope));
                        if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                            continue;
                        }
                        return_coefficient = return_coefficient / (bounce_count/2);
                        auto ray = PackagedRay(
                            ray_data.position + NEAR_THRESHOLD * results.normal * 1.1,
                            diffuse_slope,
                            return_coefficient,
                            ray_data.output,
                            ray_data.remaining_bounces - 1,
                            ray_data.remaining_monte_carlo - 1,
                            ray_data.monte_diffuse / 4,
                            max(1, ray_data.monte_shadow / 4),
                            true
                        );
                        process_ray(stats, ray);
                        stats.diffuses_cast++;
                    }
                    for (int i = 0; i < bounce_count / 2; i++) {
                        XYZ specular_slope = results.material->reflective_bounce(reflection_rot_m);
                        XYZ return_coefficient = ray_data.coefficient * results.material->specular_BRDF(results.normal, specular_slope, flipped_output);
                        if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                            continue;
                        }
                        return_coefficient = specular_slope / (bounce_count/2);
                        auto ray = PackagedRay(
                            ray_data.position + NEAR_THRESHOLD * results.normal * 1.1,
                            specular_slope,
                            return_coefficient,
                            ray_data.output,
                            ray_data.remaining_bounces - 1,
                            ray_data.remaining_monte_carlo - 1,
                            ray_data.monte_diffuse / 4,
                            max(1, ray_data.monte_shadow / 4),
                            true
                        );
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
                    Matrix3x3 light_rot_m = Matrix3x3::quatToMatrix(light_rot);
                    decimal distance = XYZ::distance(ray_data.position, PLL.origin);
                    decimal half_arc = atan(PLL.radius / distance);
                    XYZ light_falloff_coefficient = ray_data.coefficient / (4 * PI * distance * distance);
                    XYZ output = XYZ(0, 0, 0);
                    for (int i = 0; i < shadow_count; i++) {
                        XYZ cast_slope;
                        if (shadow_count > 1) {
                            
                            //cast_slope = light_slope+ XYZ(fRand(-0.1, 0.1), fRand(-0.1, 0.1), fRand(-0.1, 0.1));
                            //cast_slope.normalize();
                            
                            cast_slope = Matrix3x3::aligned_random(half_arc, light_rot_m);
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
                        XYZ start_position = ray_data.position + NEAR_THRESHOLD * results.normal * 1.1;
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
            XY scaled_position = focus_position*make_scale();
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

        }
    }
};

namespace ImageHandler {
    XYZ post_process_pixel(XYZ luminence,int pre_scaled_gain, int post_scaled_gain) {
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

    static decimal Filmic_curve(decimal t) {
        return 0.371 * (sqrt(t) + 0.28257 * log(t) + 1.69542);
    }
    static XYZ reinhard_extended_luminance(XYZ v, decimal max_white_l)
    {
        decimal l_old = luminance(v);
        decimal numerator = l_old * (1.0 + (l_old / (max_white_l * max_white_l)));
        decimal l_new = numerator / (1.0 + l_old);
        return change_luminance(v, l_new);
    }


    static XYZ RRTAndODTFit(XYZ v)
    {
        XYZ a = v * (v + 0.0245786f) - 0.000090537f;
        XYZ b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
        return a / b;
    }
    static XYZ ACESFitted(XYZ color)
    {
        color = Matrix3x3::multiply_horizontally(ACESInputMat, color);

        // Apply RRT and ODT
        color = RRTAndODTFit(color);

        color = Matrix3x3::multiply_horizontally(ACESOutputMat, color);

        return color;
    }
    static XYZ FastACES(XYZ x) {
        decimal a = 2.51;
        decimal b = 0.03;
        decimal c = 2.43;
        decimal d = 0.59;
        decimal e = 0.14;
        return (x * (a * x + b)) / (x * (c * x + d) + e);
    }
    static void postProcessRaw(vector<vector<XYZ*>>* data) {
        for (vector<XYZ*>& data_row : *data) {
            for (XYZ*& pixel_ptr : data_row) {
                *pixel_ptr = ImageHandler::post_process_pixel(*pixel_ptr,1,1);
            }
        }
    }
    static GUIHandler* openImageVector(vector<vector<XYZ*>>* data) {
        int resX = data->at(0).size();
        int resY = data->size();
        cout << resX << " " << resY << endl;
        GUIHandler* GUI = new GUIHandler(resX,resY,0);
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




class SceneManager {
public:
    Camera* camera;
    Scene* scene;
    RayEngine RE;
    Block_Manager BM;
    GUIHandler GUI;
    vector<vector<XYZ*>> raw_output; 

    int current_resolution_x = 0;
    int current_resolution_y = 0;
    int subdivision_count = 0;
    int current_samples_per_pixel = 0;
    const int block_size = 32;

    int y_increment;
    int x_increment;

    chrono::steady_clock::time_point render_start;

    SceneManager(Camera* _camera, Scene* _scene):
    camera(_camera), scene(_scene){}

    void render(int resolution_x, int resolution_y, int _subdivision_count = 1) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;

        y_increment = ceil(current_resolution_y / block_size);
        x_increment = ceil(current_resolution_x / block_size);

        subdivision_count = _subdivision_count;
        current_samples_per_pixel = subdivision_count*subdivision_count;

        BM = Block_Manager();
        GUI = GUIHandler(current_resolution_x, current_resolution_y,block_size);

        render_start = chrono::high_resolution_clock::now();
        prep();
        enqueue_rays();
        GUI.create_window();
        prepGUI();
        cout <<endl << "+RAYCASTING+" << endl;
        run_engine(true);

    }
    void hold_window() {
        GUI.hold_window();
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

        camera->prep(current_resolution_x, current_resolution_y, subdivision_count);
        RE.load_scene_objects(scene);
        RE.prep();

        
        
        auto prep_end = chrono::high_resolution_clock::now();

        cout << "Done [" << chrono::duration_cast<chrono::milliseconds>(prep_end - prep_start).count() << "ms]" << endl;

    }
    void prepGUI() {
        GUI.add_focus(); 
    }
    void enqueue_rays() {
        cout << "Queueing Rays............" << flush;
        auto queue_start = chrono::high_resolution_clock::now();

        XYZ emit_coord = camera->focal_position + camera->position;

        XYZ start_position = emit_coord;
        XYZ starting_coefficient = XYZ(1, 1, 1)/current_samples_per_pixel;
        int max_bounces = 2;
        int monte_bounce_count = 2;
        int diffuse_emit_count = 256/current_samples_per_pixel;
        int lighting_emit_count = 256/current_samples_per_pixel;

        
        for (int i_y = 0; i_y < y_increment; i_y++) {
            for (int i_x = 0; i_x < x_increment; i_x++) {
                int offset_y = i_y * block_size;
                int offset_x = i_x * block_size;
                int final_y = offset_y + block_size;
                int final_x = offset_x + block_size;
                for (int y = offset_y; y < final_y;y++) {
                    for (int x = offset_x; x < final_x; x++) {
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
                GUI.commit_canvas();
                GUI.handle_events();
                frame_time = time;
            }
        }
        refresh_canvas();
        cout << endl;

        auto end_time = chrono::high_resolution_clock::now();
        double millis = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        long long approx_rays = blocks_processed * RAY_BLOCK_SIZE;
        cout << "Ray Casting Done [" << to_string(millis) << "ms][" << blocks_processed << " Blocks Processed]" << endl;
        cout << "Approximate speed: " << intToEng((rays_cast / (millis / 1000.0))) << " rays/s";
    }
    void refresh_canvas() {
        for (int y = 0; y < current_resolution_y; y++) {
            for (int x = 0; x < current_resolution_x; x++) {
                XYZ color = ImageHandler::post_process_pixel(*raw_output[y][x], 1, 1);
                GUI.commit_pixel(color, x, y);
            }
        }
        GUI.commit_canvas();
    }
    void process_block(long long& rays_cast, bool verbose) {
        Ray_Block* block = BM.next();
        fire_casting_event(block, rays_cast, false);
        iXY::set* changes = BM.changed_pixels();
        
        for (iXY change : *changes) {
            int block_x = floor(change.X / block_size)*block_size;
            int block_y = floor(change.Y / block_size)*block_size;
            GUI.set_focus_position(0, block_x, block_y);
            XYZ color = ImageHandler::post_process_pixel(*raw_output[change.Y][change.X],1,1);
            GUI.commit_pixel(color,change.X, change.Y);
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
        ifstream file(fName, ios::in);
        return loadObj(file, mat);
    }
    static Mesh* loadObj(ifstream& file, Material* mat) {
        Mesh* mesh = new Mesh();
        vector<XYZ> verts;
        for (std::string line; getline(file, line);){
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
            if (line_type == "f") {
                int i0 = stoi(words[1])-1;
                int i1 = stoi(words[2])-1;
                int i2 = stoi(words[3])-1;
                Tri f = Tri(verts[i0], verts[i1], verts[i2], mat);
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
            if (words.size()>1) {
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

    Lens* lens = new RectLens(1, 1);
    Camera* camera = new Camera(scene, XYZ(0, -0.3, -2), lens , XYZ(0,0,-2));


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
        Material* mat = new Material(XYZ(r, g, b));
        mat->roughness = 1;
        Sphere* sphere =new Sphere(radius, origin, mat);
        spheres.push_back(sphere);

    }

    Material* sphere_mat = new Material();
    sphere_mat->color = XYZ(0.2, 0.2, 1);//XYZ(0.24725, 0.1995, 0.0745);
    sphere_mat->specular = 0;
    sphere_mat->metallic = 1;
    sphere_mat->roughness = 0.2;
    Sphere* my_sphere =new Sphere(1, XYZ(0, 0, 0), sphere_mat);

    Material* sphere_2_mat = new Material();
    sphere_2_mat->color = XYZ(1, 0.5, 0.5);//XYZ(0.24725, 0.1995, 0.0745);
    sphere_2_mat->specular = 0;
    sphere_2_mat->metallic = 0;
    sphere_2_mat->roughness = 1;
    Sphere* my_sphere_2 = new Sphere(0.4, XYZ(-0.9, -0.6, -0.7), sphere_2_mat);
    //Sphere* my_sphere_2 = new Sphere(0.4, XYZ(0, -0.6, 0.5), sphere_2_mat);

    Material* tri_mat = new Material();
    tri_mat->color = XYZ(1, 0.5, 1);//XYZ(0.24725, 0.1995, 0.0745);
    tri_mat->specular = 0.8;
    tri_mat->metallic = 0;
    tri_mat->roughness = 0.05;
    Tri* my_tri = new Tri(XYZ(3, -0.7, -3), XYZ(0, 3, 0), XYZ(-3, -0.7, 3), tri_mat);

    Material* light_mat = new Material();
    light_mat->color = XYZ(1, 1, 1);
    light_mat->emission = 5000;
    light_mat->emissive_color = XYZ(1, 1, 1);
    //light_mat.emissive_color = XYZ(1, 0.662, 0.341);
    Sphere* glow_sphere =new Sphere(1, XYZ(-3, 2, 0), light_mat);


    Material* floor_mat = new Material();
    floor_mat->color = XYZ(1, 1, 1);
    floor_mat->roughness = 0.25;
    floor_mat->metallic = 0;
    floor_mat->specular = 0.5;
    Plane* floor =new Plane(XYZ(0, 1, 0), XYZ(0, -1, 0), floor_mat);

    Material* mesh_mat = new Material();
    mesh_mat->color = XYZ(1, 0.2, 0.2);
    //mesh_mat->color = XYZ(1);
    mesh_mat->roughness = 0.2;
    mesh_mat->specular = 0.5;
    mesh_mat->metallic = 0;
    //Mesh* mesh = FileManager::loadTriFile("C:\\Users\\Charlie\\Documents\\RobotLab.tri", mesh_mat);
    Mesh* mesh = FileManager::loadObjFile("C:\\Users\\Charlie\\Documents\\models\\objs\\mike.obj", mesh_mat);
    //Mesh* mesh = FileManager::loadMeshFile("C:\\Users\\Charlie\\source\\repos\\spatialPartitioning\\spatialPartitioning\\lowfumo.obj",mesh_mat);
    Object* m_obj = new Object(XYZ(0,0,0),1*XYZ(1,1,1));
    m_obj->addMesh(mesh);
   // m_obj->applyTransformXYZ(1, 0, 2);
    m_obj->applyTransformXYZ(0, 0, 2);
    m_obj->rotateX(0);
    m_obj->rotateY(PI / 4);
    m_obj->rotateZ(0);

    //scene->register_sphere(my_sphere);
    //scene->register_sphere(my_sphere_2);
    //scene->register_tri(my_tri);
    scene->register_sphere(glow_sphere);
    scene->register_plane(floor);
    scene->register_object(m_obj);
    for (int i = 0; i < spheres.size(); i++) {
        scene->register_sphere(spheres.at(i));
    }

    SceneManager* SM = new SceneManager(camera, scene);

    return SM;
}



int main()
{
    cout << std::numeric_limits<float>::epsilon() << endl;
    srand(0);

    xe_seed[0] = rand();
    xe_seed[1] = rand();
    xe_seed[2] = rand();
    xe_seed[3] = rand();

    
    XYZ direction = XYZ(1, 1, 0);
    direction.normalize();

    Quat rot = Quat::makeRotation(XYZ(0,1,0),direction);
    Matrix3x3 rot_m = Matrix3x3::quatToMatrix(rot);

    cout << "A::-1::-1::-1::0::2::A::1::0::0::0::0;" << endl;
    cout << "A::1::1::1::0::2::A::1::0::0::0::0;" << endl;
    for (float i = -PI/2; i <= PI/2; i += PI / 100) {
        //cout << "(" << i << "," << Scos(i) << ")" << endl;
    }
    
    //exit(0);
    
    int count = pow(10,8);
    XYZ out = XYZ();
    VecLib::prep();
    int num = 512;
    float t_c_e = 0;
    float t_s_e = 0;
    float t_c_e_a = 0;
    float t_s_e_a = 0;
    for (int i = 0; i < num; i++) {
        float f = (float)i / num;
        float c_error = VecLib::Lcos(f)-cos(f);
        float s_error = VecLib::Lsin(f)-sin(f);
        t_c_e += c_error;
        t_s_e += s_error;
        t_c_e_a += abs(c_error);
        t_s_e_a += abs(s_error);
    }
    cout << "cosine: total error: " << t_c_e << ", total absolute error: " << t_c_e_a << ", average error: " << t_c_e / num << ", average absolute error: " << t_c_e_a / num << endl;
    cout << "sine: total error: " << t_s_e << ", total absolute error: " << t_s_e_a << ", average error: " << t_s_e / num << ", average absolute error: " << t_s_e_a / num << endl;
    cout << endl;
    auto start_1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count;i++) {
        //out+= Matrix3x3::aligned_random(0.1,rot_m);
    }
    auto end_1 = chrono::high_resolution_clock::now();
    auto start_2 = chrono::high_resolution_clock::now();
    decimal f = cos(0.1);
    for (int i = 0; i < count;i++) {
        //out+= (VecLib::aligned_random(0.1,rot));
    }
    auto end_2 = chrono::high_resolution_clock::now();
    auto start_3 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count*2;i+=2) {
        //out+= XYZ::test_slope(XYZ::random(100), XYZ::random(100));
    }
    auto end_3 = chrono::high_resolution_clock::now();
    auto start_4 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count * 2;i += 2) {
        //out+= XYZ::slope(XYZ::random(100), XYZ::random(100));
    }
    auto end_4 = chrono::high_resolution_clock::now();
    auto start_5 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count;i++) {
        //volatile XYZ k = VecLib::lookup_random_cone(0.1);
    }
    auto end_5 = chrono::high_resolution_clock::now();
    cout << out << endl;
    cout << intToEng((double)count/chrono::duration_cast<chrono::milliseconds>(end_1 - start_1).count()*1000) << "/s" << endl;
    cout << intToEng((double)count/chrono::duration_cast<chrono::milliseconds>(end_2 - start_2).count()*1000) << "/s" << endl;
    cout << intToEng((double)count / chrono::duration_cast<chrono::milliseconds>(end_3 - start_3).count() * 1000) << "/s" << endl;
    cout << intToEng((double)count / chrono::duration_cast<chrono::milliseconds>(end_4 - start_4).count() * 1000) << "/s" << endl;
    cout << intToEng((double)count / chrono::duration_cast<chrono::milliseconds>(end_5 - start_5).count() * 1000) << "/s" << endl;

    //exit(0);
    /**
    for (int i = 0; i < 0; i++) {
        //XYZ p = VecLib::random_cone(0.1);
        XYZ p = VecLib::aligned_random(0.1, rot);
        cout << "A::" << p.X << "::" << p.Y << "::" << p.Z << "::0::2::A::1::0::0::0::0;" << endl;
    }
    */
    //exit(0);
    

    //GUIHandler* GUI = FileManager::openRawFile("outputs.raw");
    
    //GUI->hold_window();
    PackagedBVH pbvh = PackagedBVH(new BVH());
    cout << "asdfasdfasdf   " << sizeof(pbvh) << endl;
    SceneManager* scene_manager = load_default_scene();
    scene_manager->render(640*1.6,640*1.6,8);
    /*auto l = scene_manager->RE.test_bvh->flatten(-1);
    cout << l.first->size() << endl;
    string out_s = "join(";
    int max = min((int)l.first->size(), 1000);
    for (int i = 0; i < max;i++) {
        BVH* b = l.first->at(i);
        //Tri* t = l.second->at(i);
        //XYZ M = t->AABB_max;
        //XYZ m = t->AABB_min;
        XYZ M = b->max;
        XYZ m = b->min;
        string piece = "B(";
        piece += M.to_string();
        piece += ", ";
        piece += m.to_string();
        piece += ")";
        if (i < max - 1) {
            piece += ",";
        }
        piece += " ";
        out_s += piece;
    }
    //cout << out_s << ")" << flush;
    out_s = "triangle(";
    string v1 = "[";
    string v2 = "[";
    string v3 = "[";
    max = min((int)l.second->size(), 2000);
    for (int i = 0; i < max;i++) {
        if (i > 0) {
            v1 += ",";
            v2 += ",";
            v3 += ",";
        }
        Tri* t = l.second->at(i);
        v1 += t->p1.to_string();
        v2 += t->p2.to_string();
        v3 += t->p3.to_string();
    }
    v1 += "]";
    v2 += "]";
    v3 += "]";
    out_s += v1 + ", " + v2 + ", " + v3 + ")";
    //cout << out_s << flush;
    //exit(0);
    //scene_manager->render(1920 / 10, 1080 / 10, 1);
    */
    cout << endl << "writing out raw......." << flush;
    FileManager::writeRawFile(&scene_manager->raw_output, "outputs.raw");
    cout << "done" << endl;
    scene_manager->hold_window();

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
