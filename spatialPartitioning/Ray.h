#pragma once
#include "commons.h"
#include "XYZ.h"
#include "Materials.h"


struct PackagedRay {
    XYZ position;
    XYZ slope;
    char generation = 0;
    XYZ* output;
    XYZ coefficient;
    short IOR_stack[IOR_STACK_SIZE];
    char IOR_size = 0;
    PackagedRay() {};
    PackagedRay(XYZ _position, XYZ _slope, char gen, XYZ coef = XYZ(1,1,1)) :
        position(_position),
        slope(_slope),
        generation(gen),
        coefficient(coef)
    {}
    void move(float distance) {
        position += slope * distance;
    }
};

struct MinimumRay {
    XYZ position;
    XYZ inv_slope;
    uint64_t index;
};

struct AVX_AABB_ray {
    m256_vec3 fusedposition;
    m256_vec3 inv_slope;
    unsigned int index;
};