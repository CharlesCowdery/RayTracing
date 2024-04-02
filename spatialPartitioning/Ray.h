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
    Material* IOR_stack[IOR_STACK_SIZE];
    PackagedRay() {};
    PackagedRay(XYZ _position, XYZ _slope, char gen) :
        position(_position),
        slope(_slope),
        generation(gen)
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