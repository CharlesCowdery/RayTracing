#pragma once
#include "Materials.h"
#include "XYZ.h"
#include <chrono>

struct CastResults {
    static const float default_distance;
    uint32_t material;
    float distance;
    XYZ TUV = XYZ(-1);
    XY texUV = XY(-1);
    int tri_index = 0;
    CastResults();
    CastResults( short _mat);
};

struct Casting_Diagnostics {
    long reflections_cast = 0;
    long shadows_cast = 0;
    long diffuses_cast = 0;
    long long rays_processed = 0;
    std::chrono::nanoseconds duration;
};
