#pragma once

#include "Config.h"
#include "XorRandom.h"
#include <chrono>

#define NOMINMAX 1

#undef max
#undef min

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

extern const __m256 AVX_ZEROS;
const int one = 1;
extern float float_one;
extern const __m256 AVX_FINT_ONES;

extern thread_local XorGen gen;



float clamp(float value, float low, float high);

float lerp(float a, float b, float factor);

std::string padString(std::string s, std::string padder, int size);

std::string iToFixedLength(int input, int length, std::string fill = " ");

std::string intToEng(int inum, bool force_small = false);

std::vector<std::string> split_string(std::string input, char splitter);