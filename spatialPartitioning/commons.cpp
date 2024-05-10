#pragma once

#include "Config.h"
#include "XorRandom.h"
#include <chrono>
#include "commons.h"
#include <string>


using std::string;
using std::max;
using std::min;
using namespace std::chrono;


const __m256 AVX_ZEROS = _mm256_set1_ps(0);
const __m256 AVX_ONES = _mm256_set1_ps(1);
const __m256 AVX_HALFS = _mm256_set1_ps(0.5);
const __m256 AVX_NEG_ONES = _mm256_set1_ps(-1);
const __m256 AVX_PI = _mm256_set1_ps(PI);
const __m256 AVX_HALF_PI = _mm256_set1_ps(PI*0.5);
const __m256 AVX_3O2_PI = _mm256_set1_ps(PI*1.5);
const __m256 AVX_2PI = _mm256_set1_ps(2*PI);
const __m256 AVX_SIGN_MASK = _mm256_set1_ps(std::bit_cast<float>(0x80000000));
const __m256 AVX_ABS_MASK = _mm256_set1_ps(std::bit_cast<float>(~0x80000000));
float float_one = *((float*)&one);
const __m256 AVX_FINT_ONES = _mm256_set1_ps(float_one);
const __m256 AVX_PACKED_ONES = _mm256_set1_ps(std::bit_cast<float>(0xffffffff));

thread_local XorGen gen;



float clamp(float value, float low, float high) {
    return std::max(std::min(value, high), low);
}

float lerp(float a, float b, float factor) {
    return a + factor * (b - a);
}

std::vector<std::string> eng_levels = std::vector<std::string>{
    "k",
    "M",
    "B",
    "T",
    "qu",
    "Qu"
};

std::string padString(std::string s, std::string padder, int size) {
    std::string out = s;
    while (out.size() < size) {
        out += padder;
    }
    return out;
}

std::string iToFixedLength(int input, int length, std::string fill) {
    std::string in = std::to_string(input);
    int spaces = length - in.size();
    for (int i = 0; i < spaces;i++) {
        in = fill + in;
    }
    return in;
}

std::string intToEng(int inum, bool force_small) {
    double num = inum;
    std::string prefix = "";
    if (num < 1000) {
        return iToFixedLength(inum, 4) + " ";
    }
    for (auto level : eng_levels) {
        if (num < 1000) {
            std::string number = std::to_string(num);
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

std::vector<std::string> split_string(std::string input, char splitter) {
    std::vector<std::string> output;
    std::string s;
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