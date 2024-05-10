// spatialPartitioning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#pragma once
#define _AMD64_ true
#define _WIN64 true
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>

#include <tuple>
#include <math.h>
#include <format>
#include <string>
#include <ctime>
#include <chrono>
#include <functional>
#include <fstream>


#include <thread>
#include <vector>
#include <array>
#include <unordered_set>
#include <set>
#include <queue>
#include <processthreadsapi.h>

#include <immintrin.h>


//#include "imgui-master/imgui.h"
//#include "imgui-master/backends/imgui_impl_win32.h"
//#include "imgui-master/backends/imgui_impl_dx12.h"

#include "commons.h"
#include "SceneManager.h"
#include "FileLoading.h"
//#include "imgui_gui_maker.h"

/*Xorshiro256+ pseudorandom end*/








int main(int argc, char* argv[])
{
    srand(0);

    gen.seed(rand(), rand(), rand(), rand());

    //TCHAR buffer[MAX_PATH] = { 0 };
    //GetModuleFileName(NULL, buffer, MAX_PATH);
    //std::wstring::size_type pos = std::wstring(buffer).find_last_of(L"\\/");
    //auto str = std::wstring(buffer).substr(0, pos);
    //cout << string(str.begin(), str.end()) << endl;

    cout << sizeof(PackagedTri) << endl;

    VecLib::prep(gen);

    int mode = 0;
    string fpath = "C:\\Users\\Charlie\\Documents\\models\\Scenes\\kitchen.glb";
    int scalar = 2;
    int resolution_x = 1080/scalar;
    int resolution_y = 1080/scalar;
    int subdivisions = 5;
    if (subdivisions % 2 == 0)throw exception("I am lazy and even subdivisions play hell with the noise algo");
    //GUIHandler* GUI = FileManager::openRawFile("outputs.raw");

    //GUI->hold_window();
    //SceneManager* scene_manager = load_cornell_box();

    SceneManager* scene_manager;

    int test_num = 99999999*0+1;
    __m256* rands = (__m256*) _aligned_malloc(sizeof(__m256)*test_num,64);
    __m256 out = _mm256_set1_ps(0);

    int iteration_count = 10;
    long long total = 0;
    out = _mm256_set1_ps(0);

    unsigned int NAN_mask = 0x7f800000u;

    for (int i = 0; i < iteration_count; i++) {

        for (int i = 0; i < test_num; i++) {
            gen.regen_pool();
            rands[i] = _mm256_mul_ps(gen.pool, AVX_2PI);
            //rands[i] = _mm256_sub_ps(rands[i], AVX_HALF_PI);
        }


        steady_clock::time_point start = steady_clock::now();
#pragma loop(no_vector)
        for (int i = 0; i < test_num; i++) {
            __m256& reg = rands[i]; //AVX register filled with random values between 0-1
            float* reg_ptr = (float*)&reg;
            reg_ptr[0] = VecLib::full_sin(reg_ptr[0]); //run sine on each value
            reg_ptr[1] = VecLib::full_sin(reg_ptr[1]);
            reg_ptr[2] = VecLib::full_sin(reg_ptr[2]);
            reg_ptr[3] = VecLib::full_sin(reg_ptr[3]);
            reg_ptr[4] = VecLib::full_sin(reg_ptr[4]);
            reg_ptr[5] = VecLib::full_sin(reg_ptr[5]);
            reg_ptr[6] = VecLib::full_sin(reg_ptr[6]);
            reg_ptr[7] = VecLib::full_sin(reg_ptr[7]);
            //reg = VecLib::full_sin_AVX(reg);
            //add to an out register to prevent the test from being optimized out



            //out =  _mm256_add_ps(out, reg);
        }

        steady_clock::time_point end = steady_clock::now();
        long long millis = duration_cast<milliseconds>(end - start).count();
        long long nanos = duration_cast<nanoseconds>(end - start).count();
        std::cout << "Count: " << test_num << ". Test duration: " << millis << "ms. " << (int)(test_num / (millis / 1000.0)) << "ops/s. " << std::setprecision(3) << (float)nanos / test_num << " ns/op" << endl;
        total += millis;
        for (int i = 0; i < test_num; i++) {
            out = _mm256_add_ps(out, rands[i]);
        }
    }

    for(int i = 0; i < 8;i ++) cout << ((float*)&out)[i];

    cout << "Average test time: " << total / iteration_count << "ms"<<endl;

    for (int i = 0; i < test_num; i++) {
        gen.regen_pool();
        rands[i] = _mm256_mul_ps(gen.pool, AVX_PI);
        rands[i] = _mm256_sub_ps(rands[i], AVX_HALF_PI);
        rands[i] = _mm256_mul_ps(rands[i], _mm256_set1_ps(0.5));
    }

    out = _mm256_set1_ps(0);
    float max_error = 0;
    for (int i = 0; i < test_num; i++) {
        __m256& reg = rands[i];
        __m256 reg2;
        float* reg_ptr = (float*)&reg;//m256 with each value a random 0-1
        float* reg2_ptr = (float*)&reg2;
        reg2_ptr[0] = sin(reg_ptr[0]); //run sine on each value
        reg2_ptr[1] = sin(reg_ptr[1]);
        reg2_ptr[2] = sin(reg_ptr[2]);
        reg2_ptr[3] = sin(reg_ptr[3]);
        reg2_ptr[4] = sin(reg_ptr[4]);
        reg2_ptr[5] = sin(reg_ptr[5]);
        reg2_ptr[6] = sin(reg_ptr[6]);
        reg2_ptr[7] = sin(reg_ptr[7]);
        reg_ptr[0] = VecLib::full_sin(reg_ptr[0]); //run sine on each value
        reg_ptr[1] = VecLib::full_sin(reg_ptr[1]);
        reg_ptr[2] = VecLib::full_sin(reg_ptr[2]);
        reg_ptr[3] = VecLib::full_sin(reg_ptr[3]);
        reg_ptr[4] = VecLib::full_sin(reg_ptr[4]);
        reg_ptr[5] = VecLib::full_sin(reg_ptr[5]);
        reg_ptr[6] = VecLib::full_sin(reg_ptr[6]);
        reg_ptr[7] = VecLib::full_sin(reg_ptr[7]);
        //reg = VecLib::full_sin_AVX(reg);
        __m256 diff = _mm256_and_ps(AVX_ABS_MASK,_mm256_sub_ps(reg, reg2));
        float* diff_ptr = (float*)&diff;
        max_error = std::max(max_error, diff_ptr[0]);
        max_error = std::max(max_error, diff_ptr[1]);
        max_error = std::max(max_error, diff_ptr[2]);
        max_error = std::max(max_error, diff_ptr[3]);
        max_error = std::max(max_error, diff_ptr[4]);
        max_error = std::max(max_error, diff_ptr[5]);
        max_error = std::max(max_error, diff_ptr[6]);
        max_error = std::max(max_error, diff_ptr[7]);

        //add to an out register to prevent the test from being optimized out
        out = _mm256_add_ps(out,diff);
    }

    double total_error = 0;
    for (int i = 0; i < 8; i++) {
        total_error+= ((float*)&out)[i];
    }
    cout << "Total error: " << to_string(total_error) << endl;
    cout << "Average error: " << std::setprecision(15) << total_error / test_num /8.0f << endl;
    cout << "Maximum error: " << std::setprecision(15) << max_error << endl;

    _aligned_free(rands);


    switch (mode) {
    case 0:
        scene_manager = FileManager::loadGLTFFile(fpath);
        scene_manager->render(resolution_x, resolution_y, subdivisions, 1);
        scene_manager->hold_window();
        break;
    }


    //scene_manager->render(1080/4, 1080/4, 2);

    //cout << endl << "writing out raw......." << flush;
    //FileManager::writeRawFile(&scene_manager->raw_output, "outputs.raw");
    //cout << "done" << endl;

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
