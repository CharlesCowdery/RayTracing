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

    VecLib::prep(gen);

    int mode = 0;
    string fpath = "C:\\Users\\Charlie\\Documents\\models\\Scenes\\dining.glb";
    int scalar = 2;
    int resolution_x = 1080/scalar;
    int resolution_y = 1080/scalar;
    int subdivisions = 1;
    //GUIHandler* GUI = FileManager::openRawFile("outputs.raw");

    //GUI->hold_window();
    //SceneManager* scene_manager = load_cornell_box();

    SceneManager* scene_manager;


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
