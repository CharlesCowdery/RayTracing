#pragma once
#include <chrono>
#include <vector>
#include <map>
#include <iostream>

#include "commons.h"
#include "BVH.h"
#include "Camera.h"

using namespace std::chrono;
using std::cout;
using std::endl;
using std::flush;
using std::pair;
using std::map;

class PackagedScene {
public:
    vector<Material> materials;
    vector<PackagedTri> tri_data;
    vector<PTri_AVX> avx_tri_data;
    vector<Light*> lights;
    vector<PackagedTri> emissive_tris;
    PackagedBVH* flat_bvh = nullptr;
    vector<BVH_AVX> avx_bvh;
    short monte_carlo_generations = 3;
    short max_generations = 10;
    short monte_carlo_max = 32;
    float monte_carlo_modifier = 1.0 / 4;
};

class Scene {
public:

    int object_count = 0;
    int primitive_count = 0;

    vector<Object*> objects;
    vector<Camera*> cameras;
    Camera* camera;

    vector<Material> materials;

    vector<Light*> lights;

    int current_resolution_x;
    int current_resolution_y;
    Scene() {}

    void register_object(Object* obj) {
        object_count++;
        for (Mesh* M : obj->meshes) {
            primitive_count += M->primitive_count;
        }
        objects.push_back(obj);
    }
    void register_camera(Camera* cam) {
        cameras.push_back(cam);
    }
    void register_light(Light* li) {
        lights.push_back(li);
    }
    void merge(Scene& mergee) {
        for (Object* O : mergee.objects) {
            objects.push_back(O);
        }
        for (Camera* cam : mergee.cameras) {
            cameras.push_back(cam);
        }
        for (Light* light : mergee.lights) {
            lights.push_back(light);
        }
    }

    void prep(int res_x, int res_y, int subdiv_count) {
        auto start = high_resolution_clock::now();
        cout << padString("[Scene] Prepping", ".", 100) << flush;
        for (Object* O : objects) {
            O->prep(materials);
        }
        camera = cameras[0];
        for (Camera* camera : cameras) {
            camera->prep(res_x, res_y, subdiv_count);
        }
        auto end = high_resolution_clock::now();
        cout << "[Done][" << duration_cast<milliseconds>(end - start).count() << "ms]" << endl;
    }

    PackagedScene* package() {
        vector<Tri> tris;
        BVH* bvh = new BVH();


        PackagedScene* PS = new PackagedScene();

        PS->materials = materials;

        auto data_group_start = high_resolution_clock::now();
        cout << padString("[Scene] Regrouping data", ".", 100) << flush;

        for (auto O : objects) {
            O->fetchData(tris);
        }
        for (Light* L : lights) {
            PS->lights.push_back(L);
        }

        auto data_group_end = high_resolution_clock::now();
        cout << "[Done][" << duration_cast<milliseconds>(data_group_end - data_group_start).count() << "ms]" << endl;

        //for (int i = 0; i < tris.size();i++) {
        //    Tri T = tris[i];
        //    cout << "triangle(" << T.p1.to_string() << "," << T.p2.to_string() << "," << T.p3.to_string() << ")" << endl;
        //}

        if (tris.size() > 0) {
            auto BVH_con_start = high_resolution_clock::now();
            cout << padString("[BVH] Constructing", ".", 100) << flush;

            bvh->construct_dangerous(tris); //note, this will cause the BVH to break when tris moves out of scope

            auto BVH_con_end = high_resolution_clock::now();
            cout << "[Done][" << duration_cast<milliseconds>(BVH_con_end - BVH_con_start).count() << "ms]" << endl;
            auto BVH_collapse_start = high_resolution_clock::now();
            cout << padString("[BVH] Flattening", ".", 100) << flush;

            PS->emissive_tris = bvh->get_emissive_tris(materials);
#if USE_AVX_BVH
            auto collapse_out = BVH_AVX::collapse(bvh);
            PS->avx_bvh = *collapse_out.first;
            PS->tri_data = *collapse_out.second;

#else
            auto collapse_out = PackagedBVH::collapse(bvh);
            PS->flat_bvh = collapse_out.first;
            PS->tri_data = *collapse_out.second;
#endif
            //bvh->printDesmos();
            delete bvh;

#if USE_AVX_TRI
            map<int, pair<BVH_AVX*, int>> order;
            for (int i = 0; i < PS->avx_bvh.size(); i++) {
                BVH_AVX& current = PS->avx_bvh[i];
                for (int i = 0; i < 8; i++) {
                    int leaf_size = current.leaf_size[i];
                    int index = current.indexes[i];
                    if (leaf_size > 0) {
                        order[index] = pair<BVH_AVX*, int>(&current, i);
                    }
                }
            }
            for (auto iter = order.begin(); iter != order.end(); ++iter)
            {
                int index = iter->first;
                auto datapoint = iter->second;
                BVH_AVX& current = *datapoint.first;
                int selection_index = datapoint.second;
                int leaf_size = current.leaf_size[selection_index];
                current.indexes[selection_index] = PS->avx_tri_data.size();
                current.leaf_size[selection_index] = ceil(((float)leaf_size) / 8.0);
                int i = 0;
                while (i < ceil(((float)leaf_size) / 8.0)) {
                    vector<PackagedTri> PTris;
                    for (int j = 0; j < 8 && i * 8 + j < leaf_size; j++) {
                        PTris.push_back(PS->tri_data[index + 8 * i + j]);
                    }
                    PS->avx_tri_data.push_back(PTri_AVX(PTris));
                    PS->avx_tri_data.back().index = index + 8 * i;
                    i++;
                }
            }
#endif
            auto BVH_collapse_end = high_resolution_clock::now();
            cout << "[Done][" << duration_cast<milliseconds>(BVH_collapse_end - BVH_collapse_start).count() << "ms]" << endl;
        }
        return PS;
    }
};
#define AVX_STACK_SIZE 256