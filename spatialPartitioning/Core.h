#pragma once
#include <cassert>
#include "commons.h"
#include "Scene.h"
#include "Ray.h"
#include "Coredefs.h"
#include "PDF.h"
#include <immintrin.h>
#define __SIMD_BITONIC_IMPLEMENTATION__
#include "simd_bitonic.h"

class ray_packet {
public:
    uint8_t ray_count = 0;
    m256_vec3 origins;
    m256_vec3 slopes;
    m256_vec3 coefficient;
    bool has_ray(int i) {
        return ray_count > i;
    }
    void add_ray(PackagedRay& PR) {
        origins.set(ray_count, PR.position);
        slopes.set(ray_count, PR.slope);
        coefficient.set(ray_count, PR.coefficient);
        ray_count++;
    }
};

unsigned int LOWER_BITMASK = 0x00000007u;
float MAX_NEG = std::bit_cast<float>(0xff7fffff);
const __m256 LOWER_BITMASK_AVX = _mm256_set1_ps(std::bit_cast<float>(LOWER_BITMASK));
const __m256 INDEX_ASCENDING = std::bit_cast<__m256>(_mm256_set_epi32(7,6,5,4,3,2,1,0));
const __m256 MAX_NEGATIVE = _mm256_set1_ps(MAX_NEG);


class RayEngine {
public:
    PackagedScene* data;

    static const int queue_size = 32;
    ray_packet queue[queue_size];
    int queue_bottom;
    int queue_top;

    RayEngine() {}
    void load_scene(PackagedScene* PS) {
        data = PS;
    }
    void iterativeBVH_AVX(CastResults& res, const XYZ& position, const XYZ& slope) {
        //global_counter++;
        m256_vec3 position_avx(position);
        m256_vec3 slope_avx(slope);
        m256_vec3 fusedposition(-1 * position / slope);
        m256_vec3 inv_slope(1 / slope);
        int stack[AVX_STACK_SIZE];
        float distances[AVX_STACK_SIZE];
        int tri_count[AVX_STACK_SIZE];
        int stack_index = 0;
        stack[0] = 0;
        distances[0] = 0;
        tri_count[0] = 0;
        while (stack_index >= 0) {
            const float dist = distances[stack_index];
            const int self_triangles = tri_count[stack_index];
            if (dist >= res.distance) {
                stack_index--;
                continue;
            };
            if (self_triangles == 0) {
                const BVH_AVX& current = data->avx_bvh[stack[stack_index]];
                stack_index--;
                auto m256_results = current.intersection(fusedposition, inv_slope);

                __m256 results_persist = m256_results;
                __m256 allow_mask = _mm256_cmp_ps(m256_results, AVX_ZEROS, _CMP_NEQ_OQ);
                m256_results = _mm256_andnot_ps(LOWER_BITMASK_AVX, m256_results);
                m256_results = _mm256_or_ps(m256_results, INDEX_ASCENDING);
                m256_results = _mm256_and_ps(m256_results, allow_mask);
                //float results[8];
                float* results = (float*)&m256_results;

                int sorted_size = 0;
                //m256_results = simd_sort_1V(m256_results);
                for (int i = 0; i < 8; i++) {
                    float num = results[i];
                    if (num != 0) {
                        int j = sorted_size - 1;
                        for (; j >= 0; j--) {
                            if (results[j] > num) {
                                results[j + 1] = results[j];
                            }
                            else break;
                        }
                        results[j + 1] = num;
                        sorted_size++;
                    }
                }
                for (int i =sorted_size-1; i >= 0; i--) {
                    //if (results[i] == MAX_NEG) break;
                    unsigned int selection_index = std::bit_cast<unsigned int>(results[i])&LOWER_BITMASK;
                    float dist2 = ((float*)&results_persist)[selection_index];
                    stack_index++;
                    stack[stack_index] = current.indexes[selection_index];
                    distances[stack_index] = dist2;
                    tri_count[stack_index] = current.leaf_size[selection_index];
                }


                //int sorted_index[8];
                //int sorted_index_size = 0;
                //for (int i = 0; i < 8; i++) {
                //    float& num = results[i];
                //    if (num != 0) {
                //        int j = sorted_index_size - 1;
                //        for (; j >= 0; j--) {
                //            if (results[sorted_index[j]] > num) {
                //                sorted_index[j + 1] = sorted_index[j];
                //            }
                //            else break;
                //        }
                //        sorted_index[j + 1] = i;
                //        sorted_index_size++;
                //    }
                //}
                //for (int i = sorted_index_size - 1; i >= 0; i--) {
                //    int selection_index = sorted_index[i];
                //    float dist2 = results[selection_index];
                //    stack_index++;
                //    stack[stack_index] = current.indexes[selection_index];
                //    distances[stack_index] = dist2;
                //    tri_count[stack_index] = current.leaf_size[selection_index];
                //
                //}

                //stack_index++; places the bvhs in without order. Slower overall, but the simpler method makes it still viable.
                //for (int i = 0; i < 8; i++) {
                //    bool sgnbit = results[i] != 0;
                //    stack[stack_index]      = current.indexes[i];
                //    tri_count[stack_index]  = current.leaf_size[i];
                //    distances[stack_index] = results[i];
                //    stack_index+=sgnbit;
                //}
                //stack_index--;

                //__m256 temp_m256 = simd_sort_1V(_mm256_max_ps(m256_results,AVX_ZEROS));
                //float* temp = (float*) & temp_m256;
                //
                //float values[8];
                //int pointer[8];
                //int size = 0;
                //for (int i = 0; i < 8; i++) {
                //    if (results[i] > 0) {
                //        values[size] = results[i];
                //        pointer[size] = i;
                //        size++;
                //    }
                //}
                //
                //for (int i = 7; i > 7-size; i--) {
                //    float dist2 = temp[i];
                //    int index = 0;
                //    for (int k = size-1; k >= 0; k--) {
                //        if (values[k] == dist2) {
                //            index = pointer[k];
                //            values[k] = -1;
                //            break;
                //        }
                //    }
                //    stack_index++;
                //    stack[stack_index] = current.indexes[index];
                //    distances[stack_index] = dist2;
                //    tri_count[stack_index] = current.leaf_size[index];
                //
                //}
            }
            else {
                int start_index = stack[stack_index];
                stack_index--;
                for (int i = 0; i < self_triangles; i++) {
#if !USE_AVX_TRI

                    const PackagedTri& t = (data->tri_data)[i + start_index];
                    XYZ distance = Tri::intersection_check_MoTru(t, position, slope);
                    if (distance.X >= 0) {
                        if (distance.X < res.distance) {
                            res.distance = distance.X;
                            res.normal = Tri::get_normal(t.normal, position - t.origin_offset);
                            res.material = t.material;
                            res.UV = distance;
                        }
                    }
#else

                    const PTri_AVX& PTri_AVX_pack = (data->avx_tri_data)[i + start_index];
                    m256_vec3 intersection_stats;
                    //m256_vec3 normal_stats;
                    //vector<PackagedTri> pTris;
                    //vector<XYZ> real_res;
                    //vector<XYZ> real_normal;
                    //for (int l = 0; l < 8; l++) {
                    //    PackagedTri pTri = PackagedTri();
                    //    pTri.p1 = PTri_AVX_pack.p1.at(l);
                    //    pTri.p1p2 = PTri_AVX_pack.p1p2.at(l);
                    //    pTri.p1p3 = PTri_AVX_pack.p1p3.at(l);
                    //    pTri.normal = PTri_AVX_pack.normal.at(l);
                    //    pTri.origin_offset = PTri_AVX_pack.origin_offset.at(l);
                    //    pTri.UV_basis = PTri_AVX_pack.UV_basis.at(l);
                    //    pTri.U_delta = PTri_AVX_pack.U_delta.at(l);
                    //    pTri.V_delta = PTri_AVX_pack.V_delta.at(l);
                    //    pTri.material = PTri_AVX_pack.materials[l];
                    //    real_res.push_back(Tri::intersection_check(pTri, position, slope));
                    //    real_normal.push_back(Tri::get_normal(pTri.normal, position - pTri.origin_offset));
                    //    pTris.push_back(pTri);
                    //}
                    //
                    //PTri_AVX::intersection_check(PTri_AVX_pack, position_avx, slope_avx, intersection_stats);
                    //PTri_AVX::get_normal(PTri_AVX_pack, position_avx, normal_stats);
                    //
                    //float* dists = (float*)&intersection_stats.X;
                    //for (int i = 0; i < 8; i++) {
                    //    if (dists[i] > 0 && dists[i] < res.distance) {
                    //        res.distance = dists[i];
                    //        res.normal = normal_stats.at(i);
                    //        res.material = PTri_AVX_pack.materials[i];
                    //        res.UV = intersection_stats.at(i);
                    //        res.tri_index = PTri_AVX_pack.get_index(i);
                    //    }
                    //}
                    PTri_AVX::intersection_check(PTri_AVX_pack, position_avx, slope_avx, intersection_stats);

                    float* dists = (float*)&intersection_stats.X;
                    for (int i = 0; i < 8; i++) {
                        if (dists[i] > 0 && dists[i] < res.distance) {
                            res.distance = dists[i];
                            res.material = PTri_AVX_pack.materials[i];
                            res.TUV = intersection_stats.at(i);
                            res.tri_index = PTri_AVX_pack.get_index(i);
                        }
                    }
#endif
                }
            }
        }
    }
    void iterativeBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
        int stack[64];
        float dist[64];
        int stack_index = 0;
        stack[0] = 0;
        while (stack_index >= 0) {
            const PackagedBVH& current = data->flat_bvh[stack[stack_index]];
            stack_index--;
            if (dist[stack_index + 1] > res.distance) continue;
            //if (distance <= 0 || distance >= res.distance) {
            //    continue;
            //}
            if (current.leaf_size == 0) {
                int c1_index = current.index;
                int c2_index = c1_index + 1;
                const PackagedBVH& c1 = data->flat_bvh[c1_index];
                const PackagedBVH& c2 = data->flat_bvh[c2_index];
                float t1 = BVH::intersection(c1.sMax, c1.sMin, position, inv_slope);
                float t2 = BVH::intersection(c2.sMax, c2.sMin, position, inv_slope);

                bool f1 = t1 < 0;//t1 < 0 || t1>res.distance;
                bool f2 = t2 < 0;//(t2 < 0 || t2>res.distance);

                if (f1 && f2) {
                    continue;
                }
                if (f1) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                    continue;
                }
                if (f2) {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                    continue;
                }
                if (t1 > t2 && abs(t1 - t2) > 0.1) {
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                }
                else {
                    stack_index++;
                    stack[stack_index] = c1_index;
                    dist[stack_index] = t1;
                    stack_index++;
                    stack[stack_index] = c2_index;
                    dist[stack_index] = t2;
                }
            }
            else {
                for (int i = 0; i < current.leaf_size; i++) {
                    const PackagedTri& t = (data->tri_data)[i + current.index];
                    XYZ distance = Tri::intersection_check(t, position, slope);
                    if (distance.X >= 0) {
                        if (distance.X < res.distance) {
                            res.distance = distance.X;
                            res.material = t.material;
                            res.TUV = distance;
                        }
                    }
                }
            }
        }
    }
    void navBVH(CastResults& res, const XYZ& position, const XYZ& slope, const XYZ& inv_slope) {
#if USE_AVX_BVH
        //iterativeBVH_AVX_unordered(starting_index, res, position, slope);
        iterativeBVH_AVX(res, position, slope);
#else
        iterativeBVH(res, position, slope, inv_slope);
#endif
    }
    CastResults execute_bvh_cast(XYZ& position, const XYZ& slope) {
        CastResults returner;
        XYZ inv_slope = XYZ(1) / slope;
        navBVH(returner, position, slope, inv_slope);
        position += returner.distance * slope * 0.999;
        return returner;
    }
    CastResults execute_naive_cast(XYZ& position, const XYZ& slope) {
#define default_smallest_distance 9999999;
        float smallest_distance = default_smallest_distance;
        CastResults returner = CastResults(-1);
        for (const PackagedTri& t : data->tri_data) {
            XYZ distance = Tri::intersection_check(t, position, slope);
            if (distance.X >= 0) {
                if (distance.X < returner.distance) {
                    returner.distance = distance.X;
                    returner.material = t.material;
                }
            }
        }
        position += returner.distance * slope;
        return returner;
    }
    CastResults execute_ray_cast(XYZ& position, const XYZ& slope) {
        return execute_bvh_cast(position, slope);
        //return execute_naive_cast(position,slope);

    }
    bool execute_lighting_cast(XYZ position, Light* L) {
        XYZ slope = L->vectorTo(position);
        auto res = execute_bvh_cast(position, slope);
        if (L->intersects(position, slope, res)) {
            return true;
        }
        return false;
    }
    XYZ process_ray(Casting_Diagnostics& stats, PackagedRay& ray_data) {
        stats.rays_processed++;

        XYZ orig = ray_data.position;
        CastResults results = execute_ray_cast(ray_data.position, ray_data.slope);
        PackagedTri& Tri = data->tri_data[results.tri_index];
        results.texUV.X = results.TUV.Y * Tri.U_delta.X + results.TUV.Z * Tri.V_delta.X + Tri.UV_basis.X;
        results.texUV.Y = results.TUV.Y * Tri.U_delta.Y + results.TUV.Z * Tri.V_delta.Y + Tri.UV_basis.Y;

        if (results.material == -1) {
            return XYZ(0);
        }
        Material& parent_material = data->materials[results.material];
        MaterialSample material = parent_material.sample_UV(results.texUV.X, results.texUV.Y);

        bool backfacing = XYZ::dot(Tri.normal, orig - Tri.origin_offset) < 0;
        float facing_sign = backfacing ? -1 : 1;

        XYZ geo_normal = backfacing ? -Tri.normal : Tri.normal;
        XYZ smoothed_normal_raw = Tri.get_normal(results.TUV.Y, results.TUV.Z);

        // this is kind of ass, but ensures that the tris smoothed normal is aligned to the facing of the normal
        //in theory this should be flipped by the backfacing attribute, but severe smooth shading can break that
        //this doesnt break under that circumstances
        XYZ normal = (XYZ::dot(geo_normal, smoothed_normal_raw) < 0) ? -smoothed_normal_raw : smoothed_normal_raw;
        XYZ T = Tri.get_tangent(results.TUV.Y, results.TUV.Z);
        XYZ BT = Tri.get_bitangent(results.TUV.Y, results.TUV.Z);
        if (parent_material.use_normals) {
            normal = material.normalPertubation.X * T + material.normalPertubation.Y * BT + material.normalPertubation.Z * normal;
        }

        normal = XYZ::normalize(normal);

        XYZ aggregate = XYZ();
        {
            if (DRAW_UV) {
                if (results.TUV.X != -1) {
                    XYZ color;
                    color.X = results.texUV.X;
                    color.Y = results.texUV.Y;
                    color.Z = 1 - color.X - color.Y;
                    return color * 1;
                }
            }
            if (DRAW_NORMAL) {
                XYZ res = (normal);
                return (normal / 2 + XYZ(0.5, 0.5, 0.5)) * 1;

            }
            if (DRAW_NORMAL_DELTAS) {
                XYZ res = (geo_normal - normal);
                if (res.magnitude() < 0.1) {
                    return XYZ(0, 0, 0);
                }
                return (XYZ::normalize(res) / 2 + XYZ(0.5, 0.5, 0.5)) * 1;

            }
            if (DRAW_NORMAL_FACING) return XYZ::dot(normal, orig) / abs(XYZ::dot(normal, orig));
            if (DRAW_COLOR) return material.color;
            if (false) {
                if (ray_data.generation == 1)
                    aggregate += ray_data.position / 2 + XYZ(0.5, 0.5, 0.5);
            }
            //if (DRAW_EMISSIVE) return material.calculate_emissions();
        }


        aggregate += material.calculate_emissions();

        int bounce_count = data->monte_carlo_max * pow(data->monte_carlo_modifier, ray_data.generation);
        XYZ flipped_output = XYZ::flip(ray_data.slope);
        XYZ reflection_slope = XYZ::normalize(XYZ::reflect(XYZ::flip(ray_data.slope), normal));

        float dot_NO = XYZ::dot(normal, flipped_output);

        short n1_mId = -1;
        if (ray_data.IOR_size > 0)
            n1_mId = ray_data.IOR_stack[ray_data.IOR_size - 1];
        float n1 = 1.0003;
        if (n1_mId != -1) data->materials[n1_mId].IOR.getValue(0, 0);
        float ior_u = n1 / material.IOR;

        uint16_t trans_stack[IOR_STACK_SIZE];
        int trans_stack_size = ray_data.IOR_size;
        for (int i = 0; i < ray_data.IOR_size;i++) trans_stack[i] = ray_data.IOR_stack[i];
        if (!backfacing) {
            trans_stack[ray_data.IOR_size] = results.material;
            trans_stack_size++;
        } else {
            int i = 0;
            for(; i < ray_data.IOR_size;i++) if (trans_stack[i] == results.material) break;
            for (; i < ray_data.IOR_size-1;i++) trans_stack[i] = trans_stack[i + 1];
            trans_stack_size--;
        }

        material.config(flipped_output, normal, n1, material.IOR);

#if DRAW_LIGHTS
        for (Light* L : data->lights) {
            XYZ pos = ray_data.position + 0.001 * geo_normal;
            XYZ slope = L->vectorTo(pos);
            material.set_input(slope);
            XYZ return_coefficient = material.fast_BRDF_co(true);
            if (execute_lighting_cast(pos, L)) {
                aggregate += return_coefficient * L->emission;
#if DRAW_LIGHT_EXPOSURE
                return aggregate;
            }
            return XYZ();
#else
            }
#endif
        }
#endif
        if (ray_data.generation >= data->max_generations) return aggregate;
        if (bounce_count < 1) {
            bounce_count = 1;
        }
        //if (diffuse_bounces < 4 && diffuse_bounces > 0) diffuse_bounces = 4;
        //int diffuse_Vslices = floor(log2(diffuse_bounces) - 1);
        //int diffuse_Rslices = floor(diffuse_bounces / diffuse_Vslices);
        //diffuse_bounces = diffuse_Vslices * diffuse_Rslices;
        //double diffuse_multiplier = 1.0 / bounce_count;
        //double specular_multiplier = 1.0 / bounce_count;
        //float diffuse_V_increment = 1.0 / diffuse_Vslices;
        //float diffuse_R_increment = 1.0 / diffuse_Rslices;
        //float diffuse_V_position = 0;
        //float diffuse_R_position = 0;

        int ri = 0;
        if (XYZ::dot(flipped_output, normal) < 0) {
            return aggregate;
        }

        vector<PackagedRay> r_queue;
        r_queue.reserve(bounce_count);

        XYZ throughput = XYZ(1, 1, 1);

        for (int i = 0; i < bounce_count; i++) {
            ri++;
            gen.regen_pool();
            if (ri > bounce_count * 10 && i == 0) {
                return aggregate;
            }
            //if (ri > bounce_count * 2) {
            //    cout << ">:(";
            //}
            float factor = gen.pool_rand(0);//gen.fRand(0, 1);
            float c_theta = PDF::sample(factor, material.i_roughness); //theta of one side of the cavity
            float azimuth = gen.pool_rand(1) * 2.0f * fPI;//gen.fRand(0, 2 * PI);
            float c_theta_sin = VecLib::full_sin(c_theta + fPI);
            float c_theta_cos = VecLib::full_cos(c_theta + fPI);
            XYZ MN_1 = XYZ(-c_theta_sin * VecLib::full_cos(azimuth), -c_theta_sin * VecLib::full_sin(azimuth), -c_theta_cos); //micro normal 1, and its flipped counterpart
            MN_1 = Matrix3x3::applyTBNMatrix(MN_1, T, BT, normal);
            XYZ MN_2 = XYZ::reflect(MN_1, normal);
            float choice_prob = XYZ::cdot(flipped_output, MN_2) / (XYZ::cdot(flipped_output, MN_1) + XYZ::cdot(flipped_output, MN_2));
            float choice_factor = gen.pool_rand(2);//gen.fRand(0, 1);
            XYZ MN = MN_1;
            if (choice_factor < choice_prob) {
                MN = MN_2;
            }
            XYZ outgoing;
            float ref_prob = material.fast_fresnel(MN, flipped_output);
            float ref_factor = gen.pool_rand(3);//gen.fRand(0, 1);
            XYZ slope;
            XYZ return_coefficient = throughput;
            bool is_transmissive = false;
            if (ref_factor < ref_prob) {
                slope = XYZ::reflect(flipped_output, MN);
                material.set_input(slope);
                return_coefficient *= material.PDF_specular_BRDF();
            }
            else {
                float trans_factor = gen.pool_rand(4);
                if (trans_factor < material.transmission) {
                    float dot_HI = XYZ::dot(MN, flipped_output);
                    float f1 = sqrt(1 - ior_u * ior_u * (1 - dot_HI * dot_HI));
                    slope = -(f1 * MN + ior_u * (flipped_output - dot_HI * MN));
                    is_transmissive = true;
                }
                else {
                    slope = VecLib::biased_random_hemi(gen, 0, 1, 0, 1);
                    slope = { slope.X, slope.Z, slope.Y };
                    slope = Matrix3x3::applyTBNMatrix(slope, T, BT, normal);
                    //slope = material.biased_diffuse_bounce(gen, normal_rot_m);
                    material.set_input(slope);
                    return_coefficient *= material.diffuse_BRDF();
                }
            }
            bool t_flip = is_transmissive;
            if (t_flip != (XYZ::dot(slope, normal) < 0)) {
                throughput = return_coefficient;
                i--;
                continue;
            }
            if (t_flip != (XYZ::dot(slope, geo_normal) < 0)) {
                throughput = XYZ(1, 1, 1);
                i--;
                continue;
            }
            if (return_coefficient.X < 0 || return_coefficient.Y < 0 || return_coefficient.Z < 0) {
                throughput = XYZ(1, 1, 1);
                i--;
                continue;
            }
            XYZ global_return = return_coefficient;
            float russian_p = std::max(global_return.X, std::max(global_return.Y, global_return.Z));
            float russian_factor = gen.pool_rand(5);//gen.fRand(0, 1);
            if (russian_factor > russian_p) {
                throughput = XYZ(1, 1, 1);
                continue;
            }
            return_coefficient = return_coefficient / russian_p;
            auto ray = PackagedRay(
                ray_data.position + (is_transmissive?-1:1)*0.0001 * geo_normal,
                slope,
                ray_data.generation+1,
                return_coefficient
            );
            ray.output = ray_data.output;
            if (is_transmissive) {
                 for (int i = 0; i < trans_stack_size; i++) ray.IOR_stack[i] = trans_stack[i];
                 ray.IOR_size = trans_stack_size;
            }
            else {
                for (int i = 0; i < IOR_STACK_SIZE; i++) ray.IOR_stack[i] = ray_data.IOR_stack[i];
                ray.IOR_size = ray_data.IOR_size;
            }
            r_queue.push_back(ray);
            throughput = XYZ(1, 1, 1);

        }
        for (int i = 0; i < r_queue.size(); i++) {
            PackagedRay& ray = r_queue[i];
            ray.coefficient = ray.coefficient/bounce_count;
            stats.diffuses_cast++;
            XYZ returned_light = process_ray(stats, ray);
            aggregate += ray.coefficient * returned_light;
        }

        return aggregate;
    }
    XYZ process_top_packet(Casting_Diagnostics& stats) {
        ray_packet& packet = queue[queue_bottom];
        queue_bottom = (queue_bottom+1)%queue_size;
        stats.rays_processed += packet.ray_count;
        m256_vec3 orig = packet.origins;
        CastResults results[8];
        PackagedTri* tris[8];
        Material* materials[8];
        for (int i = 0; i < packet.ray_count; i++) {
            XYZ ray_origin = packet.origins.at(i);
            XYZ ray_slope = packet.slopes.at(i);
            results[i] = execute_ray_cast(ray_origin, ray_slope);
            packet.origins.set(i, ray_origin);
            _mm_prefetch((char*)&(data->tri_data[results[i].tri_index]), _MM_HINT_T1); //get tri fetching
            if (i > 0) _mm_prefetch(
                (char*)&(data->materials[data->tri_data[results[i-1].tri_index].material]),
                _MM_HINT_T1
            ); //now that the previous cycles tri is fetched, get its material
        }
        m256_vec3 aggregate;
        aggregate.X = _mm256_set1_ps(0);
        aggregate.Y = _mm256_set1_ps(0);
        aggregate.Z = _mm256_set1_ps(0);
        __m256 UV_tU;
        __m256 UV_tV;
        __m256 UV_tK;
        __m256 tex_UV_U;
        __m256 tex_UV_V;
        m256_vec2 UV_basis;
        m256_vec2 U_delta;
        m256_vec2 V_delta;
        m256_vec3 origin_offset;
        m256_vec3 n1;
        m256_vec3 n2;
        m256_vec3 n3;
        m256_vec3 t1;
        m256_vec3 t2;
        m256_vec3 t3;
        m256_vec3 b1;
        m256_vec3 b2;
        m256_vec3 b3;
        m256_vec3 normals;
        m256_vec3 bitangents;
        m256_vec3 tangents;
        m256_vec3 geo_normals;
        for (int i = 0; i < packet.ray_count; i++) {
            tris[i] = &(data->tri_data[results[i].tri_index]);
            ((float*)&UV_tU)[i] = results[i].TUV.Y;
            ((float*)&UV_tV)[i] = results[i].TUV.Z;
            UV_basis.set(i,tris[i]->UV_basis);
            U_delta.set(i, tris[i]->UV_basis);
            V_delta.set(i, tris[i]->UV_basis);
            geo_normals.set(i, tris[i]->normal);
            origin_offset.set(i, tris[i]->origin_offset);
            n1.set(i, tris[i]->n1);
            n2.set(i, tris[i]->n2);
            n3.set(i, tris[i]->n3);
            t1.set(i, tris[i]->t1);
            t2.set(i, tris[i]->t2);
            t3.set(i, tris[i]->t3);
            b1.set(i, tris[i]->b1);
            b2.set(i, tris[i]->b2);
            b3.set(i, tris[i]->b3);
        }
        tex_UV_U = _mm256_fmadd_ps(
            UV_tU,
            U_delta.X,
            _mm256_fmadd_ps(
                UV_tV,
                V_delta.X,
                UV_basis.X
            )
        );
        tex_UV_V = _mm256_fmadd_ps(
            UV_tU,
            U_delta.Y,
            _mm256_fmadd_ps(
                UV_tV,
                V_delta.Y,
                UV_basis.Y
            )
        );
        UV_tK = _mm256_sub_ps(
            AVX_ONES,
            _mm256_sub_ps(
                UV_tU,
                UV_tV
            )
        );

        AVX_MaterialSample materialSample(materials, (float*)&tex_UV_U, (float*)&tex_UV_V);

        orig.sub(orig, origin_offset, orig);
        __m256 dot_outputs;
        VecLib::dot_avx(geo_normals, orig, dot_outputs);
        //this has a 1 in the sign bit if a vector needs to be flipped.
        __m256 flip_mask = _mm256_and_ps(AVX_SIGN_MASK,
            _mm256_cmp_ps(dot_outputs, AVX_ZEROS, _CMP_LT_OQ)); 
        geo_normals = m256_vec3::bitwise_xor(geo_normals, flip_mask);
        normals = m256_vec3::bitwise_xor( //flip according to flip mask
            m256_vec3::normalize(
                m256_vec3::muladd( //calculate smoothed normal
                    n1,
                    UV_tK,
                    m256_vec3::muladd(
                        n2,
                        UV_tU,
                        m256_vec3::mul(
                            n3,
                            UV_tV
                        )))),
            flip_mask
        );
        tangents = m256_vec3::bitwise_xor( //flip according to flip mask
            m256_vec3::normalize(
                m256_vec3::muladd( //calculate smoothed normal
                    t1,
                    UV_tK,
                    m256_vec3::muladd(
                        t2,
                        UV_tU,
                        m256_vec3::mul(
                            t3,
                            UV_tV
                        )))),
            flip_mask
        );
        bitangents = m256_vec3::bitwise_xor( //flip according to flip mask
            m256_vec3::normalize(
                m256_vec3::muladd( //calculate smoothed normal
                    b1,
                    UV_tK,
                    m256_vec3::muladd(
                        b2,
                        UV_tU,
                        m256_vec3::mul(
                            b3,
                            UV_tV
                        )))),
            flip_mask
        );
        normals = m256_vec3::bitwise_or(
            m256_vec3::bitwise_and(
                m256_vec3::muladd(
                    tangents,
                    materialSample.normalPertubation.X,
                    m256_vec3::muladd(
                        bitangents,
                        materialSample.normalPertubation.Y,
                        m256_vec3::mul(
                            normals,
                            materialSample.normalPertubation.Z
                        )
                    )
                ),
                materialSample.use_normal
            ),
            m256_vec3::bitwise_and(
                normals,
                _mm256_xor_ps(AVX_ZEROS, materialSample.use_normal)
            )
        );//and we FINALLY have normal mapped normals. Christ.
        aggregate = m256_vec3::mul(
            packet.coefficient,
            materialSample.emissive
        );
        __m256 max_coefficient = _mm256_max_ps(
            packet.coefficient.X,
            _mm256_max_ps(
                packet.coefficient.Y,
                packet.coefficient.Z
            )
        );
        __m256i bounce_count = _mm256_cvtps_epi32(
            _mm256_ceil_ps(
                _mm256_mul_ps(
                    _mm256_sqrt_ps(
                        max_coefficient
                    ),
                    _mm256_set1_ps(data->monte_carlo_max)
                )
            )
        );

        m256_vec3 flipped_output = m256_vec3::bitwise_xor(packet.slopes, AVX_SIGN_MASK);

        tangents = m256_vec3::normalize(
            m256_vec3::sub(tangents, m256_vec3::mul(normals, m256_vec3::dot(tangents, normals))));

        bitangents = m256_vec3::cross(normals, tangents);

        __m256 dot_NO = m256_vec3::dot(flipped_output, normals);

        __m256 qualified_mask = _mm256_cmp_ps(
            dot_NO,
            AVX_ZEROS,
            _CMP_GT_OQ
        );

        bounce_count = _mm256_or_epi32(
            std::bit_cast<__m256i>(_mm256_and_ps(
                AVX_ONES,
                qualified_mask
            )),
            _mm256_and_epi32(
                bounce_count,
                std::bit_cast<__m256i>(_mm256_xor_ps(
                    qualified_mask,
                    AVX_ONES
                ))
            )
        );
        __m256i bounces_remaining = bounce_count;
        int* bounce_lookup = (int*)&bounces_remaining;
        int* total_bounce_lookup = (int*)&bounce_count;
        int true_index[8] = { 0,1,2,3,4,5,6,7 };
        int rebalance_in = 0;
        m256_vec3 return_coefficient;
        __m256 allow_bounce;
        int iterations = -1;
        while (true) {
            iterations++;
            if (rebalance_in == 0) {
                int minimum_bounces = 9999;
                bool continue_processing = false;
                int max_index = 0;
                int max_bounces = 0;
                for (int i = 0; i < 8;i++) {
                    int current_remainder = bounce_lookup[true_index[i]];
                    if (current_remainder > 0 && iterations < total_bounce_lookup[i]*10) {
                        continue_processing = true;
                        minimum_bounces = std::min(minimum_bounces, current_remainder);
                        if (current_remainder > max_bounces) {
                            max_bounces = current_remainder;
                            max_index = i;
                        }
                    }
                }
                if (!continue_processing) break;
                for (int i = 0; i < 8;i++) {
                    int current_remainder = bounce_lookup[true_index[i]];
                    if (current_remainder == 0 && max_bounces > 0) {
                        true_index[i] = max_index;
                        packet.coefficient.transfer(max_index, i);
                        packet.origins.transfer(max_index, i);
                        flipped_output.transfer(max_index, i);
                        materialSample.transfer(max_index, i);
                        ((float*)&dot_NO)[i] = ((float*)&dot_NO)[max_index];
                    }
                }
            }
            gen.regen_pool();
            __m256 azimuth = _mm256_mul_ps(gen.pool, AVX_2PI);
            gen.regen_pool();
            __m256 c_thetas;
            for (int i = 0; i < 8; i++) {
                float r = ((float*)&materialSample.roughness)[i];
                float factor = gen.pool_rand(i);
                ((float*)&c_thetas)[i] = PDF::sample(factor, r*0.001);
            }
            __m256 sin_c_theta = VecLib::full_sin_AVX(c_thetas);
            m256_vec3 MN_1;
            //performing a trig identity is probably faster for cosine since we need sine anyways.
            MN_1.X = _mm256_mul_ps(sin_c_theta, VecLib::full_cos_AVX(azimuth));
            MN_1.Y = _mm256_mul_ps(sin_c_theta, VecLib::full_sin_AVX(azimuth));
            MN_1.Z = VecLib::full_cos_AVX(c_thetas);
            MN_1 = m256_vec3::apply_matrix(MN_1, tangents, bitangents, normals); //theoretically transforms into normal space

        }

    }
};
