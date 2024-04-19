#pragma once
#include <cassert>
#include "commons.h"
#include "Scene.h"
#include "Ray.h"
#include "Coredefs.h"
#include "PDF.h"


class RayEngine {
public:
    PackagedScene* data;
    vector<vector<PackagedRay>> ray_packet_queue;
    const int max_queue = pow(2, 20); //bit over a million
    const int sub_packet_size = 1024;
    RayEngine() {}
    void load_scene(PackagedScene* PS) {
        data = PS;
    }
    void traverse_packets() {
        vector<PackagedRay>& packet = ray_packet_queue.front();
        vector<MinimumRay> minimal_sub_packet;
        vector<AVX_AABB_ray> sub_packet;
        sub_packet.resize(sub_packet_size);
        minimal_sub_packet.resize(sub_packet_size);
        unsigned int packet_ray_index = 0;
        while (packet_ray_index < packet.size()) {
            for (unsigned int i = 0; i < sub_packet_size; i++) {
                auto& ray_data = packet[packet_ray_index];
                minimal_sub_packet[i] = MinimumRay{ ray_data.position,1 / ray_data.slope,packet_ray_index };
                sub_packet[i] = AVX_AABB_ray{ m256_vec3(-1 * ray_data.position / ray_data.slope),m256_vec3(1 / ray_data.slope),i };
                packet_ray_index++;
            }
        }


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


                //float results[8];
                float* results = (float*)&m256_results;
                int sorted_index[8];
                int sorted_index_size = 0;
                for (int i = 0; i < 8; i++) {
                    float& num = results[i];
                    if (num != 0) {
                        int j = sorted_index_size - 1;
                        for (; j >= 0; j--) {
                            if (results[sorted_index[j]] > num) {
                                sorted_index[j + 1] = sorted_index[j];
                            }
                            else break;
                        }
                        sorted_index[j + 1] = i;
                        sorted_index_size++;
                    }
                }
                for (int i = sorted_index_size - 1; i >= 0; i--) {
                    int selection_index = sorted_index[i];
                    float dist2 = results[selection_index];
                    stack_index++;
                    stack[stack_index] = current.indexes[selection_index];
                    distances[stack_index] = dist2;
                    tri_count[stack_index] = current.leaf_size[selection_index];

                }

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
                            res.normal = Tri::get_normal(PTri_AVX_pack.normal[i], position - PTri_AVX_pack.origin_offset[i]);
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
                            res.normal = Tri::get_normal(t.normal, position - t.origin_offset);
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
        CastResults returner = CastResults(XYZ(-1, -1, -1), nullptr);
        for (const PackagedTri& t : data->tri_data) {
            XYZ distance = Tri::intersection_check(t, position, slope);
            if (distance.X >= 0) {
                if (distance.X < returner.distance) {
                    returner.distance = distance.X;
                    returner.normal = Tri::get_normal(t.normal, position - t.origin_offset);
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

        if (results.material == nullptr) {
            return XYZ(0);
        }
        MaterialSample material = results.material->sample_UV(results.texUV.X, results.texUV.Y);
        XYZ aggregate = XYZ();
        XYZ geo_normal = Tri::get_normal(Tri.normal, orig - Tri.origin_offset);
        XYZ smoothed_normal_raw = Tri.get_normal(results.TUV.Y, results.TUV.Z);
        XYZ normal = (XYZ::dot(geo_normal, smoothed_normal_raw) < 0)? -smoothed_normal_raw : smoothed_normal_raw;
        if (results.material->use_normals) {
            XYZ T = Tri.get_tangent(results.TUV.Y, results.TUV.Z);
            XYZ BT = Tri.get_bitangent(results.TUV.Y, results.TUV.Z);
            normal = material.normalPertubation.X * T + material.normalPertubation.Y * BT + material.normalPertubation.Z * normal;

        }
        normal = XYZ::normalize(normal);
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
                XYZ res = (results.normal - normal);
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

        Quat normal_rot = Quat::makeRotationFromY(normal);
        Quat reflection_rot = Quat::makeRotationFromY(reflection_slope);
        Matrix3x3 normal_rot_m = Matrix3x3::quatToMatrix(normal_rot);
        Matrix3x3 reflection_rot_m = Matrix3x3::quatToMatrix(reflection_rot);
        //ray_data.PreLL.value = XYZ(0, 0, 100);

        float dot_NO = XYZ::dot(normal, flipped_output);
        material.config(flipped_output, normal, 1.0003, material.IOR);

#if DRAW_LIGHTS
        for (Light* L : data->lights) {
            XYZ pos = ray_data.position + 0.001 * results.normal;
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
        if (ray_data.generation < data->monte_carlo_generations) {


            int diffuse_bounces = bounce_count;
            int specular_bounces = 0;
            for (int i = 0; i < bounce_count; i++) {
                if (gen.fRand(0, 1) > material.f_spec) {
                    specular_bounces++;
                    diffuse_bounces--;
                }
            }
#if !DRAW_DIFFUSE
            diffuse_bounces = 0;
            specular_bounces = bounce_count;
#endif
#if !DRAW_SPECULAR
            diffuse_bounces = bounce_count;
            specular_bounces = 0;
#endif
#if !DRAW_SPECULAR && !DRAW_DIFFUSE
            diffuse_bounces = 0;
            specular_bounces = 0;
#endif
            if (diffuse_bounces < 4 && diffuse_bounces > 0) diffuse_bounces = 4;
            int diffuse_Vslices = floor(log2(diffuse_bounces) - 1);
            int diffuse_Rslices = floor(diffuse_bounces / diffuse_Vslices);
            diffuse_bounces = diffuse_Vslices * diffuse_Rslices;

            double diffuse_multiplier = 1.0 / bounce_count;
            double specular_multiplier = 1.0 / bounce_count;



            float diffuse_V_increment = 1.0 / diffuse_Vslices;
            float diffuse_R_increment = 1.0 / diffuse_Rslices;
            float diffuse_V_position = 0;
            float diffuse_R_position = 0;

            int ri = 0;
            if (XYZ::dot(flipped_output, results.normal) < 0) {
                cout << ">:(";
            }
            if (XYZ::dot(flipped_output, normal) < 0) {
                return aggregate;
            }
            
            vector<PackagedRay> r_queue;

            for (int i = 0; i < bounce_count; i++) {
                ri++;
                if (ri > bounce_count * 10 && i == 0) {
                    return aggregate;
                }
                //if (ri > bounce_count * 2) {
                //    cout << ">:(";
                //}
                float factor = gen.fRand(0, 1);
                float c_theta = PDF::sample(factor, material.i_roughness); //theta of one side of the cavity
                float azimuth = gen.fRand(0, 2 * PI);
                XYZ MN_1 = XYZ(sin(c_theta) * cos(azimuth), cos(c_theta), sin(c_theta) * sin(azimuth)); //micro normal 1, and its flipped counterpart
                MN_1 = Matrix3x3::applyRotationMatrix(MN_1, normal_rot_m);
                XYZ MN_2 = XYZ::reflect(MN_1, normal);
                float choice_prob = XYZ::cdot(flipped_output, MN_2) / (XYZ::cdot(flipped_output, MN_1) + XYZ::cdot(flipped_output, MN_2));
                float choice_factor = gen.fRand(0, 1);
                XYZ MN = MN_1;
                if (choice_factor < choice_prob) {
                    MN = MN_2;
                }
                XYZ outgoing;
                float ref_prob = material.fast_fresnel(MN, flipped_output);
                float ref_factor = gen.fRand(0, 1);
                XYZ slope;
                XYZ return_coefficient;
                if (ref_factor < ref_prob) {
                    slope = XYZ::reflect(flipped_output, MN);
                    material.set_input(slope);
                    return_coefficient = material.PDF_specular_BRDF();
                }
                else {
                    slope = material.biased_diffuse_bounce(gen,normal_rot_m);
                    material.set_input(slope);
                    return_coefficient = material.diffuse_BRDF();
                }

                if (XYZ::dot(slope, geo_normal) < 0 || XYZ::dot(slope, normal) < 0) {
                    i--;
                    continue;
                }
                if (return_coefficient.X < 0 || return_coefficient.Y < 0 || return_coefficient.Z < 0) {
                    i--;
                    continue;
                }
                XYZ global_return = return_coefficient;
                float russian_p = std::max(global_return.X, std::max(global_return.Y, global_return.Z));
                float russian_factor = gen.fRand(0, 1);
                if (russian_factor > russian_p) {
                    continue;
                }
                return_coefficient = return_coefficient/russian_p;
                auto ray = PackagedRay(
                    ray_data.position + 0.001 * results.normal,
                    slope,
                    ray_data.generation + 1,
                    return_coefficient
                );
                ray.output = ray_data.output;
                r_queue.push_back(ray);

            }
            for (int i = 0; i < r_queue.size(); i++) {
                PackagedRay& ray = r_queue[i];
                stats.diffuses_cast++;
                XYZ returned_light = process_ray(stats, ray);
                aggregate += ray.coefficient/bounce_count * returned_light;
            }

            /*
            for (int i = 0; i < diffuse_bounces; i++) {
                if (i % diffuse_Rslices == 0 && i != 0) {
                    diffuse_V_position += diffuse_V_increment;
                }
                diffuse_R_position += diffuse_R_increment;
                if (diffuse_R_position > 1) diffuse_R_position -= 1;
                XYZ diffuse_slope = material.biased_diffuse_bounce(
                    gen,
                    normal_rot_m
                    , diffuse_V_position
                    , diffuse_V_position + diffuse_V_increment
                    , diffuse_R_position
                    , diffuse_R_position + diffuse_R_increment
                );
                //XYZ diffuse_slope = material.biased_diffuse_bounce(
                //    normal_rot_m
                //);
                if (DRAW_BOUNCE_DIRECTION) {
                    if (ray_data.generation == 0) {
                        aggregate += diffuse_multiplier * (diffuse_slope / 2 + XYZ(0.5, 0.5, 0.5)) * 1;
                        continue;
                    }
                }

                //XYZ return_coefficient = ray_data.coefficient * material.fast_BRDF_co(normal, diffuse_slope, flipped_output);
                material.set_input(diffuse_slope);
                XYZ return_coefficient = material.fast_BRDF_co(false);
#if DRAW_DIFFUSE_COEFFICIENTS
                aggregate += return_coefficient * diffuse_multiplier;
                continue;
#endif
                if (XYZ::magnitude(return_coefficient) < 0.0001) continue;

                return_coefficient = return_coefficient * diffuse_multiplier;

                auto ray = PackagedRay(
                    ray_data.position + 0.001 * results.normal,
                    diffuse_slope,
                    ray_data.generation + 1
                );
                ray.output = ray_data.output;
                stats.diffuses_cast++;
                XYZ returned_light = process_ray(stats, ray);
                aggregate += return_coefficient * returned_light;
            }
            for (int i = 0; i < specular_bounces;i++) {
                XYZ specular_slope;
                int tries = 0;
                do {
                    if (tries > 10) throw exception();
                    if (material.roughness > 0) {
                        specular_slope = material.reflective_bounce(gen, reflection_rot_m);
                    }
                    else {
                        specular_slope = reflection_slope;
                    }
                } while (false);
                material.use_fresnel = false;
                material.set_input(specular_slope);
                XYZ return_coefficient = material.fast_BRDF_co(true);
                //XYZ return_coefficient = XYZ(1, 1, 1);
#if DRAW_REFLECTIVE_COEFFICIENTS
                aggregate += return_coefficient * specular_multiplier;
                continue;
#endif
#if DRAW_DOT_RETURNS

#endif
                
                //stringstream s;
                //float f = XYZ::dot(flipped_output, specular_slope);
                //s << left << setw(10) << ((f>0) ? "+" + to_string(f) : to_string(f)) << flipped_output.to_string() << " " << specular_slope.to_string() << " " << normal.to_string() << " " << endl;
                //cout << s.str();

                if (DRAW_BOUNCE_DIRECTION) {
                    if (ray_data.generation == 0) {
                        aggregate += specular_multiplier * (specular_slope / 2 + XYZ(0.5, 0.5, 0.5)) * 1;
                        continue;
                    }

                }

                if (XYZ::magnitude(return_coefficient) < 0.0001) continue;
                return_coefficient = return_coefficient * specular_multiplier;


                auto ray = PackagedRay(
                    ray_data.position + 0.001 * results.normal,
                    specular_slope,
                    ray_data.generation + 1
                );

                ray.output = ray_data.output;
                stats.diffuses_cast++;
                XYZ returned_light = process_ray(stats, ray);
                aggregate += return_coefficient * returned_light;

            }
            */
        }
        else {
            material.set_input(reflection_slope);
            XYZ return_coefficient = material.fast_BRDF_co(true);

            if (XYZ::equals(return_coefficient, XYZ(0, 0, 0))) {
                return aggregate;
            }
            auto ray = PackagedRay(
                ray_data.position + NEAR_THRESHOLD * results.normal * 10,
                reflection_slope,
                ray_data.generation + 1
            );
            ray.output = ray_data.output;
            aggregate += return_coefficient * process_ray(stats, ray);

            stats.reflections_cast++;
        }
        return aggregate;
}
};
