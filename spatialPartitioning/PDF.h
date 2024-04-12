#pragma once

#include <iostream>
#include "Config.h"
#include "commons.h"
#include <chrono>

using namespace std::chrono;

namespace PDF {
    bool prepped = false;
    float* lookup;
    int lookup_resolution = PDF_LOOKUP_RESOLUTION;
    int lookup_size = lookup_resolution * lookup_resolution;

    const double lookup_start = 0;
    const double lookup_end = 1 - lookup_start;
    const double lookup_range = lookup_end - lookup_start;

    double t_offset = -PI / 2;
    double segment_scalar = lookup_range / (double)lookup_resolution;
    double integral_scalar = 1.0 / (double)PDF_INTEGRAL_RESOLUTION;

    double max_integration_step = integral_scalar*2;
    double min_integration_step = integral_scalar / 10;

    double a_start = 0.001;
    double a_range = 1 - a_start;

    double normal_dist_GGXTR(const double theta, const double a) {
        double dot = cos(theta);
        double final = (a * a)
            /
            (PI * pow((dot * dot) * (a * a - 1) + 1, 2));

        return final;
    }
    float normal_dist_GGX(const double theta, const double a) {
        return pow(1 + pow(tan(theta), 2) / a*a, -2);
    }
    double NDF(double theta, double a) {
        return normal_dist_GGXTR(theta, a);
    }
    double get_pdf_weight(double theta, double a) {
        return cos(theta) * NDF(theta, a);
    }
    double find_weight(double a, double val, double start, double iteration_constant, int iterations) {
        const double epsilon = 0.00000000001;

        double check_theta = start;
        double point_value = 0;

        double high = check_theta;
        double low = start;
        while (point_value < val) {
            low = check_theta;
            check_theta += iteration_constant;
            high = check_theta;
            point_value = (check_theta - start) / PI * get_pdf_weight(check_theta, a);
            if (check_theta > PI / 2) {
                high = PI/2;
                break;
            }
        }

        double low_v = (low - start) / PI * get_pdf_weight(low, a);;
        double high_v = (high - start) / PI * get_pdf_weight(high, a);
        for (int i = 0; i < iterations; i++) {
            check_theta = (low + high) / 2;
            point_value = (check_theta - start) / PI * get_pdf_weight(check_theta, a);
            if (abs(point_value - val) < epsilon) break;
            if ((point_value < val)&&(low_v<high_v)) {
                low = check_theta;
                low_v = point_value;
            }
            else {
                high = check_theta;
                high_v = point_value;
            }
        }
        return check_theta;
    }
    void prep() {
        steady_clock::time_point start = steady_clock::now();
        lookup = (float*)_aligned_malloc(sizeof(float) * lookup_size, 64);
        for (int i_a = 0; i_a < lookup_resolution; i_a++) { //loops once per a value
            steady_clock::time_point astart = steady_clock::now();
            double a = i_a * segment_scalar*a_range+a_start;
            int lookup_offset = i_a * lookup_resolution;
            double int_total = 0;
            double gen_total = 0;
            double prev_scaled_total = 0;
            double prev_theta = -PI / 2;
            double theta_int = -PI / 2;
            while(theta_int<PI/2) {
                double int_step = pow(sin(theta_int), 2) * max_integration_step + min_integration_step;
                theta_int += PI * int_step;
                double weight = get_pdf_weight(theta_int, a);
                int_total += weight * int_step;
            }
            double v_step = int_total / (double)lookup_resolution;
            double v_pos = 0;
            double t_pos = -PI/2;
            double iteration = std::min(pow(a,1), 0.1);
            for (int i = 0; i < lookup_resolution; i++) {
                if (abs(t_pos) < 0.01) {
                    cout << "";
                }
                t_pos = find_weight(a, v_step, t_pos, sin(t_pos)*sin(t_pos)*iteration+iteration/100, 5);
                v_pos += v_step;
                lookup[lookup_offset + i] = t_pos;
            }
            //for (int i = 0; i < lookup_resolution; i++) {
            //    std::cout << "(" << lookup[lookup_offset+i] << "," << v_step * i/int_total << "), ";
            //}
            steady_clock::time_point aend = steady_clock::now();
            //cout << a << ": " << duration_cast<milliseconds>(end - start).count() << endl;
        }
        prepped = true;
        steady_clock::time_point end = steady_clock::now();
        cout << "NDF Inverse Lookup Table Generation Done [" << duration_cast<milliseconds>(end - start).count() << "ms]" << endl;

    }
    float sample(float factor, float a) {
        int offset = round(a * (lookup_resolution-1))*lookup_resolution;
        int index = round(factor * (lookup_resolution-1));
        index = std::max(0,std::min(1023, index));
        return lookup[offset + index];
    }

}
