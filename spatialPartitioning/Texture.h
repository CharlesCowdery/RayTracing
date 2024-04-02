#pragma once

#include "XYZ.h"
#include <malloc.h>
#include <immintrin.h>
#include <math.h>

#define USE_MORTON 0

using namespace std;

template <typename T>
class Texture {
public:
    int resolution_x;
    int resolution_y;
    double U_multiplier;
    double V_multiplier;
    XY scale = XY(1, 1);
    XY offset = XY(0, 0);
    T* data;//was originally a vector, but the red tape surrounding vector made texture loading very slow in debug mode


    Texture() {
        data = nullptr;
    }
    Texture(T* data_ptr) {
        data = data_ptr;
    }
    Texture(unsigned char* img, int width, int height, int components, T(*converter)(float, float, float)) {
        allocate(width, height);
        unsigned int position = 0;
        unsigned int data_pos = 0;
        for (unsigned int y = 0; y < resolution_y; y++) {
            for (unsigned int x = 0; x < resolution_x; x++) {
                position += components;
                float r = ((float)img[position + 0]) / 255.0;
                float g = ((float)img[position + 1]) / 255.0;
                float b = ((float)img[position + 2]) / 255.0;
                unsigned int address = to_memory(x, y);
                data[address] = (converter(r, g, b));
                data_pos++;
            }
        }
    }
    ~Texture() {
        if (data != nullptr) {
            free(data);
        }
    }
    Texture<T>* clone_shallow() {
        Texture<T>* tex = new Texture<T>();
        tex->resolution_x = resolution_x;
        tex->resolution_y = resolution_y;
        tex->set_transform(scale, offset);
        tex->data = data;
        return tex;
    }


    unsigned int direct_size(unsigned int rX, unsigned int rY) {
        return rX * rY;
    }
    unsigned int morton_size(unsigned int rX, unsigned rY) {
        int largest_axis = max(rX, rY);
        int power_2_axes = ceil(log2(largest_axis));
        return pow(pow(2, power_2_axes), 2);
    }
    unsigned int get_memory_size(unsigned int rX, unsigned int rY) {
#if USE_MORTON 
        return morton_size(rX, rY);
#else          
        return direct_size(rX, rY);
#endif
    }
    void allocate(int rx, int ry) {
        resolution_x = rx;
        resolution_y = ry;
        unsigned int memory_size = get_memory_size(rx, ry);
        data = (T*)_aligned_malloc(sizeof(T) * memory_size, 64);
    }

    void set_transform(XY _scale, XY _offset) {
        scale = _scale;
        offset = _offset;
    }
    int size() {
        return resolution_x * resolution_y;
    }

    Texture<float>* export_channel(int channel) {
        Texture<float>* out = new Texture<float>();
        out->allocate(resolution_x, resolution_y);
        out->set_transform(scale, offset);
        int pos = 0;
        for (int y = 0; y < resolution_y; y++) {
            for (int x = 0; x < resolution_x; x++) {
                out->data[pos] = lookup(x, y)[channel];
                pos++;
            }
        }
        return out;
    }
    void prep() {
        U_multiplier = resolution_x * scale.X;
        V_multiplier = resolution_y * scale.Y;
    }
    unsigned int to_morton(unsigned int x, unsigned int y) {
        return _pdep_u64(x, 0x5555555555555555) | _pdep_u64(y, 0xaaaaaaaaaaaaaaaa);
    }
    unsigned int to_direct_address(unsigned int x, unsigned int y) {
        return y * resolution_x + x;
    }
    unsigned int to_memory(unsigned int x, unsigned int y) {
#if USE_MORTON
        return to_morton(x, y);
#else
        return to_direct_address(x, y);
#endif
    }

    T lookup(int x, int y) {
        x = x % resolution_x;
        y = y % resolution_y;
        if (x < 0) x = resolution_x + x;
        if (y < 0) y = resolution_y + y;
        unsigned int address = to_memory(x, y);
        return data[address];
    }
    T getUV(double U, double V) {
        V = 1 - V;
        int pos_x = (int)floor(U_multiplier * U) % resolution_x;
        int pos_y = (int)floor(V_multiplier * V) % resolution_y;
        return lookup(pos_x, pos_y);
    }
    T getPixelLinear(double U, double V) {
        double x_1 = U_multiplier * U - 0.5 * U_multiplier;
        double x_2 = U_multiplier * U + 0.5 * U_multiplier;
        double y_1 = V_multiplier * V - 0.5 * V_multiplier;
        double y_2 = V_multiplier * V + 0.5 * V_multiplier;
        int pos_x_1 = max(0, std::min(resolution_x - 1, (int)floor(x_1)));
        int pos_x_2 = max(0, std::min(resolution_x - 1, (int)floor(x_2)));
        int pos_y_1 = max(0, std::min(resolution_y - 1, (int)floor(y_1)));
        int pos_y_2 = max(0, std::min(resolution_y - 1, (int)floor(y_2)));
        T c1 = lookup(pos_x_1, pos_y_1);
        T c2 = lookup(pos_x_2, pos_y_1);
        T c3 = lookup(pos_x_1, pos_y_2);
        T c4 = lookup(pos_x_2, pos_y_2);
        double f1 = pos_x_2 - x_1;
        double f2 = x_2 - pos_x_2;
        double f3 = pos_y_2 - y_1;
        double f4 = y_2 - pos_y_2;
        T c_m1 = f1 * c1 + f2 * c2;
        T c_m2 = f1 * c3 + f2 * c4;
        return f3 * c_m1 + f4 * c_m2;
    }
};

template <typename T>
class Parameter {
public:
    T static_value = 0;
    Texture<T>* texture = nullptr;
    Parameter(T value) : static_value(value) {}
    Parameter(Texture<T>* texture_ptr) : texture(texture_ptr) {}
    //Parameter(string file_name) : texture(new XYZTexture(file_name)) {}
    void set_static(T value) {
        static_value = value;
    }
    void prep() {
        if (texture != nullptr) {
            texture->prep();
        }
    }
    void set_texture(Texture<T>* texture_ptr, bool free = false) {
        if (texture != nullptr) {
            if (free) {
                delete texture;
            }
        }
        texture = texture_ptr;
    }
    void set_texture(string file_name, bool free = false) {
        //set_texture(new XYZTexture(file_name), free);
    }
    T getValue(double U = 0, double V = 0) {
        if (texture != nullptr) {
            return texture->getUV(U, V);
        }
        else {
            return static_value;
        }
    }
};