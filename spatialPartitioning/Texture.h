#pragma once

#include "XYZ.h"
#include <malloc.h>
#include <immintrin.h>
#include <math.h>
#include "commons.h"

#define USE_MORTON 0
#define USE_CHAR 1

typedef XYZ(*ConverterFunction)(float, float, float);

class TransferFunction {
public:
    float scale;
    ConverterFunction converter;
    TransferFunction(ConverterFunction f, float s = -9999) {
        converter = f;
        scale = s;
    }
    bool operator==(const TransferFunction& other) const {
        return scale == other.scale && converter == other.converter;
    }
    bool operator<(const TransferFunction& other) const {
        bool c1 = converter < other.converter;
        bool e1 = converter == other.converter;
        bool c2 = scale < other.scale;
        return c1 || (e1 && c2);
    }
};

template <typename T>
class Texture {
private:
    void load_no_scale(unsigned char* img, int width, int height, int components, ConverterFunction converter) {
        unsigned int position = 0;
        unsigned int data_pos = 0;
        for (unsigned int y = 0; y < resolution_y; y++) {
            for (unsigned int x = 0; x < resolution_x; x++) {
                position += components;
                float r = ((float)img[position + 0]) / 255.0;
                float g = ((float)img[position + 1]) / 255.0;
                float b = ((float)img[position + 2]) / 255.0;
                unsigned int address = to_memory(x, y);
                commit(address,converter(r, g, b));
                data_pos++;
            }
        }
    }
    void load_with_scale(unsigned char* img, int width, int height, int components, ConverterFunction converter, float scale) {
        use_normal_transform = true;
        unsigned int position = 0;
        unsigned int data_pos = 0;
        for (unsigned int y = 0; y < resolution_y; y++) {
            for (unsigned int x = 0; x < resolution_x; x++) {
                position += components;
                float r = ((float)img[position + 0]) / 255.0;
                float g = ((float)img[position + 1]) / 255.0;
                float b = ((float)img[position + 2]) / 255.0;
                unsigned int address = to_memory(x, y);
                XYZ in = (converter(r, g, b));
                commit(address,XYZ::normalize((in * 2 - 1) * XYZ(scale, scale, 1.0)));
                data_pos++;
            }
        }
    }
public:
    int resolution_x;
    int resolution_y;
    double U_multiplier;
    double V_multiplier;
    XY scale = XY(1, 1);
    XY offset = XY(0, 0);
    bool use_normal_transform = false;
#if USE_CHAR
    unsigned char* data = nullptr;//was originally a vector, but the red tape surrounding vectors made texture loading very slow in debug mode
#else
    T* data = nullptr;
#endif


    Texture() {
        data = nullptr;
    }
    Texture(unsigned char* img, int width, int height, int components, TransferFunction TF) {
        load(img, width, height, components, TF);
    }
    ~Texture() {
        if (data != nullptr) {
            free(data);
        }
    }
    void load(unsigned char* img, int width, int height, int components, TransferFunction TF) {
        allocate(width, height);
        if (TF.scale == -9999) {
            load_no_scale(img, width, height, components, TF.converter);
        }
        else {
            load_with_scale(img, width, height, components, TF.converter, TF.scale);
        }
    }

    Texture<T>* clone_shallow() {
        Texture<T>* tex = new Texture<T>();
        tex->resolution_x = resolution_x;
        tex->resolution_y = resolution_y;
        tex->set_transform(scale, offset);
        tex->data = data;
        tex->use_normal_transform = use_normal_transform;
        return tex;
    }
    void commit_default(unsigned int addr, const T& v) {
        data[addr] = v;
    }
    void commit_char(unsigned int addr, const XYZ& v) {
        if (use_normal_transform) {
            unsigned char r = v.X * 127+128;
            unsigned char g = v.Y * 127+128;
            unsigned char b = v.Z * 127+128;
            data[addr * 3 + 0] = r;
            data[addr * 3 + 1] = g;
            data[addr * 3 + 2] = b;
            return;
        }
        unsigned char r = v.X * 255;
        unsigned char g = v.Y * 255;
        unsigned char b = v.Z * 255;
        data[addr * 3 + 0] = r;
        data[addr * 3 + 1] = g;
        data[addr * 3 + 2] = b;
    }
    void commit_char(unsigned int addr, const float& v) {
        unsigned char r = v * 255;
        data[addr * 3 + 0] = r;
    }
    void commit(unsigned int addr, const T& v){
#if USE_CHAR
        commit_char(addr, v);
#else
        commit_default(addr, v);
#endif
    }
    T fetch_default(unsigned int addr) {
        return data[addr];
    }
    XYZ fetch_char_XYZ(unsigned int addr) {
        unsigned char r = data[addr * 3 + 0];
        unsigned char g = data[addr * 3 + 1];
        unsigned char b = data[addr * 3 + 2];
        if (use_normal_transform) {
            return XYZ((r - 128) / 127.0f, (g - 128) / 127.0f, (b - 128) / 127.0f);
        }
        return XYZ(r / 255.0f, g / 255.0f, b / 255.0f);
    }
    float fetch_char_float(unsigned int addr) {
        unsigned char r = data[addr * 3 + 0];
        return r/255.0f;
    }
    template<typename Tf> Tf fetch(unsigned int addr) {
#if USE_CHAR
        return Texture<float>::fetch_char_float(addr);
#else
        return fetch_default(addr);
#endif
    }
    template<> XYZ fetch<XYZ>(unsigned int addr) {
#if USE_CHAR
        return Texture<XYZ>::fetch_char_XYZ(addr);
#else
        return fetch_default(addr);
#endif
    }
    unsigned int direct_size(unsigned int rX, unsigned int rY) {
        return rX * rY;
    }
    unsigned int morton_size(unsigned int rX, unsigned rY) {
        int largest_axis = std::max(rX, rY);
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
#if USE_CHAR
    void allocate(int rx, int ry) {
        resolution_x = rx;
        resolution_y = ry;
        unsigned int memory_size = get_memory_size(rx, ry);
        data = (unsigned char*)_aligned_malloc(sizeof(char) * 3 * memory_size, 64);
    }
#else
    void allocate(int rx, int ry) {
        resolution_x = rx;
        resolution_y = ry;
        unsigned int memory_size = get_memory_size(rx, ry);
        data = (T*)_aligned_malloc(sizeof(T) * memory_size, 64);
    }
#endif
    

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
                out->commit(pos,lookup(x, y)[channel]);
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

    T lookup(int lx, int ly) {
        int x = lx % resolution_x;
        int y = ly % resolution_y;
        x = (x < 0)*resolution_x + x;
        y = (y < 0)*resolution_y + y;
        unsigned int address = to_memory(x, y);
        return fetch<T>(address);
    }
    T getUV(double U, double V) {
        V = 1 - V;
        int pos_x = (int) (U_multiplier * U);
        int pos_y = (int) (V_multiplier * V);
        return lookup(pos_x, pos_y);
    }
    T getPixelLinear(double U, double V) {
        double x_1 = U_multiplier * U - 0.5 * U_multiplier;
        double x_2 = U_multiplier * U + 0.5 * U_multiplier;
        double y_1 = V_multiplier * V - 0.5 * V_multiplier;
        double y_2 = V_multiplier * V + 0.5 * V_multiplier;
        int pos_x_1 = std::max(0, std::min(resolution_x - 1, (int)floor(x_1)));
        int pos_x_2 = std::max(0, std::min(resolution_x - 1, (int)floor(x_2)));
        int pos_y_1 = std::max(0, std::min(resolution_y - 1, (int)floor(y_1)));
        int pos_y_2 = std::max(0, std::min(resolution_y - 1, (int)floor(y_2)));
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
        if (texture_ptr == NULL) throw std::exception("Illegal texture address");
        texture = texture_ptr;
    }
    void set_texture(std::string file_name, bool free = false) {
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