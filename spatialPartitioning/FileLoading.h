#pragma once

#include <chrono>
#include <thread>
#include <vector>
#include <fstream>
#include <set>

#include "commons.h"
#include "XYZ.h"
#include "Primitives.h"
#include "Core.h"
#include "SceneManager.h"


#define TINYGLTF_IMPLEMENTATION
#include "tiny_gltf.h"


using std::ofstream;
using std::ios;
using std::ifstream;
using std::iostream;
using std::vector;
using std::is_same;
using std::set;

using namespace std::chrono;


class FileManager {
public:
    FileManager() {}
    class TextureLoader {
    public:
        class TextureWorker {
        public:
            tinygltf::Model* data;
            thread worker_thread;
            atomic<bool> running;
            atomic<bool> stop;
            atomic<int> job_address;
            atomic<int> converter_address;
            map<int, map<TransferFunction, Texture<XYZ>*>>* texture_lookup;
            TextureWorker(tinygltf::Model* _data, map<int, map<TransferFunction, Texture<XYZ>*>>* _texture_lookup) :
                data(_data), texture_lookup(_texture_lookup) {
                running = false;
                stop = false;
                job_address = -1;
                worker_thread = thread(&TextureWorker::main, this);
            }
            void main() {
                while (true) {
                    if (running.load(std::memory_order_relaxed)) {
                        runJob();
                    }
                    if (stop.load(std::memory_order_relaxed)) {
                        return;
                    }
                    std::this_thread::sleep_for(milliseconds(10));
                }
            }
            void runJob() {
                int address = job_address.load(std::memory_order_relaxed);
                int c_address = converter_address.load(std::memory_order_relaxed);
                auto map_iterator = (*texture_lookup)[address].begin();
                for (int i = 0; i < c_address; i++) map_iterator++;
                TransferFunction TF = map_iterator->first;
                Texture<XYZ>* texture_ptr = map_iterator->second;
                loadTexture<XYZ>(*data, address, TF, false, texture_ptr);
                running.store(false, std::memory_order_relaxed);
            }
        };
        tinygltf::Model* data;
        map<int, map<TransferFunction, Texture<XYZ>*>> texture_lookup;
        vector<TextureWorker*> workers;
        void createWorkers(int num) {
            for (int i = 0; i < num; i++) {
                TextureWorker* worker = new TextureWorker(data,&texture_lookup);
                string name = "TextureLoadingThread-" + to_string(i);
                wstring wname = wstring(name.begin(), name.end());
                const wchar_t* wcname = wname.c_str();
                workers.push_back(worker);
                SetThreadPriority(worker->worker_thread.native_handle(), 0);
                SetThreadDescription(worker->worker_thread.native_handle(), wcname);
            }
        }
        void dispatchThreads() {
            int threads_stopped = 0;
            auto tex_iterator = texture_lookup.begin();
            auto tex_end = texture_lookup.end();
            if (tex_iterator == tex_end) return;
            auto TF_iterator = 0;
            auto TF_end = tex_iterator->second.size();
            while (threads_stopped < workers.size()) {
                for (int i = 0; i < workers.size();i++) {
                    TextureWorker* worker = workers[i];
                    if (worker->stop.load(std::memory_order_relaxed)) continue;
                    if (worker->running.load(std::memory_order_relaxed)) continue;
                    if (TF_iterator == TF_end) {
                        if (tex_iterator!=tex_end)tex_iterator++;
                        if (tex_iterator == tex_end) {
                            worker->stop.store(true, std::memory_order_relaxed);
                            threads_stopped++;
                            continue;
                        }

                        TF_iterator = 0;
                        TF_end = tex_iterator->second.size();
                    }
                    worker->job_address = tex_iterator->first;
                    worker->converter_address = TF_iterator;
                    worker->running.store(true, std::memory_order_relaxed);
                    TF_iterator++;
                }
            }
        }
        void loadTextures(tinygltf::Model* _data, int thread_count) {
            data = _data;
            cout << padString("", " ", 100) << "\r";
            cout << "[Load] constructing texture lookup" << endl;
            for (auto& material_data : data->materials) {
                auto pbr = material_data.pbrMetallicRoughness;
                if (pbr.baseColorTexture.index != -1) {
                    int t_index = pbr.baseColorTexture.index;
                    TransferFunction tf = TransferFunction(converter_xyz_sRGB);
                    if (!texture_lookup.contains(t_index)) {
                        texture_lookup[t_index] = map<TransferFunction, Texture<XYZ>*>();
                    }
                    if (!texture_lookup[t_index].contains(tf)) {
                        texture_lookup[t_index][tf] = new Texture<XYZ>();
                    }
                }
                if (pbr.metallicRoughnessTexture.index != -1) {
                    int t_index = pbr.metallicRoughnessTexture.index;
                    TransferFunction tf = TransferFunction(converter_xyz);
                    if (!texture_lookup.contains(t_index)) {
                        texture_lookup[t_index] = map<TransferFunction, Texture<XYZ>*>();
                    }
                    if (!texture_lookup[t_index].contains(tf)) {
                        texture_lookup[t_index][tf] = new Texture<XYZ>();
                    }
                }
                if (material_data.normalTexture.index != -1) {
                    int t_index = material_data.normalTexture.index;
                    const float scale = material_data.normalTexture.scale;
                    TransferFunction tf = TransferFunction(converter_xyz, scale);
                    if (!texture_lookup.contains(t_index)) {
                        texture_lookup[t_index] = map<TransferFunction, Texture<XYZ>*>();
                    }
                    if (!texture_lookup[t_index].contains(tf)) {
                        texture_lookup[t_index][tf] = new Texture<XYZ>();
                    }
                }
            }
            createWorkers(thread_count);
            dispatchThreads();
            for (int i = 0; i < workers.size(); i++) {
                workers[i]->stop = true;
                workers[i]->worker_thread.join();
                delete workers[i];
            }
        }

    };
    static void writeRawFile(vector<vector<XYZ*>>* data, string fName) {
        ofstream file(fName, ios::out | ios::binary | ios::trunc);
        FileManager::writeRaw(data, file);
        file.close();
    }
    static void writeRaw(vector<vector<XYZ*>>* data, ofstream& file) {
        unsigned int resX = data->at(0).size();
        unsigned int resY = data->size();
        cout << resX << " " << resY << endl;
        file.write((char*)&resX, sizeof(unsigned int));
        file.write((char*)&resY, sizeof(unsigned int));
        for (vector<XYZ*>& data_row : *data) {
            for (XYZ* pixel_ptr : data_row) {
                file.write((char*)pixel_ptr, sizeof(XYZ));
            }
        }
    }
    static vector<vector<XYZ*>>* readRaw(ifstream& file) {
        vector<vector<XYZ*>>* data = new vector<vector<XYZ*>>();
        int resolutionX = 0;
        int resolutionY = 0;
        file.read((char*)&resolutionX, sizeof(int));
        file.read((char*)&resolutionY, sizeof(int));
        cout << resolutionX << " " << resolutionY << endl;
        for (int y = 0; y < resolutionY; y++) {
            vector<XYZ*> output_row;
            for (int x = 0; x < resolutionX; x++) {
                output_row.push_back(new XYZ(0, 0, 0));
                file.read((char*)output_row[x], sizeof(XYZ));
            }
            data->push_back(output_row);
        }
        return data;
    }
    static Mesh* loadObjFile(string fName, Material* mat) {
        cout << padString("[File] Loading resource: " + fName, ".", 100);
        ifstream file(fName, ios::in);
        Mesh* return_value = loadObj(file, mat);
        cout << "[Done]" << endl;
        return return_value;
    }
    static Mesh* loadObj(ifstream& file, Material* mat) {
        Mesh* mesh = new Mesh();
        vector<XYZ> verts;
        vector<XY> UV_coords;
        for (std::string line; getline(file, line);) {
            string line_type = line.substr(0, line.find(" "));
            vector<string> words = split_string(line, ' ');
            if (line_type == "#") {
                continue;
            }
            if (line_type == "mtllib") {
                continue;
            }
            if (line_type == "o") {
                continue;
            }
            if (line_type == "v") {
                verts.push_back(XYZ(
                    stod(words[1]),
                    stod(words[2]),
                    stod(words[3])
                ));
                //cout << verts.back() << endl;
            }
            if (line_type == "vt") {
                UV_coords.push_back(XY(
                    stod(words[1]), stod(words[2])
                ));
            }
            if (line_type == "f") {
                vector<string> sub_words_1 = split_string(words[1], '/');
                vector<string> sub_words_2 = split_string(words[2], '/');
                vector<string> sub_words_3 = split_string(words[3], '/');
                int v_i0 = stoi(sub_words_1[0]) - 1;
                int v_i1 = stoi(sub_words_2[0]) - 1;
                int v_i2 = stoi(sub_words_3[0]) - 1;
                XY  vt1 = XY(-1, -1);
                XY  vt2 = XY(-1, -1);
                XY  vt3 = XY(-1, -1);
                if (sub_words_1.size() > 1) {
                    if (sub_words_1[1] != "") {
                        int vt_i0 = stoi(sub_words_1[1]) - 1;
                        int vt_i1 = stoi(sub_words_2[1]) - 1;
                        int vt_i2 = stoi(sub_words_3[1]) - 1;
                        vt1 = UV_coords[vt_i0];
                        vt2 = UV_coords[vt_i1];
                        vt3 = UV_coords[vt_i2];
                    }
                }
                Tri f = Tri(
                    verts[v_i0], verts[v_i1], verts[v_i2],
                    vt1, vt2, vt3,
                    mat);
                mesh->addTri(f);
            }
        }
        return mesh;
    }
    static Mesh* loadTriFile(string fName, Material* mat) {
        ifstream file(fName, ios::in);
        return loadTri(file, mat);
    }
    static Mesh* loadTri(ifstream& file, Material* mat) {
        Mesh* mesh = new Mesh();
        vector<XYZ> verts;
        for (std::string line; getline(file, line);) {
            vector<string> words = split_string(line, ' ');
            if (words.size() > 1) {
                double v1_x = stod(words[0]);
                double v1_y = stod(words[1]);
                double v1_z = stod(words[2]);
                double v2_x = stod(words[3]);
                double v2_y = stod(words[4]);
                double v2_z = stod(words[5]);
                double v3_x = stod(words[6]);
                double v3_y = stod(words[7]);
                double v3_z = stod(words[8]);
                Tri f1 = Tri(
                    XYZ(v1_x, v1_y, v1_z),
                    XYZ(v2_x, v2_y, v2_z),
                    XYZ(v3_x, v3_y, v3_z),
                    mat);
                mesh->addTri(f1);
                //cout << verts.back() << endl;
            }
        }
        return mesh;
    }
    static SceneManager* loadGLTFFile(string fName) {
        tinygltf::Model model;
        tinygltf::TinyGLTF loader;
        std::string err;
        std::string warn;

        cout << padString("[File] Loading scene " + fName + "...\n", ".", 0);
        //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, argv[1]);
        bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, fName); // for binary glTF(.glb)

        if (!warn.empty()) {
            printf("Warn: %s\n", warn.c_str());
        } if (!err.empty()) {
            printf("Err: %s\n", err.c_str());
        } if (!ret) {
            throw exception("Failed to parse glTF\n");
        }

        SceneManager* out = parseGLTFData(model);

        cout << "[Load] loading done";

        return out;

    }
    template<typename T> struct iterator_pair {
        std::_Vector_iterator<std::_Vector_val<std::_Simple_types<T>>> begin;
        std::_Vector_iterator<std::_Vector_val<std::_Simple_types<T>>> end;
        T* begin_ptr;
        iterator_pair() {}
    };
    struct buffer_accessor {
        tinygltf::Model& data;
        vector<char> swizzle;
        buffer_accessor(tinygltf::Model& _data) : data(_data) {}
        iterator_pair<unsigned char> get_buffer(int accessor_index) {
            iterator_pair<unsigned char> returner;
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            int buffer_index = view.buffer;
            int start_index = view.byteOffset;
            int end_index = start_index + view.byteLength;

            auto& buffer = data.buffers[buffer_index];
            returner.begin = buffer.data.begin() + start_index;
            returner.end = buffer.data.begin() + end_index;
            returner.begin_ptr = &(*returner.begin);

            return returner;

        }
        void get_float_data(int accessor_index, vector<float>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5126);
            assert(view.byteLength % 4 == 0);

            target.clear();
            target.resize(view.byteLength / 4, 0);
            std::copy_n(reinterpret_cast<float*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_ushort_data(int accessor_index, vector<unsigned short>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5123);
            assert(view.byteLength % 2 == 0);

            target.clear();
            target.resize(view.byteLength / 2, 0);
            std::copy_n(reinterpret_cast<unsigned short*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_uint_data(int accessor_index, vector<unsigned int>& target) {
            iterator_pair<unsigned char> iterators = get_buffer(accessor_index);
            auto& accessor = data.accessors[accessor_index];
            auto& view = data.bufferViews[accessor.bufferView];

            assert(accessor.componentType == 5125);
            assert(view.byteLength % 4 == 0);

            target.clear();
            target.resize(view.byteLength / 4, 0);
            std::copy_n(reinterpret_cast<unsigned int*>(iterators.begin_ptr), target.size(), target.begin());
        }
        void get_data_any_int(int accessor_index, vector<unsigned int>& target) {
            auto& accessor = data.accessors[accessor_index];
            vector<unsigned short> intermediate;
            switch (accessor.componentType) {
            case 5125:
                get_uint_data(accessor_index, target);
                break;
            case 5123:
                get_ushort_data(accessor_index, intermediate);
                target.clear();
                target.reserve(intermediate.size());
                for (int i = 0; i < intermediate.size(); i++) {
                    target.push_back((int)intermediate[i]);
                }
                break;
            default:
                assert(0);
                break;
            }
        }
        void get_vec4_data(int accessor_index, vector<Quat>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 4);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(Quat(float_data[i * 4], float_data[i * 4 + 1], float_data[i * 4 + 2], float_data[i * 4 + 3]));
            }
        }
        void get_vec4_data_strip(int accessor_index, vector<XYZ>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 4);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XYZ(float_data[i * 4], float_data[i * 4 + 1], float_data[i * 4 + 2]));
            }
        }
        void get_vec3_data(int accessor_index, vector<XYZ>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 3);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XYZ(float_data[i * 3], float_data[i * 3 + 1], float_data[i * 3 + 2]));
            }
        }
        void get_vec2_data(int accessor_index, vector<XY>& target) {
            auto& accessor = data.accessors[accessor_index];
            assert(accessor.type == 2);

            vector<float> float_data;
            get_float_data(accessor_index, float_data);
            for (int i = 0; i < accessor.count; i++) {
                target.push_back(XY(float_data[i * 2], float_data[i * 2 + 1]));
            }
        }
    };
    struct TextureWrapper {
        Texture<XYZ>* xyz;
        Texture<float>* f;
    };
    template <int t> static float converter_float(float x, float y, float z) {
        if (t == 0) return x;
        if (t == 1) return y;
        if (t == 2) return z;
    }
    static XYZ converter_xyz(float x, float y, float z) {
        return XYZ(x, y, z);
    }
    static XYZ converter_xyz_sRGB(float x, float y, float z) {
        return ImageHandler::apply_gamma(XYZ(x, y, z), 1 / 2.2);
    }
    static XYZ converter_xzy(float x, float y, float z) {
        return XYZ(x, z, y);
    }
    template<float scale>
    class converter_normal_scaler {
    public:
        static XYZ func(float x, float y, float z) {
            return (XYZ::normalize(XYZ(x, y, z)) * 2 - 1) * XYZ(scale, scale, 1 / 2.2);
        }

    };
    template <typename T> static Texture<T>* loadTexture(tinygltf::Model& data, int index, TransferFunction TF, bool verbose = false) {
        Texture<T>* tex_ptr = new Texture<T>();
        loadTexture<T>(data, index, TF, verbose, tex_ptr);
        return tex_ptr;
    }
    template <typename T> static void loadTexture(tinygltf::Model& data, int index, TransferFunction TF, bool verbose, Texture<T>* tex_ptr) {
        auto& texture_data = data.textures[index];
        int image_index = texture_data.source;
        auto& image_data = data.images[image_index];
        int len = image_data.image.size();
        int res_x = image_data.width;
        int res_y = image_data.height;
        int component = image_data.component;

        if(verbose) cout << endl;
        if(verbose) cout << "[Load] Loading texture \"" << image_data.name << "\"\r" << flush;

        if (component != 3 && component != 4) throw exception("illegal texture!");

        int is_addressable = -1;
        if (is_same<T, XYZ>()) {
            is_addressable = 1;
        }
        if (is_same<T, float>()) {
            is_addressable = 0;
        }
        if (is_addressable == -1) throw exception("illegal texture typing!");

        tex_ptr->load(image_data.image.data(), res_x, res_y, component, TF);

        image_data.image.clear(); //freeing memory

        if(verbose) cout << padString("", " ", 100) << "\x1b[A\r";
    }

    static pair<XY, XY> fetch_transform(tinygltf::TextureInfo& Ti) {
        XY scale = XY(1, 1);
        XY offset = XY(0, 0);
        if (Ti.extensions.count("KHR_texture_transform")) {
            auto& ext = Ti.extensions["KHR_texture_transform"];
            auto& offset_object = ext.Get("offset");
            auto& scale_object = ext.Get("scale");
            offset.X = offset_object.Get(0).GetNumberAsDouble();
            offset.Y = offset_object.Get(1).GetNumberAsDouble();
            scale.X = scale_object.Get(0).GetNumberAsDouble();
            scale.Y = scale_object.Get(1).GetNumberAsDouble();
        }
        return { scale,offset };
    }
    static pair<XY, XY> fetch_transform(tinygltf::NormalTextureInfo& Ti) {
        XY scale = XY(1, 1);
        XY offset = XY(0, 0);
        if (Ti.extensions.count("KHR_texture_transform")) {
            auto& ext = Ti.extensions["KHR_texture_transform"];
            auto& offset_object = ext.Get("offset");
            auto& scale_object = ext.Get("scale");
            offset.X = offset_object.Get(0).GetNumberAsDouble();
            offset.Y = offset_object.Get(1).GetNumberAsDouble();
            scale.X = scale_object.Get(0).GetNumberAsDouble();
            scale.Y = scale_object.Get(1).GetNumberAsDouble();
        }
        return { scale,offset };
    }
    static SceneManager* parseGLTFData(tinygltf::Model data) {

        const bool print_in_place = true;

        Scene* scene = new Scene();
        SceneManager* SM = new SceneManager(scene);
        vector<Material*> materials;
        vector<Mesh*> meshes;
        vector<Camera*> cameras;
        vector<Object*> objects;
        vector<Light*> lights;
        map<int, pair<int, int>> node_lookup;
        vector<char> swizzle = { 0,-2,1 };
        vector<char> rot_swizzle = { 0,-2, 1, 3 };
        vector<char> scale_swizzle = { 0,2,1 };
        buffer_accessor buf_accessor = buffer_accessor(data);

        TextureLoader TL = TextureLoader();
        TL.loadTextures(&data,12);

        map<int,map<TransferFunction,Texture<XYZ>*>>& texture_lookup = TL.texture_lookup;
        
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Textures loaded" << endl;
        for (auto& material_data : data.materials) {
            string name = material_data.name;
            cout << padString("", " ", 100) << "\r";
            cout << "[Load] Loading material \"" << name << "\"\r" << flush;
            Material* mat = new Material();
            auto pbr = material_data.pbrMetallicRoughness;
            auto emissive = material_data.emissiveFactor;

            mat->color.set_static(XYZ(pbr.baseColorFactor));
            mat->metallic.set_static(pbr.metallicFactor);
            mat->roughness.set_static(pbr.roughnessFactor);
            mat->emissive.set_static(XYZ(emissive) * 25);

            mat->name = name;

            if (pbr.baseColorTexture.index != -1) {
                int t_index = pbr.baseColorTexture.index;
                TransferFunction tf = TransferFunction(converter_xyz_sRGB);
                mat->UV_map_index = pbr.baseColorTexture.texCoord;
                Texture<XYZ>* base_tex = texture_lookup[t_index][tf];
                Texture<XYZ>* tex = base_tex->clone_shallow();
                auto transform_data = fetch_transform(pbr.baseColorTexture);
                tex->set_transform(transform_data.first, transform_data.second);
                mat->color.set_texture(tex);
                
            }
            if (pbr.metallicRoughnessTexture.index != -1) {
                int t_index = pbr.metallicRoughnessTexture.index;
                TransferFunction tf = TransferFunction(converter_xyz);
                mat->UV_map_index = pbr.baseColorTexture.texCoord;
                Texture<XYZ>* base_tex = texture_lookup[t_index][tf];
                Texture<XYZ>* tex = base_tex->clone_shallow();
                auto& texture_data = data.textures[t_index];
                auto& image_data = data.images[texture_data.source];
                cout << endl << "[Load] Binding texture \"" << image_data.name << "\"\r" << flush;

                auto transform_data = fetch_transform(pbr.metallicRoughnessTexture);
                tex->set_transform(transform_data.first, transform_data.second);

                Texture<float>* metallic = tex->export_channel(2);
                Texture<float>* roughness = tex->export_channel(1);


                if (pbr.metallicFactor != 0) {
                    if (pbr.metallicFactor != 1) {
                        for (int i = 0; i < metallic->size(); i++) {
                            metallic->data[i] *= pbr.metallicFactor;
                        }
                    }
                    mat->metallic.set_texture(metallic);
                }
                if (pbr.roughnessFactor != 0) {
                    if (pbr.roughnessFactor != 1) throw exception("Illegal roughness factor");
                    mat->roughness.set_texture(roughness);
                }
                cout << padString("", " ", 100) << "\x1b[A\r";
            }
            if (material_data.normalTexture.index != -1) {
                int t_index = material_data.normalTexture.index;
                mat->UV_map_index = material_data.normalTexture.texCoord;
                const float scale = material_data.normalTexture.scale;
                TransferFunction tf = TransferFunction(converter_xyz,scale);
                Texture<XYZ>* base_tex = texture_lookup[t_index][tf];
                Texture<XYZ>* tex = base_tex->clone_shallow();
                auto& texture_data = data.textures[t_index];
                auto& image_data = data.images[texture_data.source];

                auto transform_data = fetch_transform(material_data.normalTexture);
                tex->set_transform(transform_data.first, transform_data.second);

                mat->normal.set_texture(tex);
                mat->use_normals = 1;
            }
            if (material_data.extensions.count("KHR_materials_specular")) {
                auto specular_iterator = material_data.extensions.find("KHR_materials_specular");
                if (specular_iterator != material_data.extensions.end()) {
                    auto specular_value = specular_iterator->second.Get("specularColorFactor");
                    if (specular_value.Type() == 0) {
                        mat->specular.set_static(0.0f);
                    }
                    else {
                        auto specular_0 = specular_value.Get(0).GetNumberAsDouble();
                        auto specular_1 = specular_value.Get(1).GetNumberAsDouble();
                        auto specular_2 = specular_value.Get(2).GetNumberAsDouble();
                        mat->specular.set_static(specular_0);
                    }
                }
            }
            else {
                mat->specular = 1;
            }
            if (material_data.extensions.count("KHR_materials_ior")) {
                auto ior_iterator = material_data.extensions.find("KHR_materials_ior");
                if (ior_iterator != material_data.extensions.end()) {
                    auto ior_value = ior_iterator->second.Get("ior");
                    auto ior = ior_value.GetNumberAsDouble();
                    mat->IOR.set_static(ior);
                }
            }
            else {
                mat->IOR = 1.5;
            }
            materials.push_back(mat);
        }
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Materials loaded" << endl;
        for (auto& model_data : data.meshes) {
            auto& primitives = model_data.primitives;
            Mesh* mesh = new Mesh();
            mesh->name = model_data.name;
            cout << padString("", " ", 100) << "\r";
            cout << "[Load] Loading mesh \"" << mesh->name << "\"\r" << flush;
            for (auto& primitive : primitives) {
                auto& attributes = primitive.attributes;
                Material* mat;
                if (primitive.material > -1) {
                    mat = materials[primitive.material];
                }
                else {
                    mat = new Material();
                }

                assert(attributes.count("POSITION") > 0);
                auto  position_accessor_index = attributes["POSITION"];
                auto texcoord_accessor_index = attributes["TEXCOORD_" + to_string(mat->UV_map_index)];
                auto normal_accessor_index = attributes["NORMAL"];
                auto tangent_accessor_index = attributes["TANGENT"];
                auto tangent_accessor = data.accessors[tangent_accessor_index];

                vector<XYZ> position_data;
                buf_accessor.get_vec3_data(position_accessor_index, position_data);

                bool UV_mapped = texcoord_accessor_index != 0;
                vector<XY> texcoord_data;
                if (UV_mapped) {
                    buf_accessor.get_vec2_data(texcoord_accessor_index, texcoord_data);
                }

                bool smooth_shading = normal_accessor_index != 0;
                smooth_shading &= 1;
                vector<XYZ> normal_data;
                if (smooth_shading) {
                    buf_accessor.get_vec3_data(normal_accessor_index, normal_data);
                }


                //vector<XYZ> tangent_data;
                //if (tangent_accessor.type == 4) {
                //    buf_accessor.get_vec4_data_strip(tangent_accessor_index, tangent_data);
                //}
                //else {
                //    buf_accessor.get_vec3_data(tangent_accessor_index, tangent_data);
                //}

                assert(primitive.indices >= 0);
                auto  indices_accessor_index = primitive.indices;
                auto& indices_accessor = data.accessors[position_accessor_index];

                vector<unsigned int> indices;
                buf_accessor.get_data_any_int(indices_accessor_index, indices);

                assert(indices.size() % 3 == 0);
                int tri_count = indices.size() / 3;

                mesh->tris.reserve(mesh->tris.size() + tri_count);
                for (int i = 0; i < tri_count; i++) {
                    int look_1 = indices[i * 3 + 0];
                    int look_2 = indices[i * 3 + 1];
                    int look_3 = indices[i * 3 + 2];
                    XYZ v1 = position_data[look_1].swizzle(swizzle);
                    XYZ v2 = position_data[look_2].swizzle(swizzle);
                    XYZ v3 = position_data[look_3].swizzle(swizzle);
                    XY uv1 = 0;
                    XY uv2 = 0;
                    XY uv3 = 0;
                    if (UV_mapped) {
                        uv1 = XY(0, 1) - XY(-1, 1) * texcoord_data[look_1];
                        uv2 = XY(0, 1) - XY(-1, 1) * texcoord_data[look_2];
                        uv3 = XY(0, 1) - XY(-1, 1) * texcoord_data[look_3];
                    }

                    //cout << v1 << " " << v2 << " " << v3 << " " << uv1 << " " << uv2 << " " << uv3 << endl;
                    Tri T = Tri(v1, v2, v3, uv1, uv2, uv3, mat);
                    if (smooth_shading) {
                        T.n1 = normal_data[look_1].swizzle(swizzle);
                        T.n2 = normal_data[look_2].swizzle(swizzle);
                        T.n3 = normal_data[look_3].swizzle(swizzle);
                    }
                    mesh->addTri(T);
                }
            }
            meshes.push_back(mesh);
        }
        cout << padString("", " ", 100) << "\r";
        cout << "[Load] Meshes loaded" << endl;
        for (auto& camera_data : data.cameras) {
            assert(camera_data.type == "perspective");
            auto& perspective = camera_data.perspective;
            auto aspect_ratio = perspective.aspectRatio;
            Lens* lens = new RectLens(aspect_ratio, 1);
            float distance = 0.5 / tan(0.5 * perspective.yfov);
            Camera* camera = new Camera(XYZ(), lens, XYZ(0, distance, 0));
            cameras.push_back(camera);
        }
        int node_index = -1;
        for (auto& node_data : data.nodes) {
            node_index++;
            if (node_data.camera > -1) {
                XYZ scale = XYZ(1);
                XYZ translation = XYZ();
                Quat rotation = Quat();
                if (node_data.scale.size() > 0) scale = XYZ(node_data.scale, swizzle);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                cameras[node_data.camera]->position += translation;
                cameras[node_data.camera]->rotation = rotation;
                node_lookup[node_index] = pair<int, int>(0, node_data.camera);
                continue;
            }
            if (node_data.mesh > -1) {
                XYZ scale = XYZ(1);
                XYZ translation = XYZ(0, 0, 0);
                Quat rotation = Quat(0, 0, 0, 1);
                if (node_data.scale.size() > 0) scale = XYZ(node_data.scale, scale_swizzle);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                Object* obj = new Object(XYZ(), scale);
                obj->name = node_data.name;
                obj->addMesh(meshes[node_data.mesh]);
                obj->_registerRotation(rotation);
                obj->_registerMove(translation);
                objects.push_back(obj);
                node_lookup[node_index] = pair<int, int>(1, objects.size() - 1);
                continue;
            }
            if (node_data.light > -1) {
                XYZ translation = XYZ();
                Quat rotation = Quat(0, 0, 0, 1);
                if (node_data.translation.size() > 0) translation = XYZ(node_data.translation, swizzle);
                if (node_data.rotation.size() > 0) rotation = Quat(node_data.rotation, rot_swizzle);
                Light* li = nullptr;
                auto& light_data = data.lights[node_data.light];
                if (light_data.type == "directional") {
                    li = new SunLight();
                    XYZ dir = XYZ(0, -1, 0);
                    dir = Quat::applyRotation(dir, rotation);
                    cout << dir << endl;
                    li->position = translation;
                    li->rotation = dir;
                    li->emission = XYZ(light_data.color) * light_data.intensity / 50000;
                }
                if (li != nullptr) {
                    lights.push_back(li);
                }
                else {
                    cout << "Unrecognized camera type: " << light_data.type;
                }
                node_lookup[node_index] = pair<int, int>(2, lights.size() - 1);
                continue;
            }
        }
        for (auto& object_lookup : node_lookup) {
            int node_index = object_lookup.first;
            int object_index = object_lookup.second.second;
            Object* obj = objects[object_index];
            auto& node_data = data.nodes[node_index];
            for (int child_node_index : node_data.children) {
                auto& lookup_entry = node_lookup[child_node_index];
                int child_object_index = lookup_entry.second;
                Object* child_obj = objects[child_object_index];
                obj->addChild(child_obj);
            }

        }
        int default_scene = data.defaultScene;
        auto scene_data = data.scenes[default_scene];
        //for (auto& node_index : scene_data.nodes) {
        for (int i = 0; i < node_lookup.size(); i++) { //this is 100% incorrect, but I was getting weird import issues
            auto& node_register = node_lookup[i];
            int node_type = node_register.first;
            int lookup_index = node_register.second;
            if (node_type == 0) {
                Camera* camera = cameras[lookup_index];
                scene->register_camera(camera);
            }
            if (node_type == 1) {
                Object* obj = objects[lookup_index];
                if (obj->parent == nullptr) {
                    scene->register_object(obj);
                }
            }
            if (node_type == 2) {
                Light* li = lights[lookup_index];
                scene->register_light(li);
            }
        }

        scene->camera = scene->cameras[0];
        return SM;
    }

};