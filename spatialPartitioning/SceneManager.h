#pragma once

#include <chrono>
#include "Core.h"
#include "GUI.h"
#include "ImageProcessing.h"
#include "PDF.h"

using std::atomic;
using namespace std::chrono;
using std::thread;
using std::wstring;

class RenderThread {
public:
    struct block {
        int size;
        int block_index_x;
        int block_index_y;
        iXY offset;
        vector<XYZ*> raws;
        atomic<int> pixels_done;
        atomic<int> iterations_done{ 0 };
        int pixels_last;
        atomic<bool> done = false;
        block(int _block_size, int _block_x, int _block_y) {
            size = _block_size;
            block_index_x = _block_x;
            block_index_y = _block_y;
            int x_offset = block_index_x * size;
            int y_offset = block_index_y * size;
            offset = iXY(x_offset, y_offset);
            pixels_done = 0;
            pixels_last = 0;
            for (int i = 0; i < size * size; i++) {
                raws.push_back(new XYZ());
            }
        }
    };
    struct processing_info {
        processing_info() {}

        Casting_Diagnostics CD;
        atomic<long long> parent_rays_cast{ 0 };
        atomic<long long> child_rays_cast{ 0 };

        int block_size = 0;
        int res_y = 0;
        int res_x = 0;
        int y_increment = 0;
        int x_increment = 0;
        int pixels_per_block = 0;
        int samples_per_pixel = 0;

        atomic<int> pixels_done{ 0 };
        atomic<int> iterations_done{ 0 };

        Camera* camera = nullptr;
        RayEngine* RE = nullptr;

        XYZ emit_coord;
        XYZ starting_coefficient;


        void add_casts() {
            parent_rays_cast++;
            child_rays_cast += CD.rays_processed - 1;
            CD.rays_processed = 0;
        }
        long long rays_cast() {
            return parent_rays_cast + child_rays_cast;
        }
        double percent() {
            return ((double)pixels_done) / (res_x * res_y);
        }

        ~processing_info() {
            //delete camera;
        }
    };
    atomic<bool> idle = true;
    atomic<bool> stop_thread = false;
    atomic<bool> halt_ops = false;
    atomic<int> mode = 0;

    thread worker;
    block* current;
    processing_info* process_info;
    vector<block*> blocks;

    RenderThread(processing_info* config) :process_info(config), worker(&RenderThread::main, this) {

    }

    void main() {
        process_info->parent_rays_cast = 0;
        process_info->child_rays_cast = 0;
        process_info->pixels_done = 0;
        srand(high_resolution_clock::now().time_since_epoch().count());
        gen.seed(rand(), rand(), rand(), rand());
        PackagedScene* PS = process_info->RE->data;
        while (true) {
            if (stop_thread.load(std::memory_order_relaxed)) {
                return;
            }
            if (!idle.load(std::memory_order_relaxed)) {
                switch (mode.load(std::memory_order_relaxed)) {
                case 0:
                    current = (blocks.back());
                    process_linearly();
                    current->done = true;
                    idle = true;
                    break;
                case 1:
                    int block_id = 0;
                    while (!halt_ops.load(std::memory_order_relaxed)) {
                        current = blocks[block_id];
                        process_iteratively();
                        block_id++;
                        block_id %= blocks.size();
                    }
                }
            }
            std::this_thread::sleep_for(milliseconds(1));
        }
    }
    void process_linearly() {
        for (int pixel_index = 0; pixel_index < process_info->pixels_per_block; pixel_index++) {
            iXY pixel_coordinate = get_coord(pixel_index);
            XYZ* output_link = get_output_link(pixel_index);
            for (int sub_index = 0; sub_index < process_info->samples_per_pixel; sub_index++) {
                if (halt_ops.load(std::memory_order_relaxed)) {
                    return;
                };
                XYZ ray_slope = slope_at(pixel_coordinate, sub_index);
                auto ray = PackagedRay(
                    process_info->emit_coord,
                    ray_slope,
                    0
                );
                (*output_link) += process_info->starting_coefficient * process_info->RE->process_ray(process_info->CD, ray);
                process_info->add_casts();
            }
            process_info->pixels_done++;
            current->pixels_done.store(pixel_index + 1, std::memory_order_relaxed);
        }
        current->done.store(true, std::memory_order_relaxed);
    }
    void process_iteratively() {
        int iterations_done = current->iterations_done.load(std::memory_order_relaxed);
        for (int pixel_index = 0; pixel_index < process_info->pixels_per_block; pixel_index++) {
            if (halt_ops.load(std::memory_order_relaxed)) {
                return;
            };
            iXY pixel_coordinate = get_coord(pixel_index);
            XYZ* output_link = get_output_link(pixel_index);
            XYZ ray_slope = random_slope_at(pixel_coordinate);
            auto ray = PackagedRay(
                process_info->emit_coord,
                ray_slope,
                0
            );
            (*output_link) *= (double)iterations_done / (iterations_done + 1);
            (*output_link) += 1.0 / (iterations_done + 1) * process_info->RE->process_ray(process_info->CD, ray);
            process_info->add_casts();
        }
        process_info->pixels_done++;
        current->iterations_done.store(iterations_done + 1, std::memory_order_relaxed);
    }
    iXY get_coord(int pixel_index) {
        int block_size = process_info->block_size;
        iXY inside_offset = iXY(pixel_index % block_size, block_size - pixel_index / block_size - 1);
        return current->offset + inside_offset;
    }
    XYZ* get_output_link(int pixel_index) {
        return current->raws[pixel_index];
    }
    XYZ slope_at(iXY coord, int sub_index) {
        return process_info->camera->slope_at(coord.X, coord.Y, sub_index);
    }
    XYZ random_slope_at(iXY coord) {
        return process_info->camera->random_slope_at(coord.X, coord.Y);
    }



};

class SceneManager {
public:
    Scene* scene;
    RayEngine RE;
    GUIHandler GUI;
    vector<vector<XYZ*>> raw_output;

    int thread_count = 12;
    vector<RenderThread*> threads;

    int current_resolution_x = 0;
    int current_resolution_y = 0;
    int subdivision_count = 0;
    int current_samples_per_pixel = 0;
    const int block_size = 8;

    int y_increment;
    int x_increment;

    steady_clock::time_point render_start;

    SceneManager(Scene* _scene) : scene(_scene) {}

    void render(int resolution_x, int resolution_y, int _subdivision_count = 1, int mode = 0) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;

        y_increment = ceil(current_resolution_y / block_size);
        x_increment = ceil(current_resolution_x / block_size);

        subdivision_count = _subdivision_count;
        current_samples_per_pixel = subdivision_count * subdivision_count;

        GUI = GUIHandler(current_resolution_x, current_resolution_y, block_size);

        render_start = high_resolution_clock::now();
        prep();
        //enqueue_rays();
        GUI.create_window();
        if (mode == 0) prepGUI();
        cout << endl << "+RAYCASTING+" << endl;
        spawn_threads();
        run_threaded_engine(mode);
        //_run_engine(true);

    }
    void hold_window() {
        GUI.hold_window();
    }

private:
    void prep() {
        cout << endl << "+PREPROCESSING+" << endl;

        auto prep_start = high_resolution_clock::now();

        for (int y = 0; y < current_resolution_y; y++) {
            vector<XYZ*> row;
            for (int x = 0; x < current_resolution_x; x++) {
                row.push_back(new XYZ(0, 0, 0));
            }
            raw_output.push_back(row);
        }

        scene->prep(current_resolution_x, current_resolution_y, subdivision_count);
        PDF::prep();
        RE.load_scene(scene->package());

        auto prep_end = high_resolution_clock::now();

        cout << "Preprocessing done [" << duration_cast<milliseconds>(prep_end - prep_start).count() << "ms]" << endl;

    }
    void spawn_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread::processing_info* config = create_thread_config();
            RenderThread* thread = new RenderThread(config);
            string name = "RenderingThread-" + to_string(i);
            wstring wname = wstring(name.begin(), name.end());
            const wchar_t* wcname = wname.c_str();
            threads.push_back(thread);
            SetThreadPriority(thread->worker.native_handle(), 0);
            SetThreadDescription(thread->worker.native_handle(), wcname);

        }
    }
    struct render_stats {
        steady_clock::time_point start_time;
        steady_clock::time_point now;
        steady_clock::time_point last_refresh = high_resolution_clock::now();
        long long parent_rays = 0;
        long long child_rays = 0;
        long long parent_rays_last = 0;
        long long child_rays_last = 0;
        double percent_last = 0;
        long long pixels_done = 0;
        double percent_done = 0;
        vector<long long> parent_rays_rate_history;
        vector<long long> child_rays_rate_history;
        vector<double>    percent_rate_history;
        int remaining_threads;
        render_stats(int history_size_rays, int history_size_rate) {
            parent_rays_rate_history = vector<long long>(history_size_rays, 0);
            child_rays_rate_history = vector<long long>(history_size_rays, 0);
            percent_rate_history = vector<double>(history_size_rate, 0);
        }
    };
    struct pixel_strip {
        int strip_start;
        iXY block;
        vector<XYZ> outputs;
        pixel_strip(int start, iXY _block) {
            strip_start = start;
            block = _block;
        }
    };
    void assign_iterative_groups(vector<RenderThread::block*>& job_stack) {
        double jobs_per = (double)job_stack.size() / threads.size();
        double amount_given = 0;
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->mode = 1;
            double num_now = amount_given + jobs_per;
            double to_do = num_now - ceil(amount_given);
            amount_given = num_now;
            if (i == thread_count - 1) num_now = ceil(num_now);
            for (int j = 0; j < to_do; j++) {
                int job_index = rand() % job_stack.size();
                RenderThread::block* job = job_stack[job_index];
                render_thread->blocks.push_back(job);
                job_stack.erase(job_stack.begin() + job_index);
            }
        }
    }
    void start_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->idle = false;
        }
    }
    void halt_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->halt_ops = true;
        }
    }
    void unhalt_threads() {
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            render_thread->halt_ops = false;
        }
    }
    void reset_blocks(vector<RenderThread::block*>& blocks) {
        for (int i = 0;i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            block->iterations_done = 0;
            for (int j = 0; j < block->raws.size(); j++) {
                *block->raws[j] = XYZ();
            }
        }
    }
    void gather_thread_stats(render_stats& rStat) {
        rStat.parent_rays = 0;
        rStat.child_rays = 0;
        rStat.pixels_done = 0;
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            rStat.parent_rays += render_thread->process_info->parent_rays_cast;
            rStat.child_rays += render_thread->process_info->child_rays_cast;
            rStat.pixels_done += render_thread->process_info->pixels_done;
        }
    }
    void manage_threads(render_stats& rStat, vector<RenderThread::block*>& job_stack) {
        steady_clock::time_point now = high_resolution_clock::now();
        for (int i = 0; i < thread_count; i++) {
            RenderThread* render_thread = threads[i];
            bool is_idle = render_thread->idle.load(std::memory_order_relaxed);
            bool is_stopped = render_thread->stop_thread.load(std::memory_order_relaxed);
            if (is_idle && !is_stopped) {
                if (job_stack.size() > 0) {
                    RenderThread::block* job_block = job_stack.back();
                    GUI.set_focus_position(i, job_block->offset.X, job_block->offset.Y);
                    render_thread->blocks.push_back(job_block);
                    job_stack.pop_back();
                    render_thread->idle = false;
                }
                else {
                    render_thread->stop_thread = true;
                    rStat.remaining_threads--;
                    GUI.set_focus_position(i, -10000, -10000);
                }
            }
        }
    }
    iXY get_pixel_offset(int index) {
        return iXY(index % block_size, block_size - index / block_size - 1);
    }
    void get_block_updates(vector<RenderThread::block*>& blocks, vector<pixel_strip>& outputs) {
        for (int i = 0; i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            int pixels_done = block->pixels_done;
            int pixels_last = block->pixels_last;
            if (pixels_last < pixels_done) {
                iXY block_offset = block->offset;
                pixel_strip p_strip(pixels_last, block_offset);
                for (int pixel_index = pixels_last; pixel_index < pixels_done; pixel_index++) {
                    iXY internal_offset = get_pixel_offset(pixel_index);
                    iXY final_position = block_offset + internal_offset;
                    (*raw_output[final_position.Y][final_position.X]) = *block->raws[pixel_index];
                    XYZ color = ImageHandler::post_process_pixel(*raw_output[final_position.Y][final_position.X]);
                    p_strip.outputs.push_back(color);
                    //GUI.commit_pixel(color, final_position.X, final_position.Y);

                }
                outputs.push_back(p_strip);
            }
            block->pixels_last = pixels_done;
        }
    }
    void get_screen_update(vector<RenderThread::block*>& blocks) {
        for (int i = 0;i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            iXY block_offset = block->offset;
            for (int pixel_index = 0; pixel_index < block_size * block_size; pixel_index++) {
                iXY internal_offset = get_pixel_offset(pixel_index);
                iXY final_coord = block_offset + internal_offset;
                (*raw_output[final_coord.Y][final_coord.X]) = *block->raws[pixel_index];
                XYZ color = ImageHandler::post_process_pixel(*raw_output[final_coord.Y][final_coord.X]);
                GUI.commit_pixel(color, final_coord.X, final_coord.Y);
            }
        }
    }
    void pre_process_stats(render_stats& rStat) {
        int seconds = duration_cast<std::chrono::seconds>(rStat.now - rStat.start_time).count();
        int minutes = seconds / 60;
        int hours = minutes / 60;
        seconds %= 60;
        minutes %= 60;
        int micro_delta = duration_cast<microseconds>(rStat.now - rStat.last_refresh).count();
        int milli_delta = micro_delta / 1000;
        double second_ratio = ((double)1000000) / micro_delta;

        long long total_pixels = current_resolution_x * current_resolution_y;
        double percent = ((double)rStat.pixels_done) / total_pixels;
        rStat.percent_done = percent * 100;
        double percent_rate = (rStat.percent_done - rStat.percent_last) * second_ratio;


        long long delta_parent = rStat.parent_rays - rStat.parent_rays_last;
        long long delta_child = rStat.child_rays - rStat.child_rays_last;

        long long parent_rate = delta_parent * second_ratio;
        long long child_rate = delta_child * second_ratio;

        rStat.percent_rate_history.push_back(percent_rate);
        rStat.percent_rate_history.erase(rStat.percent_rate_history.begin());
        rStat.parent_rays_rate_history.push_back(parent_rate);
        rStat.parent_rays_rate_history.erase(rStat.parent_rays_rate_history.begin());
        rStat.child_rays_rate_history.push_back(child_rate);
        rStat.child_rays_rate_history.erase(rStat.child_rays_rate_history.begin());
    }
    void post_process_stats(render_stats& rStat) {
        rStat.parent_rays_last = rStat.parent_rays;
        rStat.child_rays_last = rStat.child_rays;
        rStat.percent_last = rStat.percent_done;
        rStat.last_refresh = high_resolution_clock::now();
    }
    void display_stats(render_stats& rStat) {
        int seconds = duration_cast<std::chrono::seconds>(rStat.now - rStat.start_time).count();
        int minutes = seconds / 60;
        int hours = minutes / 60;
        seconds %= 60;
        minutes %= 60;
        int micro_delta = duration_cast<std::chrono::microseconds>(rStat.now - rStat.last_refresh).count();
        int milli_delta = micro_delta / 1000;
        double second_ratio = ((double)1000000) / micro_delta;

        double percent_rate_average = 0;
        long long parent_rate_average = 0;
        long long child_rate_average = 0;

        for (int i = 0; i < rStat.parent_rays_rate_history.size(); i++) {
            parent_rate_average += rStat.parent_rays_rate_history[i];
        }
        parent_rate_average /= rStat.parent_rays_rate_history.size();
        for (int i = 0; i < rStat.child_rays_rate_history.size(); i++) {
            child_rate_average += rStat.child_rays_rate_history[i];
        }
        child_rate_average /= rStat.child_rays_rate_history.size();
        for (int i = 0; i < rStat.percent_rate_history.size(); i++) {
            percent_rate_average += rStat.percent_rate_history[i];
        }
        percent_rate_average /= rStat.percent_rate_history.size();

        string primary_per_sec = intToEng(parent_rate_average);
        string secondary_per_sec = intToEng(child_rate_average);
        string total_per_sec = intToEng(parent_rate_average + child_rate_average);



        string format_string = "[%02i:%02i:%02i][%7ims] - ([%.5srays/s Primary][%.5srays/s Secondary])[%.5srays/s] - [%4.1f%% (+%4.3f/s)]\r";
        printf(format_string.c_str(), hours, minutes, seconds, milli_delta, primary_per_sec.c_str(), secondary_per_sec.c_str(), total_per_sec.c_str(), rStat.percent_done, percent_rate_average);

    }
    void process_strips(vector<pixel_strip> p_strips) {
        for (int o_index = 0; o_index < p_strips.size(); o_index++) {
            pixel_strip& p_strip = p_strips[o_index];
            iXY block_offset = p_strip.block;
            int strip_start = p_strip.strip_start;
            vector<XYZ>& outputs = p_strip.outputs;
            for (int i = 0; i < outputs.size(); i++) {
                XYZ color = outputs[i];
                iXY internal_offset = get_pixel_offset(i + strip_start);
                iXY final_coord = internal_offset + block_offset;
                GUI.commit_pixel(color, final_coord.X, final_coord.Y);
            }

        }
        p_strips.clear();
    }
    void refresh_display() {
        GUI.commit_canvas();
        GUI.handle_events();
    }
    void run_threaded_engine(int mode = 0) {
        render_stats rStat(120, 300);
        rStat.start_time = high_resolution_clock::now();
        render_start = rStat.start_time;
        rStat.remaining_threads = thread_count;
        vector<RenderThread::block*> blocks;
        vector<RenderThread::block*> job_stack;
        vector<pixel_strip> update_strips;
        for (int block_index_y = 0; block_index_y < y_increment; block_index_y++) {
            for (int block_index_x = 0; block_index_x < x_increment; block_index_x++) {
                job_stack.push_back(new RenderThread::block(block_size, block_index_x, block_index_y));
                blocks.push_back(job_stack.back());
            }
        }
        switch (mode) {
        case 0:
            while (true) {
                rStat.now = high_resolution_clock::now();
                manage_threads(rStat, job_stack);
                if (duration_cast<std::chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 16) {
                    gather_thread_stats(rStat);
                    get_block_updates(blocks, update_strips);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    process_strips(update_strips);
                    refresh_display();
                }
                if (rStat.remaining_threads <= 0) {
                    break;
                }
            }
            break;
        case 1:
            assign_iterative_groups(job_stack);
            start_threads();
            while (true) {
                rStat.now = high_resolution_clock::now();
                if (duration_cast<std::chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 16) {
                    gather_thread_stats(rStat);
                    get_screen_update(blocks);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    refresh_display();
                }

            }
            break;
        case 2:
            assign_iterative_groups(job_stack);
            start_threads();
            while (true) {
                rStat.now = high_resolution_clock::now();
                if (duration_cast<std::chrono::milliseconds>(rStat.now - rStat.last_refresh).count() > 64) {
                    halt_threads();
                    gather_thread_stats(rStat);
                    get_screen_update(blocks);
                    pre_process_stats(rStat);
                    display_stats(rStat);
                    post_process_stats(rStat);
                    refresh_display();
                    reset_blocks(blocks);
                    unhalt_threads();
                }
            }
        }
        for (int i = 0; i < thread_count; i++) {
            threads[i]->stop_thread = true;
            threads[i]->worker.join();
        }
        for (int i = 0; i < blocks.size(); i++) {
            RenderThread::block* block = blocks[i];
            iXY block_offset = block->offset;
            for (int pixel_index = 0; pixel_index < block->raws.size(); pixel_index++) {
                iXY internal_offset = iXY(pixel_index % block_size, block_size - pixel_index / block_size - 1);
                iXY final_position = block_offset + internal_offset;
                (*raw_output[final_position.Y][final_position.X]) = *block->raws[pixel_index];
                XYZ color = ImageHandler::post_process_pixel(*raw_output[final_position.Y][final_position.X]);
                GUI.commit_pixel(color, final_position.X, final_position.Y);
            }
        }
        refresh_canvas();
        cout << endl;

        auto end_time = high_resolution_clock::now();
        double millis = duration_cast<std::chrono::milliseconds>(end_time - rStat.start_time).count();
        cout << "Ray Casting Done [" << to_string(millis) << "ms]" << endl;
        //cout << "Approximate speed: " << intToEng((rays_cast / (millis / 1000.0))) << " rays/s";
    }
    void prepGUI() {
        for (int i = 0; i < thread_count; i++) {
            GUI.add_focus();
        }
    }
    RenderThread::processing_info* create_thread_config() {
        RenderThread::processing_info* process_info = new RenderThread::processing_info();

        int max_bounces = 2;

        process_info->emit_coord = scene->camera->position;
        process_info->starting_coefficient = XYZ(1, 1, 1) / current_samples_per_pixel;

        process_info->res_y = current_resolution_y;
        process_info->res_x = current_resolution_x;
        process_info->y_increment = y_increment;
        process_info->x_increment = x_increment;
        process_info->block_size = block_size;
        process_info->samples_per_pixel = current_samples_per_pixel;
        process_info->pixels_per_block = block_size * block_size;

        process_info->camera = scene->camera;//camera->clone();
        process_info->RE = new RayEngine(RE);

        return process_info;
    }


    void refresh_canvas() {
        for (int y = 0; y < current_resolution_y; y++) {
            for (int x = 0; x < current_resolution_x; x++) {
                XYZ color = ImageHandler::post_process_pixel(*raw_output[y][x]);
                GUI.commit_pixel(color, x, y);
            }
        }
        GUI.commit_canvas();
    }
};