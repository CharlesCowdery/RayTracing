#pragma once

#include <SFML/Graphics.hpp>
#include "XYZ.h"
#include <chrono>

using namespace std::chrono;
using std::string;
using std::to_string;

class GUIHandler {
public:
    struct noise_rect {
        iXY begin;
        iXY end;
        float noise;
        bool done = false;
        float max_noise = 1;
        Quat get_color() {
            float r_noise = std::max(0.0f,std::min(noise, 1.0f));
            XYZ c_color;
            if (done)c_color = XYZ(0, 255, 0);
            else c_color = XYZ(0, 0, 255);
            XYZ n_color = XYZ(255, 0, 0);
            float f = r_noise/max_noise;
            f = std::min(1.0f, std::max(f, 0.0f));
            int alpha = (done) ? 0 : 100;
            return Quat((1 - f) * c_color + (f)*n_color, alpha);
        }
        sf::RectangleShape R;
    };
    int current_resolution_x = 0;
    int current_resolution_y = 0;

    sf::RenderWindow* window;
    sf::Image canvas;
    sf::Texture render_texture;
    sf::Sprite render_sprite;

    sf::RenderWindow* noise_window;
    sf::Texture noise_texture;
    sf::Sprite  noise_sprite;

    sf::RenderWindow* denoised_window;
    sf::Texture denoised_texture;
    sf::Sprite  denoised_sprite;
    sf::Image denoised_canvas;

    map<string, sf::Texture> textures;
    map<string, sf::Sprite> sprites;
    map<string, sf::Image> images;
    map<string, sf::RectangleShape> rects;
    map<string, sf::RectangleShape> noise;

    double scalar_exponent = 0;

    XY focus_size;
    vector<XY> focuses;

    vector<noise_rect> noise_rects;

    GUIHandler() {}

    GUIHandler(int resX, int resY, int block_size) {
        current_resolution_x = resX;
        current_resolution_y = resY;
        canvas.create(current_resolution_x, current_resolution_y, sf::Color::Black);
        render_sprite.setScale(make_scale(), -make_scale());
        render_sprite.setPosition(0, make_scale() * current_resolution_y);
        focus_size = XY(block_size, block_size);
        denoised_canvas.create(current_resolution_x, current_resolution_y, sf::Color::Black);
        denoised_sprite.setScale(make_scale(), -make_scale());
        denoised_sprite.setPosition(0, make_scale() * current_resolution_y);
    }

    void hold_window() {
        auto frame_time = high_resolution_clock::now();
        while (window->isOpen())
        {
            sf::Event event;
            while (window->pollEvent(event))
            {
                handle_events(event);
            }
            partial_render();
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }
    }
    XY get_relative_position(int x, int y) {
        double relative_x = ((double)x) / current_resolution_x;
        double relative_y = ((double)y) / current_resolution_y;
        return XY(make_scale() * relative_x, make_scale() * relative_y);

    }
    void create_window() {
        cout << "Opening Window...." << flush;
        window = new sf::RenderWindow(sf::VideoMode(current_resolution_x * make_scale(), current_resolution_y * make_scale()), "Render Window!");
        cout << "Done" << endl;
        denoised_window = new sf::RenderWindow(sf::VideoMode(current_resolution_x * 1, current_resolution_y * 1), "denoised Window!");
    }
    void create_noise_window() {
        noise_window = new sf::RenderWindow(sf::VideoMode(current_resolution_x * 1, current_resolution_y * 1), "Noise Window!");
    }
    void commit_pixel(XYZ color, int x, int y) {
        auto sfcolor = sf::Color(color.X, color.Y, color.Z);
        canvas.setPixel(x, y, sfcolor);
    }
    void commit_denoised_pixel(XYZ color, int x, int y) {
        auto sfcolor = sf::Color(color.X, color.Y, color.Z);
        denoised_canvas.setPixel(x, y, sfcolor);
    }
    void partial_render() {
        window->clear();
        draw_render_sprite();
        window->display();
    }
    void render_pass() {
        window->clear();
        update_focuses();
        draw_render_sprite();
        draw_secondary_sprites();
        draw_shapes();
        window->display();
    }
    void noise_pass() {
        noise_window->clear();
        noise_window->draw(render_sprite);
        draw_noise_rects();
        noise_window->display();
    }
    void denoised_pass() {
        denoised_window->clear();
        draw_denoised_sprite();
        denoised_window->display();
    }
    void update_focuses() {
        XY scaled_size = focus_size * make_scale();
        for (int i = 0; i < focuses.size(); i++) {
            XY focus_position = focuses.at(i);
            focus_position.Y = current_resolution_y - focus_position.Y - focus_size.Y;
            XY scaled_position = focus_position * make_scale();
            string focus_id = "render_focus_" + to_string(i);
            rects[focus_id].setSize(sf::Vector2f(scaled_size.X, scaled_size.Y));
            rects[focus_id].setPosition(scaled_position.X, scaled_position.Y);
        }
    }
    void draw_render_sprite() {
        window->draw(render_sprite);
    }
    void draw_denoised_sprite() {
        denoised_window->draw(denoised_sprite);
    }
    void draw_secondary_sprites() {
        for (auto i = sprites.begin(); i != sprites.end(); i++) {
            sf::Sprite& sprite = i->second;
            window->draw(sprite);

        }
    }
    void draw_shapes() {
        for (auto i = rects.begin(); i != rects.end(); i++) {
            sf::RectangleShape& rect = i->second;
            window->draw(rect);
        }
    }
    void draw_noise_rects() {
        for (int i = 0; i < noise_rects.size(); i++) {
            noise_rect rect = noise_rects[i];
            iXY delta = rect.end - rect.begin;
            Quat color = rect.get_color();
            rect.R.setSize(sf::Vector2f(delta.X,delta.Y));
            rect.R.setPosition(rect.begin.X, (current_resolution_y-rect.begin.Y-delta.Y));
            //rect.R.setPosition(100, 100);
            rect.R.setFillColor(sf::Color(color.X, color.Y, color.Z, color.W));
            //rect.R.setFillColor(sf::Color(255, 255, 255, 255));
            //rect.R.setOutlineColor(sf::Color(color.X, color.Y, color.Z, 200));
            noise_window->draw(rect.R);
        }
    }
    void update_image_textures() {
        for (auto i = images.begin(); i != images.end(); i++) {
            string key = i->first;
            sf::Image& image = i->second;
            getTexture(key).loadFromImage(image);
        }
    }
    void commit_canvas() {
        render_texture.loadFromImage(canvas);
        render_sprite.setTexture(render_texture, false);
        render_pass();
    }
    void commit_denoised_canvas() {
        denoised_texture.loadFromImage(denoised_canvas);
        denoised_sprite.setTexture(denoised_texture, false);
        denoised_pass();
    }
    void update_noise() {
        noise_pass();
    }

    void add_texture_sprite(string name) {
        textures[name] = sf::Texture();
        sprites[name] = sf::Sprite();
        getSprite(name).setTexture(getTexture(name));
    }
    void add_texture_sprite_image(string name) {
        textures[name] = sf::Texture();
        sprites[name] = sf::Sprite();
        images[name] = sf::Image();
        getSprite(name).setTexture(getTexture(name));
    }
    sf::RectangleShape& add_rect(string name) {
        rects[name] = sf::RectangleShape();
        return rects[name];
    }
    void add_focus() {
        string focus_id = "render_focus_" + to_string(focuses.size());
        add_rect(focus_id);
        focuses.push_back(XY(0, 0));

        sf::RectangleShape& focus = rects[focus_id];
        focus.setFillColor(sf::Color::Transparent);
        focus.setOutlineThickness(1);
        focus.setOutlineColor(sf::Color(255, 211, 110));
    }
    void set_focus_position(int index, int x, int y) {
        focuses[index] = XY(x, y);
    }
    sf::Sprite& getSprite(string name) {
        return sprites[name];
    }
    sf::Texture& getTexture(string name) {
        return textures[name];
    }
    sf::Image& getImage(string name) {
        return images[name];
    }
    void handle_events() {
        sf::Event event;
        while (window->pollEvent(event)) {
            handle_events(event);
        }
    }
    void handle_events_noise() {
        sf::Event event;
        while (noise_window->pollEvent(event)) {
            if (event.type == sf::Event::Resized)
            {
                // update the view to the new size of the window
                sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
                noise_window->setView(sf::View(visibleArea));
            }
        }
    }
    void handle_events_denoised() {
        sf::Event event;
        while (denoised_window->pollEvent(event)) {
            if (event.type == sf::Event::Resized)
            {
                // update the view to the new size of the window
                sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
                denoised_window->setView(sf::View(visibleArea));
            }
        }
    }
    double make_scale() {
        return pow(2, scalar_exponent);
    }
    void handle_events(sf::Event event) {
        if (event.type == sf::Event::Closed) window->close();
        if (event.type == sf::Event::Resized)
        {
            // update the view to the new size of the window
            sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
            window->setView(sf::View(visibleArea));
        }
        if (event.type == sf::Event::MouseWheelMoved)
        {
            // display number of ticks mouse wheel has moved
            double amount = event.mouseWheel.delta;
            scalar_exponent += amount / 5;
            double scale = make_scale();
            render_sprite.setScale(scale, -scale);
            render_sprite.setPosition(0, make_scale() * current_resolution_y);
        }
    }
};
