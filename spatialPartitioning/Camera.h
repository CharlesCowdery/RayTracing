#pragma once


class Lens {
public:
    int resolution_x;
    int resolution_y;
    int subdivision_size;
    void (*_prep)(Lens* self, int resolution_x, int resolution_y, int subdivision_size);
    Lens* (*_clone)(Lens* self);
    XYZ(*_at)(Lens* self, int p_x, int p_y, int sample_index);
    XYZ(*_random_at)(Lens* self, int p_x, int p_y, int subindex);
    void prep(int resolution_x, int resolution_y, int subdivision_size) {
        _prep(this, resolution_x, resolution_y, subdivision_size);
    }
    XYZ at(int p_x, int p_y, int sample_index) {
        return _at(this, p_x, p_y, sample_index);
    }
    XYZ random_at(int p_x, int p_y, int subdiv_index) {
        return _random_at(this, p_x, p_y, subdiv_index);
    }
    Lens* clone() {
        Lens* obj = _clone(this);
        obj->resolution_x = resolution_x;
        obj->resolution_y = resolution_y;
        obj->subdivision_size = subdivision_size;
        return obj;
    }

};

class RectLens : public Lens {
public:
    float width;
    float height;
    vector<XY> subdiv_offsets;
    vector<vector<XY>> outputs;
    float pixel_width;
    float pixel_height;
    float subdiv_width;
    float subdiv_height;
    RectLens(float _width, float _height) :
        width(_width), height(_height) {
        _random_at = _random_at_function;
        _at = _at_function;
        _prep = _prep_function;
        _clone = _clone_function;
    }
private:
    static void _prep_function(Lens* self, int resolution_x, int resolution_y, int subdivision_count) {
        RectLens* self_rect = (RectLens*)self;
        float width = self_rect->width;
        float height = self_rect->height;
        float partial_width = width / resolution_x;
        float partial_height = height / resolution_y;
        self_rect->pixel_width = partial_width;
        self_rect->pixel_height = partial_height;
        float offset_x = -width / 2;
        float offset_y = -height / 2;
        float subdiv_width_x = partial_width / subdivision_count;
        float subdiv_width_y = partial_height / subdivision_count;
        self_rect->subdiv_width = subdiv_width_x;
        self_rect->subdiv_height = subdiv_width_y;
        for (int y = 0; y < subdivision_count; y++) {
            for (int x = 0; x < subdivision_count; x++) {
                self_rect->subdiv_offsets.push_back(XY(
                    subdiv_width_x * (x + 0.5),
                    subdiv_width_y * (y + 0.5)
                ));
            }
        }
        for (int y = 0; y < resolution_y; y++) {
            vector<XY> row;
            for (int x = 0; x < resolution_x; x++) {
                float final_x = offset_x + partial_width * x;
                float final_y = offset_y + partial_height * y;
                row.push_back(XY(
                    final_x,
                    final_y
                ));
            }
            self_rect->outputs.push_back(row);
        }

    }
    static XYZ _at_function(Lens* self, int p_x, int p_y, int sample_index) {
        RectLens* self_rect = (RectLens*)self;
        XY pos = self_rect->outputs[p_y][p_x] + self_rect->subdiv_offsets[sample_index];
        return XYZ(pos.X, 0, pos.Y);
    }
    static XYZ _random_at_function(Lens* self, int p_x, int p_y, int subindex) {
        RectLens* self_rect = (RectLens*)self;
        XY offset = self_rect->outputs[p_y][p_x];
        XY subdiv_offset = self_rect->subdiv_offsets[subindex]-XY(self_rect->subdiv_width,self_rect->subdiv_height)*0.5;
        XY pos = offset+ subdiv_offset + XY(gen.fRand(0, self_rect->subdiv_width), gen.fRand(0, self_rect->subdiv_height));
        return XYZ(pos.X, 0, pos.Y);
    }

    static Lens* _clone_function(Lens* self) {
        RectLens* rect_self = (RectLens*)self;
        Lens* obj = new RectLens(rect_self->width, rect_self->height);
        RectLens* rect_obj = (RectLens*)obj;
        rect_obj->outputs = rect_self->outputs;
        rect_obj->subdiv_offsets = rect_self->subdiv_offsets;
        return obj;
    }
};

class Camera {
    //Ill note this is a fairly non-direct class. the original was more streamlined, but I decided to update it to remove gunk
    //and to allow progressive output updates during monte carlo rendering
public:
    XYZ position;
    Quat rotation = Quat(0, 0, 0, 1);

    Lens* lens;
    XYZ focal_position;

    int current_resolution_x = 0;
    int current_resolution_y = 0;

    Camera(XYZ _position, Lens* _lens, float focal_distance) :
        lens(_lens), position(_position), focal_position(XYZ(0, 0, -focal_distance)) {}

    Camera(XYZ _position, Lens* _lens, XYZ _focal_position) :
        lens(_lens), position(_position), focal_position(_focal_position) {}

    void prep(int resolution_x, int resolution_y, int samples) {
        current_resolution_x = resolution_x;
        current_resolution_y = resolution_y;
        lens->prep(resolution_x, resolution_y, samples);
    }

    XYZ slope_at(int p_x, int p_y, int sample_index) {
        return Quat::applyRotation(XYZ::slope(XYZ(0, 0, 0), lens->at(p_x, p_y, sample_index) + focal_position), rotation);
    }

    XYZ random_slope_at(int p_x, int p_y, int subdiv_index) {
        return Quat::applyRotation(XYZ::slope(XYZ(0, 0, 0), lens->random_at(p_x, p_y, subdiv_index) + focal_position), rotation);
    }

    Camera* clone() {
        Lens* new_lens = lens->clone();
        Camera* obj = new Camera(position, new_lens, focal_position);
        obj->current_resolution_x = current_resolution_x;
        obj->current_resolution_y = current_resolution_y;
        return obj;
    }

    ~Camera() {
        delete lens;
    }
};