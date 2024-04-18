#pragma once

#include <queue>
#include <set>

#include "commons.h"
#include "Primitives.h"

using std::queue;
using std::set;
using std::vector;
using std::exception;
using std::cout;
using std::endl;
using std::pair;

class BVH {
public:
    XYZ max = XYZ(0, 0, 0);
    XYZ min = XYZ(0, 0, 0);
    vector<Tri*> elements;
    BVH* c1 = nullptr;
    BVH* c2 = nullptr;
    BVH* parent;
    int _count = -1;
    BVH() {    }
    BVH(vector<Tri*>* T_vec) {
        elements = vector<Tri*>(*T_vec);
        for (Tri* T_ptr : *T_vec) {
            max = XYZ::max(max, T_ptr->AABB_max);
            min = XYZ::min(min, T_ptr->AABB_min);
        }

    }
    ~BVH() {
        delete c1;
        delete c2;
    }
    void printDesmos() {
        //cout << "B(" + min.to_string() + "," + max.to_string() + ")" << endl;
        if (c1 != nullptr) c1->printDesmos();
        if (c2 != nullptr) c2->printDesmos();
        for (auto& T : elements) {
            cout << "triangle(" + T->p1.to_string() + "," + T->p2.to_string() + "," + T->p3.to_string() + ")" << endl;
        }
    }
    static float intersection(const XYZ& max, const XYZ& min, const XYZ& origin, const XYZ& inv_slope) {
        float tx1 = (min.X - origin.X) * inv_slope.X;
        float tx2 = (max.X - origin.X) * inv_slope.X;

        float tmin = std::min(tx1, tx2);
        float tmax = std::max(tx1, tx2);

        float ty1 = (min.Y - origin.Y) * inv_slope.Y;
        float ty2 = (max.Y - origin.Y) * inv_slope.Y;

        tmin = std::max(tmin, std::min(ty1, ty2));
        tmax = std::min(tmax, std::max(ty1, ty2));

        float tz1 = (min.Z - origin.Z) * inv_slope.Z;
        float tz2 = (max.Z - origin.Z) * inv_slope.Z;

        tmin = std::max(tmin, std::min(tz1, tz2));
        tmax = std::min(tmax, std::max(tz1, tz2));

        if (tmax >= tmin) { //introduces a branch, but itll probably be fine. Lets me order search checks by nearest intersection
            if (tmin >= 0 && tmin < tmax) {
                return tmin;
            }
            else {
                return tmax;
            }
        }
        else {
            return -1;
        }
    }
    float intersection(const XYZ& origin, const XYZ& inv_slope) { //https://tavianator.com/2011/ray_box.html
        return BVH::intersection(max, min, origin, inv_slope);
    }

    struct WorkPacket {
        BVH* target;
        vector<Tri*>* contents;
        WorkPacket(BVH* t, vector<Tri*>* c) : target(t), contents(c) {}
    };
    //very dangerous. only ever use if BVH is being used as an intermediate
    //will turn every tri into a hanging pointer if the passed vector moves out of scope
    int construct_dangerous(vector<Tri>& initial_geo) {
        vector<Tri*>* working_data = new vector<Tri*>();
        working_data->reserve(initial_geo.size());
        for (Tri& T : initial_geo) {
            working_data->push_back(&T);
        }
        int return_value = construct(working_data);
        return return_value;
    }
    int construct(vector<Tri*>* initial_geo) {
        queue<WorkPacket> packets;
        assertm((initial_geo->size() > 0), "BVH was creation was attempted with zero geometry. Did everything load right?");
        auto initial_packet = WorkPacket(this, initial_geo);
        packets.push(initial_packet);
        int i = 0;
        while (packets.size() > 0) {
            i++;
            WorkPacket packet = packets.front();
            vector<Tri*>* geo = packet.contents;
            BVH* target = packet.target;

            assertm(geo->size() > 0, "zero size bin");

            AABB my_bounds = getBounds(geo);

            target->max = my_bounds.max;
            target->min = my_bounds.min;

            XYZ extent = target->max - target->min;
            int axis = 0;
            if (extent.Y > extent.X) axis = 1;
            if (extent.Z > extent[axis]) axis = 2;
            float splitPos = target->min[axis] + extent[axis] * 0.5f;

            Split* split = new Split(splitPos, axis);
            get_stats(geo, *split);
            evaluate_split(*split, 1);
            Split* probe_split = probe(geo);
            get_stats(geo, *probe_split);
            evaluate_split(*probe_split, 1);


            if (probe_split->score < split->score) {
                delete split;
                split = probe_split;
            }
            auto bins = bin(split, geo);
            //if (geo->size() > 200) {
            //    auto binned_split = binned_split_probe(geo);
            //    if (binned_split.first->score < split->score) {
            //        bins = binned_split.second;
            //        split = binned_split.first;
            //    }
            //}


            if (((size_t)abs((int)(bins.p_geo->size() - bins.n_geo->size()))) == geo->size()) {
                target->elements = vector<Tri*>(*geo);
                packets.pop();
                continue;
            }

            if (bins.p_geo->size() <= LEAF_SIZE) {
                assertm(bins.p_geo->size() > 0, "attempted bin creation of zero elements");
                target->c1 = new BVH(bins.p_geo);
                target->c1->parent = target;
            }
            else {
                BVH* leaf = new BVH();
                target->c1 = leaf;
                target->c1->parent = target;
                WorkPacket k = WorkPacket(leaf, bins.p_geo);
                packets.push(k);
            }
            if (bins.n_geo->size() <= LEAF_SIZE) {
                assertm(bins.n_geo->size() > 0, "attempted bin creation of zero elements");
                target->c2 = new BVH(bins.n_geo);
                target->c2->parent = target;
            }
            else {
                BVH* leaf = new BVH();
                target->c2 = leaf;
                target->c2->parent = target;
                WorkPacket k = WorkPacket(leaf, bins.n_geo);
                packets.push(k);
            }

            delete geo;
            packets.pop();
        }
        return i;
    }
    pair<vector<BVH*>*, vector<Tri*>*> flatten(int depth = 0) {
        vector<BVH*>* l = new vector<BVH*>();
        vector<Tri*>* t = new vector<Tri*>();

        flatten(l, t, depth);
        return pair<vector<BVH*>*, vector<Tri*>*>(l, t);
    }
    void flatten(vector<BVH*>* l, vector<Tri*>* t, int depth) {
        depth = depth - 1;
        l->push_back(this);
        if (depth != 0) {
            if (c1 != nullptr) {
                c1->flatten(l, t, depth);
            }
            if (c2 != nullptr) {
                c2->flatten(l, t, depth);
            }
            if (elements.size() > 0) {
                for (auto T : elements) {
                    t->push_back(T);
                }
            }
        }
    }
    int size() {
        int s = 0;
        size(s);
        return s;
    }
    void size(int& s) {
        s++;
        if (c1 != nullptr) {
            c1->size(s);
            c2->size(s);
        }
    }
    int count() {
        if (_count == -1) {
            count(_count);
        }
        return _count;
    }
    void count(int& s) const {
        s += this->elements.size();
        if (c1 != nullptr) {
            c1->count(s);
            c2->count(s);
        }
    }
    float avg_path() const {
        if (c1 == nullptr) return 0;
        return (c1->avg_path() + c2->avg_path()) / 2.0;
    }
    vector<PackagedTri> get_emissive_tris() {
        vector<PackagedTri> out;
        get_emissive_tris(out);
        return out;
    }
    void get_emissive_tris(vector<PackagedTri>& vec) {
        if (c1 != nullptr) {
            c1->get_emissive_tris(vec);
            c2->get_emissive_tris(vec);
        }
        else {
            for (int i = 0; i < elements.size(); i++) {
                Tri* t = elements[i];
                if (t->material->emissive.getValue().magnitude() > 0) {
                    vec.push_back(t->pack());
                }
            }
        }
    }

private:
    struct BinResults {
        vector<Tri*>* p_geo;
        vector<Tri*>* n_geo;
    };
    struct Stats {
        XYZ min = XYZ(99999999999);
        XYZ max = XYZ();
        int count = 0;
        float SA() {
            return VecLib::surface_area(max, min);
        }
        float volume() {
            return VecLib::volume(max, min);
        }
    };
    struct Split {
        Stats p;
        Stats n;
        float score;
        XYZ placement;
        int facing;
        Split(XYZ place, int _facing) :placement(place), facing(_facing) {};
    };
    struct AABB {
        XYZ max;
        XYZ min;
        AABB() : max(-99999999999999), min(99999999999999) {}
        AABB(XYZ _max, XYZ _min) :max(_max), min(_min) {}
        void expand(XYZ point) {
            max = XYZ::max(max, point);
            min = XYZ::min(min, point);
        }
    };
    struct reduced_tri {
        XYZ aabb_max;
        XYZ aabb_min;
        XYZ midpoint;
        reduced_tri(Tri T) {
            aabb_max = T.AABB_max;
            aabb_min = T.AABB_min;
            midpoint = T.midpoint;
        }
        struct less_than_axis_operator {
            int axis;
            less_than_axis_operator(int _axis) {
                axis = _axis;
            }
            inline bool operator() (const reduced_tri t1, const reduced_tri t2)
            {
                return (t1.midpoint[axis] < t2.midpoint[axis]);
            }
        };
    };
    AABB getBounds(vector<Tri*>* geo) {
        AABB results;
        for (Tri* T_ptr : *geo) {
            results.expand(T_ptr->AABB_max);
            results.expand(T_ptr->AABB_min);
        }
        return results;
    }
    BinResults bin(Split* split, vector<Tri*>* geo) {
        BinResults results;
        results.p_geo = new vector<Tri*>();
        results.n_geo = new vector<Tri*>();
        int axis = split->facing;
        for (Tri* T_ptr : *geo) {
            if (T_ptr->midpoint[axis] > split->placement[axis]) {
                results.p_geo->push_back(T_ptr);
            }
            else {
                results.n_geo->push_back(T_ptr);
            }
        }
        return results;
    }
    struct relevant_value {
        XYZ midpoint;
        XYZ max;
        XYZ min;
        XYZ right_max;
        XYZ right_min;
        Tri* parent;
        int value_selector = 0;
        relevant_value(Tri* source, int v_selector = 0) {
            parent = source;
            midpoint = parent->midpoint;
            max = parent->AABB_max;
            min = parent->AABB_min;
            value_selector = v_selector;
        }
        XYZ get_value() const {
            switch (value_selector) {
            case 0:
                return min;
            case 1:
                return midpoint;
            case 2:
                return max;
            default:
                throw exception("out of bounds value selector");
            }
        }
        struct comparer {
        public:
            comparer(char compare_index, bool _use_selector = false) {
                internal_index = compare_index;
                use_selector = _use_selector;
            }
            inline bool operator()(const relevant_value& v1, const relevant_value& v2) {
                if (use_selector) {
                    XYZ v1_value = v1.get_value();
                    XYZ v2_value = v2.get_value();
                    return v1_value[internal_index] < v2_value[internal_index];
                }
                return v1.midpoint[internal_index] < v2.midpoint[internal_index];
            }
        private:
            char internal_index;
            bool use_selector;
        };
    };
    pair<Split*, BinResults> binned_split_probe(vector<Tri*>* geo) {
        int min_bins = 10;
        int max_bins = 1000;
        int bin_count = std::max(min_bins, std::min(max_bins, (int)geo->size()));
        int adaptive_bin_count = 3;
        int adpative_sweep_recusion_count = 1;

        int geo_count = geo->size();
        auto sorted = vector<relevant_value>();

        Split operator_split = Split(XYZ(), 0);
        Split best = operator_split;
        best.score = 999999999999999999999999.0;

        for (Tri* T_ptr : *geo) {
            sorted.push_back(relevant_value(T_ptr, 0));
            //sorted.push_back(relevant_value(T_ptr, 1));
            sorted.push_back(relevant_value(T_ptr, 2));
        }
        for (int i = 0; i < 3; i++) {
            operator_split.facing = i;
            sort(sorted.begin(), sorted.end(), relevant_value::comparer(i, true));

            XYZ left_min = sorted.front().min;
            XYZ left_max = sorted.front().min;
            XYZ right_min = sorted.back().min;
            XYZ right_max = sorted.back().max;
            int left_count = 0;
            int right_count = geo_count;



            for (int j = sorted.size() - 1; j >= 0; j--) {
                relevant_value& v = sorted[j];
                v.right_min = right_min;
                v.right_max = right_max;
                if (v.value_selector == 0) {
                    right_min = XYZ::min(right_min, v.min);
                    right_max = XYZ::max(right_max, v.max);
                }
            }

            double span = right_max[i] - right_min[i];
            double bin_interval = span / (bin_count);
            vector<double> bin_scores;

            operator_split.n.min = XYZ();
            operator_split.n.max = XYZ();
            operator_split.n.count = 0;
            operator_split.p.min = right_min;
            operator_split.p.max = right_max;
            operator_split.p.count = geo_count;
            evaluate_split(operator_split, 1);
            if (operator_split.score < best.score) best = operator_split;

            double pos = left_min[i];
            int index = 0;
            set<Tri*> partials;
            for (int j = 0; j < bin_count; j++) {
                while (index < sorted.size() && sorted[index].get_value()[i] < pos) {
                    relevant_value& v = sorted[index];
                    assert(v.parent != nullptr);
                    if (v.value_selector == 0) {
                        partials.insert(v.parent);
                        left_count++;
                    }
                    if (v.value_selector == 2) {
                        partials.erase(v.parent);
                        right_min = v.right_min;
                        right_max = v.right_max;
                        left_min = XYZ::min(v.min, left_min);
                        left_max = XYZ::max(v.max, left_max);
                        right_count--;
                    }
                    index++;
                }

                XYZ temp_right_min = right_min;
                XYZ temp_right_max = right_max;
                XYZ temp_left_min = left_min;
                XYZ temp_left_max = left_max;

                vector<XYZ> planar;
                for (Tri* T : partials) {
                    vector<XYZ> points;
                    vector<XYZ> left; //note, this is probably pretty inefficient, but its a simple solution
                    vector<XYZ> right;
                    points.push_back(T->p1);
                    points.push_back(T->p2);
                    points.push_back(T->p3);
                    for (int ti = 0; ti < 3; ti++) {
                        XYZ p = points[ti];
                        if (p[i] < pos)left.push_back(p);
                        if (p[i] == pos)planar.push_back(p);
                        if (p[i] > pos)right.push_back(p);
                    }

                    for (int li = 0; li < left.size(); li++) {
                        for (int ri = 0; ri < right.size(); ri++) {
                            XYZ v1 = left[li];
                            XYZ v2 = right[ri];
                            XYZ slope = XYZ::slope(v1, v2);
                            double t_span = v2[i] - v1[i];
                            double t = (pos - v1[i]) / t_span;
                            XYZ planar_point = slope * t + v1;
                            planar.push_back(planar_point);
                        }
                    }
                }
                for (XYZ& p : planar) {
                    temp_right_min = XYZ::min(p, temp_right_min);
                    temp_right_max = XYZ::max(p, temp_right_max);
                    temp_left_min = XYZ::min(p, temp_left_min);
                    temp_left_max = XYZ::max(p, temp_left_max);
                }
                operator_split.n.min = temp_left_min;
                operator_split.n.max = temp_left_max;
                operator_split.n.count = left_count;
                operator_split.p.min = temp_right_min;
                operator_split.p.max = temp_right_max;
                operator_split.p.count = right_count;
                operator_split.placement = temp_left_max;
                evaluate_split(operator_split, 1);
                if (operator_split.score < 0) {
                    assert(0);
                }
                if (operator_split.score < best.score) {
                    best = operator_split;
                }
                bin_scores.push_back(operator_split.score);
                pos += bin_interval;

            }
            for (int i = 0; i < bin_scores.size(); i++) {
                cout << bin_scores[i] << endl;
            }
        }
        exit(0);

        BinResults out_bin = BinResults();
        out_bin.n_geo = new vector<Tri*>();
        out_bin.p_geo = new vector<Tri*>();
        int facing = best.facing;
        for (Tri* T : *geo) {
            if (T->AABB_max[facing] < best.placement[facing]) {
                out_bin.n_geo->push_back(T);
                continue;
            }
            if (T->AABB_min[facing] > best.placement[facing]) {
                out_bin.p_geo->push_back(T);
                continue;
            }

            vector<XYZ> points;
            vector<XYZ> left; //note, this is probably pretty inefficient, but its a simple solution
            vector<XYZ> right;
            points.push_back(T->p1);
            points.push_back(T->p2);
            points.push_back(T->p3);
            Tri* left_T = new Tri(*T);
            Tri* right_T = new Tri(*T);
            left_T->AABB_min = T->AABB_max;
            left_T->AABB_max = T->AABB_min;
            right_T->AABB_min = T->AABB_max;
            right_T->AABB_max = T->AABB_min;
            double pos = best.placement[facing];
            for (int ti = 0; ti < 3; ti++) {
                XYZ p = points[ti];
                if (p[facing] < pos) {
                    left.push_back(p);
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                }
                if (p[facing] == pos) {
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
                if (p[facing] > pos) {
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
            }

            for (int li = 0; li < left.size(); li++) {
                for (int ri = 0; ri < right.size(); ri++) {
                    XYZ v1 = left[li];
                    XYZ v2 = right[ri];
                    XYZ slope = XYZ::slope(v1, v2);
                    double t_span = v2[facing] - v1[facing];
                    double t = (pos - v1[facing]) / t_span;
                    XYZ p = slope * t + v1;
                    left_T->AABB_min = XYZ::min(left_T->AABB_min, p);
                    left_T->AABB_max = XYZ::max(left_T->AABB_max, p);
                    right_T->AABB_min = XYZ::min(right_T->AABB_min, p);
                    right_T->AABB_max = XYZ::max(right_T->AABB_max, p);
                }
            }
            out_bin.n_geo->push_back(left_T);
            out_bin.p_geo->push_back(right_T);
        }
        return pair<Split*, BinResults>(new Split(best), out_bin);
    }
    Split* probe(vector<Tri*>* geo) {

        auto sorted = vector<relevant_value>();
        for (Tri* T_ptr : *geo) {
            sorted.push_back(relevant_value(T_ptr));
        }
        Split operator_split = Split(XYZ(), 0);
        Split best = operator_split;
        best.score = 999999999999999999999999.0;
        int geo_count = geo->size();
        int iteration = 0;
        double last_read = 0;
        for (int i = 0; i < 3; i++) {
            operator_split.facing = i;
            sort(sorted.begin(), sorted.end(), relevant_value::comparer(i));

            XYZ left_min = sorted.front().min;
            XYZ left_max = sorted.front().max;
            XYZ right_min = sorted.back().min;
            XYZ right_max = sorted.back().max;
            int count = 0;
            for (int j = sorted.size() - 1; j >= 0; j--) {
                relevant_value& v = sorted[j];
                right_min = XYZ::min(right_min, v.min);
                right_max = XYZ::max(right_max, v.max);
                v.right_min = right_min;
                v.right_max = right_max;
            }

            operator_split.n.min = XYZ();
            operator_split.n.max = XYZ();
            operator_split.n.count = 0;
            operator_split.p.min = right_min;
            operator_split.p.max = right_max;
            operator_split.p.count = geo_count;
            evaluate_split(operator_split, 1);
            if (operator_split.score < best.score) best = operator_split;

            double prev_pos = -99999999;
            int segments = 100;
            double span = right_max[i] - right_min[i];
            double segment = segments / span;


            for (int j = 0; j < sorted.size(); j++) {
                auto v = sorted[j];
                left_min = XYZ::min(left_min, v.min);
                left_max = XYZ::max(left_max, v.max);
                count++;
                prev_pos = v.midpoint[i];
                operator_split.n.min = left_min;
                operator_split.n.max = left_max;
                operator_split.n.count = count;
                operator_split.p.min = v.right_min;
                operator_split.p.max = v.right_max;
                operator_split.p.count = geo_count - count;
                operator_split.placement = v.midpoint;
                evaluate_split(operator_split, 1);
                if (operator_split.score < best.score) {
                    best = operator_split;
                }
                if (abs(operator_split.score - last_read) / operator_split.score > 0.05) {
                    //cout << iteration << endl;
                    //cout << operator_split.score << endl;
                    last_read = operator_split.score;
                }
                iteration++;
            }
        }
        //exit(0);

        return new Split(best);
    }
    static void operate_stats(Stats& stat, const XYZ& point_max, const XYZ& point_min) {
        stat.count++;
        stat.max = XYZ::max(point_max, stat.max);
        stat.min = XYZ::min(point_min, stat.min);
    }
    static void get_stats(vector<Tri*>* geo, Split& split) {
        split.p.count = 0;
        split.n.count = 0;
        split.p.max = XYZ(-99999999);
        split.n.max = XYZ(-99999999);
        split.p.min = XYZ(9999999);
        split.n.min = XYZ(9999999);
        int axis = split.facing;
        for (const Tri* T_ptr : *geo) {
            if (T_ptr->midpoint[axis] > split.placement[axis]) {
                operate_stats(split.p, T_ptr->AABB_max, T_ptr->AABB_min);
            }
            if (T_ptr->midpoint[axis] <= split.placement[axis]) {
                operate_stats(split.n, T_ptr->AABB_max, T_ptr->AABB_min);
            }
        }
    }
    static float evaluate_split(Split& split, float SA_parent) {
        float traversal_cost = 1;
        float intersect_cost = 32;
        float Sa = split.p.SA() / SA_parent;
        float Sb = split.n.SA() / SA_parent;
        int count_a = split.p.count;
        int count_b = split.n.count;
#if PENALIZE_UNFILLED_LEAFS
        count_a = ceil(count_a / (float)LEAF_SIZE);
        count_b = ceil(count_b / (float)LEAF_SIZE);
        intersect_cost = 2;
#endif
        float Ha = Sa * count_a * intersect_cost;
        float Hb = Sb * count_b * intersect_cost;
        float final_h = traversal_cost + Ha + Hb;
        split.score = final_h;
        return final_h;
        //return (split.p.count+split.n.count)/geo->size()*(split.p.count+split.n.count);
    }

};



class OctBVH {
public:
    OctBVH* children[8];
    XYZ max;
    XYZ min;
    vector<Tri*> elements;
    OctBVH* parent;

    struct sTri {
        XYZ max;
        XYZ min;
        Tri* parent;
        char x = -1;
        char y = -1;
        char z = -1;
        sTri(Tri* t) {
            max = t->AABB_max;
            min = t->AABB_min;
            parent = t;
        }
    };

    struct Packet {
        XYZ max;
        XYZ min;
        sTri* tris;
        int size = 0;
        void allocate(int num) {
            tris = (sTri*) _aligned_malloc(sizeof(sTri) * num, 64);
            size = num;
        }
    };

    void construct(vector<Tri>* tris) {
        Packet* entry_packet = new Packet();
        entry_packet->allocate(tris->size());
        for (int i = 0; i < tris->size(); i++) {
            entry_packet->tris[i] = sTri(&(*tris)[i]);
        }
        recurse_construct(entry_packet);
    }
private:
    void recurse_construct(Packet* pack) {
        int subdivs = 4;
        int total_segs = pow(2, subdivs);
        XYZ offset = pack->min;
        XYZ subdiv_mult = ((pack->max - pack->min) / total_segs);
        XYZ subdiv_mult_inv = 1/subdiv_mult;
        
        vector<vector<vector<int>>> segments;

        segments.resize(total_segs);
        for (auto& vec_y : segments) {
            vec_y.resize(total_segs);
            for (auto& vec_z : vec_y) {
                vec_z.resize(total_segs);
                for (auto& v : vec_z) {
                    v = 0;
                }
            }
        }

        vector<sTri> extras;
        for (int i = 0; i < pack->size; i++) {
            sTri& t = pack->tris[i];
            XYZ Mn_loc = XYZ::floor((t.min - offset) * subdiv_mult_inv);
            XYZ Mx_loc = XYZ::floor((t.max - offset) * subdiv_mult_inv);
            if (!XYZ::equals(Mn_loc, Mx_loc)) { //note, this may heavily duplicates tris
                vector<sTri*> current_segments;
                current_segments.push_back(&t);
                XYZ split_loc = Mx_loc*subdiv_mult+offset; 
                if (Mn_loc.X != Mx_loc.X) {
                    for (sTri* sT : current_segments) {
                        extras.push_back(sTri(*sT));
                        current_segments.push_back(&extras.back());
                        sTri* sT2 = &extras.back();
                        sT2->min.X = split_loc.X;
                        sT->max.X = split_loc.X;
                    }
                }
                if (Mn_loc.Y != Mx_loc.Y) {
                    for (sTri* sT : current_segments) {
                        extras.push_back(sTri(*sT));
                        current_segments.push_back(&extras.back());
                        sTri* sT2 = &extras.back();
                        sT2->min.Y = split_loc.Y;
                        sT->max.Y = split_loc.Y;
                    }
                }
                if (Mn_loc.Z != Mx_loc.Z) {
                    for (sTri* sT : current_segments) {
                        extras.push_back(sTri(*sT));
                        current_segments.push_back(&extras.back());
                        sTri* sT2 = &extras.back();
                        sT2->min.Y = split_loc.Y;
                        sT->max.Y = split_loc.Y;
                    }
                }

            }
            t.x = Mn_loc.X;
            t.y = Mn_loc.Y;
            t.z = Mn_loc.Z;
        }


    }
};

class PackagedBVH { //memory optimized. produced once the BVH tree is finalized
public:
    int32_t index;
    XYZ sMax;
    XYZ sMin;
    int32_t leaf_size;
    PackagedBVH(BVH* target) {
        sMax = target->max;
        sMin = target->min;
        index = -1;
    }
    static pair<PackagedBVH*, vector<PackagedTri>*> collapse(BVH* top) {
        vector<PackagedTri>* tri_vec = new vector<PackagedTri>();
        int tri_count = top->count();
        int tree_size = top->size();
        PackagedBVH* out = (PackagedBVH*)_aligned_malloc((tree_size + 1) * sizeof(PackagedBVH), 32);
        tri_vec->reserve(tri_count + 1);
        PackagedBVH::collapse(out, tri_vec, top);
        return pair<PackagedBVH*, vector<PackagedTri>*>(out, tri_vec);
    }
private:
    static void collapse(PackagedBVH* vec, vector<PackagedTri>* tri_vec, BVH* top) {
        vector<pair<BVH*, int>> current;
        vector<pair<BVH*, int>> next;

        int accumulated_index = 0;
        current.push_back(pair<BVH*, int>(top, -1));
        while (current.size() > 0) {
            for (int i = 0; i < current.size(); i++) {
                auto selection = current.at(i);
                BVH* target = selection.first;
                int parent_index = selection.second;
                auto PBVH = PackagedBVH(target);
                XYZ a = 0.5 * (PBVH.sMax + PBVH.sMin);
                //PBVH.sMax += 0.5*(PBVH.sMax-a);
                //PBVH.sMin += 0.5*(PBVH.sMin-a);
                PBVH.leaf_size = target->elements.size();
                if (PBVH.leaf_size > 0) {
                    PBVH.index = tri_vec->size();
                    for (int i = 0; i < PBVH.leaf_size; i++) {
                        tri_vec->push_back(target->elements[i]->pack());
                    }
                }
                if (target->c1 != nullptr) {
                    next.push_back(pair<BVH*, int>(target->c1, accumulated_index));
                    next.push_back(pair<BVH*, int>(target->c2, accumulated_index));
                }
                vec[accumulated_index] = PBVH;
                if (parent_index != -1) {
                    if (vec[parent_index].index == -1) {
                        vec[parent_index].index = accumulated_index;
                    }

                }
                else accumulated_index++;
                accumulated_index++;
            }
            current = next;
            next.clear();
        }
    }
};



class BVH_AVX {
public:
    m256_vec3 max;
    m256_vec3 min;
    unsigned int leaf_size[8];
    unsigned int indexes[8];
    int parent_index;
    int self_index;
    BVH_AVX(vector<BVH*>& batch) {
        vector<XYZ> max_vec;
        vector<XYZ> min_vec;
        for (int i = 0; i < batch.size(); i++) {
            max_vec.push_back(batch[i]->max);
            min_vec.push_back(batch[i]->min);
            leaf_size[i] = batch[i]->elements.size();
        }
        for (int i = batch.size(); i < 8;i++) {
            max_vec.push_back(XYZ());
            min_vec.push_back(XYZ());
            leaf_size[i] = 0;
            indexes[i] = 0;
        }
        max = m256_vec3(max_vec);
        min = m256_vec3(min_vec);
    }
    __m256 intersectionv1(const m256_vec3& fusedorigin, const m256_vec3& inv_slope) const { //holy mother of moving data. I pray for you, my cpu, I pray
        __m256 tx1 = _mm256_fmadd_ps(min.X, inv_slope.X, fusedorigin.X);
        __m256 tx2 = _mm256_fmadd_ps(max.X, inv_slope.X, fusedorigin.X);
        __m256 ty1 = _mm256_fmadd_ps(min.Y, inv_slope.Y, fusedorigin.Y);
        __m256 ty2 = _mm256_fmadd_ps(max.Y, inv_slope.Y, fusedorigin.Y);
        __m256 tz1 = _mm256_fmadd_ps(min.Z, inv_slope.Z, fusedorigin.Z);
        __m256 tz2 = _mm256_fmadd_ps(max.Z, inv_slope.Z, fusedorigin.Z);

        __m256 tmin = _mm256_max_ps(
            _mm256_min_ps(tz1, tz2),
            _mm256_max_ps(
                _mm256_min_ps(ty1, ty2),
                _mm256_min_ps(tx1, tx2)
            ));
        __m256 tmax = _mm256_min_ps(
            _mm256_max_ps(tz1, tz2),
            _mm256_min_ps(
                _mm256_max_ps(ty1, ty2),
                _mm256_max_ps(tx1, tx2)
            ));
        __m256 diff = _mm256_sub_ps(tmax, tmin);
        __m256 masked_diff = _mm256_max_ps(
            diff,
            AVX_ZEROS
        );
        __m256 return_mask = _mm256_cmp_ps(
            masked_diff,
            diff,
            _CMP_EQ_OQ
        );
        __m256 final = _mm256_and_ps(
            _mm256_add_ps(
                tmin,
                masked_diff
            ),
            return_mask
        );
        //__m256 diff_is_pos = _mm256_cmp_ps(diff, AVX_ZEROS, _CMP_GT_OQ);

        return final;
    }
    __m256 intersection(const m256_vec3& fusedorigin, const m256_vec3& inv_slope) const { //holy mother of moving data. I pray for you, my cpu, I pray
        __m256 tx1 = _mm256_fmadd_ps(min.X, inv_slope.X, fusedorigin.X);
        __m256 tx2 = _mm256_fmadd_ps(max.X, inv_slope.X, fusedorigin.X);
        __m256 ty1 = _mm256_fmadd_ps(min.Y, inv_slope.Y, fusedorigin.Y);
        __m256 ty2 = _mm256_fmadd_ps(max.Y, inv_slope.Y, fusedorigin.Y);
        __m256 tz1 = _mm256_fmadd_ps(min.Z, inv_slope.Z, fusedorigin.Z);
        __m256 tz2 = _mm256_fmadd_ps(max.Z, inv_slope.Z, fusedorigin.Z);

        __m256 tminx = _mm256_min_ps(tx1, tx2);
        __m256 tminy = _mm256_min_ps(ty1, ty2);
        __m256 tminz = _mm256_min_ps(tz1, tz2);
        __m256 tmaxx = _mm256_max_ps(tx1, tx2);
        __m256 tmaxy = _mm256_max_ps(ty1, ty2);
        __m256 tmaxz = _mm256_max_ps(tz1, tz2);

        __m256 tmin = _mm256_max_ps(
            tminz,
            _mm256_max_ps(
                tminy,
                tminx
            ));
        __m256 tmax = _mm256_min_ps(
            tmaxz,
            _mm256_min_ps(
                tmaxy,
                tmaxx
            ));
        __m256 return_mask = _mm256_and_ps(
            _mm256_cmp_ps(
                tmax,
                AVX_ZEROS,
                _CMP_GE_OQ
            ),
            _mm256_cmp_ps(
                tmin,
                tmax,
                _CMP_LT_OQ
            )
        );


        __m256 final = _mm256_and_ps(
            tmin,
            return_mask
        );
        //__m256 diff_is_pos = _mm256_cmp_ps(diff, AVX_ZEROS, _CMP_GT_OQ);

        return final;
    }
    static pair<vector<BVH_AVX>*, vector<PackagedTri>*> collapse(BVH* top) {
        vector<PackagedTri>* tri_vec = new vector<PackagedTri>();
        int tri_count = top->count();
        int tree_size = top->size();
        vector<BVH_AVX> v;
        v.reserve(top->size());
        tri_vec->reserve(tri_count + 1);
        BVH_AVX::collapse(v, tri_vec, top);
        vector<BVH_AVX>* out = new vector<BVH_AVX>();
        out->swap(v);//trimming overallocation of vector
        return pair<vector<BVH_AVX>*, vector<PackagedTri>*>(out, tri_vec);
    }
private:
    struct comparer {
        bool operator()(BVH* b1, BVH* b2) {
            //return b1->avg_path() < b2->avg_path();
            return VecLib::surface_area(b1->max, b1->min) * b1->count() > VecLib::surface_area(b2->max, b2->min) * b2->count();
        }
    };
    static void collapse(vector<BVH_AVX>& vec, vector<PackagedTri>* tri_vec, BVH* top) {
        vector<pair<BVH*, pair<int, int>>> current;
        vector<pair<BVH*, pair<int, int>>> next;

        int accumulated_index = 0;
        current.push_back(pair<BVH*, pair<int, int>>(top, pair<int, int>(-1, -1)));
        while (current.size() > 0) {
            for (int i = 0; i < current.size(); i++) {
                auto selection = current.at(i);
                vector<BVH*> batch;
                vector<BVH*> has_tris;
                batch.push_back(selection.first);
                while (true) {
                    while (batch.size() + has_tris.size() < 8) {
                        BVH* at_index = batch[0];
                        if (at_index->c1 != nullptr) {
                            batch.push_back(at_index->c1);
                            batch.push_back(at_index->c2);
                            batch.erase(batch.begin());
                        }
                        else {
                            has_tris.push_back(at_index);
                            batch.erase(batch.begin());
                            if (batch.size() == 0) {
                                break;
                            }
                        }
                        sort(batch.begin(), batch.end(), comparer());
                    }
                    if (has_tris.size() + batch.size() >= 8 || batch.size() == 0) {
                        if (has_tris.size() + batch.size() > 8) {
                            cout << "error occured in AVX BVH construction: over allocated" << endl;
                            throw exception();
                        }
                        for (int i = 0; i < has_tris.size(); i++) {
                            batch.push_back(has_tris[i]);
                        }
                        break;
                    }
                }
                BVH* target = selection.first;
                int parent_index = selection.second.first;
                int point_back_index = selection.second.second;
                BVH_AVX ABVH(batch);
                ABVH.parent_index = parent_index;
                ABVH.self_index = point_back_index;
                for (int i = 0; i < batch.size(); i++) {
                    if (ABVH.leaf_size[i] > 0) {
                        ABVH.indexes[i] = tri_vec->size();
                        for (int j = 0; j < ABVH.leaf_size[i]; j++) {
                            tri_vec->push_back(batch[i]->elements[j]->pack());
                        }
                    }
                    else {
                        next.push_back(pair<BVH*, pair<int, int>>(batch[i], pair<int, int>(accumulated_index, i)));
                    }
                }
                vec.push_back(ABVH);
                if (parent_index != -1) {
                    BVH_AVX& parent = vec[parent_index];
                    for (int k = 0; k < 8; k++) {
                        parent.indexes[point_back_index] = accumulated_index;
                    }
                }
                accumulated_index++;
            }
            current = next;
            next.clear();
        }
    }
};
