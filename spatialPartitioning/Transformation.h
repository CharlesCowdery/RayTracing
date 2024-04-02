#pragma once
#include "XYZ.h"

class Transformation {
private:
    Quat* rot_transform = nullptr;
    XYZ* XYZ_transform = nullptr;
public:

    Transformation(XYZ transformation, Quat rotation) {
        rot_transform = new Quat(rotation);
        XYZ_transform = new XYZ(transformation);
    }
    Transformation(Quat rotation) {
        rot_transform = new Quat(rotation);
    }
    Transformation(XYZ transformation) {
        XYZ_transform = new XYZ(transformation);
    }
    Transformation rotation() {
        return Transformation(*rot_transform);
    }
    void stack(XYZ& XYZ_transf, Quat& rotation) {
        if (XYZ_transform != nullptr) {
            XYZ_transf += *XYZ_transform;
        }
        if (rot_transform != nullptr) {
            rotation = Quat::multiply(rotation, *rot_transform);
            XYZ_transf = Quat::applyRotation(XYZ_transf, *rot_transform);
        }
    }
    //takes an XYZ position reference, and applies itself to it
    void apply(XYZ& position) {
        if (rot_transform != nullptr) {
            position = Quat::applyRotation(position, *rot_transform);
        }
        if (XYZ_transform != nullptr) {
            position += *XYZ_transform;
        }
    }
    static Transformation collpase(vector<Transformation>& transforms) {
        XYZ position = XYZ(0, 0, 0);
        Quat rotation = Quat(0, 0, 0, 1);
        for (Transformation& t : transforms) {
            t.stack(position, rotation);
        }
        return Transformation(position, rotation);
    }
};