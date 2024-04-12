#pragma once
#include "XYZ.h"
#include "commons.h"
#include "Coredefs.h"

CastResults::CastResults() : normal(XYZ()), material(nullptr), distance(default_distance) {};
CastResults::CastResults(XYZ _normal, Material* _mat) : normal(_normal), material(_mat), distance(default_distance) {}

const float CastResults::default_distance = 99999999999999999;

