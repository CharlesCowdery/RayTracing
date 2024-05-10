#pragma once
#include "XYZ.h"
#include "commons.h"
#include "Coredefs.h"

CastResults::CastResults() :  material(-1), distance(default_distance) {};
CastResults::CastResults(short _mat) : material(_mat), distance(default_distance) {}

const float CastResults::default_distance = 99999999999999999;

