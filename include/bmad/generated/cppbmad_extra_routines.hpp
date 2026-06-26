#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace CppBmadExtra {
extern "C" bool fortran_set_ele_misalignments(
    void *ele /* 0D_NOT_type inout */,
    double &x_offset /* 0D_NOT_real in */,
    double &y_offset /* 0D_NOT_real in */,
    double &z_offset /* 0D_NOT_real in */,
    double &x_pitch /* 0D_NOT_real in */,
    double &y_pitch /* 0D_NOT_real in */,
    double &tilt /* 0D_NOT_real in */,
    bool *check_free /* 0D_NOT_logical in */,
    bool &ok /* 0D_NOT_logical out */
);
bool set_ele_misalignments(
    EleStruct &ele,
    double x_offset,
    double y_offset,
    double z_offset,
    double x_pitch,
    double y_pitch,
    double tilt,
    std::optional<bool> check_free = std::nullopt
);
} // namespace CppBmadExtra
