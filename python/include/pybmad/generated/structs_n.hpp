#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_normal_modes_struct(nb::module_ &m, nb::class_<NormalModesStruct> &class_);
void init_nametable_struct(nb::module_ &m, nb::class_<NametableStruct> &class_);
