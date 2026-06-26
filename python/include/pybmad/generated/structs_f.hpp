#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_floor_position_struct(nb::module_ &m, nb::class_<FloorPositionStruct> &class_);
void init_fibre(nb::module_ &m, nb::class_<Fibre> &class_);
