#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_diffuse_param_struct(nb::module_ &m, nb::class_<DiffuseParamStruct> &class_);
void init_do_loop_struct(nb::module_ &m, nb::class_<DoLoopStruct> &class_);
