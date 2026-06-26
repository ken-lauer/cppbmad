#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_out_io_output_direct_struct(nb::module_ &m, nb::class_<OutIoOutputDirectStruct> &class_);
