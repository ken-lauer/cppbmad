#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_kv_beam_init_struct(nb::module_ &m, nb::class_<KvBeamInitStruct> &class_);
