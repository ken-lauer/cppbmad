#pragma once

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <functional>

#include "bmad/routines.hpp"
#include "pybmad/arrays.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

// ${forward_declarations}
