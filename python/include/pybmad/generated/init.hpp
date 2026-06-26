#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <functional>

#include "bmad/routines.hpp"
#include "pybmad/arrays.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;

namespace Pybmad {} // namespace Pybmad
