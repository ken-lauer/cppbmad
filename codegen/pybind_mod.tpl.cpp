#include <pybind11/complex.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <string>

#include "bmad/json.hpp"
#include "pybmad/common_structs.hpp"
#include "pybmad/generated/init.hpp"
#include "pybmad/util.hpp"

using namespace Bmad;
using namespace Pybmad;

namespace py = pybind11;
using namespace pybind11::literals;

// ${pybind11_routine_wrappers}

PYBIND11_MODULE(_pybmad, m, py::mod_gil_not_used()) {
  // Generated definitions:
  // ${pybind11_definitions}
}
