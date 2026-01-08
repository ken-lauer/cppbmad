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

using namespace Bmad;
using namespace Pybmad;

using json = nlohmann::json;

namespace py = pybind11;
using namespace pybind11::literals;

template <typename T>
std::string instance_to_json(const T& self) {
  json j;
  to_json(j, self);
  return to_string(j);
}

template <typename T>
std::optional<std::reference_wrapper<T>> make_opt_ref(
    std::optional<T>& opt_val) {
  if (opt_val.has_value()) {
    return std::ref(opt_val.value());
  }
  return std::nullopt;
}

// ${pybind11_routine_wrappers}

PYBIND11_MODULE(_pybmad, m) {
  // Generated definitions:
  // ${pybind11_definitions}
}
