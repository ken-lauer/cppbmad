#include "pybmad/generated/Bmad_routines_y.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_y(py::module& m) {
  m.def(
      "ylafun",
      &Bmad::ylafun,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("res"),
      R"""(Parameters
----------
x : 
y : 
z : 
res : 
)""");
}
