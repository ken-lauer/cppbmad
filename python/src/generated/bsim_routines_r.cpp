#include "pybmad/generated/bsim_routines_r.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_r(py::module &m) {
  m.def(
      "rf_cav_names",
      &bsim::rf_cav_names,
      py::arg("lat"),
      R"""(Wrapper for Fortran routine rf_cav_names

Parameters
----------
lat : 
)"""
  );
}
