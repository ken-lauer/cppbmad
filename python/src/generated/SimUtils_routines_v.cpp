#include "pybmad/generated/SimUtils_routines_v.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_v(nb::module_ &m) {
  m.def(
      "value_of_all_ptr",
      &SimUtils::value_of_all_ptr,
      nb::arg("a_ptr"),
      R"""(Wrapper for Fortran routine value_of_all_ptr

Parameters
----------
a_ptr : AllPointerStruct

Returns
-------
value : float
)"""
  );
  m.def(
      "virtual_memory_usage",
      &SimUtils::virtual_memory_usage,
      R"""(Wrapper for Fortran routine virtual_memory_usage

Parameters
----------

Returns
-------
usage : int
)"""
  );
}
