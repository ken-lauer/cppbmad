#include "pybmad/generated/structs_i.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// interval1_coef_struct
void init_interval1_coef_struct(py::module &m, py::class_<Interval1CoefStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         py::arg("c0") = py::none(),
         py::arg("c1") = py::none(),
         py::arg("n_exp") = py::none()
  )
      // Interval1CoefStruct.c0 (0D_NOT_real -
      .def_property("c0", &Interval1CoefStruct::c0, &Interval1CoefStruct::set_c0)
      // Interval1CoefStruct.c1 (0D_NOT_real -
      .def_property("c1", &Interval1CoefStruct::c1, &Interval1CoefStruct::set_c1)
      // Interval1CoefStruct.n_exp (0D_NOT_real -
      .def_property("n_exp", &Interval1CoefStruct::n_exp, &Interval1CoefStruct::set_n_exp)
      .def_static(
          "new_array1d",
          [](int sz) { return Interval1CoefStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Interval1CoefStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const Interval1CoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Interval1CoefStruct &self) {
            return Interval1CoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Interval1CoefStruct &self, py::dict &memo) { return Interval1CoefStruct(self); }
      )

      ;

  bind_FTypeArrayND<Interval1CoefStructArray1D>(m, "Interval1CoefStructArray1D");
  bind_FTypeAlloc1D<Interval1CoefStructAlloc1D>(m, "Interval1CoefStructAlloc1D");
  // 2D Interval1CoefStruct arrays are not used in structs/routines
  // 3D Interval1CoefStruct arrays are not used in structs/routines
}