#include "pybmad/generated/structs_i.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

using namespace Pybmad;
namespace nb = nanobind;

// =============================================================================
// interval1_coef_struct
void init_interval1_coef_struct(nb::module_ &m, nb::class_<Interval1CoefStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("c0") = nb::none(),
         nb::arg("c1") = nb::none(),
         nb::arg("n_exp") = nb::none()
  )
      .def_prop_rw("c0", &Interval1CoefStruct::c0, &Interval1CoefStruct::set_c0)
      .def_prop_rw("c1", &Interval1CoefStruct::c1, &Interval1CoefStruct::set_c1)
      .def_prop_rw("n_exp", &Interval1CoefStruct::n_exp, &Interval1CoefStruct::set_n_exp)
      .def_static(
          "new_array1d",
          [](int sz) { return Interval1CoefStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Interval1CoefStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const Interval1CoefStruct &self, nb::dict &memo) { return Interval1CoefStruct(self); }
      )
      .def(
          "__eq__",
          [](const Interval1CoefStruct &self, const Interval1CoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Interval1CoefStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<Interval1CoefStructArray1D, Interval1CoefStructAlloc1D>(
      m,
      "Interval1CoefStructArray1D",
      "Interval1CoefStructAlloc1D"
  );
  // 2D Interval1CoefStruct arrays are not used in structs/routines
  // 3D Interval1CoefStruct arrays are not used in structs/routines
}