#include "pybmad/generated/structs_v.hpp"

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
// var_length_string_struct
void init_var_length_string_struct(nb::module_ &m, nb::class_<VarLengthStringStruct> &cls) {
  cls.def(nb::init<std::optional<std::string>>(), nb::arg("str") = nb::none())
      .def_prop_rw("str", &VarLengthStringStruct::str, &VarLengthStringStruct::set_str)
      .def_static(
          "new_array1d",
          [](int sz) { return VarLengthStringStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = VarLengthStringStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const VarLengthStringStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const VarLengthStringStruct &self) {
            return VarLengthStringStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const VarLengthStringStruct &self, nb::dict &memo) {
            return VarLengthStringStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const VarLengthStringStruct &self, const VarLengthStringStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const VarLengthStringStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<VarLengthStringStructArray1D, VarLengthStringStructAlloc1D>(
      m,
      "VarLengthStringStructArray1D",
      "VarLengthStringStructAlloc1D"
  );
  // 2D VarLengthStringStruct arrays are not used in structs/routines
  // 3D VarLengthStringStruct arrays are not used in structs/routines
}