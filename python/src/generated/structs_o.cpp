#include "pybmad/generated/structs_o.hpp"

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
// out_io_output_direct_struct
void init_out_io_output_direct_struct(nb::module_ &m, nb::class_<OutIoOutputDirectStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::vector<bool>>, std::optional<std::vector<int>>>(),
         nb::arg("print_and_capture") = nb::none(),
         nb::arg("file_unit") = nb::none()
  )
      .def_prop_rw(
          "print_and_capture",
          &OutIoOutputDirectStruct::print_and_capture,
          &OutIoOutputDirectStruct::set_print_and_capture,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "file_unit",
          &OutIoOutputDirectStruct::file_unit,
          &OutIoOutputDirectStruct::set_file_unit,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const OutIoOutputDirectStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const OutIoOutputDirectStruct &self) {
            return OutIoOutputDirectStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const OutIoOutputDirectStruct &self, nb::dict &memo) {
            return OutIoOutputDirectStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const OutIoOutputDirectStruct &self, const OutIoOutputDirectStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const OutIoOutputDirectStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D OutIoOutputDirectStruct arrays are not used in structs/routines
  // 2D OutIoOutputDirectStruct arrays are not used in structs/routines
  // 3D OutIoOutputDirectStruct arrays are not used in structs/routines
}