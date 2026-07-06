#include "pybmad/generated/structs_x.hpp"

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
// xy_disp_struct
void init_xy_disp_struct(nb::module_ &m, nb::class_<XyDispStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("eta") = nb::none(),
         nb::arg("etap") = nb::none(),
         nb::arg("deta_ds") = nb::none(),
         nb::arg("sigma") = nb::none(),
         nb::arg("deta_dpz") = nb::none(),
         nb::arg("detap_dpz") = nb::none()
  )
      .def_prop_rw("eta", &XyDispStruct::eta, &XyDispStruct::set_eta)
      .def_prop_rw("etap", &XyDispStruct::etap, &XyDispStruct::set_etap)
      .def_prop_rw("deta_ds", &XyDispStruct::deta_ds, &XyDispStruct::set_deta_ds)
      .def_prop_rw("sigma", &XyDispStruct::sigma, &XyDispStruct::set_sigma)
      .def_prop_rw("deta_dpz", &XyDispStruct::deta_dpz, &XyDispStruct::set_deta_dpz)
      .def_prop_rw("detap_dpz", &XyDispStruct::detap_dpz, &XyDispStruct::set_detap_dpz)

      .def("__repr__", [](const XyDispStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const XyDispStruct &self) {
            return XyDispStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const XyDispStruct &self, nb::dict &memo) { return XyDispStruct(self); }
      )
      .def(
          "__eq__",
          [](const XyDispStruct &self, const XyDispStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const XyDispStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D XyDispStruct arrays are not used in structs/routines
  // 2D XyDispStruct arrays are not used in structs/routines
  // 3D XyDispStruct arrays are not used in structs/routines
}