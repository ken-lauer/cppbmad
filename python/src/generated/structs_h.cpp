#include "pybmad/generated/structs_h.hpp"

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
// high_energy_space_charge_struct
void init_high_energy_space_charge_struct(
    nb::module_ &m,
    nb::class_<HighEnergySpaceChargeStruct> &cls
) {
  cls.def(
         "__init__",
         [](HighEnergySpaceChargeStruct *self,
            const CoordStruct *closed_orb,
            std::optional<double> kick_const,
            std::optional<double> sig_x,
            std::optional<double> sig_y,
            std::optional<double> phi,
            std::optional<double> sin_phi,
            std::optional<double> cos_phi,
            std::optional<double> sig_z) {
           new (self) HighEnergySpaceChargeStruct(
               ptr_to_opt_ref(closed_orb),
               kick_const,
               sig_x,
               sig_y,
               phi,
               sin_phi,
               cos_phi,
               sig_z
           );
         },
         nb::arg("closed_orb") = nb::none(),
         nb::arg("kick_const") = nb::none(),
         nb::arg("sig_x") = nb::none(),
         nb::arg("sig_y") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("sin_phi") = nb::none(),
         nb::arg("cos_phi") = nb::none(),
         nb::arg("sig_z") = nb::none()
  )
      .def_prop_rw(
          "closed_orb",
          &HighEnergySpaceChargeStruct::closed_orb,
          &HighEnergySpaceChargeStruct::set_closed_orb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "beam orbit"
      )
      .def_prop_rw(
          "kick_const",
          &HighEnergySpaceChargeStruct::kick_const,
          &HighEnergySpaceChargeStruct::set_kick_const
      )
      .def_prop_rw(
          "sig_x",
          &HighEnergySpaceChargeStruct::sig_x,
          &HighEnergySpaceChargeStruct::set_sig_x
      )
      .def_prop_rw(
          "sig_y",
          &HighEnergySpaceChargeStruct::sig_y,
          &HighEnergySpaceChargeStruct::set_sig_y
      )
      .def_prop_rw(
          "phi",
          &HighEnergySpaceChargeStruct::phi,
          &HighEnergySpaceChargeStruct::set_phi,
          "Rotation angle to go from lab frame to rotated frame."
      )
      .def_prop_rw(
          "sin_phi",
          &HighEnergySpaceChargeStruct::sin_phi,
          &HighEnergySpaceChargeStruct::set_sin_phi
      )
      .def_prop_rw(
          "cos_phi",
          &HighEnergySpaceChargeStruct::cos_phi,
          &HighEnergySpaceChargeStruct::set_cos_phi
      )
      .def_prop_rw(
          "sig_z",
          &HighEnergySpaceChargeStruct::sig_z,
          &HighEnergySpaceChargeStruct::set_sig_z
      )

      .def("__repr__", [](const HighEnergySpaceChargeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const HighEnergySpaceChargeStruct &self) {
            return HighEnergySpaceChargeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const HighEnergySpaceChargeStruct &self, nb::dict &memo) {
            return HighEnergySpaceChargeStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const HighEnergySpaceChargeStruct &self, const HighEnergySpaceChargeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const HighEnergySpaceChargeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D HighEnergySpaceChargeStruct arrays are not used in structs/routines
  // 2D HighEnergySpaceChargeStruct arrays are not used in structs/routines
  // 3D HighEnergySpaceChargeStruct arrays are not used in structs/routines
}