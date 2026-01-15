#include "pybmad/generated/structs_h.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// high_energy_space_charge_struct
void init_high_energy_space_charge_struct(
    py::module &m,
    py::class_<HighEnergySpaceChargeStruct> &cls
) {
  cls.def(py::init<>())
      // HighEnergySpaceChargeStruct.closed_orb (0D_NOT_type - beam orbit
      .def_property(
          "closed_orb",
          &HighEnergySpaceChargeStruct::closed_orb,
          &HighEnergySpaceChargeStruct::set_closed_orb
      )
      // HighEnergySpaceChargeStruct.kick_const (0D_NOT_real -
      .def_property(
          "kick_const",
          &HighEnergySpaceChargeStruct::kick_const,
          &HighEnergySpaceChargeStruct::set_kick_const
      )
      // HighEnergySpaceChargeStruct.sig_x (0D_NOT_real -
      .def_property(
          "sig_x",
          &HighEnergySpaceChargeStruct::sig_x,
          &HighEnergySpaceChargeStruct::set_sig_x
      )
      // HighEnergySpaceChargeStruct.sig_y (0D_NOT_real -
      .def_property(
          "sig_y",
          &HighEnergySpaceChargeStruct::sig_y,
          &HighEnergySpaceChargeStruct::set_sig_y
      )
      // HighEnergySpaceChargeStruct.phi (0D_NOT_real - Rotation angle to go from lab frame to
      // rotated frame.
      .def_property("phi", &HighEnergySpaceChargeStruct::phi, &HighEnergySpaceChargeStruct::set_phi)
      // HighEnergySpaceChargeStruct.sin_phi (0D_NOT_real -
      .def_property(
          "sin_phi",
          &HighEnergySpaceChargeStruct::sin_phi,
          &HighEnergySpaceChargeStruct::set_sin_phi
      )
      // HighEnergySpaceChargeStruct.cos_phi (0D_NOT_real -
      .def_property(
          "cos_phi",
          &HighEnergySpaceChargeStruct::cos_phi,
          &HighEnergySpaceChargeStruct::set_cos_phi
      )
      // HighEnergySpaceChargeStruct.sig_z (0D_NOT_real -
      .def_property(
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
          [](const HighEnergySpaceChargeStruct &self, py::dict &memo) {
            return HighEnergySpaceChargeStruct(self);
          }
      )

      ;

  // 1D HighEnergySpaceChargeStruct arrays are not used in structs/routines
  // 2D HighEnergySpaceChargeStruct arrays are not used in structs/routines
  // 3D HighEnergySpaceChargeStruct arrays are not used in structs/routines
}