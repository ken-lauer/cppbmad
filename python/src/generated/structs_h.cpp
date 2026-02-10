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
  cls.def(
         py::init<
             optional_ref<const CoordStruct>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("closed_orb") = py::none(),
         py::arg("kick_const") = py::none(),
         py::arg("sig_x") = py::none(),
         py::arg("sig_y") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("sin_phi") = py::none(),
         py::arg("cos_phi") = py::none(),
         py::arg("sig_z") = py::none()
  )
      .def_property(
          "closed_orb",
          py::cpp_function(&HighEnergySpaceChargeStruct::closed_orb, py::keep_alive<0, 1>()),
          &HighEnergySpaceChargeStruct::set_closed_orb,
          "beam orbit"
      )
      .def_property(
          "kick_const",
          &HighEnergySpaceChargeStruct::kick_const,
          &HighEnergySpaceChargeStruct::set_kick_const
      )
      .def_property(
          "sig_x",
          &HighEnergySpaceChargeStruct::sig_x,
          &HighEnergySpaceChargeStruct::set_sig_x
      )
      .def_property(
          "sig_y",
          &HighEnergySpaceChargeStruct::sig_y,
          &HighEnergySpaceChargeStruct::set_sig_y
      )
      .def_property(
          "phi",
          &HighEnergySpaceChargeStruct::phi,
          &HighEnergySpaceChargeStruct::set_phi,
          "Rotation angle to go from lab frame to rotated frame."
      )
      .def_property(
          "sin_phi",
          &HighEnergySpaceChargeStruct::sin_phi,
          &HighEnergySpaceChargeStruct::set_sin_phi
      )
      .def_property(
          "cos_phi",
          &HighEnergySpaceChargeStruct::cos_phi,
          &HighEnergySpaceChargeStruct::set_cos_phi
      )
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