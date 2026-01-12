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
    py::module& m,
    py::class_<HighEnergySpaceChargeProxy>& cls) {
  cls.def(py::init<>())
      // HighEnergySpaceChargeProxy.closed_orb (0D_NOT_type - beam orbit
      .def_property(
          "closed_orb",
          &HighEnergySpaceChargeProxy::closed_orb,
          &HighEnergySpaceChargeProxy::set_closed_orb)
      // HighEnergySpaceChargeProxy.kick_const (0D_NOT_real -
      .def_property(
          "kick_const",
          &HighEnergySpaceChargeProxy::kick_const,
          &HighEnergySpaceChargeProxy::set_kick_const)
      // HighEnergySpaceChargeProxy.sig_x (0D_NOT_real -
      .def_property(
          "sig_x",
          &HighEnergySpaceChargeProxy::sig_x,
          &HighEnergySpaceChargeProxy::set_sig_x)
      // HighEnergySpaceChargeProxy.sig_y (0D_NOT_real -
      .def_property(
          "sig_y",
          &HighEnergySpaceChargeProxy::sig_y,
          &HighEnergySpaceChargeProxy::set_sig_y)
      // HighEnergySpaceChargeProxy.phi (0D_NOT_real - Rotation angle to go from lab frame to rotated frame.
      .def_property(
          "phi",
          &HighEnergySpaceChargeProxy::phi,
          &HighEnergySpaceChargeProxy::set_phi)
      // HighEnergySpaceChargeProxy.sin_phi (0D_NOT_real -
      .def_property(
          "sin_phi",
          &HighEnergySpaceChargeProxy::sin_phi,
          &HighEnergySpaceChargeProxy::set_sin_phi)
      // HighEnergySpaceChargeProxy.cos_phi (0D_NOT_real -
      .def_property(
          "cos_phi",
          &HighEnergySpaceChargeProxy::cos_phi,
          &HighEnergySpaceChargeProxy::set_cos_phi)
      // HighEnergySpaceChargeProxy.sig_z (0D_NOT_real -
      .def_property(
          "sig_z",
          &HighEnergySpaceChargeProxy::sig_z,
          &HighEnergySpaceChargeProxy::set_sig_z)

      .def(
          "__repr__",
          [](const HighEnergySpaceChargeProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const HighEnergySpaceChargeProxy& self) {
            return HighEnergySpaceChargeProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const HighEnergySpaceChargeProxy& self, py::dict& memo) {
            return HighEnergySpaceChargeProxy(self);
          })

      ;

  // 1D HighEnergySpaceChargeProxy arrays are not used in structs/routines
  // 2D HighEnergySpaceChargeProxy arrays are not used in structs/routines
  // 3D HighEnergySpaceChargeProxy arrays are not used in structs/routines
}