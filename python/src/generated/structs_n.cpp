#include "pybmad/generated/structs_n.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// normal_modes_struct
void init_normal_modes_struct(
    py::module& m,
    py::class_<NormalModesProxy>& cls) {
  cls.def(py::init<>())
      // NormalModesProxy.synch_int (1D_NOT_real - Synchrotron integrals I0, I1, I2, and I3
      .def_property_readonly("synch_int", &NormalModesProxy::synch_int)
      // NormalModesProxy.sigE_E (0D_NOT_real - SigmaE/E
      .def_property(
          "sigE_E", &NormalModesProxy::sigE_E, &NormalModesProxy::set_sigE_E)
      // NormalModesProxy.sig_z (0D_NOT_real - Sigma_Z
      .def_property(
          "sig_z", &NormalModesProxy::sig_z, &NormalModesProxy::set_sig_z)
      // NormalModesProxy.e_loss (0D_NOT_real - Energy loss / turn (eV)
      .def_property(
          "e_loss", &NormalModesProxy::e_loss, &NormalModesProxy::set_e_loss)
      // NormalModesProxy.rf_voltage (0D_NOT_real - Total rfcavity voltage (eV)
      .def_property(
          "rf_voltage",
          &NormalModesProxy::rf_voltage,
          &NormalModesProxy::set_rf_voltage)
      // NormalModesProxy.pz_aperture (0D_NOT_real - pz aperture limit. Used with Touschek calculations.
      .def_property(
          "pz_aperture",
          &NormalModesProxy::pz_aperture,
          &NormalModesProxy::set_pz_aperture)
      // NormalModesProxy.pz_average (0D_NOT_real - Average over branch due to damping.
      .def_property(
          "pz_average",
          &NormalModesProxy::pz_average,
          &NormalModesProxy::set_pz_average)
      // NormalModesProxy.momentum_compaction (0D_NOT_real -
      .def_property(
          "momentum_compaction",
          &NormalModesProxy::momentum_compaction,
          &NormalModesProxy::set_momentum_compaction)
      // NormalModesProxy.dpz_damp (0D_NOT_real - Change in pz without RF
      .def_property(
          "dpz_damp",
          &NormalModesProxy::dpz_damp,
          &NormalModesProxy::set_dpz_damp)
      // NormalModesProxy.a (0D_NOT_type -
      .def_property("a", &NormalModesProxy::a, &NormalModesProxy::set_a)
      // NormalModesProxy.b (0D_NOT_type -
      .def_property("b", &NormalModesProxy::b, &NormalModesProxy::set_b)
      // NormalModesProxy.z (0D_NOT_type -
      .def_property("z", &NormalModesProxy::z, &NormalModesProxy::set_z)
      // NormalModesProxy.lin (0D_NOT_type -
      .def_property("lin", &NormalModesProxy::lin, &NormalModesProxy::set_lin)

      .def(
          "__repr__",
          [](const NormalModesProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NormalModesProxy& self) {
            return NormalModesProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const NormalModesProxy& self, py::dict& memo) {
            return NormalModesProxy(self);
          })

      ;

  // 1D NormalModesProxy arrays are not used in structs/routines
  // 2D NormalModesProxy arrays are not used in structs/routines
  // 3D NormalModesProxy arrays are not used in structs/routines
}

// =============================================================================
// nametable_struct
void init_nametable_struct(py::module& m, py::class_<NametableProxy>& cls) {
  cls.def(py::init<>())
      // NametableProxy.name (1D_ALLOC_character - Array of names.
      .def_property_readonly("name", &NametableProxy::name)
      // NametableProxy.index (1D_ALLOC_integer - Sorted index for names(:) array. names(an_index(i)) is in alphabetical order.
      .def_property_readonly("index", &NametableProxy::index)
      // NametableProxy.n_min (0D_NOT_integer - Set to 0 for use in a lattice.
      .def_property("n_min", &NametableProxy::n_min, &NametableProxy::set_n_min)
      // NametableProxy.n_max (0D_NOT_integer - Use only names(n_min:n_max) part of array.
      .def_property("n_max", &NametableProxy::n_max, &NametableProxy::set_n_max)

      .def(
          "__repr__",
          [](const NametableProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NametableProxy& self) {
            return NametableProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const NametableProxy& self, py::dict& memo) {
            return NametableProxy(self);
          })

      ;

  // 1D NametableProxy arrays are not used in structs/routines
  // 2D NametableProxy arrays are not used in structs/routines
  // 3D NametableProxy arrays are not used in structs/routines
}