#include "pybmad/generated/structs_n.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// normal_modes_struct
void init_normal_modes_struct(py::module &m, py::class_<NormalModesStruct> &cls) {
  cls.def(py::init<>())
      // NormalModesStruct.synch_int (1D_NOT_real - Synchrotron integrals I0, I1, I2, and I3
      .def_property("synch_int", &NormalModesStruct::synch_int, &NormalModesStruct::set_synch_int)
      // NormalModesStruct.sigE_E (0D_NOT_real - SigmaE/E
      .def_property("sigE_E", &NormalModesStruct::sigE_E, &NormalModesStruct::set_sigE_E)
      // NormalModesStruct.sig_z (0D_NOT_real - Sigma_Z
      .def_property("sig_z", &NormalModesStruct::sig_z, &NormalModesStruct::set_sig_z)
      // NormalModesStruct.e_loss (0D_NOT_real - Energy loss / turn (eV)
      .def_property("e_loss", &NormalModesStruct::e_loss, &NormalModesStruct::set_e_loss)
      // NormalModesStruct.rf_voltage (0D_NOT_real - Total rfcavity voltage (eV)
      .def_property(
          "rf_voltage",
          &NormalModesStruct::rf_voltage,
          &NormalModesStruct::set_rf_voltage
      )
      // NormalModesStruct.pz_aperture (0D_NOT_real - pz aperture limit. Used with Touschek
      // calculations.
      .def_property(
          "pz_aperture",
          &NormalModesStruct::pz_aperture,
          &NormalModesStruct::set_pz_aperture
      )
      // NormalModesStruct.pz_average (0D_NOT_real - Average over branch due to damping.
      .def_property(
          "pz_average",
          &NormalModesStruct::pz_average,
          &NormalModesStruct::set_pz_average
      )
      // NormalModesStruct.momentum_compaction (0D_NOT_real -
      .def_property(
          "momentum_compaction",
          &NormalModesStruct::momentum_compaction,
          &NormalModesStruct::set_momentum_compaction
      )
      // NormalModesStruct.dpz_damp (0D_NOT_real - Change in pz without RF
      .def_property("dpz_damp", &NormalModesStruct::dpz_damp, &NormalModesStruct::set_dpz_damp)
      // NormalModesStruct.a (0D_NOT_type -
      .def_property("a", &NormalModesStruct::a, &NormalModesStruct::set_a)
      // NormalModesStruct.b (0D_NOT_type -
      .def_property("b", &NormalModesStruct::b, &NormalModesStruct::set_b)
      // NormalModesStruct.z (0D_NOT_type -
      .def_property("z", &NormalModesStruct::z, &NormalModesStruct::set_z)
      // NormalModesStruct.lin (0D_NOT_type -
      .def_property("lin", &NormalModesStruct::lin, &NormalModesStruct::set_lin)

      .def("__repr__", [](const NormalModesStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NormalModesStruct &self) {
            return NormalModesStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const NormalModesStruct &self, py::dict &memo) { return NormalModesStruct(self); }
      )

      ;

  // 1D NormalModesStruct arrays are not used in structs/routines
  // 2D NormalModesStruct arrays are not used in structs/routines
  // 3D NormalModesStruct arrays are not used in structs/routines
}

// =============================================================================
// nametable_struct
void init_nametable_struct(py::module &m, py::class_<NametableStruct> &cls) {
  cls.def(py::init<>())
      // NametableStruct.name (1D_ALLOC_character - Array of names.
      .def_property_readonly("name", &NametableStruct::name)
      // NametableStruct.index (1D_ALLOC_integer - Sorted index for names(:) array.
      // names(an_index(i)) is in alphabetical order.
      .def_property("index", &NametableStruct::index, &NametableStruct::set_index)
      // NametableStruct.n_min (0D_NOT_integer - Set to 0 for use in a lattice.
      .def_property("n_min", &NametableStruct::n_min, &NametableStruct::set_n_min)
      // NametableStruct.n_max (0D_NOT_integer - Use only names(n_min:n_max) part of array.
      .def_property("n_max", &NametableStruct::n_max, &NametableStruct::set_n_max)

      .def("__repr__", [](const NametableStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NametableStruct &self) {
            return NametableStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const NametableStruct &self, py::dict &memo) { return NametableStruct(self); }
      )

      ;

  // 1D NametableStruct arrays are not used in structs/routines
  // 2D NametableStruct arrays are not used in structs/routines
  // 3D NametableStruct arrays are not used in structs/routines
}