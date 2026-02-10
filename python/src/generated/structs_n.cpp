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
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const AnormalModeStruct>,
             optional_ref<const AnormalModeStruct>,
             optional_ref<const AnormalModeStruct>,
             optional_ref<const LinacNormalModeStruct>>(),
         py::arg("synch_int") = py::none(),
         py::arg("sigE_E") = py::none(),
         py::arg("sig_z") = py::none(),
         py::arg("e_loss") = py::none(),
         py::arg("rf_voltage") = py::none(),
         py::arg("pz_aperture") = py::none(),
         py::arg("pz_average") = py::none(),
         py::arg("momentum_compaction") = py::none(),
         py::arg("dpz_damp") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("z") = py::none(),
         py::arg("lin") = py::none()
  )
      .def_property(
          "synch_int",
          &NormalModesStruct::synch_int,
          &NormalModesStruct::set_synch_int,
          "Synchrotron integrals I0, I1, I2, and I3"
      )
      .def_property(
          "sigE_E",
          &NormalModesStruct::sigE_E,
          &NormalModesStruct::set_sigE_E,
          "SigmaE/E"
      )
      .def_property("sig_z", &NormalModesStruct::sig_z, &NormalModesStruct::set_sig_z, "Sigma_Z")
      .def_property(
          "e_loss",
          &NormalModesStruct::e_loss,
          &NormalModesStruct::set_e_loss,
          "Energy loss / turn (eV)"
      )
      .def_property(
          "rf_voltage",
          &NormalModesStruct::rf_voltage,
          &NormalModesStruct::set_rf_voltage,
          "Total rfcavity voltage (eV)"
      )
      .def_property(
          "pz_aperture",
          &NormalModesStruct::pz_aperture,
          &NormalModesStruct::set_pz_aperture,
          "pz aperture limit. Used with Touschek calculations."
      )
      .def_property(
          "pz_average",
          &NormalModesStruct::pz_average,
          &NormalModesStruct::set_pz_average,
          "Average over branch due to damping."
      )
      .def_property(
          "momentum_compaction",
          &NormalModesStruct::momentum_compaction,
          &NormalModesStruct::set_momentum_compaction
      )
      .def_property(
          "dpz_damp",
          &NormalModesStruct::dpz_damp,
          &NormalModesStruct::set_dpz_damp,
          "Change in pz without RF"
      )
      .def_property("a", &NormalModesStruct::a, &NormalModesStruct::set_a)
      .def_property("b", &NormalModesStruct::b, &NormalModesStruct::set_b)
      .def_property("z", &NormalModesStruct::z, &NormalModesStruct::set_z)
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
  cls.def(
         py::init<std::optional<std::vector<int>>, std::optional<int>, std::optional<int>>(),
         py::arg("index") = py::none(),
         py::arg("n_min") = py::none(),
         py::arg("n_max") = py::none()
  )
      .def_property_readonly("name", &NametableStruct::name, "Array of names.")
      .def_property(
          "index",
          &NametableStruct::index,
          &NametableStruct::set_index,
          "Sorted index for names(:) array. names(an_index(i)) is in alphabetical order."
      )
      .def_property(
          "n_min",
          &NametableStruct::n_min,
          &NametableStruct::set_n_min,
          "Set to 0 for use in a lattice."
      )
      .def_property(
          "n_max",
          &NametableStruct::n_max,
          &NametableStruct::set_n_max,
          "Use only names(n_min:n_max) part of array."
      )

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