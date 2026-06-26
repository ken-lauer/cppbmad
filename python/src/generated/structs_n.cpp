#include "pybmad/generated/structs_n.hpp"

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
// normal_modes_struct
void init_normal_modes_struct(nb::module_ &m, nb::class_<NormalModesStruct> &cls) {
  cls.def(
         "__init__",
         [](NormalModesStruct *self,
            std::optional<std::vector<double>> synch_int,
            std::optional<double> sigE_E,
            std::optional<double> sig_z,
            std::optional<double> e_loss,
            std::optional<double> rf_voltage,
            std::optional<double> pz_aperture,
            std::optional<double> pz_average,
            std::optional<double> momentum_compaction,
            std::optional<double> dpz_damp,
            const AnormalModeStruct *a,
            const AnormalModeStruct *b,
            const AnormalModeStruct *z,
            const LinacNormalModeStruct *lin) {
           new (self) NormalModesStruct(
               synch_int,
               sigE_E,
               sig_z,
               e_loss,
               rf_voltage,
               pz_aperture,
               pz_average,
               momentum_compaction,
               dpz_damp,
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(z),
               ptr_to_opt_ref(lin)
           );
         },
         nb::arg("synch_int") = nb::none(),
         nb::arg("sigE_E") = nb::none(),
         nb::arg("sig_z") = nb::none(),
         nb::arg("e_loss") = nb::none(),
         nb::arg("rf_voltage") = nb::none(),
         nb::arg("pz_aperture") = nb::none(),
         nb::arg("pz_average") = nb::none(),
         nb::arg("momentum_compaction") = nb::none(),
         nb::arg("dpz_damp") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("z") = nb::none(),
         nb::arg("lin") = nb::none()
  )
      .def_prop_rw(
          "synch_int",
          &NormalModesStruct::synch_int,
          &NormalModesStruct::set_synch_int,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Synchrotron integrals I0, I1, I2, and I3"
      )
      .def_prop_rw("sigE_E", &NormalModesStruct::sigE_E, &NormalModesStruct::set_sigE_E, "SigmaE/E")
      .def_prop_rw("sig_z", &NormalModesStruct::sig_z, &NormalModesStruct::set_sig_z, "Sigma_Z")
      .def_prop_rw(
          "e_loss",
          &NormalModesStruct::e_loss,
          &NormalModesStruct::set_e_loss,
          "Energy loss / turn (eV)"
      )
      .def_prop_rw(
          "rf_voltage",
          &NormalModesStruct::rf_voltage,
          &NormalModesStruct::set_rf_voltage,
          "Total rfcavity voltage (eV)"
      )
      .def_prop_rw(
          "pz_aperture",
          &NormalModesStruct::pz_aperture,
          &NormalModesStruct::set_pz_aperture,
          "pz aperture limit. Used with Touschek calculations."
      )
      .def_prop_rw(
          "pz_average",
          &NormalModesStruct::pz_average,
          &NormalModesStruct::set_pz_average,
          "Average over branch due to damping."
      )
      .def_prop_rw(
          "momentum_compaction",
          &NormalModesStruct::momentum_compaction,
          &NormalModesStruct::set_momentum_compaction
      )
      .def_prop_rw(
          "dpz_damp",
          &NormalModesStruct::dpz_damp,
          &NormalModesStruct::set_dpz_damp,
          "Change in pz without RF"
      )
      .def_prop_rw(
          "a",
          &NormalModesStruct::a,
          &NormalModesStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b",
          &NormalModesStruct::b,
          &NormalModesStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "z",
          &NormalModesStruct::z,
          &NormalModesStruct::set_z,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "lin",
          &NormalModesStruct::lin,
          &NormalModesStruct::set_lin,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const NormalModesStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NormalModesStruct &self) {
            return NormalModesStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const NormalModesStruct &self, nb::dict &memo) { return NormalModesStruct(self); }
      )
      .def(
          "__eq__",
          [](const NormalModesStruct &self, const NormalModesStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const NormalModesStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D NormalModesStruct arrays are not used in structs/routines
  // 2D NormalModesStruct arrays are not used in structs/routines
  // 3D NormalModesStruct arrays are not used in structs/routines
}

// =============================================================================
// nametable_struct
void init_nametable_struct(nb::module_ &m, nb::class_<NametableStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::vector<int>>, std::optional<int>, std::optional<int>>(),
         nb::arg("index") = nb::none(),
         nb::arg("n_min") = nb::none(),
         nb::arg("n_max") = nb::none()
  )
      .def_prop_ro("name", &NametableStruct::name, nb::keep_alive<0, 1>(), "Array of names.")
      .def_prop_rw(
          "index",
          &NametableStruct::index,
          &NametableStruct::set_index,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Sorted index for names(:) array. names(an_index(i)) is in alphabetical order."
      )
      .def_prop_rw(
          "n_min",
          &NametableStruct::n_min,
          &NametableStruct::set_n_min,
          "Set to 0 for use in a lattice."
      )
      .def_prop_rw(
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
          [](const NametableStruct &self, nb::dict &memo) { return NametableStruct(self); }
      )
      .def(
          "__eq__",
          [](const NametableStruct &self, const NametableStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const NametableStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D NametableStruct arrays are not used in structs/routines
  // 2D NametableStruct arrays are not used in structs/routines
  // 3D NametableStruct arrays are not used in structs/routines
}