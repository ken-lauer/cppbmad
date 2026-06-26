#include "pybmad/generated/structs_m.hpp"

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
// mode3_struct
void init_mode3_struct(nb::module_ &m, nb::class_<Mode3Struct> &cls) {
  cls.def(
         "__init__",
         [](Mode3Struct *self,
            std::optional<std::vector<std::vector<double>>> v,
            const TwissStruct *a,
            const TwissStruct *b,
            const TwissStruct *c,
            const TwissStruct *x,
            const TwissStruct *y) {
           new (self) Mode3Struct(
               v,
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(c),
               ptr_to_opt_ref(x),
               ptr_to_opt_ref(y)
           );
         },
         nb::arg("v") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("c") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none()
  )
      .def_prop_rw(
          "v",
          &Mode3Struct::v,
          &Mode3Struct::set_v,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "a",
          &Mode3Struct::a,
          &Mode3Struct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b",
          &Mode3Struct::b,
          &Mode3Struct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "c",
          &Mode3Struct::c,
          &Mode3Struct::set_c,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "x",
          &Mode3Struct::x,
          &Mode3Struct::set_x,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "y",
          &Mode3Struct::y,
          &Mode3Struct::set_y,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const Mode3Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mode3Struct &self) {
            return Mode3Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Mode3Struct &self, nb::dict &memo) { return Mode3Struct(self); }
      )
      .def(
          "__eq__",
          [](const Mode3Struct &self, const Mode3Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Mode3Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D Mode3Struct arrays are not used in structs/routines
  // 2D Mode3Struct arrays are not used in structs/routines
  // 3D Mode3Struct arrays are not used in structs/routines
}

// =============================================================================
// mode_info_struct
void init_mode_info_struct(nb::module_ &m, nb::class_<ModeInfoStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("stable") = nb::none(),
         nb::arg("tune") = nb::none(),
         nb::arg("emit") = nb::none(),
         nb::arg("chrom") = nb::none(),
         nb::arg("sigma") = nb::none(),
         nb::arg("sigmap") = nb::none()
  )
      .def_prop_rw(
          "stable",
          &ModeInfoStruct::stable,
          &ModeInfoStruct::set_stable,
          "Is the mode stable?"
      )
      .def_prop_rw(
          "tune",
          &ModeInfoStruct::tune,
          &ModeInfoStruct::set_tune,
          "'fractional' tune in radians"
      )
      .def_prop_rw(
          "emit",
          &ModeInfoStruct::emit,
          &ModeInfoStruct::set_emit,
          "Emittance (unnormalized)."
      )
      .def_prop_rw("chrom", &ModeInfoStruct::chrom, &ModeInfoStruct::set_chrom, "Chromaticity.")
      .def_prop_rw("sigma", &ModeInfoStruct::sigma, &ModeInfoStruct::set_sigma, "Beam size.")
      .def_prop_rw(
          "sigmap",
          &ModeInfoStruct::sigmap,
          &ModeInfoStruct::set_sigmap,
          "Beam divergence."
      )

      .def("__repr__", [](const ModeInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ModeInfoStruct &self) {
            return ModeInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ModeInfoStruct &self, nb::dict &memo) { return ModeInfoStruct(self); }
      )
      .def(
          "__eq__",
          [](const ModeInfoStruct &self, const ModeInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ModeInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ModeInfoStruct arrays are not used in structs/routines
  // 2D ModeInfoStruct arrays are not used in structs/routines
  // 3D ModeInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_energy_struct
void init_mad_energy_struct(nb::module_ &m, nb::class_<MadEnergyStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("total") = nb::none(),
         nb::arg("beta") = nb::none(),
         nb::arg("gamma") = nb::none(),
         nb::arg("kinetic") = nb::none(),
         nb::arg("p0c") = nb::none(),
         nb::arg("particle") = nb::none()
  )
      .def_prop_rw("total", &MadEnergyStruct::total, &MadEnergyStruct::set_total)
      .def_prop_rw(
          "beta",
          &MadEnergyStruct::beta,
          &MadEnergyStruct::set_beta,
          "normalized velocity: v/c"
      )
      .def_prop_rw(
          "gamma",
          &MadEnergyStruct::gamma,
          &MadEnergyStruct::set_gamma,
          "relativistic factor: 1/sqrt(1-beta^2)"
      )
      .def_prop_rw(
          "kinetic",
          &MadEnergyStruct::kinetic,
          &MadEnergyStruct::set_kinetic,
          "kinetic energy"
      )
      .def_prop_rw("p0c", &MadEnergyStruct::p0c, &MadEnergyStruct::set_p0c, "particle momentum")
      .def_prop_rw(
          "particle",
          &MadEnergyStruct::particle,
          &MadEnergyStruct::set_particle,
          "particle species"
      )

      .def("__repr__", [](const MadEnergyStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadEnergyStruct &self) {
            return MadEnergyStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MadEnergyStruct &self, nb::dict &memo) { return MadEnergyStruct(self); }
      )
      .def(
          "__eq__",
          [](const MadEnergyStruct &self, const MadEnergyStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MadEnergyStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MadEnergyStruct arrays are not used in structs/routines
  // 2D MadEnergyStruct arrays are not used in structs/routines
  // 3D MadEnergyStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_map_struct
void init_mad_map_struct(nb::module_ &m, nb::class_<MadMapStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>>(),
         nb::arg("k") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("t") = nb::none()
  )
      .def_prop_rw(
          "k",
          &MadMapStruct::k,
          &MadMapStruct::set_k,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th order map."
      )
      .def_prop_rw(
          "r",
          &MadMapStruct::r,
          &MadMapStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1st order map."
      )
      .def_prop_rw(
          "t",
          &MadMapStruct::t,
          &MadMapStruct::set_t,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "2nd order map."
      )

      .def("__repr__", [](const MadMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadMapStruct &self) {
            return MadMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MadMapStruct &self, nb::dict &memo) { return MadMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const MadMapStruct &self, const MadMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MadMapStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MadMapStruct arrays are not used in structs/routines
  // 2D MadMapStruct arrays are not used in structs/routines
  // 3D MadMapStruct arrays are not used in structs/routines
}