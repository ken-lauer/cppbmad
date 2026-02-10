#include "pybmad/generated/structs_m.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// mode3_struct
void init_mode3_struct(py::module &m, py::class_<Mode3Struct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<double>>>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>>(),
         py::arg("v") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("c") = py::none(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none()
  )
      .def_property("v", &Mode3Struct::v, &Mode3Struct::set_v, py::keep_alive<0, 1>())
      .def_property("a", &Mode3Struct::a, &Mode3Struct::set_a, py::keep_alive<0, 1>())
      .def_property("b", &Mode3Struct::b, &Mode3Struct::set_b, py::keep_alive<0, 1>())
      .def_property("c", &Mode3Struct::c, &Mode3Struct::set_c, py::keep_alive<0, 1>())
      .def_property("x", &Mode3Struct::x, &Mode3Struct::set_x, py::keep_alive<0, 1>())
      .def_property("y", &Mode3Struct::y, &Mode3Struct::set_y, py::keep_alive<0, 1>())

      .def("__repr__", [](const Mode3Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mode3Struct &self) {
            return Mode3Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Mode3Struct &self, py::dict &memo) { return Mode3Struct(self); }
      )

      ;

  // 1D Mode3Struct arrays are not used in structs/routines
  // 2D Mode3Struct arrays are not used in structs/routines
  // 3D Mode3Struct arrays are not used in structs/routines
}

// =============================================================================
// mode_info_struct
void init_mode_info_struct(py::module &m, py::class_<ModeInfoStruct> &cls) {
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("stable") = py::none(),
         py::arg("tune") = py::none(),
         py::arg("emit") = py::none(),
         py::arg("chrom") = py::none(),
         py::arg("sigma") = py::none(),
         py::arg("sigmap") = py::none()
  )
      .def_property(
          "stable",
          &ModeInfoStruct::stable,
          &ModeInfoStruct::set_stable,
          "Is the mode stable?"
      )
      .def_property(
          "tune",
          &ModeInfoStruct::tune,
          &ModeInfoStruct::set_tune,
          "'fractional' tune in radians"
      )
      .def_property(
          "emit",
          &ModeInfoStruct::emit,
          &ModeInfoStruct::set_emit,
          "Emittance (unnormalized)."
      )
      .def_property("chrom", &ModeInfoStruct::chrom, &ModeInfoStruct::set_chrom, "Chromaticity.")
      .def_property("sigma", &ModeInfoStruct::sigma, &ModeInfoStruct::set_sigma, "Beam size.")
      .def_property(
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
          [](const ModeInfoStruct &self, py::dict &memo) { return ModeInfoStruct(self); }
      )

      ;

  // 1D ModeInfoStruct arrays are not used in structs/routines
  // 2D ModeInfoStruct arrays are not used in structs/routines
  // 3D ModeInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_energy_struct
void init_mad_energy_struct(py::module &m, py::class_<MadEnergyStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         py::arg("total") = py::none(),
         py::arg("beta") = py::none(),
         py::arg("gamma") = py::none(),
         py::arg("kinetic") = py::none(),
         py::arg("p0c") = py::none(),
         py::arg("particle") = py::none()
  )
      .def_property("total", &MadEnergyStruct::total, &MadEnergyStruct::set_total)
      .def_property(
          "beta",
          &MadEnergyStruct::beta,
          &MadEnergyStruct::set_beta,
          "normalized velocity: v/c"
      )
      .def_property(
          "gamma",
          &MadEnergyStruct::gamma,
          &MadEnergyStruct::set_gamma,
          "relativistic factor: 1/sqrt(1-beta^2)"
      )
      .def_property(
          "kinetic",
          &MadEnergyStruct::kinetic,
          &MadEnergyStruct::set_kinetic,
          "kinetic energy"
      )
      .def_property("p0c", &MadEnergyStruct::p0c, &MadEnergyStruct::set_p0c, "particle momentum")
      .def_property(
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
          [](const MadEnergyStruct &self, py::dict &memo) { return MadEnergyStruct(self); }
      )

      ;

  // 1D MadEnergyStruct arrays are not used in structs/routines
  // 2D MadEnergyStruct arrays are not used in structs/routines
  // 3D MadEnergyStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_map_struct
void init_mad_map_struct(py::module &m, py::class_<MadMapStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>>(),
         py::arg("k") = py::none(),
         py::arg("r") = py::none(),
         py::arg("t") = py::none()
  )
      .def_property(
          "k",
          &MadMapStruct::k,
          &MadMapStruct::set_k,
          py::keep_alive<0, 1>(),
          "0th order map."
      )
      .def_property(
          "r",
          &MadMapStruct::r,
          &MadMapStruct::set_r,
          py::keep_alive<0, 1>(),
          "1st order map."
      )
      .def_property(
          "t",
          &MadMapStruct::t,
          &MadMapStruct::set_t,
          py::keep_alive<0, 1>(),
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
          [](const MadMapStruct &self, py::dict &memo) { return MadMapStruct(self); }
      )

      ;

  // 1D MadMapStruct arrays are not used in structs/routines
  // 2D MadMapStruct arrays are not used in structs/routines
  // 3D MadMapStruct arrays are not used in structs/routines
}