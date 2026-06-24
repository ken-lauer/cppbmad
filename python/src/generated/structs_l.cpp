#include "pybmad/generated/structs_l.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// lat_ele_loc_struct
void init_lat_ele_loc_struct(py::module &m, py::class_<LatEleLocStruct> &cls) {
  cls.def(
         py::init<std::optional<int>, std::optional<int>>(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_branch") = py::none()
  )
      .def_property("ix_ele", &LatEleLocStruct::ix_ele, &LatEleLocStruct::set_ix_ele)
      .def_property("ix_branch", &LatEleLocStruct::ix_branch, &LatEleLocStruct::set_ix_branch)
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleLocStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleLocStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const LatEleLocStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleLocStruct &self) {
            return LatEleLocStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatEleLocStruct &self, py::dict &memo) { return LatEleLocStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleLocStruct &self, const LatEleLocStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleLocStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<LatEleLocStructArray1D, LatEleLocStructAlloc1D>(
      m,
      "LatEleLocStructArray1D",
      "LatEleLocStructAlloc1D"
  );
  // 2D LatEleLocStruct arrays are not used in structs/routines
  // 3D LatEleLocStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order1_struct
void init_lat_ele_order1_struct(py::module &m, py::class_<LatEleOrder1Struct> &cls) {
  cls.def(
         py::init<std::optional<int>, std::optional<int>>(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_order") = py::none()
  )
      .def_property(
          "ix_branch",
          &LatEleOrder1Struct::ix_branch,
          &LatEleOrder1Struct::set_ix_branch,
          "Branch index"
      )
      .def_property(
          "ix_order",
          &LatEleOrder1Struct::ix_order,
          &LatEleOrder1Struct::set_ix_order,
          "Order index. -1 -> Unique in lattice, 0 -> unique in branch."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleOrder1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleOrder1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const LatEleOrder1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrder1Struct &self) {
            return LatEleOrder1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatEleOrder1Struct &self, py::dict &memo) { return LatEleOrder1Struct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleOrder1Struct &self, const LatEleOrder1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrder1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<LatEleOrder1StructArray1D, LatEleOrder1StructAlloc1D>(
      m,
      "LatEleOrder1StructArray1D",
      "LatEleOrder1StructAlloc1D"
  );
  // 2D LatEleOrder1Struct arrays are not used in structs/routines
  // 3D LatEleOrder1Struct arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_array_struct
void init_lat_ele_order_array_struct(py::module &m, py::class_<LatEleOrderArrayStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly(
          "ele",
          py::cpp_function(&LatEleOrderArrayStruct::ele, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleOrderArrayStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleOrderArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const LatEleOrderArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderArrayStruct &self) {
            return LatEleOrderArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatEleOrderArrayStruct &self, py::dict &memo) {
            return LatEleOrderArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const LatEleOrderArrayStruct &self, const LatEleOrderArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrderArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<LatEleOrderArrayStructArray1D, LatEleOrderArrayStructAlloc1D>(
      m,
      "LatEleOrderArrayStructArray1D",
      "LatEleOrderArrayStructAlloc1D"
  );
  // 2D LatEleOrderArrayStruct arrays are not used in structs/routines
  // 3D LatEleOrderArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_struct
void init_lat_ele_order_struct(py::module &m, py::class_<LatEleOrderStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly(
          "branch",
          py::cpp_function(&LatEleOrderStruct::branch, py::keep_alive<0, 1>())
      )

      .def("__repr__", [](const LatEleOrderStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderStruct &self) {
            return LatEleOrderStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatEleOrderStruct &self, py::dict &memo) { return LatEleOrderStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleOrderStruct &self, const LatEleOrderStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrderStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D LatEleOrderStruct arrays are not used in structs/routines
  // 2D LatEleOrderStruct arrays are not used in structs/routines
  // 3D LatEleOrderStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_param_struct
void init_lat_param_struct(py::module &m, py::class_<LatParamStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const BookkeepingStateStruct>,
             optional_ref<const BeamInitStruct>>(),
         py::arg("n_part") = py::none(),
         py::arg("total_length") = py::none(),
         py::arg("unstable_factor") = py::none(),
         py::arg("t1_with_RF") = py::none(),
         py::arg("t1_no_RF") = py::none(),
         py::arg("spin_tune") = py::none(),
         py::arg("particle") = py::none(),
         py::arg("default_tracking_species") = py::none(),
         py::arg("geometry") = py::none(),
         py::arg("ixx") = py::none(),
         py::arg("stable") = py::none(),
         py::arg("live_branch") = py::none(),
         py::arg("g1_integral") = py::none(),
         py::arg("g2_integral") = py::none(),
         py::arg("g3_integral") = py::none(),
         py::arg("bookkeeping_state") = py::none(),
         py::arg("beam_init") = py::none()
  )
      .def_property(
          "n_part",
          &LatParamStruct::n_part,
          &LatParamStruct::set_n_part,
          "Particles/bunch (for BeamBeam elements)."
      )
      .def_property(
          "total_length",
          &LatParamStruct::total_length,
          &LatParamStruct::set_total_length,
          "total_length of branch. Warning: branch may not start at s = 0."
      )
      .def_property(
          "unstable_factor",
          &LatParamStruct::unstable_factor,
          &LatParamStruct::set_unstable_factor,
          "If positive: Growth rate/turn if unstable in closed branches or "
          "|orbit-aperture|/aperture if particle hits wall. Zero otherwise."
      )
      .def_property(
          "t1_with_RF",
          py::cpp_function(&LatParamStruct::t1_with_RF, py::keep_alive<0, 1>()),
          &LatParamStruct::set_t1_with_RF,
          "Full 1-turn matrix with RF on."
      )
      .def_property(
          "t1_no_RF",
          py::cpp_function(&LatParamStruct::t1_no_RF, py::keep_alive<0, 1>()),
          &LatParamStruct::set_t1_no_RF,
          "Full 1-turn matrix with RF off."
      )
      .def_property(
          "spin_tune",
          &LatParamStruct::spin_tune,
          &LatParamStruct::set_spin_tune,
          "Closed orbit spin tune."
      )
      .def_property(
          "particle",
          &LatParamStruct::particle,
          &LatParamStruct::set_particle,
          "Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is "
          "changed."
      )
      .def_property(
          "default_tracking_species",
          &LatParamStruct::default_tracking_species,
          &LatParamStruct::set_default_tracking_species,
          "Default particle type to use in tracking."
      )
      .def_property(
          "geometry",
          &LatParamStruct::geometry,
          &LatParamStruct::set_geometry,
          "open$ or closed$"
      )
      .def_property(
          "ixx",
          &LatParamStruct::ixx,
          &LatParamStruct::set_ixx,
          "Integer for general use"
      )
      .def_property(
          "stable",
          &LatParamStruct::stable,
          &LatParamStruct::set_stable,
          "is closed lat stable?"
      )
      .def_property(
          "live_branch",
          &LatParamStruct::live_branch,
          &LatParamStruct::set_live_branch,
          "Should tracking be done on the branch?"
      )
      .def_property(
          "g1_integral",
          &LatParamStruct::g1_integral,
          &LatParamStruct::set_g1_integral,
          "Approximate |g| (bending strength) integral of branch."
      )
      .def_property(
          "g2_integral",
          &LatParamStruct::g2_integral,
          &LatParamStruct::set_g2_integral,
          "Approximate g^2 integral of branch."
      )
      .def_property(
          "g3_integral",
          &LatParamStruct::g3_integral,
          &LatParamStruct::set_g3_integral,
          "Approximate g^2 integral of branch."
      )
      .def_property(
          "bookkeeping_state",
          py::cpp_function(&LatParamStruct::bookkeeping_state, py::keep_alive<0, 1>()),
          &LatParamStruct::set_bookkeeping_state,
          "Overall status for the branch."
      )
      .def_property(
          "beam_init",
          py::cpp_function(&LatParamStruct::beam_init, py::keep_alive<0, 1>()),
          &LatParamStruct::set_beam_init,
          "For beam initialization."
      )

      .def("__repr__", [](const LatParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatParamStruct &self) {
            return LatParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatParamStruct &self, py::dict &memo) { return LatParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatParamStruct &self, const LatParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatParamStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D LatParamStruct arrays are not used in structs/routines
  // 2D LatParamStruct arrays are not used in structs/routines
  // 3D LatParamStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_struct
void init_lat_struct(py::module &m, py::class_<LatStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const LatParamStruct>,
             optional_ref<const BookkeepingStateStruct>,
             optional_ref<const EleStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const BeamInitStruct>,
             optional_ref<const PreTrackerStruct>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<int>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("use_name") = py::none(),
         py::arg("lattice") = py::none(),
         py::arg("machine") = py::none(),
         py::arg("input_file_name") = py::none(),
         py::arg("title") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("z") = py::none(),
         py::arg("param") = py::none(),
         py::arg("lord_state") = py::none(),
         py::arg("ele_init") = py::none(),
         py::arg("particle_start") = py::none(),
         py::arg("beam_init") = py::none(),
         py::arg("pre_tracker") = py::none(),
         py::arg("custom") = py::none(),
         py::arg("version") = py::none(),
         py::arg("n_ele_track") = py::none(),
         py::arg("n_ele_max") = py::none(),
         py::arg("n_control_max") = py::none(),
         py::arg("n_ic_max") = py::none(),
         py::arg("input_taylor_order") = py::none(),
         py::arg("ic") = py::none(),
         py::arg("photon_type") = py::none(),
         py::arg("creation_hash") = py::none(),
         py::arg("ramper_slave_bookkeeping") = py::none(),
         py::arg("parser_make_xfer_mats") = py::none()
  )
      .def_property(
          "use_name",
          &LatStruct::use_name,
          &LatStruct::set_use_name,
          "Name of lat given by USE statement"
      )
      .def_property("lattice", &LatStruct::lattice, &LatStruct::set_lattice, "Lattice")
      .def_property(
          "machine",
          &LatStruct::machine,
          &LatStruct::set_machine,
          "Name of the machine the lattice is for ('LHC', etc)."
      )
      .def_property(
          "input_file_name",
          &LatStruct::input_file_name,
          &LatStruct::set_input_file_name,
          "Name of the lattice input file"
      )
      .def_property("title", &LatStruct::title, &LatStruct::set_title, "General title")
      .def_property_readonly(
          "print_str",
          py::cpp_function(&LatStruct::print_str, py::keep_alive<0, 1>()),
          "Saved print statements."
      )
      .def_property_readonly(
          "constant",
          py::cpp_function(&LatStruct::constant, py::keep_alive<0, 1>()),
          "Constants defined in the lattice"
      )
      .def_property(
          "a",
          py::cpp_function(&LatStruct::a, py::keep_alive<0, 1>()),
          &LatStruct::set_a,
          "Tunes (fractional part), etc."
      )
      .def_property(
          "b",
          py::cpp_function(&LatStruct::b, py::keep_alive<0, 1>()),
          &LatStruct::set_b,
          "Tunes (fractional part), etc."
      )
      .def_property(
          "z",
          py::cpp_function(&LatStruct::z, py::keep_alive<0, 1>()),
          &LatStruct::set_z,
          "Tunes (fractional part), etc."
      )
      .def_property(
          "param",
          py::cpp_function(&LatStruct::param, py::keep_alive<0, 1>()),
          &LatStruct::set_param,
          "Parameters"
      )
      .def_property(
          "lord_state",
          py::cpp_function(&LatStruct::lord_state, py::keep_alive<0, 1>()),
          &LatStruct::set_lord_state,
          "lord bookkeeping status."
      )
      .def_property(
          "ele_init",
          py::cpp_function(&LatStruct::ele_init, py::keep_alive<0, 1>()),
          &LatStruct::set_ele_init,
          "For use by any program"
      )
      .def_property_readonly(
          "ele",
          py::cpp_function(&LatStruct::ele, py::keep_alive<0, 1>()),
          "Array of elements [=> branch(0)]."
      )
      .def_property_readonly(
          "branch",
          py::cpp_function(&LatStruct::branch, py::keep_alive<0, 1>()),
          "Branch(0:) array"
      )
      .def_property_readonly(
          "control",
          py::cpp_function(&LatStruct::control, py::keep_alive<0, 1>()),
          "Control list"
      )
      .def_property(
          "particle_start",
          py::cpp_function(&LatStruct::particle_start, py::keep_alive<0, 1>()),
          &LatStruct::set_particle_start,
          "Starting particle_coords."
      )
      .def_property(
          "beam_init",
          py::cpp_function(&LatStruct::beam_init, py::keep_alive<0, 1>()),
          &LatStruct::set_beam_init,
          "Beam initialization."
      )
      .def_property(
          "pre_tracker",
          py::cpp_function(&LatStruct::pre_tracker, py::keep_alive<0, 1>()),
          &LatStruct::set_pre_tracker,
          "For OPAL/IMPACT-T"
      )
      .def_property(
          "custom",
          py::cpp_function(&LatStruct::custom, py::keep_alive<0, 1>()),
          &LatStruct::set_custom,
          "Custom attributes."
      )
      .def_property("version", &LatStruct::version, &LatStruct::set_version, "Version number")
      .def_property(
          "n_ele_track",
          &LatStruct::n_ele_track,
          &LatStruct::set_n_ele_track,
          "Number of lat elements to track through."
      )
      .def_property(
          "n_ele_max",
          &LatStruct::n_ele_max,
          &LatStruct::set_n_ele_max,
          "Index of last valid element in %ele(:) array"
      )
      .def_property(
          "n_control_max",
          &LatStruct::n_control_max,
          &LatStruct::set_n_control_max,
          "Last index used in control_array"
      )
      .def_property(
          "n_ic_max",
          &LatStruct::n_ic_max,
          &LatStruct::set_n_ic_max,
          "Last index used in ic_array"
      )
      .def_property(
          "input_taylor_order",
          &LatStruct::input_taylor_order,
          &LatStruct::set_input_taylor_order,
          "As set in the input file"
      )
      .def_property(
          "ic",
          py::cpp_function(&LatStruct::ic, py::keep_alive<0, 1>()),
          &LatStruct::set_ic,
          "Index to %control(:) from slaves."
      )
      .def_property(
          "photon_type",
          &LatStruct::photon_type,
          &LatStruct::set_photon_type,
          "Or coherent$. For X-ray simulations."
      )
      .def_property(
          "creation_hash",
          &LatStruct::creation_hash,
          &LatStruct::set_creation_hash,
          "Set by bmad_parser. creation_hash will vary if any of the lattice files are modified."
      )
      .def_property(
          "ramper_slave_bookkeeping",
          &LatStruct::ramper_slave_bookkeeping,
          &LatStruct::set_ramper_slave_bookkeeping
      )
      .def_property(
          "parser_make_xfer_mats",
          &LatStruct::parser_make_xfer_mats,
          &LatStruct::set_parser_make_xfer_mats,
          "Is Bmad parser to make element transfer matrices?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const LatStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatStruct &self) {
            return LatStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const LatStruct &self, py::dict &memo) { return LatStruct(self); })
      .def(
          "__eq__",
          [](const LatStruct &self, const LatStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LatStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<LatStructArray1D, LatStructAlloc1D>(
      m,
      "LatStructArray1D",
      "LatStructAlloc1D"
  );
  // 2D LatStruct arrays are not used in structs/routines
  // 3D LatStruct arrays are not used in structs/routines
}

// =============================================================================
// linac_normal_mode_struct
void init_linac_normal_mode_struct(py::module &m, py::class_<LinacNormalModeStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("i2_E4") = py::none(),
         py::arg("i3_E7") = py::none(),
         py::arg("i5a_E6") = py::none(),
         py::arg("i5b_E6") = py::none(),
         py::arg("sig_E1") = py::none(),
         py::arg("a_emittance_end") = py::none(),
         py::arg("b_emittance_end") = py::none()
  )
      .def_property(
          "i2_E4",
          &LinacNormalModeStruct::i2_E4,
          &LinacNormalModeStruct::set_i2_E4,
          "Integral: g^2 * gamma^4"
      )
      .def_property(
          "i3_E7",
          &LinacNormalModeStruct::i3_E7,
          &LinacNormalModeStruct::set_i3_E7,
          "Integral: g^3 * gamma^7"
      )
      .def_property(
          "i5a_E6",
          &LinacNormalModeStruct::i5a_E6,
          &LinacNormalModeStruct::set_i5a_E6,
          "Integral: (g^3 * H_a) * gamma^6"
      )
      .def_property(
          "i5b_E6",
          &LinacNormalModeStruct::i5b_E6,
          &LinacNormalModeStruct::set_i5b_E6,
          "Integral: (g^3 * H_b) * gamma^6"
      )
      .def_property(
          "sig_E1",
          &LinacNormalModeStruct::sig_E1,
          &LinacNormalModeStruct::set_sig_E1,
          "Energy spread after 1 pass (eV)"
      )
      .def_property(
          "a_emittance_end",
          &LinacNormalModeStruct::a_emittance_end,
          &LinacNormalModeStruct::set_a_emittance_end,
          "a mode emittance at end of linac"
      )
      .def_property(
          "b_emittance_end",
          &LinacNormalModeStruct::b_emittance_end,
          &LinacNormalModeStruct::set_b_emittance_end,
          "b mode emittance at end of linac"
      )

      .def("__repr__", [](const LinacNormalModeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LinacNormalModeStruct &self) {
            return LinacNormalModeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LinacNormalModeStruct &self, py::dict &memo) {
            return LinacNormalModeStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const LinacNormalModeStruct &self, const LinacNormalModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const LinacNormalModeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D LinacNormalModeStruct arrays are not used in structs/routines
  // 2D LinacNormalModeStruct arrays are not used in structs/routines
  // 3D LinacNormalModeStruct arrays are not used in structs/routines
}

// =============================================================================
// layout
void init_layout(py::module &m, py::class_<Layout> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<int>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             optional_ref<const Layout>,
             optional_ref<const Layout>>(),
         py::arg("NAME") = py::none(),
         py::arg("INDEX") = py::none(),
         py::arg("HARMONIC_NUMBER") = py::none(),
         py::arg("CLOSED") = py::none(),
         py::arg("N") = py::none(),
         py::arg("NTHIN") = py::none(),
         py::arg("THIN") = py::none(),
         py::arg("LASTPOS") = py::none(),
         py::arg("LAST") = py::none(),
         py::arg("END") = py::none(),
         py::arg("START") = py::none(),
         py::arg("START_GROUND") = py::none(),
         py::arg("END_GROUND") = py::none(),
         py::arg("NEXT") = py::none(),
         py::arg("PREVIOUS") = py::none()
  )
      .def_property("NAME", &Layout::NAME, &Layout::set_NAME, "IDENTIFICATION")
      .def_property("INDEX", &Layout::INDEX, &Layout::set_INDEX, "IDENTIFICATION, CHARGE SIGN")
      .def_property("HARMONIC_NUMBER", &Layout::HARMONIC_NUMBER, &Layout::set_HARMONIC_NUMBER)
      .def_property("CLOSED", &Layout::CLOSED, &Layout::set_CLOSED)
      .def_property("N", &Layout::N, &Layout::set_N, "TOTAL ELEMENT IN THE CHAIN")
      .def_property(
          "NTHIN",
          &Layout::NTHIN,
          &Layout::set_NTHIN,
          "NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)"
      )
      .def_property(
          "THIN",
          &Layout::THIN,
          &Layout::set_THIN,
          "PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS POINTERS OF LINK LAYOUT"
      )
      .def_property("LASTPOS", &Layout::LASTPOS, &Layout::set_LASTPOS, "POSITION OF LAST VISITED")
      .def_property(
          "LAST",
          py::cpp_function(&Layout::LAST, py::keep_alive<0, 1>()),
          &Layout::set_LAST,
          "LAST VISITED"
      )
      .def_property("END", py::cpp_function(&Layout::END, py::keep_alive<0, 1>()), &Layout::set_END)
      .def_property(
          "START",
          py::cpp_function(&Layout::START, py::keep_alive<0, 1>()),
          &Layout::set_START
      )
      .def_property(
          "START_GROUND",
          py::cpp_function(&Layout::START_GROUND, py::keep_alive<0, 1>()),
          &Layout::set_START_GROUND,
          "STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING"
      )
      .def_property(
          "END_GROUND",
          py::cpp_function(&Layout::END_GROUND, py::keep_alive<0, 1>()),
          &Layout::set_END_GROUND,
          "STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING"
      )
      .def_property(
          "NEXT",
          py::cpp_function(&Layout::NEXT, py::keep_alive<0, 1>()),
          &Layout::set_NEXT
      )
      .def_property(
          "PREVIOUS",
          py::cpp_function(&Layout::PREVIOUS, py::keep_alive<0, 1>()),
          &Layout::set_PREVIOUS
      )

      .def("__repr__", [](const Layout &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Layout &self) {
            return Layout(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const Layout &self, py::dict &memo) { return Layout(self); })
      .def(
          "__eq__",
          [](const Layout &self, const Layout &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const Layout &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D Layout arrays are not used in structs/routines
  // 2D Layout arrays are not used in structs/routines
  // 3D Layout arrays are not used in structs/routines
}