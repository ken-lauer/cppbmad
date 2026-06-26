#include "pybmad/generated/structs_l.hpp"

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
// lat_ele_loc_struct
void init_lat_ele_loc_struct(nb::module_ &m, nb::class_<LatEleLocStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<int>>(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_branch") = nb::none()
  )
      .def_prop_rw("ix_ele", &LatEleLocStruct::ix_ele, &LatEleLocStruct::set_ix_ele)
      .def_prop_rw("ix_branch", &LatEleLocStruct::ix_branch, &LatEleLocStruct::set_ix_branch)
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleLocStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleLocStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const LatEleLocStruct &self, nb::dict &memo) { return LatEleLocStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleLocStruct &self, const LatEleLocStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleLocStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_lat_ele_order1_struct(nb::module_ &m, nb::class_<LatEleOrder1Struct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<int>>(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_order") = nb::none()
  )
      .def_prop_rw(
          "ix_branch",
          &LatEleOrder1Struct::ix_branch,
          &LatEleOrder1Struct::set_ix_branch,
          "Branch index"
      )
      .def_prop_rw(
          "ix_order",
          &LatEleOrder1Struct::ix_order,
          &LatEleOrder1Struct::set_ix_order,
          "Order index. -1 -> Unique in lattice, 0 -> unique in branch."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleOrder1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleOrder1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const LatEleOrder1Struct &self, nb::dict &memo) { return LatEleOrder1Struct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleOrder1Struct &self, const LatEleOrder1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrder1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_lat_ele_order_array_struct(nb::module_ &m, nb::class_<LatEleOrderArrayStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("ele", &LatEleOrderArrayStruct::ele, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return LatEleOrderArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatEleOrderArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const LatEleOrderArrayStruct &self, nb::dict &memo) {
            return LatEleOrderArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const LatEleOrderArrayStruct &self, const LatEleOrderArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrderArrayStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_lat_ele_order_struct(nb::module_ &m, nb::class_<LatEleOrderStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("branch", &LatEleOrderStruct::branch, nb::keep_alive<0, 1>())

      .def("__repr__", [](const LatEleOrderStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderStruct &self) {
            return LatEleOrderStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatEleOrderStruct &self, nb::dict &memo) { return LatEleOrderStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatEleOrderStruct &self, const LatEleOrderStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatEleOrderStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D LatEleOrderStruct arrays are not used in structs/routines
  // 2D LatEleOrderStruct arrays are not used in structs/routines
  // 3D LatEleOrderStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_param_struct
void init_lat_param_struct(nb::module_ &m, nb::class_<LatParamStruct> &cls) {
  cls.def(
         "__init__",
         [](LatParamStruct *self,
            std::optional<double> n_part,
            std::optional<double> total_length,
            std::optional<double> unstable_factor,
            std::optional<std::vector<std::vector<double>>> t1_with_RF,
            std::optional<std::vector<std::vector<double>>> t1_no_RF,
            std::optional<double> spin_tune,
            std::optional<int> particle,
            std::optional<int> default_tracking_species,
            std::optional<int> geometry,
            std::optional<int> ixx,
            std::optional<bool> stable,
            std::optional<bool> live_branch,
            std::optional<double> g1_integral,
            std::optional<double> g2_integral,
            std::optional<double> g3_integral,
            const BookkeepingStateStruct *bookkeeping_state,
            const BeamInitStruct *beam_init) {
           new (self) LatParamStruct(
               n_part,
               total_length,
               unstable_factor,
               t1_with_RF,
               t1_no_RF,
               spin_tune,
               particle,
               default_tracking_species,
               geometry,
               ixx,
               stable,
               live_branch,
               g1_integral,
               g2_integral,
               g3_integral,
               ptr_to_opt_ref(bookkeeping_state),
               ptr_to_opt_ref(beam_init)
           );
         },
         nb::arg("n_part") = nb::none(),
         nb::arg("total_length") = nb::none(),
         nb::arg("unstable_factor") = nb::none(),
         nb::arg("t1_with_RF") = nb::none(),
         nb::arg("t1_no_RF") = nb::none(),
         nb::arg("spin_tune") = nb::none(),
         nb::arg("particle") = nb::none(),
         nb::arg("default_tracking_species") = nb::none(),
         nb::arg("geometry") = nb::none(),
         nb::arg("ixx") = nb::none(),
         nb::arg("stable") = nb::none(),
         nb::arg("live_branch") = nb::none(),
         nb::arg("g1_integral") = nb::none(),
         nb::arg("g2_integral") = nb::none(),
         nb::arg("g3_integral") = nb::none(),
         nb::arg("bookkeeping_state") = nb::none(),
         nb::arg("beam_init") = nb::none()
  )
      .def_prop_rw(
          "n_part",
          &LatParamStruct::n_part,
          &LatParamStruct::set_n_part,
          "Particles/bunch (for BeamBeam elements)."
      )
      .def_prop_rw(
          "total_length",
          &LatParamStruct::total_length,
          &LatParamStruct::set_total_length,
          "total_length of branch. Warning: branch may not start at s = 0."
      )
      .def_prop_rw(
          "unstable_factor",
          &LatParamStruct::unstable_factor,
          &LatParamStruct::set_unstable_factor,
          "If positive: Growth rate/turn if unstable in closed branches or "
          "|orbit-aperture|/aperture if particle hits wall. Zero otherwise."
      )
      .def_prop_rw(
          "t1_with_RF",
          &LatParamStruct::t1_with_RF,
          &LatParamStruct::set_t1_with_RF,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Full 1-turn matrix with RF on."
      )
      .def_prop_rw(
          "t1_no_RF",
          &LatParamStruct::t1_no_RF,
          &LatParamStruct::set_t1_no_RF,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Full 1-turn matrix with RF off."
      )
      .def_prop_rw(
          "spin_tune",
          &LatParamStruct::spin_tune,
          &LatParamStruct::set_spin_tune,
          "Closed orbit spin tune."
      )
      .def_prop_rw(
          "particle",
          &LatParamStruct::particle,
          &LatParamStruct::set_particle,
          "Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is "
          "changed."
      )
      .def_prop_rw(
          "default_tracking_species",
          &LatParamStruct::default_tracking_species,
          &LatParamStruct::set_default_tracking_species,
          "Default particle type to use in tracking."
      )
      .def_prop_rw(
          "geometry",
          &LatParamStruct::geometry,
          &LatParamStruct::set_geometry,
          "open$ or closed$"
      )
      .def_prop_rw("ixx", &LatParamStruct::ixx, &LatParamStruct::set_ixx, "Integer for general use")
      .def_prop_rw(
          "stable",
          &LatParamStruct::stable,
          &LatParamStruct::set_stable,
          "is closed lat stable?"
      )
      .def_prop_rw(
          "live_branch",
          &LatParamStruct::live_branch,
          &LatParamStruct::set_live_branch,
          "Should tracking be done on the branch?"
      )
      .def_prop_rw(
          "g1_integral",
          &LatParamStruct::g1_integral,
          &LatParamStruct::set_g1_integral,
          "Approximate |g| (bending strength) integral of branch."
      )
      .def_prop_rw(
          "g2_integral",
          &LatParamStruct::g2_integral,
          &LatParamStruct::set_g2_integral,
          "Approximate g^2 integral of branch."
      )
      .def_prop_rw(
          "g3_integral",
          &LatParamStruct::g3_integral,
          &LatParamStruct::set_g3_integral,
          "Approximate g^2 integral of branch."
      )
      .def_prop_rw(
          "bookkeeping_state",
          &LatParamStruct::bookkeeping_state,
          &LatParamStruct::set_bookkeeping_state,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Overall status for the branch."
      )
      .def_prop_rw(
          "beam_init",
          &LatParamStruct::beam_init,
          &LatParamStruct::set_beam_init,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const LatParamStruct &self, nb::dict &memo) { return LatParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatParamStruct &self, const LatParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D LatParamStruct arrays are not used in structs/routines
  // 2D LatParamStruct arrays are not used in structs/routines
  // 3D LatParamStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_pointer_struct
void init_lat_pointer_struct(nb::module_ &m, nb::class_<LatPointerStruct> &cls) {
  cls.def(
         "__init__",
         [](LatPointerStruct *self, const LatStruct *lat) {
           new (self) LatPointerStruct(ptr_to_opt_ref(lat));
         },
         nb::arg("lat") = nb::none()
  )
      .def_prop_rw(
          "lat",
          &LatPointerStruct::lat,
          &LatPointerStruct::set_lat,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatPointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatPointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const LatPointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatPointerStruct &self) {
            return LatPointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LatPointerStruct &self, nb::dict &memo) { return LatPointerStruct(self); }
      )
      .def(
          "__eq__",
          [](const LatPointerStruct &self, const LatPointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatPointerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<LatPointerStructArray1D, LatPointerStructAlloc1D>(
      m,
      "LatPointerStructArray1D",
      "LatPointerStructAlloc1D"
  );
  // 2D LatPointerStruct arrays are not used in structs/routines
  // 3D LatPointerStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_struct
void init_lat_struct(nb::module_ &m, nb::class_<LatStruct> &cls) {
  cls.def(
         "__init__",
         [](LatStruct *self,
            std::optional<std::string> use_name,
            std::optional<std::string> lattice,
            std::optional<std::string> machine,
            std::optional<std::string> input_file_name,
            std::optional<std::string> title,
            const ModeInfoStruct *a,
            const ModeInfoStruct *b,
            const ModeInfoStruct *z,
            const LatParamStruct *param,
            const BookkeepingStateStruct *lord_state,
            const EleStruct *ele_init,
            const CoordStruct *particle_start,
            const BeamInitStruct *beam_init,
            const PreTrackerStruct *pre_tracker,
            const NametableStruct *nametable,
            std::optional<std::vector<double>> custom,
            std::optional<int> version,
            std::optional<int> n_ele_track,
            std::optional<int> n_ele_max,
            std::optional<int> n_control_max,
            std::optional<int> n_ic_max,
            std::optional<int> input_taylor_order,
            std::optional<std::vector<int>> ic,
            std::optional<int> photon_type,
            std::optional<int> creation_hash,
            std::optional<int> ramper_slave_bookkeeping,
            std::optional<bool> parser_make_xfer_mats) {
           new (self) LatStruct(
               use_name,
               lattice,
               machine,
               input_file_name,
               title,
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(z),
               ptr_to_opt_ref(param),
               ptr_to_opt_ref(lord_state),
               ptr_to_opt_ref(ele_init),
               ptr_to_opt_ref(particle_start),
               ptr_to_opt_ref(beam_init),
               ptr_to_opt_ref(pre_tracker),
               ptr_to_opt_ref(nametable),
               custom,
               version,
               n_ele_track,
               n_ele_max,
               n_control_max,
               n_ic_max,
               input_taylor_order,
               ic,
               photon_type,
               creation_hash,
               ramper_slave_bookkeeping,
               parser_make_xfer_mats
           );
         },
         nb::arg("use_name") = nb::none(),
         nb::arg("lattice") = nb::none(),
         nb::arg("machine") = nb::none(),
         nb::arg("input_file_name") = nb::none(),
         nb::arg("title") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("z") = nb::none(),
         nb::arg("param") = nb::none(),
         nb::arg("lord_state") = nb::none(),
         nb::arg("ele_init") = nb::none(),
         nb::arg("particle_start") = nb::none(),
         nb::arg("beam_init") = nb::none(),
         nb::arg("pre_tracker") = nb::none(),
         nb::arg("nametable") = nb::none(),
         nb::arg("custom") = nb::none(),
         nb::arg("version") = nb::none(),
         nb::arg("n_ele_track") = nb::none(),
         nb::arg("n_ele_max") = nb::none(),
         nb::arg("n_control_max") = nb::none(),
         nb::arg("n_ic_max") = nb::none(),
         nb::arg("input_taylor_order") = nb::none(),
         nb::arg("ic") = nb::none(),
         nb::arg("photon_type") = nb::none(),
         nb::arg("creation_hash") = nb::none(),
         nb::arg("ramper_slave_bookkeeping") = nb::none(),
         nb::arg("parser_make_xfer_mats") = nb::none()
  )
      .def_prop_rw(
          "use_name",
          &LatStruct::use_name,
          &LatStruct::set_use_name,
          "Name of lat given by USE statement"
      )
      .def_prop_rw("lattice", &LatStruct::lattice, &LatStruct::set_lattice, "Lattice")
      .def_prop_rw(
          "machine",
          &LatStruct::machine,
          &LatStruct::set_machine,
          "Name of the machine the lattice is for ('LHC', etc)."
      )
      .def_prop_rw(
          "input_file_name",
          &LatStruct::input_file_name,
          &LatStruct::set_input_file_name,
          "Name of the lattice input file"
      )
      .def_prop_rw("title", &LatStruct::title, &LatStruct::set_title, "General title")
      .def_prop_ro(
          "print_str",
          &LatStruct::print_str,
          nb::keep_alive<0, 1>(),
          "Saved print statements."
      )
      .def_prop_ro(
          "constant",
          &LatStruct::constant,
          nb::keep_alive<0, 1>(),
          "Constants defined in the lattice"
      )
      .def_prop_rw(
          "a",
          &LatStruct::a,
          &LatStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Tunes (fractional part), etc."
      )
      .def_prop_rw(
          "b",
          &LatStruct::b,
          &LatStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Tunes (fractional part), etc."
      )
      .def_prop_rw(
          "z",
          &LatStruct::z,
          &LatStruct::set_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Tunes (fractional part), etc."
      )
      .def_prop_rw(
          "param",
          &LatStruct::param,
          &LatStruct::set_param,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Parameters"
      )
      .def_prop_rw(
          "lord_state",
          &LatStruct::lord_state,
          &LatStruct::set_lord_state,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "lord bookkeeping status."
      )
      .def_prop_rw(
          "ele_init",
          &LatStruct::ele_init,
          &LatStruct::set_ele_init,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For use by any program"
      )
      .def_prop_ro(
          "ele",
          &LatStruct::ele,
          nb::keep_alive<0, 1>(),
          "Array of elements [=> branch(0)]."
      )
      .def_prop_ro("branch", &LatStruct::branch, nb::keep_alive<0, 1>(), "Branch(0:) array")
      .def_prop_ro("control", &LatStruct::control, nb::keep_alive<0, 1>(), "Control list")
      .def_prop_rw(
          "particle_start",
          &LatStruct::particle_start,
          &LatStruct::set_particle_start,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Starting particle_coords."
      )
      .def_prop_rw(
          "beam_init",
          &LatStruct::beam_init,
          &LatStruct::set_beam_init,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Beam initialization."
      )
      .def_prop_rw(
          "pre_tracker",
          &LatStruct::pre_tracker,
          &LatStruct::set_pre_tracker,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For OPAL/IMPACT-T"
      )
      .def_prop_rw(
          "nametable",
          &LatStruct::nametable,
          &LatStruct::set_nametable,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For quick searching by element name."
      )
      .def_prop_rw(
          "custom",
          &LatStruct::custom,
          &LatStruct::set_custom,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Custom attributes."
      )
      .def_prop_rw("version", &LatStruct::version, &LatStruct::set_version, "Version number")
      .def_prop_rw(
          "n_ele_track",
          &LatStruct::n_ele_track,
          &LatStruct::set_n_ele_track,
          "Number of lat elements to track through."
      )
      .def_prop_rw(
          "n_ele_max",
          &LatStruct::n_ele_max,
          &LatStruct::set_n_ele_max,
          "Index of last valid element in %ele(:) array"
      )
      .def_prop_rw(
          "n_control_max",
          &LatStruct::n_control_max,
          &LatStruct::set_n_control_max,
          "Last index used in control_array"
      )
      .def_prop_rw(
          "n_ic_max",
          &LatStruct::n_ic_max,
          &LatStruct::set_n_ic_max,
          "Last index used in ic_array"
      )
      .def_prop_rw(
          "input_taylor_order",
          &LatStruct::input_taylor_order,
          &LatStruct::set_input_taylor_order,
          "As set in the input file"
      )
      .def_prop_rw(
          "ic",
          &LatStruct::ic,
          &LatStruct::set_ic,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Index to %control(:) from slaves."
      )
      .def_prop_rw(
          "photon_type",
          &LatStruct::photon_type,
          &LatStruct::set_photon_type,
          "Or coherent$. For X-ray simulations."
      )
      .def_prop_rw(
          "creation_hash",
          &LatStruct::creation_hash,
          &LatStruct::set_creation_hash,
          "Set by bmad_parser. creation_hash will vary if any of the lattice files are modified."
      )
      .def_prop_rw(
          "ramper_slave_bookkeeping",
          &LatStruct::ramper_slave_bookkeeping,
          &LatStruct::set_ramper_slave_bookkeeping
      )
      .def_prop_rw(
          "parser_make_xfer_mats",
          &LatStruct::parser_make_xfer_mats,
          &LatStruct::set_parser_make_xfer_mats,
          "Is Bmad parser to make element transfer matrices?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LatStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LatStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const LatStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatStruct &self) {
            return LatStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const LatStruct &self, nb::dict &memo) { return LatStruct(self); })
      .def(
          "__eq__",
          [](const LatStruct &self, const LatStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LatStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_linac_normal_mode_struct(nb::module_ &m, nb::class_<LinacNormalModeStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("i2_E4") = nb::none(),
         nb::arg("i3_E7") = nb::none(),
         nb::arg("i5a_E6") = nb::none(),
         nb::arg("i5b_E6") = nb::none(),
         nb::arg("sig_E1") = nb::none(),
         nb::arg("a_emittance_end") = nb::none(),
         nb::arg("b_emittance_end") = nb::none()
  )
      .def_prop_rw(
          "i2_E4",
          &LinacNormalModeStruct::i2_E4,
          &LinacNormalModeStruct::set_i2_E4,
          "Integral: g^2 * gamma^4"
      )
      .def_prop_rw(
          "i3_E7",
          &LinacNormalModeStruct::i3_E7,
          &LinacNormalModeStruct::set_i3_E7,
          "Integral: g^3 * gamma^7"
      )
      .def_prop_rw(
          "i5a_E6",
          &LinacNormalModeStruct::i5a_E6,
          &LinacNormalModeStruct::set_i5a_E6,
          "Integral: (g^3 * H_a) * gamma^6"
      )
      .def_prop_rw(
          "i5b_E6",
          &LinacNormalModeStruct::i5b_E6,
          &LinacNormalModeStruct::set_i5b_E6,
          "Integral: (g^3 * H_b) * gamma^6"
      )
      .def_prop_rw(
          "sig_E1",
          &LinacNormalModeStruct::sig_E1,
          &LinacNormalModeStruct::set_sig_E1,
          "Energy spread after 1 pass (eV)"
      )
      .def_prop_rw(
          "a_emittance_end",
          &LinacNormalModeStruct::a_emittance_end,
          &LinacNormalModeStruct::set_a_emittance_end,
          "a mode emittance at end of linac"
      )
      .def_prop_rw(
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
          [](const LinacNormalModeStruct &self, nb::dict &memo) {
            return LinacNormalModeStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const LinacNormalModeStruct &self, const LinacNormalModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LinacNormalModeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D LinacNormalModeStruct arrays are not used in structs/routines
  // 2D LinacNormalModeStruct arrays are not used in structs/routines
  // 3D LinacNormalModeStruct arrays are not used in structs/routines
}

// =============================================================================
// linear_ele_isf_struct
void init_linear_ele_isf_struct(nb::module_ &m, nb::class_<LinearEleIsfStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro(
          "node",
          &LinearEleIsfStruct::node,
          nb::keep_alive<0, 1>(),
          "Array per PTC integration node."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LinearEleIsfStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LinearEleIsfStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const LinearEleIsfStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LinearEleIsfStruct &self) {
            return LinearEleIsfStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LinearEleIsfStruct &self, nb::dict &memo) { return LinearEleIsfStruct(self); }
      )
      .def(
          "__eq__",
          [](const LinearEleIsfStruct &self, const LinearEleIsfStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LinearEleIsfStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<LinearEleIsfStructArray1D, LinearEleIsfStructAlloc1D>(
      m,
      "LinearEleIsfStructArray1D",
      "LinearEleIsfStructAlloc1D"
  );
  // 2D LinearEleIsfStruct arrays are not used in structs/routines
  // 3D LinearEleIsfStruct arrays are not used in structs/routines
}

// =============================================================================
// linear_isf1_struct
void init_linear_isf1_struct(nb::module_ &m, nb::class_<LinearIsf1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>>(),
         nb::arg("orb0") = nb::none(),
         nb::arg("isf") = nb::none(),
         nb::arg("s") = nb::none()
  )
      .def_prop_rw(
          "orb0",
          &LinearIsf1Struct::orb0,
          &LinearIsf1Struct::set_orb0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Closed orbit."
      )
      .def_prop_rw(
          "isf",
          &LinearIsf1Struct::isf,
          &LinearIsf1Struct::set_isf,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Linear ISF map at a given point."
      )
      .def_prop_rw(
          "s",
          &LinearIsf1Struct::s,
          &LinearIsf1Struct::set_s,
          "Offset from beginning of element. !! real(rp) :: m_1turn(6,6) = 0   ! Orbital 1-turn "
          "matrix."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return LinearIsf1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = LinearIsf1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const LinearIsf1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LinearIsf1Struct &self) {
            return LinearIsf1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const LinearIsf1Struct &self, nb::dict &memo) { return LinearIsf1Struct(self); }
      )
      .def(
          "__eq__",
          [](const LinearIsf1Struct &self, const LinearIsf1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const LinearIsf1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<LinearIsf1StructArray1D, LinearIsf1StructAlloc1D>(
      m,
      "LinearIsf1StructArray1D",
      "LinearIsf1StructAlloc1D"
  );
  // 2D LinearIsf1Struct arrays are not used in structs/routines
  // 3D LinearIsf1Struct arrays are not used in structs/routines
}

// =============================================================================
// layout
void init_layout(nb::module_ &m, nb::class_<Layout> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("NAME") = nb::none(),
         nb::arg("INDEX") = nb::none(),
         nb::arg("HARMONIC_NUMBER") = nb::none(),
         nb::arg("CLOSED") = nb::none(),
         nb::arg("N") = nb::none(),
         nb::arg("NTHIN") = nb::none(),
         nb::arg("THIN") = nb::none(),
         nb::arg("LASTPOS") = nb::none()
  )
      .def_prop_rw("NAME", &Layout::NAME, &Layout::set_NAME, "IDENTIFICATION")
      .def_prop_rw("INDEX", &Layout::INDEX, &Layout::set_INDEX, "IDENTIFICATION, CHARGE SIGN")
      .def_prop_rw("HARMONIC_NUMBER", &Layout::HARMONIC_NUMBER, &Layout::set_HARMONIC_NUMBER)
      .def_prop_rw("CLOSED", &Layout::CLOSED, &Layout::set_CLOSED)
      .def_prop_rw("N", &Layout::N, &Layout::set_N, "TOTAL ELEMENT IN THE CHAIN")
      .def_prop_rw(
          "NTHIN",
          &Layout::NTHIN,
          &Layout::set_NTHIN,
          "NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)"
      )
      .def_prop_rw(
          "THIN",
          &Layout::THIN,
          &Layout::set_THIN,
          "PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS POINTERS OF LINK LAYOUT"
      )
      .def_prop_rw("LASTPOS", &Layout::LASTPOS, &Layout::set_LASTPOS, "POSITION OF LAST VISITED")

      .def("__repr__", [](const Layout &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Layout &self) {
            return Layout(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const Layout &self, nb::dict &memo) { return Layout(self); })
      .def(
          "__eq__",
          [](const Layout &self, const Layout &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Layout &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D Layout arrays are not used in structs/routines
  // 2D Layout arrays are not used in structs/routines
  // 3D Layout arrays are not used in structs/routines
}