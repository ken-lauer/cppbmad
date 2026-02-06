#include "pybmad/generated/structs_l.hpp"

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
      // LatEleLocStruct.ix_ele (0D_NOT_integer -
      .def_property("ix_ele", &LatEleLocStruct::ix_ele, &LatEleLocStruct::set_ix_ele)
      // LatEleLocStruct.ix_branch (0D_NOT_integer -
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

      ;

  bind_FTypeArrayND<LatEleLocStructArray1D>(m, "LatEleLocStructArray1D");
  bind_FTypeAlloc1D<LatEleLocStructAlloc1D>(m, "LatEleLocStructAlloc1D");
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
      // LatEleOrder1Struct.ix_branch (0D_NOT_integer - Branch index
      .def_property("ix_branch", &LatEleOrder1Struct::ix_branch, &LatEleOrder1Struct::set_ix_branch)
      // LatEleOrder1Struct.ix_order (0D_NOT_integer - Order index. -1 -> Unique in lattice, 0 ->
      // unique in branch.
      .def_property("ix_order", &LatEleOrder1Struct::ix_order, &LatEleOrder1Struct::set_ix_order)
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

      ;

  bind_FTypeArrayND<LatEleOrder1StructArray1D>(m, "LatEleOrder1StructArray1D");
  bind_FTypeAlloc1D<LatEleOrder1StructAlloc1D>(m, "LatEleOrder1StructAlloc1D");
  // 2D LatEleOrder1Struct arrays are not used in structs/routines
  // 3D LatEleOrder1Struct arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_array_struct
void init_lat_ele_order_array_struct(py::module &m, py::class_<LatEleOrderArrayStruct> &cls) {
  cls.def(py::init<>())
      // LatEleOrderArrayStruct.ele (1D_ALLOC_type -
      .def_property_readonly("ele", &LatEleOrderArrayStruct::ele)
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

      ;

  bind_FTypeArrayND<LatEleOrderArrayStructArray1D>(m, "LatEleOrderArrayStructArray1D");
  bind_FTypeAlloc1D<LatEleOrderArrayStructAlloc1D>(m, "LatEleOrderArrayStructAlloc1D");
  // 2D LatEleOrderArrayStruct arrays are not used in structs/routines
  // 3D LatEleOrderArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_struct
void init_lat_ele_order_struct(py::module &m, py::class_<LatEleOrderStruct> &cls) {
  cls.def(py::init<>())
      // LatEleOrderStruct.branch (1D_ALLOC_type -
      .def_property_readonly("branch", &LatEleOrderStruct::branch)

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
             optional_ref<const std::vector<std::vector<double>>>,
             optional_ref<const std::vector<std::vector<double>>>,
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
      // LatParamStruct.n_part (0D_NOT_real - Particles/bunch (for BeamBeam elements).
      .def_property("n_part", &LatParamStruct::n_part, &LatParamStruct::set_n_part)
      // LatParamStruct.total_length (0D_NOT_real - total_length of branch. Warning: branch may not
      // start at s = 0.
      .def_property(
          "total_length",
          &LatParamStruct::total_length,
          &LatParamStruct::set_total_length
      )
      // LatParamStruct.unstable_factor (0D_NOT_real - If positive: Growth rate/turn if unstable in
      // closed branches or |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
      .def_property(
          "unstable_factor",
          &LatParamStruct::unstable_factor,
          &LatParamStruct::set_unstable_factor
      )
      // LatParamStruct.t1_with_RF (2D_NOT_real - Full 1-turn matrix with RF on.
      .def_property("t1_with_RF", &LatParamStruct::t1_with_RF, &LatParamStruct::set_t1_with_RF)
      // LatParamStruct.t1_no_RF (2D_NOT_real - Full 1-turn matrix with RF off.
      .def_property("t1_no_RF", &LatParamStruct::t1_no_RF, &LatParamStruct::set_t1_no_RF)
      // LatParamStruct.spin_tune (0D_NOT_real - Closed orbit spin tune.
      .def_property("spin_tune", &LatParamStruct::spin_tune, &LatParamStruct::set_spin_tune)
      // LatParamStruct.particle (0D_NOT_integer - Reference particle: positron$, electron$, etc.
      // Call lattice_bookkeeper if this is changed.
      .def_property("particle", &LatParamStruct::particle, &LatParamStruct::set_particle)
      // LatParamStruct.default_tracking_species (0D_NOT_integer - Default particle type to use in
      // tracking.
      .def_property(
          "default_tracking_species",
          &LatParamStruct::default_tracking_species,
          &LatParamStruct::set_default_tracking_species
      )
      // LatParamStruct.geometry (0D_NOT_integer - open$ or closed$
      .def_property("geometry", &LatParamStruct::geometry, &LatParamStruct::set_geometry)
      // LatParamStruct.ixx (0D_NOT_integer - Integer for general use
      .def_property("ixx", &LatParamStruct::ixx, &LatParamStruct::set_ixx)
      // LatParamStruct.stable (0D_NOT_logical - is closed lat stable?
      .def_property("stable", &LatParamStruct::stable, &LatParamStruct::set_stable)
      // LatParamStruct.live_branch (0D_NOT_logical - Should tracking be done on the branch?
      .def_property("live_branch", &LatParamStruct::live_branch, &LatParamStruct::set_live_branch)
      // LatParamStruct.g1_integral (0D_NOT_real - Approximate |g| (bending strength) integral of
      // branch.
      .def_property("g1_integral", &LatParamStruct::g1_integral, &LatParamStruct::set_g1_integral)
      // LatParamStruct.g2_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property("g2_integral", &LatParamStruct::g2_integral, &LatParamStruct::set_g2_integral)
      // LatParamStruct.g3_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property("g3_integral", &LatParamStruct::g3_integral, &LatParamStruct::set_g3_integral)
      // LatParamStruct.bookkeeping_state (0D_NOT_type - Overall status for the branch.
      .def_property(
          "bookkeeping_state",
          &LatParamStruct::bookkeeping_state,
          &LatParamStruct::set_bookkeeping_state
      )
      // LatParamStruct.beam_init (0D_NOT_type - For beam initialization.
      .def_property("beam_init", &LatParamStruct::beam_init, &LatParamStruct::set_beam_init)

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
             optional_ref<const std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             optional_ref<const std::vector<int>>,
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
      // LatStruct.use_name (0D_NOT_character - Name of lat given by USE statement
      .def_property("use_name", &LatStruct::use_name, &LatStruct::set_use_name)
      // LatStruct.lattice (0D_NOT_character - Lattice
      .def_property("lattice", &LatStruct::lattice, &LatStruct::set_lattice)
      // LatStruct.machine (0D_NOT_character - Name of the machine the lattice is for ('LHC', etc).
      .def_property("machine", &LatStruct::machine, &LatStruct::set_machine)
      // LatStruct.input_file_name (0D_NOT_character - Name of the lattice input file
      .def_property("input_file_name", &LatStruct::input_file_name, &LatStruct::set_input_file_name)
      // LatStruct.title (0D_NOT_character - General title
      .def_property("title", &LatStruct::title, &LatStruct::set_title)
      // LatStruct.print_str (1D_ALLOC_character - Saved print statements.
      .def_property_readonly("print_str", &LatStruct::print_str)
      // LatStruct.constant (1D_ALLOC_type - Constants defined in the lattice
      .def_property_readonly("constant", &LatStruct::constant)
      // LatStruct.a (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("a", &LatStruct::a, &LatStruct::set_a)
      // LatStruct.b (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("b", &LatStruct::b, &LatStruct::set_b)
      // LatStruct.z (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("z", &LatStruct::z, &LatStruct::set_z)
      // LatStruct.param (0D_PTR_type - Parameters
      .def_property("param", &LatStruct::param, &LatStruct::set_param)
      // LatStruct.lord_state (0D_NOT_type - lord bookkeeping status.
      .def_property("lord_state", &LatStruct::lord_state, &LatStruct::set_lord_state)
      // LatStruct.ele_init (0D_NOT_type - For use by any program
      .def_property("ele_init", &LatStruct::ele_init, &LatStruct::set_ele_init)
      // LatStruct.ele (1D_PTR_type - Array of elements [=> branch(0)].
      .def_property_readonly("ele", &LatStruct::ele)
      // LatStruct.branch (1D_ALLOC_type - Branch(0:) array
      .def_property_readonly("branch", &LatStruct::branch)
      // LatStruct.control (1D_ALLOC_type - Control list
      .def_property_readonly("control", &LatStruct::control)
      // LatStruct.particle_start (0D_PTR_type - Starting particle_coords.
      .def_property("particle_start", &LatStruct::particle_start, &LatStruct::set_particle_start)
      // LatStruct.beam_init (0D_NOT_type - Beam initialization.
      .def_property("beam_init", &LatStruct::beam_init, &LatStruct::set_beam_init)
      // LatStruct.pre_tracker (0D_NOT_type - For OPAL/IMPACT-T
      .def_property("pre_tracker", &LatStruct::pre_tracker, &LatStruct::set_pre_tracker)
      // LatStruct.custom (1D_ALLOC_real - Custom attributes.
      .def_property("custom", &LatStruct::custom, &LatStruct::set_custom)
      // LatStruct.version (0D_NOT_integer - Version number
      .def_property("version", &LatStruct::version, &LatStruct::set_version)
      // LatStruct.n_ele_track (0D_PTR_integer - Number of lat elements to track through.
      .def_property("n_ele_track", &LatStruct::n_ele_track, &LatStruct::set_n_ele_track)
      // LatStruct.n_ele_max (0D_PTR_integer - Index of last valid element in %ele(:) array
      .def_property("n_ele_max", &LatStruct::n_ele_max, &LatStruct::set_n_ele_max)
      // LatStruct.n_control_max (0D_NOT_integer - Last index used in control_array
      .def_property("n_control_max", &LatStruct::n_control_max, &LatStruct::set_n_control_max)
      // LatStruct.n_ic_max (0D_NOT_integer - Last index used in ic_array
      .def_property("n_ic_max", &LatStruct::n_ic_max, &LatStruct::set_n_ic_max)
      // LatStruct.input_taylor_order (0D_NOT_integer - As set in the input file
      .def_property(
          "input_taylor_order",
          &LatStruct::input_taylor_order,
          &LatStruct::set_input_taylor_order
      )
      // LatStruct.ic (1D_ALLOC_integer - Index to %control(:) from slaves.
      .def_property("ic", &LatStruct::ic, &LatStruct::set_ic)
      // LatStruct.photon_type (0D_NOT_integer - Or coherent$. For X-ray simulations.
      .def_property("photon_type", &LatStruct::photon_type, &LatStruct::set_photon_type)
      // LatStruct.creation_hash (0D_NOT_integer - Set by bmad_parser. creation_hash will vary if
      // any of the lattice files are modified.
      .def_property("creation_hash", &LatStruct::creation_hash, &LatStruct::set_creation_hash)
      // LatStruct.ramper_slave_bookkeeping (0D_NOT_integer -
      .def_property(
          "ramper_slave_bookkeeping",
          &LatStruct::ramper_slave_bookkeeping,
          &LatStruct::set_ramper_slave_bookkeeping
      )
      // LatStruct.parser_make_xfer_mats (0D_NOT_logical - Is Bmad parser to make element transfer
      // matrices?
      .def_property(
          "parser_make_xfer_mats",
          &LatStruct::parser_make_xfer_mats,
          &LatStruct::set_parser_make_xfer_mats
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

      ;

  bind_FTypeArrayND<LatStructArray1D>(m, "LatStructArray1D");
  bind_FTypeAlloc1D<LatStructAlloc1D>(m, "LatStructAlloc1D");
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
      // LinacNormalModeStruct.i2_E4 (0D_NOT_real - Integral: g^2 * gamma^4
      .def_property("i2_E4", &LinacNormalModeStruct::i2_E4, &LinacNormalModeStruct::set_i2_E4)
      // LinacNormalModeStruct.i3_E7 (0D_NOT_real - Integral: g^3 * gamma^7
      .def_property("i3_E7", &LinacNormalModeStruct::i3_E7, &LinacNormalModeStruct::set_i3_E7)
      // LinacNormalModeStruct.i5a_E6 (0D_NOT_real - Integral: (g^3 * H_a) * gamma^6
      .def_property("i5a_E6", &LinacNormalModeStruct::i5a_E6, &LinacNormalModeStruct::set_i5a_E6)
      // LinacNormalModeStruct.i5b_E6 (0D_NOT_real - Integral: (g^3 * H_b) * gamma^6
      .def_property("i5b_E6", &LinacNormalModeStruct::i5b_E6, &LinacNormalModeStruct::set_i5b_E6)
      // LinacNormalModeStruct.sig_E1 (0D_NOT_real - Energy spread after 1 pass (eV)
      .def_property("sig_E1", &LinacNormalModeStruct::sig_E1, &LinacNormalModeStruct::set_sig_E1)
      // LinacNormalModeStruct.a_emittance_end (0D_NOT_real - a mode emittance at end of linac
      .def_property(
          "a_emittance_end",
          &LinacNormalModeStruct::a_emittance_end,
          &LinacNormalModeStruct::set_a_emittance_end
      )
      // LinacNormalModeStruct.b_emittance_end (0D_NOT_real - b mode emittance at end of linac
      .def_property(
          "b_emittance_end",
          &LinacNormalModeStruct::b_emittance_end,
          &LinacNormalModeStruct::set_b_emittance_end
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
      // Layout.NAME (0D_PTR_character - IDENTIFICATION
      .def_property("NAME", &Layout::NAME, &Layout::set_NAME)
      // Layout.INDEX (0D_PTR_integer - IDENTIFICATION, CHARGE SIGN
      .def_property("INDEX", &Layout::INDEX, &Layout::set_INDEX)
      // Layout.HARMONIC_NUMBER (0D_PTR_real -
      .def_property("HARMONIC_NUMBER", &Layout::HARMONIC_NUMBER, &Layout::set_HARMONIC_NUMBER)
      // Layout.CLOSED (0D_PTR_logical -
      .def_property("CLOSED", &Layout::CLOSED, &Layout::set_CLOSED)
      // Layout.N (0D_PTR_integer - TOTAL ELEMENT IN THE CHAIN
      .def_property("N", &Layout::N, &Layout::set_N)
      // Layout.NTHIN (0D_PTR_integer - NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)
      .def_property("NTHIN", &Layout::NTHIN, &Layout::set_NTHIN)
      // Layout.THIN (0D_PTR_real - PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS POINTERS OF
      // LINK LAYOUT
      .def_property("THIN", &Layout::THIN, &Layout::set_THIN)
      // Layout.LASTPOS (0D_PTR_integer - POSITION OF LAST VISITED
      .def_property("LASTPOS", &Layout::LASTPOS, &Layout::set_LASTPOS)
      // Layout.LAST (0D_PTR_type - LAST VISITED
      .def_property("LAST", &Layout::LAST, &Layout::set_LAST)
      // Layout.END (0D_PTR_type -
      .def_property("END", &Layout::END, &Layout::set_END)
      // Layout.START (0D_PTR_type -
      .def_property("START", &Layout::START, &Layout::set_START)
      // Layout.START_GROUND (0D_PTR_type - STORE THE GROUNDED VALUE OF START DURING CIRCULAR
      // SCANNING
      .def_property("START_GROUND", &Layout::START_GROUND, &Layout::set_START_GROUND)
      // Layout.END_GROUND (0D_PTR_type - STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
      .def_property("END_GROUND", &Layout::END_GROUND, &Layout::set_END_GROUND)
      // Layout.NEXT (0D_PTR_type -
      .def_property("NEXT", &Layout::NEXT, &Layout::set_NEXT)
      // Layout.PREVIOUS (0D_PTR_type -
      .def_property("PREVIOUS", &Layout::PREVIOUS, &Layout::set_PREVIOUS)

      .def("__repr__", [](const Layout &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Layout &self) {
            return Layout(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const Layout &self, py::dict &memo) { return Layout(self); })

      ;

  // 1D Layout arrays are not used in structs/routines
  // 2D Layout arrays are not used in structs/routines
  // 3D Layout arrays are not used in structs/routines
}