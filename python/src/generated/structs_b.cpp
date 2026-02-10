#include "pybmad/generated/structs_b.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// bbu_beam_struct
void init_bbu_beam_struct(py::module &m, py::class_<BbuBeamStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<int>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("ix_ele_bunch") = py::none(),
         py::arg("ix_bunch_head") = py::none(),
         py::arg("ix_bunch_end") = py::none(),
         py::arg("n_bunch_in_lat") = py::none(),
         py::arg("ix_stage_voltage_max") = py::none(),
         py::arg("hom_voltage_max") = py::none(),
         py::arg("time_now") = py::none(),
         py::arg("one_turn_time") = py::none(),
         py::arg("rf_wavelength_max") = py::none()
  )
      .def_property_readonly("bunch", &BbuBeamStruct::bunch, "Bunches in the lattice")
      .def_property_readonly("stage", &BbuBeamStruct::stage)
      .def_property(
          "ix_ele_bunch",
          &BbuBeamStruct::ix_ele_bunch,
          &BbuBeamStruct::set_ix_ele_bunch,
          "element where bunch is"
      )
      .def_property(
          "ix_bunch_head",
          &BbuBeamStruct::ix_bunch_head,
          &BbuBeamStruct::set_ix_bunch_head,
          "Index to head bunch(:)"
      )
      .def_property(
          "ix_bunch_end",
          &BbuBeamStruct::ix_bunch_end,
          &BbuBeamStruct::set_ix_bunch_end,
          "Index of the end bunch(:). -1 -> no bunches."
      )
      .def_property(
          "n_bunch_in_lat",
          &BbuBeamStruct::n_bunch_in_lat,
          &BbuBeamStruct::set_n_bunch_in_lat,
          "Number of bunches transversing the lattice."
      )
      .def_property(
          "ix_stage_voltage_max",
          &BbuBeamStruct::ix_stage_voltage_max,
          &BbuBeamStruct::set_ix_stage_voltage_max
      )
      .def_property(
          "hom_voltage_max",
          &BbuBeamStruct::hom_voltage_max,
          &BbuBeamStruct::set_hom_voltage_max
      )
      .def_property("time_now", &BbuBeamStruct::time_now, &BbuBeamStruct::set_time_now)
      .def_property(
          "one_turn_time",
          &BbuBeamStruct::one_turn_time,
          &BbuBeamStruct::set_one_turn_time
      )
      .def_property(
          "rf_wavelength_max",
          &BbuBeamStruct::rf_wavelength_max,
          &BbuBeamStruct::set_rf_wavelength_max
      )

      .def("__repr__", [](const BbuBeamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuBeamStruct &self) {
            return BbuBeamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BbuBeamStruct &self, py::dict &memo) { return BbuBeamStruct(self); }
      )

      ;

  // 1D BbuBeamStruct arrays are not used in structs/routines
  // 2D BbuBeamStruct arrays are not used in structs/routines
  // 3D BbuBeamStruct arrays are not used in structs/routines
}

// =============================================================================
// bbu_param_struct
void init_bbu_param_struct(py::module &m, py::class_<BbuParamStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("lat_filename") = py::none(),
         py::arg("lat2_filename") = py::none(),
         py::arg("bunch_by_bunch_info_file") = py::none(),
         py::arg("hybridize") = py::none(),
         py::arg("write_digested_hybrid_lat") = py::none(),
         py::arg("write_voltage_vs_time_dat") = py::none(),
         py::arg("keep_overlays_and_groups") = py::none(),
         py::arg("keep_all_lcavities") = py::none(),
         py::arg("use_taylor_for_hybrids") = py::none(),
         py::arg("stable_orbit_anal") = py::none(),
         py::arg("limit_factor") = py::none(),
         py::arg("simulation_turns_max") = py::none(),
         py::arg("bunch_freq") = py::none(),
         py::arg("init_particle_offset") = py::none(),
         py::arg("current") = py::none(),
         py::arg("rel_tol") = py::none(),
         py::arg("drscan") = py::none(),
         py::arg("use_interpolated_threshold") = py::none(),
         py::arg("write_hom_info") = py::none(),
         py::arg("elindex") = py::none(),
         py::arg("elname") = py::none(),
         py::arg("nstep") = py::none(),
         py::arg("begdr") = py::none(),
         py::arg("enddr") = py::none(),
         py::arg("nrep") = py::none(),
         py::arg("ran_seed") = py::none(),
         py::arg("hom_order_cutoff") = py::none(),
         py::arg("ran_gauss_sigma_cut") = py::none(),
         py::arg("ele_track_end") = py::none(),
         py::arg("ix_ele_track_end") = py::none(),
         py::arg("regression") = py::none(),
         py::arg("normalize_z_to_rf") = py::none(),
         py::arg("ramp_on") = py::none(),
         py::arg("ramp_pattern") = py::none(),
         py::arg("ramp_n_start") = py::none(),
         py::arg("n_ramp_pattern") = py::none()
  )
      .def_property(
          "lat_filename",
          &BbuParamStruct::lat_filename,
          &BbuParamStruct::set_lat_filename,
          "Bmad lattice file name"
      )
      .def_property(
          "lat2_filename",
          &BbuParamStruct::lat2_filename,
          &BbuParamStruct::set_lat2_filename,
          "Bmad lattice2 file name for secondary parser"
      )
      .def_property(
          "bunch_by_bunch_info_file",
          &BbuParamStruct::bunch_by_bunch_info_file,
          &BbuParamStruct::set_bunch_by_bunch_info_file,
          "For outputting bunch-by-bunch info."
      )
      .def_property(
          "hybridize",
          &BbuParamStruct::hybridize,
          &BbuParamStruct::set_hybridize,
          "Combine non-hom elements to speed up simulation?"
      )
      .def_property(
          "write_digested_hybrid_lat",
          &BbuParamStruct::write_digested_hybrid_lat,
          &BbuParamStruct::set_write_digested_hybrid_lat,
          "For debugging purposes."
      )
      .def_property(
          "write_voltage_vs_time_dat",
          &BbuParamStruct::write_voltage_vs_time_dat,
          &BbuParamStruct::set_write_voltage_vs_time_dat,
          "For debugging purposes."
      )
      .def_property(
          "keep_overlays_and_groups",
          &BbuParamStruct::keep_overlays_and_groups,
          &BbuParamStruct::set_keep_overlays_and_groups,
          "Keep when hybridizing?"
      )
      .def_property(
          "keep_all_lcavities",
          &BbuParamStruct::keep_all_lcavities,
          &BbuParamStruct::set_keep_all_lcavities,
          "Keep when hybridizing?"
      )
      .def_property(
          "use_taylor_for_hybrids",
          &BbuParamStruct::use_taylor_for_hybrids,
          &BbuParamStruct::set_use_taylor_for_hybrids,
          "Use taylor map for hybrids when true. Otherwise tracking method is linear."
      )
      .def_property(
          "stable_orbit_anal",
          &BbuParamStruct::stable_orbit_anal,
          &BbuParamStruct::set_stable_orbit_anal,
          "Write stable_orbit.out and hom_voltage.out?"
      )
      .def_property(
          "limit_factor",
          &BbuParamStruct::limit_factor,
          &BbuParamStruct::set_limit_factor,
          "Init_hom_amp * limit_factor = simulation unstable limit"
      )
      .def_property(
          "simulation_turns_max",
          &BbuParamStruct::simulation_turns_max,
          &BbuParamStruct::set_simulation_turns_max,
          "Sets the duration of the simulation."
      )
      .def_property(
          "bunch_freq",
          &BbuParamStruct::bunch_freq,
          &BbuParamStruct::set_bunch_freq,
          "Freq in Hz."
      )
      .def_property(
          "init_particle_offset",
          &BbuParamStruct::init_particle_offset,
          &BbuParamStruct::set_init_particle_offset,
          "Initial particle offset for particles born in the first turn period."
      )
      .def_property(
          "current",
          &BbuParamStruct::current,
          &BbuParamStruct::set_current,
          "Starting current (amps)"
      )
      .def_property(
          "rel_tol",
          &BbuParamStruct::rel_tol,
          &BbuParamStruct::set_rel_tol,
          "Final threshold current accuracy."
      )
      .def_property(
          "drscan",
          &BbuParamStruct::drscan,
          &BbuParamStruct::set_drscan,
          "If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3."
      )
      .def_property(
          "use_interpolated_threshold",
          &BbuParamStruct::use_interpolated_threshold,
          &BbuParamStruct::set_use_interpolated_threshold
      )
      .def_property(
          "write_hom_info",
          &BbuParamStruct::write_hom_info,
          &BbuParamStruct::set_write_hom_info,
          "Write HOM parameters to main output file?"
      )
      .def_property("elindex", &BbuParamStruct::elindex, &BbuParamStruct::set_elindex)
      .def_property(
          "elname",
          &BbuParamStruct::elname,
          &BbuParamStruct::set_elname,
          "Element to step length for DRSCAN"
      )
      .def_property(
          "nstep",
          &BbuParamStruct::nstep,
          &BbuParamStruct::set_nstep,
          "Number of steps for DRSCAN."
      )
      .def_property(
          "begdr",
          &BbuParamStruct::begdr,
          &BbuParamStruct::set_begdr,
          "Beginning DR value for DRSCAN."
      )
      .def_property(
          "enddr",
          &BbuParamStruct::enddr,
          &BbuParamStruct::set_enddr,
          "End DR value for DRSCAN."
      )
      .def_property(
          "nrep",
          &BbuParamStruct::nrep,
          &BbuParamStruct::set_nrep,
          "Number of times to repeat threshold calculation"
      )
      .def_property(
          "ran_seed",
          &BbuParamStruct::ran_seed,
          &BbuParamStruct::set_ran_seed,
          "If set to 0, the output results will vary from run to run."
      )
      .def_property(
          "hom_order_cutoff",
          &BbuParamStruct::hom_order_cutoff,
          &BbuParamStruct::set_hom_order_cutoff,
          "If positive -> ignore HOM's with order greater than this."
      )
      .def_property(
          "ran_gauss_sigma_cut",
          &BbuParamStruct::ran_gauss_sigma_cut,
          &BbuParamStruct::set_ran_gauss_sigma_cut
      )
      .def_property(
          "ele_track_end",
          &BbuParamStruct::ele_track_end,
          &BbuParamStruct::set_ele_track_end
      )
      .def_property(
          "ix_ele_track_end",
          &BbuParamStruct::ix_ele_track_end,
          &BbuParamStruct::set_ix_ele_track_end,
          "Default: set to last element with a wake"
      )
      .def_property(
          "regression",
          &BbuParamStruct::regression,
          &BbuParamStruct::set_regression,
          "Do regression test?"
      )
      .def_property(
          "normalize_z_to_rf",
          &BbuParamStruct::normalize_z_to_rf,
          &BbuParamStruct::set_normalize_z_to_rf,
          "make starting z = mod(z, rf_wavelength)? Ramp parameters"
      )
      .def_property("ramp_on", &BbuParamStruct::ramp_on, &BbuParamStruct::set_ramp_on)
      .def_property(
          "ramp_pattern",
          &BbuParamStruct::ramp_pattern,
          &BbuParamStruct::set_ramp_pattern
      )
      .def_property(
          "ramp_n_start",
          &BbuParamStruct::ramp_n_start,
          &BbuParamStruct::set_ramp_n_start,
          "Index of start of ramp Internal parameters"
      )
      .def_property(
          "n_ramp_pattern",
          &BbuParamStruct::n_ramp_pattern,
          &BbuParamStruct::set_n_ramp_pattern,
          "Number of valid ramp_pattern"
      )

      .def("__repr__", [](const BbuParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuParamStruct &self) {
            return BbuParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BbuParamStruct &self, py::dict &memo) { return BbuParamStruct(self); }
      )

      ;

  // 1D BbuParamStruct arrays are not used in structs/routines
  // 2D BbuParamStruct arrays are not used in structs/routines
  // 3D BbuParamStruct arrays are not used in structs/routines
}

// =============================================================================
// bbu_stage_struct
void init_bbu_stage_struct(py::module &m, py::class_<BbuStageStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int>>(),
         py::arg("ix_ele_lr_wake") = py::none(),
         py::arg("ix_ele_stage_end") = py::none(),
         py::arg("ix_pass") = py::none(),
         py::arg("ix_stage_pass1") = py::none(),
         py::arg("ix_head_bunch") = py::none(),
         py::arg("ix_hom_max") = py::none(),
         py::arg("hom_voltage_max") = py::none(),
         py::arg("time_at_wake_ele") = py::none(),
         py::arg("ave_orb") = py::none(),
         py::arg("rms_orb") = py::none(),
         py::arg("min_orb") = py::none(),
         py::arg("max_orb") = py::none(),
         py::arg("n_orb") = py::none()
  )
      .def_property(
          "ix_ele_lr_wake",
          &BbuStageStruct::ix_ele_lr_wake,
          &BbuStageStruct::set_ix_ele_lr_wake,
          "Element index of element with the wake"
      )
      .def_property(
          "ix_ele_stage_end",
          &BbuStageStruct::ix_ele_stage_end,
          &BbuStageStruct::set_ix_ele_stage_end,
          "Element at end of stage."
      )
      .def_property(
          "ix_pass",
          &BbuStageStruct::ix_pass,
          &BbuStageStruct::set_ix_pass,
          "Pass index when in multipass section"
      )
      .def_property(
          "ix_stage_pass1",
          &BbuStageStruct::ix_stage_pass1,
          &BbuStageStruct::set_ix_stage_pass1,
          "Index of corresponding stage on first pass"
      )
      .def_property(
          "ix_head_bunch",
          &BbuStageStruct::ix_head_bunch,
          &BbuStageStruct::set_ix_head_bunch
      )
      .def_property("ix_hom_max", &BbuStageStruct::ix_hom_max, &BbuStageStruct::set_ix_hom_max)
      .def_property(
          "hom_voltage_max",
          &BbuStageStruct::hom_voltage_max,
          &BbuStageStruct::set_hom_voltage_max
      )
      .def_property(
          "time_at_wake_ele",
          &BbuStageStruct::time_at_wake_ele,
          &BbuStageStruct::set_time_at_wake_ele
      )
      .def_property("ave_orb", &BbuStageStruct::ave_orb, &BbuStageStruct::set_ave_orb)
      .def_property("rms_orb", &BbuStageStruct::rms_orb, &BbuStageStruct::set_rms_orb)
      .def_property("min_orb", &BbuStageStruct::min_orb, &BbuStageStruct::set_min_orb)
      .def_property("max_orb", &BbuStageStruct::max_orb, &BbuStageStruct::set_max_orb)
      .def_property("n_orb", &BbuStageStruct::n_orb, &BbuStageStruct::set_n_orb)
      .def_static(
          "new_array1d",
          [](int sz) { return BbuStageStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BbuStageStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const BbuStageStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuStageStruct &self) {
            return BbuStageStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BbuStageStruct &self, py::dict &memo) { return BbuStageStruct(self); }
      )

      ;

  bind_1d_type_array_pair<BbuStageStructArray1D, BbuStageStructAlloc1D>(
      m,
      "BbuStageStructArray1D",
      "BbuStageStructAlloc1D"
  );
  // 2D BbuStageStruct arrays are not used in structs/routines
  // 3D BbuStageStruct arrays are not used in structs/routines
}

// =============================================================================
// beam_init_struct
void init_beam_init_struct(py::module &m, py::class_<BeamInitStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::vector<double>>,
             optional_ref<const KvBeamInitStruct>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>>(),
         py::arg("position_file") = py::none(),
         py::arg("spin") = py::none(),
         py::arg("KV") = py::none(),
         py::arg("center_jitter") = py::none(),
         py::arg("emit_jitter") = py::none(),
         py::arg("sig_z_jitter") = py::none(),
         py::arg("sig_pz_jitter") = py::none(),
         py::arg("n_particle") = py::none(),
         py::arg("renorm_center") = py::none(),
         py::arg("renorm_sigma") = py::none(),
         py::arg("random_engine") = py::none(),
         py::arg("random_gauss_converter") = py::none(),
         py::arg("random_sigma_cutoff") = py::none(),
         py::arg("a_norm_emit") = py::none(),
         py::arg("b_norm_emit") = py::none(),
         py::arg("a_emit") = py::none(),
         py::arg("b_emit") = py::none(),
         py::arg("dPz_dz") = py::none(),
         py::arg("center") = py::none(),
         py::arg("t_offset") = py::none(),
         py::arg("dt_bunch") = py::none(),
         py::arg("sig_z") = py::none(),
         py::arg("sig_pz") = py::none(),
         py::arg("bunch_charge") = py::none(),
         py::arg("n_bunch") = py::none(),
         py::arg("ix_turn") = py::none(),
         py::arg("species") = py::none(),
         py::arg("full_6D_coupling_calc") = py::none(),
         py::arg("use_particle_start") = py::none(),
         py::arg("use_t_coords") = py::none(),
         py::arg("use_z_as_t") = py::none(),
         py::arg("file_name") = py::none()
  )
      .def_property(
          "position_file",
          &BeamInitStruct::position_file,
          &BeamInitStruct::set_position_file,
          "File with particle positions."
      )
      .def_property_readonly(
          "distribution_type",
          &BeamInitStruct::distribution_type,
          "distribution type (in x-px, y-py, and z-pz planes) 'ELLIPSE', 'KV', 'GRID', 'FILE', "
          "'RAN_GAUSS' or '' = 'RAN_GAUSS'"
      )
      .def_property("spin", &BeamInitStruct::spin, &BeamInitStruct::set_spin, "Spin (x, y, z)")
      .def_property_readonly("ellipse", &BeamInitStruct::ellipse, "Ellipse beam distribution")
      .def_property("KV", &BeamInitStruct::KV, &BeamInitStruct::set_KV, "KV beam distribution")
      .def_property_readonly("grid", &BeamInitStruct::grid, "Grid beam distribution")
      .def_property(
          "center_jitter",
          &BeamInitStruct::center_jitter,
          &BeamInitStruct::set_center_jitter,
          "Bunch center rms jitter"
      )
      .def_property(
          "emit_jitter",
          &BeamInitStruct::emit_jitter,
          &BeamInitStruct::set_emit_jitter,
          "a and b bunch emittance rms jitter normalized to emittance"
      )
      .def_property(
          "sig_z_jitter",
          &BeamInitStruct::sig_z_jitter,
          &BeamInitStruct::set_sig_z_jitter,
          "bunch length RMS jitter"
      )
      .def_property(
          "sig_pz_jitter",
          &BeamInitStruct::sig_pz_jitter,
          &BeamInitStruct::set_sig_pz_jitter,
          "RMS pz spread jitter"
      )
      .def_property(
          "n_particle",
          &BeamInitStruct::n_particle,
          &BeamInitStruct::set_n_particle,
          "Number of particles per bunch."
      )
      .def_property(
          "renorm_center",
          &BeamInitStruct::renorm_center,
          &BeamInitStruct::set_renorm_center,
          "Renormalize centroid?"
      )
      .def_property(
          "renorm_sigma",
          &BeamInitStruct::renorm_sigma,
          &BeamInitStruct::set_renorm_sigma,
          "Renormalize sigma?"
      )
      .def_property(
          "random_engine",
          &BeamInitStruct::random_engine,
          &BeamInitStruct::set_random_engine,
          "Or 'quasi'. Random number engine to use."
      )
      .def_property(
          "random_gauss_converter",
          &BeamInitStruct::random_gauss_converter,
          &BeamInitStruct::set_random_gauss_converter,
          "Or 'quick'. Uniform to gauss conversion method."
      )
      .def_property(
          "random_sigma_cutoff",
          &BeamInitStruct::random_sigma_cutoff,
          &BeamInitStruct::set_random_sigma_cutoff,
          "Cut-off in sigmas."
      )
      .def_property(
          "a_norm_emit",
          &BeamInitStruct::a_norm_emit,
          &BeamInitStruct::set_a_norm_emit,
          "a-mode normalized emittance (emit * beta * gamma)"
      )
      .def_property(
          "b_norm_emit",
          &BeamInitStruct::b_norm_emit,
          &BeamInitStruct::set_b_norm_emit,
          "b-mode normalized emittance (emit * beta * gamma)"
      )
      .def_property(
          "a_emit",
          &BeamInitStruct::a_emit,
          &BeamInitStruct::set_a_emit,
          "a-mode emittance"
      )
      .def_property(
          "b_emit",
          &BeamInitStruct::b_emit,
          &BeamInitStruct::set_b_emit,
          "b-mode emittance"
      )
      .def_property(
          "dPz_dz",
          &BeamInitStruct::dPz_dz,
          &BeamInitStruct::set_dPz_dz,
          "Correlation of Pz with long position."
      )
      .def_property(
          "center",
          &BeamInitStruct::center,
          &BeamInitStruct::set_center,
          "Bench phase space center offset relative to reference."
      )
      .def_property(
          "t_offset",
          &BeamInitStruct::t_offset,
          &BeamInitStruct::set_t_offset,
          "Time center offset"
      )
      .def_property(
          "dt_bunch",
          &BeamInitStruct::dt_bunch,
          &BeamInitStruct::set_dt_bunch,
          "Time between bunches."
      )
      .def_property("sig_z", &BeamInitStruct::sig_z, &BeamInitStruct::set_sig_z, "Z sigma in m.")
      .def_property("sig_pz", &BeamInitStruct::sig_pz, &BeamInitStruct::set_sig_pz, "pz sigma")
      .def_property(
          "bunch_charge",
          &BeamInitStruct::bunch_charge,
          &BeamInitStruct::set_bunch_charge,
          "charge (Coul) in a bunch."
      )
      .def_property(
          "n_bunch",
          &BeamInitStruct::n_bunch,
          &BeamInitStruct::set_n_bunch,
          "Number of bunches."
      )
      .def_property(
          "ix_turn",
          &BeamInitStruct::ix_turn,
          &BeamInitStruct::set_ix_turn,
          "Turn index used to adjust particles time if needed."
      )
      .def_property(
          "species",
          &BeamInitStruct::species,
          &BeamInitStruct::set_species,
          "'positron', etc. '' => use referece particle."
      )
      .def_property(
          "full_6D_coupling_calc",
          &BeamInitStruct::full_6D_coupling_calc,
          &BeamInitStruct::set_full_6D_coupling_calc,
          "Use V from 6x6 1-turn mat to match distribution? Else use 4x4 1-turn mat used."
      )
      .def_property(
          "use_particle_start",
          &BeamInitStruct::use_particle_start,
          &BeamInitStruct::set_use_particle_start,
          "Use lat%particle_start instead of beam_init%center, %t_offset, and %spin?"
      )
      .def_property(
          "use_t_coords",
          &BeamInitStruct::use_t_coords,
          &BeamInitStruct::set_use_t_coords,
          "If true, the distributions will be taken as in t-coordinates"
      )
      .def_property(
          "use_z_as_t",
          &BeamInitStruct::use_z_as_t,
          &BeamInitStruct::set_use_z_as_t,
          "Only used if  use_t_coords = .true. If true,  z describes the t distribution If false, "
          "z describes the s distribution"
      )
      .def_property(
          "file_name",
          &BeamInitStruct::file_name,
          &BeamInitStruct::set_file_name,
          "OLD!! DO NOT USE!!"
      )

      .def("__repr__", [](const BeamInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BeamInitStruct &self) {
            return BeamInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BeamInitStruct &self, py::dict &memo) { return BeamInitStruct(self); }
      )

      ;

  // 1D BeamInitStruct arrays are not used in structs/routines
  // 2D BeamInitStruct arrays are not used in structs/routines
  // 3D BeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// beam_struct
void init_beam_struct(py::module &m, py::class_<BeamStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly("bunch", &BeamStruct::bunch)

      .def("__repr__", [](const BeamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BeamStruct &self) {
            return BeamStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const BeamStruct &self, py::dict &memo) { return BeamStruct(self); })

      ;

  // 1D BeamStruct arrays are not used in structs/routines
  // 2D BeamStruct arrays are not used in structs/routines
  // 3D BeamStruct arrays are not used in structs/routines
}

// =============================================================================
// bmad_common_struct
void init_bmad_common_struct(py::module &m, py::class_<BmadCommonStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("max_aperture_limit") = py::none(),
         py::arg("d_orb") = py::none(),
         py::arg("default_ds_step") = py::none(),
         py::arg("significant_length") = py::none(),
         py::arg("rel_tol_tracking") = py::none(),
         py::arg("abs_tol_tracking") = py::none(),
         py::arg("rel_tol_adaptive_tracking") = py::none(),
         py::arg("abs_tol_adaptive_tracking") = py::none(),
         py::arg("init_ds_adaptive_tracking") = py::none(),
         py::arg("min_ds_adaptive_tracking") = py::none(),
         py::arg("fatal_ds_adaptive_tracking") = py::none(),
         py::arg("autoscale_amp_abs_tol") = py::none(),
         py::arg("autoscale_amp_rel_tol") = py::none(),
         py::arg("autoscale_phase_tol") = py::none(),
         py::arg("electric_dipole_moment") = py::none(),
         py::arg("synch_rad_scale") = py::none(),
         py::arg("sad_eps_scale") = py::none(),
         py::arg("sad_amp_max") = py::none(),
         py::arg("sad_n_div_max") = py::none(),
         py::arg("taylor_order") = py::none(),
         py::arg("runge_kutta_order") = py::none(),
         py::arg("default_integ_order") = py::none(),
         py::arg("max_num_runge_kutta_step") = py::none(),
         py::arg("rf_phase_below_transition_ref") = py::none(),
         py::arg("sr_wakes_on") = py::none(),
         py::arg("lr_wakes_on") = py::none(),
         py::arg("auto_bookkeeper") = py::none(),
         py::arg("high_energy_space_charge_on") = py::none(),
         py::arg("csr_and_space_charge_on") = py::none(),
         py::arg("spin_tracking_on") = py::none(),
         py::arg("spin_sokolov_ternov_flipping_on") = py::none(),
         py::arg("radiation_damping_on") = py::none(),
         py::arg("radiation_zero_average") = py::none(),
         py::arg("radiation_fluctuations_on") = py::none(),
         py::arg("conserve_taylor_maps") = py::none(),
         py::arg("absolute_time_tracking") = py::none(),
         py::arg("absolute_time_ref_shift") = py::none(),
         py::arg("convert_to_kinetic_momentum") = py::none(),
         py::arg("normalize_twiss") = py::none(),
         py::arg("aperture_limit_on") = py::none(),
         py::arg("spin_n0_direction_user_set") = py::none(),
         py::arg("debug") = py::none()
  )
      .def_property(
          "max_aperture_limit",
          &BmadCommonStruct::max_aperture_limit,
          &BmadCommonStruct::set_max_aperture_limit,
          "Max Aperture."
      )
      .def_property(
          "d_orb",
          &BmadCommonStruct::d_orb,
          &BmadCommonStruct::set_d_orb,
          "Orbit deltas for the mat6 via tracking calc."
      )
      .def_property(
          "default_ds_step",
          &BmadCommonStruct::default_ds_step,
          &BmadCommonStruct::set_default_ds_step,
          "Default integration step for eles without an explicit step calc."
      )
      .def_property(
          "significant_length",
          &BmadCommonStruct::significant_length,
          &BmadCommonStruct::set_significant_length,
          "meter"
      )
      .def_property(
          "rel_tol_tracking",
          &BmadCommonStruct::rel_tol_tracking,
          &BmadCommonStruct::set_rel_tol_tracking,
          "Closed orbit relative tolerance."
      )
      .def_property(
          "abs_tol_tracking",
          &BmadCommonStruct::abs_tol_tracking,
          &BmadCommonStruct::set_abs_tol_tracking,
          "Closed orbit absolute tolerance."
      )
      .def_property(
          "rel_tol_adaptive_tracking",
          &BmadCommonStruct::rel_tol_adaptive_tracking,
          &BmadCommonStruct::set_rel_tol_adaptive_tracking,
          "Runge-Kutta tracking relative tolerance."
      )
      .def_property(
          "abs_tol_adaptive_tracking",
          &BmadCommonStruct::abs_tol_adaptive_tracking,
          &BmadCommonStruct::set_abs_tol_adaptive_tracking,
          "Runge-Kutta tracking absolute tolerance."
      )
      .def_property(
          "init_ds_adaptive_tracking",
          &BmadCommonStruct::init_ds_adaptive_tracking,
          &BmadCommonStruct::set_init_ds_adaptive_tracking,
          "Initial step size"
      )
      .def_property(
          "min_ds_adaptive_tracking",
          &BmadCommonStruct::min_ds_adaptive_tracking,
          &BmadCommonStruct::set_min_ds_adaptive_tracking,
          "Min step size to take."
      )
      .def_property(
          "fatal_ds_adaptive_tracking",
          &BmadCommonStruct::fatal_ds_adaptive_tracking,
          &BmadCommonStruct::set_fatal_ds_adaptive_tracking,
          "If actual step size is below this particle is lost."
      )
      .def_property(
          "autoscale_amp_abs_tol",
          &BmadCommonStruct::autoscale_amp_abs_tol,
          &BmadCommonStruct::set_autoscale_amp_abs_tol,
          "Autoscale absolute amplitude tolerance (eV)."
      )
      .def_property(
          "autoscale_amp_rel_tol",
          &BmadCommonStruct::autoscale_amp_rel_tol,
          &BmadCommonStruct::set_autoscale_amp_rel_tol,
          "Autoscale relative amplitude tolerance"
      )
      .def_property(
          "autoscale_phase_tol",
          &BmadCommonStruct::autoscale_phase_tol,
          &BmadCommonStruct::set_autoscale_phase_tol,
          "Autoscale phase tolerance."
      )
      .def_property(
          "electric_dipole_moment",
          &BmadCommonStruct::electric_dipole_moment,
          &BmadCommonStruct::set_electric_dipole_moment,
          "Particle's EDM. Call set_ptc to transfer value to PTC."
      )
      .def_property(
          "synch_rad_scale",
          &BmadCommonStruct::synch_rad_scale,
          &BmadCommonStruct::set_synch_rad_scale,
          "Synch radiation kick scale. 1 => normal, 0 => no kicks."
      )
      .def_property(
          "sad_eps_scale",
          &BmadCommonStruct::sad_eps_scale,
          &BmadCommonStruct::set_sad_eps_scale,
          "Used in sad_mult step length calc."
      )
      .def_property(
          "sad_amp_max",
          &BmadCommonStruct::sad_amp_max,
          &BmadCommonStruct::set_sad_amp_max,
          "Used in sad_mult step length calc."
      )
      .def_property(
          "sad_n_div_max",
          &BmadCommonStruct::sad_n_div_max,
          &BmadCommonStruct::set_sad_n_div_max,
          "Used in sad_mult step length calc."
      )
      .def_property(
          "taylor_order",
          &BmadCommonStruct::taylor_order,
          &BmadCommonStruct::set_taylor_order,
          "Taylor order to use. 0 -> default = ptc_private%taylor_order_saved."
      )
      .def_property(
          "runge_kutta_order",
          &BmadCommonStruct::runge_kutta_order,
          &BmadCommonStruct::set_runge_kutta_order,
          "Runge Kutta order."
      )
      .def_property(
          "default_integ_order",
          &BmadCommonStruct::default_integ_order,
          &BmadCommonStruct::set_default_integ_order,
          "PTC integration order."
      )
      .def_property(
          "max_num_runge_kutta_step",
          &BmadCommonStruct::max_num_runge_kutta_step,
          &BmadCommonStruct::set_max_num_runge_kutta_step,
          "Maximum number of RK steps before particle is considered lost."
      )
      .def_property(
          "rf_phase_below_transition_ref",
          &BmadCommonStruct::rf_phase_below_transition_ref,
          &BmadCommonStruct::set_rf_phase_below_transition_ref,
          "Autoscale uses below transition stable point for RFCavities?"
      )
      .def_property(
          "sr_wakes_on",
          &BmadCommonStruct::sr_wakes_on,
          &BmadCommonStruct::set_sr_wakes_on,
          "Short range wakefields?"
      )
      .def_property(
          "lr_wakes_on",
          &BmadCommonStruct::lr_wakes_on,
          &BmadCommonStruct::set_lr_wakes_on,
          "Long range wakefields"
      )
      .def_property(
          "auto_bookkeeper",
          &BmadCommonStruct::auto_bookkeeper,
          &BmadCommonStruct::set_auto_bookkeeper,
          "Deprecated and no longer used."
      )
      .def_property(
          "high_energy_space_charge_on",
          &BmadCommonStruct::high_energy_space_charge_on,
          &BmadCommonStruct::set_high_energy_space_charge_on,
          "High energy space charge effect switch."
      )
      .def_property(
          "csr_and_space_charge_on",
          &BmadCommonStruct::csr_and_space_charge_on,
          &BmadCommonStruct::set_csr_and_space_charge_on,
          "Space charge switch."
      )
      .def_property(
          "spin_tracking_on",
          &BmadCommonStruct::spin_tracking_on,
          &BmadCommonStruct::set_spin_tracking_on,
          "spin tracking?"
      )
      .def_property(
          "spin_sokolov_ternov_flipping_on",
          &BmadCommonStruct::spin_sokolov_ternov_flipping_on,
          &BmadCommonStruct::set_spin_sokolov_ternov_flipping_on,
          "Spin flipping during synchrotron radiation emission?"
      )
      .def_property(
          "radiation_damping_on",
          &BmadCommonStruct::radiation_damping_on,
          &BmadCommonStruct::set_radiation_damping_on,
          "Radiation damping toggle."
      )
      .def_property(
          "radiation_zero_average",
          &BmadCommonStruct::radiation_zero_average,
          &BmadCommonStruct::set_radiation_zero_average,
          "Shift damping to be zero on the zero orbit to get rid of sawtooth?"
      )
      .def_property(
          "radiation_fluctuations_on",
          &BmadCommonStruct::radiation_fluctuations_on,
          &BmadCommonStruct::set_radiation_fluctuations_on,
          "Radiation fluctuations toggle."
      )
      .def_property(
          "conserve_taylor_maps",
          &BmadCommonStruct::conserve_taylor_maps,
          &BmadCommonStruct::set_conserve_taylor_maps,
          "Enable bookkeeper to set ele%taylor_map_includes_offsets = F?"
      )
      .def_property(
          "absolute_time_tracking",
          &BmadCommonStruct::absolute_time_tracking,
          &BmadCommonStruct::set_absolute_time_tracking,
          "Absolute or relative time tracking?"
      )
      .def_property(
          "absolute_time_ref_shift",
          &BmadCommonStruct::absolute_time_ref_shift,
          &BmadCommonStruct::set_absolute_time_ref_shift,
          "Apply reference time shift when using absolute time tracking?"
      )
      .def_property(
          "convert_to_kinetic_momentum",
          &BmadCommonStruct::convert_to_kinetic_momentum,
          &BmadCommonStruct::set_convert_to_kinetic_momentum,
          "Cancel kicks due to finite vector potential when doing symplectic tracking? Set to True "
          "to test symp_lie_bmad against runge_kutta."
      )
      .def_property(
          "normalize_twiss",
          &BmadCommonStruct::normalize_twiss,
          &BmadCommonStruct::set_normalize_twiss,
          "Normalize matrix when computing Twiss for off-energy ref?"
      )
      .def_property(
          "aperture_limit_on",
          &BmadCommonStruct::aperture_limit_on,
          &BmadCommonStruct::set_aperture_limit_on,
          "Use apertures in tracking?"
      )
      .def_property(
          "spin_n0_direction_user_set",
          &BmadCommonStruct::spin_n0_direction_user_set,
          &BmadCommonStruct::set_spin_n0_direction_user_set,
          "User sets direction of n0 for closed geometry branches?"
      )
      .def_property(
          "debug",
          &BmadCommonStruct::debug,
          &BmadCommonStruct::set_debug,
          "Used for code debugging."
      )

      .def("__repr__", [](const BmadCommonStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BmadCommonStruct &self) {
            return BmadCommonStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BmadCommonStruct &self, py::dict &memo) { return BmadCommonStruct(self); }
      )

      ;

  // 1D BmadCommonStruct arrays are not used in structs/routines
  // 2D BmadCommonStruct arrays are not used in structs/routines
  // 3D BmadCommonStruct arrays are not used in structs/routines
}

// =============================================================================
// bmad_normal_form_struct
void init_bmad_normal_form_struct(py::module &m, py::class_<BmadNormalFormStruct> &cls) {
  cls.def(py::init<optional_ref<const EleStruct>>(), py::arg("ele_origin") = py::none())
      .def_property(
          "ele_origin",
          &BmadNormalFormStruct::ele_origin,
          &BmadNormalFormStruct::set_ele_origin,
          "Element at which the on-turn map was created."
      )
      .def_property_readonly(
          "M",
          &BmadNormalFormStruct::M,
          "One-turn taylor map: M = A o N o A_inv, N = exp(:h:)"
      )
      .def_property_readonly("A", &BmadNormalFormStruct::A, "Map from Floquet -> Lab coordinates")
      .def_property_readonly(
          "A_inv",
          &BmadNormalFormStruct::A_inv,
          "Map from Lab -> Floquet coordinates"
      )
      .def_property_readonly(
          "dhdj",
          &BmadNormalFormStruct::dhdj,
          "Nonlinear tune function operating on Floquet coordinates"
      )
      .def_property_readonly(
          "F",
          &BmadNormalFormStruct::F,
          "Vector field factorization in phasor basis:"
      )
      .def_property_readonly("L", &BmadNormalFormStruct::L, "L component")
      .def_property_readonly("h", &BmadNormalFormStruct::h)

      .def("__repr__", [](const BmadNormalFormStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BmadNormalFormStruct &self) {
            return BmadNormalFormStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BmadNormalFormStruct &self, py::dict &memo) {
            return BmadNormalFormStruct(self);
          }
      )

      ;

  // 1D BmadNormalFormStruct arrays are not used in structs/routines
  // 2D BmadNormalFormStruct arrays are not used in structs/routines
  // 3D BmadNormalFormStruct arrays are not used in structs/routines
}

// =============================================================================
// bookkeeping_state_struct
void init_bookkeeping_state_struct(py::module &m, py::class_<BookkeepingStateStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("attributes") = py::none(),
         py::arg("control") = py::none(),
         py::arg("floor_position") = py::none(),
         py::arg("s_position") = py::none(),
         py::arg("ref_energy") = py::none(),
         py::arg("mat6") = py::none(),
         py::arg("rad_int") = py::none(),
         py::arg("ptc") = py::none(),
         py::arg("has_misalign") = py::none()
  )
      .def_property(
          "attributes",
          &BookkeepingStateStruct::attributes,
          &BookkeepingStateStruct::set_attributes,
          "Element dependent attributes: super_ok$, ok$ or stale$"
      )
      .def_property(
          "control",
          &BookkeepingStateStruct::control,
          &BookkeepingStateStruct::set_control,
          "Lord/slave bookkeeping status: super_ok$, ok$ or stale$"
      )
      .def_property(
          "floor_position",
          &BookkeepingStateStruct::floor_position,
          &BookkeepingStateStruct::set_floor_position,
          "Global (floor) geometry: super_ok$, ok$ or stale$"
      )
      .def_property(
          "s_position",
          &BookkeepingStateStruct::s_position,
          &BookkeepingStateStruct::set_s_position,
          "Longitudinal position & element length: super_ok$, ok$ or stale$"
      )
      .def_property(
          "ref_energy",
          &BookkeepingStateStruct::ref_energy,
          &BookkeepingStateStruct::set_ref_energy,
          "Reference energy and ref time: super_ok$, ok$ or stale$"
      )
      .def_property(
          "mat6",
          &BookkeepingStateStruct::mat6,
          &BookkeepingStateStruct::set_mat6,
          "Linear transfer map status: super_ok$, ok$ or stale$"
      )
      .def_property(
          "rad_int",
          &BookkeepingStateStruct::rad_int,
          &BookkeepingStateStruct::set_rad_int,
          "Radiation integrals cache status"
      )
      .def_property(
          "ptc",
          &BookkeepingStateStruct::ptc,
          &BookkeepingStateStruct::set_ptc,
          "Associated PTC fibre (or layout) status."
      )
      .def_property(
          "has_misalign",
          &BookkeepingStateStruct::has_misalign,
          &BookkeepingStateStruct::set_has_misalign,
          "Used to avoid unnecessary calls to offset_particle."
      )

      .def("__repr__", [](const BookkeepingStateStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BookkeepingStateStruct &self) {
            return BookkeepingStateStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BookkeepingStateStruct &self, py::dict &memo) {
            return BookkeepingStateStruct(self);
          }
      )

      ;

  // 1D BookkeepingStateStruct arrays are not used in structs/routines
  // 2D BookkeepingStateStruct arrays are not used in structs/routines
  // 3D BookkeepingStateStruct arrays are not used in structs/routines
}

// =============================================================================
// bpm_phase_coupling_struct
void init_bpm_phase_coupling_struct(py::module &m, py::class_<BpmPhaseCouplingStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("K_22a") = py::none(),
         py::arg("K_12a") = py::none(),
         py::arg("K_11b") = py::none(),
         py::arg("K_12b") = py::none(),
         py::arg("Cbar22_a") = py::none(),
         py::arg("Cbar12_a") = py::none(),
         py::arg("Cbar11_b") = py::none(),
         py::arg("Cbar12_b") = py::none(),
         py::arg("phi_a") = py::none(),
         py::arg("phi_b") = py::none()
  )
      .def_property(
          "K_22a",
          &BpmPhaseCouplingStruct::K_22a,
          &BpmPhaseCouplingStruct::set_K_22a,
          "In-phase y/x for a-mode oscillations."
      )
      .def_property(
          "K_12a",
          &BpmPhaseCouplingStruct::K_12a,
          &BpmPhaseCouplingStruct::set_K_12a,
          "Out-of-phase y/x for a-mode oscillations."
      )
      .def_property(
          "K_11b",
          &BpmPhaseCouplingStruct::K_11b,
          &BpmPhaseCouplingStruct::set_K_11b,
          "In-phase x/y for b-mode oscillations."
      )
      .def_property(
          "K_12b",
          &BpmPhaseCouplingStruct::K_12b,
          &BpmPhaseCouplingStruct::set_K_12b,
          "Out-of-phase x/y for b-mode oscillations."
      )
      .def_property(
          "Cbar22_a",
          &BpmPhaseCouplingStruct::Cbar22_a,
          &BpmPhaseCouplingStruct::set_Cbar22_a,
          "Cbar22 as calculated from K_22a."
      )
      .def_property(
          "Cbar12_a",
          &BpmPhaseCouplingStruct::Cbar12_a,
          &BpmPhaseCouplingStruct::set_Cbar12_a,
          "Cbar12 as calculated from K_12a."
      )
      .def_property(
          "Cbar11_b",
          &BpmPhaseCouplingStruct::Cbar11_b,
          &BpmPhaseCouplingStruct::set_Cbar11_b,
          "Cbar11 as calculated from K_11b."
      )
      .def_property(
          "Cbar12_b",
          &BpmPhaseCouplingStruct::Cbar12_b,
          &BpmPhaseCouplingStruct::set_Cbar12_b,
          "Cbar12 as calculated from K_12b."
      )
      .def_property(
          "phi_a",
          &BpmPhaseCouplingStruct::phi_a,
          &BpmPhaseCouplingStruct::set_phi_a,
          "a-mode betatron phase."
      )
      .def_property(
          "phi_b",
          &BpmPhaseCouplingStruct::phi_b,
          &BpmPhaseCouplingStruct::set_phi_b,
          "b-mode betatron phase."
      )

      .def("__repr__", [](const BpmPhaseCouplingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BpmPhaseCouplingStruct &self) {
            return BpmPhaseCouplingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BpmPhaseCouplingStruct &self, py::dict &memo) {
            return BpmPhaseCouplingStruct(self);
          }
      )

      ;

  // 1D BpmPhaseCouplingStruct arrays are not used in structs/routines
  // 2D BpmPhaseCouplingStruct arrays are not used in structs/routines
  // 3D BpmPhaseCouplingStruct arrays are not used in structs/routines
}

// =============================================================================
// branch_struct
void init_branch_struct(py::module &m, py::class_<BranchStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             optional_ref<const LatStruct>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const ModeInfoStruct>,
             optional_ref<const LatParamStruct>,
             optional_ref<const CoordStruct>>(),
         py::arg("name") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_from_branch") = py::none(),
         py::arg("ix_from_ele") = py::none(),
         py::arg("ix_to_ele") = py::none(),
         py::arg("ix_fixer") = py::none(),
         py::arg("n_ele_track") = py::none(),
         py::arg("n_ele_max") = py::none(),
         py::arg("lat") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("z") = py::none(),
         py::arg("param") = py::none(),
         py::arg("particle_start") = py::none()
  )
      .def_property(
          "name",
          &BranchStruct::name,
          &BranchStruct::set_name,
          "Name of line that defines the branch."
      )
      .def_property(
          "ix_branch",
          &BranchStruct::ix_branch,
          &BranchStruct::set_ix_branch,
          "Index of this branch. 0 => Main branch"
      )
      .def_property(
          "ix_from_branch",
          &BranchStruct::ix_from_branch,
          &BranchStruct::set_ix_from_branch,
          "-1 => No creating fork element to this branch."
      )
      .def_property(
          "ix_from_ele",
          &BranchStruct::ix_from_ele,
          &BranchStruct::set_ix_from_ele,
          "Index of creating fork element which forks to this branch."
      )
      .def_property(
          "ix_to_ele",
          &BranchStruct::ix_to_ele,
          &BranchStruct::set_ix_to_ele,
          "Index of element in this branch that creating fork element forks to."
      )
      .def_property(
          "ix_fixer",
          &BranchStruct::ix_fixer,
          &BranchStruct::set_ix_fixer,
          "Index of active fixer or beginning_ele element."
      )
      .def_property("n_ele_track", &BranchStruct::n_ele_track, &BranchStruct::set_n_ele_track)
      .def_property("n_ele_max", &BranchStruct::n_ele_max, &BranchStruct::set_n_ele_max)
      .def_property("lat", &BranchStruct::lat, &BranchStruct::set_lat)
      .def_property(
          "a",
          &BranchStruct::a,
          &BranchStruct::set_a,
          "Note: Tunes are the fractional part."
      )
      .def_property(
          "b",
          &BranchStruct::b,
          &BranchStruct::set_b,
          "Note: Tunes are the fractional part."
      )
      .def_property(
          "z",
          &BranchStruct::z,
          &BranchStruct::set_z,
          "Note: Tunes are the fractional part."
      )
      .def_property_readonly("ele", &BranchStruct::ele)
      .def_property("param", &BranchStruct::param, &BranchStruct::set_param)
      .def_property(
          "particle_start",
          &BranchStruct::particle_start,
          &BranchStruct::set_particle_start
      )
      .def_property_readonly("wall3d", &BranchStruct::wall3d)
      .def_static(
          "new_array1d",
          [](int sz) { return BranchStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const BranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BranchStruct &self) {
            return BranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BranchStruct &self, py::dict &memo) { return BranchStruct(self); }
      )

      ;

  bind_1d_type_array_pair<BranchStructArray1D, BranchStructAlloc1D>(
      m,
      "BranchStructArray1D",
      "BranchStructAlloc1D"
  );
  // 2D BranchStruct arrays are not used in structs/routines
  // 3D BranchStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_params_struct
void init_bunch_params_struct(py::module &m, py::class_<BunchParamsStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const CoordStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("centroid") = py::none(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("z") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("c") = py::none(),
         py::arg("sigma") = py::none(),
         py::arg("rel_max") = py::none(),
         py::arg("rel_min") = py::none(),
         py::arg("s") = py::none(),
         py::arg("t") = py::none(),
         py::arg("sigma_t") = py::none(),
         py::arg("charge_live") = py::none(),
         py::arg("charge_tot") = py::none(),
         py::arg("n_particle_tot") = py::none(),
         py::arg("n_particle_live") = py::none(),
         py::arg("n_particle_lost_in_ele") = py::none(),
         py::arg("n_good_steps") = py::none(),
         py::arg("n_bad_steps") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("location") = py::none(),
         py::arg("twiss_valid") = py::none()
  )
      .def_property(
          "centroid",
          &BunchParamsStruct::centroid,
          &BunchParamsStruct::set_centroid,
          "Lab frame"
      )
      .def_property(
          "x",
          &BunchParamsStruct::x,
          &BunchParamsStruct::set_x,
          "Projected Twiss parameters"
      )
      .def_property(
          "y",
          &BunchParamsStruct::y,
          &BunchParamsStruct::set_y,
          "Projected Twiss parameters"
      )
      .def_property(
          "z",
          &BunchParamsStruct::z,
          &BunchParamsStruct::set_z,
          "Projected Twiss parameters"
      )
      .def_property(
          "a",
          &BunchParamsStruct::a,
          &BunchParamsStruct::set_a,
          "Normal mode twiss parameters"
      )
      .def_property(
          "b",
          &BunchParamsStruct::b,
          &BunchParamsStruct::set_b,
          "Normal mode twiss parameters"
      )
      .def_property(
          "c",
          &BunchParamsStruct::c,
          &BunchParamsStruct::set_c,
          "Normal mode twiss parameters"
      )
      .def_property(
          "sigma",
          &BunchParamsStruct::sigma,
          &BunchParamsStruct::set_sigma,
          "beam size matrix"
      )
      .def_property(
          "rel_max",
          &BunchParamsStruct::rel_max,
          &BunchParamsStruct::set_rel_max,
          "Max orbit relative to centroid. 7 -> time."
      )
      .def_property(
          "rel_min",
          &BunchParamsStruct::rel_min,
          &BunchParamsStruct::set_rel_min,
          "Min orbit relative to_centroid. 7 -> time."
      )
      .def_property("s", &BunchParamsStruct::s, &BunchParamsStruct::set_s, "Longitudinal position.")
      .def_property("t", &BunchParamsStruct::t, &BunchParamsStruct::set_t, "Time.")
      .def_property(
          "sigma_t",
          &BunchParamsStruct::sigma_t,
          &BunchParamsStruct::set_sigma_t,
          "RMS of time spread."
      )
      .def_property(
          "charge_live",
          &BunchParamsStruct::charge_live,
          &BunchParamsStruct::set_charge_live,
          "Charge of all non-lost particle"
      )
      .def_property(
          "charge_tot",
          &BunchParamsStruct::charge_tot,
          &BunchParamsStruct::set_charge_tot,
          "Charge of all particles."
      )
      .def_property(
          "n_particle_tot",
          &BunchParamsStruct::n_particle_tot,
          &BunchParamsStruct::set_n_particle_tot,
          "Total number of particles"
      )
      .def_property(
          "n_particle_live",
          &BunchParamsStruct::n_particle_live,
          &BunchParamsStruct::set_n_particle_live,
          "Number of non-lost particles"
      )
      .def_property(
          "n_particle_lost_in_ele",
          &BunchParamsStruct::n_particle_lost_in_ele,
          &BunchParamsStruct::set_n_particle_lost_in_ele,
          "Number lost in element (not calculated by Bmad)"
      )
      .def_property(
          "n_good_steps",
          &BunchParamsStruct::n_good_steps,
          &BunchParamsStruct::set_n_good_steps,
          "Number of good steps (set when tracking with space charge)"
      )
      .def_property(
          "n_bad_steps",
          &BunchParamsStruct::n_bad_steps,
          &BunchParamsStruct::set_n_bad_steps,
          "Number of bad steps (set when tracking with space charge)"
      )
      .def_property(
          "ix_ele",
          &BunchParamsStruct::ix_ele,
          &BunchParamsStruct::set_ix_ele,
          "Lattice element where params evaluated at."
      )
      .def_property(
          "location",
          &BunchParamsStruct::location,
          &BunchParamsStruct::set_location,
          "Location in element: upstream_end$, inside$, or downstream_end$"
      )
      .def_property(
          "twiss_valid",
          &BunchParamsStruct::twiss_valid,
          &BunchParamsStruct::set_twiss_valid,
          "Is the data here valid? Note: IF there is no energy variation (RF off) twiss_valid may "
          "be true but in this case the z-twiss will not be valid."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchParamsStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchParamsStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const BunchParamsStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchParamsStruct &self) {
            return BunchParamsStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BunchParamsStruct &self, py::dict &memo) { return BunchParamsStruct(self); }
      )

      ;

  bind_1d_type_array_pair<BunchParamsStructArray1D, BunchParamsStructAlloc1D>(
      m,
      "BunchParamsStructArray1D",
      "BunchParamsStructAlloc1D"
  );
  // 2D BunchParamsStruct arrays are not used in structs/routines
  // 3D BunchParamsStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_struct
void init_bunch_struct(py::module &m, py::class_<BunchStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<int>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("ix_z") = py::none(),
         py::arg("charge_tot") = py::none(),
         py::arg("charge_live") = py::none(),
         py::arg("z_center") = py::none(),
         py::arg("t_center") = py::none(),
         py::arg("t0") = py::none(),
         py::arg("drift_between_t_and_s") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_bunch") = py::none(),
         py::arg("ix_turn") = py::none(),
         py::arg("n_live") = py::none(),
         py::arg("n_good") = py::none(),
         py::arg("n_bad") = py::none()
  )
      .def_property_readonly("particle", &BunchStruct::particle)
      .def_property(
          "ix_z",
          &BunchStruct::ix_z,
          &BunchStruct::set_ix_z,
          "bunch%ix_z(1) is index of head particle, etc."
      )
      .def_property(
          "charge_tot",
          &BunchStruct::charge_tot,
          &BunchStruct::set_charge_tot,
          "Total charge in a bunch (Coul)."
      )
      .def_property(
          "charge_live",
          &BunchStruct::charge_live,
          &BunchStruct::set_charge_live,
          "Charge of live particles (Coul)."
      )
      .def_property(
          "z_center",
          &BunchStruct::z_center,
          &BunchStruct::set_z_center,
          "Longitudinal center of bunch at creation time. Note: Generally, z_center of bunch #1 is "
          "0 and z_center of the other bunches is negative."
      )
      .def_property(
          "t_center",
          &BunchStruct::t_center,
          &BunchStruct::set_t_center,
          "Center of bunch at creation time relative to head bunch."
      )
      .def_property(
          "t0",
          &BunchStruct::t0,
          &BunchStruct::set_t0,
          "Used by track1_bunch_space_charge for tracking so particles have constant t."
      )
      .def_property(
          "drift_between_t_and_s",
          &BunchStruct::drift_between_t_and_s,
          &BunchStruct::set_drift_between_t_and_s,
          "Drift (ignore any fields) instead of tracking to speed up the calculation? This can "
          "only be done under certain circumstances."
      )
      .def_property(
          "ix_ele",
          &BunchStruct::ix_ele,
          &BunchStruct::set_ix_ele,
          "Nominal element bunch is at. But, EG, dead particles can be someplace else."
      )
      .def_property(
          "ix_bunch",
          &BunchStruct::ix_bunch,
          &BunchStruct::set_ix_bunch,
          "Bunch index. Head bunch = 1, etc."
      )
      .def_property(
          "ix_turn",
          &BunchStruct::ix_turn,
          &BunchStruct::set_ix_turn,
          "Turn index for long term tracking. ix_turn = 0 before end of first turn, etc."
      )
      .def_property("n_live", &BunchStruct::n_live, &BunchStruct::set_n_live)
      .def_property(
          "n_good",
          &BunchStruct::n_good,
          &BunchStruct::set_n_good,
          "Number of accepted steps when using adaptive step size control."
      )
      .def_property(
          "n_bad",
          &BunchStruct::n_bad,
          &BunchStruct::set_n_bad,
          "Number of rejected steps when using adaptive step size control."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const BunchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchStruct &self) {
            return BunchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BunchStruct &self, py::dict &memo) { return BunchStruct(self); }
      )

      ;

  bind_1d_type_array_pair<BunchStructArray1D, BunchStructAlloc1D>(
      m,
      "BunchStructArray1D",
      "BunchStructAlloc1D"
  );
  // 2D BunchStruct arrays are not used in structs/routines
  // 3D BunchStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_track_struct
void init_bunch_track_struct(py::module &m, py::class_<BunchTrackStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<int>>(),
         py::arg("ds_save") = py::none(),
         py::arg("n_pt") = py::none()
  )
      .def_property_readonly("pt", &BunchTrackStruct::pt, "Array indexed from 0")
      .def_property(
          "ds_save",
          &BunchTrackStruct::ds_save,
          &BunchTrackStruct::set_ds_save,
          "Min distance between points."
      )
      .def_property(
          "n_pt",
          &BunchTrackStruct::n_pt,
          &BunchTrackStruct::set_n_pt,
          "Track upper bound"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchTrackStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchTrackStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const BunchTrackStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchTrackStruct &self) {
            return BunchTrackStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BunchTrackStruct &self, py::dict &memo) { return BunchTrackStruct(self); }
      )

      ;

  bind_1d_type_array_pair<BunchTrackStructArray1D, BunchTrackStructAlloc1D>(
      m,
      "BunchTrackStructArray1D",
      "BunchTrackStructAlloc1D"
  );
  // 2D BunchTrackStruct arrays are not used in structs/routines
  // 3D BunchTrackStruct arrays are not used in structs/routines
}

// =============================================================================
// bicubic_cmplx_coef_struct
void init_bicubic_cmplx_coef_struct(py::module &m, py::class_<BicubicCmplxCoefStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<std::complex<double>>>>,
             std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("i_box") = py::none()
  )
      .def_property(
          "coef",
          &BicubicCmplxCoefStruct::coef,
          &BicubicCmplxCoefStruct::set_coef,
          "Coefs"
      )
      .def_property(
          "i_box",
          &BicubicCmplxCoefStruct::i_box,
          &BicubicCmplxCoefStruct::set_i_box,
          "index at lower box corner."
      )

      .def("__repr__", [](const BicubicCmplxCoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BicubicCmplxCoefStruct &self) {
            return BicubicCmplxCoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BicubicCmplxCoefStruct &self, py::dict &memo) {
            return BicubicCmplxCoefStruct(self);
          }
      )

      ;

  // 1D BicubicCmplxCoefStruct arrays are not used in structs/routines
  // 2D BicubicCmplxCoefStruct arrays are not used in structs/routines
  bind_FTypeArrayND<BicubicCmplxCoefStructArray3D>(m, "BicubicCmplxCoefStructArray3D");
}