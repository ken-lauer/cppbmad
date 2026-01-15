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
  cls.def(py::init<>())
      // BbuBeamStruct.bunch (1D_ALLOC_type - Bunches in the lattice
      .def_property_readonly("bunch", &BbuBeamStruct::bunch)
      // BbuBeamStruct.stage (1D_ALLOC_type -
      .def_property_readonly("stage", &BbuBeamStruct::stage)
      // BbuBeamStruct.ix_ele_bunch (1D_ALLOC_integer - element where bunch is
      .def_property_readonly("ix_ele_bunch", &BbuBeamStruct::ix_ele_bunch)
      // BbuBeamStruct.ix_bunch_head (0D_NOT_integer - Index to head bunch(:)
      .def_property(
          "ix_bunch_head",
          &BbuBeamStruct::ix_bunch_head,
          &BbuBeamStruct::set_ix_bunch_head
      )
      // BbuBeamStruct.ix_bunch_end (0D_NOT_integer - Index of the end bunch(:). -1 -> no bunches.
      .def_property("ix_bunch_end", &BbuBeamStruct::ix_bunch_end, &BbuBeamStruct::set_ix_bunch_end)
      // BbuBeamStruct.n_bunch_in_lat (0D_NOT_integer - Number of bunches transversing the lattice.
      .def_property(
          "n_bunch_in_lat",
          &BbuBeamStruct::n_bunch_in_lat,
          &BbuBeamStruct::set_n_bunch_in_lat
      )
      // BbuBeamStruct.ix_stage_voltage_max (0D_NOT_integer -
      .def_property(
          "ix_stage_voltage_max",
          &BbuBeamStruct::ix_stage_voltage_max,
          &BbuBeamStruct::set_ix_stage_voltage_max
      )
      // BbuBeamStruct.hom_voltage_max (0D_NOT_real -
      .def_property(
          "hom_voltage_max",
          &BbuBeamStruct::hom_voltage_max,
          &BbuBeamStruct::set_hom_voltage_max
      )
      // BbuBeamStruct.time_now (0D_NOT_real -
      .def_property("time_now", &BbuBeamStruct::time_now, &BbuBeamStruct::set_time_now)
      // BbuBeamStruct.one_turn_time (0D_NOT_real -
      .def_property(
          "one_turn_time",
          &BbuBeamStruct::one_turn_time,
          &BbuBeamStruct::set_one_turn_time
      )
      // BbuBeamStruct.rf_wavelength_max (0D_NOT_real -
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
  cls.def(py::init<>())
      // BbuParamStruct.lat_filename (0D_NOT_character - Bmad lattice file name
      .def_property(
          "lat_filename",
          &BbuParamStruct::lat_filename,
          &BbuParamStruct::set_lat_filename
      )
      // BbuParamStruct.lat2_filename (0D_NOT_character - Bmad lattice2 file name for secondary
      // parser
      .def_property(
          "lat2_filename",
          &BbuParamStruct::lat2_filename,
          &BbuParamStruct::set_lat2_filename
      )
      // BbuParamStruct.bunch_by_bunch_info_file (0D_NOT_character - For outputting bunch-by-bunch
      // info.
      .def_property(
          "bunch_by_bunch_info_file",
          &BbuParamStruct::bunch_by_bunch_info_file,
          &BbuParamStruct::set_bunch_by_bunch_info_file
      )
      // BbuParamStruct.hybridize (0D_NOT_logical - Combine non-hom elements to speed up simulation?
      .def_property("hybridize", &BbuParamStruct::hybridize, &BbuParamStruct::set_hybridize)
      // BbuParamStruct.write_digested_hybrid_lat (0D_NOT_logical - For debugging purposes.
      .def_property(
          "write_digested_hybrid_lat",
          &BbuParamStruct::write_digested_hybrid_lat,
          &BbuParamStruct::set_write_digested_hybrid_lat
      )
      // BbuParamStruct.write_voltage_vs_time_dat (0D_NOT_logical - For debugging purposes.
      .def_property(
          "write_voltage_vs_time_dat",
          &BbuParamStruct::write_voltage_vs_time_dat,
          &BbuParamStruct::set_write_voltage_vs_time_dat
      )
      // BbuParamStruct.keep_overlays_and_groups (0D_NOT_logical - Keep when hybridizing?
      .def_property(
          "keep_overlays_and_groups",
          &BbuParamStruct::keep_overlays_and_groups,
          &BbuParamStruct::set_keep_overlays_and_groups
      )
      // BbuParamStruct.keep_all_lcavities (0D_NOT_logical - Keep when hybridizing?
      .def_property(
          "keep_all_lcavities",
          &BbuParamStruct::keep_all_lcavities,
          &BbuParamStruct::set_keep_all_lcavities
      )
      // BbuParamStruct.use_taylor_for_hybrids (0D_NOT_logical - Use taylor map for hybrids when
      // true. Otherwise tracking method is linear.
      .def_property(
          "use_taylor_for_hybrids",
          &BbuParamStruct::use_taylor_for_hybrids,
          &BbuParamStruct::set_use_taylor_for_hybrids
      )
      // BbuParamStruct.stable_orbit_anal (0D_NOT_logical - Write stable_orbit.out and
      // hom_voltage.out?
      .def_property(
          "stable_orbit_anal",
          &BbuParamStruct::stable_orbit_anal,
          &BbuParamStruct::set_stable_orbit_anal
      )
      // BbuParamStruct.limit_factor (0D_NOT_real - Init_hom_amp * limit_factor = simulation
      // unstable limit
      .def_property(
          "limit_factor",
          &BbuParamStruct::limit_factor,
          &BbuParamStruct::set_limit_factor
      )
      // BbuParamStruct.simulation_turns_max (0D_NOT_real - Sets the duration of the simulation.
      .def_property(
          "simulation_turns_max",
          &BbuParamStruct::simulation_turns_max,
          &BbuParamStruct::set_simulation_turns_max
      )
      // BbuParamStruct.bunch_freq (0D_NOT_real - Freq in Hz.
      .def_property("bunch_freq", &BbuParamStruct::bunch_freq, &BbuParamStruct::set_bunch_freq)
      // BbuParamStruct.init_particle_offset (0D_NOT_real - Initial particle offset for particles
      // born in the first turn period.
      .def_property(
          "init_particle_offset",
          &BbuParamStruct::init_particle_offset,
          &BbuParamStruct::set_init_particle_offset
      )
      // BbuParamStruct.current (0D_NOT_real - Starting current (amps)
      .def_property("current", &BbuParamStruct::current, &BbuParamStruct::set_current)
      // BbuParamStruct.rel_tol (0D_NOT_real - Final threshold current accuracy.
      .def_property("rel_tol", &BbuParamStruct::rel_tol, &BbuParamStruct::set_rel_tol)
      // BbuParamStruct.drscan (0D_NOT_logical - If true, scan DR variable as in PRSTAB 7 (2004)
      // Fig. 3.
      .def_property("drscan", &BbuParamStruct::drscan, &BbuParamStruct::set_drscan)
      // BbuParamStruct.use_interpolated_threshold (0D_NOT_logical -
      .def_property(
          "use_interpolated_threshold",
          &BbuParamStruct::use_interpolated_threshold,
          &BbuParamStruct::set_use_interpolated_threshold
      )
      // BbuParamStruct.write_hom_info (0D_NOT_logical - Write HOM parameters to main output file?
      .def_property(
          "write_hom_info",
          &BbuParamStruct::write_hom_info,
          &BbuParamStruct::set_write_hom_info
      )
      // BbuParamStruct.elindex (0D_NOT_integer -
      .def_property("elindex", &BbuParamStruct::elindex, &BbuParamStruct::set_elindex)
      // BbuParamStruct.elname (0D_NOT_character - Element to step length for DRSCAN
      .def_property("elname", &BbuParamStruct::elname, &BbuParamStruct::set_elname)
      // BbuParamStruct.nstep (0D_NOT_integer - Number of steps for DRSCAN.
      .def_property("nstep", &BbuParamStruct::nstep, &BbuParamStruct::set_nstep)
      // BbuParamStruct.begdr (0D_NOT_real - Beginning DR value for DRSCAN.
      .def_property("begdr", &BbuParamStruct::begdr, &BbuParamStruct::set_begdr)
      // BbuParamStruct.enddr (0D_NOT_real - End DR value for DRSCAN.
      .def_property("enddr", &BbuParamStruct::enddr, &BbuParamStruct::set_enddr)
      // BbuParamStruct.nrep (0D_NOT_integer - Number of times to repeat threshold calculation
      .def_property("nrep", &BbuParamStruct::nrep, &BbuParamStruct::set_nrep)
      // BbuParamStruct.ran_seed (0D_NOT_integer - If set to 0, the output results will vary from
      // run to run.
      .def_property("ran_seed", &BbuParamStruct::ran_seed, &BbuParamStruct::set_ran_seed)
      // BbuParamStruct.hom_order_cutoff (0D_NOT_integer - If positive -> ignore HOM's with order
      // greater than this.
      .def_property(
          "hom_order_cutoff",
          &BbuParamStruct::hom_order_cutoff,
          &BbuParamStruct::set_hom_order_cutoff
      )
      // BbuParamStruct.ran_gauss_sigma_cut (0D_NOT_real -
      .def_property(
          "ran_gauss_sigma_cut",
          &BbuParamStruct::ran_gauss_sigma_cut,
          &BbuParamStruct::set_ran_gauss_sigma_cut
      )
      // BbuParamStruct.ele_track_end (0D_NOT_character -
      .def_property(
          "ele_track_end",
          &BbuParamStruct::ele_track_end,
          &BbuParamStruct::set_ele_track_end
      )
      // BbuParamStruct.ix_ele_track_end (0D_NOT_integer - Default: set to last element with a wake
      .def_property(
          "ix_ele_track_end",
          &BbuParamStruct::ix_ele_track_end,
          &BbuParamStruct::set_ix_ele_track_end
      )
      // BbuParamStruct.regression (0D_NOT_logical - Do regression test?
      .def_property("regression", &BbuParamStruct::regression, &BbuParamStruct::set_regression)
      // BbuParamStruct.normalize_z_to_rf (0D_NOT_logical - make starting z = mod(z, rf_wavelength)?
      // Ramp parameters
      .def_property(
          "normalize_z_to_rf",
          &BbuParamStruct::normalize_z_to_rf,
          &BbuParamStruct::set_normalize_z_to_rf
      )
      // BbuParamStruct.ramp_on (0D_NOT_logical -
      .def_property("ramp_on", &BbuParamStruct::ramp_on, &BbuParamStruct::set_ramp_on)
      // BbuParamStruct.ramp_pattern (1D_NOT_real -
      .def_property_readonly("ramp_pattern", &BbuParamStruct::ramp_pattern)
      // BbuParamStruct.ramp_n_start (0D_NOT_integer - Index of start of ramp Internal parameters
      .def_property(
          "ramp_n_start",
          &BbuParamStruct::ramp_n_start,
          &BbuParamStruct::set_ramp_n_start
      )
      // BbuParamStruct.n_ramp_pattern (0D_NOT_integer - Number of valid ramp_pattern
      .def_property(
          "n_ramp_pattern",
          &BbuParamStruct::n_ramp_pattern,
          &BbuParamStruct::set_n_ramp_pattern
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
  cls.def(py::init<>())
      // BbuStageStruct.ix_ele_lr_wake (0D_NOT_integer - Element index of element with the wake
      .def_property(
          "ix_ele_lr_wake",
          &BbuStageStruct::ix_ele_lr_wake,
          &BbuStageStruct::set_ix_ele_lr_wake
      )
      // BbuStageStruct.ix_ele_stage_end (0D_NOT_integer - Element at end of stage.
      .def_property(
          "ix_ele_stage_end",
          &BbuStageStruct::ix_ele_stage_end,
          &BbuStageStruct::set_ix_ele_stage_end
      )
      // BbuStageStruct.ix_pass (0D_NOT_integer - Pass index when in multipass section
      .def_property("ix_pass", &BbuStageStruct::ix_pass, &BbuStageStruct::set_ix_pass)
      // BbuStageStruct.ix_stage_pass1 (0D_NOT_integer - Index of corresponding stage on first pass
      .def_property(
          "ix_stage_pass1",
          &BbuStageStruct::ix_stage_pass1,
          &BbuStageStruct::set_ix_stage_pass1
      )
      // BbuStageStruct.ix_head_bunch (0D_NOT_integer -
      .def_property(
          "ix_head_bunch",
          &BbuStageStruct::ix_head_bunch,
          &BbuStageStruct::set_ix_head_bunch
      )
      // BbuStageStruct.ix_hom_max (0D_NOT_integer -
      .def_property("ix_hom_max", &BbuStageStruct::ix_hom_max, &BbuStageStruct::set_ix_hom_max)
      // BbuStageStruct.hom_voltage_max (0D_NOT_real -
      .def_property(
          "hom_voltage_max",
          &BbuStageStruct::hom_voltage_max,
          &BbuStageStruct::set_hom_voltage_max
      )
      // BbuStageStruct.time_at_wake_ele (0D_NOT_real -
      .def_property(
          "time_at_wake_ele",
          &BbuStageStruct::time_at_wake_ele,
          &BbuStageStruct::set_time_at_wake_ele
      )
      // BbuStageStruct.ave_orb (1D_NOT_real -
      .def_property_readonly("ave_orb", &BbuStageStruct::ave_orb)
      // BbuStageStruct.rms_orb (1D_NOT_real -
      .def_property_readonly("rms_orb", &BbuStageStruct::rms_orb)
      // BbuStageStruct.min_orb (1D_NOT_real -
      .def_property_readonly("min_orb", &BbuStageStruct::min_orb)
      // BbuStageStruct.max_orb (1D_NOT_real -
      .def_property_readonly("max_orb", &BbuStageStruct::max_orb)
      // BbuStageStruct.n_orb (0D_NOT_integer -
      .def_property("n_orb", &BbuStageStruct::n_orb, &BbuStageStruct::set_n_orb)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BbuStageStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
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

  bind_FTypeArrayND<BbuStageStructArray1D>(m, "BbuStageStructArray1D");
  bind_FTypeAlloc1D<BbuStageStructAlloc1D>(m, "BbuStageStructAlloc1D");
  // 2D BbuStageStruct arrays are not used in structs/routines
  // 3D BbuStageStruct arrays are not used in structs/routines
}

// =============================================================================
// beam_init_struct
void init_beam_init_struct(py::module &m, py::class_<BeamInitStruct> &cls) {
  cls.def(py::init<>())
      // BeamInitStruct.position_file (0D_NOT_character - File with particle positions.
      .def_property(
          "position_file",
          &BeamInitStruct::position_file,
          &BeamInitStruct::set_position_file
      )
      // BeamInitStruct.distribution_type (1D_NOT_character - distribution type (in x-px, y-py, and
      // z-pz planes) 'ELLIPSE', 'KV', 'GRID', 'FILE', 'RAN_GAUSS' or '' = 'RAN_GAUSS'
      .def_property_readonly("distribution_type", &BeamInitStruct::distribution_type)
      // BeamInitStruct.spin (1D_NOT_real - Spin (x, y, z)
      .def_property_readonly("spin", &BeamInitStruct::spin)
      // BeamInitStruct.ellipse (1D_NOT_type - Ellipse beam distribution
      .def_property_readonly("ellipse", &BeamInitStruct::ellipse)
      // BeamInitStruct.KV (0D_NOT_type - KV beam distribution
      .def_property("KV", &BeamInitStruct::KV, &BeamInitStruct::set_KV)
      // BeamInitStruct.grid (1D_NOT_type - Grid beam distribution
      .def_property_readonly("grid", &BeamInitStruct::grid)
      // BeamInitStruct.center_jitter (1D_NOT_real - Bunch center rms jitter
      .def_property_readonly("center_jitter", &BeamInitStruct::center_jitter)
      // BeamInitStruct.emit_jitter (1D_NOT_real - a and b bunch emittance rms jitter normalized to
      // emittance
      .def_property_readonly("emit_jitter", &BeamInitStruct::emit_jitter)
      // BeamInitStruct.sig_z_jitter (0D_NOT_real - bunch length RMS jitter
      .def_property(
          "sig_z_jitter",
          &BeamInitStruct::sig_z_jitter,
          &BeamInitStruct::set_sig_z_jitter
      )
      // BeamInitStruct.sig_pz_jitter (0D_NOT_real - RMS pz spread jitter
      .def_property(
          "sig_pz_jitter",
          &BeamInitStruct::sig_pz_jitter,
          &BeamInitStruct::set_sig_pz_jitter
      )
      // BeamInitStruct.n_particle (0D_NOT_integer - Number of particles per bunch.
      .def_property("n_particle", &BeamInitStruct::n_particle, &BeamInitStruct::set_n_particle)
      // BeamInitStruct.renorm_center (0D_NOT_logical - Renormalize centroid?
      .def_property(
          "renorm_center",
          &BeamInitStruct::renorm_center,
          &BeamInitStruct::set_renorm_center
      )
      // BeamInitStruct.renorm_sigma (0D_NOT_logical - Renormalize sigma?
      .def_property(
          "renorm_sigma",
          &BeamInitStruct::renorm_sigma,
          &BeamInitStruct::set_renorm_sigma
      )
      // BeamInitStruct.random_engine (0D_NOT_character - Or 'quasi'. Random number engine to use.
      .def_property(
          "random_engine",
          &BeamInitStruct::random_engine,
          &BeamInitStruct::set_random_engine
      )
      // BeamInitStruct.random_gauss_converter (0D_NOT_character - Or 'quick'. Uniform to gauss
      // conversion method.
      .def_property(
          "random_gauss_converter",
          &BeamInitStruct::random_gauss_converter,
          &BeamInitStruct::set_random_gauss_converter
      )
      // BeamInitStruct.random_sigma_cutoff (0D_NOT_real - Cut-off in sigmas.
      .def_property(
          "random_sigma_cutoff",
          &BeamInitStruct::random_sigma_cutoff,
          &BeamInitStruct::set_random_sigma_cutoff
      )
      // BeamInitStruct.a_norm_emit (0D_NOT_real - a-mode normalized emittance (emit * beta * gamma)
      .def_property("a_norm_emit", &BeamInitStruct::a_norm_emit, &BeamInitStruct::set_a_norm_emit)
      // BeamInitStruct.b_norm_emit (0D_NOT_real - b-mode normalized emittance (emit * beta * gamma)
      .def_property("b_norm_emit", &BeamInitStruct::b_norm_emit, &BeamInitStruct::set_b_norm_emit)
      // BeamInitStruct.a_emit (0D_NOT_real - a-mode emittance
      .def_property("a_emit", &BeamInitStruct::a_emit, &BeamInitStruct::set_a_emit)
      // BeamInitStruct.b_emit (0D_NOT_real - b-mode emittance
      .def_property("b_emit", &BeamInitStruct::b_emit, &BeamInitStruct::set_b_emit)
      // BeamInitStruct.dPz_dz (0D_NOT_real - Correlation of Pz with long position.
      .def_property("dPz_dz", &BeamInitStruct::dPz_dz, &BeamInitStruct::set_dPz_dz)
      // BeamInitStruct.center (1D_NOT_real - Bench phase space center offset relative to reference.
      .def_property_readonly("center", &BeamInitStruct::center)
      // BeamInitStruct.t_offset (0D_NOT_real - Time center offset
      .def_property("t_offset", &BeamInitStruct::t_offset, &BeamInitStruct::set_t_offset)
      // BeamInitStruct.dt_bunch (0D_NOT_real - Time between bunches.
      .def_property("dt_bunch", &BeamInitStruct::dt_bunch, &BeamInitStruct::set_dt_bunch)
      // BeamInitStruct.sig_z (0D_NOT_real - Z sigma in m.
      .def_property("sig_z", &BeamInitStruct::sig_z, &BeamInitStruct::set_sig_z)
      // BeamInitStruct.sig_pz (0D_NOT_real - pz sigma
      .def_property("sig_pz", &BeamInitStruct::sig_pz, &BeamInitStruct::set_sig_pz)
      // BeamInitStruct.bunch_charge (0D_NOT_real - charge (Coul) in a bunch.
      .def_property(
          "bunch_charge",
          &BeamInitStruct::bunch_charge,
          &BeamInitStruct::set_bunch_charge
      )
      // BeamInitStruct.n_bunch (0D_NOT_integer - Number of bunches.
      .def_property("n_bunch", &BeamInitStruct::n_bunch, &BeamInitStruct::set_n_bunch)
      // BeamInitStruct.ix_turn (0D_NOT_integer - Turn index used to adjust particles time if
      // needed.
      .def_property("ix_turn", &BeamInitStruct::ix_turn, &BeamInitStruct::set_ix_turn)
      // BeamInitStruct.species (0D_NOT_character - 'positron', etc. '' => use referece particle.
      .def_property("species", &BeamInitStruct::species, &BeamInitStruct::set_species)
      // BeamInitStruct.full_6D_coupling_calc (0D_NOT_logical - Use V from 6x6 1-turn mat to match
      // distribution? Else use 4x4 1-turn mat used.
      .def_property(
          "full_6D_coupling_calc",
          &BeamInitStruct::full_6D_coupling_calc,
          &BeamInitStruct::set_full_6D_coupling_calc
      )
      // BeamInitStruct.use_particle_start (0D_NOT_logical - Use lat%particle_start instead of
      // beam_init%center, %t_offset, and %spin?
      .def_property(
          "use_particle_start",
          &BeamInitStruct::use_particle_start,
          &BeamInitStruct::set_use_particle_start
      )
      // BeamInitStruct.use_t_coords (0D_NOT_logical - If true, the distributions will be taken as
      // in t-coordinates
      .def_property(
          "use_t_coords",
          &BeamInitStruct::use_t_coords,
          &BeamInitStruct::set_use_t_coords
      )
      // BeamInitStruct.use_z_as_t (0D_NOT_logical - Only used if  use_t_coords = .true. If true,  z
      // describes the t distribution If false, z describes the s distribution
      .def_property("use_z_as_t", &BeamInitStruct::use_z_as_t, &BeamInitStruct::set_use_z_as_t)
      // BeamInitStruct.file_name (0D_NOT_character - OLD!! DO NOT USE!!
      .def_property("file_name", &BeamInitStruct::file_name, &BeamInitStruct::set_file_name)

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
      // BeamStruct.bunch (1D_ALLOC_type -
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
  cls.def(py::init<>())
      // BmadCommonStruct.max_aperture_limit (0D_NOT_real - Max Aperture.
      .def_property(
          "max_aperture_limit",
          &BmadCommonStruct::max_aperture_limit,
          &BmadCommonStruct::set_max_aperture_limit
      )
      // BmadCommonStruct.d_orb (1D_NOT_real - Orbit deltas for the mat6 via tracking calc.
      .def_property_readonly("d_orb", &BmadCommonStruct::d_orb)
      // BmadCommonStruct.default_ds_step (0D_NOT_real - Default integration step for eles without
      // an explicit step calc.
      .def_property(
          "default_ds_step",
          &BmadCommonStruct::default_ds_step,
          &BmadCommonStruct::set_default_ds_step
      )
      // BmadCommonStruct.significant_length (0D_NOT_real - meter
      .def_property(
          "significant_length",
          &BmadCommonStruct::significant_length,
          &BmadCommonStruct::set_significant_length
      )
      // BmadCommonStruct.rel_tol_tracking (0D_NOT_real - Closed orbit relative tolerance.
      .def_property(
          "rel_tol_tracking",
          &BmadCommonStruct::rel_tol_tracking,
          &BmadCommonStruct::set_rel_tol_tracking
      )
      // BmadCommonStruct.abs_tol_tracking (0D_NOT_real - Closed orbit absolute tolerance.
      .def_property(
          "abs_tol_tracking",
          &BmadCommonStruct::abs_tol_tracking,
          &BmadCommonStruct::set_abs_tol_tracking
      )
      // BmadCommonStruct.rel_tol_adaptive_tracking (0D_NOT_real - Runge-Kutta tracking relative
      // tolerance.
      .def_property(
          "rel_tol_adaptive_tracking",
          &BmadCommonStruct::rel_tol_adaptive_tracking,
          &BmadCommonStruct::set_rel_tol_adaptive_tracking
      )
      // BmadCommonStruct.abs_tol_adaptive_tracking (0D_NOT_real - Runge-Kutta tracking absolute
      // tolerance.
      .def_property(
          "abs_tol_adaptive_tracking",
          &BmadCommonStruct::abs_tol_adaptive_tracking,
          &BmadCommonStruct::set_abs_tol_adaptive_tracking
      )
      // BmadCommonStruct.init_ds_adaptive_tracking (0D_NOT_real - Initial step size
      .def_property(
          "init_ds_adaptive_tracking",
          &BmadCommonStruct::init_ds_adaptive_tracking,
          &BmadCommonStruct::set_init_ds_adaptive_tracking
      )
      // BmadCommonStruct.min_ds_adaptive_tracking (0D_NOT_real - Min step size to take.
      .def_property(
          "min_ds_adaptive_tracking",
          &BmadCommonStruct::min_ds_adaptive_tracking,
          &BmadCommonStruct::set_min_ds_adaptive_tracking
      )
      // BmadCommonStruct.fatal_ds_adaptive_tracking (0D_NOT_real - If actual step size is below
      // this particle is lost.
      .def_property(
          "fatal_ds_adaptive_tracking",
          &BmadCommonStruct::fatal_ds_adaptive_tracking,
          &BmadCommonStruct::set_fatal_ds_adaptive_tracking
      )
      // BmadCommonStruct.autoscale_amp_abs_tol (0D_NOT_real - Autoscale absolute amplitude
      // tolerance (eV).
      .def_property(
          "autoscale_amp_abs_tol",
          &BmadCommonStruct::autoscale_amp_abs_tol,
          &BmadCommonStruct::set_autoscale_amp_abs_tol
      )
      // BmadCommonStruct.autoscale_amp_rel_tol (0D_NOT_real - Autoscale relative amplitude
      // tolerance
      .def_property(
          "autoscale_amp_rel_tol",
          &BmadCommonStruct::autoscale_amp_rel_tol,
          &BmadCommonStruct::set_autoscale_amp_rel_tol
      )
      // BmadCommonStruct.autoscale_phase_tol (0D_NOT_real - Autoscale phase tolerance.
      .def_property(
          "autoscale_phase_tol",
          &BmadCommonStruct::autoscale_phase_tol,
          &BmadCommonStruct::set_autoscale_phase_tol
      )
      // BmadCommonStruct.electric_dipole_moment (0D_NOT_real - Particle's EDM. Call set_ptc to
      // transfer value to PTC.
      .def_property(
          "electric_dipole_moment",
          &BmadCommonStruct::electric_dipole_moment,
          &BmadCommonStruct::set_electric_dipole_moment
      )
      // BmadCommonStruct.synch_rad_scale (0D_NOT_real - Synch radiation kick scale. 1 => normal, 0
      // => no kicks.
      .def_property(
          "synch_rad_scale",
          &BmadCommonStruct::synch_rad_scale,
          &BmadCommonStruct::set_synch_rad_scale
      )
      // BmadCommonStruct.sad_eps_scale (0D_NOT_real - Used in sad_mult step length calc.
      .def_property(
          "sad_eps_scale",
          &BmadCommonStruct::sad_eps_scale,
          &BmadCommonStruct::set_sad_eps_scale
      )
      // BmadCommonStruct.sad_amp_max (0D_NOT_real - Used in sad_mult step length calc.
      .def_property(
          "sad_amp_max",
          &BmadCommonStruct::sad_amp_max,
          &BmadCommonStruct::set_sad_amp_max
      )
      // BmadCommonStruct.sad_n_div_max (0D_NOT_integer - Used in sad_mult step length calc.
      .def_property(
          "sad_n_div_max",
          &BmadCommonStruct::sad_n_div_max,
          &BmadCommonStruct::set_sad_n_div_max
      )
      // BmadCommonStruct.taylor_order (0D_NOT_integer - Taylor order to use. 0 -> default =
      // ptc_private%taylor_order_saved.
      .def_property(
          "taylor_order",
          &BmadCommonStruct::taylor_order,
          &BmadCommonStruct::set_taylor_order
      )
      // BmadCommonStruct.runge_kutta_order (0D_NOT_integer - Runge Kutta order.
      .def_property(
          "runge_kutta_order",
          &BmadCommonStruct::runge_kutta_order,
          &BmadCommonStruct::set_runge_kutta_order
      )
      // BmadCommonStruct.default_integ_order (0D_NOT_integer - PTC integration order.
      .def_property(
          "default_integ_order",
          &BmadCommonStruct::default_integ_order,
          &BmadCommonStruct::set_default_integ_order
      )
      // BmadCommonStruct.max_num_runge_kutta_step (0D_NOT_integer - Maximum number of RK steps
      // before particle is considered lost.
      .def_property(
          "max_num_runge_kutta_step",
          &BmadCommonStruct::max_num_runge_kutta_step,
          &BmadCommonStruct::set_max_num_runge_kutta_step
      )
      // BmadCommonStruct.rf_phase_below_transition_ref (0D_NOT_logical - Autoscale uses below
      // transition stable point for RFCavities?
      .def_property(
          "rf_phase_below_transition_ref",
          &BmadCommonStruct::rf_phase_below_transition_ref,
          &BmadCommonStruct::set_rf_phase_below_transition_ref
      )
      // BmadCommonStruct.sr_wakes_on (0D_NOT_logical - Short range wakefields?
      .def_property(
          "sr_wakes_on",
          &BmadCommonStruct::sr_wakes_on,
          &BmadCommonStruct::set_sr_wakes_on
      )
      // BmadCommonStruct.lr_wakes_on (0D_NOT_logical - Long range wakefields
      .def_property(
          "lr_wakes_on",
          &BmadCommonStruct::lr_wakes_on,
          &BmadCommonStruct::set_lr_wakes_on
      )
      // BmadCommonStruct.auto_bookkeeper (0D_NOT_logical - Deprecated and no longer used.
      .def_property(
          "auto_bookkeeper",
          &BmadCommonStruct::auto_bookkeeper,
          &BmadCommonStruct::set_auto_bookkeeper
      )
      // BmadCommonStruct.high_energy_space_charge_on (0D_NOT_logical - High energy space charge
      // effect switch.
      .def_property(
          "high_energy_space_charge_on",
          &BmadCommonStruct::high_energy_space_charge_on,
          &BmadCommonStruct::set_high_energy_space_charge_on
      )
      // BmadCommonStruct.csr_and_space_charge_on (0D_NOT_logical - Space charge switch.
      .def_property(
          "csr_and_space_charge_on",
          &BmadCommonStruct::csr_and_space_charge_on,
          &BmadCommonStruct::set_csr_and_space_charge_on
      )
      // BmadCommonStruct.spin_tracking_on (0D_NOT_logical - spin tracking?
      .def_property(
          "spin_tracking_on",
          &BmadCommonStruct::spin_tracking_on,
          &BmadCommonStruct::set_spin_tracking_on
      )
      // BmadCommonStruct.spin_sokolov_ternov_flipping_on (0D_NOT_logical - Spin flipping during
      // synchrotron radiation emission?
      .def_property(
          "spin_sokolov_ternov_flipping_on",
          &BmadCommonStruct::spin_sokolov_ternov_flipping_on,
          &BmadCommonStruct::set_spin_sokolov_ternov_flipping_on
      )
      // BmadCommonStruct.radiation_damping_on (0D_NOT_logical - Radiation damping toggle.
      .def_property(
          "radiation_damping_on",
          &BmadCommonStruct::radiation_damping_on,
          &BmadCommonStruct::set_radiation_damping_on
      )
      // BmadCommonStruct.radiation_zero_average (0D_NOT_logical - Shift damping to be zero on the
      // zero orbit to get rid of sawtooth?
      .def_property(
          "radiation_zero_average",
          &BmadCommonStruct::radiation_zero_average,
          &BmadCommonStruct::set_radiation_zero_average
      )
      // BmadCommonStruct.radiation_fluctuations_on (0D_NOT_logical - Radiation fluctuations toggle.
      .def_property(
          "radiation_fluctuations_on",
          &BmadCommonStruct::radiation_fluctuations_on,
          &BmadCommonStruct::set_radiation_fluctuations_on
      )
      // BmadCommonStruct.conserve_taylor_maps (0D_NOT_logical - Enable bookkeeper to set
      // ele%taylor_map_includes_offsets = F?
      .def_property(
          "conserve_taylor_maps",
          &BmadCommonStruct::conserve_taylor_maps,
          &BmadCommonStruct::set_conserve_taylor_maps
      )
      // BmadCommonStruct.absolute_time_tracking (0D_NOT_logical - Absolute or relative time
      // tracking?
      .def_property(
          "absolute_time_tracking",
          &BmadCommonStruct::absolute_time_tracking,
          &BmadCommonStruct::set_absolute_time_tracking
      )
      // BmadCommonStruct.absolute_time_ref_shift (0D_NOT_logical - Apply reference time shift when
      // using absolute time tracking?
      .def_property(
          "absolute_time_ref_shift",
          &BmadCommonStruct::absolute_time_ref_shift,
          &BmadCommonStruct::set_absolute_time_ref_shift
      )
      // BmadCommonStruct.convert_to_kinetic_momentum (0D_NOT_logical - Cancel kicks due to finite
      // vector potential when doing symplectic tracking? Set to True to test symp_lie_bmad against
      // runge_kutta.
      .def_property(
          "convert_to_kinetic_momentum",
          &BmadCommonStruct::convert_to_kinetic_momentum,
          &BmadCommonStruct::set_convert_to_kinetic_momentum
      )
      // BmadCommonStruct.normalize_twiss (0D_NOT_logical - Normalize matrix when computing Twiss
      // for off-energy ref?
      .def_property(
          "normalize_twiss",
          &BmadCommonStruct::normalize_twiss,
          &BmadCommonStruct::set_normalize_twiss
      )
      // BmadCommonStruct.aperture_limit_on (0D_NOT_logical - Use apertures in tracking?
      .def_property(
          "aperture_limit_on",
          &BmadCommonStruct::aperture_limit_on,
          &BmadCommonStruct::set_aperture_limit_on
      )
      // BmadCommonStruct.spin_n0_direction_user_set (0D_NOT_logical - User sets direction of n0 for
      // closed geometry branches?
      .def_property(
          "spin_n0_direction_user_set",
          &BmadCommonStruct::spin_n0_direction_user_set,
          &BmadCommonStruct::set_spin_n0_direction_user_set
      )
      // BmadCommonStruct.debug (0D_NOT_logical - Used for code debugging.
      .def_property("debug", &BmadCommonStruct::debug, &BmadCommonStruct::set_debug)

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
  cls.def(py::init<>())
      // BmadNormalFormStruct.ele_origin (0D_PTR_type - Element at which the on-turn map was
      // created.
      .def_property(
          "ele_origin",
          &BmadNormalFormStruct::ele_origin,
          &BmadNormalFormStruct::set_ele_origin
      )
      // BmadNormalFormStruct.M (1D_NOT_type - One-turn taylor map: M = A o N o A_inv, N = exp(:h:)
      .def_property_readonly("M", &BmadNormalFormStruct::M)
      // BmadNormalFormStruct.A (1D_NOT_type - Map from Floquet -> Lab coordinates
      .def_property_readonly("A", &BmadNormalFormStruct::A)
      // BmadNormalFormStruct.A_inv (1D_NOT_type - Map from Lab -> Floquet coordinates
      .def_property_readonly("A_inv", &BmadNormalFormStruct::A_inv)
      // BmadNormalFormStruct.dhdj (1D_NOT_type - Nonlinear tune function operating on Floquet
      // coordinates
      .def_property_readonly("dhdj", &BmadNormalFormStruct::dhdj)
      // BmadNormalFormStruct.F (1D_NOT_type - Vector field factorization in phasor basis:
      .def_property_readonly("F", &BmadNormalFormStruct::F)
      // BmadNormalFormStruct.L (1D_NOT_type - L component
      .def_property_readonly("L", &BmadNormalFormStruct::L)
      // BmadNormalFormStruct.h (1D_ALLOC_type -
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
  cls.def(py::init<>())
      // BookkeepingStateStruct.attributes (0D_NOT_integer - Element dependent attributes:
      // super_ok$, ok$ or stale$
      .def_property(
          "attributes",
          &BookkeepingStateStruct::attributes,
          &BookkeepingStateStruct::set_attributes
      )
      // BookkeepingStateStruct.control (0D_NOT_integer - Lord/slave bookkeeping status: super_ok$,
      // ok$ or stale$
      .def_property(
          "control",
          &BookkeepingStateStruct::control,
          &BookkeepingStateStruct::set_control
      )
      // BookkeepingStateStruct.floor_position (0D_NOT_integer - Global (floor) geometry: super_ok$,
      // ok$ or stale$
      .def_property(
          "floor_position",
          &BookkeepingStateStruct::floor_position,
          &BookkeepingStateStruct::set_floor_position
      )
      // BookkeepingStateStruct.s_position (0D_NOT_integer - Longitudinal position & element length:
      // super_ok$, ok$ or stale$
      .def_property(
          "s_position",
          &BookkeepingStateStruct::s_position,
          &BookkeepingStateStruct::set_s_position
      )
      // BookkeepingStateStruct.ref_energy (0D_NOT_integer - Reference energy and ref time:
      // super_ok$, ok$ or stale$
      .def_property(
          "ref_energy",
          &BookkeepingStateStruct::ref_energy,
          &BookkeepingStateStruct::set_ref_energy
      )
      // BookkeepingStateStruct.mat6 (0D_NOT_integer - Linear transfer map status: super_ok$, ok$ or
      // stale$
      .def_property("mat6", &BookkeepingStateStruct::mat6, &BookkeepingStateStruct::set_mat6)
      // BookkeepingStateStruct.rad_int (0D_NOT_integer - Radiation integrals cache status
      .def_property(
          "rad_int",
          &BookkeepingStateStruct::rad_int,
          &BookkeepingStateStruct::set_rad_int
      )
      // BookkeepingStateStruct.ptc (0D_NOT_integer - Associated PTC fibre (or layout) status.
      .def_property("ptc", &BookkeepingStateStruct::ptc, &BookkeepingStateStruct::set_ptc)
      // BookkeepingStateStruct.has_misalign (0D_NOT_logical - Used to avoid unnecessary calls to
      // offset_particle.
      .def_property(
          "has_misalign",
          &BookkeepingStateStruct::has_misalign,
          &BookkeepingStateStruct::set_has_misalign
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
  cls.def(py::init<>())
      // BpmPhaseCouplingStruct.K_22a (0D_NOT_real - In-phase y/x for a-mode oscillations.
      .def_property("K_22a", &BpmPhaseCouplingStruct::K_22a, &BpmPhaseCouplingStruct::set_K_22a)
      // BpmPhaseCouplingStruct.K_12a (0D_NOT_real - Out-of-phase y/x for a-mode oscillations.
      .def_property("K_12a", &BpmPhaseCouplingStruct::K_12a, &BpmPhaseCouplingStruct::set_K_12a)
      // BpmPhaseCouplingStruct.K_11b (0D_NOT_real - In-phase x/y for b-mode oscillations.
      .def_property("K_11b", &BpmPhaseCouplingStruct::K_11b, &BpmPhaseCouplingStruct::set_K_11b)
      // BpmPhaseCouplingStruct.K_12b (0D_NOT_real - Out-of-phase x/y for b-mode oscillations.
      .def_property("K_12b", &BpmPhaseCouplingStruct::K_12b, &BpmPhaseCouplingStruct::set_K_12b)
      // BpmPhaseCouplingStruct.Cbar22_a (0D_NOT_real - Cbar22 as calculated from K_22a.
      .def_property(
          "Cbar22_a",
          &BpmPhaseCouplingStruct::Cbar22_a,
          &BpmPhaseCouplingStruct::set_Cbar22_a
      )
      // BpmPhaseCouplingStruct.Cbar12_a (0D_NOT_real - Cbar12 as calculated from K_12a.
      .def_property(
          "Cbar12_a",
          &BpmPhaseCouplingStruct::Cbar12_a,
          &BpmPhaseCouplingStruct::set_Cbar12_a
      )
      // BpmPhaseCouplingStruct.Cbar11_b (0D_NOT_real - Cbar11 as calculated from K_11b.
      .def_property(
          "Cbar11_b",
          &BpmPhaseCouplingStruct::Cbar11_b,
          &BpmPhaseCouplingStruct::set_Cbar11_b
      )
      // BpmPhaseCouplingStruct.Cbar12_b (0D_NOT_real - Cbar12 as calculated from K_12b.
      .def_property(
          "Cbar12_b",
          &BpmPhaseCouplingStruct::Cbar12_b,
          &BpmPhaseCouplingStruct::set_Cbar12_b
      )
      // BpmPhaseCouplingStruct.phi_a (0D_NOT_real - a-mode betatron phase.
      .def_property("phi_a", &BpmPhaseCouplingStruct::phi_a, &BpmPhaseCouplingStruct::set_phi_a)
      // BpmPhaseCouplingStruct.phi_b (0D_NOT_real - b-mode betatron phase.
      .def_property("phi_b", &BpmPhaseCouplingStruct::phi_b, &BpmPhaseCouplingStruct::set_phi_b)

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
  cls.def(py::init<>())
      // BranchStruct.name (0D_NOT_character - Name of line that defines the branch.
      .def_property("name", &BranchStruct::name, &BranchStruct::set_name)
      // BranchStruct.ix_branch (0D_NOT_integer - Index of this branch. 0 => Main branch
      .def_property("ix_branch", &BranchStruct::ix_branch, &BranchStruct::set_ix_branch)
      // BranchStruct.ix_from_branch (0D_NOT_integer - -1 => No creating fork element to this
      // branch.
      .def_property(
          "ix_from_branch",
          &BranchStruct::ix_from_branch,
          &BranchStruct::set_ix_from_branch
      )
      // BranchStruct.ix_from_ele (0D_NOT_integer - Index of creating fork element which forks to
      // this branch.
      .def_property("ix_from_ele", &BranchStruct::ix_from_ele, &BranchStruct::set_ix_from_ele)
      // BranchStruct.ix_to_ele (0D_NOT_integer - Index of element in this branch that creating fork
      // element forks to.
      .def_property("ix_to_ele", &BranchStruct::ix_to_ele, &BranchStruct::set_ix_to_ele)
      // BranchStruct.ix_fixer (0D_NOT_integer - Index of active fixer or beginning_ele element.
      .def_property("ix_fixer", &BranchStruct::ix_fixer, &BranchStruct::set_ix_fixer)
      // BranchStruct.n_ele_track (0D_NOT_integer -
      .def_property("n_ele_track", &BranchStruct::n_ele_track, &BranchStruct::set_n_ele_track)
      // BranchStruct.n_ele_max (0D_NOT_integer -
      .def_property("n_ele_max", &BranchStruct::n_ele_max, &BranchStruct::set_n_ele_max)
      // BranchStruct.lat (0D_PTR_type -
      .def_property("lat", &BranchStruct::lat, &BranchStruct::set_lat)
      // BranchStruct.a (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("a", &BranchStruct::a, &BranchStruct::set_a)
      // BranchStruct.b (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("b", &BranchStruct::b, &BranchStruct::set_b)
      // BranchStruct.z (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("z", &BranchStruct::z, &BranchStruct::set_z)
      // BranchStruct.ele (1D_PTR_type -
      .def_property_readonly("ele", &BranchStruct::ele)
      // BranchStruct.param (0D_NOT_type -
      .def_property("param", &BranchStruct::param, &BranchStruct::set_param)
      // BranchStruct.particle_start (0D_NOT_type -
      .def_property(
          "particle_start",
          &BranchStruct::particle_start,
          &BranchStruct::set_particle_start
      )
      // BranchStruct.wall3d (1D_PTR_type -
      .def_property_readonly("wall3d", &BranchStruct::wall3d)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BranchStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
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

  bind_FTypeArrayND<BranchStructArray1D>(m, "BranchStructArray1D");
  bind_FTypeAlloc1D<BranchStructAlloc1D>(m, "BranchStructAlloc1D");
  // 2D BranchStruct arrays are not used in structs/routines
  // 3D BranchStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_params_struct
void init_bunch_params_struct(py::module &m, py::class_<BunchParamsStruct> &cls) {
  cls.def(py::init<>())
      // BunchParamsStruct.centroid (0D_NOT_type - Lab frame
      .def_property("centroid", &BunchParamsStruct::centroid, &BunchParamsStruct::set_centroid)
      // BunchParamsStruct.x (0D_NOT_type - Projected Twiss parameters
      .def_property("x", &BunchParamsStruct::x, &BunchParamsStruct::set_x)
      // BunchParamsStruct.y (0D_NOT_type - Projected Twiss parameters
      .def_property("y", &BunchParamsStruct::y, &BunchParamsStruct::set_y)
      // BunchParamsStruct.z (0D_NOT_type - Projected Twiss parameters
      .def_property("z", &BunchParamsStruct::z, &BunchParamsStruct::set_z)
      // BunchParamsStruct.a (0D_NOT_type - Normal mode twiss parameters
      .def_property("a", &BunchParamsStruct::a, &BunchParamsStruct::set_a)
      // BunchParamsStruct.b (0D_NOT_type - Normal mode twiss parameters
      .def_property("b", &BunchParamsStruct::b, &BunchParamsStruct::set_b)
      // BunchParamsStruct.c (0D_NOT_type - Normal mode twiss parameters
      .def_property("c", &BunchParamsStruct::c, &BunchParamsStruct::set_c)
      // BunchParamsStruct.sigma (2D_NOT_real - beam size matrix
      .def_property_readonly("sigma", &BunchParamsStruct::sigma)
      // BunchParamsStruct.rel_max (1D_NOT_real - Max orbit relative to centroid. 7 -> time.
      .def_property_readonly("rel_max", &BunchParamsStruct::rel_max)
      // BunchParamsStruct.rel_min (1D_NOT_real - Min orbit relative to_centroid. 7 -> time.
      .def_property_readonly("rel_min", &BunchParamsStruct::rel_min)
      // BunchParamsStruct.s (0D_NOT_real - Longitudinal position.
      .def_property("s", &BunchParamsStruct::s, &BunchParamsStruct::set_s)
      // BunchParamsStruct.t (0D_NOT_real - Time.
      .def_property("t", &BunchParamsStruct::t, &BunchParamsStruct::set_t)
      // BunchParamsStruct.sigma_t (0D_NOT_real - RMS of time spread.
      .def_property("sigma_t", &BunchParamsStruct::sigma_t, &BunchParamsStruct::set_sigma_t)
      // BunchParamsStruct.charge_live (0D_NOT_real - Charge of all non-lost particle
      .def_property(
          "charge_live",
          &BunchParamsStruct::charge_live,
          &BunchParamsStruct::set_charge_live
      )
      // BunchParamsStruct.charge_tot (0D_NOT_real - Charge of all particles.
      .def_property(
          "charge_tot",
          &BunchParamsStruct::charge_tot,
          &BunchParamsStruct::set_charge_tot
      )
      // BunchParamsStruct.n_particle_tot (0D_NOT_integer - Total number of particles
      .def_property(
          "n_particle_tot",
          &BunchParamsStruct::n_particle_tot,
          &BunchParamsStruct::set_n_particle_tot
      )
      // BunchParamsStruct.n_particle_live (0D_NOT_integer - Number of non-lost particles
      .def_property(
          "n_particle_live",
          &BunchParamsStruct::n_particle_live,
          &BunchParamsStruct::set_n_particle_live
      )
      // BunchParamsStruct.n_particle_lost_in_ele (0D_NOT_integer - Number lost in element (not
      // calculated by Bmad)
      .def_property(
          "n_particle_lost_in_ele",
          &BunchParamsStruct::n_particle_lost_in_ele,
          &BunchParamsStruct::set_n_particle_lost_in_ele
      )
      // BunchParamsStruct.n_good_steps (0D_NOT_integer - Number of good steps (set when tracking
      // with space charge)
      .def_property(
          "n_good_steps",
          &BunchParamsStruct::n_good_steps,
          &BunchParamsStruct::set_n_good_steps
      )
      // BunchParamsStruct.n_bad_steps (0D_NOT_integer - Number of bad steps (set when tracking with
      // space charge)
      .def_property(
          "n_bad_steps",
          &BunchParamsStruct::n_bad_steps,
          &BunchParamsStruct::set_n_bad_steps
      )
      // BunchParamsStruct.ix_ele (0D_NOT_integer - Lattice element where params evaluated at.
      .def_property("ix_ele", &BunchParamsStruct::ix_ele, &BunchParamsStruct::set_ix_ele)
      // BunchParamsStruct.location (0D_NOT_integer - Location in element: upstream_end$, inside$,
      // or downstream_end$
      .def_property("location", &BunchParamsStruct::location, &BunchParamsStruct::set_location)
      // BunchParamsStruct.twiss_valid (0D_NOT_logical - Is the data here valid? Note: IF there is
      // no energy variation (RF off) twiss_valid may be true but in this case the z-twiss will not
      // be valid.
      .def_property(
          "twiss_valid",
          &BunchParamsStruct::twiss_valid,
          &BunchParamsStruct::set_twiss_valid
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BunchParamsStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
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

  bind_FTypeArrayND<BunchParamsStructArray1D>(m, "BunchParamsStructArray1D");
  bind_FTypeAlloc1D<BunchParamsStructAlloc1D>(m, "BunchParamsStructAlloc1D");
  // 2D BunchParamsStruct arrays are not used in structs/routines
  // 3D BunchParamsStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_struct
void init_bunch_struct(py::module &m, py::class_<BunchStruct> &cls) {
  cls.def(py::init<>())
      // BunchStruct.particle (1D_ALLOC_type -
      .def_property_readonly("particle", &BunchStruct::particle)
      // BunchStruct.ix_z (1D_ALLOC_integer - bunch%ix_z(1) is index of head particle, etc.
      .def_property_readonly("ix_z", &BunchStruct::ix_z)
      // BunchStruct.charge_tot (0D_NOT_real - Total charge in a bunch (Coul).
      .def_property("charge_tot", &BunchStruct::charge_tot, &BunchStruct::set_charge_tot)
      // BunchStruct.charge_live (0D_NOT_real - Charge of live particles (Coul).
      .def_property("charge_live", &BunchStruct::charge_live, &BunchStruct::set_charge_live)
      // BunchStruct.z_center (0D_NOT_real - Longitudinal center of bunch at creation time. Note:
      // Generally, z_center of bunch #1 is 0 and z_center of the other bunches is negative.
      .def_property("z_center", &BunchStruct::z_center, &BunchStruct::set_z_center)
      // BunchStruct.t_center (0D_NOT_real - Center of bunch at creation time relative to head
      // bunch.
      .def_property("t_center", &BunchStruct::t_center, &BunchStruct::set_t_center)
      // BunchStruct.t0 (0D_NOT_real - Used by track1_bunch_space_charge for tracking so particles
      // have constant t.
      .def_property("t0", &BunchStruct::t0, &BunchStruct::set_t0)
      // BunchStruct.drift_between_t_and_s (0D_NOT_logical - Drift (ignore any fields) instead of
      // tracking to speed up the calculation? This can only be done under certain circumstances.
      .def_property(
          "drift_between_t_and_s",
          &BunchStruct::drift_between_t_and_s,
          &BunchStruct::set_drift_between_t_and_s
      )
      // BunchStruct.ix_ele (0D_NOT_integer - Nominal element bunch is at. But, EG, dead particles
      // can be someplace else.
      .def_property("ix_ele", &BunchStruct::ix_ele, &BunchStruct::set_ix_ele)
      // BunchStruct.ix_bunch (0D_NOT_integer - Bunch index. Head bunch = 1, etc.
      .def_property("ix_bunch", &BunchStruct::ix_bunch, &BunchStruct::set_ix_bunch)
      // BunchStruct.ix_turn (0D_NOT_integer - Turn index for long term tracking. ix_turn = 0 before
      // end of first turn, etc.
      .def_property("ix_turn", &BunchStruct::ix_turn, &BunchStruct::set_ix_turn)
      // BunchStruct.n_live (0D_NOT_integer -
      .def_property("n_live", &BunchStruct::n_live, &BunchStruct::set_n_live)
      // BunchStruct.n_good (0D_NOT_integer - Number of accepted steps when using adaptive step size
      // control.
      .def_property("n_good", &BunchStruct::n_good, &BunchStruct::set_n_good)
      // BunchStruct.n_bad (0D_NOT_integer - Number of rejected steps when using adaptive step size
      // control.
      .def_property("n_bad", &BunchStruct::n_bad, &BunchStruct::set_n_bad)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BunchStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
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

  bind_FTypeArrayND<BunchStructArray1D>(m, "BunchStructArray1D");
  bind_FTypeAlloc1D<BunchStructAlloc1D>(m, "BunchStructAlloc1D");
  // 2D BunchStruct arrays are not used in structs/routines
  // 3D BunchStruct arrays are not used in structs/routines
}

// =============================================================================
// bunch_track_struct
void init_bunch_track_struct(py::module &m, py::class_<BunchTrackStruct> &cls) {
  cls.def(py::init<>())
      // BunchTrackStruct.pt (1D_ALLOC_type - Array indexed from 0
      .def_property_readonly("pt", &BunchTrackStruct::pt)
      // BunchTrackStruct.ds_save (0D_NOT_real - Min distance between points.
      .def_property("ds_save", &BunchTrackStruct::ds_save, &BunchTrackStruct::set_ds_save)
      // BunchTrackStruct.n_pt (0D_NOT_integer - Track upper bound
      .def_property("n_pt", &BunchTrackStruct::n_pt, &BunchTrackStruct::set_n_pt)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BunchTrackStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
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

  bind_FTypeArrayND<BunchTrackStructArray1D>(m, "BunchTrackStructArray1D");
  bind_FTypeAlloc1D<BunchTrackStructAlloc1D>(m, "BunchTrackStructAlloc1D");
  // 2D BunchTrackStruct arrays are not used in structs/routines
  // 3D BunchTrackStruct arrays are not used in structs/routines
}

// =============================================================================
// bicubic_cmplx_coef_struct
void init_bicubic_cmplx_coef_struct(py::module &m, py::class_<BicubicCmplxCoefStruct> &cls) {
  cls.def(py::init<>())
      // BicubicCmplxCoefStruct.coef (2D_NOT_complex - Coefs
      .def_property_readonly("coef", &BicubicCmplxCoefStruct::coef)
      // BicubicCmplxCoefStruct.i_box (1D_NOT_integer - index at lower box corner.
      .def_property_readonly("i_box", &BicubicCmplxCoefStruct::i_box)

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