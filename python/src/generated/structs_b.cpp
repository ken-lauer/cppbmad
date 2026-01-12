#include "pybmad/generated/structs_b.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// bbu_beam_struct
void init_bbu_beam_struct(py::module& m, py::class_<BbuBeamProxy>& cls) {
  cls.def(py::init<>())
      // BbuBeamProxy.bunch (1D_ALLOC_type - Bunches in the lattice
      .def_property_readonly("bunch", &BbuBeamProxy::bunch)
      // BbuBeamProxy.stage (1D_ALLOC_type -
      .def_property_readonly("stage", &BbuBeamProxy::stage)
      // BbuBeamProxy.ix_ele_bunch (1D_ALLOC_integer - element where bunch is
      .def_property_readonly("ix_ele_bunch", &BbuBeamProxy::ix_ele_bunch)
      // BbuBeamProxy.ix_bunch_head (0D_NOT_integer - Index to head bunch(:)
      .def_property(
          "ix_bunch_head",
          &BbuBeamProxy::ix_bunch_head,
          &BbuBeamProxy::set_ix_bunch_head)
      // BbuBeamProxy.ix_bunch_end (0D_NOT_integer - Index of the end bunch(:). -1 -> no bunches.
      .def_property(
          "ix_bunch_end",
          &BbuBeamProxy::ix_bunch_end,
          &BbuBeamProxy::set_ix_bunch_end)
      // BbuBeamProxy.n_bunch_in_lat (0D_NOT_integer - Number of bunches transversing the lattice.
      .def_property(
          "n_bunch_in_lat",
          &BbuBeamProxy::n_bunch_in_lat,
          &BbuBeamProxy::set_n_bunch_in_lat)
      // BbuBeamProxy.ix_stage_voltage_max (0D_NOT_integer -
      .def_property(
          "ix_stage_voltage_max",
          &BbuBeamProxy::ix_stage_voltage_max,
          &BbuBeamProxy::set_ix_stage_voltage_max)
      // BbuBeamProxy.hom_voltage_max (0D_NOT_real -
      .def_property(
          "hom_voltage_max",
          &BbuBeamProxy::hom_voltage_max,
          &BbuBeamProxy::set_hom_voltage_max)
      // BbuBeamProxy.time_now (0D_NOT_real -
      .def_property(
          "time_now", &BbuBeamProxy::time_now, &BbuBeamProxy::set_time_now)
      // BbuBeamProxy.one_turn_time (0D_NOT_real -
      .def_property(
          "one_turn_time",
          &BbuBeamProxy::one_turn_time,
          &BbuBeamProxy::set_one_turn_time)
      // BbuBeamProxy.rf_wavelength_max (0D_NOT_real -
      .def_property(
          "rf_wavelength_max",
          &BbuBeamProxy::rf_wavelength_max,
          &BbuBeamProxy::set_rf_wavelength_max)

      .def("__repr__", [](const BbuBeamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuBeamProxy& self) {
            return BbuBeamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BbuBeamProxy& self, py::dict& memo) {
            return BbuBeamProxy(self);
          })

      ;

  // 1D BbuBeamProxy arrays are not used in structs/routines
  // 2D BbuBeamProxy arrays are not used in structs/routines
  // 3D BbuBeamProxy arrays are not used in structs/routines
}

// =============================================================================
// bbu_param_struct
void init_bbu_param_struct(py::module& m, py::class_<BbuParamProxy>& cls) {
  cls.def(py::init<>())
      // BbuParamProxy.lat_filename (0D_NOT_character - Bmad lattice file name
      .def_property(
          "lat_filename",
          &BbuParamProxy::lat_filename,
          &BbuParamProxy::set_lat_filename)
      // BbuParamProxy.lat2_filename (0D_NOT_character - Bmad lattice2 file name for secondary parser
      .def_property(
          "lat2_filename",
          &BbuParamProxy::lat2_filename,
          &BbuParamProxy::set_lat2_filename)
      // BbuParamProxy.bunch_by_bunch_info_file (0D_NOT_character - For outputting bunch-by-bunch info.
      .def_property(
          "bunch_by_bunch_info_file",
          &BbuParamProxy::bunch_by_bunch_info_file,
          &BbuParamProxy::set_bunch_by_bunch_info_file)
      // BbuParamProxy.hybridize (0D_NOT_logical - Combine non-hom elements to speed up simulation?
      .def_property(
          "hybridize", &BbuParamProxy::hybridize, &BbuParamProxy::set_hybridize)
      // BbuParamProxy.write_digested_hybrid_lat (0D_NOT_logical - For debugging purposes.
      .def_property(
          "write_digested_hybrid_lat",
          &BbuParamProxy::write_digested_hybrid_lat,
          &BbuParamProxy::set_write_digested_hybrid_lat)
      // BbuParamProxy.write_voltage_vs_time_dat (0D_NOT_logical - For debugging purposes.
      .def_property(
          "write_voltage_vs_time_dat",
          &BbuParamProxy::write_voltage_vs_time_dat,
          &BbuParamProxy::set_write_voltage_vs_time_dat)
      // BbuParamProxy.keep_overlays_and_groups (0D_NOT_logical - Keep when hybridizing?
      .def_property(
          "keep_overlays_and_groups",
          &BbuParamProxy::keep_overlays_and_groups,
          &BbuParamProxy::set_keep_overlays_and_groups)
      // BbuParamProxy.keep_all_lcavities (0D_NOT_logical - Keep when hybridizing?
      .def_property(
          "keep_all_lcavities",
          &BbuParamProxy::keep_all_lcavities,
          &BbuParamProxy::set_keep_all_lcavities)
      // BbuParamProxy.use_taylor_for_hybrids (0D_NOT_logical - Use taylor map for hybrids when true. Otherwise tracking method is linear.
      .def_property(
          "use_taylor_for_hybrids",
          &BbuParamProxy::use_taylor_for_hybrids,
          &BbuParamProxy::set_use_taylor_for_hybrids)
      // BbuParamProxy.stable_orbit_anal (0D_NOT_logical - Write stable_orbit.out and hom_voltage.out?
      .def_property(
          "stable_orbit_anal",
          &BbuParamProxy::stable_orbit_anal,
          &BbuParamProxy::set_stable_orbit_anal)
      // BbuParamProxy.limit_factor (0D_NOT_real - Init_hom_amp * limit_factor = simulation unstable limit
      .def_property(
          "limit_factor",
          &BbuParamProxy::limit_factor,
          &BbuParamProxy::set_limit_factor)
      // BbuParamProxy.simulation_turns_max (0D_NOT_real - Sets the duration of the simulation.
      .def_property(
          "simulation_turns_max",
          &BbuParamProxy::simulation_turns_max,
          &BbuParamProxy::set_simulation_turns_max)
      // BbuParamProxy.bunch_freq (0D_NOT_real - Freq in Hz.
      .def_property(
          "bunch_freq",
          &BbuParamProxy::bunch_freq,
          &BbuParamProxy::set_bunch_freq)
      // BbuParamProxy.init_particle_offset (0D_NOT_real - Initial particle offset for particles born in the first turn period.
      .def_property(
          "init_particle_offset",
          &BbuParamProxy::init_particle_offset,
          &BbuParamProxy::set_init_particle_offset)
      // BbuParamProxy.current (0D_NOT_real - Starting current (amps)
      .def_property(
          "current", &BbuParamProxy::current, &BbuParamProxy::set_current)
      // BbuParamProxy.rel_tol (0D_NOT_real - Final threshold current accuracy.
      .def_property(
          "rel_tol", &BbuParamProxy::rel_tol, &BbuParamProxy::set_rel_tol)
      // BbuParamProxy.drscan (0D_NOT_logical - If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3.
      .def_property(
          "drscan", &BbuParamProxy::drscan, &BbuParamProxy::set_drscan)
      // BbuParamProxy.use_interpolated_threshold (0D_NOT_logical -
      .def_property(
          "use_interpolated_threshold",
          &BbuParamProxy::use_interpolated_threshold,
          &BbuParamProxy::set_use_interpolated_threshold)
      // BbuParamProxy.write_hom_info (0D_NOT_logical - Write HOM parameters to main output file?
      .def_property(
          "write_hom_info",
          &BbuParamProxy::write_hom_info,
          &BbuParamProxy::set_write_hom_info)
      // BbuParamProxy.elindex (0D_NOT_integer -
      .def_property(
          "elindex", &BbuParamProxy::elindex, &BbuParamProxy::set_elindex)
      // BbuParamProxy.elname (0D_NOT_character - Element to step length for DRSCAN
      .def_property(
          "elname", &BbuParamProxy::elname, &BbuParamProxy::set_elname)
      // BbuParamProxy.nstep (0D_NOT_integer - Number of steps for DRSCAN.
      .def_property("nstep", &BbuParamProxy::nstep, &BbuParamProxy::set_nstep)
      // BbuParamProxy.begdr (0D_NOT_real - Beginning DR value for DRSCAN.
      .def_property("begdr", &BbuParamProxy::begdr, &BbuParamProxy::set_begdr)
      // BbuParamProxy.enddr (0D_NOT_real - End DR value for DRSCAN.
      .def_property("enddr", &BbuParamProxy::enddr, &BbuParamProxy::set_enddr)
      // BbuParamProxy.nrep (0D_NOT_integer - Number of times to repeat threshold calculation
      .def_property("nrep", &BbuParamProxy::nrep, &BbuParamProxy::set_nrep)
      // BbuParamProxy.ran_seed (0D_NOT_integer - If set to 0, the output results will vary from run to run.
      .def_property(
          "ran_seed", &BbuParamProxy::ran_seed, &BbuParamProxy::set_ran_seed)
      // BbuParamProxy.hom_order_cutoff (0D_NOT_integer - If positive -> ignore HOM's with order greater than this.
      .def_property(
          "hom_order_cutoff",
          &BbuParamProxy::hom_order_cutoff,
          &BbuParamProxy::set_hom_order_cutoff)
      // BbuParamProxy.ran_gauss_sigma_cut (0D_NOT_real -
      .def_property(
          "ran_gauss_sigma_cut",
          &BbuParamProxy::ran_gauss_sigma_cut,
          &BbuParamProxy::set_ran_gauss_sigma_cut)
      // BbuParamProxy.ele_track_end (0D_NOT_character -
      .def_property(
          "ele_track_end",
          &BbuParamProxy::ele_track_end,
          &BbuParamProxy::set_ele_track_end)
      // BbuParamProxy.ix_ele_track_end (0D_NOT_integer - Default: set to last element with a wake
      .def_property(
          "ix_ele_track_end",
          &BbuParamProxy::ix_ele_track_end,
          &BbuParamProxy::set_ix_ele_track_end)
      // BbuParamProxy.regression (0D_NOT_logical - Do regression test?
      .def_property(
          "regression",
          &BbuParamProxy::regression,
          &BbuParamProxy::set_regression)
      // BbuParamProxy.normalize_z_to_rf (0D_NOT_logical - make starting z = mod(z, rf_wavelength)? Ramp parameters
      .def_property(
          "normalize_z_to_rf",
          &BbuParamProxy::normalize_z_to_rf,
          &BbuParamProxy::set_normalize_z_to_rf)
      // BbuParamProxy.ramp_on (0D_NOT_logical -
      .def_property(
          "ramp_on", &BbuParamProxy::ramp_on, &BbuParamProxy::set_ramp_on)
      // BbuParamProxy.ramp_pattern (1D_NOT_real -
      .def_property_readonly("ramp_pattern", &BbuParamProxy::ramp_pattern)
      // BbuParamProxy.ramp_n_start (0D_NOT_integer - Index of start of ramp Internal parameters
      .def_property(
          "ramp_n_start",
          &BbuParamProxy::ramp_n_start,
          &BbuParamProxy::set_ramp_n_start)
      // BbuParamProxy.n_ramp_pattern (0D_NOT_integer - Number of valid ramp_pattern
      .def_property(
          "n_ramp_pattern",
          &BbuParamProxy::n_ramp_pattern,
          &BbuParamProxy::set_n_ramp_pattern)

      .def(
          "__repr__", [](const BbuParamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuParamProxy& self) {
            return BbuParamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BbuParamProxy& self, py::dict& memo) {
            return BbuParamProxy(self);
          })

      ;

  // 1D BbuParamProxy arrays are not used in structs/routines
  // 2D BbuParamProxy arrays are not used in structs/routines
  // 3D BbuParamProxy arrays are not used in structs/routines
}

// =============================================================================
// bbu_stage_struct
void init_bbu_stage_struct(py::module& m, py::class_<BbuStageProxy>& cls) {
  cls.def(py::init<>())
      // BbuStageProxy.ix_ele_lr_wake (0D_NOT_integer - Element index of element with the wake
      .def_property(
          "ix_ele_lr_wake",
          &BbuStageProxy::ix_ele_lr_wake,
          &BbuStageProxy::set_ix_ele_lr_wake)
      // BbuStageProxy.ix_ele_stage_end (0D_NOT_integer - Element at end of stage.
      .def_property(
          "ix_ele_stage_end",
          &BbuStageProxy::ix_ele_stage_end,
          &BbuStageProxy::set_ix_ele_stage_end)
      // BbuStageProxy.ix_pass (0D_NOT_integer - Pass index when in multipass section
      .def_property(
          "ix_pass", &BbuStageProxy::ix_pass, &BbuStageProxy::set_ix_pass)
      // BbuStageProxy.ix_stage_pass1 (0D_NOT_integer - Index of corresponding stage on first pass
      .def_property(
          "ix_stage_pass1",
          &BbuStageProxy::ix_stage_pass1,
          &BbuStageProxy::set_ix_stage_pass1)
      // BbuStageProxy.ix_head_bunch (0D_NOT_integer -
      .def_property(
          "ix_head_bunch",
          &BbuStageProxy::ix_head_bunch,
          &BbuStageProxy::set_ix_head_bunch)
      // BbuStageProxy.ix_hom_max (0D_NOT_integer -
      .def_property(
          "ix_hom_max",
          &BbuStageProxy::ix_hom_max,
          &BbuStageProxy::set_ix_hom_max)
      // BbuStageProxy.hom_voltage_max (0D_NOT_real -
      .def_property(
          "hom_voltage_max",
          &BbuStageProxy::hom_voltage_max,
          &BbuStageProxy::set_hom_voltage_max)
      // BbuStageProxy.time_at_wake_ele (0D_NOT_real -
      .def_property(
          "time_at_wake_ele",
          &BbuStageProxy::time_at_wake_ele,
          &BbuStageProxy::set_time_at_wake_ele)
      // BbuStageProxy.ave_orb (1D_NOT_real -
      .def_property_readonly("ave_orb", &BbuStageProxy::ave_orb)
      // BbuStageProxy.rms_orb (1D_NOT_real -
      .def_property_readonly("rms_orb", &BbuStageProxy::rms_orb)
      // BbuStageProxy.min_orb (1D_NOT_real -
      .def_property_readonly("min_orb", &BbuStageProxy::min_orb)
      // BbuStageProxy.max_orb (1D_NOT_real -
      .def_property_readonly("max_orb", &BbuStageProxy::max_orb)
      // BbuStageProxy.n_orb (0D_NOT_integer -
      .def_property("n_orb", &BbuStageProxy::n_orb, &BbuStageProxy::set_n_orb)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BbuStageProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const BbuStageProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BbuStageProxy& self) {
            return BbuStageProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BbuStageProxy& self, py::dict& memo) {
            return BbuStageProxy(self);
          })

      ;

  bind_FTypeArrayND<BbuStageProxyArray1D>(m, "BbuStageStructArray1D");
  bind_FTypeAlloc1D<BbuStageProxyAlloc1D>(m, "BbuStageStructAlloc1D");
  // 2D BbuStageProxy arrays are not used in structs/routines
  // 3D BbuStageProxy arrays are not used in structs/routines
}

// =============================================================================
// beam_init_struct
void init_beam_init_struct(py::module& m, py::class_<BeamInitProxy>& cls) {
  cls.def(py::init<>())
      // BeamInitProxy.position_file (0D_NOT_character - File with particle positions.
      .def_property(
          "position_file",
          &BeamInitProxy::position_file,
          &BeamInitProxy::set_position_file)
      // BeamInitProxy.distribution_type (1D_NOT_character - distribution type (in x-px, y-py, and z-pz planes) 'ELLIPSE', 'KV', 'GRID', 'FILE', 'RAN_GAUSS' or '' = 'RAN_GAUSS'
      .def_property_readonly(
          "distribution_type", &BeamInitProxy::distribution_type)
      // BeamInitProxy.spin (1D_NOT_real - Spin (x, y, z)
      .def_property_readonly("spin", &BeamInitProxy::spin)
      // BeamInitProxy.ellipse (1D_NOT_type - Ellipse beam distribution
      .def_property_readonly("ellipse", &BeamInitProxy::ellipse)
      // BeamInitProxy.KV (0D_NOT_type - KV beam distribution
      .def_property("KV", &BeamInitProxy::KV, &BeamInitProxy::set_KV)
      // BeamInitProxy.grid (1D_NOT_type - Grid beam distribution
      .def_property_readonly("grid", &BeamInitProxy::grid)
      // BeamInitProxy.center_jitter (1D_NOT_real - Bunch center rms jitter
      .def_property_readonly("center_jitter", &BeamInitProxy::center_jitter)
      // BeamInitProxy.emit_jitter (1D_NOT_real - a and b bunch emittance rms jitter normalized to emittance
      .def_property_readonly("emit_jitter", &BeamInitProxy::emit_jitter)
      // BeamInitProxy.sig_z_jitter (0D_NOT_real - bunch length RMS jitter
      .def_property(
          "sig_z_jitter",
          &BeamInitProxy::sig_z_jitter,
          &BeamInitProxy::set_sig_z_jitter)
      // BeamInitProxy.sig_pz_jitter (0D_NOT_real - RMS pz spread jitter
      .def_property(
          "sig_pz_jitter",
          &BeamInitProxy::sig_pz_jitter,
          &BeamInitProxy::set_sig_pz_jitter)
      // BeamInitProxy.n_particle (0D_NOT_integer - Number of particles per bunch.
      .def_property(
          "n_particle",
          &BeamInitProxy::n_particle,
          &BeamInitProxy::set_n_particle)
      // BeamInitProxy.renorm_center (0D_NOT_logical - Renormalize centroid?
      .def_property(
          "renorm_center",
          &BeamInitProxy::renorm_center,
          &BeamInitProxy::set_renorm_center)
      // BeamInitProxy.renorm_sigma (0D_NOT_logical - Renormalize sigma?
      .def_property(
          "renorm_sigma",
          &BeamInitProxy::renorm_sigma,
          &BeamInitProxy::set_renorm_sigma)
      // BeamInitProxy.random_engine (0D_NOT_character - Or 'quasi'. Random number engine to use.
      .def_property(
          "random_engine",
          &BeamInitProxy::random_engine,
          &BeamInitProxy::set_random_engine)
      // BeamInitProxy.random_gauss_converter (0D_NOT_character - Or 'quick'. Uniform to gauss conversion method.
      .def_property(
          "random_gauss_converter",
          &BeamInitProxy::random_gauss_converter,
          &BeamInitProxy::set_random_gauss_converter)
      // BeamInitProxy.random_sigma_cutoff (0D_NOT_real - Cut-off in sigmas.
      .def_property(
          "random_sigma_cutoff",
          &BeamInitProxy::random_sigma_cutoff,
          &BeamInitProxy::set_random_sigma_cutoff)
      // BeamInitProxy.a_norm_emit (0D_NOT_real - a-mode normalized emittance (emit * beta * gamma)
      .def_property(
          "a_norm_emit",
          &BeamInitProxy::a_norm_emit,
          &BeamInitProxy::set_a_norm_emit)
      // BeamInitProxy.b_norm_emit (0D_NOT_real - b-mode normalized emittance (emit * beta * gamma)
      .def_property(
          "b_norm_emit",
          &BeamInitProxy::b_norm_emit,
          &BeamInitProxy::set_b_norm_emit)
      // BeamInitProxy.a_emit (0D_NOT_real - a-mode emittance
      .def_property(
          "a_emit", &BeamInitProxy::a_emit, &BeamInitProxy::set_a_emit)
      // BeamInitProxy.b_emit (0D_NOT_real - b-mode emittance
      .def_property(
          "b_emit", &BeamInitProxy::b_emit, &BeamInitProxy::set_b_emit)
      // BeamInitProxy.dPz_dz (0D_NOT_real - Correlation of Pz with long position.
      .def_property(
          "dPz_dz", &BeamInitProxy::dPz_dz, &BeamInitProxy::set_dPz_dz)
      // BeamInitProxy.center (1D_NOT_real - Bench phase space center offset relative to reference.
      .def_property_readonly("center", &BeamInitProxy::center)
      // BeamInitProxy.t_offset (0D_NOT_real - Time center offset
      .def_property(
          "t_offset", &BeamInitProxy::t_offset, &BeamInitProxy::set_t_offset)
      // BeamInitProxy.dt_bunch (0D_NOT_real - Time between bunches.
      .def_property(
          "dt_bunch", &BeamInitProxy::dt_bunch, &BeamInitProxy::set_dt_bunch)
      // BeamInitProxy.sig_z (0D_NOT_real - Z sigma in m.
      .def_property("sig_z", &BeamInitProxy::sig_z, &BeamInitProxy::set_sig_z)
      // BeamInitProxy.sig_pz (0D_NOT_real - pz sigma
      .def_property(
          "sig_pz", &BeamInitProxy::sig_pz, &BeamInitProxy::set_sig_pz)
      // BeamInitProxy.bunch_charge (0D_NOT_real - charge (Coul) in a bunch.
      .def_property(
          "bunch_charge",
          &BeamInitProxy::bunch_charge,
          &BeamInitProxy::set_bunch_charge)
      // BeamInitProxy.n_bunch (0D_NOT_integer - Number of bunches.
      .def_property(
          "n_bunch", &BeamInitProxy::n_bunch, &BeamInitProxy::set_n_bunch)
      // BeamInitProxy.ix_turn (0D_NOT_integer - Turn index used to adjust particles time if needed.
      .def_property(
          "ix_turn", &BeamInitProxy::ix_turn, &BeamInitProxy::set_ix_turn)
      // BeamInitProxy.species (0D_NOT_character - 'positron', etc. '' => use referece particle.
      .def_property(
          "species", &BeamInitProxy::species, &BeamInitProxy::set_species)
      // BeamInitProxy.full_6D_coupling_calc (0D_NOT_logical - Use V from 6x6 1-turn mat to match distribution? Else use 4x4 1-turn mat used.
      .def_property(
          "full_6D_coupling_calc",
          &BeamInitProxy::full_6D_coupling_calc,
          &BeamInitProxy::set_full_6D_coupling_calc)
      // BeamInitProxy.use_particle_start (0D_NOT_logical - Use lat%particle_start instead of beam_init%center, %t_offset, and %spin?
      .def_property(
          "use_particle_start",
          &BeamInitProxy::use_particle_start,
          &BeamInitProxy::set_use_particle_start)
      // BeamInitProxy.use_t_coords (0D_NOT_logical - If true, the distributions will be taken as in t-coordinates
      .def_property(
          "use_t_coords",
          &BeamInitProxy::use_t_coords,
          &BeamInitProxy::set_use_t_coords)
      // BeamInitProxy.use_z_as_t (0D_NOT_logical - Only used if  use_t_coords = .true. If true,  z describes the t distribution If false, z describes the s distribution
      .def_property(
          "use_z_as_t",
          &BeamInitProxy::use_z_as_t,
          &BeamInitProxy::set_use_z_as_t)
      // BeamInitProxy.file_name (0D_NOT_character - OLD!! DO NOT USE!!
      .def_property(
          "file_name", &BeamInitProxy::file_name, &BeamInitProxy::set_file_name)

      .def(
          "__repr__", [](const BeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BeamInitProxy& self) {
            return BeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BeamInitProxy& self, py::dict& memo) {
            return BeamInitProxy(self);
          })

      ;

  // 1D BeamInitProxy arrays are not used in structs/routines
  // 2D BeamInitProxy arrays are not used in structs/routines
  // 3D BeamInitProxy arrays are not used in structs/routines
}

// =============================================================================
// beam_struct
void init_beam_struct(py::module& m, py::class_<BeamProxy>& cls) {
  cls.def(py::init<>())
      // BeamProxy.bunch (1D_ALLOC_type -
      .def_property_readonly("bunch", &BeamProxy::bunch)

      .def("__repr__", [](const BeamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BeamProxy& self) {
            return BeamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BeamProxy& self, py::dict& memo) { return BeamProxy(self); })

      ;

  // 1D BeamProxy arrays are not used in structs/routines
  // 2D BeamProxy arrays are not used in structs/routines
  // 3D BeamProxy arrays are not used in structs/routines
}

// =============================================================================
// bmad_common_struct
void init_bmad_common_struct(py::module& m, py::class_<BmadCommonProxy>& cls) {
  cls.def(py::init<>())
      // BmadCommonProxy.max_aperture_limit (0D_NOT_real - Max Aperture.
      .def_property(
          "max_aperture_limit",
          &BmadCommonProxy::max_aperture_limit,
          &BmadCommonProxy::set_max_aperture_limit)
      // BmadCommonProxy.d_orb (1D_NOT_real - Orbit deltas for the mat6 via tracking calc.
      .def_property_readonly("d_orb", &BmadCommonProxy::d_orb)
      // BmadCommonProxy.default_ds_step (0D_NOT_real - Default integration step for eles without an explicit step calc.
      .def_property(
          "default_ds_step",
          &BmadCommonProxy::default_ds_step,
          &BmadCommonProxy::set_default_ds_step)
      // BmadCommonProxy.significant_length (0D_NOT_real - meter
      .def_property(
          "significant_length",
          &BmadCommonProxy::significant_length,
          &BmadCommonProxy::set_significant_length)
      // BmadCommonProxy.rel_tol_tracking (0D_NOT_real - Closed orbit relative tolerance.
      .def_property(
          "rel_tol_tracking",
          &BmadCommonProxy::rel_tol_tracking,
          &BmadCommonProxy::set_rel_tol_tracking)
      // BmadCommonProxy.abs_tol_tracking (0D_NOT_real - Closed orbit absolute tolerance.
      .def_property(
          "abs_tol_tracking",
          &BmadCommonProxy::abs_tol_tracking,
          &BmadCommonProxy::set_abs_tol_tracking)
      // BmadCommonProxy.rel_tol_adaptive_tracking (0D_NOT_real - Runge-Kutta tracking relative tolerance.
      .def_property(
          "rel_tol_adaptive_tracking",
          &BmadCommonProxy::rel_tol_adaptive_tracking,
          &BmadCommonProxy::set_rel_tol_adaptive_tracking)
      // BmadCommonProxy.abs_tol_adaptive_tracking (0D_NOT_real - Runge-Kutta tracking absolute tolerance.
      .def_property(
          "abs_tol_adaptive_tracking",
          &BmadCommonProxy::abs_tol_adaptive_tracking,
          &BmadCommonProxy::set_abs_tol_adaptive_tracking)
      // BmadCommonProxy.init_ds_adaptive_tracking (0D_NOT_real - Initial step size
      .def_property(
          "init_ds_adaptive_tracking",
          &BmadCommonProxy::init_ds_adaptive_tracking,
          &BmadCommonProxy::set_init_ds_adaptive_tracking)
      // BmadCommonProxy.min_ds_adaptive_tracking (0D_NOT_real - Min step size to take.
      .def_property(
          "min_ds_adaptive_tracking",
          &BmadCommonProxy::min_ds_adaptive_tracking,
          &BmadCommonProxy::set_min_ds_adaptive_tracking)
      // BmadCommonProxy.fatal_ds_adaptive_tracking (0D_NOT_real - If actual step size is below this particle is lost.
      .def_property(
          "fatal_ds_adaptive_tracking",
          &BmadCommonProxy::fatal_ds_adaptive_tracking,
          &BmadCommonProxy::set_fatal_ds_adaptive_tracking)
      // BmadCommonProxy.autoscale_amp_abs_tol (0D_NOT_real - Autoscale absolute amplitude tolerance (eV).
      .def_property(
          "autoscale_amp_abs_tol",
          &BmadCommonProxy::autoscale_amp_abs_tol,
          &BmadCommonProxy::set_autoscale_amp_abs_tol)
      // BmadCommonProxy.autoscale_amp_rel_tol (0D_NOT_real - Autoscale relative amplitude tolerance
      .def_property(
          "autoscale_amp_rel_tol",
          &BmadCommonProxy::autoscale_amp_rel_tol,
          &BmadCommonProxy::set_autoscale_amp_rel_tol)
      // BmadCommonProxy.autoscale_phase_tol (0D_NOT_real - Autoscale phase tolerance.
      .def_property(
          "autoscale_phase_tol",
          &BmadCommonProxy::autoscale_phase_tol,
          &BmadCommonProxy::set_autoscale_phase_tol)
      // BmadCommonProxy.electric_dipole_moment (0D_NOT_real - Particle's EDM. Call set_ptc to transfer value to PTC.
      .def_property(
          "electric_dipole_moment",
          &BmadCommonProxy::electric_dipole_moment,
          &BmadCommonProxy::set_electric_dipole_moment)
      // BmadCommonProxy.synch_rad_scale (0D_NOT_real - Synch radiation kick scale. 1 => normal, 0 => no kicks.
      .def_property(
          "synch_rad_scale",
          &BmadCommonProxy::synch_rad_scale,
          &BmadCommonProxy::set_synch_rad_scale)
      // BmadCommonProxy.sad_eps_scale (0D_NOT_real - Used in sad_mult step length calc.
      .def_property(
          "sad_eps_scale",
          &BmadCommonProxy::sad_eps_scale,
          &BmadCommonProxy::set_sad_eps_scale)
      // BmadCommonProxy.sad_amp_max (0D_NOT_real - Used in sad_mult step length calc.
      .def_property(
          "sad_amp_max",
          &BmadCommonProxy::sad_amp_max,
          &BmadCommonProxy::set_sad_amp_max)
      // BmadCommonProxy.sad_n_div_max (0D_NOT_integer - Used in sad_mult step length calc.
      .def_property(
          "sad_n_div_max",
          &BmadCommonProxy::sad_n_div_max,
          &BmadCommonProxy::set_sad_n_div_max)
      // BmadCommonProxy.taylor_order (0D_NOT_integer - Taylor order to use. 0 -> default = ptc_private%taylor_order_saved.
      .def_property(
          "taylor_order",
          &BmadCommonProxy::taylor_order,
          &BmadCommonProxy::set_taylor_order)
      // BmadCommonProxy.runge_kutta_order (0D_NOT_integer - Runge Kutta order.
      .def_property(
          "runge_kutta_order",
          &BmadCommonProxy::runge_kutta_order,
          &BmadCommonProxy::set_runge_kutta_order)
      // BmadCommonProxy.default_integ_order (0D_NOT_integer - PTC integration order.
      .def_property(
          "default_integ_order",
          &BmadCommonProxy::default_integ_order,
          &BmadCommonProxy::set_default_integ_order)
      // BmadCommonProxy.max_num_runge_kutta_step (0D_NOT_integer - Maximum number of RK steps before particle is considered lost.
      .def_property(
          "max_num_runge_kutta_step",
          &BmadCommonProxy::max_num_runge_kutta_step,
          &BmadCommonProxy::set_max_num_runge_kutta_step)
      // BmadCommonProxy.rf_phase_below_transition_ref (0D_NOT_logical - Autoscale uses below transition stable point for RFCavities?
      .def_property(
          "rf_phase_below_transition_ref",
          &BmadCommonProxy::rf_phase_below_transition_ref,
          &BmadCommonProxy::set_rf_phase_below_transition_ref)
      // BmadCommonProxy.sr_wakes_on (0D_NOT_logical - Short range wakefields?
      .def_property(
          "sr_wakes_on",
          &BmadCommonProxy::sr_wakes_on,
          &BmadCommonProxy::set_sr_wakes_on)
      // BmadCommonProxy.lr_wakes_on (0D_NOT_logical - Long range wakefields
      .def_property(
          "lr_wakes_on",
          &BmadCommonProxy::lr_wakes_on,
          &BmadCommonProxy::set_lr_wakes_on)
      // BmadCommonProxy.auto_bookkeeper (0D_NOT_logical - Deprecated and no longer used.
      .def_property(
          "auto_bookkeeper",
          &BmadCommonProxy::auto_bookkeeper,
          &BmadCommonProxy::set_auto_bookkeeper)
      // BmadCommonProxy.high_energy_space_charge_on (0D_NOT_logical - High energy space charge effect switch.
      .def_property(
          "high_energy_space_charge_on",
          &BmadCommonProxy::high_energy_space_charge_on,
          &BmadCommonProxy::set_high_energy_space_charge_on)
      // BmadCommonProxy.csr_and_space_charge_on (0D_NOT_logical - Space charge switch.
      .def_property(
          "csr_and_space_charge_on",
          &BmadCommonProxy::csr_and_space_charge_on,
          &BmadCommonProxy::set_csr_and_space_charge_on)
      // BmadCommonProxy.spin_tracking_on (0D_NOT_logical - spin tracking?
      .def_property(
          "spin_tracking_on",
          &BmadCommonProxy::spin_tracking_on,
          &BmadCommonProxy::set_spin_tracking_on)
      // BmadCommonProxy.spin_sokolov_ternov_flipping_on (0D_NOT_logical - Spin flipping during synchrotron radiation emission?
      .def_property(
          "spin_sokolov_ternov_flipping_on",
          &BmadCommonProxy::spin_sokolov_ternov_flipping_on,
          &BmadCommonProxy::set_spin_sokolov_ternov_flipping_on)
      // BmadCommonProxy.radiation_damping_on (0D_NOT_logical - Radiation damping toggle.
      .def_property(
          "radiation_damping_on",
          &BmadCommonProxy::radiation_damping_on,
          &BmadCommonProxy::set_radiation_damping_on)
      // BmadCommonProxy.radiation_zero_average (0D_NOT_logical - Shift damping to be zero on the zero orbit to get rid of sawtooth?
      .def_property(
          "radiation_zero_average",
          &BmadCommonProxy::radiation_zero_average,
          &BmadCommonProxy::set_radiation_zero_average)
      // BmadCommonProxy.radiation_fluctuations_on (0D_NOT_logical - Radiation fluctuations toggle.
      .def_property(
          "radiation_fluctuations_on",
          &BmadCommonProxy::radiation_fluctuations_on,
          &BmadCommonProxy::set_radiation_fluctuations_on)
      // BmadCommonProxy.conserve_taylor_maps (0D_NOT_logical - Enable bookkeeper to set ele%taylor_map_includes_offsets = F?
      .def_property(
          "conserve_taylor_maps",
          &BmadCommonProxy::conserve_taylor_maps,
          &BmadCommonProxy::set_conserve_taylor_maps)
      // BmadCommonProxy.absolute_time_tracking (0D_NOT_logical - Absolute or relative time tracking?
      .def_property(
          "absolute_time_tracking",
          &BmadCommonProxy::absolute_time_tracking,
          &BmadCommonProxy::set_absolute_time_tracking)
      // BmadCommonProxy.absolute_time_ref_shift (0D_NOT_logical - Apply reference time shift when using absolute time tracking?
      .def_property(
          "absolute_time_ref_shift",
          &BmadCommonProxy::absolute_time_ref_shift,
          &BmadCommonProxy::set_absolute_time_ref_shift)
      // BmadCommonProxy.convert_to_kinetic_momentum (0D_NOT_logical - Cancel kicks due to finite vector potential when doing symplectic tracking? Set to True to test symp_lie_bmad against runge_kutta.
      .def_property(
          "convert_to_kinetic_momentum",
          &BmadCommonProxy::convert_to_kinetic_momentum,
          &BmadCommonProxy::set_convert_to_kinetic_momentum)
      // BmadCommonProxy.normalize_twiss (0D_NOT_logical - Normalize matrix when computing Twiss for off-energy ref?
      .def_property(
          "normalize_twiss",
          &BmadCommonProxy::normalize_twiss,
          &BmadCommonProxy::set_normalize_twiss)
      // BmadCommonProxy.aperture_limit_on (0D_NOT_logical - Use apertures in tracking?
      .def_property(
          "aperture_limit_on",
          &BmadCommonProxy::aperture_limit_on,
          &BmadCommonProxy::set_aperture_limit_on)
      // BmadCommonProxy.spin_n0_direction_user_set (0D_NOT_logical - User sets direction of n0 for closed geometry branches?
      .def_property(
          "spin_n0_direction_user_set",
          &BmadCommonProxy::spin_n0_direction_user_set,
          &BmadCommonProxy::set_spin_n0_direction_user_set)
      // BmadCommonProxy.debug (0D_NOT_logical - Used for code debugging.
      .def_property(
          "debug", &BmadCommonProxy::debug, &BmadCommonProxy::set_debug)

      .def(
          "__repr__",
          [](const BmadCommonProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BmadCommonProxy& self) {
            return BmadCommonProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BmadCommonProxy& self, py::dict& memo) {
            return BmadCommonProxy(self);
          })

      ;

  // 1D BmadCommonProxy arrays are not used in structs/routines
  // 2D BmadCommonProxy arrays are not used in structs/routines
  // 3D BmadCommonProxy arrays are not used in structs/routines
}

// =============================================================================
// bmad_normal_form_struct
void init_bmad_normal_form_struct(
    py::module& m,
    py::class_<BmadNormalFormProxy>& cls) {
  cls.def(py::init<>())
      // BmadNormalFormProxy.ele_origin (0D_PTR_type - Element at which the on-turn map was created.
      .def_property(
          "ele_origin",
          &BmadNormalFormProxy::ele_origin,
          &BmadNormalFormProxy::set_ele_origin)
      // BmadNormalFormProxy.M (1D_NOT_type - One-turn taylor map: M = A o N o A_inv, N = exp(:h:)
      .def_property_readonly("M", &BmadNormalFormProxy::M)
      // BmadNormalFormProxy.A (1D_NOT_type - Map from Floquet -> Lab coordinates
      .def_property_readonly("A", &BmadNormalFormProxy::A)
      // BmadNormalFormProxy.A_inv (1D_NOT_type - Map from Lab -> Floquet coordinates
      .def_property_readonly("A_inv", &BmadNormalFormProxy::A_inv)
      // BmadNormalFormProxy.dhdj (1D_NOT_type - Nonlinear tune function operating on Floquet coordinates
      .def_property_readonly("dhdj", &BmadNormalFormProxy::dhdj)
      // BmadNormalFormProxy.F (1D_NOT_type - Vector field factorization in phasor basis:
      .def_property_readonly("F", &BmadNormalFormProxy::F)
      // BmadNormalFormProxy.L (1D_NOT_type - L component
      .def_property_readonly("L", &BmadNormalFormProxy::L)
      // BmadNormalFormProxy.h (1D_ALLOC_type -
      .def_property_readonly("h", &BmadNormalFormProxy::h)

      .def(
          "__repr__",
          [](const BmadNormalFormProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BmadNormalFormProxy& self) {
            return BmadNormalFormProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BmadNormalFormProxy& self, py::dict& memo) {
            return BmadNormalFormProxy(self);
          })

      ;

  // 1D BmadNormalFormProxy arrays are not used in structs/routines
  // 2D BmadNormalFormProxy arrays are not used in structs/routines
  // 3D BmadNormalFormProxy arrays are not used in structs/routines
}

// =============================================================================
// bookkeeping_state_struct
void init_bookkeeping_state_struct(
    py::module& m,
    py::class_<BookkeepingStateProxy>& cls) {
  cls.def(py::init<>())
      // BookkeepingStateProxy.attributes (0D_NOT_integer - Element dependent attributes: super_ok$, ok$ or stale$
      .def_property(
          "attributes",
          &BookkeepingStateProxy::attributes,
          &BookkeepingStateProxy::set_attributes)
      // BookkeepingStateProxy.control (0D_NOT_integer - Lord/slave bookkeeping status: super_ok$, ok$ or stale$
      .def_property(
          "control",
          &BookkeepingStateProxy::control,
          &BookkeepingStateProxy::set_control)
      // BookkeepingStateProxy.floor_position (0D_NOT_integer - Global (floor) geometry: super_ok$, ok$ or stale$
      .def_property(
          "floor_position",
          &BookkeepingStateProxy::floor_position,
          &BookkeepingStateProxy::set_floor_position)
      // BookkeepingStateProxy.s_position (0D_NOT_integer - Longitudinal position & element length: super_ok$, ok$ or stale$
      .def_property(
          "s_position",
          &BookkeepingStateProxy::s_position,
          &BookkeepingStateProxy::set_s_position)
      // BookkeepingStateProxy.ref_energy (0D_NOT_integer - Reference energy and ref time: super_ok$, ok$ or stale$
      .def_property(
          "ref_energy",
          &BookkeepingStateProxy::ref_energy,
          &BookkeepingStateProxy::set_ref_energy)
      // BookkeepingStateProxy.mat6 (0D_NOT_integer - Linear transfer map status: super_ok$, ok$ or stale$
      .def_property(
          "mat6",
          &BookkeepingStateProxy::mat6,
          &BookkeepingStateProxy::set_mat6)
      // BookkeepingStateProxy.rad_int (0D_NOT_integer - Radiation integrals cache status
      .def_property(
          "rad_int",
          &BookkeepingStateProxy::rad_int,
          &BookkeepingStateProxy::set_rad_int)
      // BookkeepingStateProxy.ptc (0D_NOT_integer - Associated PTC fibre (or layout) status.
      .def_property(
          "ptc", &BookkeepingStateProxy::ptc, &BookkeepingStateProxy::set_ptc)
      // BookkeepingStateProxy.has_misalign (0D_NOT_logical - Used to avoid unnecessary calls to offset_particle.
      .def_property(
          "has_misalign",
          &BookkeepingStateProxy::has_misalign,
          &BookkeepingStateProxy::set_has_misalign)

      .def(
          "__repr__",
          [](const BookkeepingStateProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BookkeepingStateProxy& self) {
            return BookkeepingStateProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BookkeepingStateProxy& self, py::dict& memo) {
            return BookkeepingStateProxy(self);
          })

      ;

  // 1D BookkeepingStateProxy arrays are not used in structs/routines
  // 2D BookkeepingStateProxy arrays are not used in structs/routines
  // 3D BookkeepingStateProxy arrays are not used in structs/routines
}

// =============================================================================
// bpm_phase_coupling_struct
void init_bpm_phase_coupling_struct(
    py::module& m,
    py::class_<BpmPhaseCouplingProxy>& cls) {
  cls.def(py::init<>())
      // BpmPhaseCouplingProxy.K_22a (0D_NOT_real - In-phase y/x for a-mode oscillations.
      .def_property(
          "K_22a",
          &BpmPhaseCouplingProxy::K_22a,
          &BpmPhaseCouplingProxy::set_K_22a)
      // BpmPhaseCouplingProxy.K_12a (0D_NOT_real - Out-of-phase y/x for a-mode oscillations.
      .def_property(
          "K_12a",
          &BpmPhaseCouplingProxy::K_12a,
          &BpmPhaseCouplingProxy::set_K_12a)
      // BpmPhaseCouplingProxy.K_11b (0D_NOT_real - In-phase x/y for b-mode oscillations.
      .def_property(
          "K_11b",
          &BpmPhaseCouplingProxy::K_11b,
          &BpmPhaseCouplingProxy::set_K_11b)
      // BpmPhaseCouplingProxy.K_12b (0D_NOT_real - Out-of-phase x/y for b-mode oscillations.
      .def_property(
          "K_12b",
          &BpmPhaseCouplingProxy::K_12b,
          &BpmPhaseCouplingProxy::set_K_12b)
      // BpmPhaseCouplingProxy.Cbar22_a (0D_NOT_real - Cbar22 as calculated from K_22a.
      .def_property(
          "Cbar22_a",
          &BpmPhaseCouplingProxy::Cbar22_a,
          &BpmPhaseCouplingProxy::set_Cbar22_a)
      // BpmPhaseCouplingProxy.Cbar12_a (0D_NOT_real - Cbar12 as calculated from K_12a.
      .def_property(
          "Cbar12_a",
          &BpmPhaseCouplingProxy::Cbar12_a,
          &BpmPhaseCouplingProxy::set_Cbar12_a)
      // BpmPhaseCouplingProxy.Cbar11_b (0D_NOT_real - Cbar11 as calculated from K_11b.
      .def_property(
          "Cbar11_b",
          &BpmPhaseCouplingProxy::Cbar11_b,
          &BpmPhaseCouplingProxy::set_Cbar11_b)
      // BpmPhaseCouplingProxy.Cbar12_b (0D_NOT_real - Cbar12 as calculated from K_12b.
      .def_property(
          "Cbar12_b",
          &BpmPhaseCouplingProxy::Cbar12_b,
          &BpmPhaseCouplingProxy::set_Cbar12_b)
      // BpmPhaseCouplingProxy.phi_a (0D_NOT_real - a-mode betatron phase.
      .def_property(
          "phi_a",
          &BpmPhaseCouplingProxy::phi_a,
          &BpmPhaseCouplingProxy::set_phi_a)
      // BpmPhaseCouplingProxy.phi_b (0D_NOT_real - b-mode betatron phase.
      .def_property(
          "phi_b",
          &BpmPhaseCouplingProxy::phi_b,
          &BpmPhaseCouplingProxy::set_phi_b)

      .def(
          "__repr__",
          [](const BpmPhaseCouplingProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BpmPhaseCouplingProxy& self) {
            return BpmPhaseCouplingProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BpmPhaseCouplingProxy& self, py::dict& memo) {
            return BpmPhaseCouplingProxy(self);
          })

      ;

  // 1D BpmPhaseCouplingProxy arrays are not used in structs/routines
  // 2D BpmPhaseCouplingProxy arrays are not used in structs/routines
  // 3D BpmPhaseCouplingProxy arrays are not used in structs/routines
}

// =============================================================================
// branch_struct
void init_branch_struct(py::module& m, py::class_<BranchProxy>& cls) {
  cls.def(py::init<>())
      // BranchProxy.name (0D_NOT_character - Name of line that defines the branch.
      .def_property("name", &BranchProxy::name, &BranchProxy::set_name)
      // BranchProxy.ix_branch (0D_NOT_integer - Index of this branch. 0 => Main branch
      .def_property(
          "ix_branch", &BranchProxy::ix_branch, &BranchProxy::set_ix_branch)
      // BranchProxy.ix_from_branch (0D_NOT_integer - -1 => No creating fork element to this branch.
      .def_property(
          "ix_from_branch",
          &BranchProxy::ix_from_branch,
          &BranchProxy::set_ix_from_branch)
      // BranchProxy.ix_from_ele (0D_NOT_integer - Index of creating fork element which forks to this branch.
      .def_property(
          "ix_from_ele",
          &BranchProxy::ix_from_ele,
          &BranchProxy::set_ix_from_ele)
      // BranchProxy.ix_to_ele (0D_NOT_integer - Index of element in this branch that creating fork element forks to.
      .def_property(
          "ix_to_ele", &BranchProxy::ix_to_ele, &BranchProxy::set_ix_to_ele)
      // BranchProxy.ix_fixer (0D_NOT_integer - Index of active fixer or beginning_ele element.
      .def_property(
          "ix_fixer", &BranchProxy::ix_fixer, &BranchProxy::set_ix_fixer)
      // BranchProxy.n_ele_track (0D_NOT_integer -
      .def_property(
          "n_ele_track",
          &BranchProxy::n_ele_track,
          &BranchProxy::set_n_ele_track)
      // BranchProxy.n_ele_max (0D_NOT_integer -
      .def_property(
          "n_ele_max", &BranchProxy::n_ele_max, &BranchProxy::set_n_ele_max)
      // BranchProxy.lat (0D_PTR_type -
      .def_property("lat", &BranchProxy::lat, &BranchProxy::set_lat)
      // BranchProxy.a (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("a", &BranchProxy::a, &BranchProxy::set_a)
      // BranchProxy.b (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("b", &BranchProxy::b, &BranchProxy::set_b)
      // BranchProxy.z (0D_NOT_type - Note: Tunes are the fractional part.
      .def_property("z", &BranchProxy::z, &BranchProxy::set_z)
      // BranchProxy.ele (1D_PTR_type -
      .def_property_readonly("ele", &BranchProxy::ele)
      // BranchProxy.param (0D_NOT_type -
      .def_property("param", &BranchProxy::param, &BranchProxy::set_param)
      // BranchProxy.particle_start (0D_NOT_type -
      .def_property(
          "particle_start",
          &BranchProxy::particle_start,
          &BranchProxy::set_particle_start)
      // BranchProxy.wall3d (1D_PTR_type -
      .def_property_readonly("wall3d", &BranchProxy::wall3d)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BranchProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const BranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BranchProxy& self) {
            return BranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BranchProxy& self, py::dict& memo) {
            return BranchProxy(self);
          })

      ;

  bind_FTypeArrayND<BranchProxyArray1D>(m, "BranchStructArray1D");
  bind_FTypeAlloc1D<BranchProxyAlloc1D>(m, "BranchStructAlloc1D");
  // 2D BranchProxy arrays are not used in structs/routines
  // 3D BranchProxy arrays are not used in structs/routines
}

// =============================================================================
// bunch_params_struct
void init_bunch_params_struct(
    py::module& m,
    py::class_<BunchParamsProxy>& cls) {
  cls.def(py::init<>())
      // BunchParamsProxy.centroid (0D_NOT_type - Lab frame
      .def_property(
          "centroid",
          &BunchParamsProxy::centroid,
          &BunchParamsProxy::set_centroid)
      // BunchParamsProxy.x (0D_NOT_type - Projected Twiss parameters
      .def_property("x", &BunchParamsProxy::x, &BunchParamsProxy::set_x)
      // BunchParamsProxy.y (0D_NOT_type - Projected Twiss parameters
      .def_property("y", &BunchParamsProxy::y, &BunchParamsProxy::set_y)
      // BunchParamsProxy.z (0D_NOT_type - Projected Twiss parameters
      .def_property("z", &BunchParamsProxy::z, &BunchParamsProxy::set_z)
      // BunchParamsProxy.a (0D_NOT_type - Normal mode twiss parameters
      .def_property("a", &BunchParamsProxy::a, &BunchParamsProxy::set_a)
      // BunchParamsProxy.b (0D_NOT_type - Normal mode twiss parameters
      .def_property("b", &BunchParamsProxy::b, &BunchParamsProxy::set_b)
      // BunchParamsProxy.c (0D_NOT_type - Normal mode twiss parameters
      .def_property("c", &BunchParamsProxy::c, &BunchParamsProxy::set_c)
      // BunchParamsProxy.sigma (2D_NOT_real - beam size matrix
      .def_property_readonly("sigma", &BunchParamsProxy::sigma)
      // BunchParamsProxy.rel_max (1D_NOT_real - Max orbit relative to centroid. 7 -> time.
      .def_property_readonly("rel_max", &BunchParamsProxy::rel_max)
      // BunchParamsProxy.rel_min (1D_NOT_real - Min orbit relative to_centroid. 7 -> time.
      .def_property_readonly("rel_min", &BunchParamsProxy::rel_min)
      // BunchParamsProxy.s (0D_NOT_real - Longitudinal position.
      .def_property("s", &BunchParamsProxy::s, &BunchParamsProxy::set_s)
      // BunchParamsProxy.t (0D_NOT_real - Time.
      .def_property("t", &BunchParamsProxy::t, &BunchParamsProxy::set_t)
      // BunchParamsProxy.sigma_t (0D_NOT_real - RMS of time spread.
      .def_property(
          "sigma_t", &BunchParamsProxy::sigma_t, &BunchParamsProxy::set_sigma_t)
      // BunchParamsProxy.charge_live (0D_NOT_real - Charge of all non-lost particle
      .def_property(
          "charge_live",
          &BunchParamsProxy::charge_live,
          &BunchParamsProxy::set_charge_live)
      // BunchParamsProxy.charge_tot (0D_NOT_real - Charge of all particles.
      .def_property(
          "charge_tot",
          &BunchParamsProxy::charge_tot,
          &BunchParamsProxy::set_charge_tot)
      // BunchParamsProxy.n_particle_tot (0D_NOT_integer - Total number of particles
      .def_property(
          "n_particle_tot",
          &BunchParamsProxy::n_particle_tot,
          &BunchParamsProxy::set_n_particle_tot)
      // BunchParamsProxy.n_particle_live (0D_NOT_integer - Number of non-lost particles
      .def_property(
          "n_particle_live",
          &BunchParamsProxy::n_particle_live,
          &BunchParamsProxy::set_n_particle_live)
      // BunchParamsProxy.n_particle_lost_in_ele (0D_NOT_integer - Number lost in element (not calculated by Bmad)
      .def_property(
          "n_particle_lost_in_ele",
          &BunchParamsProxy::n_particle_lost_in_ele,
          &BunchParamsProxy::set_n_particle_lost_in_ele)
      // BunchParamsProxy.n_good_steps (0D_NOT_integer - Number of good steps (set when tracking with space charge)
      .def_property(
          "n_good_steps",
          &BunchParamsProxy::n_good_steps,
          &BunchParamsProxy::set_n_good_steps)
      // BunchParamsProxy.n_bad_steps (0D_NOT_integer - Number of bad steps (set when tracking with space charge)
      .def_property(
          "n_bad_steps",
          &BunchParamsProxy::n_bad_steps,
          &BunchParamsProxy::set_n_bad_steps)
      // BunchParamsProxy.ix_ele (0D_NOT_integer - Lattice element where params evaluated at.
      .def_property(
          "ix_ele", &BunchParamsProxy::ix_ele, &BunchParamsProxy::set_ix_ele)
      // BunchParamsProxy.location (0D_NOT_integer - Location in element: upstream_end$, inside$, or downstream_end$
      .def_property(
          "location",
          &BunchParamsProxy::location,
          &BunchParamsProxy::set_location)
      // BunchParamsProxy.twiss_valid (0D_NOT_logical - Is the data here valid? Note: IF there is no energy variation (RF off) twiss_valid may be true but in this case the z-twiss will not be valid.
      .def_property(
          "twiss_valid",
          &BunchParamsProxy::twiss_valid,
          &BunchParamsProxy::set_twiss_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return BunchParamsProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const BunchParamsProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchParamsProxy& self) {
            return BunchParamsProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BunchParamsProxy& self, py::dict& memo) {
            return BunchParamsProxy(self);
          })

      ;

  bind_FTypeArrayND<BunchParamsProxyArray1D>(m, "BunchParamsStructArray1D");
  bind_FTypeAlloc1D<BunchParamsProxyAlloc1D>(m, "BunchParamsStructAlloc1D");
  // 2D BunchParamsProxy arrays are not used in structs/routines
  // 3D BunchParamsProxy arrays are not used in structs/routines
}

// =============================================================================
// bunch_struct
void init_bunch_struct(py::module& m, py::class_<BunchProxy>& cls) {
  cls.def(py::init<>())
      // BunchProxy.particle (1D_ALLOC_type -
      .def_property_readonly("particle", &BunchProxy::particle)
      // BunchProxy.ix_z (1D_ALLOC_integer - bunch%ix_z(1) is index of head particle, etc.
      .def_property_readonly("ix_z", &BunchProxy::ix_z)
      // BunchProxy.charge_tot (0D_NOT_real - Total charge in a bunch (Coul).
      .def_property(
          "charge_tot", &BunchProxy::charge_tot, &BunchProxy::set_charge_tot)
      // BunchProxy.charge_live (0D_NOT_real - Charge of live particles (Coul).
      .def_property(
          "charge_live", &BunchProxy::charge_live, &BunchProxy::set_charge_live)
      // BunchProxy.z_center (0D_NOT_real - Longitudinal center of bunch at creation time. Note: Generally, z_center of bunch #1 is 0 and z_center of the other bunches is negative.
      .def_property(
          "z_center", &BunchProxy::z_center, &BunchProxy::set_z_center)
      // BunchProxy.t_center (0D_NOT_real - Center of bunch at creation time relative to head bunch.
      .def_property(
          "t_center", &BunchProxy::t_center, &BunchProxy::set_t_center)
      // BunchProxy.t0 (0D_NOT_real - Used by track1_bunch_space_charge for tracking so particles have constant t.
      .def_property("t0", &BunchProxy::t0, &BunchProxy::set_t0)
      // BunchProxy.drift_between_t_and_s (0D_NOT_logical - Drift (ignore any fields) instead of tracking to speed up the calculation? This can only be done under certain circumstances.
      .def_property(
          "drift_between_t_and_s",
          &BunchProxy::drift_between_t_and_s,
          &BunchProxy::set_drift_between_t_and_s)
      // BunchProxy.ix_ele (0D_NOT_integer - Nominal element bunch is at. But, EG, dead particles can be someplace else.
      .def_property("ix_ele", &BunchProxy::ix_ele, &BunchProxy::set_ix_ele)
      // BunchProxy.ix_bunch (0D_NOT_integer - Bunch index. Head bunch = 1, etc.
      .def_property(
          "ix_bunch", &BunchProxy::ix_bunch, &BunchProxy::set_ix_bunch)
      // BunchProxy.ix_turn (0D_NOT_integer - Turn index for long term tracking. ix_turn = 0 before end of first turn, etc.
      .def_property("ix_turn", &BunchProxy::ix_turn, &BunchProxy::set_ix_turn)
      // BunchProxy.n_live (0D_NOT_integer -
      .def_property("n_live", &BunchProxy::n_live, &BunchProxy::set_n_live)
      // BunchProxy.n_good (0D_NOT_integer - Number of accepted steps when using adaptive step size control.
      .def_property("n_good", &BunchProxy::n_good, &BunchProxy::set_n_good)
      // BunchProxy.n_bad (0D_NOT_integer - Number of rejected steps when using adaptive step size control.
      .def_property("n_bad", &BunchProxy::n_bad, &BunchProxy::set_n_bad)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BunchProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const BunchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchProxy& self) {
            return BunchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BunchProxy& self, py::dict& memo) {
            return BunchProxy(self);
          })

      ;

  bind_FTypeArrayND<BunchProxyArray1D>(m, "BunchStructArray1D");
  bind_FTypeAlloc1D<BunchProxyAlloc1D>(m, "BunchStructAlloc1D");
  // 2D BunchProxy arrays are not used in structs/routines
  // 3D BunchProxy arrays are not used in structs/routines
}

// =============================================================================
// bunch_track_struct
void init_bunch_track_struct(py::module& m, py::class_<BunchTrackProxy>& cls) {
  cls.def(py::init<>())
      // BunchTrackProxy.pt (1D_ALLOC_type - Array indexed from 0
      .def_property_readonly("pt", &BunchTrackProxy::pt)
      // BunchTrackProxy.ds_save (0D_NOT_real - Min distance between points.
      .def_property(
          "ds_save", &BunchTrackProxy::ds_save, &BunchTrackProxy::set_ds_save)
      // BunchTrackProxy.n_pt (0D_NOT_integer - Track upper bound
      .def_property("n_pt", &BunchTrackProxy::n_pt, &BunchTrackProxy::set_n_pt)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return BunchTrackProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const BunchTrackProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BunchTrackProxy& self) {
            return BunchTrackProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BunchTrackProxy& self, py::dict& memo) {
            return BunchTrackProxy(self);
          })

      ;

  bind_FTypeArrayND<BunchTrackProxyArray1D>(m, "BunchTrackStructArray1D");
  bind_FTypeAlloc1D<BunchTrackProxyAlloc1D>(m, "BunchTrackStructAlloc1D");
  // 2D BunchTrackProxy arrays are not used in structs/routines
  // 3D BunchTrackProxy arrays are not used in structs/routines
}

// =============================================================================
// bicubic_cmplx_coef_struct
void init_bicubic_cmplx_coef_struct(
    py::module& m,
    py::class_<BicubicCmplxCoefProxy>& cls) {
  cls.def(py::init<>())
      // BicubicCmplxCoefProxy.coef (2D_NOT_complex - Coefs
      .def_property_readonly("coef", &BicubicCmplxCoefProxy::coef)
      // BicubicCmplxCoefProxy.i_box (1D_NOT_integer - index at lower box corner.
      .def_property_readonly("i_box", &BicubicCmplxCoefProxy::i_box)

      .def(
          "__repr__",
          [](const BicubicCmplxCoefProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BicubicCmplxCoefProxy& self) {
            return BicubicCmplxCoefProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const BicubicCmplxCoefProxy& self, py::dict& memo) {
            return BicubicCmplxCoefProxy(self);
          })

      ;

  // 1D BicubicCmplxCoefProxy arrays are not used in structs/routines
  // 2D BicubicCmplxCoefProxy arrays are not used in structs/routines
  bind_FTypeArrayND<BicubicCmplxCoefProxyArray3D>(
      m, "BicubicCmplxCoefStructArray3D");
}