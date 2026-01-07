#include <pybind11/pybind11.h>
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

namespace Pybmad {

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
// ac_kicker_freq_struct
void init_ac_kicker_freq_struct(
    py::module& m,
    py::class_<AcKickerFreqProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerFreqProxy.f (0D_NOT_real -
      .def_property("f", &AcKickerFreqProxy::f, &AcKickerFreqProxy::set_f)
      // AcKickerFreqProxy.amp (0D_NOT_real -
      .def_property("amp", &AcKickerFreqProxy::amp, &AcKickerFreqProxy::set_amp)
      // AcKickerFreqProxy.phi (0D_NOT_real -
      .def_property("phi", &AcKickerFreqProxy::phi, &AcKickerFreqProxy::set_phi)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AcKickerFreqProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AcKickerFreqProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerFreqProxy& self) {
            return AcKickerFreqProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerFreqProxy& self, py::dict& memo) {
            return AcKickerFreqProxy(self);
          })

      ;

  bind_FTypeArrayND<AcKickerFreqProxyArray1D>(m, "AcKickerFreqStructArray1D");
  bind_FTypeAlloc1D<AcKickerFreqProxyAlloc1D>(m, "AcKickerFreqStructAlloc1D");
  // 2D AcKickerFreqProxy arrays are not used in structs/routines
  // 3D AcKickerFreqProxy arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_struct
void init_ac_kicker_struct(py::module& m, py::class_<AcKickerProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerProxy.amp_vs_time (1D_ALLOC_type -
      .def_property_readonly("amp_vs_time", &AcKickerProxy::amp_vs_time)
      // AcKickerProxy.frequency (1D_ALLOC_type -
      .def_property_readonly("frequency", &AcKickerProxy::frequency)

      .def(
          "__repr__", [](const AcKickerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerProxy& self) {
            return AcKickerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerProxy& self, py::dict& memo) {
            return AcKickerProxy(self);
          })

      ;

  // 1D AcKickerProxy arrays are not used in structs/routines
  // 2D AcKickerProxy arrays are not used in structs/routines
  // 3D AcKickerProxy arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_time_struct
void init_ac_kicker_time_struct(
    py::module& m,
    py::class_<AcKickerTimeProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerTimeProxy.amp (0D_NOT_real -
      .def_property("amp", &AcKickerTimeProxy::amp, &AcKickerTimeProxy::set_amp)
      // AcKickerTimeProxy.time (0D_NOT_real -
      .def_property(
          "time", &AcKickerTimeProxy::time, &AcKickerTimeProxy::set_time)
      // AcKickerTimeProxy.spline (0D_NOT_type -
      .def_property(
          "spline", &AcKickerTimeProxy::spline, &AcKickerTimeProxy::set_spline)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AcKickerTimeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AcKickerTimeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerTimeProxy& self) {
            return AcKickerTimeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerTimeProxy& self, py::dict& memo) {
            return AcKickerTimeProxy(self);
          })

      ;

  bind_FTypeArrayND<AcKickerTimeProxyArray1D>(m, "AcKickerTimeStructArray1D");
  bind_FTypeAlloc1D<AcKickerTimeProxyAlloc1D>(m, "AcKickerTimeStructAlloc1D");
  // 2D AcKickerTimeProxy arrays are not used in structs/routines
  // 3D AcKickerTimeProxy arrays are not used in structs/routines
}

// =============================================================================
// anormal_mode_struct
void init_anormal_mode_struct(
    py::module& m,
    py::class_<AnormalModeProxy>& cls) {
  cls.def(py::init<>())
      // AnormalModeProxy.emittance (0D_NOT_real - Beam emittance (unnormalized). Includes vertical photon opening angle.
      .def_property(
          "emittance",
          &AnormalModeProxy::emittance,
          &AnormalModeProxy::set_emittance)
      // AnormalModeProxy.emittance_no_vert (0D_NOT_real - Unnormalized beam emittance without the vertical photon opening angle taken into account.
      .def_property(
          "emittance_no_vert",
          &AnormalModeProxy::emittance_no_vert,
          &AnormalModeProxy::set_emittance_no_vert)
      // AnormalModeProxy.synch_int (1D_NOT_real - Synchrotron integrals
      .def_property_readonly("synch_int", &AnormalModeProxy::synch_int)
      // AnormalModeProxy.j_damp (0D_NOT_real - damping partition number
      .def_property(
          "j_damp", &AnormalModeProxy::j_damp, &AnormalModeProxy::set_j_damp)
      // AnormalModeProxy.alpha_damp (0D_NOT_real - damping per turn
      .def_property(
          "alpha_damp",
          &AnormalModeProxy::alpha_damp,
          &AnormalModeProxy::set_alpha_damp)
      // AnormalModeProxy.chrom (0D_NOT_real - Chromaticity
      .def_property(
          "chrom", &AnormalModeProxy::chrom, &AnormalModeProxy::set_chrom)
      // AnormalModeProxy.tune (0D_NOT_real - 'Fractional' tune in radians
      .def_property(
          "tune", &AnormalModeProxy::tune, &AnormalModeProxy::set_tune)

      .def(
          "__repr__",
          [](const AnormalModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AnormalModeProxy& self) {
            return AnormalModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AnormalModeProxy& self, py::dict& memo) {
            return AnormalModeProxy(self);
          })

      ;

  // 1D AnormalModeProxy arrays are not used in structs/routines
  // 2D AnormalModeProxy arrays are not used in structs/routines
  // 3D AnormalModeProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_param_struct
void init_aperture_param_struct(
    py::module& m,
    py::class_<ApertureParamProxy>& cls) {
  cls.def(py::init<>())
      // ApertureParamProxy.min_angle (0D_NOT_real -
      .def_property(
          "min_angle",
          &ApertureParamProxy::min_angle,
          &ApertureParamProxy::set_min_angle)
      // ApertureParamProxy.max_angle (0D_NOT_real -
      .def_property(
          "max_angle",
          &ApertureParamProxy::max_angle,
          &ApertureParamProxy::set_max_angle)
      // ApertureParamProxy.n_angle (0D_NOT_integer -
      .def_property(
          "n_angle",
          &ApertureParamProxy::n_angle,
          &ApertureParamProxy::set_n_angle)
      // ApertureParamProxy.n_turn (0D_NOT_integer - Number of turns a particle must survive.
      .def_property(
          "n_turn",
          &ApertureParamProxy::n_turn,
          &ApertureParamProxy::set_n_turn)
      // ApertureParamProxy.x_init (0D_NOT_real - Initial x coordinate to start with for theta_xy = 0.
      .def_property(
          "x_init",
          &ApertureParamProxy::x_init,
          &ApertureParamProxy::set_x_init)
      // ApertureParamProxy.y_init (0D_NOT_real - Initial y coordinate to start with for theta_xy = pi/2.
      .def_property(
          "y_init",
          &ApertureParamProxy::y_init,
          &ApertureParamProxy::set_y_init)
      // ApertureParamProxy.rel_accuracy (0D_NOT_real - Relative resolution of bracketed aperture.
      .def_property(
          "rel_accuracy",
          &ApertureParamProxy::rel_accuracy,
          &ApertureParamProxy::set_rel_accuracy)
      // ApertureParamProxy.abs_accuracy (0D_NOT_real - Absolute resolution of bracketed aperture (meters).
      .def_property(
          "abs_accuracy",
          &ApertureParamProxy::abs_accuracy,
          &ApertureParamProxy::set_abs_accuracy)
      // ApertureParamProxy.start_ele (0D_NOT_character - Element to start tracking at.
      .def_property(
          "start_ele",
          &ApertureParamProxy::start_ele,
          &ApertureParamProxy::set_start_ele)

      .def(
          "__repr__",
          [](const ApertureParamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureParamProxy& self) {
            return ApertureParamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ApertureParamProxy& self, py::dict& memo) {
            return ApertureParamProxy(self);
          })

      ;

  // 1D ApertureParamProxy arrays are not used in structs/routines
  // 2D ApertureParamProxy arrays are not used in structs/routines
  // 3D ApertureParamProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_point_struct
void init_aperture_point_struct(
    py::module& m,
    py::class_<AperturePointProxy>& cls) {
  cls.def(py::init<>())
      // AperturePointProxy.x (0D_NOT_real - (x,y) aperture point with respect to the reference orbit.
      .def_property("x", &AperturePointProxy::x, &AperturePointProxy::set_x)
      // AperturePointProxy.y (0D_NOT_real - (x,y) aperture point with respect to the reference orbit.
      .def_property("y", &AperturePointProxy::y, &AperturePointProxy::set_y)
      // AperturePointProxy.plane (0D_NOT_integer - plane determining loss
      .def_property(
          "plane", &AperturePointProxy::plane, &AperturePointProxy::set_plane)
      // AperturePointProxy.ix_ele (0D_NOT_integer - ele index particle lost at
      .def_property(
          "ix_ele",
          &AperturePointProxy::ix_ele,
          &AperturePointProxy::set_ix_ele)
      // AperturePointProxy.i_turn (0D_NOT_integer - turn particle lost at
      .def_property(
          "i_turn",
          &AperturePointProxy::i_turn,
          &AperturePointProxy::set_i_turn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AperturePointProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AperturePointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AperturePointProxy& self) {
            return AperturePointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AperturePointProxy& self, py::dict& memo) {
            return AperturePointProxy(self);
          })

      ;

  bind_FTypeArrayND<AperturePointProxyArray1D>(m, "AperturePointStructArray1D");
  bind_FTypeAlloc1D<AperturePointProxyAlloc1D>(m, "AperturePointStructAlloc1D");
  // 2D AperturePointProxy arrays are not used in structs/routines
  // 3D AperturePointProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_scan_struct
void init_aperture_scan_struct(
    py::module& m,
    py::class_<ApertureScanProxy>& cls) {
  cls.def(py::init<>())
      // ApertureScanProxy.point (1D_ALLOC_type - Set of aperture points at different angles.
      .def_property_readonly("point", &ApertureScanProxy::point)
      // ApertureScanProxy.ref_orb (0D_NOT_type - Ref orbit around which the scan is made.
      .def_property(
          "ref_orb",
          &ApertureScanProxy::ref_orb,
          &ApertureScanProxy::set_ref_orb)
      // ApertureScanProxy.pz_start (0D_NOT_real - Starting pz.
      .def_property(
          "pz_start",
          &ApertureScanProxy::pz_start,
          &ApertureScanProxy::set_pz_start)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ApertureScanProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ApertureScanProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureScanProxy& self) {
            return ApertureScanProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ApertureScanProxy& self, py::dict& memo) {
            return ApertureScanProxy(self);
          })

      ;

  bind_FTypeArrayND<ApertureScanProxyArray1D>(m, "ApertureScanStructArray1D");
  bind_FTypeAlloc1D<ApertureScanProxyAlloc1D>(m, "ApertureScanStructAlloc1D");
  // 2D ApertureScanProxy arrays are not used in structs/routines
  // 3D ApertureScanProxy arrays are not used in structs/routines
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
// cartesian_map_struct
void init_cartesian_map_struct(
    py::module& m,
    py::class_<CartesianMapProxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CartesianMapProxy::field_scale,
          &CartesianMapProxy::set_field_scale)
      // CartesianMapProxy.r0 (1D_NOT_real - Field origin offset.
      .def_property_readonly("r0", &CartesianMapProxy::r0)
      // CartesianMapProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CartesianMapProxy::master_parameter,
          &CartesianMapProxy::set_master_parameter)
      // CartesianMapProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CartesianMapProxy::ele_anchor_pt,
          &CartesianMapProxy::set_ele_anchor_pt)
      // CartesianMapProxy.field_type (0D_NOT_integer - or electric$
      .def_property(
          "field_type",
          &CartesianMapProxy::field_type,
          &CartesianMapProxy::set_field_type)
      // CartesianMapProxy.ptr (0D_PTR_type -
      .def_property("ptr", &CartesianMapProxy::ptr, &CartesianMapProxy::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CartesianMapProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CartesianMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapProxy& self) {
            return CartesianMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapProxy& self, py::dict& memo) {
            return CartesianMapProxy(self);
          })

      ;

  bind_FTypeArrayND<CartesianMapProxyArray1D>(m, "CartesianMapStructArray1D");
  bind_FTypeAlloc1D<CartesianMapProxyAlloc1D>(m, "CartesianMapStructAlloc1D");
  // 2D CartesianMapProxy arrays are not used in structs/routines
  // 3D CartesianMapProxy arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term1_struct
void init_cartesian_map_term1_struct(
    py::module& m,
    py::class_<CartesianMapTerm1Proxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapTerm1Proxy.coef (0D_NOT_real -
      .def_property(
          "coef",
          &CartesianMapTerm1Proxy::coef,
          &CartesianMapTerm1Proxy::set_coef)
      // CartesianMapTerm1Proxy.kx (0D_NOT_real -
      .def_property(
          "kx", &CartesianMapTerm1Proxy::kx, &CartesianMapTerm1Proxy::set_kx)
      // CartesianMapTerm1Proxy.ky (0D_NOT_real -
      .def_property(
          "ky", &CartesianMapTerm1Proxy::ky, &CartesianMapTerm1Proxy::set_ky)
      // CartesianMapTerm1Proxy.kz (0D_NOT_real -
      .def_property(
          "kz", &CartesianMapTerm1Proxy::kz, &CartesianMapTerm1Proxy::set_kz)
      // CartesianMapTerm1Proxy.x0 (0D_NOT_real -
      .def_property(
          "x0", &CartesianMapTerm1Proxy::x0, &CartesianMapTerm1Proxy::set_x0)
      // CartesianMapTerm1Proxy.y0 (0D_NOT_real -
      .def_property(
          "y0", &CartesianMapTerm1Proxy::y0, &CartesianMapTerm1Proxy::set_y0)
      // CartesianMapTerm1Proxy.phi_z (0D_NOT_real -
      .def_property(
          "phi_z",
          &CartesianMapTerm1Proxy::phi_z,
          &CartesianMapTerm1Proxy::set_phi_z)
      // CartesianMapTerm1Proxy.family (0D_NOT_integer - family_x$, etc.
      .def_property(
          "family",
          &CartesianMapTerm1Proxy::family,
          &CartesianMapTerm1Proxy::set_family)
      // CartesianMapTerm1Proxy.form (0D_NOT_integer - hyper_y$, etc.
      .def_property(
          "form",
          &CartesianMapTerm1Proxy::form,
          &CartesianMapTerm1Proxy::set_form)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CartesianMapTerm1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CartesianMapTerm1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTerm1Proxy& self) {
            return CartesianMapTerm1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapTerm1Proxy& self, py::dict& memo) {
            return CartesianMapTerm1Proxy(self);
          })

      ;

  bind_FTypeArrayND<CartesianMapTerm1ProxyArray1D>(
      m, "CartesianMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CartesianMapTerm1ProxyAlloc1D>(
      m, "CartesianMapTerm1StructAlloc1D");
  // 2D CartesianMapTerm1Proxy arrays are not used in structs/routines
  // 3D CartesianMapTerm1Proxy arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term_struct
void init_cartesian_map_term_struct(
    py::module& m,
    py::class_<CartesianMapTermProxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapTermProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property(
          "file",
          &CartesianMapTermProxy::file,
          &CartesianMapTermProxy::set_file)
      // CartesianMapTermProxy.n_link (0D_NOT_integer - For memory management of %term
      .def_property(
          "n_link",
          &CartesianMapTermProxy::n_link,
          &CartesianMapTermProxy::set_n_link)
      // CartesianMapTermProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &CartesianMapTermProxy::term)

      .def(
          "__repr__",
          [](const CartesianMapTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTermProxy& self) {
            return CartesianMapTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapTermProxy& self, py::dict& memo) {
            return CartesianMapTermProxy(self);
          })

      ;

  // 1D CartesianMapTermProxy arrays are not used in structs/routines
  // 2D CartesianMapTermProxy arrays are not used in structs/routines
  // 3D CartesianMapTermProxy arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_struct
void init_complex_taylor_struct(
    py::module& m,
    py::class_<ComplexTaylorProxy>& cls) {
  cls.def(py::init<>())
      // ComplexTaylorProxy.ref (0D_NOT_complex -
      .def_property(
          "ref", &ComplexTaylorProxy::ref, &ComplexTaylorProxy::set_ref)
      // ComplexTaylorProxy.term (1D_PTR_type -
      .def_property_readonly("term", &ComplexTaylorProxy::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ComplexTaylorProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ComplexTaylorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorProxy& self) {
            return ComplexTaylorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ComplexTaylorProxy& self, py::dict& memo) {
            return ComplexTaylorProxy(self);
          })

      ;

  bind_FTypeArrayND<ComplexTaylorProxyArray1D>(m, "ComplexTaylorStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorProxyAlloc1D>(m, "ComplexTaylorStructAlloc1D");
  // 2D ComplexTaylorProxy arrays are not used in structs/routines
  // 3D ComplexTaylorProxy arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_term_struct
void init_complex_taylor_term_struct(
    py::module& m,
    py::class_<ComplexTaylorTermProxy>& cls) {
  cls.def(py::init<>())
      // ComplexTaylorTermProxy.coef (0D_NOT_complex -
      .def_property(
          "coef",
          &ComplexTaylorTermProxy::coef,
          &ComplexTaylorTermProxy::set_coef)
      // ComplexTaylorTermProxy.expn (1D_NOT_integer -
      .def_property_readonly("expn", &ComplexTaylorTermProxy::expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ComplexTaylorTermProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ComplexTaylorTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorTermProxy& self) {
            return ComplexTaylorTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ComplexTaylorTermProxy& self, py::dict& memo) {
            return ComplexTaylorTermProxy(self);
          })

      ;

  bind_FTypeArrayND<ComplexTaylorTermProxyArray1D>(
      m, "ComplexTaylorTermStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorTermProxyAlloc1D>(
      m, "ComplexTaylorTermStructAlloc1D");
  // 2D ComplexTaylorTermProxy arrays are not used in structs/routines
  // 3D ComplexTaylorTermProxy arrays are not used in structs/routines
}

// =============================================================================
// control_ramp1_struct
void init_control_ramp1_struct(
    py::module& m,
    py::class_<ControlRamp1Proxy>& cls) {
  cls.def(py::init<>())
      // ControlRamp1Proxy.y_knot (1D_ALLOC_real -
      .def_property_readonly("y_knot", &ControlRamp1Proxy::y_knot)
      // ControlRamp1Proxy.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlRamp1Proxy::stack)
      // ControlRamp1Proxy.attribute (0D_NOT_character - Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps.
      .def_property(
          "attribute",
          &ControlRamp1Proxy::attribute,
          &ControlRamp1Proxy::set_attribute)
      // ControlRamp1Proxy.slave_name (0D_NOT_character - Name of slave.
      .def_property(
          "slave_name",
          &ControlRamp1Proxy::slave_name,
          &ControlRamp1Proxy::set_slave_name)
      // ControlRamp1Proxy.is_controller (0D_NOT_logical - Is the slave a controller? If so bookkeeping is different.
      .def_property(
          "is_controller",
          &ControlRamp1Proxy::is_controller,
          &ControlRamp1Proxy::set_is_controller)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ControlRamp1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ControlRamp1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlRamp1Proxy& self) {
            return ControlRamp1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlRamp1Proxy& self, py::dict& memo) {
            return ControlRamp1Proxy(self);
          })

      ;

  bind_FTypeArrayND<ControlRamp1ProxyArray1D>(m, "ControlRamp1StructArray1D");
  bind_FTypeAlloc1D<ControlRamp1ProxyAlloc1D>(m, "ControlRamp1StructAlloc1D");
  // 2D ControlRamp1Proxy arrays are not used in structs/routines
  // 3D ControlRamp1Proxy arrays are not used in structs/routines
}

// =============================================================================
// control_struct
void init_control_struct(py::module& m, py::class_<ControlProxy>& cls) {
  cls.def(py::init<>())
      // ControlProxy.value (0D_NOT_real - Used by group, and overlay elements.
      .def_property("value", &ControlProxy::value, &ControlProxy::set_value)
      // ControlProxy.y_knot (1D_ALLOC_real -
      .def_property_readonly("y_knot", &ControlProxy::y_knot)
      // ControlProxy.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlProxy::stack)
      // ControlProxy.slave (0D_NOT_type -
      .def_property("slave", &ControlProxy::slave, &ControlProxy::set_slave)
      // ControlProxy.lord (0D_NOT_type -
      .def_property("lord", &ControlProxy::lord, &ControlProxy::set_lord)
      // ControlProxy.slave_name (0D_NOT_character - Name of slave.
      .def_property(
          "slave_name",
          &ControlProxy::slave_name,
          &ControlProxy::set_slave_name)
      // ControlProxy.attribute (0D_NOT_character - Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps. Set to 'INPUT' or 'OUTPUT' for feedback slaves.
      .def_property(
          "attribute", &ControlProxy::attribute, &ControlProxy::set_attribute)
      // ControlProxy.ix_attrib (0D_NOT_integer - Index of attribute controlled. See note above!
      .def_property(
          "ix_attrib", &ControlProxy::ix_attrib, &ControlProxy::set_ix_attrib)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ControlProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const ControlProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlProxy& self) {
            return ControlProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlProxy& self, py::dict& memo) {
            return ControlProxy(self);
          })

      ;

  bind_FTypeArrayND<ControlProxyArray1D>(m, "ControlStructArray1D");
  bind_FTypeAlloc1D<ControlProxyAlloc1D>(m, "ControlStructAlloc1D");
  // 2D ControlProxy arrays are not used in structs/routines
  // 3D ControlProxy arrays are not used in structs/routines
}

// =============================================================================
// control_var1_struct
void init_control_var1_struct(
    py::module& m,
    py::class_<ControlVar1Proxy>& cls) {
  cls.def(py::init<>())
      // ControlVar1Proxy.name (0D_NOT_character -
      .def_property(
          "name", &ControlVar1Proxy::name, &ControlVar1Proxy::set_name)
      // ControlVar1Proxy.value (0D_NOT_real -
      .def_property(
          "value", &ControlVar1Proxy::value, &ControlVar1Proxy::set_value)
      // ControlVar1Proxy.old_value (0D_NOT_real -
      .def_property(
          "old_value",
          &ControlVar1Proxy::old_value,
          &ControlVar1Proxy::set_old_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ControlVar1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ControlVar1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlVar1Proxy& self) {
            return ControlVar1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlVar1Proxy& self, py::dict& memo) {
            return ControlVar1Proxy(self);
          })

      ;

  bind_FTypeArrayND<ControlVar1ProxyArray1D>(m, "ControlVar1StructArray1D");
  bind_FTypeAlloc1D<ControlVar1ProxyAlloc1D>(m, "ControlVar1StructAlloc1D");
  // 2D ControlVar1Proxy arrays are not used in structs/routines
  // 3D ControlVar1Proxy arrays are not used in structs/routines
}

// =============================================================================
// controller_struct
void init_controller_struct(py::module& m, py::class_<ControllerProxy>& cls) {
  cls.def(py::init<>())
      // ControllerProxy.var (1D_ALLOC_type -
      .def_property_readonly("var", &ControllerProxy::var)
      // ControllerProxy.ramp (1D_ALLOC_type - For ramper lord elements
      .def_property_readonly("ramp", &ControllerProxy::ramp)
      // ControllerProxy.ramper_lord (1D_ALLOC_type - Ramper lord info for this slave
      .def_property_readonly("ramper_lord", &ControllerProxy::ramper_lord)
      // ControllerProxy.x_knot (1D_ALLOC_real -
      .def_property_readonly("x_knot", &ControllerProxy::x_knot)

      .def(
          "__repr__",
          [](const ControllerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControllerProxy& self) {
            return ControllerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControllerProxy& self, py::dict& memo) {
            return ControllerProxy(self);
          })

      ;

  // 1D ControllerProxy arrays are not used in structs/routines
  // 2D ControllerProxy arrays are not used in structs/routines
  // 3D ControllerProxy arrays are not used in structs/routines
}

// =============================================================================
// coord_array_struct
void init_coord_array_struct(py::module& m, py::class_<CoordArrayProxy>& cls) {
  cls.def(py::init<>())
      // CoordArrayProxy.orbit (1D_ALLOC_type -
      .def_property_readonly("orbit", &CoordArrayProxy::orbit)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordArrayProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CoordArrayProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordArrayProxy& self) {
            return CoordArrayProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CoordArrayProxy& self, py::dict& memo) {
            return CoordArrayProxy(self);
          })

      ;

  bind_FTypeArrayND<CoordArrayProxyArray1D>(m, "CoordArrayStructArray1D");
  bind_FTypeAlloc1D<CoordArrayProxyAlloc1D>(m, "CoordArrayStructAlloc1D");
  // 2D CoordArrayProxy arrays are not used in structs/routines
  // 3D CoordArrayProxy arrays are not used in structs/routines
}

// =============================================================================
// coord_struct
void init_coord_struct(py::module& m, py::class_<CoordProxy>& cls) {
  cls.def(py::init<>())
      // CoordProxy.vec (1D_NOT_real - (x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual.
      .def_property_readonly("vec", &CoordProxy::vec)
      // CoordProxy.s (0D_NOT_real - Longitudinal position
      .def_property("s", &CoordProxy::s, &CoordProxy::set_s)
      // CoordProxy.t (0D_NOT_real16 - Absolute time (not relative to reference). Note: Quad precision!
      .def_property("t", &CoordProxy::t, &CoordProxy::set_t)
      // CoordProxy.spin (1D_NOT_real - Spin.
      .def_property_readonly("spin", &CoordProxy::spin)
      // CoordProxy.field (1D_NOT_real - Photon E-field intensity (x,y).
      .def_property_readonly("field", &CoordProxy::field)
      // CoordProxy.phase (1D_NOT_real - Photon E-field phase (x,y). For charged particles, phase(1) is RF phase.
      .def_property_readonly("phase", &CoordProxy::phase)
      // CoordProxy.charge (0D_NOT_real - Macroparticle weight (which is different from particle species charge). For some space charge calcs the weight is in Coulombs.
      .def_property("charge", &CoordProxy::charge, &CoordProxy::set_charge)
      // CoordProxy.dt_ref (0D_NOT_real - Used in: * time tracking for computing z. * by coherent photons = path_length/c_light.
      .def_property("dt_ref", &CoordProxy::dt_ref, &CoordProxy::set_dt_ref)
      // CoordProxy.r (0D_NOT_real - For general use. Not used by Bmad.
      .def_property("r", &CoordProxy::r, &CoordProxy::set_r)
      // CoordProxy.p0c (0D_NOT_real - For non-photons: Reference momentum. For photons: Photon momentum (not reference).
      .def_property("p0c", &CoordProxy::p0c, &CoordProxy::set_p0c)
      // CoordProxy.E_potential (0D_NOT_real - Potential energy.
      .def_property(
          "E_potential", &CoordProxy::E_potential, &CoordProxy::set_E_potential)
      // CoordProxy.beta (0D_NOT_real - Velocity / c_light.
      .def_property("beta", &CoordProxy::beta, &CoordProxy::set_beta)
      // CoordProxy.ix_ele (0D_NOT_integer - Index of the lattice element the particle is in. May be -1 if element is not associated with a lattice.
      .def_property("ix_ele", &CoordProxy::ix_ele, &CoordProxy::set_ix_ele)
      // CoordProxy.ix_branch (0D_NOT_integer - Index of the lattice branch the particle is in.
      .def_property(
          "ix_branch", &CoordProxy::ix_branch, &CoordProxy::set_ix_branch)
      // CoordProxy.ix_turn (0D_NOT_integer - Turn index for multiturn tracking.
      .def_property("ix_turn", &CoordProxy::ix_turn, &CoordProxy::set_ix_turn)
      // CoordProxy.ix_user (0D_NOT_integer - For general use, not used by Bmad.
      .def_property("ix_user", &CoordProxy::ix_user, &CoordProxy::set_ix_user)
      // CoordProxy.state (0D_NOT_integer - alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc.
      .def_property("state", &CoordProxy::state, &CoordProxy::set_state)
      // CoordProxy.direction (0D_NOT_integer - +1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the element orientation.
      .def_property(
          "direction", &CoordProxy::direction, &CoordProxy::set_direction)
      // CoordProxy.time_dir (0D_NOT_integer - +1 or -1. Time direction. -1 => Traveling backwards in time.
      .def_property(
          "time_dir", &CoordProxy::time_dir, &CoordProxy::set_time_dir)
      // CoordProxy.species (0D_NOT_integer - positron$, proton$, etc.
      .def_property("species", &CoordProxy::species, &CoordProxy::set_species)
      // CoordProxy.location (0D_NOT_integer - upstream_end$, inside$, or downstream_end$
      .def_property(
          "location", &CoordProxy::location, &CoordProxy::set_location)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const CoordProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordProxy& self) {
            return CoordProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CoordProxy& self, py::dict& memo) {
            return CoordProxy(self);
          })

      ;

  bind_FTypeArrayND<CoordProxyArray1D>(m, "CoordStructArray1D");
  bind_FTypeAlloc1D<CoordProxyAlloc1D>(m, "CoordStructAlloc1D");
  // 2D CoordProxy arrays are not used in structs/routines
  // 3D CoordProxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_struct
void init_cylindrical_map_struct(
    py::module& m,
    py::class_<CylindricalMapProxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapProxy.m (0D_NOT_integer - Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)
      .def_property("m", &CylindricalMapProxy::m, &CylindricalMapProxy::set_m)
      // CylindricalMapProxy.harmonic (0D_NOT_integer - Harmonic of fundamental
      .def_property(
          "harmonic",
          &CylindricalMapProxy::harmonic,
          &CylindricalMapProxy::set_harmonic)
      // CylindricalMapProxy.phi0_fieldmap (0D_NOT_real - Mode oscillates as: twopi * (f * t + phi0_fieldmap)
      .def_property(
          "phi0_fieldmap",
          &CylindricalMapProxy::phi0_fieldmap,
          &CylindricalMapProxy::set_phi0_fieldmap)
      // CylindricalMapProxy.theta0_azimuth (0D_NOT_real - Azimuthal ((x, y) plane) orientation of mode.
      .def_property(
          "theta0_azimuth",
          &CylindricalMapProxy::theta0_azimuth,
          &CylindricalMapProxy::set_theta0_azimuth)
      // CylindricalMapProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CylindricalMapProxy::field_scale,
          &CylindricalMapProxy::set_field_scale)
      // CylindricalMapProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CylindricalMapProxy::master_parameter,
          &CylindricalMapProxy::set_master_parameter)
      // CylindricalMapProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CylindricalMapProxy::ele_anchor_pt,
          &CylindricalMapProxy::set_ele_anchor_pt)
      // CylindricalMapProxy.dz (0D_NOT_real - Distance between sampled field points.
      .def_property(
          "dz", &CylindricalMapProxy::dz, &CylindricalMapProxy::set_dz)
      // CylindricalMapProxy.r0 (1D_NOT_real - Field origin offset.
      .def_property_readonly("r0", &CylindricalMapProxy::r0)
      // CylindricalMapProxy.ptr (0D_PTR_type -
      .def_property(
          "ptr", &CylindricalMapProxy::ptr, &CylindricalMapProxy::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CylindricalMapProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CylindricalMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapProxy& self) {
            return CylindricalMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapProxy& self, py::dict& memo) {
            return CylindricalMapProxy(self);
          })

      ;

  bind_FTypeArrayND<CylindricalMapProxyArray1D>(
      m, "CylindricalMapStructArray1D");
  bind_FTypeAlloc1D<CylindricalMapProxyAlloc1D>(
      m, "CylindricalMapStructAlloc1D");
  // 2D CylindricalMapProxy arrays are not used in structs/routines
  // 3D CylindricalMapProxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term1_struct
void init_cylindrical_map_term1_struct(
    py::module& m,
    py::class_<CylindricalMapTerm1Proxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapTerm1Proxy.e_coef (0D_NOT_complex -
      .def_property(
          "e_coef",
          &CylindricalMapTerm1Proxy::e_coef,
          &CylindricalMapTerm1Proxy::set_e_coef)
      // CylindricalMapTerm1Proxy.b_coef (0D_NOT_complex -
      .def_property(
          "b_coef",
          &CylindricalMapTerm1Proxy::b_coef,
          &CylindricalMapTerm1Proxy::set_b_coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CylindricalMapTerm1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CylindricalMapTerm1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTerm1Proxy& self) {
            return CylindricalMapTerm1Proxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapTerm1Proxy& self, py::dict& memo) {
            return CylindricalMapTerm1Proxy(self);
          })

      ;

  bind_FTypeArrayND<CylindricalMapTerm1ProxyArray1D>(
      m, "CylindricalMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CylindricalMapTerm1ProxyAlloc1D>(
      m, "CylindricalMapTerm1StructAlloc1D");
  // 2D CylindricalMapTerm1Proxy arrays are not used in structs/routines
  // 3D CylindricalMapTerm1Proxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term_struct
void init_cylindrical_map_term_struct(
    py::module& m,
    py::class_<CylindricalMapTermProxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapTermProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property(
          "file",
          &CylindricalMapTermProxy::file,
          &CylindricalMapTermProxy::set_file)
      // CylindricalMapTermProxy.n_link (0D_NOT_integer - For memory management of this structure
      .def_property(
          "n_link",
          &CylindricalMapTermProxy::n_link,
          &CylindricalMapTermProxy::set_n_link)
      // CylindricalMapTermProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &CylindricalMapTermProxy::term)

      .def(
          "__repr__",
          [](const CylindricalMapTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTermProxy& self) {
            return CylindricalMapTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapTermProxy& self, py::dict& memo) {
            return CylindricalMapTermProxy(self);
          })

      ;

  // 1D CylindricalMapTermProxy arrays are not used in structs/routines
  // 2D CylindricalMapTermProxy arrays are not used in structs/routines
  // 3D CylindricalMapTermProxy arrays are not used in structs/routines
}

// =============================================================================
// ele_pointer_struct
void init_ele_pointer_struct(py::module& m, py::class_<ElePointerProxy>& cls) {
  cls.def(py::init<>())
      // ElePointerProxy.ele (0D_PTR_type -
      .def_property("ele", &ElePointerProxy::ele, &ElePointerProxy::set_ele)
      // ElePointerProxy.loc (0D_NOT_type -
      .def_property("loc", &ElePointerProxy::loc, &ElePointerProxy::set_loc)
      // ElePointerProxy.id (0D_NOT_integer - For general use. Not used by Bmad. In particular, used by Tao to designate universe ele is in.
      .def_property("id", &ElePointerProxy::id, &ElePointerProxy::set_id)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ElePointerProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ElePointerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ElePointerProxy& self) {
            return ElePointerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ElePointerProxy& self, py::dict& memo) {
            return ElePointerProxy(self);
          })

      ;

  bind_FTypeArrayND<ElePointerProxyArray1D>(m, "ElePointerStructArray1D");
  bind_FTypeAlloc1D<ElePointerProxyAlloc1D>(m, "ElePointerStructAlloc1D");
  // 2D ElePointerProxy arrays are not used in structs/routines
  // 3D ElePointerProxy arrays are not used in structs/routines
}

// =============================================================================
// ele_struct
void init_ele_struct(py::module& m, py::class_<EleProxy>& cls) {
  cls.def(py::init<>())
      // EleProxy.name (0D_NOT_character - name of element.
      .def_property("name", &EleProxy::name, &EleProxy::set_name)
      // EleProxy.type (0D_NOT_character - type name.
      .def_property("type", &EleProxy::type, &EleProxy::set_type)
      // EleProxy.alias (0D_NOT_character - Another name.
      .def_property("alias", &EleProxy::alias, &EleProxy::set_alias)
      // EleProxy.component_name (0D_NOT_character - Used by overlays, multipass patch, etc.
      .def_property(
          "component_name",
          &EleProxy::component_name,
          &EleProxy::set_component_name)
      // EleProxy.descrip (0D_PTR_character - Description string.
      .def_property("descrip", &EleProxy::descrip, &EleProxy::set_descrip)
      // EleProxy.a (0D_NOT_type - Twiss parameters at end of element
      .def_property("a", &EleProxy::a, &EleProxy::set_a)
      // EleProxy.b (0D_NOT_type - Twiss parameters at end of element
      .def_property("b", &EleProxy::b, &EleProxy::set_b)
      // EleProxy.z (0D_NOT_type - Twiss parameters at end of element
      .def_property("z", &EleProxy::z, &EleProxy::set_z)
      // EleProxy.x (0D_NOT_type - Projected dispersions.
      .def_property("x", &EleProxy::x, &EleProxy::set_x)
      // EleProxy.y (0D_NOT_type - Projected dispersions.
      .def_property("y", &EleProxy::y, &EleProxy::set_y)
      // EleProxy.ac_kick (0D_PTR_type - ac_kicker element parameters.
      .def_property("ac_kick", &EleProxy::ac_kick, &EleProxy::set_ac_kick)
      // EleProxy.bookkeeping_state (0D_NOT_type - Attribute bookkeeping
      .def_property(
          "bookkeeping_state",
          &EleProxy::bookkeeping_state,
          &EleProxy::set_bookkeeping_state)
      // EleProxy.branch (0D_PTR_type - Pointer to branch containing element.
      .def_property("branch", &EleProxy::branch, &EleProxy::set_branch)
      // EleProxy.control (0D_PTR_type - group & overlay variables.
      .def_property("control", &EleProxy::control, &EleProxy::set_control)
      // EleProxy.rf (0D_PTR_type - RF parameters.
      .def_property("rf", &EleProxy::rf, &EleProxy::set_rf)
      // EleProxy.lord (0D_PTR_type - Pointer to a slice lord.
      .def_property("lord", &EleProxy::lord, &EleProxy::set_lord)
      // EleProxy.floor (0D_NOT_type -
      .def_property("floor", &EleProxy::floor, &EleProxy::set_floor)
      // EleProxy.high_energy_space_charge (0D_PTR_type -
      .def_property(
          "high_energy_space_charge",
          &EleProxy::high_energy_space_charge,
          &EleProxy::set_high_energy_space_charge)
      // EleProxy.mode3 (0D_PTR_type - 6D normal mode structure.
      .def_property("mode3", &EleProxy::mode3, &EleProxy::set_mode3)
      // EleProxy.photon (0D_PTR_type -
      .def_property("photon", &EleProxy::photon, &EleProxy::set_photon)
      // EleProxy.rad_map (0D_PTR_type - Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map.
      .def_property("rad_map", &EleProxy::rad_map, &EleProxy::set_rad_map)
      // EleProxy.taylor (1D_NOT_type - Phase space Taylor map.
      .def_property_readonly("taylor", &EleProxy::taylor)
      // EleProxy.spin_taylor_ref_orb_in (1D_NOT_real -
      .def_property_readonly(
          "spin_taylor_ref_orb_in", &EleProxy::spin_taylor_ref_orb_in)
      // EleProxy.spin_taylor (1D_NOT_type - Quaternion Spin Taylor map.
      .def_property_readonly("spin_taylor", &EleProxy::spin_taylor)
      // EleProxy.wake (0D_PTR_type - Wakes
      .def_property("wake", &EleProxy::wake, &EleProxy::set_wake)
      // EleProxy.wall3d (1D_PTR_type - Chamber or capillary wall E/M field structs.
      .def_property_readonly("wall3d", &EleProxy::wall3d)
      // EleProxy.cartesian_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cartesian_map", &EleProxy::cartesian_map)
      // EleProxy.cylindrical_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cylindrical_map", &EleProxy::cylindrical_map)
      // EleProxy.gen_grad_map (1D_PTR_type - Used to define E/M fields.
      .def_property_readonly("gen_grad_map", &EleProxy::gen_grad_map)
      // EleProxy.grid_field (1D_PTR_type - Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference particle which is generally the zero orbit (non-zero, for example, in the second slice of a sliced wiggler).
      .def_property_readonly("grid_field", &EleProxy::grid_field)
      // EleProxy.map_ref_orb_in (0D_NOT_type - Entrance end transfer map ref orbit
      .def_property(
          "map_ref_orb_in",
          &EleProxy::map_ref_orb_in,
          &EleProxy::set_map_ref_orb_in)
      // EleProxy.map_ref_orb_out (0D_NOT_type - Exit end transfer map ref orbit
      .def_property(
          "map_ref_orb_out",
          &EleProxy::map_ref_orb_out,
          &EleProxy::set_map_ref_orb_out)
      // EleProxy.time_ref_orb_in (0D_NOT_type - Reference orbit at entrance end for ref_time calc.
      .def_property(
          "time_ref_orb_in",
          &EleProxy::time_ref_orb_in,
          &EleProxy::set_time_ref_orb_in)
      // EleProxy.time_ref_orb_out (0D_NOT_type - Reference orbit at exit end for ref_time calc.
      .def_property(
          "time_ref_orb_out",
          &EleProxy::time_ref_orb_out,
          &EleProxy::set_time_ref_orb_out)
      // EleProxy.value (1D_NOT_real - attribute values.
      .def_property_readonly("value", &EleProxy::value)
      // EleProxy.old_value (1D_NOT_real - Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit matrices is %map_ref_orb_in/out
      .def_property_readonly("old_value", &EleProxy::old_value)
      // EleProxy.spin_q (2D_NOT_real - 0th and 1st order Spin transport quaternion.
      .def_property_readonly("spin_q", &EleProxy::spin_q)
      // EleProxy.vec0 (1D_NOT_real - 0th order transport vector.
      .def_property_readonly("vec0", &EleProxy::vec0)
      // EleProxy.mat6 (2D_NOT_real - 1st order transport matrix.
      .def_property_readonly("mat6", &EleProxy::mat6)
      // EleProxy.c_mat (2D_NOT_real - 2x2 C coupling matrix
      .def_property_readonly("c_mat", &EleProxy::c_mat)
      // EleProxy.dc_mat_dpz (2D_NOT_real - d(c_mat)/dpz variation.
      .def_property_readonly("dc_mat_dpz", &EleProxy::dc_mat_dpz)
      // EleProxy.gamma_c (0D_NOT_real - gamma associated with C matrix
      .def_property("gamma_c", &EleProxy::gamma_c, &EleProxy::set_gamma_c)
      // EleProxy.s_start (0D_NOT_real - longitudinal ref position at entrance_end
      .def_property("s_start", &EleProxy::s_start, &EleProxy::set_s_start)
      // EleProxy.s (0D_NOT_real - longitudinal ref position at the exit end.
      .def_property("s", &EleProxy::s, &EleProxy::set_s)
      // EleProxy.ref_time (0D_NOT_real - Time ref particle passes exit end.
      .def_property("ref_time", &EleProxy::ref_time, &EleProxy::set_ref_time)
      // EleProxy.a_pole (1D_PTR_real - knl for multipole elements.
      .def_property_readonly("a_pole", &EleProxy::a_pole)
      // EleProxy.b_pole (1D_PTR_real - tilt for multipole elements.
      .def_property_readonly("b_pole", &EleProxy::b_pole)
      // EleProxy.a_pole_elec (1D_PTR_real - Electrostatic multipoles. ksnl for multipole elements.
      .def_property_readonly("a_pole_elec", &EleProxy::a_pole_elec)
      // EleProxy.b_pole_elec (1D_PTR_real - Electrostatic multipoles.
      .def_property_readonly("b_pole_elec", &EleProxy::b_pole_elec)
      // EleProxy.custom (1D_PTR_real - Custom attributes.
      .def_property_readonly("custom", &EleProxy::custom)
      // EleProxy.r (3D_PTR_real - For general use. Not used by Bmad.
      .def_property_readonly("r", &EleProxy::r)
      // EleProxy.key (0D_NOT_integer - Element class (quadrupole, etc.).
      .def_property("key", &EleProxy::key, &EleProxy::set_key)
      // EleProxy.sub_key (0D_NOT_integer - Records bend input type.
      .def_property("sub_key", &EleProxy::sub_key, &EleProxy::set_sub_key)
      // EleProxy.ix_ele (0D_NOT_integer - Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements.
      .def_property("ix_ele", &EleProxy::ix_ele, &EleProxy::set_ix_ele)
      // EleProxy.ix_branch (0D_NOT_integer - Index in lat%branch(:) array. Note: lat%ele => lat%branch(0).
      .def_property("ix_branch", &EleProxy::ix_branch, &EleProxy::set_ix_branch)
      // EleProxy.lord_status (0D_NOT_integer - Type of lord element this is. overlay_lord$, etc.
      .def_property(
          "lord_status", &EleProxy::lord_status, &EleProxy::set_lord_status)
      // EleProxy.n_slave (0D_NOT_integer - Number of slaves (except field overlap slaves) of this element.
      .def_property("n_slave", &EleProxy::n_slave, &EleProxy::set_n_slave)
      // EleProxy.n_slave_field (0D_NOT_integer - Number of field slaves of this element.
      .def_property(
          "n_slave_field",
          &EleProxy::n_slave_field,
          &EleProxy::set_n_slave_field)
      // EleProxy.ix1_slave (0D_NOT_integer - Pointer index to this element's slaves.
      .def_property("ix1_slave", &EleProxy::ix1_slave, &EleProxy::set_ix1_slave)
      // EleProxy.slave_status (0D_NOT_integer - Type of slave element this is. multipass_slave$, slice_slave$, etc.
      .def_property(
          "slave_status", &EleProxy::slave_status, &EleProxy::set_slave_status)
      // EleProxy.n_lord (0D_NOT_integer - Number of lords (except field overlap and ramper lords).
      .def_property("n_lord", &EleProxy::n_lord, &EleProxy::set_n_lord)
      // EleProxy.n_lord_field (0D_NOT_integer - Number of field lords of this element.
      .def_property(
          "n_lord_field", &EleProxy::n_lord_field, &EleProxy::set_n_lord_field)
      // EleProxy.n_lord_ramper (0D_NOT_integer - Number of ramper lords.
      .def_property(
          "n_lord_ramper",
          &EleProxy::n_lord_ramper,
          &EleProxy::set_n_lord_ramper)
      // EleProxy.ic1_lord (0D_NOT_integer - Pointer index to this element's lords.
      .def_property("ic1_lord", &EleProxy::ic1_lord, &EleProxy::set_ic1_lord)
      // EleProxy.ix_pointer (0D_NOT_integer - For general use. Not used by Bmad.
      .def_property(
          "ix_pointer", &EleProxy::ix_pointer, &EleProxy::set_ix_pointer)
      // EleProxy.ixx (0D_NOT_integer - Index for Bmad internal use.
      .def_property("ixx", &EleProxy::ixx, &EleProxy::set_ixx)
      // EleProxy.iyy (0D_NOT_integer - Index for Bmad internal use.
      .def_property("iyy", &EleProxy::iyy, &EleProxy::set_iyy)
      // EleProxy.izz (0D_NOT_integer - Index for Bmad internal use.
      .def_property("izz", &EleProxy::izz, &EleProxy::set_izz)
      // EleProxy.mat6_calc_method (0D_NOT_integer - taylor$, symp_lie_ptc$, etc.
      .def_property(
          "mat6_calc_method",
          &EleProxy::mat6_calc_method,
          &EleProxy::set_mat6_calc_method)
      // EleProxy.tracking_method (0D_NOT_integer - taylor$, linear$, etc.
      .def_property(
          "tracking_method",
          &EleProxy::tracking_method,
          &EleProxy::set_tracking_method)
      // EleProxy.spin_tracking_method (0D_NOT_integer - symp_lie_ptc$, etc.
      .def_property(
          "spin_tracking_method",
          &EleProxy::spin_tracking_method,
          &EleProxy::set_spin_tracking_method)
      // EleProxy.csr_method (0D_NOT_integer - or one_dim$ ('1_dim'), steady_state_3d$
      .def_property(
          "csr_method", &EleProxy::csr_method, &EleProxy::set_csr_method)
      // EleProxy.space_charge_method (0D_NOT_integer - slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$
      .def_property(
          "space_charge_method",
          &EleProxy::space_charge_method,
          &EleProxy::set_space_charge_method)
      // EleProxy.ptc_integration_type (0D_NOT_integer - drift_kick$, matrix_kick$, or ripken_kick$
      .def_property(
          "ptc_integration_type",
          &EleProxy::ptc_integration_type,
          &EleProxy::set_ptc_integration_type)
      // EleProxy.field_calc (0D_NOT_integer - no_field$, fieldmap$, refer_to_lords$, or custom$
      .def_property(
          "field_calc", &EleProxy::field_calc, &EleProxy::set_field_calc)
      // EleProxy.aperture_at (0D_NOT_integer - Aperture location: entrance_end$, ...
      .def_property(
          "aperture_at", &EleProxy::aperture_at, &EleProxy::set_aperture_at)
      // EleProxy.aperture_type (0D_NOT_integer - rectangular$, elliptical$, auto_aperture$, ...
      .def_property(
          "aperture_type",
          &EleProxy::aperture_type,
          &EleProxy::set_aperture_type)
      // EleProxy.ref_species (0D_NOT_integer - Reference species
      .def_property(
          "ref_species", &EleProxy::ref_species, &EleProxy::set_ref_species)
      // EleProxy.orientation (0D_NOT_integer - -1 -> Element is longitudinally reversed. +1 -> Normal.
      .def_property(
          "orientation", &EleProxy::orientation, &EleProxy::set_orientation)
      // EleProxy.symplectify (0D_NOT_logical - Symplectify mat6 matrices.
      .def_property(
          "symplectify", &EleProxy::symplectify, &EleProxy::set_symplectify)
      // EleProxy.mode_flip (0D_NOT_logical - Have the normal modes traded places?
      .def_property("mode_flip", &EleProxy::mode_flip, &EleProxy::set_mode_flip)
      // EleProxy.multipoles_on (0D_NOT_logical - For turning multipoles on/off
      .def_property(
          "multipoles_on",
          &EleProxy::multipoles_on,
          &EleProxy::set_multipoles_on)
      // EleProxy.scale_multipoles (0D_NOT_logical - Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the element?
      .def_property(
          "scale_multipoles",
          &EleProxy::scale_multipoles,
          &EleProxy::set_scale_multipoles)
      // EleProxy.taylor_map_includes_offsets (0D_NOT_logical - Taylor map calculated with element misalignments?
      .def_property(
          "taylor_map_includes_offsets",
          &EleProxy::taylor_map_includes_offsets,
          &EleProxy::set_taylor_map_includes_offsets)
      // EleProxy.field_master (0D_NOT_logical - Calculate strength from the field value?
      .def_property(
          "field_master", &EleProxy::field_master, &EleProxy::set_field_master)
      // EleProxy.is_on (0D_NOT_logical - For turning element on/off.
      .def_property("is_on", &EleProxy::is_on, &EleProxy::set_is_on)
      // EleProxy.logic (0D_NOT_logical - For general use. Not used by Bmad (except during lattice parsing).
      .def_property("logic", &EleProxy::logic, &EleProxy::set_logic)
      // EleProxy.bmad_logic (0D_NOT_logical - For Bmad internal use only.
      .def_property(
          "bmad_logic", &EleProxy::bmad_logic, &EleProxy::set_bmad_logic)
      // EleProxy.select (0D_NOT_logical - For Bmad internal use only.
      .def_property("select", &EleProxy::select, &EleProxy::set_select)
      // EleProxy.offset_moves_aperture (0D_NOT_logical - element offsets affects aperture? ! final :: ele_finalizer
      .def_property(
          "offset_moves_aperture",
          &EleProxy::offset_moves_aperture,
          &EleProxy::set_offset_moves_aperture)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EleProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const EleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EleProxy& self) {
            return EleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EleProxy& self, py::dict& memo) { return EleProxy(self); })

      ;

  bind_FTypeArrayND<EleProxyArray1D>(m, "EleStructArray1D");
  bind_FTypeAlloc1D<EleProxyAlloc1D>(m, "EleStructAlloc1D");
  // 2D EleProxy arrays are not used in structs/routines
  // 3D EleProxy arrays are not used in structs/routines
}

// =============================================================================
// ellipse_beam_init_struct
void init_ellipse_beam_init_struct(
    py::module& m,
    py::class_<EllipseBeamInitProxy>& cls) {
  cls.def(py::init<>())
      // EllipseBeamInitProxy.part_per_ellipse (0D_NOT_integer - number of particles per ellipse
      .def_property(
          "part_per_ellipse",
          &EllipseBeamInitProxy::part_per_ellipse,
          &EllipseBeamInitProxy::set_part_per_ellipse)
      // EllipseBeamInitProxy.n_ellipse (0D_NOT_integer - number of ellipses (>= 1)
      .def_property(
          "n_ellipse",
          &EllipseBeamInitProxy::n_ellipse,
          &EllipseBeamInitProxy::set_n_ellipse)
      // EllipseBeamInitProxy.sigma_cutoff (0D_NOT_real - sigma cutoff of the representation
      .def_property(
          "sigma_cutoff",
          &EllipseBeamInitProxy::sigma_cutoff,
          &EllipseBeamInitProxy::set_sigma_cutoff)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return EllipseBeamInitProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const EllipseBeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EllipseBeamInitProxy& self) {
            return EllipseBeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EllipseBeamInitProxy& self, py::dict& memo) {
            return EllipseBeamInitProxy(self);
          })

      ;

  bind_FTypeArrayND<EllipseBeamInitProxyArray1D>(
      m, "EllipseBeamInitStructArray1D");
  bind_FTypeAlloc1D<EllipseBeamInitProxyAlloc1D>(
      m, "EllipseBeamInitStructAlloc1D");
  // 2D EllipseBeamInitProxy arrays are not used in structs/routines
  // 3D EllipseBeamInitProxy arrays are not used in structs/routines
}

// =============================================================================
// em_field_struct
void init_em_field_struct(py::module& m, py::class_<EmFieldProxy>& cls) {
  cls.def(py::init<>())
      // EmFieldProxy.E (1D_NOT_real - electric field.
      .def_property_readonly("E", &EmFieldProxy::E)
      // EmFieldProxy.B (1D_NOT_real - magnetic field.
      .def_property_readonly("B", &EmFieldProxy::B)
      // EmFieldProxy.dE (2D_NOT_real - electric field gradient.
      .def_property_readonly("dE", &EmFieldProxy::dE)
      // EmFieldProxy.dB (2D_NOT_real - magnetic field gradient.
      .def_property_readonly("dB", &EmFieldProxy::dB)
      // EmFieldProxy.phi (0D_NOT_real - Electric scalar potential.
      .def_property("phi", &EmFieldProxy::phi, &EmFieldProxy::set_phi)
      // EmFieldProxy.phi_B (0D_NOT_real - Magnetic scalar potential.
      .def_property("phi_B", &EmFieldProxy::phi_B, &EmFieldProxy::set_phi_B)
      // EmFieldProxy.A (1D_NOT_real - Magnetic vector potential.
      .def_property_readonly("A", &EmFieldProxy::A)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmFieldProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const EmFieldProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmFieldProxy& self) {
            return EmFieldProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmFieldProxy& self, py::dict& memo) {
            return EmFieldProxy(self);
          })

      ;

  bind_FTypeArrayND<EmFieldProxyArray1D>(m, "EmFieldStructArray1D");
  bind_FTypeAlloc1D<EmFieldProxyAlloc1D>(m, "EmFieldStructAlloc1D");
  // 2D EmFieldProxy arrays are not used in structs/routines
  // 3D EmFieldProxy arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_struct
void init_em_taylor_struct(py::module& m, py::class_<EmTaylorProxy>& cls) {
  cls.def(py::init<>())
      // EmTaylorProxy.ref (0D_NOT_real -
      .def_property("ref", &EmTaylorProxy::ref, &EmTaylorProxy::set_ref)
      // EmTaylorProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &EmTaylorProxy::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmTaylorProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const EmTaylorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorProxy& self) {
            return EmTaylorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmTaylorProxy& self, py::dict& memo) {
            return EmTaylorProxy(self);
          })

      ;

  bind_FTypeArrayND<EmTaylorProxyArray1D>(m, "EmTaylorStructArray1D");
  bind_FTypeAlloc1D<EmTaylorProxyAlloc1D>(m, "EmTaylorStructAlloc1D");
  // 2D EmTaylorProxy arrays are not used in structs/routines
  // 3D EmTaylorProxy arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_term_struct
void init_em_taylor_term_struct(
    py::module& m,
    py::class_<EmTaylorTermProxy>& cls) {
  cls.def(py::init<>())
      // EmTaylorTermProxy.coef (0D_NOT_real -
      .def_property(
          "coef", &EmTaylorTermProxy::coef, &EmTaylorTermProxy::set_coef)
      // EmTaylorTermProxy.expn (1D_NOT_integer -
      .def_property_readonly("expn", &EmTaylorTermProxy::expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return EmTaylorTermProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const EmTaylorTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorTermProxy& self) {
            return EmTaylorTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmTaylorTermProxy& self, py::dict& memo) {
            return EmTaylorTermProxy(self);
          })

      ;

  bind_FTypeArrayND<EmTaylorTermProxyArray1D>(m, "EmTaylorTermStructArray1D");
  bind_FTypeAlloc1D<EmTaylorTermProxyAlloc1D>(m, "EmTaylorTermStructAlloc1D");
  // 2D EmTaylorTermProxy arrays are not used in structs/routines
  // 3D EmTaylorTermProxy arrays are not used in structs/routines
}

// =============================================================================
// expression_atom_struct
void init_expression_atom_struct(
    py::module& m,
    py::class_<ExpressionAtomProxy>& cls) {
  cls.def(py::init<>())
      // ExpressionAtomProxy.name (0D_NOT_character -
      .def_property(
          "name", &ExpressionAtomProxy::name, &ExpressionAtomProxy::set_name)
      // ExpressionAtomProxy.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name
      .def_property(
          "type", &ExpressionAtomProxy::type, &ExpressionAtomProxy::set_type)
      // ExpressionAtomProxy.value (0D_NOT_real -
      .def_property(
          "value", &ExpressionAtomProxy::value, &ExpressionAtomProxy::set_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ExpressionAtomProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ExpressionAtomProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionAtomProxy& self) {
            return ExpressionAtomProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ExpressionAtomProxy& self, py::dict& memo) {
            return ExpressionAtomProxy(self);
          })

      ;

  bind_FTypeArrayND<ExpressionAtomProxyArray1D>(
      m, "ExpressionAtomStructArray1D");
  bind_FTypeAlloc1D<ExpressionAtomProxyAlloc1D>(
      m, "ExpressionAtomStructAlloc1D");
  // 2D ExpressionAtomProxy arrays are not used in structs/routines
  // 3D ExpressionAtomProxy arrays are not used in structs/routines
}

// =============================================================================
// expression_tree_struct
void init_expression_tree_struct(
    py::module& m,
    py::class_<ExpressionTreeProxy>& cls) {
  cls.def(py::init<>())
      // ExpressionTreeProxy.name (0D_NOT_character -
      .def_property(
          "name", &ExpressionTreeProxy::name, &ExpressionTreeProxy::set_name)
      // ExpressionTreeProxy.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc.
      .def_property(
          "type", &ExpressionTreeProxy::type, &ExpressionTreeProxy::set_type)
      // ExpressionTreeProxy.value (0D_NOT_real -
      .def_property(
          "value", &ExpressionTreeProxy::value, &ExpressionTreeProxy::set_value)
      // ExpressionTreeProxy.node (1D_PTR_type - Child nodes. Note: Pointer used here since Ifort does not support allocatable.
      .def_property_readonly("node", &ExpressionTreeProxy::node)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ExpressionTreeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ExpressionTreeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionTreeProxy& self) {
            return ExpressionTreeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ExpressionTreeProxy& self, py::dict& memo) {
            return ExpressionTreeProxy(self);
          })

      ;

  bind_FTypeArrayND<ExpressionTreeProxyArray1D>(
      m, "ExpressionTreeStructArray1D");
  bind_FTypeAlloc1D<ExpressionTreeProxyAlloc1D>(
      m, "ExpressionTreeStructAlloc1D");
  // 2D ExpressionTreeProxy arrays are not used in structs/routines
  // 3D ExpressionTreeProxy arrays are not used in structs/routines
}

// =============================================================================
// floor_position_struct
void init_floor_position_struct(
    py::module& m,
    py::class_<FloorPositionProxy>& cls) {
  cls.def(py::init<>())
      // FloorPositionProxy.r (1D_NOT_real - (x, y, z) offset from origin
      .def_property_readonly("r", &FloorPositionProxy::r)
      // FloorPositionProxy.w (2D_NOT_real - W matrix. Columns are unit vectors of the frame axes.
      .def_property_readonly("w", &FloorPositionProxy::w)
      // FloorPositionProxy.theta (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "theta", &FloorPositionProxy::theta, &FloorPositionProxy::set_theta)
      // FloorPositionProxy.phi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "phi", &FloorPositionProxy::phi, &FloorPositionProxy::set_phi)
      // FloorPositionProxy.psi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "psi", &FloorPositionProxy::psi, &FloorPositionProxy::set_psi)

      .def(
          "__repr__",
          [](const FloorPositionProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FloorPositionProxy& self) {
            return FloorPositionProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const FloorPositionProxy& self, py::dict& memo) {
            return FloorPositionProxy(self);
          })

      ;

  // 1D FloorPositionProxy arrays are not used in structs/routines
  // 2D FloorPositionProxy arrays are not used in structs/routines
  // 3D FloorPositionProxy arrays are not used in structs/routines
}

// =============================================================================
// gen_grad1_struct
void init_gen_grad1_struct(py::module& m, py::class_<GenGrad1Proxy>& cls) {
  cls.def(py::init<>())
      // GenGrad1Proxy.m (0D_NOT_integer - Azimuthal index
      .def_property("m", &GenGrad1Proxy::m, &GenGrad1Proxy::set_m)
      // GenGrad1Proxy.sincos (0D_NOT_integer - sin$ or cos$
      .def_property(
          "sincos", &GenGrad1Proxy::sincos, &GenGrad1Proxy::set_sincos)
      // GenGrad1Proxy.n_deriv_max (0D_NOT_integer - Max GG derivative The derivative matrix is extended to include the interpolating spline polynomial.
      .def_property(
          "n_deriv_max",
          &GenGrad1Proxy::n_deriv_max,
          &GenGrad1Proxy::set_n_deriv_max)
      // GenGrad1Proxy.deriv (2D_ALLOC_real - Range: (iz0:iz1, 0:2*n_deriv_max+1)
      .def_property_readonly("deriv", &GenGrad1Proxy::deriv)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return GenGrad1ProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const GenGrad1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GenGrad1Proxy& self) {
            return GenGrad1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GenGrad1Proxy& self, py::dict& memo) {
            return GenGrad1Proxy(self);
          })

      ;

  bind_FTypeArrayND<GenGrad1ProxyArray1D>(m, "GenGrad1StructArray1D");
  bind_FTypeAlloc1D<GenGrad1ProxyAlloc1D>(m, "GenGrad1StructAlloc1D");
  // 2D GenGrad1Proxy arrays are not used in structs/routines
  // 3D GenGrad1Proxy arrays are not used in structs/routines
}

// =============================================================================
// gen_grad_map_struct
void init_gen_grad_map_struct(py::module& m, py::class_<GenGradMapProxy>& cls) {
  cls.def(py::init<>())
      // GenGradMapProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property("file", &GenGradMapProxy::file, &GenGradMapProxy::set_file)
      // GenGradMapProxy.gg (1D_ALLOC_type -
      .def_property_readonly("gg", &GenGradMapProxy::gg)
      // GenGradMapProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &GenGradMapProxy::ele_anchor_pt,
          &GenGradMapProxy::set_ele_anchor_pt)
      // GenGradMapProxy.field_type (0D_NOT_integer - or electric$
      .def_property(
          "field_type",
          &GenGradMapProxy::field_type,
          &GenGradMapProxy::set_field_type)
      // GenGradMapProxy.iz0 (0D_NOT_integer - gg%deriv(iz0:iz1, :) lower bound.
      .def_property("iz0", &GenGradMapProxy::iz0, &GenGradMapProxy::set_iz0)
      // GenGradMapProxy.iz1 (0D_NOT_integer - gg%deriv(iz0:iz1, :) upper bound.
      .def_property("iz1", &GenGradMapProxy::iz1, &GenGradMapProxy::set_iz1)
      // GenGradMapProxy.dz (0D_NOT_real - Point spacing.
      .def_property("dz", &GenGradMapProxy::dz, &GenGradMapProxy::set_dz)
      // GenGradMapProxy.r0 (1D_NOT_real - field origin relative to ele_anchor_pt.
      .def_property_readonly("r0", &GenGradMapProxy::r0)
      // GenGradMapProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &GenGradMapProxy::field_scale,
          &GenGradMapProxy::set_field_scale)
      // GenGradMapProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &GenGradMapProxy::master_parameter,
          &GenGradMapProxy::set_master_parameter)
      // GenGradMapProxy.curved_ref_frame (0D_NOT_logical -
      .def_property(
          "curved_ref_frame",
          &GenGradMapProxy::curved_ref_frame,
          &GenGradMapProxy::set_curved_ref_frame)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return GenGradMapProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const GenGradMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GenGradMapProxy& self) {
            return GenGradMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GenGradMapProxy& self, py::dict& memo) {
            return GenGradMapProxy(self);
          })

      ;

  bind_FTypeArrayND<GenGradMapProxyArray1D>(m, "GenGradMapStructArray1D");
  bind_FTypeAlloc1D<GenGradMapProxyAlloc1D>(m, "GenGradMapStructAlloc1D");
  // 2D GenGradMapProxy arrays are not used in structs/routines
  // 3D GenGradMapProxy arrays are not used in structs/routines
}

// =============================================================================
// grid_beam_init_struct
void init_grid_beam_init_struct(
    py::module& m,
    py::class_<GridBeamInitProxy>& cls) {
  cls.def(py::init<>())
      // GridBeamInitProxy.n_x (0D_NOT_integer - Number of columns.
      .def_property("n_x", &GridBeamInitProxy::n_x, &GridBeamInitProxy::set_n_x)
      // GridBeamInitProxy.n_px (0D_NOT_integer - Number of rows.
      .def_property(
          "n_px", &GridBeamInitProxy::n_px, &GridBeamInitProxy::set_n_px)
      // GridBeamInitProxy.x_min (0D_NOT_real - Lower x limit.
      .def_property(
          "x_min", &GridBeamInitProxy::x_min, &GridBeamInitProxy::set_x_min)
      // GridBeamInitProxy.x_max (0D_NOT_real - Upper x limit.
      .def_property(
          "x_max", &GridBeamInitProxy::x_max, &GridBeamInitProxy::set_x_max)
      // GridBeamInitProxy.px_min (0D_NOT_real - Lower px limit.
      .def_property(
          "px_min", &GridBeamInitProxy::px_min, &GridBeamInitProxy::set_px_min)
      // GridBeamInitProxy.px_max (0D_NOT_real - Upper px limit.
      .def_property(
          "px_max", &GridBeamInitProxy::px_max, &GridBeamInitProxy::set_px_max)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return GridBeamInitProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const GridBeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridBeamInitProxy& self) {
            return GridBeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GridBeamInitProxy& self, py::dict& memo) {
            return GridBeamInitProxy(self);
          })

      ;

  bind_FTypeArrayND<GridBeamInitProxyArray1D>(m, "GridBeamInitStructArray1D");
  bind_FTypeAlloc1D<GridBeamInitProxyAlloc1D>(m, "GridBeamInitStructAlloc1D");
  // 2D GridBeamInitProxy arrays are not used in structs/routines
  // 3D GridBeamInitProxy arrays are not used in structs/routines
}

// =============================================================================
// grid_field_pt1_struct
void init_grid_field_pt1_struct(
    py::module& m,
    py::class_<GridFieldPt1Proxy>& cls) {
  cls.def(py::init<>())
      // GridFieldPt1Proxy.E (1D_NOT_complex -
      .def_property_readonly("E", &GridFieldPt1Proxy::E)
      // GridFieldPt1Proxy.B (1D_NOT_complex -
      .def_property_readonly("B", &GridFieldPt1Proxy::B)

      .def(
          "__repr__",
          [](const GridFieldPt1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldPt1Proxy& self) {
            return GridFieldPt1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GridFieldPt1Proxy& self, py::dict& memo) {
            return GridFieldPt1Proxy(self);
          })

      ;

  // 1D GridFieldPt1Proxy arrays are not used in structs/routines
  // 2D GridFieldPt1Proxy arrays are not used in structs/routines
  bind_FTypeArrayND<GridFieldPt1ProxyArray3D>(m, "GridFieldPt1StructArray3D");
}

// =============================================================================
// grid_field_pt_struct
void init_grid_field_pt_struct(
    py::module& m,
    py::class_<GridFieldPtProxy>& cls) {
  cls.def(py::init<>())
      // GridFieldPtProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property(
          "file", &GridFieldPtProxy::file, &GridFieldPtProxy::set_file)
      // GridFieldPtProxy.n_link (0D_NOT_integer - For memory management of this structure
      .def_property(
          "n_link", &GridFieldPtProxy::n_link, &GridFieldPtProxy::set_n_link)
      // GridFieldPtProxy.pt (3D_ALLOC_type -
      .def_property_readonly("pt", &GridFieldPtProxy::pt)

      .def(
          "__repr__",
          [](const GridFieldPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldPtProxy& self) {
            return GridFieldPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GridFieldPtProxy& self, py::dict& memo) {
            return GridFieldPtProxy(self);
          })

      ;

  // 1D GridFieldPtProxy arrays are not used in structs/routines
  // 2D GridFieldPtProxy arrays are not used in structs/routines
  // 3D GridFieldPtProxy arrays are not used in structs/routines
}

// =============================================================================
// grid_field_struct
void init_grid_field_struct(py::module& m, py::class_<GridFieldProxy>& cls) {
  cls.def(py::init<>())
      // GridFieldProxy.geometry (0D_NOT_integer - Type of grid: xyz$, or rotationally_symmetric_rz$
      .def_property(
          "geometry", &GridFieldProxy::geometry, &GridFieldProxy::set_geometry)
      // GridFieldProxy.harmonic (0D_NOT_integer - Harmonic of fundamental for AC fields.
      .def_property(
          "harmonic", &GridFieldProxy::harmonic, &GridFieldProxy::set_harmonic)
      // GridFieldProxy.phi0_fieldmap (0D_NOT_real - Mode oscillates as: twopi * (f * t + phi0_fieldmap)
      .def_property(
          "phi0_fieldmap",
          &GridFieldProxy::phi0_fieldmap,
          &GridFieldProxy::set_phi0_fieldmap)
      // GridFieldProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &GridFieldProxy::field_scale,
          &GridFieldProxy::set_field_scale)
      // GridFieldProxy.field_type (0D_NOT_integer - or magnetic$ or electric$
      .def_property(
          "field_type",
          &GridFieldProxy::field_type,
          &GridFieldProxy::set_field_type)
      // GridFieldProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &GridFieldProxy::master_parameter,
          &GridFieldProxy::set_master_parameter)
      // GridFieldProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &GridFieldProxy::ele_anchor_pt,
          &GridFieldProxy::set_ele_anchor_pt)
      // GridFieldProxy.interpolation_order (0D_NOT_integer - Possibilities are 1 or 3.
      .def_property(
          "interpolation_order",
          &GridFieldProxy::interpolation_order,
          &GridFieldProxy::set_interpolation_order)
      // GridFieldProxy.dr (1D_NOT_real - Grid spacing.
      .def_property_readonly("dr", &GridFieldProxy::dr)
      // GridFieldProxy.r0 (1D_NOT_real - Field origin relative to ele_anchor_pt.
      .def_property_readonly("r0", &GridFieldProxy::r0)
      // GridFieldProxy.curved_ref_frame (0D_NOT_logical -
      .def_property(
          "curved_ref_frame",
          &GridFieldProxy::curved_ref_frame,
          &GridFieldProxy::set_curved_ref_frame)
      // GridFieldProxy.ptr (0D_PTR_type -
      .def_property("ptr", &GridFieldProxy::ptr, &GridFieldProxy::set_ptr)
      // GridFieldProxy.bi_coef (3D_NOT_type - Save computed coefs for faster tracking
      .def_property_readonly("bi_coef", &GridFieldProxy::bi_coef)
      // GridFieldProxy.tri_coef (3D_NOT_type - Save computed coefs for faster tracking
      .def_property_readonly("tri_coef", &GridFieldProxy::tri_coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return GridFieldProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const GridFieldProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldProxy& self) {
            return GridFieldProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const GridFieldProxy& self, py::dict& memo) {
            return GridFieldProxy(self);
          })

      ;

  bind_FTypeArrayND<GridFieldProxyArray1D>(m, "GridFieldStructArray1D");
  bind_FTypeAlloc1D<GridFieldProxyAlloc1D>(m, "GridFieldStructAlloc1D");
  // 2D GridFieldProxy arrays are not used in structs/routines
  // 3D GridFieldProxy arrays are not used in structs/routines
}

// =============================================================================
// high_energy_space_charge_struct
void init_high_energy_space_charge_struct(
    py::module& m,
    py::class_<HighEnergySpaceChargeProxy>& cls) {
  cls.def(py::init<>())
      // HighEnergySpaceChargeProxy.closed_orb (0D_NOT_type - beam orbit
      .def_property(
          "closed_orb",
          &HighEnergySpaceChargeProxy::closed_orb,
          &HighEnergySpaceChargeProxy::set_closed_orb)
      // HighEnergySpaceChargeProxy.kick_const (0D_NOT_real -
      .def_property(
          "kick_const",
          &HighEnergySpaceChargeProxy::kick_const,
          &HighEnergySpaceChargeProxy::set_kick_const)
      // HighEnergySpaceChargeProxy.sig_x (0D_NOT_real -
      .def_property(
          "sig_x",
          &HighEnergySpaceChargeProxy::sig_x,
          &HighEnergySpaceChargeProxy::set_sig_x)
      // HighEnergySpaceChargeProxy.sig_y (0D_NOT_real -
      .def_property(
          "sig_y",
          &HighEnergySpaceChargeProxy::sig_y,
          &HighEnergySpaceChargeProxy::set_sig_y)
      // HighEnergySpaceChargeProxy.phi (0D_NOT_real - Rotation angle to go from lab frame to rotated frame.
      .def_property(
          "phi",
          &HighEnergySpaceChargeProxy::phi,
          &HighEnergySpaceChargeProxy::set_phi)
      // HighEnergySpaceChargeProxy.sin_phi (0D_NOT_real -
      .def_property(
          "sin_phi",
          &HighEnergySpaceChargeProxy::sin_phi,
          &HighEnergySpaceChargeProxy::set_sin_phi)
      // HighEnergySpaceChargeProxy.cos_phi (0D_NOT_real -
      .def_property(
          "cos_phi",
          &HighEnergySpaceChargeProxy::cos_phi,
          &HighEnergySpaceChargeProxy::set_cos_phi)
      // HighEnergySpaceChargeProxy.sig_z (0D_NOT_real -
      .def_property(
          "sig_z",
          &HighEnergySpaceChargeProxy::sig_z,
          &HighEnergySpaceChargeProxy::set_sig_z)

      .def(
          "__repr__",
          [](const HighEnergySpaceChargeProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const HighEnergySpaceChargeProxy& self) {
            return HighEnergySpaceChargeProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const HighEnergySpaceChargeProxy& self, py::dict& memo) {
            return HighEnergySpaceChargeProxy(self);
          })

      ;

  // 1D HighEnergySpaceChargeProxy arrays are not used in structs/routines
  // 2D HighEnergySpaceChargeProxy arrays are not used in structs/routines
  // 3D HighEnergySpaceChargeProxy arrays are not used in structs/routines
}

// =============================================================================
// interval1_coef_struct
void init_interval1_coef_struct(
    py::module& m,
    py::class_<Interval1CoefProxy>& cls) {
  cls.def(py::init<>())
      // Interval1CoefProxy.c0 (0D_NOT_real -
      .def_property("c0", &Interval1CoefProxy::c0, &Interval1CoefProxy::set_c0)
      // Interval1CoefProxy.c1 (0D_NOT_real -
      .def_property("c1", &Interval1CoefProxy::c1, &Interval1CoefProxy::set_c1)
      // Interval1CoefProxy.n_exp (0D_NOT_real -
      .def_property(
          "n_exp", &Interval1CoefProxy::n_exp, &Interval1CoefProxy::set_n_exp)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return Interval1CoefProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const Interval1CoefProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Interval1CoefProxy& self) {
            return Interval1CoefProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Interval1CoefProxy& self, py::dict& memo) {
            return Interval1CoefProxy(self);
          })

      ;

  bind_FTypeArrayND<Interval1CoefProxyArray1D>(m, "Interval1CoefStructArray1D");
  bind_FTypeAlloc1D<Interval1CoefProxyAlloc1D>(m, "Interval1CoefStructAlloc1D");
  // 2D Interval1CoefProxy arrays are not used in structs/routines
  // 3D Interval1CoefProxy arrays are not used in structs/routines
}

// =============================================================================
// kv_beam_init_struct
void init_kv_beam_init_struct(py::module& m, py::class_<KvBeamInitProxy>& cls) {
  cls.def(py::init<>())
      // KvBeamInitProxy.part_per_phi (1D_NOT_integer - number of particles per angle variable.
      .def_property_readonly("part_per_phi", &KvBeamInitProxy::part_per_phi)
      // KvBeamInitProxy.n_I2 (0D_NOT_integer - number of I2
      .def_property("n_I2", &KvBeamInitProxy::n_I2, &KvBeamInitProxy::set_n_I2)
      // KvBeamInitProxy.A (0D_NOT_real - A = I1/e
      .def_property("A", &KvBeamInitProxy::A, &KvBeamInitProxy::set_A)

      .def(
          "__repr__",
          [](const KvBeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const KvBeamInitProxy& self) {
            return KvBeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const KvBeamInitProxy& self, py::dict& memo) {
            return KvBeamInitProxy(self);
          })

      ;

  // 1D KvBeamInitProxy arrays are not used in structs/routines
  // 2D KvBeamInitProxy arrays are not used in structs/routines
  // 3D KvBeamInitProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_loc_struct
void init_lat_ele_loc_struct(py::module& m, py::class_<LatEleLocProxy>& cls) {
  cls.def(py::init<>())
      // LatEleLocProxy.ix_ele (0D_NOT_integer -
      .def_property(
          "ix_ele", &LatEleLocProxy::ix_ele, &LatEleLocProxy::set_ix_ele)
      // LatEleLocProxy.ix_branch (0D_NOT_integer -
      .def_property(
          "ix_branch",
          &LatEleLocProxy::ix_branch,
          &LatEleLocProxy::set_ix_branch)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return LatEleLocProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleLocProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleLocProxy& self) {
            return LatEleLocProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleLocProxy& self, py::dict& memo) {
            return LatEleLocProxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleLocProxyArray1D>(m, "LatEleLocStructArray1D");
  bind_FTypeAlloc1D<LatEleLocProxyAlloc1D>(m, "LatEleLocStructAlloc1D");
  // 2D LatEleLocProxy arrays are not used in structs/routines
  // 3D LatEleLocProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order1_struct
void init_lat_ele_order1_struct(
    py::module& m,
    py::class_<LatEleOrder1Proxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrder1Proxy.ix_branch (0D_NOT_integer - Branch index
      .def_property(
          "ix_branch",
          &LatEleOrder1Proxy::ix_branch,
          &LatEleOrder1Proxy::set_ix_branch)
      // LatEleOrder1Proxy.ix_order (0D_NOT_integer - Order index. -1 -> Unique in lattice, 0 -> unique in branch.
      .def_property(
          "ix_order",
          &LatEleOrder1Proxy::ix_order,
          &LatEleOrder1Proxy::set_ix_order)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return LatEleOrder1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleOrder1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrder1Proxy& self) {
            return LatEleOrder1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrder1Proxy& self, py::dict& memo) {
            return LatEleOrder1Proxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleOrder1ProxyArray1D>(m, "LatEleOrder1StructArray1D");
  bind_FTypeAlloc1D<LatEleOrder1ProxyAlloc1D>(m, "LatEleOrder1StructAlloc1D");
  // 2D LatEleOrder1Proxy arrays are not used in structs/routines
  // 3D LatEleOrder1Proxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_array_struct
void init_lat_ele_order_array_struct(
    py::module& m,
    py::class_<LatEleOrderArrayProxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrderArrayProxy.ele (1D_ALLOC_type -
      .def_property_readonly("ele", &LatEleOrderArrayProxy::ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return LatEleOrderArrayProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleOrderArrayProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderArrayProxy& self) {
            return LatEleOrderArrayProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrderArrayProxy& self, py::dict& memo) {
            return LatEleOrderArrayProxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleOrderArrayProxyArray1D>(
      m, "LatEleOrderArrayStructArray1D");
  bind_FTypeAlloc1D<LatEleOrderArrayProxyAlloc1D>(
      m, "LatEleOrderArrayStructAlloc1D");
  // 2D LatEleOrderArrayProxy arrays are not used in structs/routines
  // 3D LatEleOrderArrayProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_struct
void init_lat_ele_order_struct(
    py::module& m,
    py::class_<LatEleOrderProxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrderProxy.branch (1D_ALLOC_type -
      .def_property_readonly("branch", &LatEleOrderProxy::branch)

      .def(
          "__repr__",
          [](const LatEleOrderProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderProxy& self) {
            return LatEleOrderProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrderProxy& self, py::dict& memo) {
            return LatEleOrderProxy(self);
          })

      ;

  // 1D LatEleOrderProxy arrays are not used in structs/routines
  // 2D LatEleOrderProxy arrays are not used in structs/routines
  // 3D LatEleOrderProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_param_struct
void init_lat_param_struct(py::module& m, py::class_<LatParamProxy>& cls) {
  cls.def(py::init<>())
      // LatParamProxy.n_part (0D_NOT_real - Particles/bunch (for BeamBeam elements).
      .def_property(
          "n_part", &LatParamProxy::n_part, &LatParamProxy::set_n_part)
      // LatParamProxy.total_length (0D_NOT_real - total_length of branch. Warning: branch may not start at s = 0.
      .def_property(
          "total_length",
          &LatParamProxy::total_length,
          &LatParamProxy::set_total_length)
      // LatParamProxy.unstable_factor (0D_NOT_real - If positive: Growth rate/turn if unstable in closed branches or |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
      .def_property(
          "unstable_factor",
          &LatParamProxy::unstable_factor,
          &LatParamProxy::set_unstable_factor)
      // LatParamProxy.t1_with_RF (2D_NOT_real - Full 1-turn matrix with RF on.
      .def_property_readonly("t1_with_RF", &LatParamProxy::t1_with_RF)
      // LatParamProxy.t1_no_RF (2D_NOT_real - Full 1-turn matrix with RF off.
      .def_property_readonly("t1_no_RF", &LatParamProxy::t1_no_RF)
      // LatParamProxy.spin_tune (0D_NOT_real - Closed orbit spin tune.
      .def_property(
          "spin_tune", &LatParamProxy::spin_tune, &LatParamProxy::set_spin_tune)
      // LatParamProxy.particle (0D_NOT_integer - Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is changed.
      .def_property(
          "particle", &LatParamProxy::particle, &LatParamProxy::set_particle)
      // LatParamProxy.default_tracking_species (0D_NOT_integer - Default particle type to use in tracking.
      .def_property(
          "default_tracking_species",
          &LatParamProxy::default_tracking_species,
          &LatParamProxy::set_default_tracking_species)
      // LatParamProxy.geometry (0D_NOT_integer - open$ or closed$
      .def_property(
          "geometry", &LatParamProxy::geometry, &LatParamProxy::set_geometry)
      // LatParamProxy.ixx (0D_NOT_integer - Integer for general use
      .def_property("ixx", &LatParamProxy::ixx, &LatParamProxy::set_ixx)
      // LatParamProxy.stable (0D_NOT_logical - is closed lat stable?
      .def_property(
          "stable", &LatParamProxy::stable, &LatParamProxy::set_stable)
      // LatParamProxy.live_branch (0D_NOT_logical - Should tracking be done on the branch?
      .def_property(
          "live_branch",
          &LatParamProxy::live_branch,
          &LatParamProxy::set_live_branch)
      // LatParamProxy.g1_integral (0D_NOT_real - Approximate |g| (bending strength) integral of branch.
      .def_property(
          "g1_integral",
          &LatParamProxy::g1_integral,
          &LatParamProxy::set_g1_integral)
      // LatParamProxy.g2_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property(
          "g2_integral",
          &LatParamProxy::g2_integral,
          &LatParamProxy::set_g2_integral)
      // LatParamProxy.g3_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property(
          "g3_integral",
          &LatParamProxy::g3_integral,
          &LatParamProxy::set_g3_integral)
      // LatParamProxy.bookkeeping_state (0D_NOT_type - Overall status for the branch.
      .def_property(
          "bookkeeping_state",
          &LatParamProxy::bookkeeping_state,
          &LatParamProxy::set_bookkeeping_state)
      // LatParamProxy.beam_init (0D_NOT_type - For beam initialization.
      .def_property(
          "beam_init", &LatParamProxy::beam_init, &LatParamProxy::set_beam_init)

      .def(
          "__repr__", [](const LatParamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatParamProxy& self) {
            return LatParamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatParamProxy& self, py::dict& memo) {
            return LatParamProxy(self);
          })

      ;

  // 1D LatParamProxy arrays are not used in structs/routines
  // 2D LatParamProxy arrays are not used in structs/routines
  // 3D LatParamProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_struct
void init_lat_struct(py::module& m, py::class_<LatProxy>& cls) {
  cls.def(py::init<>())
      // LatProxy.use_name (0D_NOT_character - Name of lat given by USE statement
      .def_property("use_name", &LatProxy::use_name, &LatProxy::set_use_name)
      // LatProxy.lattice (0D_NOT_character - Lattice
      .def_property("lattice", &LatProxy::lattice, &LatProxy::set_lattice)
      // LatProxy.machine (0D_NOT_character - Name of the machine the lattice is for ('LHC', etc).
      .def_property("machine", &LatProxy::machine, &LatProxy::set_machine)
      // LatProxy.input_file_name (0D_NOT_character - Name of the lattice input file
      .def_property(
          "input_file_name",
          &LatProxy::input_file_name,
          &LatProxy::set_input_file_name)
      // LatProxy.title (0D_NOT_character - General title
      .def_property("title", &LatProxy::title, &LatProxy::set_title)
      // LatProxy.print_str (1D_ALLOC_character - Saved print statements.
      .def_property_readonly("print_str", &LatProxy::print_str)
      // LatProxy.constant (1D_ALLOC_type - Constants defined in the lattice
      .def_property_readonly("constant", &LatProxy::constant)
      // LatProxy.a (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("a", &LatProxy::a, &LatProxy::set_a)
      // LatProxy.b (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("b", &LatProxy::b, &LatProxy::set_b)
      // LatProxy.z (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("z", &LatProxy::z, &LatProxy::set_z)
      // LatProxy.param (0D_PTR_type - Parameters
      .def_property("param", &LatProxy::param, &LatProxy::set_param)
      // LatProxy.lord_state (0D_NOT_type - lord bookkeeping status.
      .def_property(
          "lord_state", &LatProxy::lord_state, &LatProxy::set_lord_state)
      // LatProxy.ele_init (0D_NOT_type - For use by any program
      .def_property("ele_init", &LatProxy::ele_init, &LatProxy::set_ele_init)
      // LatProxy.ele (1D_PTR_type - Array of elements [=> branch(0)].
      .def_property_readonly("ele", &LatProxy::ele)
      // LatProxy.branch (1D_ALLOC_type - Branch(0:) array
      .def_property_readonly("branch", &LatProxy::branch)
      // LatProxy.control (1D_ALLOC_type - Control list
      .def_property_readonly("control", &LatProxy::control)
      // LatProxy.particle_start (0D_PTR_type - Starting particle_coords.
      .def_property(
          "particle_start",
          &LatProxy::particle_start,
          &LatProxy::set_particle_start)
      // LatProxy.beam_init (0D_NOT_type - Beam initialization.
      .def_property("beam_init", &LatProxy::beam_init, &LatProxy::set_beam_init)
      // LatProxy.pre_tracker (0D_NOT_type - For OPAL/IMPACT-T
      .def_property(
          "pre_tracker", &LatProxy::pre_tracker, &LatProxy::set_pre_tracker)
      // LatProxy.custom (1D_ALLOC_real - Custom attributes.
      .def_property_readonly("custom", &LatProxy::custom)
      // LatProxy.version (0D_NOT_integer - Version number
      .def_property("version", &LatProxy::version, &LatProxy::set_version)
      // LatProxy.n_ele_track (0D_PTR_integer - Number of lat elements to track through.
      .def_property(
          "n_ele_track", &LatProxy::n_ele_track, &LatProxy::set_n_ele_track)
      // LatProxy.n_ele_max (0D_PTR_integer - Index of last valid element in %ele(:) array
      .def_property("n_ele_max", &LatProxy::n_ele_max, &LatProxy::set_n_ele_max)
      // LatProxy.n_control_max (0D_NOT_integer - Last index used in control_array
      .def_property(
          "n_control_max",
          &LatProxy::n_control_max,
          &LatProxy::set_n_control_max)
      // LatProxy.n_ic_max (0D_NOT_integer - Last index used in ic_array
      .def_property("n_ic_max", &LatProxy::n_ic_max, &LatProxy::set_n_ic_max)
      // LatProxy.input_taylor_order (0D_NOT_integer - As set in the input file
      .def_property(
          "input_taylor_order",
          &LatProxy::input_taylor_order,
          &LatProxy::set_input_taylor_order)
      // LatProxy.ic (1D_ALLOC_integer - Index to %control(:) from slaves.
      .def_property_readonly("ic", &LatProxy::ic)
      // LatProxy.photon_type (0D_NOT_integer - Or coherent$. For X-ray simulations.
      .def_property(
          "photon_type", &LatProxy::photon_type, &LatProxy::set_photon_type)
      // LatProxy.creation_hash (0D_NOT_integer - Set by bmad_parser. creation_hash will vary if any of the lattice files are modified.
      .def_property(
          "creation_hash",
          &LatProxy::creation_hash,
          &LatProxy::set_creation_hash)
      // LatProxy.ramper_slave_bookkeeping (0D_NOT_integer -
      .def_property(
          "ramper_slave_bookkeeping",
          &LatProxy::ramper_slave_bookkeeping,
          &LatProxy::set_ramper_slave_bookkeeping)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return LatProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const LatProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatProxy& self) {
            return LatProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatProxy& self, py::dict& memo) { return LatProxy(self); })

      ;

  bind_FTypeArrayND<LatProxyArray1D>(m, "LatStructArray1D");
  bind_FTypeAlloc1D<LatProxyAlloc1D>(m, "LatStructAlloc1D");
  // 2D LatProxy arrays are not used in structs/routines
  // 3D LatProxy arrays are not used in structs/routines
}

// =============================================================================
// linac_normal_mode_struct
void init_linac_normal_mode_struct(
    py::module& m,
    py::class_<LinacNormalModeProxy>& cls) {
  cls.def(py::init<>())
      // LinacNormalModeProxy.i2_E4 (0D_NOT_real - Integral: g^2 * gamma^4
      .def_property(
          "i2_E4",
          &LinacNormalModeProxy::i2_E4,
          &LinacNormalModeProxy::set_i2_E4)
      // LinacNormalModeProxy.i3_E7 (0D_NOT_real - Integral: g^3 * gamma^7
      .def_property(
          "i3_E7",
          &LinacNormalModeProxy::i3_E7,
          &LinacNormalModeProxy::set_i3_E7)
      // LinacNormalModeProxy.i5a_E6 (0D_NOT_real - Integral: (g^3 * H_a) * gamma^6
      .def_property(
          "i5a_E6",
          &LinacNormalModeProxy::i5a_E6,
          &LinacNormalModeProxy::set_i5a_E6)
      // LinacNormalModeProxy.i5b_E6 (0D_NOT_real - Integral: (g^3 * H_b) * gamma^6
      .def_property(
          "i5b_E6",
          &LinacNormalModeProxy::i5b_E6,
          &LinacNormalModeProxy::set_i5b_E6)
      // LinacNormalModeProxy.sig_E1 (0D_NOT_real - Energy spread after 1 pass (eV)
      .def_property(
          "sig_E1",
          &LinacNormalModeProxy::sig_E1,
          &LinacNormalModeProxy::set_sig_E1)
      // LinacNormalModeProxy.a_emittance_end (0D_NOT_real - a mode emittance at end of linac
      .def_property(
          "a_emittance_end",
          &LinacNormalModeProxy::a_emittance_end,
          &LinacNormalModeProxy::set_a_emittance_end)
      // LinacNormalModeProxy.b_emittance_end (0D_NOT_real - b mode emittance at end of linac
      .def_property(
          "b_emittance_end",
          &LinacNormalModeProxy::b_emittance_end,
          &LinacNormalModeProxy::set_b_emittance_end)

      .def(
          "__repr__",
          [](const LinacNormalModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LinacNormalModeProxy& self) {
            return LinacNormalModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LinacNormalModeProxy& self, py::dict& memo) {
            return LinacNormalModeProxy(self);
          })

      ;

  // 1D LinacNormalModeProxy arrays are not used in structs/routines
  // 2D LinacNormalModeProxy arrays are not used in structs/routines
  // 3D LinacNormalModeProxy arrays are not used in structs/routines
}

// =============================================================================
// mode3_struct
void init_mode3_struct(py::module& m, py::class_<Mode3Proxy>& cls) {
  cls.def(py::init<>())
      // Mode3Proxy.v (2D_NOT_real -
      .def_property_readonly("v", &Mode3Proxy::v)
      // Mode3Proxy.a (0D_NOT_type -
      .def_property("a", &Mode3Proxy::a, &Mode3Proxy::set_a)
      // Mode3Proxy.b (0D_NOT_type -
      .def_property("b", &Mode3Proxy::b, &Mode3Proxy::set_b)
      // Mode3Proxy.c (0D_NOT_type -
      .def_property("c", &Mode3Proxy::c, &Mode3Proxy::set_c)
      // Mode3Proxy.x (0D_NOT_type -
      .def_property("x", &Mode3Proxy::x, &Mode3Proxy::set_x)
      // Mode3Proxy.y (0D_NOT_type -
      .def_property("y", &Mode3Proxy::y, &Mode3Proxy::set_y)

      .def("__repr__", [](const Mode3Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mode3Proxy& self) {
            return Mode3Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Mode3Proxy& self, py::dict& memo) {
            return Mode3Proxy(self);
          })

      ;

  // 1D Mode3Proxy arrays are not used in structs/routines
  // 2D Mode3Proxy arrays are not used in structs/routines
  // 3D Mode3Proxy arrays are not used in structs/routines
}

// =============================================================================
// mode_info_struct
void init_mode_info_struct(py::module& m, py::class_<ModeInfoProxy>& cls) {
  cls.def(py::init<>())
      // ModeInfoProxy.stable (0D_NOT_logical - Is the mode stable?
      .def_property(
          "stable", &ModeInfoProxy::stable, &ModeInfoProxy::set_stable)
      // ModeInfoProxy.tune (0D_NOT_real - 'fractional' tune in radians
      .def_property("tune", &ModeInfoProxy::tune, &ModeInfoProxy::set_tune)
      // ModeInfoProxy.emit (0D_NOT_real - Emittance (unnormalized).
      .def_property("emit", &ModeInfoProxy::emit, &ModeInfoProxy::set_emit)
      // ModeInfoProxy.chrom (0D_NOT_real - Chromaticity.
      .def_property("chrom", &ModeInfoProxy::chrom, &ModeInfoProxy::set_chrom)
      // ModeInfoProxy.sigma (0D_NOT_real - Beam size.
      .def_property("sigma", &ModeInfoProxy::sigma, &ModeInfoProxy::set_sigma)
      // ModeInfoProxy.sigmap (0D_NOT_real - Beam divergence.
      .def_property(
          "sigmap", &ModeInfoProxy::sigmap, &ModeInfoProxy::set_sigmap)

      .def(
          "__repr__", [](const ModeInfoProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ModeInfoProxy& self) {
            return ModeInfoProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ModeInfoProxy& self, py::dict& memo) {
            return ModeInfoProxy(self);
          })

      ;

  // 1D ModeInfoProxy arrays are not used in structs/routines
  // 2D ModeInfoProxy arrays are not used in structs/routines
  // 3D ModeInfoProxy arrays are not used in structs/routines
}

// =============================================================================
// normal_modes_struct
void init_normal_modes_struct(
    py::module& m,
    py::class_<NormalModesProxy>& cls) {
  cls.def(py::init<>())
      // NormalModesProxy.synch_int (1D_NOT_real - Synchrotron integrals I0, I1, I2, and I3
      .def_property_readonly("synch_int", &NormalModesProxy::synch_int)
      // NormalModesProxy.sigE_E (0D_NOT_real - SigmaE/E
      .def_property(
          "sigE_E", &NormalModesProxy::sigE_E, &NormalModesProxy::set_sigE_E)
      // NormalModesProxy.sig_z (0D_NOT_real - Sigma_Z
      .def_property(
          "sig_z", &NormalModesProxy::sig_z, &NormalModesProxy::set_sig_z)
      // NormalModesProxy.e_loss (0D_NOT_real - Energy loss / turn (eV)
      .def_property(
          "e_loss", &NormalModesProxy::e_loss, &NormalModesProxy::set_e_loss)
      // NormalModesProxy.rf_voltage (0D_NOT_real - Total rfcavity voltage (eV)
      .def_property(
          "rf_voltage",
          &NormalModesProxy::rf_voltage,
          &NormalModesProxy::set_rf_voltage)
      // NormalModesProxy.pz_aperture (0D_NOT_real - pz aperture limit. Used with Touschek calculations.
      .def_property(
          "pz_aperture",
          &NormalModesProxy::pz_aperture,
          &NormalModesProxy::set_pz_aperture)
      // NormalModesProxy.pz_average (0D_NOT_real - Average over branch due to damping.
      .def_property(
          "pz_average",
          &NormalModesProxy::pz_average,
          &NormalModesProxy::set_pz_average)
      // NormalModesProxy.momentum_compaction (0D_NOT_real -
      .def_property(
          "momentum_compaction",
          &NormalModesProxy::momentum_compaction,
          &NormalModesProxy::set_momentum_compaction)
      // NormalModesProxy.dpz_damp (0D_NOT_real - Change in pz without RF
      .def_property(
          "dpz_damp",
          &NormalModesProxy::dpz_damp,
          &NormalModesProxy::set_dpz_damp)
      // NormalModesProxy.a (0D_NOT_type -
      .def_property("a", &NormalModesProxy::a, &NormalModesProxy::set_a)
      // NormalModesProxy.b (0D_NOT_type -
      .def_property("b", &NormalModesProxy::b, &NormalModesProxy::set_b)
      // NormalModesProxy.z (0D_NOT_type -
      .def_property("z", &NormalModesProxy::z, &NormalModesProxy::set_z)
      // NormalModesProxy.lin (0D_NOT_type -
      .def_property("lin", &NormalModesProxy::lin, &NormalModesProxy::set_lin)

      .def(
          "__repr__",
          [](const NormalModesProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NormalModesProxy& self) {
            return NormalModesProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const NormalModesProxy& self, py::dict& memo) {
            return NormalModesProxy(self);
          })

      ;

  // 1D NormalModesProxy arrays are not used in structs/routines
  // 2D NormalModesProxy arrays are not used in structs/routines
  // 3D NormalModesProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_element_struct
void init_photon_element_struct(
    py::module& m,
    py::class_<PhotonElementProxy>& cls) {
  cls.def(py::init<>())
      // PhotonElementProxy.curvature (0D_NOT_type -
      .def_property(
          "curvature",
          &PhotonElementProxy::curvature,
          &PhotonElementProxy::set_curvature)
      // PhotonElementProxy.target (0D_NOT_type -
      .def_property(
          "target",
          &PhotonElementProxy::target,
          &PhotonElementProxy::set_target)
      // PhotonElementProxy.material (0D_NOT_type -
      .def_property(
          "material",
          &PhotonElementProxy::material,
          &PhotonElementProxy::set_material)
      // PhotonElementProxy.segmented (0D_NOT_type -
      .def_property(
          "segmented",
          &PhotonElementProxy::segmented,
          &PhotonElementProxy::set_segmented)
      // PhotonElementProxy.h_misalign (0D_NOT_type -
      .def_property(
          "h_misalign",
          &PhotonElementProxy::h_misalign,
          &PhotonElementProxy::set_h_misalign)
      // PhotonElementProxy.displacement (0D_NOT_type -
      .def_property(
          "displacement",
          &PhotonElementProxy::displacement,
          &PhotonElementProxy::set_displacement)
      // PhotonElementProxy.pixel (0D_NOT_type -
      .def_property(
          "pixel", &PhotonElementProxy::pixel, &PhotonElementProxy::set_pixel)
      // PhotonElementProxy.reflectivity_table_type (0D_NOT_integer -
      .def_property(
          "reflectivity_table_type",
          &PhotonElementProxy::reflectivity_table_type,
          &PhotonElementProxy::set_reflectivity_table_type)
      // PhotonElementProxy.reflectivity_table_sigma (0D_NOT_type - If polarization is ignored use sigma table.
      .def_property(
          "reflectivity_table_sigma",
          &PhotonElementProxy::reflectivity_table_sigma,
          &PhotonElementProxy::set_reflectivity_table_sigma)
      // PhotonElementProxy.reflectivity_table_pi (0D_NOT_type -
      .def_property(
          "reflectivity_table_pi",
          &PhotonElementProxy::reflectivity_table_pi,
          &PhotonElementProxy::set_reflectivity_table_pi)
      // PhotonElementProxy.init_energy_prob (1D_ALLOC_type - Initial energy probability density
      .def_property_readonly(
          "init_energy_prob", &PhotonElementProxy::init_energy_prob)
      // PhotonElementProxy.integrated_init_energy_prob (1D_ALLOC_real -
      .def_property_readonly(
          "integrated_init_energy_prob",
          &PhotonElementProxy::integrated_init_energy_prob)

      .def(
          "__repr__",
          [](const PhotonElementProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonElementProxy& self) {
            return PhotonElementProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonElementProxy& self, py::dict& memo) {
            return PhotonElementProxy(self);
          })

      ;

  // 1D PhotonElementProxy arrays are not used in structs/routines
  // 2D PhotonElementProxy arrays are not used in structs/routines
  // 3D PhotonElementProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_material_struct
void init_photon_material_struct(
    py::module& m,
    py::class_<PhotonMaterialProxy>& cls) {
  cls.def(py::init<>())
      // PhotonMaterialProxy.f0_m1 (0D_NOT_complex - For multilayer_mirror only.
      .def_property(
          "f0_m1", &PhotonMaterialProxy::f0_m1, &PhotonMaterialProxy::set_f0_m1)
      // PhotonMaterialProxy.f0_m2 (0D_NOT_complex - For multilayer_mirror only.
      .def_property(
          "f0_m2", &PhotonMaterialProxy::f0_m2, &PhotonMaterialProxy::set_f0_m2)
      // PhotonMaterialProxy.f_0 (0D_NOT_complex -
      .def_property(
          "f_0", &PhotonMaterialProxy::f_0, &PhotonMaterialProxy::set_f_0)
      // PhotonMaterialProxy.f_h (0D_NOT_complex - Structure factor for H direction.
      .def_property(
          "f_h", &PhotonMaterialProxy::f_h, &PhotonMaterialProxy::set_f_h)
      // PhotonMaterialProxy.f_hbar (0D_NOT_complex - Structure factor for -H direction.
      .def_property(
          "f_hbar",
          &PhotonMaterialProxy::f_hbar,
          &PhotonMaterialProxy::set_f_hbar)
      // PhotonMaterialProxy.f_hkl (0D_NOT_complex - = sqrt(f_h * f_hbar)
      .def_property(
          "f_hkl", &PhotonMaterialProxy::f_hkl, &PhotonMaterialProxy::set_f_hkl)
      // PhotonMaterialProxy.h_norm (1D_NOT_real - Normalized H vector for crystals.
      .def_property_readonly("h_norm", &PhotonMaterialProxy::h_norm)
      // PhotonMaterialProxy.l_ref (1D_NOT_real - Crystal reference orbit displacement vector in element coords.
      .def_property_readonly("l_ref", &PhotonMaterialProxy::l_ref)

      .def(
          "__repr__",
          [](const PhotonMaterialProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonMaterialProxy& self) {
            return PhotonMaterialProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonMaterialProxy& self, py::dict& memo) {
            return PhotonMaterialProxy(self);
          })

      ;

  // 1D PhotonMaterialProxy arrays are not used in structs/routines
  // 2D PhotonMaterialProxy arrays are not used in structs/routines
  // 3D PhotonMaterialProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_surface_struct
void init_photon_reflect_surface_struct(
    py::module& m,
    py::class_<PhotonReflectSurfaceProxy>& cls) {
  cls.def(py::init<>())
      // PhotonReflectSurfaceProxy.name (0D_NOT_character -
      .def_property(
          "name",
          &PhotonReflectSurfaceProxy::name,
          &PhotonReflectSurfaceProxy::set_name)
      // PhotonReflectSurfaceProxy.description (0D_NOT_character - Descriptive name
      .def_property(
          "description",
          &PhotonReflectSurfaceProxy::description,
          &PhotonReflectSurfaceProxy::set_description)
      // PhotonReflectSurfaceProxy.reflectivity_file (0D_NOT_character -
      .def_property(
          "reflectivity_file",
          &PhotonReflectSurfaceProxy::reflectivity_file,
          &PhotonReflectSurfaceProxy::set_reflectivity_file)
      // PhotonReflectSurfaceProxy.table (1D_ALLOC_type -
      .def_property_readonly("table", &PhotonReflectSurfaceProxy::table)
      // PhotonReflectSurfaceProxy.surface_roughness_rms (0D_NOT_real - sigma in Dugan's notation
      .def_property(
          "surface_roughness_rms",
          &PhotonReflectSurfaceProxy::surface_roughness_rms,
          &PhotonReflectSurfaceProxy::set_surface_roughness_rms)
      // PhotonReflectSurfaceProxy.roughness_correlation_len (0D_NOT_real - T in Dugan's notation
      .def_property(
          "roughness_correlation_len",
          &PhotonReflectSurfaceProxy::roughness_correlation_len,
          &PhotonReflectSurfaceProxy::set_roughness_correlation_len)
      // PhotonReflectSurfaceProxy.ix_surface (0D_NOT_integer -
      .def_property(
          "ix_surface",
          &PhotonReflectSurfaceProxy::ix_surface,
          &PhotonReflectSurfaceProxy::set_ix_surface)

      .def(
          "__repr__",
          [](const PhotonReflectSurfaceProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectSurfaceProxy& self) {
            return PhotonReflectSurfaceProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonReflectSurfaceProxy& self, py::dict& memo) {
            return PhotonReflectSurfaceProxy(self);
          })

      ;

  // 1D PhotonReflectSurfaceProxy arrays are not used in structs/routines
  // 2D PhotonReflectSurfaceProxy arrays are not used in structs/routines
  // 3D PhotonReflectSurfaceProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_table_struct
void init_photon_reflect_table_struct(
    py::module& m,
    py::class_<PhotonReflectTableProxy>& cls) {
  cls.def(py::init<>())
      // PhotonReflectTableProxy.angle (1D_ALLOC_real - Vector of angle values for %p_reflect
      .def_property_readonly("angle", &PhotonReflectTableProxy::angle)
      // PhotonReflectTableProxy.energy (1D_ALLOC_real - Vector of energy values for %p_reflect
      .def_property_readonly("energy", &PhotonReflectTableProxy::energy)
      // PhotonReflectTableProxy.int1 (1D_ALLOC_type -
      .def_property_readonly("int1", &PhotonReflectTableProxy::int1)
      // PhotonReflectTableProxy.p_reflect (2D_ALLOC_real - (angle, ev) probability. Log used for smooth surface reflection
      .def_property_readonly("p_reflect", &PhotonReflectTableProxy::p_reflect)
      // PhotonReflectTableProxy.max_energy (0D_NOT_real - maximum energy for this table
      .def_property(
          "max_energy",
          &PhotonReflectTableProxy::max_energy,
          &PhotonReflectTableProxy::set_max_energy)
      // PhotonReflectTableProxy.p_reflect_scratch (1D_ALLOC_real - Scratch space
      .def_property_readonly(
          "p_reflect_scratch", &PhotonReflectTableProxy::p_reflect_scratch)
      // PhotonReflectTableProxy.bragg_angle (1D_ALLOC_real - Bragg angle at energy values.
      .def_property_readonly(
          "bragg_angle", &PhotonReflectTableProxy::bragg_angle)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return PhotonReflectTableProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const PhotonReflectTableProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectTableProxy& self) {
            return PhotonReflectTableProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonReflectTableProxy& self, py::dict& memo) {
            return PhotonReflectTableProxy(self);
          })

      ;

  bind_FTypeArrayND<PhotonReflectTableProxyArray1D>(
      m, "PhotonReflectTableStructArray1D");
  bind_FTypeAlloc1D<PhotonReflectTableProxyAlloc1D>(
      m, "PhotonReflectTableStructAlloc1D");
  // 2D PhotonReflectTableProxy arrays are not used in structs/routines
  // 3D PhotonReflectTableProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_target_struct
void init_photon_target_struct(
    py::module& m,
    py::class_<PhotonTargetProxy>& cls) {
  cls.def(py::init<>())
      // PhotonTargetProxy.type (0D_NOT_integer - or rectangular$
      .def_property(
          "type", &PhotonTargetProxy::type, &PhotonTargetProxy::set_type)
      // PhotonTargetProxy.n_corner (0D_NOT_integer -
      .def_property(
          "n_corner",
          &PhotonTargetProxy::n_corner,
          &PhotonTargetProxy::set_n_corner)
      // PhotonTargetProxy.ele_loc (0D_NOT_type -
      .def_property(
          "ele_loc",
          &PhotonTargetProxy::ele_loc,
          &PhotonTargetProxy::set_ele_loc)
      // PhotonTargetProxy.corner (1D_NOT_type -
      .def_property_readonly("corner", &PhotonTargetProxy::corner)
      // PhotonTargetProxy.center (0D_NOT_type -
      .def_property(
          "center", &PhotonTargetProxy::center, &PhotonTargetProxy::set_center)

      .def(
          "__repr__",
          [](const PhotonTargetProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonTargetProxy& self) {
            return PhotonTargetProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonTargetProxy& self, py::dict& memo) {
            return PhotonTargetProxy(self);
          })

      ;

  // 1D PhotonTargetProxy arrays are not used in structs/routines
  // 2D PhotonTargetProxy arrays are not used in structs/routines
  // 3D PhotonTargetProxy arrays are not used in structs/routines
}

// =============================================================================
// pixel_detec_struct
void init_pixel_detec_struct(py::module& m, py::class_<PixelDetecProxy>& cls) {
  cls.def(py::init<>())
      // PixelDetecProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &PixelDetecProxy::dr)
      // PixelDetecProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &PixelDetecProxy::r0)
      // PixelDetecProxy.n_track_tot (0D_NOT_integer8 - How many photons were launched from source element.
      .def_property(
          "n_track_tot",
          &PixelDetecProxy::n_track_tot,
          &PixelDetecProxy::set_n_track_tot)
      // PixelDetecProxy.n_hit_detec (0D_NOT_integer8 - How many photons hit the detector.
      .def_property(
          "n_hit_detec",
          &PixelDetecProxy::n_hit_detec,
          &PixelDetecProxy::set_n_hit_detec)
      // PixelDetecProxy.n_hit_pixel (0D_NOT_integer8 - How many photons hit the pixel grid of the detector.
      .def_property(
          "n_hit_pixel",
          &PixelDetecProxy::n_hit_pixel,
          &PixelDetecProxy::set_n_hit_pixel)
      // PixelDetecProxy.pt (2D_ALLOC_type - Grid of pixels
      .def_property_readonly("pt", &PixelDetecProxy::pt)

      .def(
          "__repr__",
          [](const PixelDetecProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelDetecProxy& self) {
            return PixelDetecProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PixelDetecProxy& self, py::dict& memo) {
            return PixelDetecProxy(self);
          })

      ;

  // 1D PixelDetecProxy arrays are not used in structs/routines
  // 2D PixelDetecProxy arrays are not used in structs/routines
  // 3D PixelDetecProxy arrays are not used in structs/routines
}

// =============================================================================
// pixel_pt_struct
void init_pixel_pt_struct(py::module& m, py::class_<PixelPtProxy>& cls) {
  cls.def(py::init<>())
      // PixelPtProxy.n_photon (0D_NOT_integer8 -
      .def_property(
          "n_photon", &PixelPtProxy::n_photon, &PixelPtProxy::set_n_photon)
      // PixelPtProxy.E_x (0D_NOT_complex -
      .def_property("E_x", &PixelPtProxy::E_x, &PixelPtProxy::set_E_x)
      // PixelPtProxy.E_y (0D_NOT_complex -
      .def_property("E_y", &PixelPtProxy::E_y, &PixelPtProxy::set_E_y)
      // PixelPtProxy.intensity_x (0D_NOT_real -
      .def_property(
          "intensity_x",
          &PixelPtProxy::intensity_x,
          &PixelPtProxy::set_intensity_x)
      // PixelPtProxy.intensity_y (0D_NOT_real -
      .def_property(
          "intensity_y",
          &PixelPtProxy::intensity_y,
          &PixelPtProxy::set_intensity_y)
      // PixelPtProxy.intensity (0D_NOT_real -
      .def_property(
          "intensity", &PixelPtProxy::intensity, &PixelPtProxy::set_intensity)
      // PixelPtProxy.orbit (1D_NOT_real - x, Vx/c, y, Vy/c, dummy, E - E_ref.
      .def_property_readonly("orbit", &PixelPtProxy::orbit)
      // PixelPtProxy.orbit_rms (1D_NOT_real - RMS statistics.
      .def_property_readonly("orbit_rms", &PixelPtProxy::orbit_rms)
      // PixelPtProxy.init_orbit (1D_NOT_real - Initial orbit at start of lattice statistics.
      .def_property_readonly("init_orbit", &PixelPtProxy::init_orbit)
      // PixelPtProxy.init_orbit_rms (1D_NOT_real - Initial orbit at start of lattice RMS statistics.
      .def_property_readonly("init_orbit_rms", &PixelPtProxy::init_orbit_rms)

      .def("__repr__", [](const PixelPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelPtProxy& self) {
            return PixelPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PixelPtProxy& self, py::dict& memo) {
            return PixelPtProxy(self);
          })

      ;

  // 1D PixelPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<PixelPtProxyArray2D>(m, "PixelPtStructArray2D");
  // 3D PixelPtProxy arrays are not used in structs/routines
}

// =============================================================================
// pre_tracker_struct
void init_pre_tracker_struct(py::module& m, py::class_<PreTrackerProxy>& cls) {
  cls.def(py::init<>())
      // PreTrackerProxy.who (0D_NOT_integer - Can be opal$, or impactt$
      .def_property("who", &PreTrackerProxy::who, &PreTrackerProxy::set_who)
      // PreTrackerProxy.ix_ele_start (0D_NOT_integer -
      .def_property(
          "ix_ele_start",
          &PreTrackerProxy::ix_ele_start,
          &PreTrackerProxy::set_ix_ele_start)
      // PreTrackerProxy.ix_ele_end (0D_NOT_integer -
      .def_property(
          "ix_ele_end",
          &PreTrackerProxy::ix_ele_end,
          &PreTrackerProxy::set_ix_ele_end)
      // PreTrackerProxy.input_file (0D_NOT_character -
      .def_property(
          "input_file",
          &PreTrackerProxy::input_file,
          &PreTrackerProxy::set_input_file)

      .def(
          "__repr__",
          [](const PreTrackerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PreTrackerProxy& self) {
            return PreTrackerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PreTrackerProxy& self, py::dict& memo) {
            return PreTrackerProxy(self);
          })

      ;

  // 1D PreTrackerProxy arrays are not used in structs/routines
  // 2D PreTrackerProxy arrays are not used in structs/routines
  // 3D PreTrackerProxy arrays are not used in structs/routines
}

// =============================================================================
// ptc_normal_form_struct
void init_ptc_normal_form_struct(
    py::module& m,
    py::class_<PtcNormalFormProxy>& cls) {
  cls.def(py::init<>())
      // PtcNormalFormProxy.ele_origin (0D_PTR_type - Element at which the on-turn map was created.
      .def_property(
          "ele_origin",
          &PtcNormalFormProxy::ele_origin,
          &PtcNormalFormProxy::set_ele_origin)
      // PtcNormalFormProxy.orb0 (1D_NOT_real - Closed orbit at element.
      .def_property_readonly("orb0", &PtcNormalFormProxy::orb0)
      // PtcNormalFormProxy.valid_map (0D_NOT_logical -
      .def_property(
          "valid_map",
          &PtcNormalFormProxy::valid_map,
          &PtcNormalFormProxy::set_valid_map)

      .def(
          "__repr__",
          [](const PtcNormalFormProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PtcNormalFormProxy& self) {
            return PtcNormalFormProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PtcNormalFormProxy& self, py::dict& memo) {
            return PtcNormalFormProxy(self);
          })

      ;

  // 1D PtcNormalFormProxy arrays are not used in structs/routines
  // 2D PtcNormalFormProxy arrays are not used in structs/routines
  // 3D PtcNormalFormProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_int1_struct
void init_rad_int1_struct(py::module& m, py::class_<RadInt1Proxy>& cls) {
  cls.def(py::init<>())
      // RadInt1Proxy.i0 (0D_NOT_real -
      .def_property("i0", &RadInt1Proxy::i0, &RadInt1Proxy::set_i0)
      // RadInt1Proxy.i1 (0D_NOT_real -
      .def_property("i1", &RadInt1Proxy::i1, &RadInt1Proxy::set_i1)
      // RadInt1Proxy.i2 (0D_NOT_real -
      .def_property("i2", &RadInt1Proxy::i2, &RadInt1Proxy::set_i2)
      // RadInt1Proxy.i3 (0D_NOT_real -
      .def_property("i3", &RadInt1Proxy::i3, &RadInt1Proxy::set_i3)
      // RadInt1Proxy.i4a (0D_NOT_real -
      .def_property("i4a", &RadInt1Proxy::i4a, &RadInt1Proxy::set_i4a)
      // RadInt1Proxy.i4b (0D_NOT_real -
      .def_property("i4b", &RadInt1Proxy::i4b, &RadInt1Proxy::set_i4b)
      // RadInt1Proxy.i4z (0D_NOT_real -
      .def_property("i4z", &RadInt1Proxy::i4z, &RadInt1Proxy::set_i4z)
      // RadInt1Proxy.i5a (0D_NOT_real -
      .def_property("i5a", &RadInt1Proxy::i5a, &RadInt1Proxy::set_i5a)
      // RadInt1Proxy.i5b (0D_NOT_real -
      .def_property("i5b", &RadInt1Proxy::i5b, &RadInt1Proxy::set_i5b)
      // RadInt1Proxy.i6b (0D_NOT_real -
      .def_property("i6b", &RadInt1Proxy::i6b, &RadInt1Proxy::set_i6b)
      // RadInt1Proxy.lin_i2_E4 (0D_NOT_real -
      .def_property(
          "lin_i2_E4", &RadInt1Proxy::lin_i2_E4, &RadInt1Proxy::set_lin_i2_E4)
      // RadInt1Proxy.lin_i3_E7 (0D_NOT_real -
      .def_property(
          "lin_i3_E7", &RadInt1Proxy::lin_i3_E7, &RadInt1Proxy::set_lin_i3_E7)
      // RadInt1Proxy.lin_i5a_E6 (0D_NOT_real -
      .def_property(
          "lin_i5a_E6",
          &RadInt1Proxy::lin_i5a_E6,
          &RadInt1Proxy::set_lin_i5a_E6)
      // RadInt1Proxy.lin_i5b_E6 (0D_NOT_real -
      .def_property(
          "lin_i5b_E6",
          &RadInt1Proxy::lin_i5b_E6,
          &RadInt1Proxy::set_lin_i5b_E6)
      // RadInt1Proxy.lin_norm_emit_a (0D_NOT_real - Running sum
      .def_property(
          "lin_norm_emit_a",
          &RadInt1Proxy::lin_norm_emit_a,
          &RadInt1Proxy::set_lin_norm_emit_a)
      // RadInt1Proxy.lin_norm_emit_b (0D_NOT_real - Running sum
      .def_property(
          "lin_norm_emit_b",
          &RadInt1Proxy::lin_norm_emit_b,
          &RadInt1Proxy::set_lin_norm_emit_b)
      // RadInt1Proxy.lin_sig_E (0D_NOT_real - Running sum
      .def_property(
          "lin_sig_E", &RadInt1Proxy::lin_sig_E, &RadInt1Proxy::set_lin_sig_E)
      // RadInt1Proxy.n_steps (0D_NOT_real - number of qromb steps needed
      .def_property(
          "n_steps", &RadInt1Proxy::n_steps, &RadInt1Proxy::set_n_steps)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return RadInt1ProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const RadInt1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadInt1Proxy& self) {
            return RadInt1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadInt1Proxy& self, py::dict& memo) {
            return RadInt1Proxy(self);
          })

      ;

  bind_FTypeArrayND<RadInt1ProxyArray1D>(m, "RadInt1StructArray1D");
  bind_FTypeAlloc1D<RadInt1ProxyAlloc1D>(m, "RadInt1StructAlloc1D");
  // 2D RadInt1Proxy arrays are not used in structs/routines
  // 3D RadInt1Proxy arrays are not used in structs/routines
}

// =============================================================================
// rad_int_all_ele_struct
void init_rad_int_all_ele_struct(
    py::module& m,
    py::class_<RadIntAllEleProxy>& cls) {
  cls.def(py::init<>())
      // RadIntAllEleProxy.branch (1D_ALLOC_type - Array is indexed from 0
      .def_property_readonly("branch", &RadIntAllEleProxy::branch)

      .def(
          "__repr__",
          [](const RadIntAllEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntAllEleProxy& self) {
            return RadIntAllEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadIntAllEleProxy& self, py::dict& memo) {
            return RadIntAllEleProxy(self);
          })

      ;

  // 1D RadIntAllEleProxy arrays are not used in structs/routines
  // 2D RadIntAllEleProxy arrays are not used in structs/routines
  // 3D RadIntAllEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_int_branch_struct
void init_rad_int_branch_struct(
    py::module& m,
    py::class_<RadIntBranchProxy>& cls) {
  cls.def(py::init<>())
      // RadIntBranchProxy.ele (1D_ALLOC_type - Array is indexed from 0
      .def_property_readonly("ele", &RadIntBranchProxy::ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return RadIntBranchProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RadIntBranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntBranchProxy& self) {
            return RadIntBranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadIntBranchProxy& self, py::dict& memo) {
            return RadIntBranchProxy(self);
          })

      ;

  bind_FTypeArrayND<RadIntBranchProxyArray1D>(m, "RadIntBranchStructArray1D");
  bind_FTypeAlloc1D<RadIntBranchProxyAlloc1D>(m, "RadIntBranchStructAlloc1D");
  // 2D RadIntBranchProxy arrays are not used in structs/routines
  // 3D RadIntBranchProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_map_ele_struct
void init_rad_map_ele_struct(py::module& m, py::class_<RadMapEleProxy>& cls) {
  cls.def(py::init<>())
      // RadMapEleProxy.rm0 (0D_NOT_type - Upstream half and downstream half matrices for an element.
      .def_property("rm0", &RadMapEleProxy::rm0, &RadMapEleProxy::set_rm0)
      // RadMapEleProxy.rm1 (0D_NOT_type - Upstream half and downstream half matrices for an element.
      .def_property("rm1", &RadMapEleProxy::rm1, &RadMapEleProxy::set_rm1)
      // RadMapEleProxy.stale (0D_NOT_logical -
      .def_property("stale", &RadMapEleProxy::stale, &RadMapEleProxy::set_stale)

      .def(
          "__repr__",
          [](const RadMapEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapEleProxy& self) {
            return RadMapEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadMapEleProxy& self, py::dict& memo) {
            return RadMapEleProxy(self);
          })

      ;

  // 1D RadMapEleProxy arrays are not used in structs/routines
  // 2D RadMapEleProxy arrays are not used in structs/routines
  // 3D RadMapEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_map_struct
void init_rad_map_struct(py::module& m, py::class_<RadMapProxy>& cls) {
  cls.def(py::init<>())
      // RadMapProxy.ref_orb (1D_NOT_real - Reference point around which damp_mat is calculated.
      .def_property_readonly("ref_orb", &RadMapProxy::ref_orb)
      // RadMapProxy.damp_dmat (2D_NOT_real - damp_correction = xfer_mat_with_damping - xfer_mat_without_damping.
      .def_property_readonly("damp_dmat", &RadMapProxy::damp_dmat)
      // RadMapProxy.xfer_damp_vec (1D_NOT_real - Transfer map with damping 0th order vector.
      .def_property_readonly("xfer_damp_vec", &RadMapProxy::xfer_damp_vec)
      // RadMapProxy.xfer_damp_mat (2D_NOT_real - 1st order matrix: xfer_no_damp_mat + xfer_damp_correction.
      .def_property_readonly("xfer_damp_mat", &RadMapProxy::xfer_damp_mat)
      // RadMapProxy.stoc_mat (2D_NOT_real - Stochastic variance or 'kick' (Cholesky decomposed) matrix.
      .def_property_readonly("stoc_mat", &RadMapProxy::stoc_mat)

      .def("__repr__", [](const RadMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapProxy& self) {
            return RadMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadMapProxy& self, py::dict& memo) {
            return RadMapProxy(self);
          })

      ;

  // 1D RadMapProxy arrays are not used in structs/routines
  // 2D RadMapProxy arrays are not used in structs/routines
  // 3D RadMapProxy arrays are not used in structs/routines
}

// =============================================================================
// ramper_lord_struct
void init_ramper_lord_struct(py::module& m, py::class_<RamperLordProxy>& cls) {
  cls.def(py::init<>())
      // RamperLordProxy.ix_ele (0D_NOT_integer - Lord index
      .def_property(
          "ix_ele", &RamperLordProxy::ix_ele, &RamperLordProxy::set_ix_ele)
      // RamperLordProxy.ix_con (0D_NOT_integer - Index in lord%control%ramp(:) array
      .def_property(
          "ix_con", &RamperLordProxy::ix_con, &RamperLordProxy::set_ix_con)
      // RamperLordProxy.attrib_ptr (0D_PTR_real - Pointer to attribute in this element.
      .def_property(
          "attrib_ptr",
          &RamperLordProxy::attrib_ptr,
          &RamperLordProxy::set_attrib_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return RamperLordProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RamperLordProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RamperLordProxy& self) {
            return RamperLordProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RamperLordProxy& self, py::dict& memo) {
            return RamperLordProxy(self);
          })

      ;

  bind_FTypeArrayND<RamperLordProxyArray1D>(m, "RamperLordStructArray1D");
  bind_FTypeAlloc1D<RamperLordProxyAlloc1D>(m, "RamperLordStructAlloc1D");
  // 2D RamperLordProxy arrays are not used in structs/routines
  // 3D RamperLordProxy arrays are not used in structs/routines
}

// =============================================================================
// resonance_h_struct
void init_resonance_h_struct(py::module& m, py::class_<ResonanceHProxy>& cls) {
  cls.def(py::init<>())
      // ResonanceHProxy.id (0D_NOT_character - 6 digit ID. EG: '003100'
      .def_property("id", &ResonanceHProxy::id, &ResonanceHProxy::set_id)
      // ResonanceHProxy.c_val (0D_NOT_complex - Resonance value
      .def_property(
          "c_val", &ResonanceHProxy::c_val, &ResonanceHProxy::set_c_val)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ResonanceHProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ResonanceHProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ResonanceHProxy& self) {
            return ResonanceHProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ResonanceHProxy& self, py::dict& memo) {
            return ResonanceHProxy(self);
          })

      ;

  bind_FTypeArrayND<ResonanceHProxyArray1D>(m, "ResonanceHStructArray1D");
  bind_FTypeAlloc1D<ResonanceHProxyAlloc1D>(m, "ResonanceHStructAlloc1D");
  // 2D ResonanceHProxy arrays are not used in structs/routines
  // 3D ResonanceHProxy arrays are not used in structs/routines
}

// =============================================================================
// rf_ele_struct
void init_rf_ele_struct(py::module& m, py::class_<RfEleProxy>& cls) {
  cls.def(py::init<>())
      // RfEleProxy.steps (1D_ALLOC_type - Energy stair step array indexed from zero.
      .def_property_readonly("steps", &RfEleProxy::steps)
      // RfEleProxy.ds_step (0D_NOT_real - length of a stair step.
      .def_property("ds_step", &RfEleProxy::ds_step, &RfEleProxy::set_ds_step)

      .def("__repr__", [](const RfEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfEleProxy& self) {
            return RfEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RfEleProxy& self, py::dict& memo) {
            return RfEleProxy(self);
          })

      ;

  // 1D RfEleProxy arrays are not used in structs/routines
  // 2D RfEleProxy arrays are not used in structs/routines
  // 3D RfEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rf_stair_step_struct
void init_rf_stair_step_struct(
    py::module& m,
    py::class_<RfStairStepProxy>& cls) {
  cls.def(py::init<>())
      // RfStairStepProxy.E_tot0 (0D_NOT_real - Reference energy in the drift region (before the kick point).
      .def_property(
          "E_tot0", &RfStairStepProxy::E_tot0, &RfStairStepProxy::set_E_tot0)
      // RfStairStepProxy.E_tot1 (0D_NOT_real - Reference energy after the kick point.
      .def_property(
          "E_tot1", &RfStairStepProxy::E_tot1, &RfStairStepProxy::set_E_tot1)
      // RfStairStepProxy.p0c (0D_NOT_real - Reference momentum in the drift region (before the kick point).
      .def_property("p0c", &RfStairStepProxy::p0c, &RfStairStepProxy::set_p0c)
      // RfStairStepProxy.p1c (0D_NOT_real - Reference momentum after the kick point.
      .def_property("p1c", &RfStairStepProxy::p1c, &RfStairStepProxy::set_p1c)
      // RfStairStepProxy.scale (0D_NOT_real - Scale for multipole kick at the kick point. Sum over all steps will be 1.
      .def_property(
          "scale", &RfStairStepProxy::scale, &RfStairStepProxy::set_scale)
      // RfStairStepProxy.time (0D_NOT_real - Reference particle time at the kick point with respect to beginning of element.
      .def_property(
          "time", &RfStairStepProxy::time, &RfStairStepProxy::set_time)
      // RfStairStepProxy.s0 (0D_NOT_real - S-position at beginning of drift region relative to the beginning of the element.
      .def_property("s0", &RfStairStepProxy::s0, &RfStairStepProxy::set_s0)
      // RfStairStepProxy.s (0D_NOT_real - S-position at the kick point relative to the beginning of the element.
      .def_property("s", &RfStairStepProxy::s, &RfStairStepProxy::set_s)
      // RfStairStepProxy.ix_step (0D_NOT_integer - Step index in ele%rf%steps(:) array
      .def_property(
          "ix_step", &RfStairStepProxy::ix_step, &RfStairStepProxy::set_ix_step)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return RfStairStepProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RfStairStepProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfStairStepProxy& self) {
            return RfStairStepProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RfStairStepProxy& self, py::dict& memo) {
            return RfStairStepProxy(self);
          })

      ;

  bind_FTypeArrayND<RfStairStepProxyArray1D>(m, "RfStairStepStructArray1D");
  bind_FTypeAlloc1D<RfStairStepProxyAlloc1D>(m, "RfStairStepStructAlloc1D");
  // 2D RfStairStepProxy arrays are not used in structs/routines
  // 3D RfStairStepProxy arrays are not used in structs/routines
}

// =============================================================================
// space_charge_common_struct
void init_space_charge_common_struct(
    py::module& m,
    py::class_<SpaceChargeCommonProxy>& cls) {
  cls.def(py::init<>())
      // SpaceChargeCommonProxy.ds_track_step (0D_NOT_real - CSR tracking step size
      .def_property(
          "ds_track_step",
          &SpaceChargeCommonProxy::ds_track_step,
          &SpaceChargeCommonProxy::set_ds_track_step)
      // SpaceChargeCommonProxy.dt_track_step (0D_NOT_real - Time Runge kutta initial step.
      .def_property(
          "dt_track_step",
          &SpaceChargeCommonProxy::dt_track_step,
          &SpaceChargeCommonProxy::set_dt_track_step)
      // SpaceChargeCommonProxy.cathode_strength_cutoff (0D_NOT_real - Cutoff for the cathode field calc.
      .def_property(
          "cathode_strength_cutoff",
          &SpaceChargeCommonProxy::cathode_strength_cutoff,
          &SpaceChargeCommonProxy::set_cathode_strength_cutoff)
      // SpaceChargeCommonProxy.rel_tol_tracking (0D_NOT_real - Relative tolerance for tracking.
      .def_property(
          "rel_tol_tracking",
          &SpaceChargeCommonProxy::rel_tol_tracking,
          &SpaceChargeCommonProxy::set_rel_tol_tracking)
      // SpaceChargeCommonProxy.abs_tol_tracking (0D_NOT_real - Absolute tolerance for tracking.
      .def_property(
          "abs_tol_tracking",
          &SpaceChargeCommonProxy::abs_tol_tracking,
          &SpaceChargeCommonProxy::set_abs_tol_tracking)
      // SpaceChargeCommonProxy.beam_chamber_height (0D_NOT_real - Used in shielding calculation.
      .def_property(
          "beam_chamber_height",
          &SpaceChargeCommonProxy::beam_chamber_height,
          &SpaceChargeCommonProxy::set_beam_chamber_height)
      // SpaceChargeCommonProxy.lsc_sigma_cutoff (0D_NOT_real - Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then ignore.
      .def_property(
          "lsc_sigma_cutoff",
          &SpaceChargeCommonProxy::lsc_sigma_cutoff,
          &SpaceChargeCommonProxy::set_lsc_sigma_cutoff)
      // SpaceChargeCommonProxy.particle_sigma_cutoff (0D_NOT_real - 3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or zero means ignore.
      .def_property(
          "particle_sigma_cutoff",
          &SpaceChargeCommonProxy::particle_sigma_cutoff,
          &SpaceChargeCommonProxy::set_particle_sigma_cutoff)
      // SpaceChargeCommonProxy.space_charge_mesh_size (1D_NOT_integer - Gird size for fft_3d space charge calc.
      .def_property_readonly(
          "space_charge_mesh_size",
          &SpaceChargeCommonProxy::space_charge_mesh_size)
      // SpaceChargeCommonProxy.csr3d_mesh_size (1D_NOT_integer - Gird size for CSR.
      .def_property_readonly(
          "csr3d_mesh_size", &SpaceChargeCommonProxy::csr3d_mesh_size)
      // SpaceChargeCommonProxy.n_bin (0D_NOT_integer - Number of bins used
      .def_property(
          "n_bin",
          &SpaceChargeCommonProxy::n_bin,
          &SpaceChargeCommonProxy::set_n_bin)
      // SpaceChargeCommonProxy.particle_bin_span (0D_NOT_integer - Longitudinal particle length / dz_bin
      .def_property(
          "particle_bin_span",
          &SpaceChargeCommonProxy::particle_bin_span,
          &SpaceChargeCommonProxy::set_particle_bin_span)
      // SpaceChargeCommonProxy.n_shield_images (0D_NOT_integer - Chamber wall shielding. 0 = no shielding.
      .def_property(
          "n_shield_images",
          &SpaceChargeCommonProxy::n_shield_images,
          &SpaceChargeCommonProxy::set_n_shield_images)
      // SpaceChargeCommonProxy.sc_min_in_bin (0D_NOT_integer - Minimum number of particles in a bin for sigmas to be valid.
      .def_property(
          "sc_min_in_bin",
          &SpaceChargeCommonProxy::sc_min_in_bin,
          &SpaceChargeCommonProxy::set_sc_min_in_bin)
      // SpaceChargeCommonProxy.lsc_kick_transverse_dependence (0D_NOT_logical -
      .def_property(
          "lsc_kick_transverse_dependence",
          &SpaceChargeCommonProxy::lsc_kick_transverse_dependence,
          &SpaceChargeCommonProxy::set_lsc_kick_transverse_dependence)
      // SpaceChargeCommonProxy.debug (0D_NOT_logical -
      .def_property(
          "debug",
          &SpaceChargeCommonProxy::debug,
          &SpaceChargeCommonProxy::set_debug)
      // SpaceChargeCommonProxy.diagnostic_output_file (0D_NOT_character - If non-blank write a diagnostic (EG wake) file
      .def_property(
          "diagnostic_output_file",
          &SpaceChargeCommonProxy::diagnostic_output_file,
          &SpaceChargeCommonProxy::set_diagnostic_output_file)

      .def(
          "__repr__",
          [](const SpaceChargeCommonProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpaceChargeCommonProxy& self) {
            return SpaceChargeCommonProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpaceChargeCommonProxy& self, py::dict& memo) {
            return SpaceChargeCommonProxy(self);
          })

      ;

  // 1D SpaceChargeCommonProxy arrays are not used in structs/routines
  // 2D SpaceChargeCommonProxy arrays are not used in structs/routines
  // 3D SpaceChargeCommonProxy arrays are not used in structs/routines
}

// =============================================================================
// spin_axis_struct
void init_spin_axis_struct(py::module& m, py::class_<SpinAxisProxy>& cls) {
  cls.def(py::init<>())
      // SpinAxisProxy.l (1D_NOT_real - Transverse axis.
      .def_property_readonly("l", &SpinAxisProxy::l)
      // SpinAxisProxy.n0 (1D_NOT_real - Invariant spin axis on closed orbit.
      .def_property_readonly("n0", &SpinAxisProxy::n0)
      // SpinAxisProxy.m (1D_NOT_real - Transverse axis.
      .def_property_readonly("m", &SpinAxisProxy::m)

      .def(
          "__repr__", [](const SpinAxisProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinAxisProxy& self) {
            return SpinAxisProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinAxisProxy& self, py::dict& memo) {
            return SpinAxisProxy(self);
          })

      ;

  // 1D SpinAxisProxy arrays are not used in structs/routines
  // 2D SpinAxisProxy arrays are not used in structs/routines
  // 3D SpinAxisProxy arrays are not used in structs/routines
}

// =============================================================================
// spin_orbit_map1_struct
void init_spin_orbit_map1_struct(
    py::module& m,
    py::class_<SpinOrbitMap1Proxy>& cls) {
  cls.def(py::init<>())
      // SpinOrbitMap1Proxy.orb_mat (2D_NOT_real - Orbital matrix
      .def_property_readonly("orb_mat", &SpinOrbitMap1Proxy::orb_mat)
      // SpinOrbitMap1Proxy.vec0 (1D_NOT_real - Orbital 0th order map: r_out = mat6 * r_in + vec0
      .def_property_readonly("vec0", &SpinOrbitMap1Proxy::vec0)
      // SpinOrbitMap1Proxy.spin_q (2D_NOT_real - 0th and 1st order quaternion spin map
      .def_property_readonly("spin_q", &SpinOrbitMap1Proxy::spin_q)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return SpinOrbitMap1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const SpinOrbitMap1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinOrbitMap1Proxy& self) {
            return SpinOrbitMap1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinOrbitMap1Proxy& self, py::dict& memo) {
            return SpinOrbitMap1Proxy(self);
          })

      ;

  bind_FTypeArrayND<SpinOrbitMap1ProxyArray1D>(m, "SpinOrbitMap1StructArray1D");
  bind_FTypeAlloc1D<SpinOrbitMap1ProxyAlloc1D>(m, "SpinOrbitMap1StructAlloc1D");
  // 2D SpinOrbitMap1Proxy arrays are not used in structs/routines
  // 3D SpinOrbitMap1Proxy arrays are not used in structs/routines
}

// =============================================================================
// spin_polar_struct
void init_spin_polar_struct(py::module& m, py::class_<SpinPolarProxy>& cls) {
  cls.def(py::init<>())
      // SpinPolarProxy.polarization (0D_NOT_real -
      .def_property(
          "polarization",
          &SpinPolarProxy::polarization,
          &SpinPolarProxy::set_polarization)
      // SpinPolarProxy.theta (0D_NOT_real - Spherical coords: Angle from z-axis.
      .def_property("theta", &SpinPolarProxy::theta, &SpinPolarProxy::set_theta)
      // SpinPolarProxy.phi (0D_NOT_real - Spherical coords: Angle in (x,y) plane.
      .def_property("phi", &SpinPolarProxy::phi, &SpinPolarProxy::set_phi)
      // SpinPolarProxy.xi (0D_NOT_real - Spinor phase angle (See Bmad manual).
      .def_property("xi", &SpinPolarProxy::xi, &SpinPolarProxy::set_xi)

      .def(
          "__repr__",
          [](const SpinPolarProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinPolarProxy& self) {
            return SpinPolarProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinPolarProxy& self, py::dict& memo) {
            return SpinPolarProxy(self);
          })

      ;

  // 1D SpinPolarProxy arrays are not used in structs/routines
  // 2D SpinPolarProxy arrays are not used in structs/routines
  // 3D SpinPolarProxy arrays are not used in structs/routines
}

// =============================================================================
// strong_beam_struct
void init_strong_beam_struct(py::module& m, py::class_<StrongBeamProxy>& cls) {
  cls.def(py::init<>())
      // StrongBeamProxy.ix_slice (0D_NOT_integer - 0 -> at element center and not at slice.
      .def_property(
          "ix_slice",
          &StrongBeamProxy::ix_slice,
          &StrongBeamProxy::set_ix_slice)
      // StrongBeamProxy.x_center (0D_NOT_real - Strong beam slice center.
      .def_property(
          "x_center",
          &StrongBeamProxy::x_center,
          &StrongBeamProxy::set_x_center)
      // StrongBeamProxy.y_center (0D_NOT_real - Strong beam slice center.
      .def_property(
          "y_center",
          &StrongBeamProxy::y_center,
          &StrongBeamProxy::set_y_center)
      // StrongBeamProxy.x_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property(
          "x_sigma", &StrongBeamProxy::x_sigma, &StrongBeamProxy::set_x_sigma)
      // StrongBeamProxy.y_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property(
          "y_sigma", &StrongBeamProxy::y_sigma, &StrongBeamProxy::set_y_sigma)
      // StrongBeamProxy.dx (0D_NOT_real - Particle - beam slice distance.
      .def_property("dx", &StrongBeamProxy::dx, &StrongBeamProxy::set_dx)
      // StrongBeamProxy.dy (0D_NOT_real - Particle - beam slice distance.
      .def_property("dy", &StrongBeamProxy::dy, &StrongBeamProxy::set_dy)

      .def(
          "__repr__",
          [](const StrongBeamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const StrongBeamProxy& self) {
            return StrongBeamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const StrongBeamProxy& self, py::dict& memo) {
            return StrongBeamProxy(self);
          })

      ;

  // 1D StrongBeamProxy arrays are not used in structs/routines
  // 2D StrongBeamProxy arrays are not used in structs/routines
  // 3D StrongBeamProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_curvature_struct
void init_surface_curvature_struct(
    py::module& m,
    py::class_<SurfaceCurvatureProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceCurvatureProxy.xy (2D_NOT_real -
      .def_property_readonly("xy", &SurfaceCurvatureProxy::xy)
      // SurfaceCurvatureProxy.spherical (0D_NOT_real -
      .def_property(
          "spherical",
          &SurfaceCurvatureProxy::spherical,
          &SurfaceCurvatureProxy::set_spherical)
      // SurfaceCurvatureProxy.elliptical (1D_NOT_real - Total curvature = elliptical + spherical
      .def_property_readonly("elliptical", &SurfaceCurvatureProxy::elliptical)
      // SurfaceCurvatureProxy.has_curvature (0D_NOT_logical - Dependent var. Will be set by Bmad
      .def_property(
          "has_curvature",
          &SurfaceCurvatureProxy::has_curvature,
          &SurfaceCurvatureProxy::set_has_curvature)

      .def(
          "__repr__",
          [](const SurfaceCurvatureProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceCurvatureProxy& self) {
            return SurfaceCurvatureProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceCurvatureProxy& self, py::dict& memo) {
            return SurfaceCurvatureProxy(self);
          })

      ;

  // 1D SurfaceCurvatureProxy arrays are not used in structs/routines
  // 2D SurfaceCurvatureProxy arrays are not used in structs/routines
  // 3D SurfaceCurvatureProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_pt_struct
void init_surface_displacement_pt_struct(
    py::module& m,
    py::class_<SurfaceDisplacementPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceDisplacementPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0",
          &SurfaceDisplacementPtProxy::x0,
          &SurfaceDisplacementPtProxy::set_x0)
      // SurfaceDisplacementPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0",
          &SurfaceDisplacementPtProxy::y0,
          &SurfaceDisplacementPtProxy::set_y0)
      // SurfaceDisplacementPtProxy.z0 (0D_NOT_real -
      .def_property(
          "z0",
          &SurfaceDisplacementPtProxy::z0,
          &SurfaceDisplacementPtProxy::set_z0)
      // SurfaceDisplacementPtProxy.dz_dx (0D_NOT_real -
      .def_property(
          "dz_dx",
          &SurfaceDisplacementPtProxy::dz_dx,
          &SurfaceDisplacementPtProxy::set_dz_dx)
      // SurfaceDisplacementPtProxy.dz_dy (0D_NOT_real -
      .def_property(
          "dz_dy",
          &SurfaceDisplacementPtProxy::dz_dy,
          &SurfaceDisplacementPtProxy::set_dz_dy)
      // SurfaceDisplacementPtProxy.d2z_dxdy (0D_NOT_real -
      .def_property(
          "d2z_dxdy",
          &SurfaceDisplacementPtProxy::d2z_dxdy,
          &SurfaceDisplacementPtProxy::set_d2z_dxdy)

      .def(
          "__repr__",
          [](const SurfaceDisplacementPtProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const SurfaceDisplacementPtProxy& self) {
            return SurfaceDisplacementPtProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementPtProxy& self, py::dict& memo) {
            return SurfaceDisplacementPtProxy(self);
          })

      ;

  // 1D SurfaceDisplacementPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceDisplacementPtProxyArray2D>(
      m, "SurfaceDisplacementPtStructArray2D");
  // 3D SurfaceDisplacementPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_struct
void init_surface_displacement_struct(
    py::module& m,
    py::class_<SurfaceDisplacementProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceDisplacementProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceDisplacementProxy::active,
          &SurfaceDisplacementProxy::set_active)
      // SurfaceDisplacementProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceDisplacementProxy::dr)
      // SurfaceDisplacementProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceDisplacementProxy::r0)
      // SurfaceDisplacementProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceDisplacementProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceDisplacementProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceDisplacementProxy& self) {
            return SurfaceDisplacementProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementProxy& self, py::dict& memo) {
            return SurfaceDisplacementProxy(self);
          })

      ;

  // 1D SurfaceDisplacementProxy arrays are not used in structs/routines
  // 2D SurfaceDisplacementProxy arrays are not used in structs/routines
  // 3D SurfaceDisplacementProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_pt_struct
void init_surface_h_misalign_pt_struct(
    py::module& m,
    py::class_<SurfaceHMisalignPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0", &SurfaceHMisalignPtProxy::x0, &SurfaceHMisalignPtProxy::set_x0)
      // SurfaceHMisalignPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0", &SurfaceHMisalignPtProxy::y0, &SurfaceHMisalignPtProxy::set_y0)
      // SurfaceHMisalignPtProxy.rot_y (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_y",
          &SurfaceHMisalignPtProxy::rot_y,
          &SurfaceHMisalignPtProxy::set_rot_y)
      // SurfaceHMisalignPtProxy.rot_t (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_t",
          &SurfaceHMisalignPtProxy::rot_t,
          &SurfaceHMisalignPtProxy::set_rot_t)
      // SurfaceHMisalignPtProxy.rot_y_rms (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_y_rms",
          &SurfaceHMisalignPtProxy::rot_y_rms,
          &SurfaceHMisalignPtProxy::set_rot_y_rms)
      // SurfaceHMisalignPtProxy.rot_t_rms (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_t_rms",
          &SurfaceHMisalignPtProxy::rot_t_rms,
          &SurfaceHMisalignPtProxy::set_rot_t_rms)

      .def(
          "__repr__",
          [](const SurfaceHMisalignPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignPtProxy& self) {
            return SurfaceHMisalignPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignPtProxy& self, py::dict& memo) {
            return SurfaceHMisalignPtProxy(self);
          })

      ;

  // 1D SurfaceHMisalignPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceHMisalignPtProxyArray2D>(
      m, "SurfaceHMisalignPtStructArray2D");
  // 3D SurfaceHMisalignPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_struct
void init_surface_h_misalign_struct(
    py::module& m,
    py::class_<SurfaceHMisalignProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceHMisalignProxy::active,
          &SurfaceHMisalignProxy::set_active)
      // SurfaceHMisalignProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceHMisalignProxy::dr)
      // SurfaceHMisalignProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceHMisalignProxy::r0)
      // SurfaceHMisalignProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceHMisalignProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceHMisalignProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignProxy& self) {
            return SurfaceHMisalignProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignProxy& self, py::dict& memo) {
            return SurfaceHMisalignProxy(self);
          })

      ;

  // 1D SurfaceHMisalignProxy arrays are not used in structs/routines
  // 2D SurfaceHMisalignProxy arrays are not used in structs/routines
  // 3D SurfaceHMisalignProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_pt_struct
void init_surface_segmented_pt_struct(
    py::module& m,
    py::class_<SurfaceSegmentedPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0", &SurfaceSegmentedPtProxy::x0, &SurfaceSegmentedPtProxy::set_x0)
      // SurfaceSegmentedPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0", &SurfaceSegmentedPtProxy::y0, &SurfaceSegmentedPtProxy::set_y0)
      // SurfaceSegmentedPtProxy.z0 (0D_NOT_real - Position at center
      .def_property(
          "z0", &SurfaceSegmentedPtProxy::z0, &SurfaceSegmentedPtProxy::set_z0)
      // SurfaceSegmentedPtProxy.dz_dx (0D_NOT_real - Slope at center
      .def_property(
          "dz_dx",
          &SurfaceSegmentedPtProxy::dz_dx,
          &SurfaceSegmentedPtProxy::set_dz_dx)
      // SurfaceSegmentedPtProxy.dz_dy (0D_NOT_real - Slope at center
      .def_property(
          "dz_dy",
          &SurfaceSegmentedPtProxy::dz_dy,
          &SurfaceSegmentedPtProxy::set_dz_dy)

      .def(
          "__repr__",
          [](const SurfaceSegmentedPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedPtProxy& self) {
            return SurfaceSegmentedPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedPtProxy& self, py::dict& memo) {
            return SurfaceSegmentedPtProxy(self);
          })

      ;

  // 1D SurfaceSegmentedPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceSegmentedPtProxyArray2D>(
      m, "SurfaceSegmentedPtStructArray2D");
  // 3D SurfaceSegmentedPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_struct
void init_surface_segmented_struct(
    py::module& m,
    py::class_<SurfaceSegmentedProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceSegmentedProxy::active,
          &SurfaceSegmentedProxy::set_active)
      // SurfaceSegmentedProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceSegmentedProxy::dr)
      // SurfaceSegmentedProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceSegmentedProxy::r0)
      // SurfaceSegmentedProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceSegmentedProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceSegmentedProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedProxy& self) {
            return SurfaceSegmentedProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedProxy& self, py::dict& memo) {
            return SurfaceSegmentedProxy(self);
          })

      ;

  // 1D SurfaceSegmentedProxy arrays are not used in structs/routines
  // 2D SurfaceSegmentedProxy arrays are not used in structs/routines
  // 3D SurfaceSegmentedProxy arrays are not used in structs/routines
}

// =============================================================================
// target_point_struct
void init_target_point_struct(
    py::module& m,
    py::class_<TargetPointProxy>& cls) {
  cls.def(py::init<>())
      // TargetPointProxy.r (1D_NOT_real - (x, y, z)
      .def_property_readonly("r", &TargetPointProxy::r)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TargetPointProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TargetPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TargetPointProxy& self) {
            return TargetPointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TargetPointProxy& self, py::dict& memo) {
            return TargetPointProxy(self);
          })

      ;

  bind_FTypeArrayND<TargetPointProxyArray1D>(m, "TargetPointStructArray1D");
  bind_FTypeAlloc1D<TargetPointProxyAlloc1D>(m, "TargetPointStructAlloc1D");
  // 2D TargetPointProxy arrays are not used in structs/routines
  // 3D TargetPointProxy arrays are not used in structs/routines
}

// =============================================================================
// taylor_struct
void init_taylor_struct(py::module& m, py::class_<TaylorProxy>& cls) {
  cls.def(py::init<>())
      // TaylorProxy.ref (0D_NOT_real -
      .def_property("ref", &TaylorProxy::ref, &TaylorProxy::set_ref)
      // TaylorProxy.term (1D_PTR_type -
      .def_property_readonly("term", &TaylorProxy::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaylorProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const TaylorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaylorProxy& self) {
            return TaylorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaylorProxy& self, py::dict& memo) {
            return TaylorProxy(self);
          })

      ;

  bind_FTypeArrayND<TaylorProxyArray1D>(m, "TaylorStructArray1D");
  bind_FTypeAlloc1D<TaylorProxyAlloc1D>(m, "TaylorStructAlloc1D");
  // 2D TaylorProxy arrays are not used in structs/routines
  // 3D TaylorProxy arrays are not used in structs/routines
}

// =============================================================================
// taylor_term_struct
void init_taylor_term_struct(py::module& m, py::class_<TaylorTermProxy>& cls) {
  cls.def(py::init<>())
      // TaylorTermProxy.coef (0D_NOT_real -
      .def_property("coef", &TaylorTermProxy::coef, &TaylorTermProxy::set_coef)
      // TaylorTermProxy.expn (1D_NOT_integer -
      .def_property_readonly("expn", &TaylorTermProxy::expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaylorTermProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaylorTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaylorTermProxy& self) {
            return TaylorTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaylorTermProxy& self, py::dict& memo) {
            return TaylorTermProxy(self);
          })

      ;

  bind_FTypeArrayND<TaylorTermProxyArray1D>(m, "TaylorTermStructArray1D");
  bind_FTypeAlloc1D<TaylorTermProxyAlloc1D>(m, "TaylorTermStructAlloc1D");
  // 2D TaylorTermProxy arrays are not used in structs/routines
  // 3D TaylorTermProxy arrays are not used in structs/routines
}

// =============================================================================
// track_point_struct
void init_track_point_struct(py::module& m, py::class_<TrackPointProxy>& cls) {
  cls.def(py::init<>())
      // TrackPointProxy.s_lab (0D_NOT_real - Longitudinal lab coord with respect to the upstream end.
      .def_property(
          "s_lab", &TrackPointProxy::s_lab, &TrackPointProxy::set_s_lab)
      // TrackPointProxy.s_body (0D_NOT_real - Longitudinal body coord with respect to the entrance end.
      .def_property(
          "s_body", &TrackPointProxy::s_body, &TrackPointProxy::set_s_body)
      // TrackPointProxy.orb (0D_NOT_type - Particle position in lab coords.
      .def_property("orb", &TrackPointProxy::orb, &TrackPointProxy::set_orb)
      // TrackPointProxy.field (0D_NOT_type - E&M fields in lab coordinates.
      .def_property(
          "field", &TrackPointProxy::field, &TrackPointProxy::set_field)
      // TrackPointProxy.strong_beam (0D_NOT_type - Strong beam info for beambeam element.
      .def_property(
          "strong_beam",
          &TrackPointProxy::strong_beam,
          &TrackPointProxy::set_strong_beam)
      // TrackPointProxy.vec0 (1D_NOT_real - 0th order part of xfer map from the beginning.
      .def_property_readonly("vec0", &TrackPointProxy::vec0)
      // TrackPointProxy.mat6 (2D_NOT_real - 1st order part of xfer map (transfer matrix).
      .def_property_readonly("mat6", &TrackPointProxy::mat6)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TrackPointProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TrackPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TrackPointProxy& self) {
            return TrackPointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TrackPointProxy& self, py::dict& memo) {
            return TrackPointProxy(self);
          })

      ;

  bind_FTypeArrayND<TrackPointProxyArray1D>(m, "TrackPointStructArray1D");
  bind_FTypeAlloc1D<TrackPointProxyAlloc1D>(m, "TrackPointStructAlloc1D");
  // 2D TrackPointProxy arrays are not used in structs/routines
  // 3D TrackPointProxy arrays are not used in structs/routines
}

// =============================================================================
// track_struct
void init_track_struct(py::module& m, py::class_<TrackProxy>& cls) {
  cls.def(py::init<>())
      // TrackProxy.pt (1D_ALLOC_type - Array of track points indexed from 0.
      .def_property_readonly("pt", &TrackProxy::pt)
      // TrackProxy.ds_save (0D_NOT_real - Min distance between points. Not positive => Save at all points.
      .def_property("ds_save", &TrackProxy::ds_save, &TrackProxy::set_ds_save)
      // TrackProxy.n_pt (0D_NOT_integer - Track upper bound for %pt(0:) array. n_bad and n_ok are used by adaptive trackers to record the number of times the step length had to be shortened.
      .def_property("n_pt", &TrackProxy::n_pt, &TrackProxy::set_n_pt)
      // TrackProxy.n_bad (0D_NOT_integer - Number of 'bad' steps where the step length was shortened.
      .def_property("n_bad", &TrackProxy::n_bad, &TrackProxy::set_n_bad)
      // TrackProxy.n_ok (0D_NOT_integer - Number of 'good' steps where the step length was not shortened.
      .def_property("n_ok", &TrackProxy::n_ok, &TrackProxy::set_n_ok)

      .def("__repr__", [](const TrackProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TrackProxy& self) {
            return TrackProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TrackProxy& self, py::dict& memo) {
            return TrackProxy(self);
          })

      ;

  // 1D TrackProxy arrays are not used in structs/routines
  // 2D TrackProxy arrays are not used in structs/routines
  // 3D TrackProxy arrays are not used in structs/routines
}

// =============================================================================
// twiss_struct
void init_twiss_struct(py::module& m, py::class_<TwissProxy>& cls) {
  cls.def(py::init<>())
      // TwissProxy.beta (0D_NOT_real -
      .def_property("beta", &TwissProxy::beta, &TwissProxy::set_beta)
      // TwissProxy.alpha (0D_NOT_real -
      .def_property("alpha", &TwissProxy::alpha, &TwissProxy::set_alpha)
      // TwissProxy.gamma (0D_NOT_real -
      .def_property("gamma", &TwissProxy::gamma, &TwissProxy::set_gamma)
      // TwissProxy.phi (0D_NOT_real -
      .def_property("phi", &TwissProxy::phi, &TwissProxy::set_phi)
      // TwissProxy.eta (0D_NOT_real -
      .def_property("eta", &TwissProxy::eta, &TwissProxy::set_eta)
      // TwissProxy.etap (0D_NOT_real -
      .def_property("etap", &TwissProxy::etap, &TwissProxy::set_etap)
      // TwissProxy.deta_ds (0D_NOT_real -
      .def_property("deta_ds", &TwissProxy::deta_ds, &TwissProxy::set_deta_ds)
      // TwissProxy.sigma (0D_NOT_real -
      .def_property("sigma", &TwissProxy::sigma, &TwissProxy::set_sigma)
      // TwissProxy.sigma_p (0D_NOT_real -
      .def_property("sigma_p", &TwissProxy::sigma_p, &TwissProxy::set_sigma_p)
      // TwissProxy.emit (0D_NOT_real -
      .def_property("emit", &TwissProxy::emit, &TwissProxy::set_emit)
      // TwissProxy.norm_emit (0D_NOT_real -
      .def_property(
          "norm_emit", &TwissProxy::norm_emit, &TwissProxy::set_norm_emit)
      // TwissProxy.chrom (0D_NOT_real -
      .def_property("chrom", &TwissProxy::chrom, &TwissProxy::set_chrom)
      // TwissProxy.dbeta_dpz (0D_NOT_real -
      .def_property(
          "dbeta_dpz", &TwissProxy::dbeta_dpz, &TwissProxy::set_dbeta_dpz)
      // TwissProxy.dalpha_dpz (0D_NOT_real -
      .def_property(
          "dalpha_dpz", &TwissProxy::dalpha_dpz, &TwissProxy::set_dalpha_dpz)
      // TwissProxy.deta_dpz (0D_NOT_real -
      .def_property(
          "deta_dpz", &TwissProxy::deta_dpz, &TwissProxy::set_deta_dpz)
      // TwissProxy.detap_dpz (0D_NOT_real -
      .def_property(
          "detap_dpz", &TwissProxy::detap_dpz, &TwissProxy::set_detap_dpz)

      .def("__repr__", [](const TwissProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TwissProxy& self) {
            return TwissProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TwissProxy& self, py::dict& memo) {
            return TwissProxy(self);
          })

      ;

  // 1D TwissProxy arrays are not used in structs/routines
  // 2D TwissProxy arrays are not used in structs/routines
  // 3D TwissProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_lr_mode_struct
void init_wake_lr_mode_struct(py::module& m, py::class_<WakeLrModeProxy>& cls) {
  cls.def(py::init<>())
      // WakeLrModeProxy.freq (0D_NOT_real - Actual Frequency in Hz.
      .def_property("freq", &WakeLrModeProxy::freq, &WakeLrModeProxy::set_freq)
      // WakeLrModeProxy.freq_in (0D_NOT_real - Input frequency in Hz.
      .def_property(
          "freq_in", &WakeLrModeProxy::freq_in, &WakeLrModeProxy::set_freq_in)
      // WakeLrModeProxy.R_over_Q (0D_NOT_real - Strength in V/C/m^(2*m_mode).
      .def_property(
          "R_over_Q",
          &WakeLrModeProxy::R_over_Q,
          &WakeLrModeProxy::set_R_over_Q)
      // WakeLrModeProxy.Q (0D_NOT_real - Used for backwards compatability.
      .def_property("Q", &WakeLrModeProxy::Q, &WakeLrModeProxy::set_Q)
      // WakeLrModeProxy.damp (0D_NOT_real - Damping factor = omega / 2 * Q = pi * freq / Q
      .def_property("damp", &WakeLrModeProxy::damp, &WakeLrModeProxy::set_damp)
      // WakeLrModeProxy.phi (0D_NOT_real - Phase in radians/2pi.
      .def_property("phi", &WakeLrModeProxy::phi, &WakeLrModeProxy::set_phi)
      // WakeLrModeProxy.angle (0D_NOT_real - polarization angle (radians/2pi).
      .def_property(
          "angle", &WakeLrModeProxy::angle, &WakeLrModeProxy::set_angle)
      // WakeLrModeProxy.b_sin (0D_NOT_real - non-skew sin-like component of the wake.
      .def_property(
          "b_sin", &WakeLrModeProxy::b_sin, &WakeLrModeProxy::set_b_sin)
      // WakeLrModeProxy.b_cos (0D_NOT_real - non-skew cos-like component of the wake.
      .def_property(
          "b_cos", &WakeLrModeProxy::b_cos, &WakeLrModeProxy::set_b_cos)
      // WakeLrModeProxy.a_sin (0D_NOT_real - skew sin-like component of the wake.
      .def_property(
          "a_sin", &WakeLrModeProxy::a_sin, &WakeLrModeProxy::set_a_sin)
      // WakeLrModeProxy.a_cos (0D_NOT_real - skew cos-like component of the wake.
      .def_property(
          "a_cos", &WakeLrModeProxy::a_cos, &WakeLrModeProxy::set_a_cos)
      // WakeLrModeProxy.m (0D_NOT_integer - Mode order (1 = dipole, 2 = quad, etc.)
      .def_property("m", &WakeLrModeProxy::m, &WakeLrModeProxy::set_m)
      // WakeLrModeProxy.polarized (0D_NOT_logical - Polaraized mode?
      .def_property(
          "polarized",
          &WakeLrModeProxy::polarized,
          &WakeLrModeProxy::set_polarized)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return WakeLrModeProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const WakeLrModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeLrModeProxy& self) {
            return WakeLrModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeLrModeProxy& self, py::dict& memo) {
            return WakeLrModeProxy(self);
          })

      ;

  bind_FTypeArrayND<WakeLrModeProxyArray1D>(m, "WakeLrModeStructArray1D");
  bind_FTypeAlloc1D<WakeLrModeProxyAlloc1D>(m, "WakeLrModeStructAlloc1D");
  // 2D WakeLrModeProxy arrays are not used in structs/routines
  // 3D WakeLrModeProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_lr_struct
void init_wake_lr_struct(py::module& m, py::class_<WakeLrProxy>& cls) {
  cls.def(py::init<>())
      // WakeLrProxy.file (0D_NOT_character -
      .def_property("file", &WakeLrProxy::file, &WakeLrProxy::set_file)
      // WakeLrProxy.mode (1D_ALLOC_type -
      .def_property_readonly("mode", &WakeLrProxy::mode)
      // WakeLrProxy.t_ref (0D_NOT_real - time reference value for computing the wake amplitude. This is used to prevent value overflow with long trains.
      .def_property("t_ref", &WakeLrProxy::t_ref, &WakeLrProxy::set_t_ref)
      // WakeLrProxy.freq_spread (0D_NOT_real - Random frequency spread of long range modes.
      .def_property(
          "freq_spread",
          &WakeLrProxy::freq_spread,
          &WakeLrProxy::set_freq_spread)
      // WakeLrProxy.amp_scale (0D_NOT_real - Wake amplitude scale factor.
      .def_property(
          "amp_scale", &WakeLrProxy::amp_scale, &WakeLrProxy::set_amp_scale)
      // WakeLrProxy.time_scale (0D_NOT_real - time scale factor.
      .def_property(
          "time_scale", &WakeLrProxy::time_scale, &WakeLrProxy::set_time_scale)
      // WakeLrProxy.self_wake_on (0D_NOT_logical - Long range self-wake used in tracking?
      .def_property(
          "self_wake_on",
          &WakeLrProxy::self_wake_on,
          &WakeLrProxy::set_self_wake_on)

      .def("__repr__", [](const WakeLrProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeLrProxy& self) {
            return WakeLrProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeLrProxy& self, py::dict& memo) {
            return WakeLrProxy(self);
          })

      ;

  // 1D WakeLrProxy arrays are not used in structs/routines
  // 2D WakeLrProxy arrays are not used in structs/routines
  // 3D WakeLrProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_mode_struct
void init_wake_sr_mode_struct(py::module& m, py::class_<WakeSrModeProxy>& cls) {
  cls.def(py::init<>())
      // WakeSrModeProxy.amp (0D_NOT_real - Amplitude
      .def_property("amp", &WakeSrModeProxy::amp, &WakeSrModeProxy::set_amp)
      // WakeSrModeProxy.damp (0D_NOT_real - Dampling factor.
      .def_property("damp", &WakeSrModeProxy::damp, &WakeSrModeProxy::set_damp)
      // WakeSrModeProxy.k (0D_NOT_real - k factor
      .def_property("k", &WakeSrModeProxy::k, &WakeSrModeProxy::set_k)
      // WakeSrModeProxy.phi (0D_NOT_real - Phase in radians/2pi
      .def_property("phi", &WakeSrModeProxy::phi, &WakeSrModeProxy::set_phi)
      // WakeSrModeProxy.b_sin (0D_NOT_real - non-skew (x) sin-like component of the wake
      .def_property(
          "b_sin", &WakeSrModeProxy::b_sin, &WakeSrModeProxy::set_b_sin)
      // WakeSrModeProxy.b_cos (0D_NOT_real - non-skew (x) cos-like component of the wake
      .def_property(
          "b_cos", &WakeSrModeProxy::b_cos, &WakeSrModeProxy::set_b_cos)
      // WakeSrModeProxy.a_sin (0D_NOT_real - skew (y) sin-like component of the wake
      .def_property(
          "a_sin", &WakeSrModeProxy::a_sin, &WakeSrModeProxy::set_a_sin)
      // WakeSrModeProxy.a_cos (0D_NOT_real - skew (y) cos-like component of the wake
      .def_property(
          "a_cos", &WakeSrModeProxy::a_cos, &WakeSrModeProxy::set_a_cos)
      // WakeSrModeProxy.polarization (0D_NOT_integer - Transverse: none$, x_axis$, y_axis$. Not used for longitudinal.
      .def_property(
          "polarization",
          &WakeSrModeProxy::polarization,
          &WakeSrModeProxy::set_polarization)
      // WakeSrModeProxy.position_dependence (0D_NOT_integer - Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$
      .def_property(
          "position_dependence",
          &WakeSrModeProxy::position_dependence,
          &WakeSrModeProxy::set_position_dependence)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return WakeSrModeProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const WakeSrModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrModeProxy& self) {
            return WakeSrModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeSrModeProxy& self, py::dict& memo) {
            return WakeSrModeProxy(self);
          })

      ;

  bind_FTypeArrayND<WakeSrModeProxyArray1D>(m, "WakeSrModeStructArray1D");
  bind_FTypeAlloc1D<WakeSrModeProxyAlloc1D>(m, "WakeSrModeStructAlloc1D");
  // 2D WakeSrModeProxy arrays are not used in structs/routines
  // 3D WakeSrModeProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_struct
void init_wake_sr_struct(py::module& m, py::class_<WakeSrProxy>& cls) {
  cls.def(py::init<>())
      // WakeSrProxy.file (0D_NOT_character -
      .def_property("file", &WakeSrProxy::file, &WakeSrProxy::set_file)
      // WakeSrProxy.z_long (0D_NOT_type -
      .def_property("z_long", &WakeSrProxy::z_long, &WakeSrProxy::set_z_long)
      // WakeSrProxy.long_wake (1D_ALLOC_type -
      .def_property_readonly("long_wake", &WakeSrProxy::long_wake)
      // WakeSrProxy.trans_wake (1D_ALLOC_type -
      .def_property_readonly("trans_wake", &WakeSrProxy::trans_wake)
      // WakeSrProxy.z_ref_long (0D_NOT_real - z reference value for computing the wake amplitude.
      .def_property(
          "z_ref_long", &WakeSrProxy::z_ref_long, &WakeSrProxy::set_z_ref_long)
      // WakeSrProxy.z_ref_trans (0D_NOT_real - This is used to prevent value overflow with long bunches.
      .def_property(
          "z_ref_trans",
          &WakeSrProxy::z_ref_trans,
          &WakeSrProxy::set_z_ref_trans)
      // WakeSrProxy.z_max (0D_NOT_real - Max allowable z value. 0-> ignore
      .def_property("z_max", &WakeSrProxy::z_max, &WakeSrProxy::set_z_max)
      // WakeSrProxy.amp_scale (0D_NOT_real - Wake amplitude scale factor.
      .def_property(
          "amp_scale", &WakeSrProxy::amp_scale, &WakeSrProxy::set_amp_scale)
      // WakeSrProxy.z_scale (0D_NOT_real - z-distance scale factor.
      .def_property("z_scale", &WakeSrProxy::z_scale, &WakeSrProxy::set_z_scale)
      // WakeSrProxy.scale_with_length (0D_NOT_logical - Scale wake with element length?
      .def_property(
          "scale_with_length",
          &WakeSrProxy::scale_with_length,
          &WakeSrProxy::set_scale_with_length)

      .def("__repr__", [](const WakeSrProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrProxy& self) {
            return WakeSrProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeSrProxy& self, py::dict& memo) {
            return WakeSrProxy(self);
          })

      ;

  // 1D WakeSrProxy arrays are not used in structs/routines
  // 2D WakeSrProxy arrays are not used in structs/routines
  // 3D WakeSrProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_z_long_struct
void init_wake_sr_z_long_struct(
    py::module& m,
    py::class_<WakeSrZLongProxy>& cls) {
  cls.def(py::init<>())
      // WakeSrZLongProxy.w (1D_ALLOC_real - Input single particle Wake. Indexed from 1.
      .def_property_readonly("w", &WakeSrZLongProxy::w)
      // WakeSrZLongProxy.fw (1D_ALLOC_complex - Fourier transform of w.
      .def_property_readonly("fw", &WakeSrZLongProxy::fw)
      // WakeSrZLongProxy.fbunch (1D_ALLOC_complex - Scratch space.
      .def_property_readonly("fbunch", &WakeSrZLongProxy::fbunch)
      // WakeSrZLongProxy.w_out (1D_ALLOC_complex - Scratch space.
      .def_property_readonly("w_out", &WakeSrZLongProxy::w_out)
      // WakeSrZLongProxy.dz (0D_NOT_real - Distance between points. If zero there is no wake.
      .def_property("dz", &WakeSrZLongProxy::dz, &WakeSrZLongProxy::set_dz)
      // WakeSrZLongProxy.z0 (0D_NOT_real - Wake extent is [-z0, z0].
      .def_property("z0", &WakeSrZLongProxy::z0, &WakeSrZLongProxy::set_z0)
      // WakeSrZLongProxy.smoothing_sigma (0D_NOT_real - 0 => No smoothing.
      .def_property(
          "smoothing_sigma",
          &WakeSrZLongProxy::smoothing_sigma,
          &WakeSrZLongProxy::set_smoothing_sigma)
      // WakeSrZLongProxy.position_dependence (0D_NOT_integer - Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$
      .def_property(
          "position_dependence",
          &WakeSrZLongProxy::position_dependence,
          &WakeSrZLongProxy::set_position_dependence)
      // WakeSrZLongProxy.time_based (0D_NOT_logical - Was input time based?
      .def_property(
          "time_based",
          &WakeSrZLongProxy::time_based,
          &WakeSrZLongProxy::set_time_based)

      .def(
          "__repr__",
          [](const WakeSrZLongProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrZLongProxy& self) {
            return WakeSrZLongProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeSrZLongProxy& self, py::dict& memo) {
            return WakeSrZLongProxy(self);
          })

      ;

  // 1D WakeSrZLongProxy arrays are not used in structs/routines
  // 2D WakeSrZLongProxy arrays are not used in structs/routines
  // 3D WakeSrZLongProxy arrays are not used in structs/routines
}

// =============================================================================
// wake_struct
void init_wake_struct(py::module& m, py::class_<WakeProxy>& cls) {
  cls.def(py::init<>())
      // WakeProxy.sr (0D_NOT_type - Short-range wake
      .def_property("sr", &WakeProxy::sr, &WakeProxy::set_sr)
      // WakeProxy.lr (0D_NOT_type - Long-range wake
      .def_property("lr", &WakeProxy::lr, &WakeProxy::set_lr)

      .def("__repr__", [](const WakeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeProxy& self) {
            return WakeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const WakeProxy& self, py::dict& memo) { return WakeProxy(self); })

      ;

  // 1D WakeProxy arrays are not used in structs/routines
  // 2D WakeProxy arrays are not used in structs/routines
  // 3D WakeProxy arrays are not used in structs/routines
}

// =============================================================================
// wall3d_section_struct
void init_wall3d_section_struct(
    py::module& m,
    py::class_<Wall3dSectionProxy>& cls) {
  cls.def(py::init<>())
      // Wall3dSectionProxy.name (0D_NOT_character - Identifying name
      .def_property(
          "name", &Wall3dSectionProxy::name, &Wall3dSectionProxy::set_name)
      // Wall3dSectionProxy.material (0D_NOT_character - Material.
      .def_property(
          "material",
          &Wall3dSectionProxy::material,
          &Wall3dSectionProxy::set_material)
      // Wall3dSectionProxy.v (1D_ALLOC_type - Array of vertices. Always stored relative.
      .def_property_readonly("v", &Wall3dSectionProxy::v)
      // Wall3dSectionProxy.surface (0D_PTR_type - Surface reflectivity tables.
      .def_property(
          "surface",
          &Wall3dSectionProxy::surface,
          &Wall3dSectionProxy::set_surface)
      // Wall3dSectionProxy.type (0D_NOT_integer - normal$, clear$, opaque$, wall_start$, wall_end$
      .def_property(
          "type", &Wall3dSectionProxy::type, &Wall3dSectionProxy::set_type)
      // Wall3dSectionProxy.n_vertex_input (0D_NOT_integer - Number of vertices specified by the user.
      .def_property(
          "n_vertex_input",
          &Wall3dSectionProxy::n_vertex_input,
          &Wall3dSectionProxy::set_n_vertex_input)
      // Wall3dSectionProxy.ix_ele (0D_NOT_integer - index of lattice element containing section
      .def_property(
          "ix_ele",
          &Wall3dSectionProxy::ix_ele,
          &Wall3dSectionProxy::set_ix_ele)
      // Wall3dSectionProxy.ix_branch (0D_NOT_integer - Index of branch lattice element is in.
      .def_property(
          "ix_branch",
          &Wall3dSectionProxy::ix_branch,
          &Wall3dSectionProxy::set_ix_branch)
      // Wall3dSectionProxy.vertices_state (0D_NOT_integer - absolute$, or shifted_to_relative$. If set to absolute$ on input, will be changed to shifted_to_relative$ by section initalizer.
      .def_property(
          "vertices_state",
          &Wall3dSectionProxy::vertices_state,
          &Wall3dSectionProxy::set_vertices_state)
      // Wall3dSectionProxy.patch_in_region (0D_NOT_logical - Patch element exists between this section and previous one?
      .def_property(
          "patch_in_region",
          &Wall3dSectionProxy::patch_in_region,
          &Wall3dSectionProxy::set_patch_in_region)
      // Wall3dSectionProxy.thickness (0D_NOT_real - Material thickness.
      .def_property(
          "thickness",
          &Wall3dSectionProxy::thickness,
          &Wall3dSectionProxy::set_thickness)
      // Wall3dSectionProxy.s (0D_NOT_real - Longitudinal position
      .def_property("s", &Wall3dSectionProxy::s, &Wall3dSectionProxy::set_s)
      // Wall3dSectionProxy.r0 (1D_NOT_real - Center of section Section-to-section spline interpolation of the center of the section
      .def_property_readonly("r0", &Wall3dSectionProxy::r0)
      // Wall3dSectionProxy.dx0_ds (0D_NOT_real - Center of wall derivative
      .def_property(
          "dx0_ds",
          &Wall3dSectionProxy::dx0_ds,
          &Wall3dSectionProxy::set_dx0_ds)
      // Wall3dSectionProxy.dy0_ds (0D_NOT_real - Center of wall derivative
      .def_property(
          "dy0_ds",
          &Wall3dSectionProxy::dy0_ds,
          &Wall3dSectionProxy::set_dy0_ds)
      // Wall3dSectionProxy.x0_coef (1D_NOT_real - Spline coefs for x-center
      .def_property_readonly("x0_coef", &Wall3dSectionProxy::x0_coef)
      // Wall3dSectionProxy.y0_coef (1D_NOT_real - Spline coefs for y-center Section-to_section spline interpolation of the wall.
      .def_property_readonly("y0_coef", &Wall3dSectionProxy::y0_coef)
      // Wall3dSectionProxy.dr_ds (0D_NOT_real - derivative of wall radius
      .def_property(
          "dr_ds", &Wall3dSectionProxy::dr_ds, &Wall3dSectionProxy::set_dr_ds)
      // Wall3dSectionProxy.p1_coef (1D_NOT_real - Spline coefs for p0 function
      .def_property_readonly("p1_coef", &Wall3dSectionProxy::p1_coef)
      // Wall3dSectionProxy.p2_coef (1D_NOT_real - Spline coefs for p1 function
      .def_property_readonly("p2_coef", &Wall3dSectionProxy::p2_coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return Wall3dSectionProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const Wall3dSectionProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dSectionProxy& self) {
            return Wall3dSectionProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Wall3dSectionProxy& self, py::dict& memo) {
            return Wall3dSectionProxy(self);
          })

      ;

  bind_FTypeArrayND<Wall3dSectionProxyArray1D>(m, "Wall3DSectionStructArray1D");
  bind_FTypeAlloc1D<Wall3dSectionProxyAlloc1D>(m, "Wall3DSectionStructAlloc1D");
  // 2D Wall3dSectionProxy arrays are not used in structs/routines
  // 3D Wall3dSectionProxy arrays are not used in structs/routines
}

// =============================================================================
// wall3d_struct
void init_wall3d_struct(py::module& m, py::class_<Wall3dProxy>& cls) {
  cls.def(py::init<>())
      // Wall3dProxy.name (0D_NOT_character -
      .def_property("name", &Wall3dProxy::name, &Wall3dProxy::set_name)
      // Wall3dProxy.type (0D_NOT_integer - or mask_plate$
      .def_property("type", &Wall3dProxy::type, &Wall3dProxy::set_type)
      // Wall3dProxy.ix_wall3d (0D_NOT_integer - Index in branch%wall3d(:) array.
      .def_property(
          "ix_wall3d", &Wall3dProxy::ix_wall3d, &Wall3dProxy::set_ix_wall3d)
      // Wall3dProxy.n_link (0D_NOT_integer - For memory management of ele%wall3d
      .def_property("n_link", &Wall3dProxy::n_link, &Wall3dProxy::set_n_link)
      // Wall3dProxy.thickness (0D_NOT_real - For diffraction_plate elements
      .def_property(
          "thickness", &Wall3dProxy::thickness, &Wall3dProxy::set_thickness)
      // Wall3dProxy.clear_material (0D_NOT_character -
      .def_property(
          "clear_material",
          &Wall3dProxy::clear_material,
          &Wall3dProxy::set_clear_material)
      // Wall3dProxy.opaque_material (0D_NOT_character -
      .def_property(
          "opaque_material",
          &Wall3dProxy::opaque_material,
          &Wall3dProxy::set_opaque_material)
      // Wall3dProxy.superimpose (0D_NOT_logical - Can overlap another wall
      .def_property(
          "superimpose",
          &Wall3dProxy::superimpose,
          &Wall3dProxy::set_superimpose)
      // Wall3dProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &Wall3dProxy::ele_anchor_pt,
          &Wall3dProxy::set_ele_anchor_pt)
      // Wall3dProxy.section (1D_ALLOC_type - Indexed from 1.
      .def_property_readonly("section", &Wall3dProxy::section)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return Wall3dProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const Wall3dProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dProxy& self) {
            return Wall3dProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Wall3dProxy& self, py::dict& memo) {
            return Wall3dProxy(self);
          })

      ;

  bind_FTypeArrayND<Wall3dProxyArray1D>(m, "Wall3DStructArray1D");
  bind_FTypeAlloc1D<Wall3dProxyAlloc1D>(m, "Wall3DStructAlloc1D");
  // 2D Wall3dProxy arrays are not used in structs/routines
  // 3D Wall3dProxy arrays are not used in structs/routines
}

// =============================================================================
// wall3d_vertex_struct
void init_wall3d_vertex_struct(
    py::module& m,
    py::class_<Wall3dVertexProxy>& cls) {
  cls.def(py::init<>())
      // Wall3dVertexProxy.x (0D_NOT_real - Coordinates of the vertex.
      .def_property("x", &Wall3dVertexProxy::x, &Wall3dVertexProxy::set_x)
      // Wall3dVertexProxy.y (0D_NOT_real - Coordinates of the vertex.
      .def_property("y", &Wall3dVertexProxy::y, &Wall3dVertexProxy::set_y)
      // Wall3dVertexProxy.radius_x (0D_NOT_real - Radius of arc or ellipse x-axis half width. 0 => Straight line.
      .def_property(
          "radius_x",
          &Wall3dVertexProxy::radius_x,
          &Wall3dVertexProxy::set_radius_x)
      // Wall3dVertexProxy.radius_y (0D_NOT_real - Ellipse y-axis half height.
      .def_property(
          "radius_y",
          &Wall3dVertexProxy::radius_y,
          &Wall3dVertexProxy::set_radius_y)
      // Wall3dVertexProxy.tilt (0D_NOT_real - Tilt of ellipse
      .def_property(
          "tilt", &Wall3dVertexProxy::tilt, &Wall3dVertexProxy::set_tilt)
      // Wall3dVertexProxy.angle (0D_NOT_real - Angle of (x, y) point.
      .def_property(
          "angle", &Wall3dVertexProxy::angle, &Wall3dVertexProxy::set_angle)
      // Wall3dVertexProxy.x0 (0D_NOT_real - Center of ellipse
      .def_property("x0", &Wall3dVertexProxy::x0, &Wall3dVertexProxy::set_x0)
      // Wall3dVertexProxy.y0 (0D_NOT_real - Center of ellipse
      .def_property("y0", &Wall3dVertexProxy::y0, &Wall3dVertexProxy::set_y0)
      // Wall3dVertexProxy.type (0D_NOT_integer - No longer used.
      .def_property(
          "type", &Wall3dVertexProxy::type, &Wall3dVertexProxy::set_type)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return Wall3dVertexProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const Wall3dVertexProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dVertexProxy& self) {
            return Wall3dVertexProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Wall3dVertexProxy& self, py::dict& memo) {
            return Wall3dVertexProxy(self);
          })

      ;

  bind_FTypeArrayND<Wall3dVertexProxyArray1D>(m, "Wall3DVertexStructArray1D");
  bind_FTypeAlloc1D<Wall3dVertexProxyAlloc1D>(m, "Wall3DVertexStructAlloc1D");
  // 2D Wall3dVertexProxy arrays are not used in structs/routines
  // 3D Wall3dVertexProxy arrays are not used in structs/routines
}

// =============================================================================
// xy_disp_struct
void init_xy_disp_struct(py::module& m, py::class_<XyDispProxy>& cls) {
  cls.def(py::init<>())
      // XyDispProxy.eta (0D_NOT_real -
      .def_property("eta", &XyDispProxy::eta, &XyDispProxy::set_eta)
      // XyDispProxy.etap (0D_NOT_real -
      .def_property("etap", &XyDispProxy::etap, &XyDispProxy::set_etap)
      // XyDispProxy.deta_ds (0D_NOT_real -
      .def_property("deta_ds", &XyDispProxy::deta_ds, &XyDispProxy::set_deta_ds)
      // XyDispProxy.sigma (0D_NOT_real -
      .def_property("sigma", &XyDispProxy::sigma, &XyDispProxy::set_sigma)
      // XyDispProxy.deta_dpz (0D_NOT_real -
      .def_property(
          "deta_dpz", &XyDispProxy::deta_dpz, &XyDispProxy::set_deta_dpz)
      // XyDispProxy.detap_dpz (0D_NOT_real -
      .def_property(
          "detap_dpz", &XyDispProxy::detap_dpz, &XyDispProxy::set_detap_dpz)

      .def("__repr__", [](const XyDispProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const XyDispProxy& self) {
            return XyDispProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const XyDispProxy& self, py::dict& memo) {
            return XyDispProxy(self);
          })

      ;

  // 1D XyDispProxy arrays are not used in structs/routines
  // 2D XyDispProxy arrays are not used in structs/routines
  // 3D XyDispProxy arrays are not used in structs/routines
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

// =============================================================================
// tricubic_cmplx_coef_struct
void init_tricubic_cmplx_coef_struct(
    py::module& m,
    py::class_<TricubicCmplxCoefProxy>& cls) {
  cls.def(py::init<>())
      // TricubicCmplxCoefProxy.coef (3D_NOT_complex - Coefs
      .def_property_readonly("coef", &TricubicCmplxCoefProxy::coef)
      // TricubicCmplxCoefProxy.i_box (1D_NOT_integer - index at lower box corner.
      .def_property_readonly("i_box", &TricubicCmplxCoefProxy::i_box)

      .def(
          "__repr__",
          [](const TricubicCmplxCoefProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TricubicCmplxCoefProxy& self) {
            return TricubicCmplxCoefProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TricubicCmplxCoefProxy& self, py::dict& memo) {
            return TricubicCmplxCoefProxy(self);
          })

      ;

  // 1D TricubicCmplxCoefProxy arrays are not used in structs/routines
  // 2D TricubicCmplxCoefProxy arrays are not used in structs/routines
  bind_FTypeArrayND<TricubicCmplxCoefProxyArray3D>(
      m, "TricubicCmplxCoefStructArray3D");
}

// =============================================================================
// mad_energy_struct
void init_mad_energy_struct(py::module& m, py::class_<MadEnergyProxy>& cls) {
  cls.def(py::init<>())
      // MadEnergyProxy.total (0D_NOT_real -
      .def_property("total", &MadEnergyProxy::total, &MadEnergyProxy::set_total)
      // MadEnergyProxy.beta (0D_NOT_real - normalized velocity: v/c
      .def_property("beta", &MadEnergyProxy::beta, &MadEnergyProxy::set_beta)
      // MadEnergyProxy.gamma (0D_NOT_real - relativistic factor: 1/sqrt(1-beta^2)
      .def_property("gamma", &MadEnergyProxy::gamma, &MadEnergyProxy::set_gamma)
      // MadEnergyProxy.kinetic (0D_NOT_real - kinetic energy
      .def_property(
          "kinetic", &MadEnergyProxy::kinetic, &MadEnergyProxy::set_kinetic)
      // MadEnergyProxy.p0c (0D_NOT_real - particle momentum
      .def_property("p0c", &MadEnergyProxy::p0c, &MadEnergyProxy::set_p0c)
      // MadEnergyProxy.particle (0D_NOT_integer - particle species
      .def_property(
          "particle", &MadEnergyProxy::particle, &MadEnergyProxy::set_particle)

      .def(
          "__repr__",
          [](const MadEnergyProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadEnergyProxy& self) {
            return MadEnergyProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const MadEnergyProxy& self, py::dict& memo) {
            return MadEnergyProxy(self);
          })

      ;

  // 1D MadEnergyProxy arrays are not used in structs/routines
  // 2D MadEnergyProxy arrays are not used in structs/routines
  // 3D MadEnergyProxy arrays are not used in structs/routines
}

// =============================================================================
// mad_map_struct
void init_mad_map_struct(py::module& m, py::class_<MadMapProxy>& cls) {
  cls.def(py::init<>())
      // MadMapProxy.k (1D_NOT_real - 0th order map.
      .def_property_readonly("k", &MadMapProxy::k)
      // MadMapProxy.r (2D_NOT_real - 1st order map.
      .def_property_readonly("r", &MadMapProxy::r)
      // MadMapProxy.t (3D_NOT_real - 2nd order map.
      .def_property_readonly("t", &MadMapProxy::t)

      .def("__repr__", [](const MadMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadMapProxy& self) {
            return MadMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const MadMapProxy& self, py::dict& memo) {
            return MadMapProxy(self);
          })

      ;

  // 1D MadMapProxy arrays are not used in structs/routines
  // 2D MadMapProxy arrays are not used in structs/routines
  // 3D MadMapProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_axis_struct
void init_qp_axis_struct(py::module& m, py::class_<QpAxisProxy>& cls) {
  cls.def(py::init<>())
      // QpAxisProxy.label (0D_NOT_character -
      .def_property("label", &QpAxisProxy::label, &QpAxisProxy::set_label)
      // QpAxisProxy.min (0D_NOT_real - Axis min/max in data units.
      .def_property("min", &QpAxisProxy::min, &QpAxisProxy::set_min)
      // QpAxisProxy.max (0D_NOT_real - Axis min/max in data units.
      .def_property("max", &QpAxisProxy::max, &QpAxisProxy::set_max)
      // QpAxisProxy.tick_min (0D_NOT_real - Min tick location along axis in data units.
      .def_property(
          "tick_min", &QpAxisProxy::tick_min, &QpAxisProxy::set_tick_min)
      // QpAxisProxy.tick_max (0D_NOT_real - Max tick location along axis in data units.
      .def_property(
          "tick_max", &QpAxisProxy::tick_max, &QpAxisProxy::set_tick_max)
      // QpAxisProxy.eval_min (0D_NOT_real - For general use. Not set by quick_plot.
      .def_property(
          "eval_min", &QpAxisProxy::eval_min, &QpAxisProxy::set_eval_min)
      // QpAxisProxy.eval_max (0D_NOT_real - For general use. Not set by quick_plot.
      .def_property(
          "eval_max", &QpAxisProxy::eval_max, &QpAxisProxy::set_eval_max)
      // QpAxisProxy.dtick (0D_NOT_real - Distance between ticks. In data units. Ticks will be drawn between %min and %max.
      .def_property("dtick", &QpAxisProxy::dtick, &QpAxisProxy::set_dtick)
      // QpAxisProxy.number_offset (0D_NOT_real - Offset from axis line in inches.
      .def_property(
          "number_offset",
          &QpAxisProxy::number_offset,
          &QpAxisProxy::set_number_offset)
      // QpAxisProxy.label_offset (0D_NOT_real - Offset from numbers in inches.
      .def_property(
          "label_offset",
          &QpAxisProxy::label_offset,
          &QpAxisProxy::set_label_offset)
      // QpAxisProxy.major_tick_len (0D_NOT_real - In inches.
      .def_property(
          "major_tick_len",
          &QpAxisProxy::major_tick_len,
          &QpAxisProxy::set_major_tick_len)
      // QpAxisProxy.minor_tick_len (0D_NOT_real - In inches.
      .def_property(
          "minor_tick_len",
          &QpAxisProxy::minor_tick_len,
          &QpAxisProxy::set_minor_tick_len)
      // QpAxisProxy.label_color (0D_NOT_character - Color of the label.
      .def_property(
          "label_color",
          &QpAxisProxy::label_color,
          &QpAxisProxy::set_label_color)
      // QpAxisProxy.major_div (0D_NOT_integer - Actual numbrer of major divisions
      .def_property(
          "major_div", &QpAxisProxy::major_div, &QpAxisProxy::set_major_div)
      // QpAxisProxy.major_div_nominal (0D_NOT_integer - Nominal value.
      .def_property(
          "major_div_nominal",
          &QpAxisProxy::major_div_nominal,
          &QpAxisProxy::set_major_div_nominal)
      // QpAxisProxy.minor_div (0D_NOT_integer - 0 = auto choose.
      .def_property(
          "minor_div", &QpAxisProxy::minor_div, &QpAxisProxy::set_minor_div)
      // QpAxisProxy.minor_div_max (0D_NOT_integer - Max number for auto choose.
      .def_property(
          "minor_div_max",
          &QpAxisProxy::minor_div_max,
          &QpAxisProxy::set_minor_div_max)
      // QpAxisProxy.places (0D_NOT_integer - Number of places after the decimal point to print.
      .def_property("places", &QpAxisProxy::places, &QpAxisProxy::set_places)
      // QpAxisProxy.type (0D_NOT_character - Or 'LOG', or 'CUSTOM'
      .def_property("type", &QpAxisProxy::type, &QpAxisProxy::set_type)
      // QpAxisProxy.bounds (0D_NOT_character - Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'
      .def_property("bounds", &QpAxisProxy::bounds, &QpAxisProxy::set_bounds)
      // QpAxisProxy.tick_side (0D_NOT_integer - +1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside.
      .def_property(
          "tick_side", &QpAxisProxy::tick_side, &QpAxisProxy::set_tick_side)
      // QpAxisProxy.number_side (0D_NOT_integer - +1 = Draw to the side inside the graph, -1 = outside.
      .def_property(
          "number_side",
          &QpAxisProxy::number_side,
          &QpAxisProxy::set_number_side)
      // QpAxisProxy.draw_label (0D_NOT_logical -
      .def_property(
          "draw_label", &QpAxisProxy::draw_label, &QpAxisProxy::set_draw_label)
      // QpAxisProxy.draw_numbers (0D_NOT_logical -
      .def_property(
          "draw_numbers",
          &QpAxisProxy::draw_numbers,
          &QpAxisProxy::set_draw_numbers)

      .def("__repr__", [](const QpAxisProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpAxisProxy& self) {
            return QpAxisProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpAxisProxy& self, py::dict& memo) {
            return QpAxisProxy(self);
          })

      ;

  // 1D QpAxisProxy arrays are not used in structs/routines
  // 2D QpAxisProxy arrays are not used in structs/routines
  // 3D QpAxisProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_legend_struct
void init_qp_legend_struct(py::module& m, py::class_<QpLegendProxy>& cls) {
  cls.def(py::init<>())
      // QpLegendProxy.row_spacing (0D_NOT_real - Spacing between rows.
      .def_property(
          "row_spacing",
          &QpLegendProxy::row_spacing,
          &QpLegendProxy::set_row_spacing)
      // QpLegendProxy.line_length (0D_NOT_real - Length of the line in points.
      .def_property(
          "line_length",
          &QpLegendProxy::line_length,
          &QpLegendProxy::set_line_length)
      // QpLegendProxy.text_offset (0D_NOT_real - Horizontal offset in points between the line and the text.
      .def_property(
          "text_offset",
          &QpLegendProxy::text_offset,
          &QpLegendProxy::set_text_offset)
      // QpLegendProxy.draw_line (0D_NOT_logical - Draw lines?
      .def_property(
          "draw_line", &QpLegendProxy::draw_line, &QpLegendProxy::set_draw_line)
      // QpLegendProxy.draw_symbol (0D_NOT_logical - Draw symbols?
      .def_property(
          "draw_symbol",
          &QpLegendProxy::draw_symbol,
          &QpLegendProxy::set_draw_symbol)
      // QpLegendProxy.draw_text (0D_NOT_logical - Draw text?
      .def_property(
          "draw_text", &QpLegendProxy::draw_text, &QpLegendProxy::set_draw_text)

      .def(
          "__repr__", [](const QpLegendProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLegendProxy& self) {
            return QpLegendProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpLegendProxy& self, py::dict& memo) {
            return QpLegendProxy(self);
          })

      ;

  // 1D QpLegendProxy arrays are not used in structs/routines
  // 2D QpLegendProxy arrays are not used in structs/routines
  // 3D QpLegendProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_line_struct
void init_qp_line_struct(py::module& m, py::class_<QpLineProxy>& cls) {
  cls.def(py::init<>())
      // QpLineProxy.width (0D_NOT_integer -
      .def_property("width", &QpLineProxy::width, &QpLineProxy::set_width)
      // QpLineProxy.color (0D_NOT_character -
      .def_property("color", &QpLineProxy::color, &QpLineProxy::set_color)
      // QpLineProxy.pattern (0D_NOT_character -
      .def_property("pattern", &QpLineProxy::pattern, &QpLineProxy::set_pattern)

      .def("__repr__", [](const QpLineProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLineProxy& self) {
            return QpLineProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpLineProxy& self, py::dict& memo) {
            return QpLineProxy(self);
          })

      ;

  // 1D QpLineProxy arrays are not used in structs/routines
  // 2D QpLineProxy arrays are not used in structs/routines
  // 3D QpLineProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_point_struct
void init_qp_point_struct(py::module& m, py::class_<QpPointProxy>& cls) {
  cls.def(py::init<>())
      // QpPointProxy.x (0D_NOT_real -
      .def_property("x", &QpPointProxy::x, &QpPointProxy::set_x)
      // QpPointProxy.y (0D_NOT_real -
      .def_property("y", &QpPointProxy::y, &QpPointProxy::set_y)
      // QpPointProxy.units (0D_NOT_character -
      .def_property("units", &QpPointProxy::units, &QpPointProxy::set_units)

      .def("__repr__", [](const QpPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpPointProxy& self) {
            return QpPointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpPointProxy& self, py::dict& memo) {
            return QpPointProxy(self);
          })

      ;

  // 1D QpPointProxy arrays are not used in structs/routines
  // 2D QpPointProxy arrays are not used in structs/routines
  // 3D QpPointProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_rect_struct
void init_qp_rect_struct(py::module& m, py::class_<QpRectProxy>& cls) {
  cls.def(py::init<>())
      // QpRectProxy.x1 (0D_NOT_real -
      .def_property("x1", &QpRectProxy::x1, &QpRectProxy::set_x1)
      // QpRectProxy.x2 (0D_NOT_real -
      .def_property("x2", &QpRectProxy::x2, &QpRectProxy::set_x2)
      // QpRectProxy.y1 (0D_NOT_real -
      .def_property("y1", &QpRectProxy::y1, &QpRectProxy::set_y1)
      // QpRectProxy.y2 (0D_NOT_real -
      .def_property("y2", &QpRectProxy::y2, &QpRectProxy::set_y2)
      // QpRectProxy.units (0D_NOT_character -
      .def_property("units", &QpRectProxy::units, &QpRectProxy::set_units)

      .def("__repr__", [](const QpRectProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpRectProxy& self) {
            return QpRectProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpRectProxy& self, py::dict& memo) {
            return QpRectProxy(self);
          })

      ;

  // 1D QpRectProxy arrays are not used in structs/routines
  // 2D QpRectProxy arrays are not used in structs/routines
  // 3D QpRectProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_symbol_struct
void init_qp_symbol_struct(py::module& m, py::class_<QpSymbolProxy>& cls) {
  cls.def(py::init<>())
      // QpSymbolProxy.type (0D_NOT_character -
      .def_property("type", &QpSymbolProxy::type, &QpSymbolProxy::set_type)
      // QpSymbolProxy.height (0D_NOT_real - in points (same as text height)
      .def_property(
          "height", &QpSymbolProxy::height, &QpSymbolProxy::set_height)
      // QpSymbolProxy.color (0D_NOT_character -
      .def_property("color", &QpSymbolProxy::color, &QpSymbolProxy::set_color)
      // QpSymbolProxy.fill_pattern (0D_NOT_character -
      .def_property(
          "fill_pattern",
          &QpSymbolProxy::fill_pattern,
          &QpSymbolProxy::set_fill_pattern)
      // QpSymbolProxy.line_width (0D_NOT_integer -
      .def_property(
          "line_width",
          &QpSymbolProxy::line_width,
          &QpSymbolProxy::set_line_width)

      .def(
          "__repr__", [](const QpSymbolProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpSymbolProxy& self) {
            return QpSymbolProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpSymbolProxy& self, py::dict& memo) {
            return QpSymbolProxy(self);
          })

      ;

  // 1D QpSymbolProxy arrays are not used in structs/routines
  // 2D QpSymbolProxy arrays are not used in structs/routines
  // 3D QpSymbolProxy arrays are not used in structs/routines
}

// =============================================================================
// random_state_struct
void init_random_state_struct(
    py::module& m,
    py::class_<RandomStateProxy>& cls) {
  cls.def(py::init<>())
      // RandomStateProxy.ix (0D_NOT_integer8 -
      .def_property("ix", &RandomStateProxy::ix, &RandomStateProxy::set_ix)
      // RandomStateProxy.iy (0D_NOT_integer8 -
      .def_property("iy", &RandomStateProxy::iy, &RandomStateProxy::set_iy)
      // RandomStateProxy.number_stored (0D_NOT_logical -
      .def_property(
          "number_stored",
          &RandomStateProxy::number_stored,
          &RandomStateProxy::set_number_stored)
      // RandomStateProxy.h_saved (0D_NOT_real -
      .def_property(
          "h_saved", &RandomStateProxy::h_saved, &RandomStateProxy::set_h_saved)
      // RandomStateProxy.engine (0D_NOT_integer - Params
      .def_property(
          "engine", &RandomStateProxy::engine, &RandomStateProxy::set_engine)
      // RandomStateProxy.seed (0D_NOT_integer -
      .def_property(
          "seed", &RandomStateProxy::seed, &RandomStateProxy::set_seed)
      // RandomStateProxy.am (0D_NOT_real -
      .def_property("am", &RandomStateProxy::am, &RandomStateProxy::set_am)
      // RandomStateProxy.gauss_converter (0D_NOT_integer -
      .def_property(
          "gauss_converter",
          &RandomStateProxy::gauss_converter,
          &RandomStateProxy::set_gauss_converter)
      // RandomStateProxy.gauss_sigma_cut (0D_NOT_real - Only used if positive.
      .def_property(
          "gauss_sigma_cut",
          &RandomStateProxy::gauss_sigma_cut,
          &RandomStateProxy::set_gauss_sigma_cut)
      // RandomStateProxy.in_sobseq (0D_NOT_integer8 -
      .def_property(
          "in_sobseq",
          &RandomStateProxy::in_sobseq,
          &RandomStateProxy::set_in_sobseq)
      // 1D_NOT_integer8 ix_sobseq proxy support missing
      // RandomStateProxy.x_sobseq (1D_NOT_real -
      .def_property_readonly("x_sobseq", &RandomStateProxy::x_sobseq)

      .def(
          "__repr__",
          [](const RandomStateProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RandomStateProxy& self) {
            return RandomStateProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RandomStateProxy& self, py::dict& memo) {
            return RandomStateProxy(self);
          })

      ;

  // 1D RandomStateProxy arrays are not used in structs/routines
  // 2D RandomStateProxy arrays are not used in structs/routines
  // 3D RandomStateProxy arrays are not used in structs/routines
}

// =============================================================================
// nametable_struct
void init_nametable_struct(py::module& m, py::class_<NametableProxy>& cls) {
  cls.def(py::init<>())
      // NametableProxy.name (1D_ALLOC_character - Array of names.
      .def_property_readonly("name", &NametableProxy::name)
      // NametableProxy.index (1D_ALLOC_integer - Sorted index for names(:) array. names(an_index(i)) is in alphabetical order.
      .def_property_readonly("index", &NametableProxy::index)
      // NametableProxy.n_min (0D_NOT_integer - Set to 0 for use in a lattice.
      .def_property("n_min", &NametableProxy::n_min, &NametableProxy::set_n_min)
      // NametableProxy.n_max (0D_NOT_integer - Use only names(n_min:n_max) part of array.
      .def_property("n_max", &NametableProxy::n_max, &NametableProxy::set_n_max)

      .def(
          "__repr__",
          [](const NametableProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const NametableProxy& self) {
            return NametableProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const NametableProxy& self, py::dict& memo) {
            return NametableProxy(self);
          })

      ;

  // 1D NametableProxy arrays are not used in structs/routines
  // 2D NametableProxy arrays are not used in structs/routines
  // 3D NametableProxy arrays are not used in structs/routines
}

// =============================================================================
// spline_struct
void init_spline_struct(py::module& m, py::class_<SplineProxy>& cls) {
  cls.def(py::init<>())
      // SplineProxy.x0 (0D_NOT_real - Point at start of spline
      .def_property("x0", &SplineProxy::x0, &SplineProxy::set_x0)
      // SplineProxy.y0 (0D_NOT_real - Point at start of spline
      .def_property("y0", &SplineProxy::y0, &SplineProxy::set_y0)
      // SplineProxy.x1 (0D_NOT_real - Point at end of spline
      .def_property("x1", &SplineProxy::x1, &SplineProxy::set_x1)
      // SplineProxy.coef (1D_NOT_real - coefficients for cubic spline
      .def_property_readonly("coef", &SplineProxy::coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return SplineProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const SplineProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SplineProxy& self) {
            return SplineProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SplineProxy& self, py::dict& memo) {
            return SplineProxy(self);
          })

      ;

  bind_FTypeArrayND<SplineProxyArray1D>(m, "SplineStructArray1D");
  bind_FTypeAlloc1D<SplineProxyAlloc1D>(m, "SplineStructAlloc1D");
  // 2D SplineProxy arrays are not used in structs/routines
  // 3D SplineProxy arrays are not used in structs/routines
}

// =============================================================================
// summation_rdt_struct
void init_summation_rdt_struct(
    py::module& m,
    py::class_<SummationRdtProxy>& cls) {
  cls.def(py::init<>())
      // SummationRdtProxy.h11001 (0D_NOT_complex -
      .def_property(
          "h11001", &SummationRdtProxy::h11001, &SummationRdtProxy::set_h11001)
      // SummationRdtProxy.h00111 (0D_NOT_complex -
      .def_property(
          "h00111", &SummationRdtProxy::h00111, &SummationRdtProxy::set_h00111)
      // SummationRdtProxy.h20001 (0D_NOT_complex -
      .def_property(
          "h20001", &SummationRdtProxy::h20001, &SummationRdtProxy::set_h20001)
      // SummationRdtProxy.h00201 (0D_NOT_complex -
      .def_property(
          "h00201", &SummationRdtProxy::h00201, &SummationRdtProxy::set_h00201)
      // SummationRdtProxy.h10002 (0D_NOT_complex -
      .def_property(
          "h10002", &SummationRdtProxy::h10002, &SummationRdtProxy::set_h10002)
      // SummationRdtProxy.h21000 (0D_NOT_complex -
      .def_property(
          "h21000", &SummationRdtProxy::h21000, &SummationRdtProxy::set_h21000)
      // SummationRdtProxy.h30000 (0D_NOT_complex -
      .def_property(
          "h30000", &SummationRdtProxy::h30000, &SummationRdtProxy::set_h30000)
      // SummationRdtProxy.h10110 (0D_NOT_complex -
      .def_property(
          "h10110", &SummationRdtProxy::h10110, &SummationRdtProxy::set_h10110)
      // SummationRdtProxy.h10020 (0D_NOT_complex -
      .def_property(
          "h10020", &SummationRdtProxy::h10020, &SummationRdtProxy::set_h10020)
      // SummationRdtProxy.h10200 (0D_NOT_complex - 2nd order in K2 moments
      .def_property(
          "h10200", &SummationRdtProxy::h10200, &SummationRdtProxy::set_h10200)
      // SummationRdtProxy.h31000 (0D_NOT_complex -
      .def_property(
          "h31000", &SummationRdtProxy::h31000, &SummationRdtProxy::set_h31000)
      // SummationRdtProxy.h40000 (0D_NOT_complex -
      .def_property(
          "h40000", &SummationRdtProxy::h40000, &SummationRdtProxy::set_h40000)
      // SummationRdtProxy.h20110 (0D_NOT_complex -
      .def_property(
          "h20110", &SummationRdtProxy::h20110, &SummationRdtProxy::set_h20110)
      // SummationRdtProxy.h11200 (0D_NOT_complex -
      .def_property(
          "h11200", &SummationRdtProxy::h11200, &SummationRdtProxy::set_h11200)
      // SummationRdtProxy.h20020 (0D_NOT_complex -
      .def_property(
          "h20020", &SummationRdtProxy::h20020, &SummationRdtProxy::set_h20020)
      // SummationRdtProxy.h20200 (0D_NOT_complex -
      .def_property(
          "h20200", &SummationRdtProxy::h20200, &SummationRdtProxy::set_h20200)
      // SummationRdtProxy.h00310 (0D_NOT_complex -
      .def_property(
          "h00310", &SummationRdtProxy::h00310, &SummationRdtProxy::set_h00310)
      // SummationRdtProxy.h00400 (0D_NOT_complex -
      .def_property(
          "h00400", &SummationRdtProxy::h00400, &SummationRdtProxy::set_h00400)
      // SummationRdtProxy.h22000 (0D_NOT_complex -
      .def_property(
          "h22000", &SummationRdtProxy::h22000, &SummationRdtProxy::set_h22000)
      // SummationRdtProxy.h00220 (0D_NOT_complex -
      .def_property(
          "h00220", &SummationRdtProxy::h00220, &SummationRdtProxy::set_h00220)
      // SummationRdtProxy.h11110 (0D_NOT_complex -
      .def_property(
          "h11110", &SummationRdtProxy::h11110, &SummationRdtProxy::set_h11110)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return SummationRdtProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const SummationRdtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SummationRdtProxy& self) {
            return SummationRdtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SummationRdtProxy& self, py::dict& memo) {
            return SummationRdtProxy(self);
          })

      ;

  bind_FTypeArrayND<SummationRdtProxyArray1D>(m, "SummationRdtStructArray1D");
  bind_FTypeAlloc1D<SummationRdtProxyAlloc1D>(m, "SummationRdtStructAlloc1D");
  // 2D SummationRdtProxy arrays are not used in structs/routines
  // 3D SummationRdtProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_beam_branch_struct
void init_tao_beam_branch_struct(
    py::module& m,
    py::class_<TaoBeamBranchProxy>& cls) {
  cls.def(py::init<>())
      // TaoBeamBranchProxy.beam_at_start (0D_NOT_type - Initial beam
      .def_property(
          "beam_at_start",
          &TaoBeamBranchProxy::beam_at_start,
          &TaoBeamBranchProxy::set_beam_at_start)
      // TaoBeamBranchProxy.beam_init (0D_NOT_type - User set beam distrubution at track start.
      .def_property(
          "beam_init",
          &TaoBeamBranchProxy::beam_init,
          &TaoBeamBranchProxy::set_beam_init)
      // TaoBeamBranchProxy.beam_init_used (0D_NOT_type - beam distribution with emit values set.
      .def_property(
          "beam_init_used",
          &TaoBeamBranchProxy::beam_init_used,
          &TaoBeamBranchProxy::set_beam_init_used)
      // TaoBeamBranchProxy.init_starting_distribution (0D_NOT_logical - Init beam
      .def_property(
          "init_starting_distribution",
          &TaoBeamBranchProxy::init_starting_distribution,
          &TaoBeamBranchProxy::set_init_starting_distribution)
      // TaoBeamBranchProxy.track_start (0D_NOT_character - Tracking start element.
      .def_property(
          "track_start",
          &TaoBeamBranchProxy::track_start,
          &TaoBeamBranchProxy::set_track_start)
      // TaoBeamBranchProxy.track_end (0D_NOT_character -
      .def_property(
          "track_end",
          &TaoBeamBranchProxy::track_end,
          &TaoBeamBranchProxy::set_track_end)
      // TaoBeamBranchProxy.ix_branch (0D_NOT_integer - Branch tracked. If track_start or track_end is a lord, ix_track_start/end index will be a index of slave.
      .def_property(
          "ix_branch",
          &TaoBeamBranchProxy::ix_branch,
          &TaoBeamBranchProxy::set_ix_branch)
      // TaoBeamBranchProxy.ix_track_start (0D_NOT_integer - Element track start index.
      .def_property(
          "ix_track_start",
          &TaoBeamBranchProxy::ix_track_start,
          &TaoBeamBranchProxy::set_ix_track_start)
      // TaoBeamBranchProxy.ix_track_end (0D_NOT_integer - Element track end index
      .def_property(
          "ix_track_end",
          &TaoBeamBranchProxy::ix_track_end,
          &TaoBeamBranchProxy::set_ix_track_end)

      .def(
          "__repr__",
          [](const TaoBeamBranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBeamBranchProxy& self) {
            return TaoBeamBranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBeamBranchProxy& self, py::dict& memo) {
            return TaoBeamBranchProxy(self);
          })

      ;

  // 1D TaoBeamBranchProxy arrays are not used in structs/routines
  // 2D TaoBeamBranchProxy arrays are not used in structs/routines
  // 3D TaoBeamBranchProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_beam_uni_struct
void init_tao_beam_uni_struct(py::module& m, py::class_<TaoBeamUniProxy>& cls) {
  cls.def(py::init<>())
      // TaoBeamUniProxy.saved_at (0D_NOT_character -
      .def_property(
          "saved_at",
          &TaoBeamUniProxy::saved_at,
          &TaoBeamUniProxy::set_saved_at)
      // TaoBeamUniProxy.dump_file (0D_NOT_character -
      .def_property(
          "dump_file",
          &TaoBeamUniProxy::dump_file,
          &TaoBeamUniProxy::set_dump_file)
      // TaoBeamUniProxy.dump_at (0D_NOT_character -
      .def_property(
          "dump_at", &TaoBeamUniProxy::dump_at, &TaoBeamUniProxy::set_dump_at)
      // TaoBeamUniProxy.track_beam_in_universe (0D_NOT_logical - Beam tracking enabled in this universe?
      .def_property(
          "track_beam_in_universe",
          &TaoBeamUniProxy::track_beam_in_universe,
          &TaoBeamUniProxy::set_track_beam_in_universe)
      // TaoBeamUniProxy.always_reinit (0D_NOT_logical -
      .def_property(
          "always_reinit",
          &TaoBeamUniProxy::always_reinit,
          &TaoBeamUniProxy::set_always_reinit)

      .def(
          "__repr__",
          [](const TaoBeamUniProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBeamUniProxy& self) {
            return TaoBeamUniProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBeamUniProxy& self, py::dict& memo) {
            return TaoBeamUniProxy(self);
          })

      ;

  // 1D TaoBeamUniProxy arrays are not used in structs/routines
  // 2D TaoBeamUniProxy arrays are not used in structs/routines
  // 3D TaoBeamUniProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_orientation_struct
void init_tao_building_wall_orientation_struct(
    py::module& m,
    py::class_<TaoBuildingWallOrientationProxy>& cls) {
  cls.def(py::init<>())
      // TaoBuildingWallOrientationProxy.theta (0D_NOT_real -
      .def_property(
          "theta",
          &TaoBuildingWallOrientationProxy::theta,
          &TaoBuildingWallOrientationProxy::set_theta)
      // TaoBuildingWallOrientationProxy.x_offset (0D_NOT_real -
      .def_property(
          "x_offset",
          &TaoBuildingWallOrientationProxy::x_offset,
          &TaoBuildingWallOrientationProxy::set_x_offset)
      // TaoBuildingWallOrientationProxy.z_offset (0D_NOT_real -
      .def_property(
          "z_offset",
          &TaoBuildingWallOrientationProxy::z_offset,
          &TaoBuildingWallOrientationProxy::set_z_offset)

      .def(
          "__repr__",
          [](const TaoBuildingWallOrientationProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const TaoBuildingWallOrientationProxy& self) {
            return TaoBuildingWallOrientationProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallOrientationProxy& self, py::dict& memo) {
            return TaoBuildingWallOrientationProxy(self);
          })

      ;

  // 1D TaoBuildingWallOrientationProxy arrays are not used in structs/routines
  // 2D TaoBuildingWallOrientationProxy arrays are not used in structs/routines
  // 3D TaoBuildingWallOrientationProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_point_struct
void init_tao_building_wall_point_struct(
    py::module& m,
    py::class_<TaoBuildingWallPointProxy>& cls) {
  cls.def(py::init<>())
      // TaoBuildingWallPointProxy.z (0D_NOT_real - Global floor position
      .def_property(
          "z", &TaoBuildingWallPointProxy::z, &TaoBuildingWallPointProxy::set_z)
      // TaoBuildingWallPointProxy.x (0D_NOT_real - Global floor position
      .def_property(
          "x", &TaoBuildingWallPointProxy::x, &TaoBuildingWallPointProxy::set_x)
      // TaoBuildingWallPointProxy.radius (0D_NOT_real - Arc radius. +r -> CW rotation, same as bends.
      .def_property(
          "radius",
          &TaoBuildingWallPointProxy::radius,
          &TaoBuildingWallPointProxy::set_radius)
      // TaoBuildingWallPointProxy.z_center (0D_NOT_real - Arc center.
      .def_property(
          "z_center",
          &TaoBuildingWallPointProxy::z_center,
          &TaoBuildingWallPointProxy::set_z_center)
      // TaoBuildingWallPointProxy.x_center (0D_NOT_real - Arc center.
      .def_property(
          "x_center",
          &TaoBuildingWallPointProxy::x_center,
          &TaoBuildingWallPointProxy::set_x_center)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoBuildingWallPointProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoBuildingWallPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallPointProxy& self) {
            return TaoBuildingWallPointProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallPointProxy& self, py::dict& memo) {
            return TaoBuildingWallPointProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoBuildingWallPointProxyArray1D>(
      m, "TaoBuildingWallPointStructArray1D");
  bind_FTypeAlloc1D<TaoBuildingWallPointProxyAlloc1D>(
      m, "TaoBuildingWallPointStructAlloc1D");
  // 2D TaoBuildingWallPointProxy arrays are not used in structs/routines
  // 3D TaoBuildingWallPointProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_section_struct
void init_tao_building_wall_section_struct(
    py::module& m,
    py::class_<TaoBuildingWallSectionProxy>& cls) {
  cls.def(py::init<>())
      // TaoBuildingWallSectionProxy.name (0D_NOT_character -
      .def_property(
          "name",
          &TaoBuildingWallSectionProxy::name,
          &TaoBuildingWallSectionProxy::set_name)
      // TaoBuildingWallSectionProxy.constraint (0D_NOT_character - 'left_side' or 'right_side' constraint.
      .def_property(
          "constraint",
          &TaoBuildingWallSectionProxy::constraint,
          &TaoBuildingWallSectionProxy::set_constraint)
      // TaoBuildingWallSectionProxy.point (1D_ALLOC_type -
      .def_property_readonly("point", &TaoBuildingWallSectionProxy::point)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoBuildingWallSectionProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoBuildingWallSectionProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const TaoBuildingWallSectionProxy& self) {
            return TaoBuildingWallSectionProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallSectionProxy& self, py::dict& memo) {
            return TaoBuildingWallSectionProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoBuildingWallSectionProxyArray1D>(
      m, "TaoBuildingWallSectionStructArray1D");
  bind_FTypeAlloc1D<TaoBuildingWallSectionProxyAlloc1D>(
      m, "TaoBuildingWallSectionStructAlloc1D");
  // 2D TaoBuildingWallSectionProxy arrays are not used in structs/routines
  // 3D TaoBuildingWallSectionProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_struct
void init_tao_building_wall_struct(
    py::module& m,
    py::class_<TaoBuildingWallProxy>& cls) {
  cls.def(py::init<>())
      // TaoBuildingWallProxy.orientation (0D_NOT_type -
      .def_property(
          "orientation",
          &TaoBuildingWallProxy::orientation,
          &TaoBuildingWallProxy::set_orientation)
      // TaoBuildingWallProxy.section (1D_ALLOC_type -
      .def_property_readonly("section", &TaoBuildingWallProxy::section)

      .def(
          "__repr__",
          [](const TaoBuildingWallProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallProxy& self) {
            return TaoBuildingWallProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallProxy& self, py::dict& memo) {
            return TaoBuildingWallProxy(self);
          })

      ;

  // 1D TaoBuildingWallProxy arrays are not used in structs/routines
  // 2D TaoBuildingWallProxy arrays are not used in structs/routines
  // 3D TaoBuildingWallProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_cmd_history_struct
void init_tao_cmd_history_struct(
    py::module& m,
    py::class_<TaoCmdHistoryProxy>& cls) {
  cls.def(py::init<>())
      // TaoCmdHistoryProxy.cmd (0D_ALLOC_character - The command
      .def_property(
          "cmd", &TaoCmdHistoryProxy::cmd, &TaoCmdHistoryProxy::set_cmd)
      // TaoCmdHistoryProxy.ix (0D_NOT_integer - Command index (1st command has ix = 1, etc.) Note: Commands from command files will be assigned an index.
      .def_property("ix", &TaoCmdHistoryProxy::ix, &TaoCmdHistoryProxy::set_ix)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoCmdHistoryProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoCmdHistoryProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCmdHistoryProxy& self) {
            return TaoCmdHistoryProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoCmdHistoryProxy& self, py::dict& memo) {
            return TaoCmdHistoryProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoCmdHistoryProxyArray1D>(m, "TaoCmdHistoryStructArray1D");
  bind_FTypeAlloc1D<TaoCmdHistoryProxyAlloc1D>(m, "TaoCmdHistoryStructAlloc1D");
  // 2D TaoCmdHistoryProxy arrays are not used in structs/routines
  // 3D TaoCmdHistoryProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_common_struct
void init_tao_common_struct(py::module& m, py::class_<TaoCommonProxy>& cls) {
  cls.def(py::init<>())
      // TaoCommonProxy.plot_place_buffer (1D_ALLOC_type - Used when %external_plotting is on.
      .def_property_readonly(
          "plot_place_buffer", &TaoCommonProxy::plot_place_buffer)
      // TaoCommonProxy.covar (2D_ALLOC_real -
      .def_property_readonly("covar", &TaoCommonProxy::covar)
      // TaoCommonProxy.alpha (2D_ALLOC_real -
      .def_property_readonly("alpha", &TaoCommonProxy::alpha)
      // TaoCommonProxy.dummy_target (0D_NOT_real - Dummy varaible
      .def_property(
          "dummy_target",
          &TaoCommonProxy::dummy_target,
          &TaoCommonProxy::set_dummy_target)
      // TaoCommonProxy.n_alias (0D_NOT_integer -
      .def_property(
          "n_alias", &TaoCommonProxy::n_alias, &TaoCommonProxy::set_n_alias)
      // TaoCommonProxy.cmd_file_level (0D_NOT_integer - For nested command files. 0 -> no command file.
      .def_property(
          "cmd_file_level",
          &TaoCommonProxy::cmd_file_level,
          &TaoCommonProxy::set_cmd_file_level)
      // TaoCommonProxy.ix_key_bank (0D_NOT_integer - For single mode.
      .def_property(
          "ix_key_bank",
          &TaoCommonProxy::ix_key_bank,
          &TaoCommonProxy::set_ix_key_bank)
      // TaoCommonProxy.ix_history (0D_NOT_integer - Index to latest command in the history circular buffer.
      .def_property(
          "ix_history",
          &TaoCommonProxy::ix_history,
          &TaoCommonProxy::set_ix_history)
      // TaoCommonProxy.n_history (0D_NOT_integer - Number of commands issued from beginning of starting Tao.
      .def_property(
          "n_history",
          &TaoCommonProxy::n_history,
          &TaoCommonProxy::set_n_history)
      // TaoCommonProxy.lev_loop (0D_NOT_integer - in do loop nest level
      .def_property(
          "lev_loop", &TaoCommonProxy::lev_loop, &TaoCommonProxy::set_lev_loop)
      // TaoCommonProxy.n_err_messages_printed (0D_NOT_integer - Used by tao_set_invalid to limit number of messages.
      .def_property(
          "n_err_messages_printed",
          &TaoCommonProxy::n_err_messages_printed,
          &TaoCommonProxy::set_n_err_messages_printed)
      // TaoCommonProxy.n_universes (0D_NOT_integer -
      .def_property(
          "n_universes",
          &TaoCommonProxy::n_universes,
          &TaoCommonProxy::set_n_universes)
      // TaoCommonProxy.ix_beam_track_active_element (0D_NOT_integer - Element being tracked through `tao_beam_track`.
      .def_property(
          "ix_beam_track_active_element",
          &TaoCommonProxy::ix_beam_track_active_element,
          &TaoCommonProxy::set_ix_beam_track_active_element)
      // TaoCommonProxy.cmd_file_paused (0D_NOT_logical -
      .def_property(
          "cmd_file_paused",
          &TaoCommonProxy::cmd_file_paused,
          &TaoCommonProxy::set_cmd_file_paused)
      // TaoCommonProxy.use_cmd_here (0D_NOT_logical - Used for commands recalled from the cmd history stack
      .def_property(
          "use_cmd_here",
          &TaoCommonProxy::use_cmd_here,
          &TaoCommonProxy::set_use_cmd_here)
      // TaoCommonProxy.cmd_from_cmd_file (0D_NOT_logical - was command from a command file?
      .def_property(
          "cmd_from_cmd_file",
          &TaoCommonProxy::cmd_from_cmd_file,
          &TaoCommonProxy::set_cmd_from_cmd_file)
      // TaoCommonProxy.use_saved_beam_in_tracking (0D_NOT_logical -
      .def_property(
          "use_saved_beam_in_tracking",
          &TaoCommonProxy::use_saved_beam_in_tracking,
          &TaoCommonProxy::set_use_saved_beam_in_tracking)
      // TaoCommonProxy.single_mode (0D_NOT_logical -
      .def_property(
          "single_mode",
          &TaoCommonProxy::single_mode,
          &TaoCommonProxy::set_single_mode)
      // TaoCommonProxy.combine_consecutive_elements_of_like_name (0D_NOT_logical -
      .def_property(
          "combine_consecutive_elements_of_like_name",
          &TaoCommonProxy::combine_consecutive_elements_of_like_name,
          &TaoCommonProxy::set_combine_consecutive_elements_of_like_name)
      // TaoCommonProxy.have_tracked_beam (0D_NOT_logical - Used to catch error when beam plotting without having tracked a beam.
      .def_property(
          "have_tracked_beam",
          &TaoCommonProxy::have_tracked_beam,
          &TaoCommonProxy::set_have_tracked_beam)
      // TaoCommonProxy.init_plot_needed (0D_NOT_logical - reinitialize plotting?
      .def_property(
          "init_plot_needed",
          &TaoCommonProxy::init_plot_needed,
          &TaoCommonProxy::set_init_plot_needed)
      // TaoCommonProxy.init_beam (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property(
          "init_beam",
          &TaoCommonProxy::init_beam,
          &TaoCommonProxy::set_init_beam)
      // TaoCommonProxy.init_var (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property(
          "init_var", &TaoCommonProxy::init_var, &TaoCommonProxy::set_init_var)
      // TaoCommonProxy.init_read_lat_info (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property(
          "init_read_lat_info",
          &TaoCommonProxy::init_read_lat_info,
          &TaoCommonProxy::set_init_read_lat_info)
      // TaoCommonProxy.optimizer_running (0D_NOT_logical -
      .def_property(
          "optimizer_running",
          &TaoCommonProxy::optimizer_running,
          &TaoCommonProxy::set_optimizer_running)
      // TaoCommonProxy.have_datums_using_expressions (0D_NOT_logical -
      .def_property(
          "have_datums_using_expressions",
          &TaoCommonProxy::have_datums_using_expressions,
          &TaoCommonProxy::set_have_datums_using_expressions)
      // TaoCommonProxy.print_to_terminal (0D_NOT_logical - Print command prompt to the terminal? For use with GUIs.
      .def_property(
          "print_to_terminal",
          &TaoCommonProxy::print_to_terminal,
          &TaoCommonProxy::set_print_to_terminal)
      // TaoCommonProxy.lattice_calc_done (0D_NOT_logical - Used by GUI for deciding when to refresh.
      .def_property(
          "lattice_calc_done",
          &TaoCommonProxy::lattice_calc_done,
          &TaoCommonProxy::set_lattice_calc_done)
      // TaoCommonProxy.add_measurement_noise (0D_NOT_logical - Turn off to take data derivatives.
      .def_property(
          "add_measurement_noise",
          &TaoCommonProxy::add_measurement_noise,
          &TaoCommonProxy::set_add_measurement_noise)
      // 1D_NOT_logical is_err_message_printed proxy support missing
      // TaoCommonProxy.command_arg_has_been_executed (0D_NOT_logical - Has the -command command line argument been executed?
      .def_property(
          "command_arg_has_been_executed",
          &TaoCommonProxy::command_arg_has_been_executed,
          &TaoCommonProxy::set_command_arg_has_been_executed)
      // TaoCommonProxy.all_merit_weights_positive (0D_NOT_logical -
      .def_property(
          "all_merit_weights_positive",
          &TaoCommonProxy::all_merit_weights_positive,
          &TaoCommonProxy::set_all_merit_weights_positive)
      // TaoCommonProxy.multi_turn_orbit_is_plotted (0D_NOT_logical - Is a multi_turn_orbit being plotted?
      .def_property(
          "multi_turn_orbit_is_plotted",
          &TaoCommonProxy::multi_turn_orbit_is_plotted,
          &TaoCommonProxy::set_multi_turn_orbit_is_plotted)
      // TaoCommonProxy.force_chrom_calc (0D_NOT_logical - Used by a routine to force a single chromaticity calculation.
      .def_property(
          "force_chrom_calc",
          &TaoCommonProxy::force_chrom_calc,
          &TaoCommonProxy::set_force_chrom_calc)
      // TaoCommonProxy.force_rad_int_calc (0D_NOT_logical - Used by a routine to force a single radiation integrals calculation
      .def_property(
          "force_rad_int_calc",
          &TaoCommonProxy::force_rad_int_calc,
          &TaoCommonProxy::set_force_rad_int_calc)
      // TaoCommonProxy.rad_int_ri_calc_on (0D_NOT_logical - 'Classical' radiation integrals calculation on/off.
      .def_property(
          "rad_int_ri_calc_on",
          &TaoCommonProxy::rad_int_ri_calc_on,
          &TaoCommonProxy::set_rad_int_ri_calc_on)
      // TaoCommonProxy.rad_int_6d_calc_on (0D_NOT_logical - 6D Radiation integrals calculation on/off.
      .def_property(
          "rad_int_6d_calc_on",
          &TaoCommonProxy::rad_int_6d_calc_on,
          &TaoCommonProxy::set_rad_int_6d_calc_on)
      // TaoCommonProxy.valid_plot_who (1D_NOT_character - model, base, ref etc...
      .def_property_readonly("valid_plot_who", &TaoCommonProxy::valid_plot_who)
      // TaoCommonProxy.single_mode_buffer (0D_NOT_character -
      .def_property(
          "single_mode_buffer",
          &TaoCommonProxy::single_mode_buffer,
          &TaoCommonProxy::set_single_mode_buffer)
      // TaoCommonProxy.cmd (0D_NOT_character - Used for the cmd history
      .def_property("cmd", &TaoCommonProxy::cmd, &TaoCommonProxy::set_cmd)
      // TaoCommonProxy.saved_cmd_line (0D_NOT_character - Saved part of command line when there are mulitple commands on a line
      .def_property(
          "saved_cmd_line",
          &TaoCommonProxy::saved_cmd_line,
          &TaoCommonProxy::set_saved_cmd_line)

      .def(
          "__repr__",
          [](const TaoCommonProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCommonProxy& self) {
            return TaoCommonProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoCommonProxy& self, py::dict& memo) {
            return TaoCommonProxy(self);
          })

      ;

  // 1D TaoCommonProxy arrays are not used in structs/routines
  // 2D TaoCommonProxy arrays are not used in structs/routines
  // 3D TaoCommonProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_color_struct
void init_tao_curve_color_struct(
    py::module& m,
    py::class_<TaoCurveColorProxy>& cls) {
  cls.def(py::init<>())
      // TaoCurveColorProxy.data_type (0D_NOT_character - Datum type to use for z-axis.
      .def_property(
          "data_type",
          &TaoCurveColorProxy::data_type,
          &TaoCurveColorProxy::set_data_type)
      // TaoCurveColorProxy.is_on (0D_NOT_logical - On/Off
      .def_property(
          "is_on", &TaoCurveColorProxy::is_on, &TaoCurveColorProxy::set_is_on)
      // TaoCurveColorProxy.min (0D_NOT_real - Min and max values for mapping z-axis to color.
      .def_property(
          "min", &TaoCurveColorProxy::min, &TaoCurveColorProxy::set_min)
      // TaoCurveColorProxy.max (0D_NOT_real - Min and max values for mapping z-axis to color.
      .def_property(
          "max", &TaoCurveColorProxy::max, &TaoCurveColorProxy::set_max)
      // TaoCurveColorProxy.autoscale (0D_NOT_logical - Set %min, %max automatically to the limits of %data_type
      .def_property(
          "autoscale",
          &TaoCurveColorProxy::autoscale,
          &TaoCurveColorProxy::set_autoscale)

      .def(
          "__repr__",
          [](const TaoCurveColorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveColorProxy& self) {
            return TaoCurveColorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoCurveColorProxy& self, py::dict& memo) {
            return TaoCurveColorProxy(self);
          })

      ;

  // 1D TaoCurveColorProxy arrays are not used in structs/routines
  // 2D TaoCurveColorProxy arrays are not used in structs/routines
  // 3D TaoCurveColorProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_orbit_struct
void init_tao_curve_orbit_struct(
    py::module& m,
    py::class_<TaoCurveOrbitProxy>& cls) {
  cls.def(py::init<>())
      // TaoCurveOrbitProxy.x (0D_NOT_real - Transverse offset
      .def_property("x", &TaoCurveOrbitProxy::x, &TaoCurveOrbitProxy::set_x)
      // TaoCurveOrbitProxy.y (0D_NOT_real - Transverse offset
      .def_property("y", &TaoCurveOrbitProxy::y, &TaoCurveOrbitProxy::set_y)
      // TaoCurveOrbitProxy.t (0D_NOT_real - Time
      .def_property("t", &TaoCurveOrbitProxy::t, &TaoCurveOrbitProxy::set_t)

      .def(
          "__repr__",
          [](const TaoCurveOrbitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveOrbitProxy& self) {
            return TaoCurveOrbitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoCurveOrbitProxy& self, py::dict& memo) {
            return TaoCurveOrbitProxy(self);
          })

      ;

  // 1D TaoCurveOrbitProxy arrays are not used in structs/routines
  // 2D TaoCurveOrbitProxy arrays are not used in structs/routines
  // 3D TaoCurveOrbitProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_struct
void init_tao_curve_struct(py::module& m, py::class_<TaoCurveProxy>& cls) {
  cls.def(py::init<>())
      // TaoCurveProxy.name (0D_NOT_character - Name identifying the curve.
      .def_property("name", &TaoCurveProxy::name, &TaoCurveProxy::set_name)
      // TaoCurveProxy.data_source (0D_NOT_character - 'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit'
      .def_property(
          "data_source",
          &TaoCurveProxy::data_source,
          &TaoCurveProxy::set_data_source)
      // TaoCurveProxy.data_index (0D_NOT_character - Used for calculating %ix_symb(:).
      .def_property(
          "data_index",
          &TaoCurveProxy::data_index,
          &TaoCurveProxy::set_data_index)
      // TaoCurveProxy.data_type_x (0D_NOT_character - Used for data slices and phase space plots.
      .def_property(
          "data_type_x",
          &TaoCurveProxy::data_type_x,
          &TaoCurveProxy::set_data_type_x)
      // TaoCurveProxy.data_type (0D_ALLOC_character - 'orbit.x', etc.
      .def_property(
          "data_type", &TaoCurveProxy::data_type, &TaoCurveProxy::set_data_type)
      // TaoCurveProxy.ele_ref_name (0D_NOT_character - Reference element.
      .def_property(
          "ele_ref_name",
          &TaoCurveProxy::ele_ref_name,
          &TaoCurveProxy::set_ele_ref_name)
      // TaoCurveProxy.legend_text (0D_NOT_character - String to draw in a curve legend.
      .def_property(
          "legend_text",
          &TaoCurveProxy::legend_text,
          &TaoCurveProxy::set_legend_text)
      // TaoCurveProxy.message_text (0D_NOT_character - Informational message to draw with graph.
      .def_property(
          "message_text",
          &TaoCurveProxy::message_text,
          &TaoCurveProxy::set_message_text)
      // TaoCurveProxy.component (0D_NOT_character - Who to plot. Eg: 'meas - design'
      .def_property(
          "component", &TaoCurveProxy::component, &TaoCurveProxy::set_component)
      // TaoCurveProxy.why_invalid (0D_NOT_character - Informative string to print.
      .def_property(
          "why_invalid",
          &TaoCurveProxy::why_invalid,
          &TaoCurveProxy::set_why_invalid)
      // TaoCurveProxy.g (0D_PTR_type - pointer to parent graph
      .def_property("g", &TaoCurveProxy::g, &TaoCurveProxy::set_g)
      // TaoCurveProxy.hist (0D_NOT_type -
      .def_property("hist", &TaoCurveProxy::hist, &TaoCurveProxy::set_hist)
      // TaoCurveProxy.z_color (0D_NOT_type -
      .def_property(
          "z_color", &TaoCurveProxy::z_color, &TaoCurveProxy::set_z_color)
      // TaoCurveProxy.x_line (1D_ALLOC_real - Coords for drawing a curve
      .def_property_readonly("x_line", &TaoCurveProxy::x_line)
      // TaoCurveProxy.y_line (1D_ALLOC_real -
      .def_property_readonly("y_line", &TaoCurveProxy::y_line)
      // TaoCurveProxy.y2_line (1D_ALLOC_real - Second array needed for beam chamber curve.
      .def_property_readonly("y2_line", &TaoCurveProxy::y2_line)
      // TaoCurveProxy.ix_line (1D_ALLOC_integer - Used by wave and aperture curves.
      .def_property_readonly("ix_line", &TaoCurveProxy::ix_line)
      // TaoCurveProxy.x_symb (1D_ALLOC_real - Coords for drawing the symbols
      .def_property_readonly("x_symb", &TaoCurveProxy::x_symb)
      // TaoCurveProxy.y_symb (1D_ALLOC_real -
      .def_property_readonly("y_symb", &TaoCurveProxy::y_symb)
      // TaoCurveProxy.z_symb (1D_ALLOC_real - Symbol color
      .def_property_readonly("z_symb", &TaoCurveProxy::z_symb)
      // TaoCurveProxy.err_symb (1D_ALLOC_real - Error bars
      .def_property_readonly("err_symb", &TaoCurveProxy::err_symb)
      // TaoCurveProxy.symb_size (1D_ALLOC_real - Symbol size. Used with symbol_size_scale.
      .def_property_readonly("symb_size", &TaoCurveProxy::symb_size)
      // TaoCurveProxy.ix_symb (1D_ALLOC_integer - Corresponding index in d1_data%d(:) array.
      .def_property_readonly("ix_symb", &TaoCurveProxy::ix_symb)
      // TaoCurveProxy.y_axis_scale_factor (0D_NOT_real - y-axis conversion from internal to plotting units.
      .def_property(
          "y_axis_scale_factor",
          &TaoCurveProxy::y_axis_scale_factor,
          &TaoCurveProxy::set_y_axis_scale_factor)
      // TaoCurveProxy.line (0D_NOT_type - Line attributes
      .def_property("line", &TaoCurveProxy::line, &TaoCurveProxy::set_line)
      // TaoCurveProxy.symbol (0D_NOT_type - Symbol attributes
      .def_property(
          "symbol", &TaoCurveProxy::symbol, &TaoCurveProxy::set_symbol)
      // TaoCurveProxy.orbit (0D_NOT_type - Used for E/B field plotting.
      .def_property("orbit", &TaoCurveProxy::orbit, &TaoCurveProxy::set_orbit)
      // TaoCurveProxy.ix_universe (0D_NOT_integer - Universe where data is. -1 => use s%global%default_universe
      .def_property(
          "ix_universe",
          &TaoCurveProxy::ix_universe,
          &TaoCurveProxy::set_ix_universe)
      // TaoCurveProxy.symbol_every (0D_NOT_integer - Symbol every how many points.
      .def_property(
          "symbol_every",
          &TaoCurveProxy::symbol_every,
          &TaoCurveProxy::set_symbol_every)
      // TaoCurveProxy.ix_branch (0D_NOT_integer -
      .def_property(
          "ix_branch", &TaoCurveProxy::ix_branch, &TaoCurveProxy::set_ix_branch)
      // TaoCurveProxy.ix_bunch (0D_NOT_integer - Bunch to plot.
      .def_property(
          "ix_bunch", &TaoCurveProxy::ix_bunch, &TaoCurveProxy::set_ix_bunch)
      // TaoCurveProxy.n_turn (0D_NOT_integer - Used for multi_turn_orbit plotting
      .def_property(
          "n_turn", &TaoCurveProxy::n_turn, &TaoCurveProxy::set_n_turn)
      // TaoCurveProxy.use_y2 (0D_NOT_logical - Use y2 axis?
      .def_property(
          "use_y2", &TaoCurveProxy::use_y2, &TaoCurveProxy::set_use_y2)
      // TaoCurveProxy.draw_line (0D_NOT_logical - Draw a line through the data points?
      .def_property(
          "draw_line", &TaoCurveProxy::draw_line, &TaoCurveProxy::set_draw_line)
      // TaoCurveProxy.draw_symbols (0D_NOT_logical - Draw a symbol at the data points?
      .def_property(
          "draw_symbols",
          &TaoCurveProxy::draw_symbols,
          &TaoCurveProxy::set_draw_symbols)
      // TaoCurveProxy.draw_symbol_index (0D_NOT_logical - Draw the symbol index number curve%ix_symb?
      .def_property(
          "draw_symbol_index",
          &TaoCurveProxy::draw_symbol_index,
          &TaoCurveProxy::set_draw_symbol_index)
      // TaoCurveProxy.draw_error_bars (0D_NOT_logical - Draw error bars based upon data%error_rms if drawing data? !! logical :: draw_rms = .false.          ! Show mean and RMS values with legend?
      .def_property(
          "draw_error_bars",
          &TaoCurveProxy::draw_error_bars,
          &TaoCurveProxy::set_draw_error_bars)
      // TaoCurveProxy.smooth_line_calc (0D_NOT_logical - Calculate data between element edge points?
      .def_property(
          "smooth_line_calc",
          &TaoCurveProxy::smooth_line_calc,
          &TaoCurveProxy::set_smooth_line_calc)
      // TaoCurveProxy.valid (0D_NOT_logical - valid data?
      .def_property("valid", &TaoCurveProxy::valid, &TaoCurveProxy::set_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoCurveProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const TaoCurveProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveProxy& self) {
            return TaoCurveProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoCurveProxy& self, py::dict& memo) {
            return TaoCurveProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoCurveProxyArray1D>(m, "TaoCurveStructArray1D");
  bind_FTypeAlloc1D<TaoCurveProxyAlloc1D>(m, "TaoCurveStructAlloc1D");
  // 2D TaoCurveProxy arrays are not used in structs/routines
  // 3D TaoCurveProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_d1_data_struct
void init_tao_d1_data_struct(py::module& m, py::class_<TaoD1DataProxy>& cls) {
  cls.def(py::init<>())
      // TaoD1DataProxy.name (0D_NOT_character - Eg: 'x', etc.
      .def_property("name", &TaoD1DataProxy::name, &TaoD1DataProxy::set_name)
      // TaoD1DataProxy.d2 (0D_PTR_type - ptr to parent d2_data
      .def_property("d2", &TaoD1DataProxy::d2, &TaoD1DataProxy::set_d2)
      // TaoD1DataProxy.d (1D_PTR_type - Pointer to the appropriate section in u%data
      .def_property_readonly("d", &TaoD1DataProxy::d)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoD1DataProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoD1DataProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD1DataProxy& self) {
            return TaoD1DataProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoD1DataProxy& self, py::dict& memo) {
            return TaoD1DataProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoD1DataProxyArray1D>(m, "TaoD1DataStructArray1D");
  bind_FTypeAlloc1D<TaoD1DataProxyAlloc1D>(m, "TaoD1DataStructAlloc1D");
  // 2D TaoD1DataProxy arrays are not used in structs/routines
  // 3D TaoD1DataProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_d2_data_struct
void init_tao_d2_data_struct(py::module& m, py::class_<TaoD2DataProxy>& cls) {
  cls.def(py::init<>())
      // TaoD2DataProxy.name (0D_NOT_character - Name to be used with commands.
      .def_property("name", &TaoD2DataProxy::name, &TaoD2DataProxy::set_name)
      // TaoD2DataProxy.data_file_name (0D_NOT_character - Data file name .
      .def_property(
          "data_file_name",
          &TaoD2DataProxy::data_file_name,
          &TaoD2DataProxy::set_data_file_name)
      // TaoD2DataProxy.ref_file_name (0D_NOT_character - Reference file name.
      .def_property(
          "ref_file_name",
          &TaoD2DataProxy::ref_file_name,
          &TaoD2DataProxy::set_ref_file_name)
      // TaoD2DataProxy.data_date (0D_NOT_character - Data measurement date.
      .def_property(
          "data_date",
          &TaoD2DataProxy::data_date,
          &TaoD2DataProxy::set_data_date)
      // TaoD2DataProxy.ref_date (0D_NOT_character - Reference data measurement date.
      .def_property(
          "ref_date", &TaoD2DataProxy::ref_date, &TaoD2DataProxy::set_ref_date)
      // TaoD2DataProxy.descrip (1D_NOT_character - Array for descriptive information.
      .def_property_readonly("descrip", &TaoD2DataProxy::descrip)
      // TaoD2DataProxy.d1 (1D_ALLOC_type - Points to children
      .def_property_readonly("d1", &TaoD2DataProxy::d1)
      // TaoD2DataProxy.ix_universe (0D_NOT_integer - Index of universe this is in.
      .def_property(
          "ix_universe",
          &TaoD2DataProxy::ix_universe,
          &TaoD2DataProxy::set_ix_universe)
      // TaoD2DataProxy.ix_d2_data (0D_NOT_integer - Index in u%d2_data(:) array.
      .def_property(
          "ix_d2_data",
          &TaoD2DataProxy::ix_d2_data,
          &TaoD2DataProxy::set_ix_d2_data)
      // TaoD2DataProxy.ix_ref (0D_NOT_integer - Index of the reference data set.
      .def_property(
          "ix_ref", &TaoD2DataProxy::ix_ref, &TaoD2DataProxy::set_ix_ref)
      // TaoD2DataProxy.data_read_in (0D_NOT_logical - A data set has been read in?
      .def_property(
          "data_read_in",
          &TaoD2DataProxy::data_read_in,
          &TaoD2DataProxy::set_data_read_in)
      // TaoD2DataProxy.ref_read_in (0D_NOT_logical - A reference data set has been read in?
      .def_property(
          "ref_read_in",
          &TaoD2DataProxy::ref_read_in,
          &TaoD2DataProxy::set_ref_read_in)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoD2DataProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoD2DataProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD2DataProxy& self) {
            return TaoD2DataProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoD2DataProxy& self, py::dict& memo) {
            return TaoD2DataProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoD2DataProxyArray1D>(m, "TaoD2DataStructArray1D");
  bind_FTypeAlloc1D<TaoD2DataProxyAlloc1D>(m, "TaoD2DataStructAlloc1D");
  // 2D TaoD2DataProxy arrays are not used in structs/routines
  // 3D TaoD2DataProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_data_struct
void init_tao_data_struct(py::module& m, py::class_<TaoDataProxy>& cls) {
  cls.def(py::init<>())
      // TaoDataProxy.ele_name (0D_NOT_character - Name of the lattice element where datum is evaluated.
      .def_property(
          "ele_name", &TaoDataProxy::ele_name, &TaoDataProxy::set_ele_name)
      // TaoDataProxy.ele_start_name (0D_NOT_character - Name of starting lattice element when there is a range
      .def_property(
          "ele_start_name",
          &TaoDataProxy::ele_start_name,
          &TaoDataProxy::set_ele_start_name)
      // TaoDataProxy.ele_ref_name (0D_NOT_character - Name of reference lattice element
      .def_property(
          "ele_ref_name",
          &TaoDataProxy::ele_ref_name,
          &TaoDataProxy::set_ele_ref_name)
      // TaoDataProxy.data_type (0D_ALLOC_character - Type of data: 'orbit.x', etc.
      .def_property(
          "data_type", &TaoDataProxy::data_type, &TaoDataProxy::set_data_type)
      // TaoDataProxy.merit_type (0D_NOT_character - Type of constraint: 'target', 'max', 'min', etc.
      .def_property(
          "merit_type",
          &TaoDataProxy::merit_type,
          &TaoDataProxy::set_merit_type)
      // TaoDataProxy.id (0D_NOT_character - Used by Tao extension code. Not used by Tao directly.
      .def_property("id", &TaoDataProxy::id, &TaoDataProxy::set_id)
      // TaoDataProxy.data_source (0D_NOT_character - 'lat', 'beam', 'data' or 'var'. Last two used for expressions.
      .def_property(
          "data_source",
          &TaoDataProxy::data_source,
          &TaoDataProxy::set_data_source)
      // TaoDataProxy.why_invalid (0D_NOT_character - Informational string if there is a problem.
      .def_property(
          "why_invalid",
          &TaoDataProxy::why_invalid,
          &TaoDataProxy::set_why_invalid)
      // TaoDataProxy.ix_uni (0D_NOT_integer - Universe index of datum.
      .def_property("ix_uni", &TaoDataProxy::ix_uni, &TaoDataProxy::set_ix_uni)
      // TaoDataProxy.ix_bunch (0D_NOT_integer - Bunch number to get the data from.
      .def_property(
          "ix_bunch", &TaoDataProxy::ix_bunch, &TaoDataProxy::set_ix_bunch)
      // TaoDataProxy.ix_branch (0D_NOT_integer - Index of the associated lattice branch.
      .def_property(
          "ix_branch", &TaoDataProxy::ix_branch, &TaoDataProxy::set_ix_branch)
      // TaoDataProxy.ix_ele (0D_NOT_integer - Index of the lattice element corresponding to ele_name
      .def_property("ix_ele", &TaoDataProxy::ix_ele, &TaoDataProxy::set_ix_ele)
      // TaoDataProxy.ix_ele_start (0D_NOT_integer - Index of lattice elment when there is a range
      .def_property(
          "ix_ele_start",
          &TaoDataProxy::ix_ele_start,
          &TaoDataProxy::set_ix_ele_start)
      // TaoDataProxy.ix_ele_ref (0D_NOT_integer - Index of lattice elment when there is a reference.
      .def_property(
          "ix_ele_ref",
          &TaoDataProxy::ix_ele_ref,
          &TaoDataProxy::set_ix_ele_ref)
      // TaoDataProxy.ix_ele_merit (0D_NOT_integer - Index of lattice elment where merit is evaluated.
      .def_property(
          "ix_ele_merit",
          &TaoDataProxy::ix_ele_merit,
          &TaoDataProxy::set_ix_ele_merit)
      // TaoDataProxy.ix_d1 (0D_NOT_integer - Index number in u%d2_data(i)%d1_data(j)%d(:) array.
      .def_property("ix_d1", &TaoDataProxy::ix_d1, &TaoDataProxy::set_ix_d1)
      // TaoDataProxy.ix_data (0D_NOT_integer - Index of this datum in the u%data(:) array of data_structs.
      .def_property(
          "ix_data", &TaoDataProxy::ix_data, &TaoDataProxy::set_ix_data)
      // TaoDataProxy.ix_dModel (0D_NOT_integer - Row number in the dModel_dVar derivative matrix.
      .def_property(
          "ix_dModel", &TaoDataProxy::ix_dModel, &TaoDataProxy::set_ix_dModel)
      // TaoDataProxy.eval_point (0D_NOT_integer - or anchor_center$, anchor_beginning$. Where to evaluate data relative to the element.
      .def_property(
          "eval_point",
          &TaoDataProxy::eval_point,
          &TaoDataProxy::set_eval_point)
      // TaoDataProxy.meas_value (0D_NOT_real - Measured datum value.
      .def_property(
          "meas_value",
          &TaoDataProxy::meas_value,
          &TaoDataProxy::set_meas_value)
      // TaoDataProxy.ref_value (0D_NOT_real - Measured datum value from the reference data set.
      .def_property(
          "ref_value", &TaoDataProxy::ref_value, &TaoDataProxy::set_ref_value)
      // TaoDataProxy.model_value (0D_NOT_real - Datum value as calculated from the model.
      .def_property(
          "model_value",
          &TaoDataProxy::model_value,
          &TaoDataProxy::set_model_value)
      // TaoDataProxy.design_value (0D_NOT_real - What the datum value is in the design lattice.
      .def_property(
          "design_value",
          &TaoDataProxy::design_value,
          &TaoDataProxy::set_design_value)
      // TaoDataProxy.old_value (0D_NOT_real - The model_value at some previous time.
      .def_property(
          "old_value", &TaoDataProxy::old_value, &TaoDataProxy::set_old_value)
      // TaoDataProxy.base_value (0D_NOT_real - The value as calculated from the base model.
      .def_property(
          "base_value",
          &TaoDataProxy::base_value,
          &TaoDataProxy::set_base_value)
      // TaoDataProxy.error_rms (0D_NOT_real - Measurement error RMS. Used in plotting.
      .def_property(
          "error_rms", &TaoDataProxy::error_rms, &TaoDataProxy::set_error_rms)
      // TaoDataProxy.delta_merit (0D_NOT_real - Diff used to calculate the merit function term.
      .def_property(
          "delta_merit",
          &TaoDataProxy::delta_merit,
          &TaoDataProxy::set_delta_merit)
      // TaoDataProxy.weight (0D_NOT_real - Weight for the merit function term.
      .def_property("weight", &TaoDataProxy::weight, &TaoDataProxy::set_weight)
      // TaoDataProxy.invalid_value (0D_NOT_real - Value used in merit calc if good_model = F (or possibly good_design & good_base).
      .def_property(
          "invalid_value",
          &TaoDataProxy::invalid_value,
          &TaoDataProxy::set_invalid_value)
      // TaoDataProxy.merit (0D_NOT_real - Merit function term value: weight * delta_merit^2
      .def_property("merit", &TaoDataProxy::merit, &TaoDataProxy::set_merit)
      // TaoDataProxy.s (0D_NOT_real - longitudinal position of ele.
      .def_property("s", &TaoDataProxy::s, &TaoDataProxy::set_s)
      // TaoDataProxy.s_offset (0D_NOT_real - Offset of the evaluation point.
      .def_property(
          "s_offset", &TaoDataProxy::s_offset, &TaoDataProxy::set_s_offset)
      // TaoDataProxy.ref_s_offset (0D_NOT_real - Offset of the reference point. In development.
      .def_property(
          "ref_s_offset",
          &TaoDataProxy::ref_s_offset,
          &TaoDataProxy::set_ref_s_offset)
      // TaoDataProxy.err_message_printed (0D_NOT_logical - Used to prevent zillions of error messages being generated
      .def_property(
          "err_message_printed",
          &TaoDataProxy::err_message_printed,
          &TaoDataProxy::set_err_message_printed)
      // TaoDataProxy.exists (0D_NOT_logical - See above
      .def_property("exists", &TaoDataProxy::exists, &TaoDataProxy::set_exists)
      // TaoDataProxy.good_model (0D_NOT_logical - See above
      .def_property(
          "good_model",
          &TaoDataProxy::good_model,
          &TaoDataProxy::set_good_model)
      // TaoDataProxy.good_base (0D_NOT_logical - See above
      .def_property(
          "good_base", &TaoDataProxy::good_base, &TaoDataProxy::set_good_base)
      // TaoDataProxy.good_design (0D_NOT_logical - See above
      .def_property(
          "good_design",
          &TaoDataProxy::good_design,
          &TaoDataProxy::set_good_design)
      // TaoDataProxy.good_meas (0D_NOT_logical - See above
      .def_property(
          "good_meas", &TaoDataProxy::good_meas, &TaoDataProxy::set_good_meas)
      // TaoDataProxy.good_ref (0D_NOT_logical - See above
      .def_property(
          "good_ref", &TaoDataProxy::good_ref, &TaoDataProxy::set_good_ref)
      // TaoDataProxy.good_user (0D_NOT_logical - See above
      .def_property(
          "good_user", &TaoDataProxy::good_user, &TaoDataProxy::set_good_user)
      // TaoDataProxy.good_opt (0D_NOT_logical - See above
      .def_property(
          "good_opt", &TaoDataProxy::good_opt, &TaoDataProxy::set_good_opt)
      // TaoDataProxy.good_plot (0D_NOT_logical - See above
      .def_property(
          "good_plot", &TaoDataProxy::good_plot, &TaoDataProxy::set_good_plot)
      // TaoDataProxy.useit_plot (0D_NOT_logical - See above
      .def_property(
          "useit_plot",
          &TaoDataProxy::useit_plot,
          &TaoDataProxy::set_useit_plot)
      // TaoDataProxy.useit_opt (0D_NOT_logical - See above
      .def_property(
          "useit_opt", &TaoDataProxy::useit_opt, &TaoDataProxy::set_useit_opt)
      // TaoDataProxy.spin_map (0D_NOT_type -
      .def_property(
          "spin_map", &TaoDataProxy::spin_map, &TaoDataProxy::set_spin_map)
      // TaoDataProxy.d1 (0D_PTR_type - Pointer to the parent d1_data_struct
      .def_property("d1", &TaoDataProxy::d1, &TaoDataProxy::set_d1)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoDataProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const TaoDataProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDataProxy& self) {
            return TaoDataProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoDataProxy& self, py::dict& memo) {
            return TaoDataProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoDataProxyArray1D>(m, "TaoDataStructArray1D");
  bind_FTypeAlloc1D<TaoDataProxyAlloc1D>(m, "TaoDataStructAlloc1D");
  // 2D TaoDataProxy arrays are not used in structs/routines
  // 3D TaoDataProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_data_var_component_struct
void init_tao_data_var_component_struct(
    py::module& m,
    py::class_<TaoDataVarComponentProxy>& cls) {
  cls.def(py::init<>())
      // TaoDataVarComponentProxy.name (0D_NOT_character - Eg: 'meas', 'ref', 'model', etc.
      .def_property(
          "name",
          &TaoDataVarComponentProxy::name,
          &TaoDataVarComponentProxy::set_name)
      // TaoDataVarComponentProxy.sign (0D_NOT_real - +1 or -1
      .def_property(
          "sign",
          &TaoDataVarComponentProxy::sign,
          &TaoDataVarComponentProxy::set_sign)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoDataVarComponentProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoDataVarComponentProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDataVarComponentProxy& self) {
            return TaoDataVarComponentProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoDataVarComponentProxy& self, py::dict& memo) {
            return TaoDataVarComponentProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoDataVarComponentProxyArray1D>(
      m, "TaoDataVarComponentStructArray1D");
  bind_FTypeAlloc1D<TaoDataVarComponentProxyAlloc1D>(
      m, "TaoDataVarComponentStructAlloc1D");
  // 2D TaoDataVarComponentProxy arrays are not used in structs/routines
  // 3D TaoDataVarComponentProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_drawing_struct
void init_tao_drawing_struct(py::module& m, py::class_<TaoDrawingProxy>& cls) {
  cls.def(py::init<>())
      // TaoDrawingProxy.ele_shape (1D_ALLOC_type -
      .def_property_readonly("ele_shape", &TaoDrawingProxy::ele_shape)

      .def(
          "__repr__",
          [](const TaoDrawingProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDrawingProxy& self) {
            return TaoDrawingProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoDrawingProxy& self, py::dict& memo) {
            return TaoDrawingProxy(self);
          })

      ;

  // 1D TaoDrawingProxy arrays are not used in structs/routines
  // 2D TaoDrawingProxy arrays are not used in structs/routines
  // 3D TaoDrawingProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_dynamic_aperture_struct
void init_tao_dynamic_aperture_struct(
    py::module& m,
    py::class_<TaoDynamicApertureProxy>& cls) {
  cls.def(py::init<>())
      // TaoDynamicApertureProxy.param (0D_NOT_type -
      .def_property(
          "param",
          &TaoDynamicApertureProxy::param,
          &TaoDynamicApertureProxy::set_param)
      // TaoDynamicApertureProxy.scan (1D_ALLOC_type - One scan for each pz.
      .def_property_readonly("scan", &TaoDynamicApertureProxy::scan)
      // TaoDynamicApertureProxy.pz (1D_ALLOC_real -
      .def_property_readonly("pz", &TaoDynamicApertureProxy::pz)
      // TaoDynamicApertureProxy.ellipse_scale (0D_NOT_real -
      .def_property(
          "ellipse_scale",
          &TaoDynamicApertureProxy::ellipse_scale,
          &TaoDynamicApertureProxy::set_ellipse_scale)
      // TaoDynamicApertureProxy.a_emit (0D_NOT_real -
      .def_property(
          "a_emit",
          &TaoDynamicApertureProxy::a_emit,
          &TaoDynamicApertureProxy::set_a_emit)
      // TaoDynamicApertureProxy.b_emit (0D_NOT_real -
      .def_property(
          "b_emit",
          &TaoDynamicApertureProxy::b_emit,
          &TaoDynamicApertureProxy::set_b_emit)

      .def(
          "__repr__",
          [](const TaoDynamicApertureProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDynamicApertureProxy& self) {
            return TaoDynamicApertureProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoDynamicApertureProxy& self, py::dict& memo) {
            return TaoDynamicApertureProxy(self);
          })

      ;

  // 1D TaoDynamicApertureProxy arrays are not used in structs/routines
  // 2D TaoDynamicApertureProxy arrays are not used in structs/routines
  // 3D TaoDynamicApertureProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_pointer_struct
void init_tao_ele_pointer_struct(
    py::module& m,
    py::class_<TaoElePointerProxy>& cls) {
  cls.def(py::init<>())
      // TaoElePointerProxy.eles (1D_ALLOC_type -
      .def_property_readonly("eles", &TaoElePointerProxy::eles)
      // TaoElePointerProxy.n_loc (0D_NOT_integer -
      .def_property(
          "n_loc", &TaoElePointerProxy::n_loc, &TaoElePointerProxy::set_n_loc)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoElePointerProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoElePointerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoElePointerProxy& self) {
            return TaoElePointerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoElePointerProxy& self, py::dict& memo) {
            return TaoElePointerProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoElePointerProxyArray1D>(m, "TaoElePointerStructArray1D");
  bind_FTypeAlloc1D<TaoElePointerProxyAlloc1D>(m, "TaoElePointerStructAlloc1D");
  // 2D TaoElePointerProxy arrays are not used in structs/routines
  // 3D TaoElePointerProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_shape_struct
void init_tao_ele_shape_struct(
    py::module& m,
    py::class_<TaoEleShapeProxy>& cls) {
  cls.def(py::init<>())
      // TaoEleShapeProxy.ele_id (0D_NOT_character - element 'key::name' to match to.
      .def_property(
          "ele_id", &TaoEleShapeProxy::ele_id, &TaoEleShapeProxy::set_ele_id)
      // TaoEleShapeProxy.shape (0D_NOT_character - Shape to draw
      .def_property(
          "shape", &TaoEleShapeProxy::shape, &TaoEleShapeProxy::set_shape)
      // TaoEleShapeProxy.color (0D_NOT_character - Color of shape
      .def_property(
          "color", &TaoEleShapeProxy::color, &TaoEleShapeProxy::set_color)
      // TaoEleShapeProxy.size (0D_NOT_real - plot vertical height
      .def_property(
          "size", &TaoEleShapeProxy::size, &TaoEleShapeProxy::set_size)
      // TaoEleShapeProxy.label (0D_NOT_character - Can be: 'name', 's', 'none'
      .def_property(
          "label", &TaoEleShapeProxy::label, &TaoEleShapeProxy::set_label)
      // TaoEleShapeProxy.draw (0D_NOT_logical - Draw the shape?
      .def_property(
          "draw", &TaoEleShapeProxy::draw, &TaoEleShapeProxy::set_draw)
      // TaoEleShapeProxy.multi (0D_NOT_logical - Can be part of a multi-shape.
      .def_property(
          "multi", &TaoEleShapeProxy::multi, &TaoEleShapeProxy::set_multi)
      // TaoEleShapeProxy.line_width (0D_NOT_integer - Width of lines used to draw the shape.
      .def_property(
          "line_width",
          &TaoEleShapeProxy::line_width,
          &TaoEleShapeProxy::set_line_width)
      // TaoEleShapeProxy.offset (0D_NOT_real - Vertical offset.
      .def_property(
          "offset", &TaoEleShapeProxy::offset, &TaoEleShapeProxy::set_offset)
      // TaoEleShapeProxy.ix_key (0D_NOT_integer - Extracted from ele_id. 0 => all classes (quadrupole, etc.)
      .def_property(
          "ix_key", &TaoEleShapeProxy::ix_key, &TaoEleShapeProxy::set_ix_key)
      // TaoEleShapeProxy.name_ele (0D_NOT_character - Name of element.
      .def_property(
          "name_ele",
          &TaoEleShapeProxy::name_ele,
          &TaoEleShapeProxy::set_name_ele)
      // TaoEleShapeProxy.uni (1D_ALLOC_type -
      .def_property_readonly("uni", &TaoEleShapeProxy::uni)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoEleShapeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoEleShapeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoEleShapeProxy& self) {
            return TaoEleShapeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoEleShapeProxy& self, py::dict& memo) {
            return TaoEleShapeProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoEleShapeProxyArray1D>(m, "TaoEleShapeStructArray1D");
  bind_FTypeAlloc1D<TaoEleShapeProxyAlloc1D>(m, "TaoEleShapeStructAlloc1D");
  // 2D TaoEleShapeProxy arrays are not used in structs/routines
  // 3D TaoEleShapeProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_floor_plan_struct
void init_tao_floor_plan_struct(
    py::module& m,
    py::class_<TaoFloorPlanProxy>& cls) {
  cls.def(py::init<>())
      // TaoFloorPlanProxy.view (0D_NOT_character - or 'xz'.
      .def_property(
          "view", &TaoFloorPlanProxy::view, &TaoFloorPlanProxy::set_view)
      // TaoFloorPlanProxy.rotation (0D_NOT_real - Rotation of floor plan plot: 1.0 -> 360^deg
      .def_property(
          "rotation",
          &TaoFloorPlanProxy::rotation,
          &TaoFloorPlanProxy::set_rotation)
      // TaoFloorPlanProxy.correct_distortion (0D_NOT_logical - T -> Shrink one axis so x-scale = y-scale.
      .def_property(
          "correct_distortion",
          &TaoFloorPlanProxy::correct_distortion,
          &TaoFloorPlanProxy::set_correct_distortion)
      // TaoFloorPlanProxy.flip_label_side (0D_NOT_logical - Draw element label on other side of element?
      .def_property(
          "flip_label_side",
          &TaoFloorPlanProxy::flip_label_side,
          &TaoFloorPlanProxy::set_flip_label_side)
      // TaoFloorPlanProxy.size_is_absolute (0D_NOT_logical - Are shape sizes in meters or window pixels?
      .def_property(
          "size_is_absolute",
          &TaoFloorPlanProxy::size_is_absolute,
          &TaoFloorPlanProxy::set_size_is_absolute)
      // TaoFloorPlanProxy.draw_only_first_pass (0D_NOT_logical - Draw only first pass with multipass elements?
      .def_property(
          "draw_only_first_pass",
          &TaoFloorPlanProxy::draw_only_first_pass,
          &TaoFloorPlanProxy::set_draw_only_first_pass)
      // TaoFloorPlanProxy.draw_building_wall (0D_NOT_logical - Draw the building wall?
      .def_property(
          "draw_building_wall",
          &TaoFloorPlanProxy::draw_building_wall,
          &TaoFloorPlanProxy::set_draw_building_wall)
      // TaoFloorPlanProxy.orbit_scale (0D_NOT_real - Scale factor for drawing orbits. 0 -> Do not draw.
      .def_property(
          "orbit_scale",
          &TaoFloorPlanProxy::orbit_scale,
          &TaoFloorPlanProxy::set_orbit_scale)
      // TaoFloorPlanProxy.orbit_color (0D_NOT_character -
      .def_property(
          "orbit_color",
          &TaoFloorPlanProxy::orbit_color,
          &TaoFloorPlanProxy::set_orbit_color)
      // TaoFloorPlanProxy.orbit_pattern (0D_NOT_character -
      .def_property(
          "orbit_pattern",
          &TaoFloorPlanProxy::orbit_pattern,
          &TaoFloorPlanProxy::set_orbit_pattern)
      // TaoFloorPlanProxy.orbit_lattice (0D_NOT_character - Or 'design' or 'base'
      .def_property(
          "orbit_lattice",
          &TaoFloorPlanProxy::orbit_lattice,
          &TaoFloorPlanProxy::set_orbit_lattice)
      // TaoFloorPlanProxy.orbit_width (0D_NOT_integer -
      .def_property(
          "orbit_width",
          &TaoFloorPlanProxy::orbit_width,
          &TaoFloorPlanProxy::set_orbit_width)

      .def(
          "__repr__",
          [](const TaoFloorPlanProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoFloorPlanProxy& self) {
            return TaoFloorPlanProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoFloorPlanProxy& self, py::dict& memo) {
            return TaoFloorPlanProxy(self);
          })

      ;

  // 1D TaoFloorPlanProxy arrays are not used in structs/routines
  // 2D TaoFloorPlanProxy arrays are not used in structs/routines
  // 3D TaoFloorPlanProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_global_struct
void init_tao_global_struct(py::module& m, py::class_<TaoGlobalProxy>& cls) {
  cls.def(py::init<>())
      // TaoGlobalProxy.beam_dead_cutoff (0D_NOT_real - Percentage of dead particles at which beam tracking is stopped.
      .def_property(
          "beam_dead_cutoff",
          &TaoGlobalProxy::beam_dead_cutoff,
          &TaoGlobalProxy::set_beam_dead_cutoff)
      // TaoGlobalProxy.lm_opt_deriv_reinit (0D_NOT_real - Reinit derivative matrix cutoff
      .def_property(
          "lm_opt_deriv_reinit",
          &TaoGlobalProxy::lm_opt_deriv_reinit,
          &TaoGlobalProxy::set_lm_opt_deriv_reinit)
      // TaoGlobalProxy.de_lm_step_ratio (0D_NOT_real - Scaling for step sizes between DE and LM optimizers.
      .def_property(
          "de_lm_step_ratio",
          &TaoGlobalProxy::de_lm_step_ratio,
          &TaoGlobalProxy::set_de_lm_step_ratio)
      // TaoGlobalProxy.de_var_to_population_factor (0D_NOT_real - DE population = max(n_var*factor, 20)
      .def_property(
          "de_var_to_population_factor",
          &TaoGlobalProxy::de_var_to_population_factor,
          &TaoGlobalProxy::set_de_var_to_population_factor)
      // TaoGlobalProxy.lmdif_eps (0D_NOT_real - Tollerance for lmdif optimizer.
      .def_property(
          "lmdif_eps",
          &TaoGlobalProxy::lmdif_eps,
          &TaoGlobalProxy::set_lmdif_eps)
      // TaoGlobalProxy.lmdif_negligible_merit (0D_NOT_real -
      .def_property(
          "lmdif_negligible_merit",
          &TaoGlobalProxy::lmdif_negligible_merit,
          &TaoGlobalProxy::set_lmdif_negligible_merit)
      // TaoGlobalProxy.svd_cutoff (0D_NOT_real - SVD singular value cutoff.
      .def_property(
          "svd_cutoff",
          &TaoGlobalProxy::svd_cutoff,
          &TaoGlobalProxy::set_svd_cutoff)
      // TaoGlobalProxy.unstable_penalty (0D_NOT_real - Used in unstable_ring datum merit calculation.
      .def_property(
          "unstable_penalty",
          &TaoGlobalProxy::unstable_penalty,
          &TaoGlobalProxy::set_unstable_penalty)
      // TaoGlobalProxy.merit_stop_value (0D_NOT_real - Merit value below which an optimizer will stop.
      .def_property(
          "merit_stop_value",
          &TaoGlobalProxy::merit_stop_value,
          &TaoGlobalProxy::set_merit_stop_value)
      // TaoGlobalProxy.dmerit_stop_value (0D_NOT_real - Fractional Merit change below which an optimizer will stop.
      .def_property(
          "dmerit_stop_value",
          &TaoGlobalProxy::dmerit_stop_value,
          &TaoGlobalProxy::set_dmerit_stop_value)
      // TaoGlobalProxy.random_sigma_cutoff (0D_NOT_real - Cut-off in sigmas.
      .def_property(
          "random_sigma_cutoff",
          &TaoGlobalProxy::random_sigma_cutoff,
          &TaoGlobalProxy::set_random_sigma_cutoff)
      // TaoGlobalProxy.delta_e_chrom (0D_NOT_real - Delta E used from chrom calc.
      .def_property(
          "delta_e_chrom",
          &TaoGlobalProxy::delta_e_chrom,
          &TaoGlobalProxy::set_delta_e_chrom)
      // TaoGlobalProxy.max_plot_time (0D_NOT_real - If plotting time (seconds) exceeds this than a message is generated.
      .def_property(
          "max_plot_time",
          &TaoGlobalProxy::max_plot_time,
          &TaoGlobalProxy::set_max_plot_time)
      // TaoGlobalProxy.default_universe (0D_NOT_integer - Default universe to work with.
      .def_property(
          "default_universe",
          &TaoGlobalProxy::default_universe,
          &TaoGlobalProxy::set_default_universe)
      // TaoGlobalProxy.default_branch (0D_NOT_integer - Default lattice branch to work with.
      .def_property(
          "default_branch",
          &TaoGlobalProxy::default_branch,
          &TaoGlobalProxy::set_default_branch)
      // TaoGlobalProxy.n_opti_cycles (0D_NOT_integer - Number of optimization cycles
      .def_property(
          "n_opti_cycles",
          &TaoGlobalProxy::n_opti_cycles,
          &TaoGlobalProxy::set_n_opti_cycles)
      // TaoGlobalProxy.n_opti_loops (0D_NOT_integer - Number of optimization loops
      .def_property(
          "n_opti_loops",
          &TaoGlobalProxy::n_opti_loops,
          &TaoGlobalProxy::set_n_opti_loops)
      // TaoGlobalProxy.n_threads (0D_NOT_integer - Number of OpenMP threads for parallel calculations.
      .def_property(
          "n_threads",
          &TaoGlobalProxy::n_threads,
          &TaoGlobalProxy::set_n_threads)
      // TaoGlobalProxy.phase_units (0D_NOT_integer - Phase units on output.
      .def_property(
          "phase_units",
          &TaoGlobalProxy::phase_units,
          &TaoGlobalProxy::set_phase_units)
      // TaoGlobalProxy.bunch_to_plot (0D_NOT_integer - Which bunch to plot
      .def_property(
          "bunch_to_plot",
          &TaoGlobalProxy::bunch_to_plot,
          &TaoGlobalProxy::set_bunch_to_plot)
      // TaoGlobalProxy.random_seed (0D_NOT_integer - Use system clock by default
      .def_property(
          "random_seed",
          &TaoGlobalProxy::random_seed,
          &TaoGlobalProxy::set_random_seed)
      // TaoGlobalProxy.n_top10_merit (0D_NOT_integer - Number of top merit constraints to print.
      .def_property(
          "n_top10_merit",
          &TaoGlobalProxy::n_top10_merit,
          &TaoGlobalProxy::set_n_top10_merit)
      // TaoGlobalProxy.srdt_gen_n_slices (0D_NOT_integer - Number times to slice elements for summation RDT calculation
      .def_property(
          "srdt_gen_n_slices",
          &TaoGlobalProxy::srdt_gen_n_slices,
          &TaoGlobalProxy::set_srdt_gen_n_slices)
      // TaoGlobalProxy.datum_err_messages_max (0D_NOT_integer - Maximum number of error messages per call to lattice_calc.
      .def_property(
          "datum_err_messages_max",
          &TaoGlobalProxy::datum_err_messages_max,
          &TaoGlobalProxy::set_datum_err_messages_max)
      // TaoGlobalProxy.srdt_sxt_n_slices (0D_NOT_integer - Number times to slice sextupoles for summation RDT calculation
      .def_property(
          "srdt_sxt_n_slices",
          &TaoGlobalProxy::srdt_sxt_n_slices,
          &TaoGlobalProxy::set_srdt_sxt_n_slices)
      // TaoGlobalProxy.srdt_use_cache (0D_NOT_logical - Create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large.
      .def_property(
          "srdt_use_cache",
          &TaoGlobalProxy::srdt_use_cache,
          &TaoGlobalProxy::set_srdt_use_cache)
      // TaoGlobalProxy.quiet (0D_NOT_character - Print I/O when running a command file?
      .def_property("quiet", &TaoGlobalProxy::quiet, &TaoGlobalProxy::set_quiet)
      // TaoGlobalProxy.random_engine (0D_NOT_character - Non-beam random number engine
      .def_property(
          "random_engine",
          &TaoGlobalProxy::random_engine,
          &TaoGlobalProxy::set_random_engine)
      // TaoGlobalProxy.random_gauss_converter (0D_NOT_character - Non-beam
      .def_property(
          "random_gauss_converter",
          &TaoGlobalProxy::random_gauss_converter,
          &TaoGlobalProxy::set_random_gauss_converter)
      // TaoGlobalProxy.track_type (0D_NOT_character - or 'beam'
      .def_property(
          "track_type",
          &TaoGlobalProxy::track_type,
          &TaoGlobalProxy::set_track_type)
      // TaoGlobalProxy.lat_sigma_calc_uses_emit_from (0D_NOT_character - Lattice derived sigma matrix uses emit values from where? Other possibilities: 'beam', 'beam_init'.
      .def_property(
          "lat_sigma_calc_uses_emit_from",
          &TaoGlobalProxy::lat_sigma_calc_uses_emit_from,
          &TaoGlobalProxy::set_lat_sigma_calc_uses_emit_from)
      // TaoGlobalProxy.prompt_string (0D_NOT_character -
      .def_property(
          "prompt_string",
          &TaoGlobalProxy::prompt_string,
          &TaoGlobalProxy::set_prompt_string)
      // TaoGlobalProxy.prompt_color (0D_NOT_character - See read_a_line routine for possible settings.
      .def_property(
          "prompt_color",
          &TaoGlobalProxy::prompt_color,
          &TaoGlobalProxy::set_prompt_color)
      // TaoGlobalProxy.optimizer (0D_NOT_character - optimizer to use.
      .def_property(
          "optimizer",
          &TaoGlobalProxy::optimizer,
          &TaoGlobalProxy::set_optimizer)
      // TaoGlobalProxy.print_command (0D_NOT_character -
      .def_property(
          "print_command",
          &TaoGlobalProxy::print_command,
          &TaoGlobalProxy::set_print_command)
      // TaoGlobalProxy.var_out_file (0D_NOT_character -
      .def_property(
          "var_out_file",
          &TaoGlobalProxy::var_out_file,
          &TaoGlobalProxy::set_var_out_file)
      // TaoGlobalProxy.history_file (0D_NOT_character -
      .def_property(
          "history_file",
          &TaoGlobalProxy::history_file,
          &TaoGlobalProxy::set_history_file)
      // TaoGlobalProxy.beam_timer_on (0D_NOT_logical - For timing the beam tracking calculation.
      .def_property(
          "beam_timer_on",
          &TaoGlobalProxy::beam_timer_on,
          &TaoGlobalProxy::set_beam_timer_on)
      // TaoGlobalProxy.box_plots (0D_NOT_logical - For debugging plot layout issues.
      .def_property(
          "box_plots",
          &TaoGlobalProxy::box_plots,
          &TaoGlobalProxy::set_box_plots)
      // TaoGlobalProxy.blank_line_between_commands (0D_NOT_logical - Add a blank line between command output?
      .def_property(
          "blank_line_between_commands",
          &TaoGlobalProxy::blank_line_between_commands,
          &TaoGlobalProxy::set_blank_line_between_commands)
      // TaoGlobalProxy.cmd_file_abort_on_error (0D_NOT_logical - Abort open command files if there is an error?
      .def_property(
          "cmd_file_abort_on_error",
          &TaoGlobalProxy::cmd_file_abort_on_error,
          &TaoGlobalProxy::set_cmd_file_abort_on_error)
      // TaoGlobalProxy.concatenate_maps (0D_NOT_logical - False => tracking using DA.
      .def_property(
          "concatenate_maps",
          &TaoGlobalProxy::concatenate_maps,
          &TaoGlobalProxy::set_concatenate_maps)
      // TaoGlobalProxy.derivative_recalc (0D_NOT_logical - Recalc before each optimizer run?
      .def_property(
          "derivative_recalc",
          &TaoGlobalProxy::derivative_recalc,
          &TaoGlobalProxy::set_derivative_recalc)
      // TaoGlobalProxy.derivative_uses_design (0D_NOT_logical - Derivative calc uses design lattice instead of model?
      .def_property(
          "derivative_uses_design",
          &TaoGlobalProxy::derivative_uses_design,
          &TaoGlobalProxy::set_derivative_uses_design)
      // TaoGlobalProxy.disable_smooth_line_calc (0D_NOT_logical - Global disable of the smooth line calculation.
      .def_property(
          "disable_smooth_line_calc",
          &TaoGlobalProxy::disable_smooth_line_calc,
          &TaoGlobalProxy::set_disable_smooth_line_calc)
      // TaoGlobalProxy.draw_curve_off_scale_warn (0D_NOT_logical - Display warning on graphs?
      .def_property(
          "draw_curve_off_scale_warn",
          &TaoGlobalProxy::draw_curve_off_scale_warn,
          &TaoGlobalProxy::set_draw_curve_off_scale_warn)
      // TaoGlobalProxy.external_plotting (0D_NOT_logical - Used with matplotlib and gui.
      .def_property(
          "external_plotting",
          &TaoGlobalProxy::external_plotting,
          &TaoGlobalProxy::set_external_plotting)
      // TaoGlobalProxy.label_lattice_elements (0D_NOT_logical - For lat_layout plots
      .def_property(
          "label_lattice_elements",
          &TaoGlobalProxy::label_lattice_elements,
          &TaoGlobalProxy::set_label_lattice_elements)
      // TaoGlobalProxy.label_keys (0D_NOT_logical - For lat_layout plots
      .def_property(
          "label_keys",
          &TaoGlobalProxy::label_keys,
          &TaoGlobalProxy::set_label_keys)
      // TaoGlobalProxy.lattice_calc_on (0D_NOT_logical - Turn on/off beam and single particle calculations.
      .def_property(
          "lattice_calc_on",
          &TaoGlobalProxy::lattice_calc_on,
          &TaoGlobalProxy::set_lattice_calc_on)
      // TaoGlobalProxy.only_limit_opt_vars (0D_NOT_logical - Only apply limits to variables used in optimization.
      .def_property(
          "only_limit_opt_vars",
          &TaoGlobalProxy::only_limit_opt_vars,
          &TaoGlobalProxy::set_only_limit_opt_vars)
      // TaoGlobalProxy.opt_with_ref (0D_NOT_logical - Use reference data in optimization?
      .def_property(
          "opt_with_ref",
          &TaoGlobalProxy::opt_with_ref,
          &TaoGlobalProxy::set_opt_with_ref)
      // TaoGlobalProxy.opt_with_base (0D_NOT_logical - Use base data in optimization?
      .def_property(
          "opt_with_base",
          &TaoGlobalProxy::opt_with_base,
          &TaoGlobalProxy::set_opt_with_base)
      // TaoGlobalProxy.opt_match_auto_recalc (0D_NOT_logical - Set recalc = True for match elements before each cycle?
      .def_property(
          "opt_match_auto_recalc",
          &TaoGlobalProxy::opt_match_auto_recalc,
          &TaoGlobalProxy::set_opt_match_auto_recalc)
      // TaoGlobalProxy.opti_write_var_file (0D_NOT_logical - 'run' command writes var_out_file
      .def_property(
          "opti_write_var_file",
          &TaoGlobalProxy::opti_write_var_file,
          &TaoGlobalProxy::set_opti_write_var_file)
      // TaoGlobalProxy.optimizer_allow_user_abort (0D_NOT_logical - See Tao manual for more details.
      .def_property(
          "optimizer_allow_user_abort",
          &TaoGlobalProxy::optimizer_allow_user_abort,
          &TaoGlobalProxy::set_optimizer_allow_user_abort)
      // TaoGlobalProxy.optimizer_var_limit_warn (0D_NOT_logical - Warn when vars reach a limit with optimization.
      .def_property(
          "optimizer_var_limit_warn",
          &TaoGlobalProxy::optimizer_var_limit_warn,
          &TaoGlobalProxy::set_optimizer_var_limit_warn)
      // TaoGlobalProxy.plot_on (0D_NOT_logical - Do plotting?
      .def_property(
          "plot_on", &TaoGlobalProxy::plot_on, &TaoGlobalProxy::set_plot_on)
      // TaoGlobalProxy.rad_int_user_calc_on (0D_NOT_logical - User set radiation integrals calculation on/off.
      .def_property(
          "rad_int_user_calc_on",
          &TaoGlobalProxy::rad_int_user_calc_on,
          &TaoGlobalProxy::set_rad_int_user_calc_on)
      // TaoGlobalProxy.rf_on (0D_NOT_logical - RFcavities on or off? Does not affect lcavities.
      .def_property("rf_on", &TaoGlobalProxy::rf_on, &TaoGlobalProxy::set_rf_on)
      // TaoGlobalProxy.single_step (0D_NOT_logical - For debugging and demonstrations: Single step through a command file?
      .def_property(
          "single_step",
          &TaoGlobalProxy::single_step,
          &TaoGlobalProxy::set_single_step)
      // TaoGlobalProxy.stop_on_error (0D_NOT_logical - For debugging: False prevents tao from exiting on an error.
      .def_property(
          "stop_on_error",
          &TaoGlobalProxy::stop_on_error,
          &TaoGlobalProxy::set_stop_on_error)
      // TaoGlobalProxy.svd_retreat_on_merit_increase (0D_NOT_logical -
      .def_property(
          "svd_retreat_on_merit_increase",
          &TaoGlobalProxy::svd_retreat_on_merit_increase,
          &TaoGlobalProxy::set_svd_retreat_on_merit_increase)
      // TaoGlobalProxy.var_limits_on (0D_NOT_logical - Respect the variable limits?
      .def_property(
          "var_limits_on",
          &TaoGlobalProxy::var_limits_on,
          &TaoGlobalProxy::set_var_limits_on)
      // TaoGlobalProxy.wait_for_CR_in_single_mode (0D_NOT_logical - For use with a python GUI.
      .def_property(
          "wait_for_CR_in_single_mode",
          &TaoGlobalProxy::wait_for_CR_in_single_mode,
          &TaoGlobalProxy::set_wait_for_CR_in_single_mode)
      // TaoGlobalProxy.symbol_import (0D_NOT_logical - Import symbols from lattice file(s)? Internal stuff
      .def_property(
          "symbol_import",
          &TaoGlobalProxy::symbol_import,
          &TaoGlobalProxy::set_symbol_import)
      // TaoGlobalProxy.debug_on (0D_NOT_logical - For debugging.
      .def_property(
          "debug_on", &TaoGlobalProxy::debug_on, &TaoGlobalProxy::set_debug_on)
      // TaoGlobalProxy.expression_tree_on (0D_NOT_logical - Use an expression tree instead of a stack?
      .def_property(
          "expression_tree_on",
          &TaoGlobalProxy::expression_tree_on,
          &TaoGlobalProxy::set_expression_tree_on)
      // TaoGlobalProxy.verbose_on (0D_NOT_logical - For verbose output. Used with debugging.
      .def_property(
          "verbose_on",
          &TaoGlobalProxy::verbose_on,
          &TaoGlobalProxy::set_verbose_on)

      .def(
          "__repr__",
          [](const TaoGlobalProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoGlobalProxy& self) {
            return TaoGlobalProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoGlobalProxy& self, py::dict& memo) {
            return TaoGlobalProxy(self);
          })

      ;

  // 1D TaoGlobalProxy arrays are not used in structs/routines
  // 2D TaoGlobalProxy arrays are not used in structs/routines
  // 3D TaoGlobalProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_graph_struct
void init_tao_graph_struct(py::module& m, py::class_<TaoGraphProxy>& cls) {
  cls.def(py::init<>())
      // TaoGraphProxy.name (0D_NOT_character - Name identifying the graph
      .def_property("name", &TaoGraphProxy::name, &TaoGraphProxy::set_name)
      // TaoGraphProxy.type (0D_NOT_character - 'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture'
      .def_property("type", &TaoGraphProxy::type, &TaoGraphProxy::set_type)
      // TaoGraphProxy.title (0D_NOT_character -
      .def_property("title", &TaoGraphProxy::title, &TaoGraphProxy::set_title)
      // TaoGraphProxy.title_suffix (0D_NOT_character -
      .def_property(
          "title_suffix",
          &TaoGraphProxy::title_suffix,
          &TaoGraphProxy::set_title_suffix)
      // TaoGraphProxy.text_legend (1D_NOT_character - Array for holding descriptive info.
      .def_property_readonly("text_legend", &TaoGraphProxy::text_legend)
      // TaoGraphProxy.text_legend_out (1D_NOT_character - Array for holding descriptive info.
      .def_property_readonly("text_legend_out", &TaoGraphProxy::text_legend_out)
      // TaoGraphProxy.why_invalid (0D_NOT_character - Informative string to print.
      .def_property(
          "why_invalid",
          &TaoGraphProxy::why_invalid,
          &TaoGraphProxy::set_why_invalid)
      // TaoGraphProxy.curve (1D_ALLOC_type -
      .def_property_readonly("curve", &TaoGraphProxy::curve)
      // TaoGraphProxy.p (0D_PTR_type - pointer to parent plot
      .def_property("p", &TaoGraphProxy::p, &TaoGraphProxy::set_p)
      // TaoGraphProxy.floor_plan (0D_NOT_type -
      .def_property(
          "floor_plan",
          &TaoGraphProxy::floor_plan,
          &TaoGraphProxy::set_floor_plan)
      // TaoGraphProxy.text_legend_origin (0D_NOT_type -
      .def_property(
          "text_legend_origin",
          &TaoGraphProxy::text_legend_origin,
          &TaoGraphProxy::set_text_legend_origin)
      // TaoGraphProxy.curve_legend_origin (0D_NOT_type -
      .def_property(
          "curve_legend_origin",
          &TaoGraphProxy::curve_legend_origin,
          &TaoGraphProxy::set_curve_legend_origin)
      // TaoGraphProxy.curve_legend (0D_NOT_type -
      .def_property(
          "curve_legend",
          &TaoGraphProxy::curve_legend,
          &TaoGraphProxy::set_curve_legend)
      // TaoGraphProxy.x (0D_NOT_type - X-axis parameters.
      .def_property("x", &TaoGraphProxy::x, &TaoGraphProxy::set_x)
      // TaoGraphProxy.y (0D_NOT_type - Y-axis attributes.
      .def_property("y", &TaoGraphProxy::y, &TaoGraphProxy::set_y)
      // TaoGraphProxy.x2 (0D_NOT_type - X2-axis attributes (Not currently used).
      .def_property("x2", &TaoGraphProxy::x2, &TaoGraphProxy::set_x2)
      // TaoGraphProxy.y2 (0D_NOT_type - Y2-axis attributes.
      .def_property("y2", &TaoGraphProxy::y2, &TaoGraphProxy::set_y2)
      // TaoGraphProxy.margin (0D_NOT_type - Margin around the graph.
      .def_property(
          "margin", &TaoGraphProxy::margin, &TaoGraphProxy::set_margin)
      // TaoGraphProxy.scale_margin (0D_NOT_type - Margin for scaling
      .def_property(
          "scale_margin",
          &TaoGraphProxy::scale_margin,
          &TaoGraphProxy::set_scale_margin)
      // TaoGraphProxy.x_axis_scale_factor (0D_NOT_real - x-axis conversion from internal to plotting units.
      .def_property(
          "x_axis_scale_factor",
          &TaoGraphProxy::x_axis_scale_factor,
          &TaoGraphProxy::set_x_axis_scale_factor)
      // TaoGraphProxy.symbol_size_scale (0D_NOT_real - Symbol size scale factor for phase_space plots.
      .def_property(
          "symbol_size_scale",
          &TaoGraphProxy::symbol_size_scale,
          &TaoGraphProxy::set_symbol_size_scale)
      // TaoGraphProxy.box (1D_NOT_integer - Defines which box the plot is put in.
      .def_property_readonly("box", &TaoGraphProxy::box)
      // TaoGraphProxy.ix_branch (0D_NOT_integer - Branch in lattice. Used when there are no associated curves.
      .def_property(
          "ix_branch", &TaoGraphProxy::ix_branch, &TaoGraphProxy::set_ix_branch)
      // TaoGraphProxy.ix_universe (0D_NOT_integer - Used for lat_layout plots.
      .def_property(
          "ix_universe",
          &TaoGraphProxy::ix_universe,
          &TaoGraphProxy::set_ix_universe)
      // TaoGraphProxy.clip (0D_NOT_logical - Clip plot at graph boundary.
      .def_property("clip", &TaoGraphProxy::clip, &TaoGraphProxy::set_clip)
      // TaoGraphProxy.y2_mirrors_y (0D_NOT_logical - Y2-axis same as Y-axis?
      .def_property(
          "y2_mirrors_y",
          &TaoGraphProxy::y2_mirrors_y,
          &TaoGraphProxy::set_y2_mirrors_y)
      // TaoGraphProxy.limited (0D_NOT_logical - True if at least one data point past graph bounds.
      .def_property(
          "limited", &TaoGraphProxy::limited, &TaoGraphProxy::set_limited)
      // TaoGraphProxy.draw_axes (0D_NOT_logical - Draw axes, labels, etc?
      .def_property(
          "draw_axes", &TaoGraphProxy::draw_axes, &TaoGraphProxy::set_draw_axes)
      // TaoGraphProxy.draw_curve_legend (0D_NOT_logical - Legend for displaying curve info.
      .def_property(
          "draw_curve_legend",
          &TaoGraphProxy::draw_curve_legend,
          &TaoGraphProxy::set_draw_curve_legend)
      // TaoGraphProxy.draw_grid (0D_NOT_logical - Draw a grid?
      .def_property(
          "draw_grid", &TaoGraphProxy::draw_grid, &TaoGraphProxy::set_draw_grid)
      // TaoGraphProxy.draw_title (0D_NOT_logical -
      .def_property(
          "draw_title",
          &TaoGraphProxy::draw_title,
          &TaoGraphProxy::set_draw_title)
      // TaoGraphProxy.draw_only_good_user_data_or_vars (0D_NOT_logical -
      .def_property(
          "draw_only_good_user_data_or_vars",
          &TaoGraphProxy::draw_only_good_user_data_or_vars,
          &TaoGraphProxy::set_draw_only_good_user_data_or_vars)
      // TaoGraphProxy.allow_wrap_around (0D_NOT_logical - 'Wrap' curves to extend past lattice boundaries?
      .def_property(
          "allow_wrap_around",
          &TaoGraphProxy::allow_wrap_around,
          &TaoGraphProxy::set_allow_wrap_around)
      // TaoGraphProxy.is_valid (0D_NOT_logical - EG: Bad x_axis_type.
      .def_property(
          "is_valid", &TaoGraphProxy::is_valid, &TaoGraphProxy::set_is_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoGraphProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const TaoGraphProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoGraphProxy& self) {
            return TaoGraphProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoGraphProxy& self, py::dict& memo) {
            return TaoGraphProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoGraphProxyArray1D>(m, "TaoGraphStructArray1D");
  bind_FTypeAlloc1D<TaoGraphProxyAlloc1D>(m, "TaoGraphStructAlloc1D");
  // 2D TaoGraphProxy arrays are not used in structs/routines
  // 3D TaoGraphProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_histogram_struct
void init_tao_histogram_struct(
    py::module& m,
    py::class_<TaoHistogramProxy>& cls) {
  cls.def(py::init<>())
      // TaoHistogramProxy.density_normalized (0D_NOT_logical -
      .def_property(
          "density_normalized",
          &TaoHistogramProxy::density_normalized,
          &TaoHistogramProxy::set_density_normalized)
      // TaoHistogramProxy.weight_by_charge (0D_NOT_logical -
      .def_property(
          "weight_by_charge",
          &TaoHistogramProxy::weight_by_charge,
          &TaoHistogramProxy::set_weight_by_charge)
      // TaoHistogramProxy.minimum (0D_NOT_real - Computed by Tao. Not User settable.
      .def_property(
          "minimum",
          &TaoHistogramProxy::minimum,
          &TaoHistogramProxy::set_minimum)
      // TaoHistogramProxy.maximum (0D_NOT_real - Computed by Tao. Not User settable.
      .def_property(
          "maximum",
          &TaoHistogramProxy::maximum,
          &TaoHistogramProxy::set_maximum)
      // TaoHistogramProxy.width (0D_NOT_real -
      .def_property(
          "width", &TaoHistogramProxy::width, &TaoHistogramProxy::set_width)
      // TaoHistogramProxy.center (0D_NOT_real -
      .def_property(
          "center", &TaoHistogramProxy::center, &TaoHistogramProxy::set_center)
      // TaoHistogramProxy.number (0D_NOT_integer -
      .def_property(
          "number", &TaoHistogramProxy::number, &TaoHistogramProxy::set_number)

      .def(
          "__repr__",
          [](const TaoHistogramProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoHistogramProxy& self) {
            return TaoHistogramProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoHistogramProxy& self, py::dict& memo) {
            return TaoHistogramProxy(self);
          })

      ;

  // 1D TaoHistogramProxy arrays are not used in structs/routines
  // 2D TaoHistogramProxy arrays are not used in structs/routines
  // 3D TaoHistogramProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_init_struct
void init_tao_init_struct(py::module& m, py::class_<TaoInitProxy>& cls) {
  cls.def(py::init<>())
      // TaoInitProxy.parse_cmd_args (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property(
          "parse_cmd_args",
          &TaoInitProxy::parse_cmd_args,
          &TaoInitProxy::set_parse_cmd_args)
      // TaoInitProxy.debug_switch (0D_NOT_logical - Is the '-debug' switch present?
      .def_property(
          "debug_switch",
          &TaoInitProxy::debug_switch,
          &TaoInitProxy::set_debug_switch)
      // TaoInitProxy.external_plotting_switch (0D_NOT_logical - Is '-external_plotting' switch present?
      .def_property(
          "external_plotting_switch",
          &TaoInitProxy::external_plotting_switch,
          &TaoInitProxy::set_external_plotting_switch)
      // TaoInitProxy.init_name (0D_NOT_character - label for initialization
      .def_property(
          "init_name", &TaoInitProxy::init_name, &TaoInitProxy::set_init_name)
      // TaoInitProxy.hook_init_file (0D_NOT_character -
      .def_property(
          "hook_init_file",
          &TaoInitProxy::hook_init_file,
          &TaoInitProxy::set_hook_init_file)
      // TaoInitProxy.hook_lat_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_lat_file",
          &TaoInitProxy::hook_lat_file,
          &TaoInitProxy::set_hook_lat_file)
      // TaoInitProxy.hook_beam_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_beam_file",
          &TaoInitProxy::hook_beam_file,
          &TaoInitProxy::set_hook_beam_file)
      // TaoInitProxy.hook_data_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_data_file",
          &TaoInitProxy::hook_data_file,
          &TaoInitProxy::set_hook_data_file)
      // TaoInitProxy.hook_plot_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_plot_file",
          &TaoInitProxy::hook_plot_file,
          &TaoInitProxy::set_hook_plot_file)
      // TaoInitProxy.hook_startup_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_startup_file",
          &TaoInitProxy::hook_startup_file,
          &TaoInitProxy::set_hook_startup_file)
      // TaoInitProxy.hook_var_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_var_file",
          &TaoInitProxy::hook_var_file,
          &TaoInitProxy::set_hook_var_file)
      // TaoInitProxy.hook_building_wall_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_building_wall_file",
          &TaoInitProxy::hook_building_wall_file,
          &TaoInitProxy::set_hook_building_wall_file)
      // TaoInitProxy.init_file_arg_path (0D_NOT_character - Path part of init_tao_file
      .def_property(
          "init_file_arg_path",
          &TaoInitProxy::init_file_arg_path,
          &TaoInitProxy::set_init_file_arg_path)
      // TaoInitProxy.lattice_file_arg (0D_NOT_character - -lattice_file        command line argument.
      .def_property(
          "lattice_file_arg",
          &TaoInitProxy::lattice_file_arg,
          &TaoInitProxy::set_lattice_file_arg)
      // TaoInitProxy.hook_init_file_arg (0D_NOT_character - -hook_init_file      command line argument
      .def_property(
          "hook_init_file_arg",
          &TaoInitProxy::hook_init_file_arg,
          &TaoInitProxy::set_hook_init_file_arg)
      // TaoInitProxy.init_file_arg (0D_NOT_character - -init_file           command line argument.
      .def_property(
          "init_file_arg",
          &TaoInitProxy::init_file_arg,
          &TaoInitProxy::set_init_file_arg)
      // TaoInitProxy.beam_file_arg (0D_NOT_character - -beam_file           command line argument.
      .def_property(
          "beam_file_arg",
          &TaoInitProxy::beam_file_arg,
          &TaoInitProxy::set_beam_file_arg)
      // TaoInitProxy.beam_init_position_file_arg (0D_NOT_character - -beam_init_position_file command line argument.
      .def_property(
          "beam_init_position_file_arg",
          &TaoInitProxy::beam_init_position_file_arg,
          &TaoInitProxy::set_beam_init_position_file_arg)
      // TaoInitProxy.command_arg (0D_NOT_character - -command             command line argument.
      .def_property(
          "command_arg",
          &TaoInitProxy::command_arg,
          &TaoInitProxy::set_command_arg)
      // TaoInitProxy.data_file_arg (0D_NOT_character - -data_file           command line argument.
      .def_property(
          "data_file_arg",
          &TaoInitProxy::data_file_arg,
          &TaoInitProxy::set_data_file_arg)
      // TaoInitProxy.plot_file_arg (0D_NOT_character - -plot_file           command line argument.
      .def_property(
          "plot_file_arg",
          &TaoInitProxy::plot_file_arg,
          &TaoInitProxy::set_plot_file_arg)
      // TaoInitProxy.startup_file_arg (0D_NOT_character - -startup_file        command line argument.
      .def_property(
          "startup_file_arg",
          &TaoInitProxy::startup_file_arg,
          &TaoInitProxy::set_startup_file_arg)
      // TaoInitProxy.var_file_arg (0D_NOT_character - -var_file            command line argument.
      .def_property(
          "var_file_arg",
          &TaoInitProxy::var_file_arg,
          &TaoInitProxy::set_var_file_arg)
      // TaoInitProxy.building_wall_file_arg (0D_NOT_character - -building_wall_file  command line argument.
      .def_property(
          "building_wall_file_arg",
          &TaoInitProxy::building_wall_file_arg,
          &TaoInitProxy::set_building_wall_file_arg)
      // TaoInitProxy.geometry_arg (0D_NOT_character - -geometry            command line argument.
      .def_property(
          "geometry_arg",
          &TaoInitProxy::geometry_arg,
          &TaoInitProxy::set_geometry_arg)
      // TaoInitProxy.slice_lattice_arg (0D_NOT_character - -slice_lattice       command line argument.
      .def_property(
          "slice_lattice_arg",
          &TaoInitProxy::slice_lattice_arg,
          &TaoInitProxy::set_slice_lattice_arg)
      // TaoInitProxy.start_branch_at_arg (0D_NOT_character - -start_branch_at     command line argument.
      .def_property(
          "start_branch_at_arg",
          &TaoInitProxy::start_branch_at_arg,
          &TaoInitProxy::set_start_branch_at_arg)
      // TaoInitProxy.log_startup_arg (0D_NOT_character - -log_startup         command line argument
      .def_property(
          "log_startup_arg",
          &TaoInitProxy::log_startup_arg,
          &TaoInitProxy::set_log_startup_arg)
      // TaoInitProxy.no_stopping_arg (0D_NOT_character - -no_stopping         command line argument
      .def_property(
          "no_stopping_arg",
          &TaoInitProxy::no_stopping_arg,
          &TaoInitProxy::set_no_stopping_arg)
      // TaoInitProxy.noplot_arg (0D_NOT_character - -noplot              command line argument
      .def_property(
          "noplot_arg",
          &TaoInitProxy::noplot_arg,
          &TaoInitProxy::set_noplot_arg)
      // TaoInitProxy.no_rad_int_arg (0D_NOT_character - -no_rad_int          command line argument
      .def_property(
          "no_rad_int_arg",
          &TaoInitProxy::no_rad_int_arg,
          &TaoInitProxy::set_no_rad_int_arg)
      // TaoInitProxy.reverse_arg (0D_NOT_character - -reverse             command line argument
      .def_property(
          "reverse_arg",
          &TaoInitProxy::reverse_arg,
          &TaoInitProxy::set_reverse_arg)
      // TaoInitProxy.debug_arg (0D_NOT_character - -debug               command line argument
      .def_property(
          "debug_arg", &TaoInitProxy::debug_arg, &TaoInitProxy::set_debug_arg)
      // TaoInitProxy.disable_smooth_line_calc_arg (0D_NOT_character - -disable_smooth_line_calc
      .def_property(
          "disable_smooth_line_calc_arg",
          &TaoInitProxy::disable_smooth_line_calc_arg,
          &TaoInitProxy::set_disable_smooth_line_calc_arg)
      // TaoInitProxy.rf_on_arg (0D_NOT_character - -rf_on               command line argument
      .def_property(
          "rf_on_arg", &TaoInitProxy::rf_on_arg, &TaoInitProxy::set_rf_on_arg)
      // TaoInitProxy.prompt_color_arg (0D_NOT_character - -prompt_color        command line argument
      .def_property(
          "prompt_color_arg",
          &TaoInitProxy::prompt_color_arg,
          &TaoInitProxy::set_prompt_color_arg)
      // TaoInitProxy.quiet_arg (0D_NOT_character - -quiet               command line argument
      .def_property(
          "quiet_arg", &TaoInitProxy::quiet_arg, &TaoInitProxy::set_quiet_arg)
      // TaoInitProxy.noinit_arg (0D_NOT_character - -noinit              command line argument
      .def_property(
          "noinit_arg",
          &TaoInitProxy::noinit_arg,
          &TaoInitProxy::set_noinit_arg)
      // TaoInitProxy.nostartup_arg (0D_NOT_character - -nostartup           command line argument
      .def_property(
          "nostartup_arg",
          &TaoInitProxy::nostartup_arg,
          &TaoInitProxy::set_nostartup_arg)
      // TaoInitProxy.symbol_import_arg (0D_NOT_character - -symbol_import       command line argument
      .def_property(
          "symbol_import_arg",
          &TaoInitProxy::symbol_import_arg,
          &TaoInitProxy::set_symbol_import_arg)
      // TaoInitProxy.unique_name_suffix (0D_NOT_character -
      .def_property(
          "unique_name_suffix",
          &TaoInitProxy::unique_name_suffix,
          &TaoInitProxy::set_unique_name_suffix)

      .def("__repr__", [](const TaoInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoInitProxy& self) {
            return TaoInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoInitProxy& self, py::dict& memo) {
            return TaoInitProxy(self);
          })

      ;

  // 1D TaoInitProxy arrays are not used in structs/routines
  // 2D TaoInitProxy arrays are not used in structs/routines
  // 3D TaoInitProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_lat_sigma_struct
void init_tao_lat_sigma_struct(
    py::module& m,
    py::class_<TaoLatSigmaProxy>& cls) {
  cls.def(py::init<>())
      // TaoLatSigmaProxy.mat (2D_NOT_real -
      .def_property_readonly("mat", &TaoLatSigmaProxy::mat)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoLatSigmaProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoLatSigmaProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatSigmaProxy& self) {
            return TaoLatSigmaProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoLatSigmaProxy& self, py::dict& memo) {
            return TaoLatSigmaProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoLatSigmaProxyArray1D>(m, "TaoLatSigmaStructArray1D");
  bind_FTypeAlloc1D<TaoLatSigmaProxyAlloc1D>(m, "TaoLatSigmaStructAlloc1D");
  // 2D TaoLatSigmaProxy arrays are not used in structs/routines
  // 3D TaoLatSigmaProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_branch_struct
void init_tao_lattice_branch_struct(
    py::module& m,
    py::class_<TaoLatticeBranchProxy>& cls) {
  cls.def(py::init<>())
      // TaoLatticeBranchProxy.tao_lat (0D_PTR_type - Parent tao_lat
      .def_property(
          "tao_lat",
          &TaoLatticeBranchProxy::tao_lat,
          &TaoLatticeBranchProxy::set_tao_lat)
      // TaoLatticeBranchProxy.lat_sigma (1D_ALLOC_type - Sigma matrix derived from lattice (not beam).
      .def_property_readonly("lat_sigma", &TaoLatticeBranchProxy::lat_sigma)
      // TaoLatticeBranchProxy.spin_ele (1D_ALLOC_type - Spin stuff
      .def_property_readonly("spin_ele", &TaoLatticeBranchProxy::spin_ele)
      // TaoLatticeBranchProxy.bunch_params (1D_ALLOC_type - Per element
      .def_property_readonly(
          "bunch_params", &TaoLatticeBranchProxy::bunch_params)
      // TaoLatticeBranchProxy.bunch_params_comb (1D_ALLOC_type - A comb for each bunch in beam.
      .def_property_readonly(
          "bunch_params_comb", &TaoLatticeBranchProxy::bunch_params_comb)
      // TaoLatticeBranchProxy.orbit (1D_ALLOC_type -
      .def_property_readonly("orbit", &TaoLatticeBranchProxy::orbit)
      // TaoLatticeBranchProxy.plot_cache (1D_ALLOC_type - Plotting data cache
      .def_property_readonly("plot_cache", &TaoLatticeBranchProxy::plot_cache)
      // TaoLatticeBranchProxy.spin (0D_NOT_type -
      .def_property(
          "spin",
          &TaoLatticeBranchProxy::spin,
          &TaoLatticeBranchProxy::set_spin)
      // TaoLatticeBranchProxy.srdt (0D_NOT_type -
      .def_property(
          "srdt",
          &TaoLatticeBranchProxy::srdt,
          &TaoLatticeBranchProxy::set_srdt)
      // TaoLatticeBranchProxy.orb0 (0D_NOT_type - For saving beginning orbit
      .def_property(
          "orb0",
          &TaoLatticeBranchProxy::orb0,
          &TaoLatticeBranchProxy::set_orb0)
      // TaoLatticeBranchProxy.modes_ri (0D_NOT_type - Synchrotron integrals stuff
      .def_property(
          "modes_ri",
          &TaoLatticeBranchProxy::modes_ri,
          &TaoLatticeBranchProxy::set_modes_ri)
      // TaoLatticeBranchProxy.modes_6d (0D_NOT_type - 6D radiation matrices.
      .def_property(
          "modes_6d",
          &TaoLatticeBranchProxy::modes_6d,
          &TaoLatticeBranchProxy::set_modes_6d)
      // TaoLatticeBranchProxy.ptc_normal_form (0D_NOT_type - Collection of normal form structures defined in PTC
      .def_property(
          "ptc_normal_form",
          &TaoLatticeBranchProxy::ptc_normal_form,
          &TaoLatticeBranchProxy::set_ptc_normal_form)
      // TaoLatticeBranchProxy.bmad_normal_form (0D_NOT_type - Collection of normal form structures defined in Bmad
      .def_property(
          "bmad_normal_form",
          &TaoLatticeBranchProxy::bmad_normal_form,
          &TaoLatticeBranchProxy::set_bmad_normal_form)
      // TaoLatticeBranchProxy.high_E_orb (1D_ALLOC_type -
      .def_property_readonly("high_E_orb", &TaoLatticeBranchProxy::high_E_orb)
      // TaoLatticeBranchProxy.low_E_orb (1D_ALLOC_type -
      .def_property_readonly("low_E_orb", &TaoLatticeBranchProxy::low_E_orb)
      // TaoLatticeBranchProxy.taylor_save (1D_NOT_type - Save to reduce computation time.
      .def_property_readonly("taylor_save", &TaoLatticeBranchProxy::taylor_save)
      // TaoLatticeBranchProxy.cache_x_min (0D_NOT_real -
      .def_property(
          "cache_x_min",
          &TaoLatticeBranchProxy::cache_x_min,
          &TaoLatticeBranchProxy::set_cache_x_min)
      // TaoLatticeBranchProxy.cache_x_max (0D_NOT_real -
      .def_property(
          "cache_x_max",
          &TaoLatticeBranchProxy::cache_x_max,
          &TaoLatticeBranchProxy::set_cache_x_max)
      // TaoLatticeBranchProxy.comb_ds_save (0D_NOT_real - Master parameter for %bunch_params_comb(:)%ds_save
      .def_property(
          "comb_ds_save",
          &TaoLatticeBranchProxy::comb_ds_save,
          &TaoLatticeBranchProxy::set_comb_ds_save)
      // TaoLatticeBranchProxy.ix_ref_taylor (0D_NOT_integer -
      .def_property(
          "ix_ref_taylor",
          &TaoLatticeBranchProxy::ix_ref_taylor,
          &TaoLatticeBranchProxy::set_ix_ref_taylor)
      // TaoLatticeBranchProxy.ix_ele_taylor (0D_NOT_integer -
      .def_property(
          "ix_ele_taylor",
          &TaoLatticeBranchProxy::ix_ele_taylor,
          &TaoLatticeBranchProxy::set_ix_ele_taylor)
      // TaoLatticeBranchProxy.track_state (0D_NOT_integer -
      .def_property(
          "track_state",
          &TaoLatticeBranchProxy::track_state,
          &TaoLatticeBranchProxy::set_track_state)
      // TaoLatticeBranchProxy.cache_n_pts (0D_NOT_integer -
      .def_property(
          "cache_n_pts",
          &TaoLatticeBranchProxy::cache_n_pts,
          &TaoLatticeBranchProxy::set_cache_n_pts)
      // TaoLatticeBranchProxy.ix_rad_int_cache (0D_NOT_integer - Radiation integrals cache index.
      .def_property(
          "ix_rad_int_cache",
          &TaoLatticeBranchProxy::ix_rad_int_cache,
          &TaoLatticeBranchProxy::set_ix_rad_int_cache)
      // TaoLatticeBranchProxy.has_open_match_element (0D_NOT_logical -
      .def_property(
          "has_open_match_element",
          &TaoLatticeBranchProxy::has_open_match_element,
          &TaoLatticeBranchProxy::set_has_open_match_element)
      // TaoLatticeBranchProxy.plot_cache_valid (0D_NOT_logical - Valid plotting data cache?
      .def_property(
          "plot_cache_valid",
          &TaoLatticeBranchProxy::plot_cache_valid,
          &TaoLatticeBranchProxy::set_plot_cache_valid)
      // TaoLatticeBranchProxy.spin_map_valid (0D_NOT_logical -
      .def_property(
          "spin_map_valid",
          &TaoLatticeBranchProxy::spin_map_valid,
          &TaoLatticeBranchProxy::set_spin_map_valid)
      // TaoLatticeBranchProxy.twiss_valid (0D_NOT_logical - Invalid EG with unstable 1-turn matrix with a closed branch. With open branch: twiss_valid = T even if some Twiss (and orbit) is invalid.
      .def_property(
          "twiss_valid",
          &TaoLatticeBranchProxy::twiss_valid,
          &TaoLatticeBranchProxy::set_twiss_valid)
      // TaoLatticeBranchProxy.mode_flip_here (0D_NOT_logical - Twiss parameter mode flip seen?
      .def_property(
          "mode_flip_here",
          &TaoLatticeBranchProxy::mode_flip_here,
          &TaoLatticeBranchProxy::set_mode_flip_here)
      // TaoLatticeBranchProxy.chrom_calc_ok (0D_NOT_logical -
      .def_property(
          "chrom_calc_ok",
          &TaoLatticeBranchProxy::chrom_calc_ok,
          &TaoLatticeBranchProxy::set_chrom_calc_ok)
      // TaoLatticeBranchProxy.rad_int_calc_ok (0D_NOT_logical -
      .def_property(
          "rad_int_calc_ok",
          &TaoLatticeBranchProxy::rad_int_calc_ok,
          &TaoLatticeBranchProxy::set_rad_int_calc_ok)
      // TaoLatticeBranchProxy.emit_6d_calc_ok (0D_NOT_logical -
      .def_property(
          "emit_6d_calc_ok",
          &TaoLatticeBranchProxy::emit_6d_calc_ok,
          &TaoLatticeBranchProxy::set_emit_6d_calc_ok)
      // TaoLatticeBranchProxy.sigma_track_ok (0D_NOT_logical -
      .def_property(
          "sigma_track_ok",
          &TaoLatticeBranchProxy::sigma_track_ok,
          &TaoLatticeBranchProxy::set_sigma_track_ok)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoLatticeBranchProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoLatticeBranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatticeBranchProxy& self) {
            return TaoLatticeBranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoLatticeBranchProxy& self, py::dict& memo) {
            return TaoLatticeBranchProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoLatticeBranchProxyArray1D>(
      m, "TaoLatticeBranchStructArray1D");
  bind_FTypeAlloc1D<TaoLatticeBranchProxyAlloc1D>(
      m, "TaoLatticeBranchStructAlloc1D");
  // 2D TaoLatticeBranchProxy arrays are not used in structs/routines
  // 3D TaoLatticeBranchProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_struct
void init_tao_lattice_struct(py::module& m, py::class_<TaoLatticeProxy>& cls) {
  cls.def(py::init<>())
      // TaoLatticeProxy.name (0D_NOT_character - 'model', 'base', or 'design'.
      .def_property("name", &TaoLatticeProxy::name, &TaoLatticeProxy::set_name)
      // TaoLatticeProxy.lat (0D_NOT_type - lattice structures
      .def_property("lat", &TaoLatticeProxy::lat, &TaoLatticeProxy::set_lat)
      // TaoLatticeProxy.high_E_lat (0D_NOT_type - For chrom calc.
      .def_property(
          "high_E_lat",
          &TaoLatticeProxy::high_E_lat,
          &TaoLatticeProxy::set_high_E_lat)
      // TaoLatticeProxy.low_E_lat (0D_NOT_type - For chrom calc.
      .def_property(
          "low_E_lat",
          &TaoLatticeProxy::low_E_lat,
          &TaoLatticeProxy::set_low_E_lat)
      // TaoLatticeProxy.rad_int_by_ele_ri (0D_NOT_type -
      .def_property(
          "rad_int_by_ele_ri",
          &TaoLatticeProxy::rad_int_by_ele_ri,
          &TaoLatticeProxy::set_rad_int_by_ele_ri)
      // TaoLatticeProxy.rad_int_by_ele_6d (0D_NOT_type -
      .def_property(
          "rad_int_by_ele_6d",
          &TaoLatticeProxy::rad_int_by_ele_6d,
          &TaoLatticeProxy::set_rad_int_by_ele_6d)
      // TaoLatticeProxy.tao_branch (1D_ALLOC_type -
      .def_property_readonly("tao_branch", &TaoLatticeProxy::tao_branch)

      .def(
          "__repr__",
          [](const TaoLatticeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatticeProxy& self) {
            return TaoLatticeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoLatticeProxy& self, py::dict& memo) {
            return TaoLatticeProxy(self);
          })

      ;

  // 1D TaoLatticeProxy arrays are not used in structs/routines
  // 2D TaoLatticeProxy arrays are not used in structs/routines
  // 3D TaoLatticeProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_model_branch_struct
void init_tao_model_branch_struct(
    py::module& m,
    py::class_<TaoModelBranchProxy>& cls) {
  cls.def(py::init<>())
      // TaoModelBranchProxy.ele (1D_ALLOC_type - Per element information
      .def_property_readonly("ele", &TaoModelBranchProxy::ele)
      // TaoModelBranchProxy.beam (0D_NOT_type -
      .def_property(
          "beam", &TaoModelBranchProxy::beam, &TaoModelBranchProxy::set_beam)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoModelBranchProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoModelBranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoModelBranchProxy& self) {
            return TaoModelBranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoModelBranchProxy& self, py::dict& memo) {
            return TaoModelBranchProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoModelBranchProxyArray1D>(
      m, "TaoModelBranchStructArray1D");
  bind_FTypeAlloc1D<TaoModelBranchProxyAlloc1D>(
      m, "TaoModelBranchStructAlloc1D");
  // 2D TaoModelBranchProxy arrays are not used in structs/routines
  // 3D TaoModelBranchProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_model_element_struct
void init_tao_model_element_struct(
    py::module& m,
    py::class_<TaoModelElementProxy>& cls) {
  cls.def(py::init<>())
      // TaoModelElementProxy.beam (0D_NOT_type - Beam distribution at element.
      .def_property(
          "beam", &TaoModelElementProxy::beam, &TaoModelElementProxy::set_beam)
      // TaoModelElementProxy.save_beam_internally (0D_NOT_logical - Save beam here? Beam also saved at fork elements and at track ends.
      .def_property(
          "save_beam_internally",
          &TaoModelElementProxy::save_beam_internally,
          &TaoModelElementProxy::set_save_beam_internally)
      // TaoModelElementProxy.save_beam_to_file (0D_NOT_logical - Save beam to a file? Beam also saved at fork elements and at track ends.
      .def_property(
          "save_beam_to_file",
          &TaoModelElementProxy::save_beam_to_file,
          &TaoModelElementProxy::set_save_beam_to_file)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoModelElementProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoModelElementProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoModelElementProxy& self) {
            return TaoModelElementProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoModelElementProxy& self, py::dict& memo) {
            return TaoModelElementProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoModelElementProxyArray1D>(
      m, "TaoModelElementStructArray1D");
  bind_FTypeAlloc1D<TaoModelElementProxyAlloc1D>(
      m, "TaoModelElementStructAlloc1D");
  // 2D TaoModelElementProxy arrays are not used in structs/routines
  // 3D TaoModelElementProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_ping_scale_struct
void init_tao_ping_scale_struct(
    py::module& m,
    py::class_<TaoPingScaleProxy>& cls) {
  cls.def(py::init<>())
      // TaoPingScaleProxy.a_mode_meas (0D_NOT_real -
      .def_property(
          "a_mode_meas",
          &TaoPingScaleProxy::a_mode_meas,
          &TaoPingScaleProxy::set_a_mode_meas)
      // TaoPingScaleProxy.a_mode_ref (0D_NOT_real -
      .def_property(
          "a_mode_ref",
          &TaoPingScaleProxy::a_mode_ref,
          &TaoPingScaleProxy::set_a_mode_ref)
      // TaoPingScaleProxy.b_mode_meas (0D_NOT_real -
      .def_property(
          "b_mode_meas",
          &TaoPingScaleProxy::b_mode_meas,
          &TaoPingScaleProxy::set_b_mode_meas)
      // TaoPingScaleProxy.b_mode_ref (0D_NOT_real -
      .def_property(
          "b_mode_ref",
          &TaoPingScaleProxy::b_mode_ref,
          &TaoPingScaleProxy::set_b_mode_ref)

      .def(
          "__repr__",
          [](const TaoPingScaleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPingScaleProxy& self) {
            return TaoPingScaleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoPingScaleProxy& self, py::dict& memo) {
            return TaoPingScaleProxy(self);
          })

      ;

  // 1D TaoPingScaleProxy arrays are not used in structs/routines
  // 2D TaoPingScaleProxy arrays are not used in structs/routines
  // 3D TaoPingScaleProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_cache_struct
void init_tao_plot_cache_struct(
    py::module& m,
    py::class_<TaoPlotCacheProxy>& cls) {
  cls.def(py::init<>())
      // TaoPlotCacheProxy.ele_to_s (0D_NOT_type - Integrated element from branch beginning. Will be marked as a hybrid element.
      .def_property(
          "ele_to_s",
          &TaoPlotCacheProxy::ele_to_s,
          &TaoPlotCacheProxy::set_ele_to_s)
      // TaoPlotCacheProxy.orbit (0D_NOT_type -
      .def_property(
          "orbit", &TaoPlotCacheProxy::orbit, &TaoPlotCacheProxy::set_orbit)
      // TaoPlotCacheProxy.err (0D_NOT_logical -
      .def_property("err", &TaoPlotCacheProxy::err, &TaoPlotCacheProxy::set_err)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoPlotCacheProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoPlotCacheProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotCacheProxy& self) {
            return TaoPlotCacheProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoPlotCacheProxy& self, py::dict& memo) {
            return TaoPlotCacheProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoPlotCacheProxyArray1D>(m, "TaoPlotCacheStructArray1D");
  bind_FTypeAlloc1D<TaoPlotCacheProxyAlloc1D>(m, "TaoPlotCacheStructAlloc1D");
  // 2D TaoPlotCacheProxy arrays are not used in structs/routines
  // 3D TaoPlotCacheProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_page_struct
void init_tao_plot_page_struct(
    py::module& m,
    py::class_<TaoPlotPageProxy>& cls) {
  cls.def(py::init<>())
      // TaoPlotPageProxy.title (0D_NOT_type - Title  at top of page.
      .def_property(
          "title", &TaoPlotPageProxy::title, &TaoPlotPageProxy::set_title)
      // TaoPlotPageProxy.subtitle (0D_NOT_type - Subtitle below title at top of page.
      .def_property(
          "subtitle",
          &TaoPlotPageProxy::subtitle,
          &TaoPlotPageProxy::set_subtitle)
      // TaoPlotPageProxy.border (0D_NOT_type - Border around plots edge of page.
      .def_property(
          "border", &TaoPlotPageProxy::border, &TaoPlotPageProxy::set_border)
      // TaoPlotPageProxy.floor_plan (0D_NOT_type -
      .def_property(
          "floor_plan",
          &TaoPlotPageProxy::floor_plan,
          &TaoPlotPageProxy::set_floor_plan)
      // TaoPlotPageProxy.lat_layout (0D_NOT_type -
      .def_property(
          "lat_layout",
          &TaoPlotPageProxy::lat_layout,
          &TaoPlotPageProxy::set_lat_layout)
      // TaoPlotPageProxy.pattern (1D_ALLOC_type -
      .def_property_readonly("pattern", &TaoPlotPageProxy::pattern)
      // TaoPlotPageProxy.template_ (1D_ALLOC_type - Templates for the plots.
      .def_property_readonly("template_", &TaoPlotPageProxy::template_)
      // TaoPlotPageProxy.region (1D_ALLOC_type -
      .def_property_readonly("region", &TaoPlotPageProxy::region)
      // TaoPlotPageProxy.plot_display_type (0D_NOT_character - 'X' or 'TK'
      .def_property(
          "plot_display_type",
          &TaoPlotPageProxy::plot_display_type,
          &TaoPlotPageProxy::set_plot_display_type)
      // TaoPlotPageProxy.size (1D_NOT_real - width and height of plot window in pixels.
      .def_property_readonly("size", &TaoPlotPageProxy::size)
      // TaoPlotPageProxy.text_height (0D_NOT_real - In points. Scales the height of all text
      .def_property(
          "text_height",
          &TaoPlotPageProxy::text_height,
          &TaoPlotPageProxy::set_text_height)
      // TaoPlotPageProxy.main_title_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "main_title_text_scale",
          &TaoPlotPageProxy::main_title_text_scale,
          &TaoPlotPageProxy::set_main_title_text_scale)
      // TaoPlotPageProxy.graph_title_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "graph_title_text_scale",
          &TaoPlotPageProxy::graph_title_text_scale,
          &TaoPlotPageProxy::set_graph_title_text_scale)
      // TaoPlotPageProxy.axis_number_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "axis_number_text_scale",
          &TaoPlotPageProxy::axis_number_text_scale,
          &TaoPlotPageProxy::set_axis_number_text_scale)
      // TaoPlotPageProxy.axis_label_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "axis_label_text_scale",
          &TaoPlotPageProxy::axis_label_text_scale,
          &TaoPlotPageProxy::set_axis_label_text_scale)
      // TaoPlotPageProxy.legend_text_scale (0D_NOT_real - Relative to text_height. For legends, plot_page, and lat_layout
      .def_property(
          "legend_text_scale",
          &TaoPlotPageProxy::legend_text_scale,
          &TaoPlotPageProxy::set_legend_text_scale)
      // TaoPlotPageProxy.key_table_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "key_table_text_scale",
          &TaoPlotPageProxy::key_table_text_scale,
          &TaoPlotPageProxy::set_key_table_text_scale)
      // TaoPlotPageProxy.floor_plan_shape_scale (0D_NOT_real -
      .def_property(
          "floor_plan_shape_scale",
          &TaoPlotPageProxy::floor_plan_shape_scale,
          &TaoPlotPageProxy::set_floor_plan_shape_scale)
      // TaoPlotPageProxy.floor_plan_text_scale (0D_NOT_real - Scale used = floor_plan_text_scale * legend_text_scale
      .def_property(
          "floor_plan_text_scale",
          &TaoPlotPageProxy::floor_plan_text_scale,
          &TaoPlotPageProxy::set_floor_plan_text_scale)
      // TaoPlotPageProxy.lat_layout_shape_scale (0D_NOT_real -
      .def_property(
          "lat_layout_shape_scale",
          &TaoPlotPageProxy::lat_layout_shape_scale,
          &TaoPlotPageProxy::set_lat_layout_shape_scale)
      // TaoPlotPageProxy.lat_layout_text_scale (0D_NOT_real - Scale used = lat_layout_text_scale * legend_text_scale
      .def_property(
          "lat_layout_text_scale",
          &TaoPlotPageProxy::lat_layout_text_scale,
          &TaoPlotPageProxy::set_lat_layout_text_scale)
      // TaoPlotPageProxy.n_curve_pts (0D_NOT_integer - Default number of points for plotting a smooth curve.
      .def_property(
          "n_curve_pts",
          &TaoPlotPageProxy::n_curve_pts,
          &TaoPlotPageProxy::set_n_curve_pts)
      // TaoPlotPageProxy.id_window (0D_NOT_integer - X window id number.
      .def_property(
          "id_window",
          &TaoPlotPageProxy::id_window,
          &TaoPlotPageProxy::set_id_window)
      // TaoPlotPageProxy.delete_overlapping_plots (0D_NOT_logical - Delete overlapping plots when a plot is placed?
      .def_property(
          "delete_overlapping_plots",
          &TaoPlotPageProxy::delete_overlapping_plots,
          &TaoPlotPageProxy::set_delete_overlapping_plots)
      // TaoPlotPageProxy.draw_graph_title_suffix (0D_NOT_logical - Draw the graph title suffix?
      .def_property(
          "draw_graph_title_suffix",
          &TaoPlotPageProxy::draw_graph_title_suffix,
          &TaoPlotPageProxy::set_draw_graph_title_suffix)

      .def(
          "__repr__",
          [](const TaoPlotPageProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotPageProxy& self) {
            return TaoPlotPageProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoPlotPageProxy& self, py::dict& memo) {
            return TaoPlotPageProxy(self);
          })

      ;

  // 1D TaoPlotPageProxy arrays are not used in structs/routines
  // 2D TaoPlotPageProxy arrays are not used in structs/routines
  // 3D TaoPlotPageProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_region_struct
void init_tao_plot_region_struct(
    py::module& m,
    py::class_<TaoPlotRegionProxy>& cls) {
  cls.def(py::init<>())
      // TaoPlotRegionProxy.name (0D_NOT_character - Region name. Eg: 'r13', etc.
      .def_property(
          "name", &TaoPlotRegionProxy::name, &TaoPlotRegionProxy::set_name)
      // TaoPlotRegionProxy.plot (0D_NOT_type - Plot associated with this region
      .def_property(
          "plot", &TaoPlotRegionProxy::plot, &TaoPlotRegionProxy::set_plot)
      // TaoPlotRegionProxy.location (1D_NOT_real - [x1, x2, y1, y2] location on page.
      .def_property_readonly("location", &TaoPlotRegionProxy::location)
      // TaoPlotRegionProxy.visible (0D_NOT_logical - To draw or not to draw.
      .def_property(
          "visible",
          &TaoPlotRegionProxy::visible,
          &TaoPlotRegionProxy::set_visible)
      // TaoPlotRegionProxy.list_with_show_plot_command (0D_NOT_logical - False used for default plots to shorten the output of 'show plot'
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotRegionProxy::list_with_show_plot_command,
          &TaoPlotRegionProxy::set_list_with_show_plot_command)
      // TaoPlotRegionProxy.setup_done (0D_NOT_logical - Used for plot bookkeeping.
      .def_property(
          "setup_done",
          &TaoPlotRegionProxy::setup_done,
          &TaoPlotRegionProxy::set_setup_done)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoPlotRegionProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoPlotRegionProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotRegionProxy& self) {
            return TaoPlotRegionProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoPlotRegionProxy& self, py::dict& memo) {
            return TaoPlotRegionProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoPlotRegionProxyArray1D>(m, "TaoPlotRegionStructArray1D");
  bind_FTypeAlloc1D<TaoPlotRegionProxyAlloc1D>(m, "TaoPlotRegionStructAlloc1D");
  // 2D TaoPlotRegionProxy arrays are not used in structs/routines
  // 3D TaoPlotRegionProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_struct
void init_tao_plot_struct(py::module& m, py::class_<TaoPlotProxy>& cls) {
  cls.def(py::init<>())
      // TaoPlotProxy.name (0D_NOT_character - Identifying name. Rule: If name is blank, plot is not valid.
      .def_property("name", &TaoPlotProxy::name, &TaoPlotProxy::set_name)
      // TaoPlotProxy.description (0D_NOT_character - Descriptive string.
      .def_property(
          "description",
          &TaoPlotProxy::description,
          &TaoPlotProxy::set_description)
      // TaoPlotProxy.graph (1D_ALLOC_type - individual graphs of a plot
      .def_property_readonly("graph", &TaoPlotProxy::graph)
      // TaoPlotProxy.r (0D_PTR_type - pointer to parent.
      .def_property("r", &TaoPlotProxy::r, &TaoPlotProxy::set_r)
      // TaoPlotProxy.ix_plot (0D_NOT_integer - Index in s%plot_page%template(:) or %region(:) arrays.
      .def_property(
          "ix_plot", &TaoPlotProxy::ix_plot, &TaoPlotProxy::set_ix_plot)
      // TaoPlotProxy.n_curve_pts (0D_NOT_integer - Overrides s%plot_page%n_curve_pts.
      .def_property(
          "n_curve_pts",
          &TaoPlotProxy::n_curve_pts,
          &TaoPlotProxy::set_n_curve_pts)
      // TaoPlotProxy.type (0D_NOT_character - or 'wave'
      .def_property("type", &TaoPlotProxy::type, &TaoPlotProxy::set_type)
      // TaoPlotProxy.x_axis_type (0D_NOT_character - 'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc.
      .def_property(
          "x_axis_type",
          &TaoPlotProxy::x_axis_type,
          &TaoPlotProxy::set_x_axis_type)
      // TaoPlotProxy.autoscale_x (0D_NOT_logical - Horizontal autoscale.
      .def_property(
          "autoscale_x",
          &TaoPlotProxy::autoscale_x,
          &TaoPlotProxy::set_autoscale_x)
      // TaoPlotProxy.autoscale_y (0D_NOT_logical - Vertical autoscale.
      .def_property(
          "autoscale_y",
          &TaoPlotProxy::autoscale_y,
          &TaoPlotProxy::set_autoscale_y)
      // TaoPlotProxy.autoscale_gang_x (0D_NOT_logical - scale cmd scales graphs together?
      .def_property(
          "autoscale_gang_x",
          &TaoPlotProxy::autoscale_gang_x,
          &TaoPlotProxy::set_autoscale_gang_x)
      // TaoPlotProxy.autoscale_gang_y (0D_NOT_logical - scale cmd scales graphs together?
      .def_property(
          "autoscale_gang_y",
          &TaoPlotProxy::autoscale_gang_y,
          &TaoPlotProxy::set_autoscale_gang_y)
      // TaoPlotProxy.list_with_show_plot_command (0D_NOT_logical - False used for default plots to shorten the output of 'show plot'
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotProxy::list_with_show_plot_command,
          &TaoPlotProxy::set_list_with_show_plot_command)
      // TaoPlotProxy.phantom (0D_NOT_logical - Used by tao_plot_init to add info lines to 'show plot -templates'
      .def_property(
          "phantom", &TaoPlotProxy::phantom, &TaoPlotProxy::set_phantom)
      // TaoPlotProxy.default_plot (0D_NOT_logical - One of Tao's default plots?
      .def_property(
          "default_plot",
          &TaoPlotProxy::default_plot,
          &TaoPlotProxy::set_default_plot)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoPlotProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const TaoPlotProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotProxy& self) {
            return TaoPlotProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoPlotProxy& self, py::dict& memo) {
            return TaoPlotProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoPlotProxyArray1D>(m, "TaoPlotStructArray1D");
  bind_FTypeAlloc1D<TaoPlotProxyAlloc1D>(m, "TaoPlotStructAlloc1D");
  // 2D TaoPlotProxy arrays are not used in structs/routines
  // 3D TaoPlotProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_point_struct
void init_tao_shape_pattern_point_struct(
    py::module& m,
    py::class_<TaoShapePatternPointProxy>& cls) {
  cls.def(py::init<>())
      // TaoShapePatternPointProxy.s (0D_NOT_real -
      .def_property(
          "s", &TaoShapePatternPointProxy::s, &TaoShapePatternPointProxy::set_s)
      // TaoShapePatternPointProxy.y (0D_NOT_real -
      .def_property(
          "y", &TaoShapePatternPointProxy::y, &TaoShapePatternPointProxy::set_y)
      // TaoShapePatternPointProxy.radius (0D_NOT_real -
      .def_property(
          "radius",
          &TaoShapePatternPointProxy::radius,
          &TaoShapePatternPointProxy::set_radius)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoShapePatternPointProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoShapePatternPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoShapePatternPointProxy& self) {
            return TaoShapePatternPointProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoShapePatternPointProxy& self, py::dict& memo) {
            return TaoShapePatternPointProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoShapePatternPointProxyArray1D>(
      m, "TaoShapePatternPointStructArray1D");
  bind_FTypeAlloc1D<TaoShapePatternPointProxyAlloc1D>(
      m, "TaoShapePatternPointStructAlloc1D");
  // 2D TaoShapePatternPointProxy arrays are not used in structs/routines
  // 3D TaoShapePatternPointProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_struct
void init_tao_shape_pattern_struct(
    py::module& m,
    py::class_<TaoShapePatternProxy>& cls) {
  cls.def(py::init<>())
      // TaoShapePatternProxy.name (0D_NOT_character -
      .def_property(
          "name", &TaoShapePatternProxy::name, &TaoShapePatternProxy::set_name)
      // TaoShapePatternProxy.line (0D_NOT_type - Line color and pattern set by shape using this pattern.
      .def_property(
          "line", &TaoShapePatternProxy::line, &TaoShapePatternProxy::set_line)
      // TaoShapePatternProxy.pt (1D_ALLOC_type -
      .def_property_readonly("pt", &TaoShapePatternProxy::pt)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoShapePatternProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoShapePatternProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoShapePatternProxy& self) {
            return TaoShapePatternProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoShapePatternProxy& self, py::dict& memo) {
            return TaoShapePatternProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoShapePatternProxyArray1D>(
      m, "TaoShapePatternStructArray1D");
  bind_FTypeAlloc1D<TaoShapePatternProxyAlloc1D>(
      m, "TaoShapePatternStructAlloc1D");
  // 2D TaoShapePatternProxy arrays are not used in structs/routines
  // 3D TaoShapePatternProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_dn_dpz_struct
void init_tao_spin_dn_dpz_struct(
    py::module& m,
    py::class_<TaoSpinDnDpzProxy>& cls) {
  cls.def(py::init<>())
      // TaoSpinDnDpzProxy.vec (1D_NOT_real - n0 derivative wrt pz.
      .def_property_readonly("vec", &TaoSpinDnDpzProxy::vec)
      // TaoSpinDnDpzProxy.partial (2D_NOT_real - partial(i:) is spin n0 derivative wrt pz for i^th oscillation mode (1 => a-mode, etc.)
      .def_property_readonly("partial", &TaoSpinDnDpzProxy::partial)
      // TaoSpinDnDpzProxy.partial2 (2D_NOT_real - partial(i:) is spin n0 derivative wrt pz with i^th oscillation mode missing (1 => a-mode, etc.)
      .def_property_readonly("partial2", &TaoSpinDnDpzProxy::partial2)

      .def(
          "__repr__",
          [](const TaoSpinDnDpzProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinDnDpzProxy& self) {
            return TaoSpinDnDpzProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoSpinDnDpzProxy& self, py::dict& memo) {
            return TaoSpinDnDpzProxy(self);
          })

      ;

  // 1D TaoSpinDnDpzProxy arrays are not used in structs/routines
  // 2D TaoSpinDnDpzProxy arrays are not used in structs/routines
  // 3D TaoSpinDnDpzProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_ele_struct
void init_tao_spin_ele_struct(py::module& m, py::class_<TaoSpinEleProxy>& cls) {
  cls.def(py::init<>())
      // TaoSpinEleProxy.dn_dpz (0D_NOT_type -
      .def_property(
          "dn_dpz", &TaoSpinEleProxy::dn_dpz, &TaoSpinEleProxy::set_dn_dpz)
      // TaoSpinEleProxy.orb_eigen_val (1D_NOT_real -
      .def_property_readonly("orb_eigen_val", &TaoSpinEleProxy::orb_eigen_val)
      // TaoSpinEleProxy.orb_eigen_vec (2D_NOT_real - (j,:) is j^th vector
      .def_property_readonly("orb_eigen_vec", &TaoSpinEleProxy::orb_eigen_vec)
      // TaoSpinEleProxy.spin_eigen_vec (2D_NOT_real - (j,:) is j^th vector
      .def_property_readonly("spin_eigen_vec", &TaoSpinEleProxy::spin_eigen_vec)
      // TaoSpinEleProxy.valid (0D_NOT_logical -
      .def_property(
          "valid", &TaoSpinEleProxy::valid, &TaoSpinEleProxy::set_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoSpinEleProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoSpinEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinEleProxy& self) {
            return TaoSpinEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoSpinEleProxy& self, py::dict& memo) {
            return TaoSpinEleProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoSpinEleProxyArray1D>(m, "TaoSpinEleStructArray1D");
  bind_FTypeAlloc1D<TaoSpinEleProxyAlloc1D>(m, "TaoSpinEleStructAlloc1D");
  // 2D TaoSpinEleProxy arrays are not used in structs/routines
  // 3D TaoSpinEleProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_map_struct
void init_tao_spin_map_struct(py::module& m, py::class_<TaoSpinMapProxy>& cls) {
  cls.def(py::init<>())
      // TaoSpinMapProxy.valid (0D_NOT_logical -
      .def_property(
          "valid", &TaoSpinMapProxy::valid, &TaoSpinMapProxy::set_valid)
      // TaoSpinMapProxy.map1 (0D_NOT_type -
      .def_property("map1", &TaoSpinMapProxy::map1, &TaoSpinMapProxy::set_map1)
      // TaoSpinMapProxy.axis_input (0D_NOT_type - Input axes.
      .def_property(
          "axis_input",
          &TaoSpinMapProxy::axis_input,
          &TaoSpinMapProxy::set_axis_input)
      // TaoSpinMapProxy.axis0 (0D_NOT_type - Initial axes.
      .def_property(
          "axis0", &TaoSpinMapProxy::axis0, &TaoSpinMapProxy::set_axis0)
      // TaoSpinMapProxy.axis1 (0D_NOT_type - Final axes.
      .def_property(
          "axis1", &TaoSpinMapProxy::axis1, &TaoSpinMapProxy::set_axis1)
      // TaoSpinMapProxy.ix_ele (0D_NOT_integer -
      .def_property(
          "ix_ele", &TaoSpinMapProxy::ix_ele, &TaoSpinMapProxy::set_ix_ele)
      // TaoSpinMapProxy.ix_ref (0D_NOT_integer -
      .def_property(
          "ix_ref", &TaoSpinMapProxy::ix_ref, &TaoSpinMapProxy::set_ix_ref)
      // TaoSpinMapProxy.ix_uni (0D_NOT_integer -
      .def_property(
          "ix_uni", &TaoSpinMapProxy::ix_uni, &TaoSpinMapProxy::set_ix_uni)
      // TaoSpinMapProxy.ix_branch (0D_NOT_integer -
      .def_property(
          "ix_branch",
          &TaoSpinMapProxy::ix_branch,
          &TaoSpinMapProxy::set_ix_branch)
      // TaoSpinMapProxy.mat8 (2D_NOT_real -
      .def_property_readonly("mat8", &TaoSpinMapProxy::mat8)

      .def(
          "__repr__",
          [](const TaoSpinMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinMapProxy& self) {
            return TaoSpinMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoSpinMapProxy& self, py::dict& memo) {
            return TaoSpinMapProxy(self);
          })

      ;

  // 1D TaoSpinMapProxy arrays are not used in structs/routines
  // 2D TaoSpinMapProxy arrays are not used in structs/routines
  // 3D TaoSpinMapProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_polarization_struct
void init_tao_spin_polarization_struct(
    py::module& m,
    py::class_<TaoSpinPolarizationProxy>& cls) {
  cls.def(py::init<>())
      // TaoSpinPolarizationProxy.tune (0D_NOT_real -
      .def_property(
          "tune",
          &TaoSpinPolarizationProxy::tune,
          &TaoSpinPolarizationProxy::set_tune)
      // TaoSpinPolarizationProxy.pol_limit_st (0D_NOT_real - Polarization calculated using Sokolov-Ternov formula.
      .def_property(
          "pol_limit_st",
          &TaoSpinPolarizationProxy::pol_limit_st,
          &TaoSpinPolarizationProxy::set_pol_limit_st)
      // TaoSpinPolarizationProxy.pol_limit_dk (0D_NOT_real - Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula.
      .def_property(
          "pol_limit_dk",
          &TaoSpinPolarizationProxy::pol_limit_dk,
          &TaoSpinPolarizationProxy::set_pol_limit_dk)
      // TaoSpinPolarizationProxy.pol_limit_dk_partial (1D_NOT_real - Limit using only single mode to calc dn_dpz
      .def_property_readonly(
          "pol_limit_dk_partial",
          &TaoSpinPolarizationProxy::pol_limit_dk_partial)
      // TaoSpinPolarizationProxy.pol_limit_dk_partial2 (1D_NOT_real - Limit using only single mode to calc dn_dpz
      .def_property_readonly(
          "pol_limit_dk_partial2",
          &TaoSpinPolarizationProxy::pol_limit_dk_partial2)
      // TaoSpinPolarizationProxy.pol_rate_bks (0D_NOT_real - BKS Polarization rate (1/sec).
      .def_property(
          "pol_rate_bks",
          &TaoSpinPolarizationProxy::pol_rate_bks,
          &TaoSpinPolarizationProxy::set_pol_rate_bks)
      // TaoSpinPolarizationProxy.depol_rate (0D_NOT_real - Depolarization rate (1/sec).
      .def_property(
          "depol_rate",
          &TaoSpinPolarizationProxy::depol_rate,
          &TaoSpinPolarizationProxy::set_depol_rate)
      // TaoSpinPolarizationProxy.depol_rate_partial (1D_NOT_real - Depolarization rate (1/sec) using only single mode to calc dn_dpz.
      .def_property_readonly(
          "depol_rate_partial", &TaoSpinPolarizationProxy::depol_rate_partial)
      // TaoSpinPolarizationProxy.depol_rate_partial2 (1D_NOT_real - Depolarization rate (1/sec) using only two modes to calc dn_dpz.
      .def_property_readonly(
          "depol_rate_partial2", &TaoSpinPolarizationProxy::depol_rate_partial2)
      // TaoSpinPolarizationProxy.integral_bn (0D_NOT_real - Integral of g^3 * b_hat * n_0
      .def_property(
          "integral_bn",
          &TaoSpinPolarizationProxy::integral_bn,
          &TaoSpinPolarizationProxy::set_integral_bn)
      // TaoSpinPolarizationProxy.integral_bdn (0D_NOT_real - Integral of g^3 * b_hat * dn/ddelta
      .def_property(
          "integral_bdn",
          &TaoSpinPolarizationProxy::integral_bdn,
          &TaoSpinPolarizationProxy::set_integral_bdn)
      // TaoSpinPolarizationProxy.integral_1ns (0D_NOT_real - Integral of g^3 (1 - 2(n * s_hat)/9)
      .def_property(
          "integral_1ns",
          &TaoSpinPolarizationProxy::integral_1ns,
          &TaoSpinPolarizationProxy::set_integral_1ns)
      // TaoSpinPolarizationProxy.integral_dn2 (0D_NOT_real - Integral of g^3 * 11 (dn/ddelta)^2 / 9
      .def_property(
          "integral_dn2",
          &TaoSpinPolarizationProxy::integral_dn2,
          &TaoSpinPolarizationProxy::set_integral_dn2)
      // TaoSpinPolarizationProxy.valid (0D_NOT_logical -
      .def_property(
          "valid",
          &TaoSpinPolarizationProxy::valid,
          &TaoSpinPolarizationProxy::set_valid)
      // TaoSpinPolarizationProxy.q_1turn (0D_NOT_type - Save results from spin_concat_linear_maps in tao_spin_polarization.
      .def_property(
          "q_1turn",
          &TaoSpinPolarizationProxy::q_1turn,
          &TaoSpinPolarizationProxy::set_q_1turn)
      // TaoSpinPolarizationProxy.q_ele (1D_ALLOC_type - Save results from spin_concat_linear_maps in tao_spin_polarization.
      .def_property_readonly("q_ele", &TaoSpinPolarizationProxy::q_ele)

      .def(
          "__repr__",
          [](const TaoSpinPolarizationProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinPolarizationProxy& self) {
            return TaoSpinPolarizationProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoSpinPolarizationProxy& self, py::dict& memo) {
            return TaoSpinPolarizationProxy(self);
          })

      ;

  // 1D TaoSpinPolarizationProxy arrays are not used in structs/routines
  // 2D TaoSpinPolarizationProxy arrays are not used in structs/routines
  // 3D TaoSpinPolarizationProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_super_universe_struct
void init_tao_super_universe_struct(
    py::module& m,
    py::class_<TaoSuperUniverseProxy>& cls) {
  cls.def(py::init<>())
      // TaoSuperUniverseProxy.global (0D_NOT_type - User accessible global variables.
      .def_property(
          "global_",
          &TaoSuperUniverseProxy::global,
          &TaoSuperUniverseProxy::set_global)
      // TaoSuperUniverseProxy.init (0D_NOT_type - Initialization parameters
      .def_property(
          "init",
          &TaoSuperUniverseProxy::init,
          &TaoSuperUniverseProxy::set_init)
      // TaoSuperUniverseProxy.com (0D_NOT_type - Non-initialization common parameters
      .def_property(
          "com", &TaoSuperUniverseProxy::com, &TaoSuperUniverseProxy::set_com)
      // TaoSuperUniverseProxy.plot_page (0D_NOT_type - Defines the plot window.
      .def_property(
          "plot_page",
          &TaoSuperUniverseProxy::plot_page,
          &TaoSuperUniverseProxy::set_plot_page)
      // TaoSuperUniverseProxy.v1_var (1D_ALLOC_type - The variable types
      .def_property_readonly("v1_var", &TaoSuperUniverseProxy::v1_var)
      // TaoSuperUniverseProxy.var (1D_ALLOC_type - array of all variables.
      .def_property_readonly("var", &TaoSuperUniverseProxy::var)
      // TaoSuperUniverseProxy.u (1D_ALLOC_type - array of universes.
      .def_property_readonly("u", &TaoSuperUniverseProxy::u)
      // TaoSuperUniverseProxy.key (1D_ALLOC_integer -
      .def_property_readonly("key", &TaoSuperUniverseProxy::key)
      // TaoSuperUniverseProxy.building_wall (0D_NOT_type -
      .def_property(
          "building_wall",
          &TaoSuperUniverseProxy::building_wall,
          &TaoSuperUniverseProxy::set_building_wall)
      // TaoSuperUniverseProxy.wave (0D_NOT_type -
      .def_property(
          "wave",
          &TaoSuperUniverseProxy::wave,
          &TaoSuperUniverseProxy::set_wave)
      // TaoSuperUniverseProxy.n_var_used (0D_NOT_integer -
      .def_property(
          "n_var_used",
          &TaoSuperUniverseProxy::n_var_used,
          &TaoSuperUniverseProxy::set_n_var_used)
      // TaoSuperUniverseProxy.n_v1_var_used (0D_NOT_integer -
      .def_property(
          "n_v1_var_used",
          &TaoSuperUniverseProxy::n_v1_var_used,
          &TaoSuperUniverseProxy::set_n_v1_var_used)
      // TaoSuperUniverseProxy.history (1D_NOT_type - command history
      .def_property_readonly("history", &TaoSuperUniverseProxy::history)
      // TaoSuperUniverseProxy.initialized (0D_NOT_logical - Does tao_init() need to be called?
      .def_property(
          "initialized",
          &TaoSuperUniverseProxy::initialized,
          &TaoSuperUniverseProxy::set_initialized)

      .def(
          "__repr__",
          [](const TaoSuperUniverseProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSuperUniverseProxy& self) {
            return TaoSuperUniverseProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoSuperUniverseProxy& self, py::dict& memo) {
            return TaoSuperUniverseProxy(self);
          })

      ;

  // 1D TaoSuperUniverseProxy arrays are not used in structs/routines
  // 2D TaoSuperUniverseProxy arrays are not used in structs/routines
  // 3D TaoSuperUniverseProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_title_struct
void init_tao_title_struct(py::module& m, py::class_<TaoTitleProxy>& cls) {
  cls.def(py::init<>())
      // TaoTitleProxy.string (0D_NOT_character - title character string.
      .def_property(
          "string", &TaoTitleProxy::string, &TaoTitleProxy::set_string)
      // TaoTitleProxy.x (0D_NOT_real - x, y rwt lower left corner
      .def_property("x", &TaoTitleProxy::x, &TaoTitleProxy::set_x)
      // TaoTitleProxy.y (0D_NOT_real - x, y rwt lower left corner
      .def_property("y", &TaoTitleProxy::y, &TaoTitleProxy::set_y)
      // TaoTitleProxy.units (0D_NOT_character - %BOX, POINTS, etc...
      .def_property("units", &TaoTitleProxy::units, &TaoTitleProxy::set_units)
      // TaoTitleProxy.justify (0D_NOT_character - Left, Center, or Right justification.
      .def_property(
          "justify", &TaoTitleProxy::justify, &TaoTitleProxy::set_justify)
      // TaoTitleProxy.draw_it (0D_NOT_logical - draw the title?
      .def_property(
          "draw_it", &TaoTitleProxy::draw_it, &TaoTitleProxy::set_draw_it)

      .def(
          "__repr__", [](const TaoTitleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoTitleProxy& self) {
            return TaoTitleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoTitleProxy& self, py::dict& memo) {
            return TaoTitleProxy(self);
          })

      ;

  // 1D TaoTitleProxy arrays are not used in structs/routines
  // 2D TaoTitleProxy arrays are not used in structs/routines
  // 3D TaoTitleProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_calc_struct
void init_tao_universe_calc_struct(
    py::module& m,
    py::class_<TaoUniverseCalcProxy>& cls) {
  cls.def(py::init<>())
      // TaoUniverseCalcProxy.srdt_for_data (0D_NOT_integer - 0 = false, 1 = 1st order, 2 = 1st & 2nd order
      .def_property(
          "srdt_for_data",
          &TaoUniverseCalcProxy::srdt_for_data,
          &TaoUniverseCalcProxy::set_srdt_for_data)
      // TaoUniverseCalcProxy.rad_int_for_data (0D_NOT_logical - Do the radiation integrals need to be computed for
      .def_property(
          "rad_int_for_data",
          &TaoUniverseCalcProxy::rad_int_for_data,
          &TaoUniverseCalcProxy::set_rad_int_for_data)
      // TaoUniverseCalcProxy.rad_int_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "rad_int_for_plotting",
          &TaoUniverseCalcProxy::rad_int_for_plotting,
          &TaoUniverseCalcProxy::set_rad_int_for_plotting)
      // TaoUniverseCalcProxy.chrom_for_data (0D_NOT_logical - Does the chromaticity need to be computed for
      .def_property(
          "chrom_for_data",
          &TaoUniverseCalcProxy::chrom_for_data,
          &TaoUniverseCalcProxy::set_chrom_for_data)
      // TaoUniverseCalcProxy.chrom_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "chrom_for_plotting",
          &TaoUniverseCalcProxy::chrom_for_plotting,
          &TaoUniverseCalcProxy::set_chrom_for_plotting)
      // TaoUniverseCalcProxy.lat_sigma_for_data (0D_NOT_logical - Do the beam sigmas need to be computed for
      .def_property(
          "lat_sigma_for_data",
          &TaoUniverseCalcProxy::lat_sigma_for_data,
          &TaoUniverseCalcProxy::set_lat_sigma_for_data)
      // TaoUniverseCalcProxy.lat_sigma_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "lat_sigma_for_plotting",
          &TaoUniverseCalcProxy::lat_sigma_for_plotting,
          &TaoUniverseCalcProxy::set_lat_sigma_for_plotting)
      // TaoUniverseCalcProxy.dynamic_aperture (0D_NOT_logical - Do the dynamic_aperture calc?
      .def_property(
          "dynamic_aperture",
          &TaoUniverseCalcProxy::dynamic_aperture,
          &TaoUniverseCalcProxy::set_dynamic_aperture)
      // TaoUniverseCalcProxy.one_turn_map (0D_NOT_logical - Compute the one turn map?
      .def_property(
          "one_turn_map",
          &TaoUniverseCalcProxy::one_turn_map,
          &TaoUniverseCalcProxy::set_one_turn_map)
      // TaoUniverseCalcProxy.lattice (0D_NOT_logical - Used to indicate which lattices need tracking done.
      .def_property(
          "lattice",
          &TaoUniverseCalcProxy::lattice,
          &TaoUniverseCalcProxy::set_lattice)
      // TaoUniverseCalcProxy.twiss (0D_NOT_logical - calc linear transfer matrix?
      .def_property(
          "twiss",
          &TaoUniverseCalcProxy::twiss,
          &TaoUniverseCalcProxy::set_twiss)
      // TaoUniverseCalcProxy.track (0D_NOT_logical - tracking needs to be done?
      .def_property(
          "track",
          &TaoUniverseCalcProxy::track,
          &TaoUniverseCalcProxy::set_track)
      // TaoUniverseCalcProxy.spin_matrices (0D_NOT_logical - Calculate G and D spin matrices?
      .def_property(
          "spin_matrices",
          &TaoUniverseCalcProxy::spin_matrices,
          &TaoUniverseCalcProxy::set_spin_matrices)

      .def(
          "__repr__",
          [](const TaoUniverseCalcProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniverseCalcProxy& self) {
            return TaoUniverseCalcProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoUniverseCalcProxy& self, py::dict& memo) {
            return TaoUniverseCalcProxy(self);
          })

      ;

  // 1D TaoUniverseCalcProxy arrays are not used in structs/routines
  // 2D TaoUniverseCalcProxy arrays are not used in structs/routines
  // 3D TaoUniverseCalcProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_pointer_struct
void init_tao_universe_pointer_struct(
    py::module& m,
    py::class_<TaoUniversePointerProxy>& cls) {
  cls.def(py::init<>())
      // TaoUniversePointerProxy.u (0D_PTR_type -
      .def_property(
          "u", &TaoUniversePointerProxy::u, &TaoUniversePointerProxy::set_u)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoUniversePointerProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoUniversePointerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniversePointerProxy& self) {
            return TaoUniversePointerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoUniversePointerProxy& self, py::dict& memo) {
            return TaoUniversePointerProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoUniversePointerProxyArray1D>(
      m, "TaoUniversePointerStructArray1D");
  bind_FTypeAlloc1D<TaoUniversePointerProxyAlloc1D>(
      m, "TaoUniversePointerStructAlloc1D");
  // 2D TaoUniversePointerProxy arrays are not used in structs/routines
  // 3D TaoUniversePointerProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_struct
void init_tao_universe_struct(
    py::module& m,
    py::class_<TaoUniverseProxy>& cls) {
  cls.def(py::init<>())
      // TaoUniverseProxy.model (0D_PTR_type -
      .def_property(
          "model", &TaoUniverseProxy::model, &TaoUniverseProxy::set_model)
      // TaoUniverseProxy.design (0D_PTR_type -
      .def_property(
          "design", &TaoUniverseProxy::design, &TaoUniverseProxy::set_design)
      // TaoUniverseProxy.base (0D_PTR_type -
      .def_property(
          "base", &TaoUniverseProxy::base, &TaoUniverseProxy::set_base)
      // TaoUniverseProxy.beam (0D_NOT_type -
      .def_property(
          "beam", &TaoUniverseProxy::beam, &TaoUniverseProxy::set_beam)
      // TaoUniverseProxy.dynamic_aperture (0D_NOT_type -
      .def_property(
          "dynamic_aperture",
          &TaoUniverseProxy::dynamic_aperture,
          &TaoUniverseProxy::set_dynamic_aperture)
      // TaoUniverseProxy.model_branch (1D_PTR_type - model specific information
      .def_property_readonly("model_branch", &TaoUniverseProxy::model_branch)
      // TaoUniverseProxy.d2_data (1D_ALLOC_type - The data types
      .def_property_readonly("d2_data", &TaoUniverseProxy::d2_data)
      // TaoUniverseProxy.data (1D_ALLOC_type - Array of all data.
      .def_property_readonly("data", &TaoUniverseProxy::data)
      // TaoUniverseProxy.ping_scale (0D_NOT_type -
      .def_property(
          "ping_scale",
          &TaoUniverseProxy::ping_scale,
          &TaoUniverseProxy::set_ping_scale)
      // TaoUniverseProxy.scratch_lat (0D_NOT_type - Scratch area.
      .def_property(
          "scratch_lat",
          &TaoUniverseProxy::scratch_lat,
          &TaoUniverseProxy::set_scratch_lat)
      // TaoUniverseProxy.calc (0D_NOT_type - What needs to be calculated?
      .def_property(
          "calc", &TaoUniverseProxy::calc, &TaoUniverseProxy::set_calc)
      // TaoUniverseProxy.ele_order (0D_NOT_type - Order of elements with same name.
      .def_property(
          "ele_order",
          &TaoUniverseProxy::ele_order,
          &TaoUniverseProxy::set_ele_order)
      // TaoUniverseProxy.spin_map (0D_NOT_type -
      .def_property(
          "spin_map",
          &TaoUniverseProxy::spin_map,
          &TaoUniverseProxy::set_spin_map)
      // TaoUniverseProxy.dModel_dVar (2D_ALLOC_real - Derivative matrix.
      .def_property_readonly("dModel_dVar", &TaoUniverseProxy::dModel_dVar)
      // TaoUniverseProxy.ix_uni (0D_NOT_integer - Universe index.
      .def_property(
          "ix_uni", &TaoUniverseProxy::ix_uni, &TaoUniverseProxy::set_ix_uni)
      // TaoUniverseProxy.n_d2_data_used (0D_NOT_integer - Number of used %d2_data(:) components.
      .def_property(
          "n_d2_data_used",
          &TaoUniverseProxy::n_d2_data_used,
          &TaoUniverseProxy::set_n_d2_data_used)
      // TaoUniverseProxy.n_data_used (0D_NOT_integer - Number of used %data(:) components.
      .def_property(
          "n_data_used",
          &TaoUniverseProxy::n_data_used,
          &TaoUniverseProxy::set_n_data_used)
      // TaoUniverseProxy.is_on (0D_NOT_logical - universe turned on
      .def_property(
          "is_on", &TaoUniverseProxy::is_on, &TaoUniverseProxy::set_is_on)
      // TaoUniverseProxy.design_same_as_previous (0D_NOT_logical - Design lat same as the previous uni?
      .def_property(
          "design_same_as_previous",
          &TaoUniverseProxy::design_same_as_previous,
          &TaoUniverseProxy::set_design_same_as_previous)
      // TaoUniverseProxy.picked_uni (0D_NOT_logical - Scratch logical.
      .def_property(
          "picked_uni",
          &TaoUniverseProxy::picked_uni,
          &TaoUniverseProxy::set_picked_uni)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoUniverseProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoUniverseProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniverseProxy& self) {
            return TaoUniverseProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoUniverseProxy& self, py::dict& memo) {
            return TaoUniverseProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoUniverseProxyArray1D>(m, "TaoUniverseStructArray1D");
  bind_FTypeAlloc1D<TaoUniverseProxyAlloc1D>(m, "TaoUniverseStructAlloc1D");
  // 2D TaoUniverseProxy arrays are not used in structs/routines
  // 3D TaoUniverseProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_v1_var_struct
void init_tao_v1_var_struct(py::module& m, py::class_<TaoV1VarProxy>& cls) {
  cls.def(py::init<>())
      // TaoV1VarProxy.name (0D_NOT_character - V1 variable name. Eg: 'quad_k1'.
      .def_property("name", &TaoV1VarProxy::name, &TaoV1VarProxy::set_name)
      // TaoV1VarProxy.ix_v1_var (0D_NOT_integer - Index to s%v1_var(:) array
      .def_property(
          "ix_v1_var", &TaoV1VarProxy::ix_v1_var, &TaoV1VarProxy::set_ix_v1_var)
      // TaoV1VarProxy.v (1D_PTR_type - Pointer to the appropriate section in s%var.
      .def_property_readonly("v", &TaoV1VarProxy::v)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoV1VarProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const TaoV1VarProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoV1VarProxy& self) {
            return TaoV1VarProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoV1VarProxy& self, py::dict& memo) {
            return TaoV1VarProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoV1VarProxyArray1D>(m, "TaoV1VarStructArray1D");
  bind_FTypeAlloc1D<TaoV1VarProxyAlloc1D>(m, "TaoV1VarStructAlloc1D");
  // 2D TaoV1VarProxy arrays are not used in structs/routines
  // 3D TaoV1VarProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_var_slave_struct
void init_tao_var_slave_struct(
    py::module& m,
    py::class_<TaoVarSlaveProxy>& cls) {
  cls.def(py::init<>())
      // TaoVarSlaveProxy.ix_uni (0D_NOT_integer - universe index.
      .def_property(
          "ix_uni", &TaoVarSlaveProxy::ix_uni, &TaoVarSlaveProxy::set_ix_uni)
      // TaoVarSlaveProxy.ix_branch (0D_NOT_integer -
      .def_property(
          "ix_branch",
          &TaoVarSlaveProxy::ix_branch,
          &TaoVarSlaveProxy::set_ix_branch)
      // TaoVarSlaveProxy.ix_ele (0D_NOT_integer - Index of element in the u%lattice%ele(:) array.
      .def_property(
          "ix_ele", &TaoVarSlaveProxy::ix_ele, &TaoVarSlaveProxy::set_ix_ele)
      // TaoVarSlaveProxy.model_value (0D_PTR_real - Pointer to the variable in the model lat.
      .def_property(
          "model_value",
          &TaoVarSlaveProxy::model_value,
          &TaoVarSlaveProxy::set_model_value)
      // TaoVarSlaveProxy.base_value (0D_PTR_real - Pointer to the variable in the base lat.
      .def_property(
          "base_value",
          &TaoVarSlaveProxy::base_value,
          &TaoVarSlaveProxy::set_base_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoVarSlaveProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoVarSlaveProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoVarSlaveProxy& self) {
            return TaoVarSlaveProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoVarSlaveProxy& self, py::dict& memo) {
            return TaoVarSlaveProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoVarSlaveProxyArray1D>(m, "TaoVarSlaveStructArray1D");
  bind_FTypeAlloc1D<TaoVarSlaveProxyAlloc1D>(m, "TaoVarSlaveStructAlloc1D");
  // 2D TaoVarSlaveProxy arrays are not used in structs/routines
  // 3D TaoVarSlaveProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_var_struct
void init_tao_var_struct(py::module& m, py::class_<TaoVarProxy>& cls) {
  cls.def(py::init<>())
      // TaoVarProxy.ele_name (0D_NOT_character - Associated lattice element name.
      .def_property(
          "ele_name", &TaoVarProxy::ele_name, &TaoVarProxy::set_ele_name)
      // TaoVarProxy.attrib_name (0D_NOT_character - Name of the attribute to vary.
      .def_property(
          "attrib_name",
          &TaoVarProxy::attrib_name,
          &TaoVarProxy::set_attrib_name)
      // TaoVarProxy.id (0D_NOT_character - Used by Tao extension code. Not used by Tao directly.
      .def_property("id", &TaoVarProxy::id, &TaoVarProxy::set_id)
      // TaoVarProxy.slave (1D_ALLOC_type -
      .def_property_readonly("slave", &TaoVarProxy::slave)
      // TaoVarProxy.ix_v1 (0D_NOT_integer - Index of this var in the s%v1_var(i)%v(:) array.
      .def_property("ix_v1", &TaoVarProxy::ix_v1, &TaoVarProxy::set_ix_v1)
      // TaoVarProxy.ix_var (0D_NOT_integer - Index number of this var in the s%var(:) array.
      .def_property("ix_var", &TaoVarProxy::ix_var, &TaoVarProxy::set_ix_var)
      // TaoVarProxy.ix_dvar (0D_NOT_integer - Column in the dData_dVar derivative matrix.
      .def_property("ix_dvar", &TaoVarProxy::ix_dvar, &TaoVarProxy::set_ix_dvar)
      // TaoVarProxy.ix_attrib (0D_NOT_integer - Index in ele%value(:) array if appropriate.
      .def_property(
          "ix_attrib", &TaoVarProxy::ix_attrib, &TaoVarProxy::set_ix_attrib)
      // TaoVarProxy.ix_key_table (0D_NOT_integer - Has a key binding?
      .def_property(
          "ix_key_table",
          &TaoVarProxy::ix_key_table,
          &TaoVarProxy::set_ix_key_table)
      // TaoVarProxy.model_value (0D_PTR_real - Model value.
      .def_property(
          "model_value",
          &TaoVarProxy::model_value,
          &TaoVarProxy::set_model_value)
      // TaoVarProxy.base_value (0D_PTR_real - Base value.
      .def_property(
          "base_value", &TaoVarProxy::base_value, &TaoVarProxy::set_base_value)
      // TaoVarProxy.design_value (0D_NOT_real - Design value from the design lattice.
      .def_property(
          "design_value",
          &TaoVarProxy::design_value,
          &TaoVarProxy::set_design_value)
      // TaoVarProxy.scratch_value (0D_NOT_real - Scratch space used by Tao.
      .def_property(
          "scratch_value",
          &TaoVarProxy::scratch_value,
          &TaoVarProxy::set_scratch_value)
      // TaoVarProxy.old_value (0D_NOT_real - Scratch space used by Tao.
      .def_property(
          "old_value", &TaoVarProxy::old_value, &TaoVarProxy::set_old_value)
      // TaoVarProxy.meas_value (0D_NOT_real - The value when the data measurement was taken.
      .def_property(
          "meas_value", &TaoVarProxy::meas_value, &TaoVarProxy::set_meas_value)
      // TaoVarProxy.ref_value (0D_NOT_real - Value when the reference measurement was taken.
      .def_property(
          "ref_value", &TaoVarProxy::ref_value, &TaoVarProxy::set_ref_value)
      // TaoVarProxy.correction_value (0D_NOT_real - Value determined by a fit to correct the lattice.
      .def_property(
          "correction_value",
          &TaoVarProxy::correction_value,
          &TaoVarProxy::set_correction_value)
      // TaoVarProxy.high_lim (0D_NOT_real - High limit for the model_value.
      .def_property(
          "high_lim", &TaoVarProxy::high_lim, &TaoVarProxy::set_high_lim)
      // TaoVarProxy.low_lim (0D_NOT_real - Low limit for the model_value.
      .def_property("low_lim", &TaoVarProxy::low_lim, &TaoVarProxy::set_low_lim)
      // TaoVarProxy.step (0D_NOT_real - Sets what is a small step for varying this var.
      .def_property("step", &TaoVarProxy::step, &TaoVarProxy::set_step)
      // TaoVarProxy.weight (0D_NOT_real - Weight for the merit function term.
      .def_property("weight", &TaoVarProxy::weight, &TaoVarProxy::set_weight)
      // TaoVarProxy.delta_merit (0D_NOT_real - Diff used to calculate the merit function term.
      .def_property(
          "delta_merit",
          &TaoVarProxy::delta_merit,
          &TaoVarProxy::set_delta_merit)
      // TaoVarProxy.merit (0D_NOT_real - merit_term = weight * delta^2.
      .def_property("merit", &TaoVarProxy::merit, &TaoVarProxy::set_merit)
      // TaoVarProxy.dMerit_dVar (0D_NOT_real - Merit derivative.
      .def_property(
          "dMerit_dVar",
          &TaoVarProxy::dMerit_dVar,
          &TaoVarProxy::set_dMerit_dVar)
      // TaoVarProxy.key_val0 (0D_NOT_real - Key base value
      .def_property(
          "key_val0", &TaoVarProxy::key_val0, &TaoVarProxy::set_key_val0)
      // TaoVarProxy.key_delta (0D_NOT_real - Change in value when a key is pressed.
      .def_property(
          "key_delta", &TaoVarProxy::key_delta, &TaoVarProxy::set_key_delta)
      // TaoVarProxy.s (0D_NOT_real - longitudinal position of ele.
      .def_property("s", &TaoVarProxy::s, &TaoVarProxy::set_s)
      // TaoVarProxy.extend_val (0D_NOT_real - For extension code. Not used by Tao.
      .def_property(
          "extend_val", &TaoVarProxy::extend_val, &TaoVarProxy::set_extend_val)
      // TaoVarProxy.merit_type (0D_NOT_character - 'target' or 'limit'
      .def_property(
          "merit_type", &TaoVarProxy::merit_type, &TaoVarProxy::set_merit_type)
      // TaoVarProxy.exists (0D_NOT_logical - See above
      .def_property("exists", &TaoVarProxy::exists, &TaoVarProxy::set_exists)
      // TaoVarProxy.good_var (0D_NOT_logical - See above
      .def_property(
          "good_var", &TaoVarProxy::good_var, &TaoVarProxy::set_good_var)
      // TaoVarProxy.good_user (0D_NOT_logical - See above
      .def_property(
          "good_user", &TaoVarProxy::good_user, &TaoVarProxy::set_good_user)
      // TaoVarProxy.good_opt (0D_NOT_logical - See above
      .def_property(
          "good_opt", &TaoVarProxy::good_opt, &TaoVarProxy::set_good_opt)
      // TaoVarProxy.good_plot (0D_NOT_logical - See above
      .def_property(
          "good_plot", &TaoVarProxy::good_plot, &TaoVarProxy::set_good_plot)
      // TaoVarProxy.useit_opt (0D_NOT_logical - See above
      .def_property(
          "useit_opt", &TaoVarProxy::useit_opt, &TaoVarProxy::set_useit_opt)
      // TaoVarProxy.useit_plot (0D_NOT_logical - See above
      .def_property(
          "useit_plot", &TaoVarProxy::useit_plot, &TaoVarProxy::set_useit_plot)
      // TaoVarProxy.key_bound (0D_NOT_logical - Variable bound to keyboard key?
      .def_property(
          "key_bound", &TaoVarProxy::key_bound, &TaoVarProxy::set_key_bound)
      // TaoVarProxy.v1 (0D_PTR_type - Pointer to the parent.
      .def_property("v1", &TaoVarProxy::v1, &TaoVarProxy::set_v1)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoVarProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const TaoVarProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoVarProxy& self) {
            return TaoVarProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoVarProxy& self, py::dict& memo) {
            return TaoVarProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoVarProxyArray1D>(m, "TaoVarStructArray1D");
  bind_FTypeAlloc1D<TaoVarProxyAlloc1D>(m, "TaoVarStructAlloc1D");
  // 2D TaoVarProxy arrays are not used in structs/routines
  // 3D TaoVarProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_kick_pt_struct
void init_tao_wave_kick_pt_struct(
    py::module& m,
    py::class_<TaoWaveKickPtProxy>& cls) {
  cls.def(py::init<>())
      // TaoWaveKickPtProxy.phi_s (0D_NOT_real -
      .def_property(
          "phi_s", &TaoWaveKickPtProxy::phi_s, &TaoWaveKickPtProxy::set_phi_s)
      // TaoWaveKickPtProxy.phi_r (0D_NOT_real -
      .def_property(
          "phi_r", &TaoWaveKickPtProxy::phi_r, &TaoWaveKickPtProxy::set_phi_r)
      // TaoWaveKickPtProxy.phi (0D_NOT_real -
      .def_property(
          "phi", &TaoWaveKickPtProxy::phi, &TaoWaveKickPtProxy::set_phi)
      // TaoWaveKickPtProxy.amp (0D_NOT_real -
      .def_property(
          "amp", &TaoWaveKickPtProxy::amp, &TaoWaveKickPtProxy::set_amp)
      // TaoWaveKickPtProxy.s (0D_NOT_real - s-position of kick
      .def_property("s", &TaoWaveKickPtProxy::s, &TaoWaveKickPtProxy::set_s)
      // TaoWaveKickPtProxy.ix_dat_before_kick (0D_NOT_integer - Index of datum in data array just before the kick.
      .def_property(
          "ix_dat_before_kick",
          &TaoWaveKickPtProxy::ix_dat_before_kick,
          &TaoWaveKickPtProxy::set_ix_dat_before_kick)
      // TaoWaveKickPtProxy.ele (0D_PTR_type - lattice element at position of kick.
      .def_property(
          "ele", &TaoWaveKickPtProxy::ele, &TaoWaveKickPtProxy::set_ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoWaveKickPtProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoWaveKickPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoWaveKickPtProxy& self) {
            return TaoWaveKickPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoWaveKickPtProxy& self, py::dict& memo) {
            return TaoWaveKickPtProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoWaveKickPtProxyArray1D>(m, "TaoWaveKickPtStructArray1D");
  bind_FTypeAlloc1D<TaoWaveKickPtProxyAlloc1D>(m, "TaoWaveKickPtStructAlloc1D");
  // 2D TaoWaveKickPtProxy arrays are not used in structs/routines
  // 3D TaoWaveKickPtProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_struct
void init_tao_wave_struct(py::module& m, py::class_<TaoWaveProxy>& cls) {
  cls.def(py::init<>())
      // TaoWaveProxy.data_type (0D_NOT_character -
      .def_property(
          "data_type", &TaoWaveProxy::data_type, &TaoWaveProxy::set_data_type)
      // TaoWaveProxy.rms_rel_a (0D_NOT_real -
      .def_property(
          "rms_rel_a", &TaoWaveProxy::rms_rel_a, &TaoWaveProxy::set_rms_rel_a)
      // TaoWaveProxy.rms_rel_b (0D_NOT_real -
      .def_property(
          "rms_rel_b", &TaoWaveProxy::rms_rel_b, &TaoWaveProxy::set_rms_rel_b)
      // TaoWaveProxy.rms_rel_as (0D_NOT_real -
      .def_property(
          "rms_rel_as",
          &TaoWaveProxy::rms_rel_as,
          &TaoWaveProxy::set_rms_rel_as)
      // TaoWaveProxy.rms_rel_bs (0D_NOT_real -
      .def_property(
          "rms_rel_bs",
          &TaoWaveProxy::rms_rel_bs,
          &TaoWaveProxy::set_rms_rel_bs)
      // TaoWaveProxy.rms_rel_ar (0D_NOT_real -
      .def_property(
          "rms_rel_ar",
          &TaoWaveProxy::rms_rel_ar,
          &TaoWaveProxy::set_rms_rel_ar)
      // TaoWaveProxy.rms_rel_br (0D_NOT_real -
      .def_property(
          "rms_rel_br",
          &TaoWaveProxy::rms_rel_br,
          &TaoWaveProxy::set_rms_rel_br)
      // TaoWaveProxy.rms_rel_k (0D_NOT_real -
      .def_property(
          "rms_rel_k", &TaoWaveProxy::rms_rel_k, &TaoWaveProxy::set_rms_rel_k)
      // TaoWaveProxy.rms_rel_ks (0D_NOT_real -
      .def_property(
          "rms_rel_ks",
          &TaoWaveProxy::rms_rel_ks,
          &TaoWaveProxy::set_rms_rel_ks)
      // TaoWaveProxy.rms_rel_kr (0D_NOT_real -
      .def_property(
          "rms_rel_kr",
          &TaoWaveProxy::rms_rel_kr,
          &TaoWaveProxy::set_rms_rel_kr)
      // TaoWaveProxy.rms_phi (0D_NOT_real -
      .def_property(
          "rms_phi", &TaoWaveProxy::rms_phi, &TaoWaveProxy::set_rms_phi)
      // TaoWaveProxy.rms_phi_s (0D_NOT_real -
      .def_property(
          "rms_phi_s", &TaoWaveProxy::rms_phi_s, &TaoWaveProxy::set_rms_phi_s)
      // TaoWaveProxy.rms_phi_r (0D_NOT_real -
      .def_property(
          "rms_phi_r", &TaoWaveProxy::rms_phi_r, &TaoWaveProxy::set_rms_phi_r)
      // TaoWaveProxy.amp_ba_s (0D_NOT_real -
      .def_property(
          "amp_ba_s", &TaoWaveProxy::amp_ba_s, &TaoWaveProxy::set_amp_ba_s)
      // TaoWaveProxy.amp_ba_r (0D_NOT_real -
      .def_property(
          "amp_ba_r", &TaoWaveProxy::amp_ba_r, &TaoWaveProxy::set_amp_ba_r)
      // TaoWaveProxy.chi_a (0D_NOT_real -
      .def_property("chi_a", &TaoWaveProxy::chi_a, &TaoWaveProxy::set_chi_a)
      // TaoWaveProxy.chi_c (0D_NOT_real -
      .def_property("chi_c", &TaoWaveProxy::chi_c, &TaoWaveProxy::set_chi_c)
      // TaoWaveProxy.chi_ba (0D_NOT_real -
      .def_property("chi_ba", &TaoWaveProxy::chi_ba, &TaoWaveProxy::set_chi_ba)
      // TaoWaveProxy.amp_a (1D_NOT_real -
      .def_property_readonly("amp_a", &TaoWaveProxy::amp_a)
      // TaoWaveProxy.amp_b (1D_NOT_real -
      .def_property_readonly("amp_b", &TaoWaveProxy::amp_b)
      // TaoWaveProxy.amp_ba (1D_NOT_real -
      .def_property_readonly("amp_ba", &TaoWaveProxy::amp_ba)
      // TaoWaveProxy.coef_a (1D_NOT_real -
      .def_property_readonly("coef_a", &TaoWaveProxy::coef_a)
      // TaoWaveProxy.coef_b (1D_NOT_real -
      .def_property_readonly("coef_b", &TaoWaveProxy::coef_b)
      // TaoWaveProxy.coef_ba (1D_NOT_real -
      .def_property_readonly("coef_ba", &TaoWaveProxy::coef_ba)
      // TaoWaveProxy.n_func (0D_NOT_integer - Number of functions used in the fit.
      .def_property("n_func", &TaoWaveProxy::n_func, &TaoWaveProxy::set_n_func)
      // TaoWaveProxy.ix_a1 (0D_NOT_integer -
      .def_property("ix_a1", &TaoWaveProxy::ix_a1, &TaoWaveProxy::set_ix_a1)
      // TaoWaveProxy.ix_a2 (0D_NOT_integer -
      .def_property("ix_a2", &TaoWaveProxy::ix_a2, &TaoWaveProxy::set_ix_a2)
      // TaoWaveProxy.ix_b1 (0D_NOT_integer -
      .def_property("ix_b1", &TaoWaveProxy::ix_b1, &TaoWaveProxy::set_ix_b1)
      // TaoWaveProxy.ix_b2 (0D_NOT_integer -
      .def_property("ix_b2", &TaoWaveProxy::ix_b2, &TaoWaveProxy::set_ix_b2)
      // TaoWaveProxy.i_a1 (0D_NOT_integer -
      .def_property("i_a1", &TaoWaveProxy::i_a1, &TaoWaveProxy::set_i_a1)
      // TaoWaveProxy.i_a2 (0D_NOT_integer -
      .def_property("i_a2", &TaoWaveProxy::i_a2, &TaoWaveProxy::set_i_a2)
      // TaoWaveProxy.i_b1 (0D_NOT_integer -
      .def_property("i_b1", &TaoWaveProxy::i_b1, &TaoWaveProxy::set_i_b1)
      // TaoWaveProxy.i_b2 (0D_NOT_integer -
      .def_property("i_b2", &TaoWaveProxy::i_b2, &TaoWaveProxy::set_i_b2)
      // TaoWaveProxy.n_a (0D_NOT_integer -
      .def_property("n_a", &TaoWaveProxy::n_a, &TaoWaveProxy::set_n_a)
      // TaoWaveProxy.n_b (0D_NOT_integer -
      .def_property("n_b", &TaoWaveProxy::n_b, &TaoWaveProxy::set_n_b)
      // TaoWaveProxy.i_curve_wrap_pt (0D_NOT_integer - Index of last point before wrap in curve array.
      .def_property(
          "i_curve_wrap_pt",
          &TaoWaveProxy::i_curve_wrap_pt,
          &TaoWaveProxy::set_i_curve_wrap_pt)
      // TaoWaveProxy.ix_data (1D_ALLOC_integer - Translates from plot point to datum index
      .def_property_readonly("ix_data", &TaoWaveProxy::ix_data)
      // TaoWaveProxy.n_kick (0D_NOT_integer -
      .def_property("n_kick", &TaoWaveProxy::n_kick, &TaoWaveProxy::set_n_kick)
      // TaoWaveProxy.kick (1D_ALLOC_type -
      .def_property_readonly("kick", &TaoWaveProxy::kick)
      // TaoWaveProxy.base_graph (0D_NOT_type - Graph before curves extended to 1.5 periods.
      .def_property(
          "base_graph",
          &TaoWaveProxy::base_graph,
          &TaoWaveProxy::set_base_graph)
      // TaoWaveProxy.region (0D_PTR_type - Where the wave plot is
      .def_property("region", &TaoWaveProxy::region, &TaoWaveProxy::set_region)
      // TaoWaveProxy.d1_dat (0D_PTR_type - D1 data for analysis
      .def_property("d1_dat", &TaoWaveProxy::d1_dat, &TaoWaveProxy::set_d1_dat)

      .def("__repr__", [](const TaoWaveProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoWaveProxy& self) {
            return TaoWaveProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoWaveProxy& self, py::dict& memo) {
            return TaoWaveProxy(self);
          })

      ;

  // 1D TaoWaveProxy arrays are not used in structs/routines
  // 2D TaoWaveProxy arrays are not used in structs/routines
  // 3D TaoWaveProxy arrays are not used in structs/routines
}

// =============================================================================
// all_encompassing_struct
void init_all_encompassing_struct(
    py::module& m,
    py::class_<AllEncompassingProxy>& cls) {
  cls.def(py::init<>())
      // AllEncompassingProxy.real_rp_0d (0D_NOT_real -
      .def_property(
          "real_rp_0d",
          &AllEncompassingProxy::real_rp_0d,
          &AllEncompassingProxy::set_real_rp_0d)
      // AllEncompassingProxy.real_rp_1d (1D_NOT_real -
      .def_property_readonly("real_rp_1d", &AllEncompassingProxy::real_rp_1d)
      // AllEncompassingProxy.real_rp_2d (2D_NOT_real -
      .def_property_readonly("real_rp_2d", &AllEncompassingProxy::real_rp_2d)
      // AllEncompassingProxy.real_rp_3d (3D_NOT_real -
      .def_property_readonly("real_rp_3d", &AllEncompassingProxy::real_rp_3d)
      // AllEncompassingProxy.real_rp_0d_ptr (0D_PTR_real -
      .def_property(
          "real_rp_0d_ptr",
          &AllEncompassingProxy::real_rp_0d_ptr,
          &AllEncompassingProxy::set_real_rp_0d_ptr)
      // AllEncompassingProxy.real_rp_1d_ptr (1D_PTR_real -
      .def_property_readonly(
          "real_rp_1d_ptr", &AllEncompassingProxy::real_rp_1d_ptr)
      // AllEncompassingProxy.real_rp_2d_ptr (2D_PTR_real -
      .def_property_readonly(
          "real_rp_2d_ptr", &AllEncompassingProxy::real_rp_2d_ptr)
      // AllEncompassingProxy.real_rp_3d_ptr (3D_PTR_real -
      .def_property_readonly(
          "real_rp_3d_ptr", &AllEncompassingProxy::real_rp_3d_ptr)
      // AllEncompassingProxy.real_rp_1d_alloc (1D_ALLOC_real -
      .def_property_readonly(
          "real_rp_1d_alloc", &AllEncompassingProxy::real_rp_1d_alloc)
      // AllEncompassingProxy.real_rp_2d_alloc (2D_ALLOC_real -
      .def_property_readonly(
          "real_rp_2d_alloc", &AllEncompassingProxy::real_rp_2d_alloc)
      // AllEncompassingProxy.real_rp_3d_alloc (3D_ALLOC_real - Real(dp)
      .def_property_readonly(
          "real_rp_3d_alloc", &AllEncompassingProxy::real_rp_3d_alloc)
      // AllEncompassingProxy.real_dp_0d (0D_NOT_real -
      .def_property(
          "real_dp_0d",
          &AllEncompassingProxy::real_dp_0d,
          &AllEncompassingProxy::set_real_dp_0d)
      // AllEncompassingProxy.real_dp_1d (1D_NOT_real -
      .def_property_readonly("real_dp_1d", &AllEncompassingProxy::real_dp_1d)
      // AllEncompassingProxy.real_dp_2d (2D_NOT_real -
      .def_property_readonly("real_dp_2d", &AllEncompassingProxy::real_dp_2d)
      // AllEncompassingProxy.real_dp_3d (3D_NOT_real -
      .def_property_readonly("real_dp_3d", &AllEncompassingProxy::real_dp_3d)
      // AllEncompassingProxy.real_dp_0d_ptr (0D_PTR_real -
      .def_property(
          "real_dp_0d_ptr",
          &AllEncompassingProxy::real_dp_0d_ptr,
          &AllEncompassingProxy::set_real_dp_0d_ptr)
      // AllEncompassingProxy.real_dp_1d_ptr (1D_PTR_real -
      .def_property_readonly(
          "real_dp_1d_ptr", &AllEncompassingProxy::real_dp_1d_ptr)
      // AllEncompassingProxy.real_dp_2d_ptr (2D_PTR_real -
      .def_property_readonly(
          "real_dp_2d_ptr", &AllEncompassingProxy::real_dp_2d_ptr)
      // AllEncompassingProxy.real_dp_3d_ptr (3D_PTR_real -
      .def_property_readonly(
          "real_dp_3d_ptr", &AllEncompassingProxy::real_dp_3d_ptr)
      // AllEncompassingProxy.real_dp_1d_alloc (1D_ALLOC_real -
      .def_property_readonly(
          "real_dp_1d_alloc", &AllEncompassingProxy::real_dp_1d_alloc)
      // AllEncompassingProxy.real_dp_2d_alloc (2D_ALLOC_real -
      .def_property_readonly(
          "real_dp_2d_alloc", &AllEncompassingProxy::real_dp_2d_alloc)
      // AllEncompassingProxy.real_dp_3d_alloc (3D_ALLOC_real - complex(dp)
      .def_property_readonly(
          "real_dp_3d_alloc", &AllEncompassingProxy::real_dp_3d_alloc)
      // AllEncompassingProxy.complex_dp_0d (0D_NOT_complex -
      .def_property(
          "complex_dp_0d",
          &AllEncompassingProxy::complex_dp_0d,
          &AllEncompassingProxy::set_complex_dp_0d)
      // AllEncompassingProxy.complex_dp_1d (1D_NOT_complex -
      .def_property_readonly(
          "complex_dp_1d", &AllEncompassingProxy::complex_dp_1d)
      // AllEncompassingProxy.complex_dp_2d (2D_NOT_complex -
      .def_property_readonly(
          "complex_dp_2d", &AllEncompassingProxy::complex_dp_2d)
      // AllEncompassingProxy.complex_dp_3d (3D_NOT_complex - TODO complex(dp), pointer :: complex_dp_0d_ptr
      .def_property_readonly(
          "complex_dp_3d", &AllEncompassingProxy::complex_dp_3d)
      // AllEncompassingProxy.complex_dp_1d_ptr (1D_PTR_complex -
      .def_property_readonly(
          "complex_dp_1d_ptr", &AllEncompassingProxy::complex_dp_1d_ptr)
      // AllEncompassingProxy.complex_dp_2d_ptr (2D_PTR_complex -
      .def_property_readonly(
          "complex_dp_2d_ptr", &AllEncompassingProxy::complex_dp_2d_ptr)
      // AllEncompassingProxy.complex_dp_3d_ptr (3D_PTR_complex -
      .def_property_readonly(
          "complex_dp_3d_ptr", &AllEncompassingProxy::complex_dp_3d_ptr)
      // AllEncompassingProxy.complex_dp_1d_alloc (1D_ALLOC_complex -
      .def_property_readonly(
          "complex_dp_1d_alloc", &AllEncompassingProxy::complex_dp_1d_alloc)
      // AllEncompassingProxy.complex_dp_2d_alloc (2D_ALLOC_complex -
      .def_property_readonly(
          "complex_dp_2d_alloc", &AllEncompassingProxy::complex_dp_2d_alloc)
      // AllEncompassingProxy.complex_dp_3d_alloc (3D_ALLOC_complex - Integer
      .def_property_readonly(
          "complex_dp_3d_alloc", &AllEncompassingProxy::complex_dp_3d_alloc)
      // AllEncompassingProxy.int_0d (0D_NOT_integer -
      .def_property(
          "int_0d",
          &AllEncompassingProxy::int_0d,
          &AllEncompassingProxy::set_int_0d)
      // AllEncompassingProxy.int_1d (1D_NOT_integer -
      .def_property_readonly("int_1d", &AllEncompassingProxy::int_1d)
      // AllEncompassingProxy.int_2d (2D_NOT_integer -
      .def_property_readonly("int_2d", &AllEncompassingProxy::int_2d)
      // AllEncompassingProxy.int_3d (3D_NOT_integer -
      .def_property_readonly("int_3d", &AllEncompassingProxy::int_3d)
      // AllEncompassingProxy.int_0d_ptr (0D_PTR_integer -
      .def_property(
          "int_0d_ptr",
          &AllEncompassingProxy::int_0d_ptr,
          &AllEncompassingProxy::set_int_0d_ptr)
      // AllEncompassingProxy.int_1d_ptr (1D_PTR_integer -
      .def_property_readonly("int_1d_ptr", &AllEncompassingProxy::int_1d_ptr)
      // AllEncompassingProxy.int_2d_ptr (2D_PTR_integer -
      .def_property_readonly("int_2d_ptr", &AllEncompassingProxy::int_2d_ptr)
      // AllEncompassingProxy.int_3d_ptr (3D_PTR_integer -
      .def_property_readonly("int_3d_ptr", &AllEncompassingProxy::int_3d_ptr)
      // AllEncompassingProxy.int_1d_alloc (1D_ALLOC_integer -
      .def_property_readonly(
          "int_1d_alloc", &AllEncompassingProxy::int_1d_alloc)
      // AllEncompassingProxy.int_2d_alloc (2D_ALLOC_integer -
      .def_property_readonly(
          "int_2d_alloc", &AllEncompassingProxy::int_2d_alloc)
      // AllEncompassingProxy.int_3d_alloc (3D_ALLOC_integer - Integer8
      .def_property_readonly(
          "int_3d_alloc", &AllEncompassingProxy::int_3d_alloc)
      // AllEncompassingProxy.int8_0d (0D_NOT_integer8 -
      .def_property(
          "int8_0d",
          &AllEncompassingProxy::int8_0d,
          &AllEncompassingProxy::set_int8_0d)
      // 1D_NOT_integer8 int8_1d proxy support missing
      // 2D_NOT_integer8 int8_2d proxy support missing
      // 3D_NOT_integer8 int8_3d proxy support missing
      // AllEncompassingProxy.int8_0d_ptr (0D_PTR_integer8 -
      .def_property(
          "int8_0d_ptr",
          &AllEncompassingProxy::int8_0d_ptr,
          &AllEncompassingProxy::set_int8_0d_ptr)
      // 1D_PTR_integer8 int8_1d_ptr proxy support missing
      // 2D_PTR_integer8 int8_2d_ptr proxy support missing
      // 3D_PTR_integer8 int8_3d_ptr proxy support missing
      // 1D_ALLOC_integer8 int8_1d_alloc proxy support missing
      // 2D_ALLOC_integer8 int8_2d_alloc proxy support missing
      // 3D_ALLOC_integer8 int8_3d_alloc proxy support missing
      // AllEncompassingProxy.logical_0d (0D_NOT_logical -
      .def_property(
          "logical_0d",
          &AllEncompassingProxy::logical_0d,
          &AllEncompassingProxy::set_logical_0d)
      // 1D_NOT_logical logical_1d proxy support missing
      // 2D_NOT_logical logical_2d proxy support missing
      // 3D_NOT_logical logical_3d proxy support missing
      // AllEncompassingProxy.logical_0d_ptr (0D_PTR_logical - logical, pointer :: logical_1d_ptr(:) logical, pointer :: logical_2d_ptr(:,:) logical, pointer :: logical_3d_ptr(:,:,:) logical, allocatable :: logical_1d_alloc(:) logical, allocatable :: logical_2d_alloc(:,:) logical, allocatable :: logical_3d_alloc(:,:,:) type
      .def_property(
          "logical_0d_ptr",
          &AllEncompassingProxy::logical_0d_ptr,
          &AllEncompassingProxy::set_logical_0d_ptr)
      // AllEncompassingProxy.type_0d (0D_NOT_type -
      .def_property(
          "type_0d",
          &AllEncompassingProxy::type_0d,
          &AllEncompassingProxy::set_type_0d)
      // AllEncompassingProxy.type_1d (1D_NOT_type -
      .def_property_readonly("type_1d", &AllEncompassingProxy::type_1d)
      // AllEncompassingProxy.type_2d (2D_NOT_type -
      .def_property_readonly("type_2d", &AllEncompassingProxy::type_2d)
      // AllEncompassingProxy.type_3d (3D_NOT_type -
      .def_property_readonly("type_3d", &AllEncompassingProxy::type_3d)
      // AllEncompassingProxy.type_0d_ptr (0D_PTR_type -
      .def_property(
          "type_0d_ptr",
          &AllEncompassingProxy::type_0d_ptr,
          &AllEncompassingProxy::set_type_0d_ptr)
      // AllEncompassingProxy.type_1d_ptr (1D_PTR_type -
      .def_property_readonly("type_1d_ptr", &AllEncompassingProxy::type_1d_ptr)
      // AllEncompassingProxy.type_2d_ptr (2D_PTR_type -
      .def_property_readonly("type_2d_ptr", &AllEncompassingProxy::type_2d_ptr)
      // AllEncompassingProxy.type_3d_ptr (3D_PTR_type -
      .def_property_readonly("type_3d_ptr", &AllEncompassingProxy::type_3d_ptr)
      // AllEncompassingProxy.type_1d_alloc (1D_ALLOC_type -
      .def_property_readonly(
          "type_1d_alloc", &AllEncompassingProxy::type_1d_alloc)
      // AllEncompassingProxy.type_2d_alloc (2D_ALLOC_type -
      .def_property_readonly(
          "type_2d_alloc", &AllEncompassingProxy::type_2d_alloc)
      // AllEncompassingProxy.type_3d_alloc (3D_ALLOC_type -
      .def_property_readonly(
          "type_3d_alloc", &AllEncompassingProxy::type_3d_alloc)

      .def(
          "__repr__",
          [](const AllEncompassingProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AllEncompassingProxy& self) {
            return AllEncompassingProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AllEncompassingProxy& self, py::dict& memo) {
            return AllEncompassingProxy(self);
          })

      ;

  // 1D AllEncompassingProxy arrays are not used in structs/routines
  // 2D AllEncompassingProxy arrays are not used in structs/routines
  // 3D AllEncompassingProxy arrays are not used in structs/routines
}

// =============================================================================
// test_sub_struct
void init_test_sub_struct(py::module& m, py::class_<TestSubProxy>& cls) {
  cls.def(py::init<>())
      // TestSubProxy.sr (0D_NOT_type -
      .def_property("sr", &TestSubProxy::sr, &TestSubProxy::set_sr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TestSubProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const TestSubProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TestSubProxy& self) {
            return TestSubProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TestSubProxy& self, py::dict& memo) {
            return TestSubProxy(self);
          })

      ;

  bind_FTypeArrayND<TestSubProxyArray1D>(m, "TestSubStructArray1D");
  bind_FTypeAlloc1D<TestSubProxyAlloc1D>(m, "TestSubStructAlloc1D");
  bind_FTypeArrayND<TestSubProxyArray2D>(m, "TestSubStructArray2D");
  bind_FTypeArrayND<TestSubProxyArray3D>(m, "TestSubStructArray3D");
}

// =============================================================================
// test_sub_sub_struct
void init_test_sub_sub_struct(py::module& m, py::class_<TestSubSubProxy>& cls) {
  cls.def(py::init<>())
      // TestSubSubProxy.aaa (0D_NOT_integer8 -
      .def_property("aaa", &TestSubSubProxy::aaa, &TestSubSubProxy::set_aaa)
      // TestSubSubProxy.bbb (0D_NOT_integer -
      .def_property("bbb", &TestSubSubProxy::bbb, &TestSubSubProxy::set_bbb)
      // TestSubSubProxy.file (0D_NOT_character -
      .def_property("file", &TestSubSubProxy::file, &TestSubSubProxy::set_file)
      // TestSubSubProxy.t_ref (0D_NOT_real - time reference value for computing the wake amplitude. This is used to prevent value overflow with long trains.
      .def_property(
          "t_ref", &TestSubSubProxy::t_ref, &TestSubSubProxy::set_t_ref)
      // TestSubSubProxy.freq_spread (0D_NOT_real - Random frequency spread of long range modes.
      .def_property(
          "freq_spread",
          &TestSubSubProxy::freq_spread,
          &TestSubSubProxy::set_freq_spread)

      .def(
          "__repr__",
          [](const TestSubSubProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TestSubSubProxy& self) {
            return TestSubSubProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TestSubSubProxy& self, py::dict& memo) {
            return TestSubSubProxy(self);
          })

      ;

  // 1D TestSubSubProxy arrays are not used in structs/routines
  // 2D TestSubSubProxy arrays are not used in structs/routines
  // 3D TestSubSubProxy arrays are not used in structs/routines
}
} // namespace Pybmad