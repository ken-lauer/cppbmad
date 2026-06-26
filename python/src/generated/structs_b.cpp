#include "pybmad/generated/structs_b.hpp"

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
// bbu_beam_struct
void init_bbu_beam_struct(nb::module_ &m, nb::class_<BbuBeamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<int>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("ix_ele_bunch") = nb::none(),
         nb::arg("ix_bunch_head") = nb::none(),
         nb::arg("ix_bunch_end") = nb::none(),
         nb::arg("n_bunch_in_lat") = nb::none(),
         nb::arg("ix_stage_voltage_max") = nb::none(),
         nb::arg("hom_voltage_max") = nb::none(),
         nb::arg("time_now") = nb::none(),
         nb::arg("one_turn_time") = nb::none(),
         nb::arg("rf_wavelength_max") = nb::none()
  )
      .def_prop_ro("bunch", &BbuBeamStruct::bunch, nb::keep_alive<0, 1>(), "Bunches in the lattice")
      .def_prop_ro("stage", &BbuBeamStruct::stage, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "ix_ele_bunch",
          &BbuBeamStruct::ix_ele_bunch,
          &BbuBeamStruct::set_ix_ele_bunch,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "element where bunch is"
      )
      .def_prop_rw(
          "ix_bunch_head",
          &BbuBeamStruct::ix_bunch_head,
          &BbuBeamStruct::set_ix_bunch_head,
          "Index to head bunch(:)"
      )
      .def_prop_rw(
          "ix_bunch_end",
          &BbuBeamStruct::ix_bunch_end,
          &BbuBeamStruct::set_ix_bunch_end,
          "Index of the end bunch(:). -1 -> no bunches."
      )
      .def_prop_rw(
          "n_bunch_in_lat",
          &BbuBeamStruct::n_bunch_in_lat,
          &BbuBeamStruct::set_n_bunch_in_lat,
          "Number of bunches transversing the lattice."
      )
      .def_prop_rw(
          "ix_stage_voltage_max",
          &BbuBeamStruct::ix_stage_voltage_max,
          &BbuBeamStruct::set_ix_stage_voltage_max
      )
      .def_prop_rw(
          "hom_voltage_max",
          &BbuBeamStruct::hom_voltage_max,
          &BbuBeamStruct::set_hom_voltage_max
      )
      .def_prop_rw("time_now", &BbuBeamStruct::time_now, &BbuBeamStruct::set_time_now)
      .def_prop_rw(
          "one_turn_time",
          &BbuBeamStruct::one_turn_time,
          &BbuBeamStruct::set_one_turn_time
      )
      .def_prop_rw(
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
          [](const BbuBeamStruct &self, nb::dict &memo) { return BbuBeamStruct(self); }
      )
      .def(
          "__eq__",
          [](const BbuBeamStruct &self, const BbuBeamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BbuBeamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BbuBeamStruct arrays are not used in structs/routines
  // 2D BbuBeamStruct arrays are not used in structs/routines
  // 3D BbuBeamStruct arrays are not used in structs/routines
}

// =============================================================================
// bbu_param_struct
void init_bbu_param_struct(nb::module_ &m, nb::class_<BbuParamStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("lat_filename") = nb::none(),
         nb::arg("lat2_filename") = nb::none(),
         nb::arg("bunch_by_bunch_info_file") = nb::none(),
         nb::arg("hybridize") = nb::none(),
         nb::arg("write_digested_hybrid_lat") = nb::none(),
         nb::arg("write_voltage_vs_time_dat") = nb::none(),
         nb::arg("keep_overlays_and_groups") = nb::none(),
         nb::arg("keep_all_lcavities") = nb::none(),
         nb::arg("use_taylor_for_hybrids") = nb::none(),
         nb::arg("stable_orbit_anal") = nb::none(),
         nb::arg("limit_factor") = nb::none(),
         nb::arg("simulation_turns_max") = nb::none(),
         nb::arg("bunch_freq") = nb::none(),
         nb::arg("init_particle_offset") = nb::none(),
         nb::arg("current") = nb::none(),
         nb::arg("rel_tol") = nb::none(),
         nb::arg("drscan") = nb::none(),
         nb::arg("use_interpolated_threshold") = nb::none(),
         nb::arg("write_hom_info") = nb::none(),
         nb::arg("elindex") = nb::none(),
         nb::arg("elname") = nb::none(),
         nb::arg("nstep") = nb::none(),
         nb::arg("begdr") = nb::none(),
         nb::arg("enddr") = nb::none(),
         nb::arg("nrep") = nb::none(),
         nb::arg("ran_seed") = nb::none(),
         nb::arg("hom_order_cutoff") = nb::none(),
         nb::arg("ran_gauss_sigma_cut") = nb::none(),
         nb::arg("ele_track_end") = nb::none(),
         nb::arg("ix_ele_track_end") = nb::none(),
         nb::arg("regression") = nb::none(),
         nb::arg("normalize_z_to_rf") = nb::none(),
         nb::arg("ramp_on") = nb::none(),
         nb::arg("ramp_pattern") = nb::none(),
         nb::arg("ramp_n_start") = nb::none(),
         nb::arg("n_ramp_pattern") = nb::none()
  )
      .def_prop_rw(
          "lat_filename",
          &BbuParamStruct::lat_filename,
          &BbuParamStruct::set_lat_filename,
          "Bmad lattice file name"
      )
      .def_prop_rw(
          "lat2_filename",
          &BbuParamStruct::lat2_filename,
          &BbuParamStruct::set_lat2_filename,
          "Bmad lattice2 file name for secondary parser"
      )
      .def_prop_rw(
          "bunch_by_bunch_info_file",
          &BbuParamStruct::bunch_by_bunch_info_file,
          &BbuParamStruct::set_bunch_by_bunch_info_file,
          "For outputting bunch-by-bunch info."
      )
      .def_prop_rw(
          "hybridize",
          &BbuParamStruct::hybridize,
          &BbuParamStruct::set_hybridize,
          "Combine non-hom elements to speed up simulation?"
      )
      .def_prop_rw(
          "write_digested_hybrid_lat",
          &BbuParamStruct::write_digested_hybrid_lat,
          &BbuParamStruct::set_write_digested_hybrid_lat,
          "For debugging purposes."
      )
      .def_prop_rw(
          "write_voltage_vs_time_dat",
          &BbuParamStruct::write_voltage_vs_time_dat,
          &BbuParamStruct::set_write_voltage_vs_time_dat,
          "For debugging purposes."
      )
      .def_prop_rw(
          "keep_overlays_and_groups",
          &BbuParamStruct::keep_overlays_and_groups,
          &BbuParamStruct::set_keep_overlays_and_groups,
          "Keep when hybridizing?"
      )
      .def_prop_rw(
          "keep_all_lcavities",
          &BbuParamStruct::keep_all_lcavities,
          &BbuParamStruct::set_keep_all_lcavities,
          "Keep when hybridizing?"
      )
      .def_prop_rw(
          "use_taylor_for_hybrids",
          &BbuParamStruct::use_taylor_for_hybrids,
          &BbuParamStruct::set_use_taylor_for_hybrids,
          "Use taylor map for hybrids when true. Otherwise tracking method is linear."
      )
      .def_prop_rw(
          "stable_orbit_anal",
          &BbuParamStruct::stable_orbit_anal,
          &BbuParamStruct::set_stable_orbit_anal,
          "Write stable_orbit.out and hom_voltage.out?"
      )
      .def_prop_rw(
          "limit_factor",
          &BbuParamStruct::limit_factor,
          &BbuParamStruct::set_limit_factor,
          "Init_hom_amp * limit_factor = simulation unstable limit"
      )
      .def_prop_rw(
          "simulation_turns_max",
          &BbuParamStruct::simulation_turns_max,
          &BbuParamStruct::set_simulation_turns_max,
          "Sets the duration of the simulation."
      )
      .def_prop_rw(
          "bunch_freq",
          &BbuParamStruct::bunch_freq,
          &BbuParamStruct::set_bunch_freq,
          "Freq in Hz."
      )
      .def_prop_rw(
          "init_particle_offset",
          &BbuParamStruct::init_particle_offset,
          &BbuParamStruct::set_init_particle_offset,
          "Initial particle offset for particles born in the first turn period."
      )
      .def_prop_rw(
          "current",
          &BbuParamStruct::current,
          &BbuParamStruct::set_current,
          "Starting current (amps)"
      )
      .def_prop_rw(
          "rel_tol",
          &BbuParamStruct::rel_tol,
          &BbuParamStruct::set_rel_tol,
          "Final threshold current accuracy."
      )
      .def_prop_rw(
          "drscan",
          &BbuParamStruct::drscan,
          &BbuParamStruct::set_drscan,
          "If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3."
      )
      .def_prop_rw(
          "use_interpolated_threshold",
          &BbuParamStruct::use_interpolated_threshold,
          &BbuParamStruct::set_use_interpolated_threshold
      )
      .def_prop_rw(
          "write_hom_info",
          &BbuParamStruct::write_hom_info,
          &BbuParamStruct::set_write_hom_info,
          "Write HOM parameters to main output file?"
      )
      .def_prop_rw("elindex", &BbuParamStruct::elindex, &BbuParamStruct::set_elindex)
      .def_prop_rw(
          "elname",
          &BbuParamStruct::elname,
          &BbuParamStruct::set_elname,
          "Element to step length for DRSCAN"
      )
      .def_prop_rw(
          "nstep",
          &BbuParamStruct::nstep,
          &BbuParamStruct::set_nstep,
          "Number of steps for DRSCAN."
      )
      .def_prop_rw(
          "begdr",
          &BbuParamStruct::begdr,
          &BbuParamStruct::set_begdr,
          "Beginning DR value for DRSCAN."
      )
      .def_prop_rw(
          "enddr",
          &BbuParamStruct::enddr,
          &BbuParamStruct::set_enddr,
          "End DR value for DRSCAN."
      )
      .def_prop_rw(
          "nrep",
          &BbuParamStruct::nrep,
          &BbuParamStruct::set_nrep,
          "Number of times to repeat threshold calculation"
      )
      .def_prop_rw(
          "ran_seed",
          &BbuParamStruct::ran_seed,
          &BbuParamStruct::set_ran_seed,
          "If set to 0, the output results will vary from run to run."
      )
      .def_prop_rw(
          "hom_order_cutoff",
          &BbuParamStruct::hom_order_cutoff,
          &BbuParamStruct::set_hom_order_cutoff,
          "If positive -> ignore HOM's with order greater than this."
      )
      .def_prop_rw(
          "ran_gauss_sigma_cut",
          &BbuParamStruct::ran_gauss_sigma_cut,
          &BbuParamStruct::set_ran_gauss_sigma_cut
      )
      .def_prop_rw(
          "ele_track_end",
          &BbuParamStruct::ele_track_end,
          &BbuParamStruct::set_ele_track_end
      )
      .def_prop_rw(
          "ix_ele_track_end",
          &BbuParamStruct::ix_ele_track_end,
          &BbuParamStruct::set_ix_ele_track_end,
          "Default: set to last element with a wake"
      )
      .def_prop_rw(
          "regression",
          &BbuParamStruct::regression,
          &BbuParamStruct::set_regression,
          "Do regression test?"
      )
      .def_prop_rw(
          "normalize_z_to_rf",
          &BbuParamStruct::normalize_z_to_rf,
          &BbuParamStruct::set_normalize_z_to_rf,
          "make starting z = mod(z, rf_wavelength)? Ramp parameters"
      )
      .def_prop_rw("ramp_on", &BbuParamStruct::ramp_on, &BbuParamStruct::set_ramp_on)
      .def_prop_rw(
          "ramp_pattern",
          &BbuParamStruct::ramp_pattern,
          &BbuParamStruct::set_ramp_pattern,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ramp_n_start",
          &BbuParamStruct::ramp_n_start,
          &BbuParamStruct::set_ramp_n_start,
          "Index of start of ramp Internal parameters"
      )
      .def_prop_rw(
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
          [](const BbuParamStruct &self, nb::dict &memo) { return BbuParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const BbuParamStruct &self, const BbuParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BbuParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BbuParamStruct arrays are not used in structs/routines
  // 2D BbuParamStruct arrays are not used in structs/routines
  // 3D BbuParamStruct arrays are not used in structs/routines
}

// =============================================================================
// bbu_stage_struct
void init_bbu_stage_struct(nb::module_ &m, nb::class_<BbuStageStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("ix_ele_lr_wake") = nb::none(),
         nb::arg("ix_ele_stage_end") = nb::none(),
         nb::arg("ix_pass") = nb::none(),
         nb::arg("ix_stage_pass1") = nb::none(),
         nb::arg("ix_head_bunch") = nb::none(),
         nb::arg("ix_hom_max") = nb::none(),
         nb::arg("hom_voltage_max") = nb::none(),
         nb::arg("time_at_wake_ele") = nb::none(),
         nb::arg("ave_orb") = nb::none(),
         nb::arg("rms_orb") = nb::none(),
         nb::arg("min_orb") = nb::none(),
         nb::arg("max_orb") = nb::none(),
         nb::arg("n_orb") = nb::none()
  )
      .def_prop_rw(
          "ix_ele_lr_wake",
          &BbuStageStruct::ix_ele_lr_wake,
          &BbuStageStruct::set_ix_ele_lr_wake,
          "Element index of element with the wake"
      )
      .def_prop_rw(
          "ix_ele_stage_end",
          &BbuStageStruct::ix_ele_stage_end,
          &BbuStageStruct::set_ix_ele_stage_end,
          "Element at end of stage."
      )
      .def_prop_rw(
          "ix_pass",
          &BbuStageStruct::ix_pass,
          &BbuStageStruct::set_ix_pass,
          "Pass index when in multipass section"
      )
      .def_prop_rw(
          "ix_stage_pass1",
          &BbuStageStruct::ix_stage_pass1,
          &BbuStageStruct::set_ix_stage_pass1,
          "Index of corresponding stage on first pass"
      )
      .def_prop_rw(
          "ix_head_bunch",
          &BbuStageStruct::ix_head_bunch,
          &BbuStageStruct::set_ix_head_bunch
      )
      .def_prop_rw("ix_hom_max", &BbuStageStruct::ix_hom_max, &BbuStageStruct::set_ix_hom_max)
      .def_prop_rw(
          "hom_voltage_max",
          &BbuStageStruct::hom_voltage_max,
          &BbuStageStruct::set_hom_voltage_max
      )
      .def_prop_rw(
          "time_at_wake_ele",
          &BbuStageStruct::time_at_wake_ele,
          &BbuStageStruct::set_time_at_wake_ele
      )
      .def_prop_rw(
          "ave_orb",
          &BbuStageStruct::ave_orb,
          &BbuStageStruct::set_ave_orb,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "rms_orb",
          &BbuStageStruct::rms_orb,
          &BbuStageStruct::set_rms_orb,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "min_orb",
          &BbuStageStruct::min_orb,
          &BbuStageStruct::set_min_orb,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "max_orb",
          &BbuStageStruct::max_orb,
          &BbuStageStruct::set_max_orb,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("n_orb", &BbuStageStruct::n_orb, &BbuStageStruct::set_n_orb)
      .def_static(
          "new_array1d",
          [](int sz) { return BbuStageStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BbuStageStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const BbuStageStruct &self, nb::dict &memo) { return BbuStageStruct(self); }
      )
      .def(
          "__eq__",
          [](const BbuStageStruct &self, const BbuStageStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BbuStageStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_beam_init_struct(nb::module_ &m, nb::class_<BeamInitStruct> &cls) {
  cls.def(
         "__init__",
         [](BeamInitStruct *self,
            std::optional<std::string> position_file,
            std::optional<std::vector<double>> spin,
            const KvBeamInitStruct *KV,
            std::optional<std::vector<double>> center_jitter,
            std::optional<std::vector<double>> emit_jitter,
            std::optional<double> sig_z_jitter,
            std::optional<double> sig_pz_jitter,
            std::optional<int> n_particle,
            std::optional<bool> renorm_center,
            std::optional<bool> renorm_sigma,
            std::optional<std::string> random_engine,
            std::optional<std::string> random_gauss_converter,
            std::optional<double> random_sigma_cutoff,
            std::optional<double> a_norm_emit,
            std::optional<double> b_norm_emit,
            std::optional<double> a_emit,
            std::optional<double> b_emit,
            std::optional<double> dPz_dz,
            std::optional<std::vector<double>> center,
            std::optional<double> t_offset,
            std::optional<double> dt_bunch,
            std::optional<double> sig_z,
            std::optional<double> sig_pz,
            std::optional<double> bunch_charge,
            std::optional<int> n_bunch,
            std::optional<int> ix_turn,
            std::optional<std::string> species,
            std::optional<bool> full_6D_coupling_calc,
            std::optional<bool> use_particle_start,
            std::optional<bool> use_t_coords,
            std::optional<bool> use_z_as_t,
            std::optional<std::string> file_name) {
           new (self) BeamInitStruct(
               position_file,
               spin,
               ptr_to_opt_ref(KV),
               center_jitter,
               emit_jitter,
               sig_z_jitter,
               sig_pz_jitter,
               n_particle,
               renorm_center,
               renorm_sigma,
               random_engine,
               random_gauss_converter,
               random_sigma_cutoff,
               a_norm_emit,
               b_norm_emit,
               a_emit,
               b_emit,
               dPz_dz,
               center,
               t_offset,
               dt_bunch,
               sig_z,
               sig_pz,
               bunch_charge,
               n_bunch,
               ix_turn,
               species,
               full_6D_coupling_calc,
               use_particle_start,
               use_t_coords,
               use_z_as_t,
               file_name
           );
         },
         nb::arg("position_file") = nb::none(),
         nb::arg("spin") = nb::none(),
         nb::arg("KV") = nb::none(),
         nb::arg("center_jitter") = nb::none(),
         nb::arg("emit_jitter") = nb::none(),
         nb::arg("sig_z_jitter") = nb::none(),
         nb::arg("sig_pz_jitter") = nb::none(),
         nb::arg("n_particle") = nb::none(),
         nb::arg("renorm_center") = nb::none(),
         nb::arg("renorm_sigma") = nb::none(),
         nb::arg("random_engine") = nb::none(),
         nb::arg("random_gauss_converter") = nb::none(),
         nb::arg("random_sigma_cutoff") = nb::none(),
         nb::arg("a_norm_emit") = nb::none(),
         nb::arg("b_norm_emit") = nb::none(),
         nb::arg("a_emit") = nb::none(),
         nb::arg("b_emit") = nb::none(),
         nb::arg("dPz_dz") = nb::none(),
         nb::arg("center") = nb::none(),
         nb::arg("t_offset") = nb::none(),
         nb::arg("dt_bunch") = nb::none(),
         nb::arg("sig_z") = nb::none(),
         nb::arg("sig_pz") = nb::none(),
         nb::arg("bunch_charge") = nb::none(),
         nb::arg("n_bunch") = nb::none(),
         nb::arg("ix_turn") = nb::none(),
         nb::arg("species") = nb::none(),
         nb::arg("full_6D_coupling_calc") = nb::none(),
         nb::arg("use_particle_start") = nb::none(),
         nb::arg("use_t_coords") = nb::none(),
         nb::arg("use_z_as_t") = nb::none(),
         nb::arg("file_name") = nb::none()
  )
      .def_prop_rw(
          "position_file",
          &BeamInitStruct::position_file,
          &BeamInitStruct::set_position_file,
          "File with particle positions."
      )
      .def_prop_ro(
          "distribution_type",
          &BeamInitStruct::distribution_type,
          nb::keep_alive<0, 1>(),
          "distribution type (in x-px, y-py, and z-pz planes) 'ELLIPSE', 'KV', 'GRID', 'FILE', "
          "'RAN_GAUSS' or '' = 'RAN_GAUSS'"
      )
      .def_prop_rw(
          "spin",
          &BeamInitStruct::spin,
          &BeamInitStruct::set_spin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spin (x, y, z)"
      )
      .def_prop_ro(
          "ellipse",
          &BeamInitStruct::ellipse,
          nb::keep_alive<0, 1>(),
          "Ellipse beam distribution"
      )
      .def_prop_rw(
          "KV",
          &BeamInitStruct::KV,
          &BeamInitStruct::set_KV,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "KV beam distribution"
      )
      .def_prop_ro("grid", &BeamInitStruct::grid, nb::keep_alive<0, 1>(), "Grid beam distribution")
      .def_prop_rw(
          "center_jitter",
          &BeamInitStruct::center_jitter,
          &BeamInitStruct::set_center_jitter,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Bunch center rms jitter"
      )
      .def_prop_rw(
          "emit_jitter",
          &BeamInitStruct::emit_jitter,
          &BeamInitStruct::set_emit_jitter,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "a and b bunch emittance rms jitter normalized to emittance"
      )
      .def_prop_rw(
          "sig_z_jitter",
          &BeamInitStruct::sig_z_jitter,
          &BeamInitStruct::set_sig_z_jitter,
          "bunch length RMS jitter"
      )
      .def_prop_rw(
          "sig_pz_jitter",
          &BeamInitStruct::sig_pz_jitter,
          &BeamInitStruct::set_sig_pz_jitter,
          "RMS pz spread jitter"
      )
      .def_prop_rw(
          "n_particle",
          &BeamInitStruct::n_particle,
          &BeamInitStruct::set_n_particle,
          "Number of particles per bunch."
      )
      .def_prop_rw(
          "renorm_center",
          &BeamInitStruct::renorm_center,
          &BeamInitStruct::set_renorm_center,
          "Renormalize centroid?"
      )
      .def_prop_rw(
          "renorm_sigma",
          &BeamInitStruct::renorm_sigma,
          &BeamInitStruct::set_renorm_sigma,
          "Renormalize sigma?"
      )
      .def_prop_rw(
          "random_engine",
          &BeamInitStruct::random_engine,
          &BeamInitStruct::set_random_engine,
          "Or 'quasi'. Random number engine to use."
      )
      .def_prop_rw(
          "random_gauss_converter",
          &BeamInitStruct::random_gauss_converter,
          &BeamInitStruct::set_random_gauss_converter,
          "Or 'quick' or 'exact'. Uniform to gauss conversion method."
      )
      .def_prop_rw(
          "random_sigma_cutoff",
          &BeamInitStruct::random_sigma_cutoff,
          &BeamInitStruct::set_random_sigma_cutoff,
          "Cut-off in sigmas."
      )
      .def_prop_rw(
          "a_norm_emit",
          &BeamInitStruct::a_norm_emit,
          &BeamInitStruct::set_a_norm_emit,
          "a-mode normalized emittance (emit * beta * gamma)"
      )
      .def_prop_rw(
          "b_norm_emit",
          &BeamInitStruct::b_norm_emit,
          &BeamInitStruct::set_b_norm_emit,
          "b-mode normalized emittance (emit * beta * gamma)"
      )
      .def_prop_rw(
          "a_emit",
          &BeamInitStruct::a_emit,
          &BeamInitStruct::set_a_emit,
          "a-mode emittance"
      )
      .def_prop_rw(
          "b_emit",
          &BeamInitStruct::b_emit,
          &BeamInitStruct::set_b_emit,
          "b-mode emittance"
      )
      .def_prop_rw(
          "dPz_dz",
          &BeamInitStruct::dPz_dz,
          &BeamInitStruct::set_dPz_dz,
          "Correlation of Pz with long position."
      )
      .def_prop_rw(
          "center",
          &BeamInitStruct::center,
          &BeamInitStruct::set_center,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Bench phase space center offset relative to reference."
      )
      .def_prop_rw(
          "t_offset",
          &BeamInitStruct::t_offset,
          &BeamInitStruct::set_t_offset,
          "Time center offset"
      )
      .def_prop_rw(
          "dt_bunch",
          &BeamInitStruct::dt_bunch,
          &BeamInitStruct::set_dt_bunch,
          "Time between bunches."
      )
      .def_prop_rw("sig_z", &BeamInitStruct::sig_z, &BeamInitStruct::set_sig_z, "Z sigma in m.")
      .def_prop_rw("sig_pz", &BeamInitStruct::sig_pz, &BeamInitStruct::set_sig_pz, "pz sigma")
      .def_prop_rw(
          "bunch_charge",
          &BeamInitStruct::bunch_charge,
          &BeamInitStruct::set_bunch_charge,
          "charge (Coul) in a bunch."
      )
      .def_prop_rw(
          "n_bunch",
          &BeamInitStruct::n_bunch,
          &BeamInitStruct::set_n_bunch,
          "Number of bunches."
      )
      .def_prop_rw(
          "ix_turn",
          &BeamInitStruct::ix_turn,
          &BeamInitStruct::set_ix_turn,
          "Turn index used to adjust particles time if needed."
      )
      .def_prop_rw(
          "species",
          &BeamInitStruct::species,
          &BeamInitStruct::set_species,
          "'positron', etc. '' => use referece particle."
      )
      .def_prop_rw(
          "full_6D_coupling_calc",
          &BeamInitStruct::full_6D_coupling_calc,
          &BeamInitStruct::set_full_6D_coupling_calc,
          "Use V from 6x6 1-turn mat to match distribution? Else use 4x4 1-turn mat used."
      )
      .def_prop_rw(
          "use_particle_start",
          &BeamInitStruct::use_particle_start,
          &BeamInitStruct::set_use_particle_start,
          "Use lat%particle_start instead of beam_init%center, %t_offset, and %spin?"
      )
      .def_prop_rw(
          "use_t_coords",
          &BeamInitStruct::use_t_coords,
          &BeamInitStruct::set_use_t_coords,
          "If true, the distributions will be taken as in t-coordinates"
      )
      .def_prop_rw(
          "use_z_as_t",
          &BeamInitStruct::use_z_as_t,
          &BeamInitStruct::set_use_z_as_t,
          "Only used if  use_t_coords = .true. If true,  z describes the t distribution If false, "
          "z describes the s distribution"
      )
      .def_prop_rw(
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
          [](const BeamInitStruct &self, nb::dict &memo) { return BeamInitStruct(self); }
      )
      .def(
          "__eq__",
          [](const BeamInitStruct &self, const BeamInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BeamInitStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BeamInitStruct arrays are not used in structs/routines
  // 2D BeamInitStruct arrays are not used in structs/routines
  // 3D BeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// beam_struct
void init_beam_struct(nb::module_ &m, nb::class_<BeamStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("bunch", &BeamStruct::bunch, nb::keep_alive<0, 1>())

      .def("__repr__", [](const BeamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BeamStruct &self) {
            return BeamStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const BeamStruct &self, nb::dict &memo) { return BeamStruct(self); })
      .def(
          "__eq__",
          [](const BeamStruct &self, const BeamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BeamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BeamStruct arrays are not used in structs/routines
  // 2D BeamStruct arrays are not used in structs/routines
  // 3D BeamStruct arrays are not used in structs/routines
}

// =============================================================================
// bmad_common_struct
void init_bmad_common_struct(nb::module_ &m, nb::class_<BmadCommonStruct> &cls) {
  cls.def(
         nb::init<
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
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("max_aperture_limit") = nb::none(),
         nb::arg("d_orb") = nb::none(),
         nb::arg("default_ds_step") = nb::none(),
         nb::arg("significant_length") = nb::none(),
         nb::arg("rel_tol_tracking") = nb::none(),
         nb::arg("abs_tol_tracking") = nb::none(),
         nb::arg("rel_tol_adaptive_tracking") = nb::none(),
         nb::arg("abs_tol_adaptive_tracking") = nb::none(),
         nb::arg("init_ds_adaptive_tracking") = nb::none(),
         nb::arg("min_ds_adaptive_tracking") = nb::none(),
         nb::arg("fatal_ds_adaptive_tracking") = nb::none(),
         nb::arg("autoscale_amp_abs_tol") = nb::none(),
         nb::arg("autoscale_amp_rel_tol") = nb::none(),
         nb::arg("autoscale_phase_tol") = nb::none(),
         nb::arg("electric_dipole_moment") = nb::none(),
         nb::arg("synch_rad_scale") = nb::none(),
         nb::arg("sad_eps_scale") = nb::none(),
         nb::arg("sad_amp_max") = nb::none(),
         nb::arg("sad_n_div_max") = nb::none(),
         nb::arg("taylor_order") = nb::none(),
         nb::arg("runge_kutta_order") = nb::none(),
         nb::arg("default_integ_order") = nb::none(),
         nb::arg("max_num_runge_kutta_step") = nb::none(),
         nb::arg("rf_phase_below_transition_ref") = nb::none(),
         nb::arg("sr_wakes_on") = nb::none(),
         nb::arg("lr_wakes_on") = nb::none(),
         nb::arg("auto_bookkeeper") = nb::none(),
         nb::arg("high_energy_space_charge_on") = nb::none(),
         nb::arg("high_energy_space_charge_linear") = nb::none(),
         nb::arg("csr_and_space_charge_on") = nb::none(),
         nb::arg("spin_tracking_on") = nb::none(),
         nb::arg("spin_sokolov_ternov_flipping_on") = nb::none(),
         nb::arg("radiation_damping_on") = nb::none(),
         nb::arg("radiation_zero_average") = nb::none(),
         nb::arg("radiation_fluctuations_on") = nb::none(),
         nb::arg("conserve_taylor_maps") = nb::none(),
         nb::arg("absolute_time_tracking") = nb::none(),
         nb::arg("absolute_time_ref_shift") = nb::none(),
         nb::arg("convert_to_kinetic_momentum") = nb::none(),
         nb::arg("normalize_twiss") = nb::none(),
         nb::arg("aperture_limit_on") = nb::none(),
         nb::arg("spin_n0_direction_user_set") = nb::none(),
         nb::arg("debug") = nb::none()
  )
      .def_prop_rw(
          "max_aperture_limit",
          &BmadCommonStruct::max_aperture_limit,
          &BmadCommonStruct::set_max_aperture_limit,
          "Max Aperture."
      )
      .def_prop_rw(
          "d_orb",
          &BmadCommonStruct::d_orb,
          &BmadCommonStruct::set_d_orb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Orbit deltas for the mat6 via tracking calc."
      )
      .def_prop_rw(
          "default_ds_step",
          &BmadCommonStruct::default_ds_step,
          &BmadCommonStruct::set_default_ds_step,
          "Default integration step for eles without an explicit step calc."
      )
      .def_prop_rw(
          "significant_length",
          &BmadCommonStruct::significant_length,
          &BmadCommonStruct::set_significant_length,
          "meter"
      )
      .def_prop_rw(
          "rel_tol_tracking",
          &BmadCommonStruct::rel_tol_tracking,
          &BmadCommonStruct::set_rel_tol_tracking,
          "Closed orbit relative tolerance."
      )
      .def_prop_rw(
          "abs_tol_tracking",
          &BmadCommonStruct::abs_tol_tracking,
          &BmadCommonStruct::set_abs_tol_tracking,
          "Closed orbit absolute tolerance."
      )
      .def_prop_rw(
          "rel_tol_adaptive_tracking",
          &BmadCommonStruct::rel_tol_adaptive_tracking,
          &BmadCommonStruct::set_rel_tol_adaptive_tracking,
          "Runge-Kutta tracking relative tolerance."
      )
      .def_prop_rw(
          "abs_tol_adaptive_tracking",
          &BmadCommonStruct::abs_tol_adaptive_tracking,
          &BmadCommonStruct::set_abs_tol_adaptive_tracking,
          "Runge-Kutta tracking absolute tolerance."
      )
      .def_prop_rw(
          "init_ds_adaptive_tracking",
          &BmadCommonStruct::init_ds_adaptive_tracking,
          &BmadCommonStruct::set_init_ds_adaptive_tracking,
          "Initial step size"
      )
      .def_prop_rw(
          "min_ds_adaptive_tracking",
          &BmadCommonStruct::min_ds_adaptive_tracking,
          &BmadCommonStruct::set_min_ds_adaptive_tracking,
          "Min step size to take."
      )
      .def_prop_rw(
          "fatal_ds_adaptive_tracking",
          &BmadCommonStruct::fatal_ds_adaptive_tracking,
          &BmadCommonStruct::set_fatal_ds_adaptive_tracking,
          "If actual step size is below this particle is lost."
      )
      .def_prop_rw(
          "autoscale_amp_abs_tol",
          &BmadCommonStruct::autoscale_amp_abs_tol,
          &BmadCommonStruct::set_autoscale_amp_abs_tol,
          "Autoscale absolute amplitude tolerance (eV)."
      )
      .def_prop_rw(
          "autoscale_amp_rel_tol",
          &BmadCommonStruct::autoscale_amp_rel_tol,
          &BmadCommonStruct::set_autoscale_amp_rel_tol,
          "Autoscale relative amplitude tolerance"
      )
      .def_prop_rw(
          "autoscale_phase_tol",
          &BmadCommonStruct::autoscale_phase_tol,
          &BmadCommonStruct::set_autoscale_phase_tol,
          "Autoscale phase tolerance."
      )
      .def_prop_rw(
          "electric_dipole_moment",
          &BmadCommonStruct::electric_dipole_moment,
          &BmadCommonStruct::set_electric_dipole_moment,
          "Particle's EDM. Call set_ptc to transfer value to PTC."
      )
      .def_prop_rw(
          "synch_rad_scale",
          &BmadCommonStruct::synch_rad_scale,
          &BmadCommonStruct::set_synch_rad_scale,
          "Synch radiation kick scale. 1 => normal, 0 => no kicks."
      )
      .def_prop_rw(
          "sad_eps_scale",
          &BmadCommonStruct::sad_eps_scale,
          &BmadCommonStruct::set_sad_eps_scale,
          "Used in sad_mult step length calc."
      )
      .def_prop_rw(
          "sad_amp_max",
          &BmadCommonStruct::sad_amp_max,
          &BmadCommonStruct::set_sad_amp_max,
          "Used in sad_mult step length calc."
      )
      .def_prop_rw(
          "sad_n_div_max",
          &BmadCommonStruct::sad_n_div_max,
          &BmadCommonStruct::set_sad_n_div_max,
          "Used in sad_mult step length calc."
      )
      .def_prop_rw(
          "taylor_order",
          &BmadCommonStruct::taylor_order,
          &BmadCommonStruct::set_taylor_order,
          "Taylor order to use. 0 -> default = ptc_private%taylor_order_saved."
      )
      .def_prop_rw(
          "runge_kutta_order",
          &BmadCommonStruct::runge_kutta_order,
          &BmadCommonStruct::set_runge_kutta_order,
          "Runge Kutta order."
      )
      .def_prop_rw(
          "default_integ_order",
          &BmadCommonStruct::default_integ_order,
          &BmadCommonStruct::set_default_integ_order,
          "PTC integration order."
      )
      .def_prop_rw(
          "max_num_runge_kutta_step",
          &BmadCommonStruct::max_num_runge_kutta_step,
          &BmadCommonStruct::set_max_num_runge_kutta_step,
          "Maximum number of RK steps before particle is considered lost."
      )
      .def_prop_rw(
          "rf_phase_below_transition_ref",
          &BmadCommonStruct::rf_phase_below_transition_ref,
          &BmadCommonStruct::set_rf_phase_below_transition_ref,
          "Autoscale uses below transition stable point for RFCavities?"
      )
      .def_prop_rw(
          "sr_wakes_on",
          &BmadCommonStruct::sr_wakes_on,
          &BmadCommonStruct::set_sr_wakes_on,
          "Short range wakefields?"
      )
      .def_prop_rw(
          "lr_wakes_on",
          &BmadCommonStruct::lr_wakes_on,
          &BmadCommonStruct::set_lr_wakes_on,
          "Long range wakefields"
      )
      .def_prop_rw(
          "auto_bookkeeper",
          &BmadCommonStruct::auto_bookkeeper,
          &BmadCommonStruct::set_auto_bookkeeper,
          "Deprecated and no longer used."
      )
      .def_prop_rw(
          "high_energy_space_charge_on",
          &BmadCommonStruct::high_energy_space_charge_on,
          &BmadCommonStruct::set_high_energy_space_charge_on,
          "High energy space charge effect switch."
      )
      .def_prop_rw(
          "high_energy_space_charge_linear",
          &BmadCommonStruct::high_energy_space_charge_linear,
          &BmadCommonStruct::set_high_energy_space_charge_linear,
          "High energy space charge effect switch."
      )
      .def_prop_rw(
          "csr_and_space_charge_on",
          &BmadCommonStruct::csr_and_space_charge_on,
          &BmadCommonStruct::set_csr_and_space_charge_on,
          "Space charge switch."
      )
      .def_prop_rw(
          "spin_tracking_on",
          &BmadCommonStruct::spin_tracking_on,
          &BmadCommonStruct::set_spin_tracking_on,
          "spin tracking?"
      )
      .def_prop_rw(
          "spin_sokolov_ternov_flipping_on",
          &BmadCommonStruct::spin_sokolov_ternov_flipping_on,
          &BmadCommonStruct::set_spin_sokolov_ternov_flipping_on,
          "Spin flipping during synchrotron radiation emission?"
      )
      .def_prop_rw(
          "radiation_damping_on",
          &BmadCommonStruct::radiation_damping_on,
          &BmadCommonStruct::set_radiation_damping_on,
          "Radiation damping toggle."
      )
      .def_prop_rw(
          "radiation_zero_average",
          &BmadCommonStruct::radiation_zero_average,
          &BmadCommonStruct::set_radiation_zero_average,
          "Shift damping to be zero on the zero orbit to get rid of sawtooth?"
      )
      .def_prop_rw(
          "radiation_fluctuations_on",
          &BmadCommonStruct::radiation_fluctuations_on,
          &BmadCommonStruct::set_radiation_fluctuations_on,
          "Radiation fluctuations toggle."
      )
      .def_prop_rw(
          "conserve_taylor_maps",
          &BmadCommonStruct::conserve_taylor_maps,
          &BmadCommonStruct::set_conserve_taylor_maps,
          "Enable bookkeeper to set ele%taylor_map_includes_offsets = F?"
      )
      .def_prop_rw(
          "absolute_time_tracking",
          &BmadCommonStruct::absolute_time_tracking,
          &BmadCommonStruct::set_absolute_time_tracking,
          "Absolute or relative time tracking?"
      )
      .def_prop_rw(
          "absolute_time_ref_shift",
          &BmadCommonStruct::absolute_time_ref_shift,
          &BmadCommonStruct::set_absolute_time_ref_shift,
          "Apply reference time shift when using absolute time tracking?"
      )
      .def_prop_rw(
          "convert_to_kinetic_momentum",
          &BmadCommonStruct::convert_to_kinetic_momentum,
          &BmadCommonStruct::set_convert_to_kinetic_momentum,
          "Cancel kicks due to finite vector potential when doing symplectic tracking? Set to True "
          "to test symp_lie_bmad against runge_kutta."
      )
      .def_prop_rw(
          "normalize_twiss",
          &BmadCommonStruct::normalize_twiss,
          &BmadCommonStruct::set_normalize_twiss,
          "Normalize matrix when computing Twiss for off-energy ref?"
      )
      .def_prop_rw(
          "aperture_limit_on",
          &BmadCommonStruct::aperture_limit_on,
          &BmadCommonStruct::set_aperture_limit_on,
          "Use apertures in tracking?"
      )
      .def_prop_rw(
          "spin_n0_direction_user_set",
          &BmadCommonStruct::spin_n0_direction_user_set,
          &BmadCommonStruct::set_spin_n0_direction_user_set,
          "User sets direction of n0 for closed geometry branches?"
      )
      .def_prop_rw(
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
          [](const BmadCommonStruct &self, nb::dict &memo) { return BmadCommonStruct(self); }
      )
      .def(
          "__eq__",
          [](const BmadCommonStruct &self, const BmadCommonStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BmadCommonStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BmadCommonStruct arrays are not used in structs/routines
  // 2D BmadCommonStruct arrays are not used in structs/routines
  // 3D BmadCommonStruct arrays are not used in structs/routines
}

// =============================================================================
// bmad_normal_form_struct
void init_bmad_normal_form_struct(nb::module_ &m, nb::class_<BmadNormalFormStruct> &cls) {
  cls.def(
         "__init__",
         [](BmadNormalFormStruct *self, const EleStruct *ele_origin) {
           new (self) BmadNormalFormStruct(ptr_to_opt_ref(ele_origin));
         },
         nb::arg("ele_origin") = nb::none()
  )
      .def_prop_rw(
          "ele_origin",
          &BmadNormalFormStruct::ele_origin,
          &BmadNormalFormStruct::set_ele_origin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Element at which the on-turn map was created."
      )
      .def_prop_ro(
          "M",
          &BmadNormalFormStruct::M,
          nb::keep_alive<0, 1>(),
          "One-turn taylor map: M = A o N o A_inv, N = exp(:h:)"
      )
      .def_prop_ro(
          "A",
          &BmadNormalFormStruct::A,
          nb::keep_alive<0, 1>(),
          "Map from Floquet -> Lab coordinates"
      )
      .def_prop_ro(
          "A_inv",
          &BmadNormalFormStruct::A_inv,
          nb::keep_alive<0, 1>(),
          "Map from Lab -> Floquet coordinates"
      )
      .def_prop_ro(
          "dhdj",
          &BmadNormalFormStruct::dhdj,
          nb::keep_alive<0, 1>(),
          "Nonlinear tune function operating on Floquet coordinates"
      )
      .def_prop_ro(
          "F",
          &BmadNormalFormStruct::F,
          nb::keep_alive<0, 1>(),
          "Vector field factorization in phasor basis:"
      )
      .def_prop_ro("L", &BmadNormalFormStruct::L, nb::keep_alive<0, 1>(), "L component")
      .def_prop_ro("h", &BmadNormalFormStruct::h, nb::keep_alive<0, 1>())

      .def("__repr__", [](const BmadNormalFormStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const BmadNormalFormStruct &self) {
            return BmadNormalFormStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const BmadNormalFormStruct &self, nb::dict &memo) {
            return BmadNormalFormStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const BmadNormalFormStruct &self, const BmadNormalFormStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BmadNormalFormStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BmadNormalFormStruct arrays are not used in structs/routines
  // 2D BmadNormalFormStruct arrays are not used in structs/routines
  // 3D BmadNormalFormStruct arrays are not used in structs/routines
}

// =============================================================================
// bookkeeping_state_struct
void init_bookkeeping_state_struct(nb::module_ &m, nb::class_<BookkeepingStateStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("attributes") = nb::none(),
         nb::arg("control") = nb::none(),
         nb::arg("floor_position") = nb::none(),
         nb::arg("s_position") = nb::none(),
         nb::arg("ref_energy") = nb::none(),
         nb::arg("mat6") = nb::none(),
         nb::arg("rad_int") = nb::none(),
         nb::arg("ptc") = nb::none(),
         nb::arg("has_misalign") = nb::none()
  )
      .def_prop_rw(
          "attributes",
          &BookkeepingStateStruct::attributes,
          &BookkeepingStateStruct::set_attributes,
          "Element dependent attributes: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "control",
          &BookkeepingStateStruct::control,
          &BookkeepingStateStruct::set_control,
          "Lord/slave bookkeeping status: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "floor_position",
          &BookkeepingStateStruct::floor_position,
          &BookkeepingStateStruct::set_floor_position,
          "Global (floor) geometry: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "s_position",
          &BookkeepingStateStruct::s_position,
          &BookkeepingStateStruct::set_s_position,
          "Longitudinal position & element length: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "ref_energy",
          &BookkeepingStateStruct::ref_energy,
          &BookkeepingStateStruct::set_ref_energy,
          "Reference energy and ref time: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "mat6",
          &BookkeepingStateStruct::mat6,
          &BookkeepingStateStruct::set_mat6,
          "Linear transfer map status: super_ok$, ok$ or stale$"
      )
      .def_prop_rw(
          "rad_int",
          &BookkeepingStateStruct::rad_int,
          &BookkeepingStateStruct::set_rad_int,
          "Radiation integrals cache status"
      )
      .def_prop_rw(
          "ptc",
          &BookkeepingStateStruct::ptc,
          &BookkeepingStateStruct::set_ptc,
          "Associated PTC fibre (or layout) status."
      )
      .def_prop_rw(
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
          [](const BookkeepingStateStruct &self, nb::dict &memo) {
            return BookkeepingStateStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const BookkeepingStateStruct &self, const BookkeepingStateStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BookkeepingStateStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BookkeepingStateStruct arrays are not used in structs/routines
  // 2D BookkeepingStateStruct arrays are not used in structs/routines
  // 3D BookkeepingStateStruct arrays are not used in structs/routines
}

// =============================================================================
// bpm_phase_coupling_struct
void init_bpm_phase_coupling_struct(nb::module_ &m, nb::class_<BpmPhaseCouplingStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("K_22a") = nb::none(),
         nb::arg("K_12a") = nb::none(),
         nb::arg("K_11b") = nb::none(),
         nb::arg("K_12b") = nb::none(),
         nb::arg("Cbar22_a") = nb::none(),
         nb::arg("Cbar12_a") = nb::none(),
         nb::arg("Cbar11_b") = nb::none(),
         nb::arg("Cbar12_b") = nb::none(),
         nb::arg("phi_a") = nb::none(),
         nb::arg("phi_b") = nb::none()
  )
      .def_prop_rw(
          "K_22a",
          &BpmPhaseCouplingStruct::K_22a,
          &BpmPhaseCouplingStruct::set_K_22a,
          "In-phase y/x for a-mode oscillations."
      )
      .def_prop_rw(
          "K_12a",
          &BpmPhaseCouplingStruct::K_12a,
          &BpmPhaseCouplingStruct::set_K_12a,
          "Out-of-phase y/x for a-mode oscillations."
      )
      .def_prop_rw(
          "K_11b",
          &BpmPhaseCouplingStruct::K_11b,
          &BpmPhaseCouplingStruct::set_K_11b,
          "In-phase x/y for b-mode oscillations."
      )
      .def_prop_rw(
          "K_12b",
          &BpmPhaseCouplingStruct::K_12b,
          &BpmPhaseCouplingStruct::set_K_12b,
          "Out-of-phase x/y for b-mode oscillations."
      )
      .def_prop_rw(
          "Cbar22_a",
          &BpmPhaseCouplingStruct::Cbar22_a,
          &BpmPhaseCouplingStruct::set_Cbar22_a,
          "Cbar22 as calculated from K_22a."
      )
      .def_prop_rw(
          "Cbar12_a",
          &BpmPhaseCouplingStruct::Cbar12_a,
          &BpmPhaseCouplingStruct::set_Cbar12_a,
          "Cbar12 as calculated from K_12a."
      )
      .def_prop_rw(
          "Cbar11_b",
          &BpmPhaseCouplingStruct::Cbar11_b,
          &BpmPhaseCouplingStruct::set_Cbar11_b,
          "Cbar11 as calculated from K_11b."
      )
      .def_prop_rw(
          "Cbar12_b",
          &BpmPhaseCouplingStruct::Cbar12_b,
          &BpmPhaseCouplingStruct::set_Cbar12_b,
          "Cbar12 as calculated from K_12b."
      )
      .def_prop_rw(
          "phi_a",
          &BpmPhaseCouplingStruct::phi_a,
          &BpmPhaseCouplingStruct::set_phi_a,
          "a-mode betatron phase."
      )
      .def_prop_rw(
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
          [](const BpmPhaseCouplingStruct &self, nb::dict &memo) {
            return BpmPhaseCouplingStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const BpmPhaseCouplingStruct &self, const BpmPhaseCouplingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BpmPhaseCouplingStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BpmPhaseCouplingStruct arrays are not used in structs/routines
  // 2D BpmPhaseCouplingStruct arrays are not used in structs/routines
  // 3D BpmPhaseCouplingStruct arrays are not used in structs/routines
}

// =============================================================================
// branch_struct
void init_branch_struct(nb::module_ &m, nb::class_<BranchStruct> &cls) {
  cls.def(
         "__init__",
         [](BranchStruct *self,
            std::optional<std::string> name,
            std::optional<int> ix_branch,
            std::optional<int> ix_from_branch,
            std::optional<int> ix_from_ele,
            std::optional<int> ix_to_ele,
            std::optional<int> ix_fixer,
            std::optional<int> n_ele_track,
            std::optional<int> n_ele_max,
            const LatStruct *lat,
            const ModeInfoStruct *a,
            const ModeInfoStruct *b,
            const ModeInfoStruct *z,
            const LatParamStruct *param,
            const CoordStruct *particle_start) {
           new (self) BranchStruct(
               name,
               ix_branch,
               ix_from_branch,
               ix_from_ele,
               ix_to_ele,
               ix_fixer,
               n_ele_track,
               n_ele_max,
               ptr_to_opt_ref(lat),
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(z),
               ptr_to_opt_ref(param),
               ptr_to_opt_ref(particle_start)
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_from_branch") = nb::none(),
         nb::arg("ix_from_ele") = nb::none(),
         nb::arg("ix_to_ele") = nb::none(),
         nb::arg("ix_fixer") = nb::none(),
         nb::arg("n_ele_track") = nb::none(),
         nb::arg("n_ele_max") = nb::none(),
         nb::arg("lat") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("z") = nb::none(),
         nb::arg("param") = nb::none(),
         nb::arg("particle_start") = nb::none()
  )
      .def_prop_rw(
          "name",
          &BranchStruct::name,
          &BranchStruct::set_name,
          "Name of line that defines the branch."
      )
      .def_prop_rw(
          "ix_branch",
          &BranchStruct::ix_branch,
          &BranchStruct::set_ix_branch,
          "Index of this branch. 0 => Main branch"
      )
      .def_prop_rw(
          "ix_from_branch",
          &BranchStruct::ix_from_branch,
          &BranchStruct::set_ix_from_branch,
          "-1 => No creating fork element to this branch."
      )
      .def_prop_rw(
          "ix_from_ele",
          &BranchStruct::ix_from_ele,
          &BranchStruct::set_ix_from_ele,
          "Index of creating fork element which forks to this branch."
      )
      .def_prop_rw(
          "ix_to_ele",
          &BranchStruct::ix_to_ele,
          &BranchStruct::set_ix_to_ele,
          "Index of element in this branch that creating fork element forks to."
      )
      .def_prop_rw(
          "ix_fixer",
          &BranchStruct::ix_fixer,
          &BranchStruct::set_ix_fixer,
          "Index of active fixer or beginning_ele element."
      )
      .def_prop_rw("n_ele_track", &BranchStruct::n_ele_track, &BranchStruct::set_n_ele_track)
      .def_prop_rw("n_ele_max", &BranchStruct::n_ele_max, &BranchStruct::set_n_ele_max)
      .def_prop_rw(
          "lat",
          &BranchStruct::lat,
          &BranchStruct::set_lat,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "a",
          &BranchStruct::a,
          &BranchStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Note: Tunes are the fractional part."
      )
      .def_prop_rw(
          "b",
          &BranchStruct::b,
          &BranchStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Note: Tunes are the fractional part."
      )
      .def_prop_rw(
          "z",
          &BranchStruct::z,
          &BranchStruct::set_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Note: Tunes are the fractional part."
      )
      .def_prop_ro("ele", &BranchStruct::ele, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "param",
          &BranchStruct::param,
          &BranchStruct::set_param,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "particle_start",
          &BranchStruct::particle_start,
          &BranchStruct::set_particle_start,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("wall3d", &BranchStruct::wall3d, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return BranchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const BranchStruct &self, nb::dict &memo) { return BranchStruct(self); }
      )
      .def(
          "__eq__",
          [](const BranchStruct &self, const BranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BranchStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_bunch_params_struct(nb::module_ &m, nb::class_<BunchParamsStruct> &cls) {
  cls.def(
         "__init__",
         [](BunchParamsStruct *self,
            const CoordStruct *centroid,
            const TwissStruct *x,
            const TwissStruct *y,
            const TwissStruct *z,
            const TwissStruct *a,
            const TwissStruct *b,
            const TwissStruct *c,
            std::optional<std::vector<std::vector<double>>> sigma,
            std::optional<std::vector<double>> rel_max,
            std::optional<std::vector<double>> rel_min,
            std::optional<double> s,
            std::optional<double> t,
            std::optional<double> sigma_t,
            std::optional<double> charge_live,
            std::optional<double> charge_tot,
            std::optional<int> n_particle_tot,
            std::optional<int> n_particle_live,
            std::optional<int> n_particle_lost_in_ele,
            std::optional<int> n_good_steps,
            std::optional<int> n_bad_steps,
            std::optional<int> ix_ele,
            std::optional<int> location,
            std::optional<bool> twiss_valid) {
           new (self) BunchParamsStruct(
               ptr_to_opt_ref(centroid),
               ptr_to_opt_ref(x),
               ptr_to_opt_ref(y),
               ptr_to_opt_ref(z),
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(c),
               sigma,
               rel_max,
               rel_min,
               s,
               t,
               sigma_t,
               charge_live,
               charge_tot,
               n_particle_tot,
               n_particle_live,
               n_particle_lost_in_ele,
               n_good_steps,
               n_bad_steps,
               ix_ele,
               location,
               twiss_valid
           );
         },
         nb::arg("centroid") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("z") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("c") = nb::none(),
         nb::arg("sigma") = nb::none(),
         nb::arg("rel_max") = nb::none(),
         nb::arg("rel_min") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("t") = nb::none(),
         nb::arg("sigma_t") = nb::none(),
         nb::arg("charge_live") = nb::none(),
         nb::arg("charge_tot") = nb::none(),
         nb::arg("n_particle_tot") = nb::none(),
         nb::arg("n_particle_live") = nb::none(),
         nb::arg("n_particle_lost_in_ele") = nb::none(),
         nb::arg("n_good_steps") = nb::none(),
         nb::arg("n_bad_steps") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("location") = nb::none(),
         nb::arg("twiss_valid") = nb::none()
  )
      .def_prop_rw(
          "centroid",
          &BunchParamsStruct::centroid,
          &BunchParamsStruct::set_centroid,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Lab frame"
      )
      .def_prop_rw(
          "x",
          &BunchParamsStruct::x,
          &BunchParamsStruct::set_x,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Projected Twiss parameters"
      )
      .def_prop_rw(
          "y",
          &BunchParamsStruct::y,
          &BunchParamsStruct::set_y,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Projected Twiss parameters"
      )
      .def_prop_rw(
          "z",
          &BunchParamsStruct::z,
          &BunchParamsStruct::set_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Projected Twiss parameters"
      )
      .def_prop_rw(
          "a",
          &BunchParamsStruct::a,
          &BunchParamsStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normal mode twiss parameters"
      )
      .def_prop_rw(
          "b",
          &BunchParamsStruct::b,
          &BunchParamsStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normal mode twiss parameters"
      )
      .def_prop_rw(
          "c",
          &BunchParamsStruct::c,
          &BunchParamsStruct::set_c,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normal mode twiss parameters"
      )
      .def_prop_rw(
          "sigma",
          &BunchParamsStruct::sigma,
          &BunchParamsStruct::set_sigma,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "beam size matrix"
      )
      .def_prop_rw(
          "rel_max",
          &BunchParamsStruct::rel_max,
          &BunchParamsStruct::set_rel_max,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Max orbit relative to centroid. 7 -> time."
      )
      .def_prop_rw(
          "rel_min",
          &BunchParamsStruct::rel_min,
          &BunchParamsStruct::set_rel_min,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Min orbit relative to_centroid. 7 -> time."
      )
      .def_prop_rw("s", &BunchParamsStruct::s, &BunchParamsStruct::set_s, "Longitudinal position.")
      .def_prop_rw("t", &BunchParamsStruct::t, &BunchParamsStruct::set_t, "Time.")
      .def_prop_rw(
          "sigma_t",
          &BunchParamsStruct::sigma_t,
          &BunchParamsStruct::set_sigma_t,
          "RMS of time spread."
      )
      .def_prop_rw(
          "charge_live",
          &BunchParamsStruct::charge_live,
          &BunchParamsStruct::set_charge_live,
          "Charge of all non-lost particle"
      )
      .def_prop_rw(
          "charge_tot",
          &BunchParamsStruct::charge_tot,
          &BunchParamsStruct::set_charge_tot,
          "Charge of all particles."
      )
      .def_prop_rw(
          "n_particle_tot",
          &BunchParamsStruct::n_particle_tot,
          &BunchParamsStruct::set_n_particle_tot,
          "Total number of particles"
      )
      .def_prop_rw(
          "n_particle_live",
          &BunchParamsStruct::n_particle_live,
          &BunchParamsStruct::set_n_particle_live,
          "Number of non-lost particles"
      )
      .def_prop_rw(
          "n_particle_lost_in_ele",
          &BunchParamsStruct::n_particle_lost_in_ele,
          &BunchParamsStruct::set_n_particle_lost_in_ele,
          "Number lost in element (not calculated by Bmad)"
      )
      .def_prop_rw(
          "n_good_steps",
          &BunchParamsStruct::n_good_steps,
          &BunchParamsStruct::set_n_good_steps,
          "Number of good steps (set when tracking with space charge)"
      )
      .def_prop_rw(
          "n_bad_steps",
          &BunchParamsStruct::n_bad_steps,
          &BunchParamsStruct::set_n_bad_steps,
          "Number of bad steps (set when tracking with space charge)"
      )
      .def_prop_rw(
          "ix_ele",
          &BunchParamsStruct::ix_ele,
          &BunchParamsStruct::set_ix_ele,
          "Lattice element where params evaluated at."
      )
      .def_prop_rw(
          "location",
          &BunchParamsStruct::location,
          &BunchParamsStruct::set_location,
          "Location in element: upstream_end$, inside$, or downstream_end$"
      )
      .def_prop_rw(
          "twiss_valid",
          &BunchParamsStruct::twiss_valid,
          &BunchParamsStruct::set_twiss_valid,
          "Is the data here valid? Note: IF there is no energy variation (RF off) twiss_valid may "
          "be true but in this case the z-twiss will not be valid."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchParamsStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchParamsStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const BunchParamsStruct &self, nb::dict &memo) { return BunchParamsStruct(self); }
      )
      .def(
          "__eq__",
          [](const BunchParamsStruct &self, const BunchParamsStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BunchParamsStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_bunch_struct(nb::module_ &m, nb::class_<BunchStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("ix_z") = nb::none(),
         nb::arg("charge_tot") = nb::none(),
         nb::arg("charge_live") = nb::none(),
         nb::arg("z_center") = nb::none(),
         nb::arg("t_center") = nb::none(),
         nb::arg("t0") = nb::none(),
         nb::arg("drift_between_t_and_s") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_bunch") = nb::none(),
         nb::arg("ix_turn") = nb::none(),
         nb::arg("n_live") = nb::none(),
         nb::arg("n_good") = nb::none(),
         nb::arg("n_bad") = nb::none()
  )
      .def_prop_ro("particle", &BunchStruct::particle, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "ix_z",
          &BunchStruct::ix_z,
          &BunchStruct::set_ix_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "bunch%ix_z(1) is index of head particle, etc."
      )
      .def_prop_rw(
          "charge_tot",
          &BunchStruct::charge_tot,
          &BunchStruct::set_charge_tot,
          "Total charge in a bunch (Coul)."
      )
      .def_prop_rw(
          "charge_live",
          &BunchStruct::charge_live,
          &BunchStruct::set_charge_live,
          "Charge of live particles (Coul)."
      )
      .def_prop_rw(
          "z_center",
          &BunchStruct::z_center,
          &BunchStruct::set_z_center,
          "Longitudinal center of bunch at creation time. Note: Generally, z_center of bunch #1 is "
          "0 and z_center of the other bunches is negative."
      )
      .def_prop_rw(
          "t_center",
          &BunchStruct::t_center,
          &BunchStruct::set_t_center,
          "Center of bunch at creation time relative to head bunch."
      )
      .def_prop_rw(
          "t0",
          &BunchStruct::t0,
          &BunchStruct::set_t0,
          "Used by track1_bunch_space_charge for tracking so particles have constant t."
      )
      .def_prop_rw(
          "drift_between_t_and_s",
          &BunchStruct::drift_between_t_and_s,
          &BunchStruct::set_drift_between_t_and_s,
          "Drift (ignore any fields) instead of tracking to speed up the calculation? This can "
          "only be done under certain circumstances."
      )
      .def_prop_rw(
          "ix_ele",
          &BunchStruct::ix_ele,
          &BunchStruct::set_ix_ele,
          "Nominal element bunch is at. But, EG, dead particles can be someplace else."
      )
      .def_prop_rw(
          "ix_bunch",
          &BunchStruct::ix_bunch,
          &BunchStruct::set_ix_bunch,
          "Bunch index. Head bunch = 1, etc."
      )
      .def_prop_rw(
          "ix_turn",
          &BunchStruct::ix_turn,
          &BunchStruct::set_ix_turn,
          "Turn index for long term tracking. ix_turn = 0 before end of first turn, etc."
      )
      .def_prop_rw("n_live", &BunchStruct::n_live, &BunchStruct::set_n_live)
      .def_prop_rw(
          "n_good",
          &BunchStruct::n_good,
          &BunchStruct::set_n_good,
          "Number of accepted steps when using adaptive step size control."
      )
      .def_prop_rw(
          "n_bad",
          &BunchStruct::n_bad,
          &BunchStruct::set_n_bad,
          "Number of rejected steps when using adaptive step size control."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const BunchStruct &self, nb::dict &memo) { return BunchStruct(self); }
      )
      .def(
          "__eq__",
          [](const BunchStruct &self, const BunchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BunchStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_bunch_track_struct(nb::module_ &m, nb::class_<BunchTrackStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<int>>(),
         nb::arg("ds_save") = nb::none(),
         nb::arg("n_pt") = nb::none()
  )
      .def_prop_ro("pt", &BunchTrackStruct::pt, nb::keep_alive<0, 1>(), "Array indexed from 0")
      .def_prop_rw(
          "ds_save",
          &BunchTrackStruct::ds_save,
          &BunchTrackStruct::set_ds_save,
          "Min distance between points."
      )
      .def_prop_rw(
          "n_pt",
          &BunchTrackStruct::n_pt,
          &BunchTrackStruct::set_n_pt,
          "Track upper bound"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return BunchTrackStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = BunchTrackStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const BunchTrackStruct &self, nb::dict &memo) { return BunchTrackStruct(self); }
      )
      .def(
          "__eq__",
          [](const BunchTrackStruct &self, const BunchTrackStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BunchTrackStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_bicubic_cmplx_coef_struct(nb::module_ &m, nb::class_<BicubicCmplxCoefStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::vector<std::complex<double>>>>,
             std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("i_box") = nb::none()
  )
      .def_prop_rw(
          "coef",
          &BicubicCmplxCoefStruct::coef,
          &BicubicCmplxCoefStruct::set_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Coefs"
      )
      .def_prop_rw(
          "i_box",
          &BicubicCmplxCoefStruct::i_box,
          &BicubicCmplxCoefStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const BicubicCmplxCoefStruct &self, nb::dict &memo) {
            return BicubicCmplxCoefStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const BicubicCmplxCoefStruct &self, const BicubicCmplxCoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const BicubicCmplxCoefStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D BicubicCmplxCoefStruct arrays are not used in structs/routines
  // 2D BicubicCmplxCoefStruct arrays are not used in structs/routines
  bind_FTypeArrayND<BicubicCmplxCoefStructArray3D>(m, "BicubicCmplxCoefStructArray3D");
}