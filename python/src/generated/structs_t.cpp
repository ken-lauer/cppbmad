#include "pybmad/generated/structs_t.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// target_point_struct
void init_target_point_struct(py::module &m, py::class_<TargetPointStruct> &cls) {
  cls.def(py::init<>())
      // TargetPointStruct.r (1D_NOT_real - (x, y, z)
      .def_property("r", &TargetPointStruct::r, &TargetPointStruct::set_r)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TargetPointStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TargetPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TargetPointStruct &self) {
            return TargetPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TargetPointStruct &self, py::dict &memo) { return TargetPointStruct(self); }
      )

      ;

  bind_FTypeArrayND<TargetPointStructArray1D>(m, "TargetPointStructArray1D");
  bind_FTypeAlloc1D<TargetPointStructAlloc1D>(m, "TargetPointStructAlloc1D");
  // 2D TargetPointStruct arrays are not used in structs/routines
  // 3D TargetPointStruct arrays are not used in structs/routines
}

// =============================================================================
// taylor_struct
void init_taylor_struct(py::module &m, py::class_<TaylorStruct> &cls) {
  cls.def(py::init<>())
      // TaylorStruct.ref (0D_NOT_real -
      .def_property("ref", &TaylorStruct::ref, &TaylorStruct::set_ref)
      // TaylorStruct.term (1D_PTR_type -
      .def_property_readonly("term", &TaylorStruct::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaylorStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaylorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaylorStruct &self) {
            return TaylorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaylorStruct &self, py::dict &memo) { return TaylorStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaylorStructArray1D>(m, "TaylorStructArray1D");
  bind_FTypeAlloc1D<TaylorStructAlloc1D>(m, "TaylorStructAlloc1D");
  // 2D TaylorStruct arrays are not used in structs/routines
  // 3D TaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// taylor_term_struct
void init_taylor_term_struct(py::module &m, py::class_<TaylorTermStruct> &cls) {
  cls.def(py::init<>())
      // TaylorTermStruct.coef (0D_NOT_real -
      .def_property("coef", &TaylorTermStruct::coef, &TaylorTermStruct::set_coef)
      // TaylorTermStruct.expn (1D_NOT_integer -
      .def_property("expn", &TaylorTermStruct::expn, &TaylorTermStruct::set_expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaylorTermStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaylorTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaylorTermStruct &self) {
            return TaylorTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaylorTermStruct &self, py::dict &memo) { return TaylorTermStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaylorTermStructArray1D>(m, "TaylorTermStructArray1D");
  bind_FTypeAlloc1D<TaylorTermStructAlloc1D>(m, "TaylorTermStructAlloc1D");
  // 2D TaylorTermStruct arrays are not used in structs/routines
  // 3D TaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// track_point_struct
void init_track_point_struct(py::module &m, py::class_<TrackPointStruct> &cls) {
  cls.def(py::init<>())
      // TrackPointStruct.s_lab (0D_NOT_real - Longitudinal lab coord with respect to the upstream
      // end.
      .def_property("s_lab", &TrackPointStruct::s_lab, &TrackPointStruct::set_s_lab)
      // TrackPointStruct.s_body (0D_NOT_real - Longitudinal body coord with respect to the entrance
      // end.
      .def_property("s_body", &TrackPointStruct::s_body, &TrackPointStruct::set_s_body)
      // TrackPointStruct.orb (0D_NOT_type - Particle position in lab coords.
      .def_property("orb", &TrackPointStruct::orb, &TrackPointStruct::set_orb)
      // TrackPointStruct.field (0D_NOT_type - E&M fields in lab coordinates.
      .def_property("field", &TrackPointStruct::field, &TrackPointStruct::set_field)
      // TrackPointStruct.strong_beam (0D_NOT_type - Strong beam info for beambeam element.
      .def_property(
          "strong_beam",
          &TrackPointStruct::strong_beam,
          &TrackPointStruct::set_strong_beam
      )
      // TrackPointStruct.vec0 (1D_NOT_real - 0th order part of xfer map from the beginning.
      .def_property("vec0", &TrackPointStruct::vec0, &TrackPointStruct::set_vec0)
      // TrackPointStruct.mat6 (2D_NOT_real - 1st order part of xfer map (transfer matrix).
      .def_property("mat6", &TrackPointStruct::mat6, &TrackPointStruct::set_mat6)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TrackPointStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TrackPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TrackPointStruct &self) {
            return TrackPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TrackPointStruct &self, py::dict &memo) { return TrackPointStruct(self); }
      )

      ;

  bind_FTypeArrayND<TrackPointStructArray1D>(m, "TrackPointStructArray1D");
  bind_FTypeAlloc1D<TrackPointStructAlloc1D>(m, "TrackPointStructAlloc1D");
  // 2D TrackPointStruct arrays are not used in structs/routines
  // 3D TrackPointStruct arrays are not used in structs/routines
}

// =============================================================================
// track_struct
void init_track_struct(py::module &m, py::class_<TrackStruct> &cls) {
  cls.def(py::init<>())
      // TrackStruct.pt (1D_ALLOC_type - Array of track points indexed from 0.
      .def_property_readonly("pt", &TrackStruct::pt)
      // TrackStruct.ds_save (0D_NOT_real - Min distance between points. Not positive => Save at all
      // points.
      .def_property("ds_save", &TrackStruct::ds_save, &TrackStruct::set_ds_save)
      // TrackStruct.n_pt (0D_NOT_integer - Track upper bound for %pt(0:) array. n_bad and n_ok are
      // used by adaptive trackers to record the number of times the step length had to be
      // shortened.
      .def_property("n_pt", &TrackStruct::n_pt, &TrackStruct::set_n_pt)
      // TrackStruct.n_bad (0D_NOT_integer - Number of 'bad' steps where the step length was
      // shortened.
      .def_property("n_bad", &TrackStruct::n_bad, &TrackStruct::set_n_bad)
      // TrackStruct.n_ok (0D_NOT_integer - Number of 'good' steps where the step length was not
      // shortened.
      .def_property("n_ok", &TrackStruct::n_ok, &TrackStruct::set_n_ok)

      .def("__repr__", [](const TrackStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TrackStruct &self) {
            return TrackStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TrackStruct &self, py::dict &memo) { return TrackStruct(self); }
      )

      ;

  // 1D TrackStruct arrays are not used in structs/routines
  // 2D TrackStruct arrays are not used in structs/routines
  // 3D TrackStruct arrays are not used in structs/routines
}

// =============================================================================
// twiss_struct
void init_twiss_struct(py::module &m, py::class_<TwissStruct> &cls) {
  cls.def(py::init<>())
      // TwissStruct.beta (0D_NOT_real -
      .def_property("beta", &TwissStruct::beta, &TwissStruct::set_beta)
      // TwissStruct.alpha (0D_NOT_real -
      .def_property("alpha", &TwissStruct::alpha, &TwissStruct::set_alpha)
      // TwissStruct.gamma (0D_NOT_real -
      .def_property("gamma", &TwissStruct::gamma, &TwissStruct::set_gamma)
      // TwissStruct.phi (0D_NOT_real -
      .def_property("phi", &TwissStruct::phi, &TwissStruct::set_phi)
      // TwissStruct.eta (0D_NOT_real -
      .def_property("eta", &TwissStruct::eta, &TwissStruct::set_eta)
      // TwissStruct.etap (0D_NOT_real -
      .def_property("etap", &TwissStruct::etap, &TwissStruct::set_etap)
      // TwissStruct.deta_ds (0D_NOT_real -
      .def_property("deta_ds", &TwissStruct::deta_ds, &TwissStruct::set_deta_ds)
      // TwissStruct.sigma (0D_NOT_real -
      .def_property("sigma", &TwissStruct::sigma, &TwissStruct::set_sigma)
      // TwissStruct.sigma_p (0D_NOT_real -
      .def_property("sigma_p", &TwissStruct::sigma_p, &TwissStruct::set_sigma_p)
      // TwissStruct.emit (0D_NOT_real -
      .def_property("emit", &TwissStruct::emit, &TwissStruct::set_emit)
      // TwissStruct.norm_emit (0D_NOT_real -
      .def_property("norm_emit", &TwissStruct::norm_emit, &TwissStruct::set_norm_emit)
      // TwissStruct.chrom (0D_NOT_real -
      .def_property("chrom", &TwissStruct::chrom, &TwissStruct::set_chrom)
      // TwissStruct.dbeta_dpz (0D_NOT_real -
      .def_property("dbeta_dpz", &TwissStruct::dbeta_dpz, &TwissStruct::set_dbeta_dpz)
      // TwissStruct.dalpha_dpz (0D_NOT_real -
      .def_property("dalpha_dpz", &TwissStruct::dalpha_dpz, &TwissStruct::set_dalpha_dpz)
      // TwissStruct.deta_dpz (0D_NOT_real -
      .def_property("deta_dpz", &TwissStruct::deta_dpz, &TwissStruct::set_deta_dpz)
      // TwissStruct.detap_dpz (0D_NOT_real -
      .def_property("detap_dpz", &TwissStruct::detap_dpz, &TwissStruct::set_detap_dpz)

      .def("__repr__", [](const TwissStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TwissStruct &self) {
            return TwissStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TwissStruct &self, py::dict &memo) { return TwissStruct(self); }
      )

      ;

  // 1D TwissStruct arrays are not used in structs/routines
  // 2D TwissStruct arrays are not used in structs/routines
  // 3D TwissStruct arrays are not used in structs/routines
}

// =============================================================================
// tricubic_cmplx_coef_struct
void init_tricubic_cmplx_coef_struct(py::module &m, py::class_<TricubicCmplxCoefStruct> &cls) {
  cls.def(py::init<>())
      // TricubicCmplxCoefStruct.coef (3D_NOT_complex - Coefs
      .def_property("coef", &TricubicCmplxCoefStruct::coef, &TricubicCmplxCoefStruct::set_coef)
      // TricubicCmplxCoefStruct.i_box (1D_NOT_integer - index at lower box corner.
      .def_property("i_box", &TricubicCmplxCoefStruct::i_box, &TricubicCmplxCoefStruct::set_i_box)

      .def("__repr__", [](const TricubicCmplxCoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TricubicCmplxCoefStruct &self) {
            return TricubicCmplxCoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TricubicCmplxCoefStruct &self, py::dict &memo) {
            return TricubicCmplxCoefStruct(self);
          }
      )

      ;

  // 1D TricubicCmplxCoefStruct arrays are not used in structs/routines
  // 2D TricubicCmplxCoefStruct arrays are not used in structs/routines
  bind_FTypeArrayND<TricubicCmplxCoefStructArray3D>(m, "TricubicCmplxCoefStructArray3D");
}

// =============================================================================
// tao_beam_branch_struct
void init_tao_beam_branch_struct(py::module &m, py::class_<TaoBeamBranchStruct> &cls) {
  cls.def(py::init<>())
      // TaoBeamBranchStruct.beam_at_start (0D_NOT_type - Initial beam
      .def_property(
          "beam_at_start",
          &TaoBeamBranchStruct::beam_at_start,
          &TaoBeamBranchStruct::set_beam_at_start
      )
      // TaoBeamBranchStruct.beam_init (0D_NOT_type - User set beam distrubution at track start.
      .def_property(
          "beam_init",
          &TaoBeamBranchStruct::beam_init,
          &TaoBeamBranchStruct::set_beam_init
      )
      // TaoBeamBranchStruct.beam_init_used (0D_NOT_type - beam distribution with emit values set.
      .def_property(
          "beam_init_used",
          &TaoBeamBranchStruct::beam_init_used,
          &TaoBeamBranchStruct::set_beam_init_used
      )
      // TaoBeamBranchStruct.init_starting_distribution (0D_NOT_logical - Init beam
      .def_property(
          "init_starting_distribution",
          &TaoBeamBranchStruct::init_starting_distribution,
          &TaoBeamBranchStruct::set_init_starting_distribution
      )
      // TaoBeamBranchStruct.track_start (0D_NOT_character - Tracking start element.
      .def_property(
          "track_start",
          &TaoBeamBranchStruct::track_start,
          &TaoBeamBranchStruct::set_track_start
      )
      // TaoBeamBranchStruct.track_end (0D_NOT_character -
      .def_property(
          "track_end",
          &TaoBeamBranchStruct::track_end,
          &TaoBeamBranchStruct::set_track_end
      )
      // TaoBeamBranchStruct.ix_branch (0D_NOT_integer - Branch tracked. If track_start or track_end
      // is a lord, ix_track_start/end index will be a index of slave.
      .def_property(
          "ix_branch",
          &TaoBeamBranchStruct::ix_branch,
          &TaoBeamBranchStruct::set_ix_branch
      )
      // TaoBeamBranchStruct.ix_track_start (0D_NOT_integer - Element track start index.
      .def_property(
          "ix_track_start",
          &TaoBeamBranchStruct::ix_track_start,
          &TaoBeamBranchStruct::set_ix_track_start
      )
      // TaoBeamBranchStruct.ix_track_end (0D_NOT_integer - Element track end index
      .def_property(
          "ix_track_end",
          &TaoBeamBranchStruct::ix_track_end,
          &TaoBeamBranchStruct::set_ix_track_end
      )

      .def("__repr__", [](const TaoBeamBranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBeamBranchStruct &self) {
            return TaoBeamBranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBeamBranchStruct &self, py::dict &memo) { return TaoBeamBranchStruct(self); }
      )

      ;

  // 1D TaoBeamBranchStruct arrays are not used in structs/routines
  // 2D TaoBeamBranchStruct arrays are not used in structs/routines
  // 3D TaoBeamBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_beam_uni_struct
void init_tao_beam_uni_struct(py::module &m, py::class_<TaoBeamUniStruct> &cls) {
  cls.def(py::init<>())
      // TaoBeamUniStruct.saved_at (0D_NOT_character -
      .def_property("saved_at", &TaoBeamUniStruct::saved_at, &TaoBeamUniStruct::set_saved_at)
      // TaoBeamUniStruct.dump_file (0D_NOT_character -
      .def_property("dump_file", &TaoBeamUniStruct::dump_file, &TaoBeamUniStruct::set_dump_file)
      // TaoBeamUniStruct.dump_at (0D_NOT_character -
      .def_property("dump_at", &TaoBeamUniStruct::dump_at, &TaoBeamUniStruct::set_dump_at)
      // TaoBeamUniStruct.track_beam_in_universe (0D_NOT_logical - Beam tracking enabled in this
      // universe?
      .def_property(
          "track_beam_in_universe",
          &TaoBeamUniStruct::track_beam_in_universe,
          &TaoBeamUniStruct::set_track_beam_in_universe
      )
      // TaoBeamUniStruct.always_reinit (0D_NOT_logical -
      .def_property(
          "always_reinit",
          &TaoBeamUniStruct::always_reinit,
          &TaoBeamUniStruct::set_always_reinit
      )

      .def("__repr__", [](const TaoBeamUniStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBeamUniStruct &self) {
            return TaoBeamUniStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBeamUniStruct &self, py::dict &memo) { return TaoBeamUniStruct(self); }
      )

      ;

  // 1D TaoBeamUniStruct arrays are not used in structs/routines
  // 2D TaoBeamUniStruct arrays are not used in structs/routines
  // 3D TaoBeamUniStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_orientation_struct
void init_tao_building_wall_orientation_struct(
    py::module &m,
    py::class_<TaoBuildingWallOrientationStruct> &cls
) {
  cls.def(py::init<>())
      // TaoBuildingWallOrientationStruct.theta (0D_NOT_real -
      .def_property(
          "theta",
          &TaoBuildingWallOrientationStruct::theta,
          &TaoBuildingWallOrientationStruct::set_theta
      )
      // TaoBuildingWallOrientationStruct.x_offset (0D_NOT_real -
      .def_property(
          "x_offset",
          &TaoBuildingWallOrientationStruct::x_offset,
          &TaoBuildingWallOrientationStruct::set_x_offset
      )
      // TaoBuildingWallOrientationStruct.z_offset (0D_NOT_real -
      .def_property(
          "z_offset",
          &TaoBuildingWallOrientationStruct::z_offset,
          &TaoBuildingWallOrientationStruct::set_z_offset
      )

      .def("__repr__", [](const TaoBuildingWallOrientationStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallOrientationStruct &self) {
            return TaoBuildingWallOrientationStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallOrientationStruct &self, py::dict &memo) {
            return TaoBuildingWallOrientationStruct(self);
          }
      )

      ;

  // 1D TaoBuildingWallOrientationStruct arrays are not used in structs/routines
  // 2D TaoBuildingWallOrientationStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallOrientationStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_point_struct
void init_tao_building_wall_point_struct(
    py::module &m,
    py::class_<TaoBuildingWallPointStruct> &cls
) {
  cls.def(py::init<>())
      // TaoBuildingWallPointStruct.z (0D_NOT_real - Global floor position
      .def_property("z", &TaoBuildingWallPointStruct::z, &TaoBuildingWallPointStruct::set_z)
      // TaoBuildingWallPointStruct.x (0D_NOT_real - Global floor position
      .def_property("x", &TaoBuildingWallPointStruct::x, &TaoBuildingWallPointStruct::set_x)
      // TaoBuildingWallPointStruct.radius (0D_NOT_real - Arc radius. +r -> CW rotation, same as
      // bends.
      .def_property(
          "radius",
          &TaoBuildingWallPointStruct::radius,
          &TaoBuildingWallPointStruct::set_radius
      )
      // TaoBuildingWallPointStruct.z_center (0D_NOT_real - Arc center.
      .def_property(
          "z_center",
          &TaoBuildingWallPointStruct::z_center,
          &TaoBuildingWallPointStruct::set_z_center
      )
      // TaoBuildingWallPointStruct.x_center (0D_NOT_real - Arc center.
      .def_property(
          "x_center",
          &TaoBuildingWallPointStruct::x_center,
          &TaoBuildingWallPointStruct::set_x_center
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoBuildingWallPointStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoBuildingWallPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallPointStruct &self) {
            return TaoBuildingWallPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallPointStruct &self, py::dict &memo) {
            return TaoBuildingWallPointStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoBuildingWallPointStructArray1D>(m, "TaoBuildingWallPointStructArray1D");
  bind_FTypeAlloc1D<TaoBuildingWallPointStructAlloc1D>(m, "TaoBuildingWallPointStructAlloc1D");
  // 2D TaoBuildingWallPointStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallPointStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_section_struct
void init_tao_building_wall_section_struct(
    py::module &m,
    py::class_<TaoBuildingWallSectionStruct> &cls
) {
  cls.def(py::init<>())
      // TaoBuildingWallSectionStruct.name (0D_NOT_character -
      .def_property(
          "name",
          &TaoBuildingWallSectionStruct::name,
          &TaoBuildingWallSectionStruct::set_name
      )
      // TaoBuildingWallSectionStruct.constraint (0D_NOT_character - 'left_side' or 'right_side'
      // constraint.
      .def_property(
          "constraint",
          &TaoBuildingWallSectionStruct::constraint,
          &TaoBuildingWallSectionStruct::set_constraint
      )
      // TaoBuildingWallSectionStruct.point (1D_ALLOC_type -
      .def_property_readonly("point", &TaoBuildingWallSectionStruct::point)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoBuildingWallSectionStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoBuildingWallSectionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallSectionStruct &self) {
            return TaoBuildingWallSectionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallSectionStruct &self, py::dict &memo) {
            return TaoBuildingWallSectionStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoBuildingWallSectionStructArray1D>(m, "TaoBuildingWallSectionStructArray1D");
  bind_FTypeAlloc1D<TaoBuildingWallSectionStructAlloc1D>(m, "TaoBuildingWallSectionStructAlloc1D");
  // 2D TaoBuildingWallSectionStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallSectionStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_struct
void init_tao_building_wall_struct(py::module &m, py::class_<TaoBuildingWallStruct> &cls) {
  cls.def(py::init<>())
      // TaoBuildingWallStruct.orientation (0D_NOT_type -
      .def_property(
          "orientation",
          &TaoBuildingWallStruct::orientation,
          &TaoBuildingWallStruct::set_orientation
      )
      // TaoBuildingWallStruct.section (1D_ALLOC_type -
      .def_property_readonly("section", &TaoBuildingWallStruct::section)

      .def("__repr__", [](const TaoBuildingWallStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallStruct &self) {
            return TaoBuildingWallStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallStruct &self, py::dict &memo) {
            return TaoBuildingWallStruct(self);
          }
      )

      ;

  // 1D TaoBuildingWallStruct arrays are not used in structs/routines
  // 2D TaoBuildingWallStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_cmd_history_struct
void init_tao_cmd_history_struct(py::module &m, py::class_<TaoCmdHistoryStruct> &cls) {
  cls.def(py::init<>())
      // TaoCmdHistoryStruct.cmd (0D_ALLOC_character - The command
      .def_property("cmd", &TaoCmdHistoryStruct::cmd, &TaoCmdHistoryStruct::set_cmd)
      // TaoCmdHistoryStruct.ix (0D_NOT_integer - Command index (1st command has ix = 1, etc.) Note:
      // Commands from command files will be assigned an index.
      .def_property("ix", &TaoCmdHistoryStruct::ix, &TaoCmdHistoryStruct::set_ix)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoCmdHistoryStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoCmdHistoryStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCmdHistoryStruct &self) {
            return TaoCmdHistoryStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCmdHistoryStruct &self, py::dict &memo) { return TaoCmdHistoryStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoCmdHistoryStructArray1D>(m, "TaoCmdHistoryStructArray1D");
  bind_FTypeAlloc1D<TaoCmdHistoryStructAlloc1D>(m, "TaoCmdHistoryStructAlloc1D");
  // 2D TaoCmdHistoryStruct arrays are not used in structs/routines
  // 3D TaoCmdHistoryStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_common_struct
void init_tao_common_struct(py::module &m, py::class_<TaoCommonStruct> &cls) {
  cls.def(py::init<>())
      // TaoCommonStruct.plot_place_buffer (1D_ALLOC_type - Used when %external_plotting is on.
      .def_property_readonly("plot_place_buffer", &TaoCommonStruct::plot_place_buffer)
      // TaoCommonStruct.covar (2D_ALLOC_real -
      .def_property("covar", &TaoCommonStruct::covar, &TaoCommonStruct::set_covar)
      // TaoCommonStruct.alpha (2D_ALLOC_real -
      .def_property("alpha", &TaoCommonStruct::alpha, &TaoCommonStruct::set_alpha)
      // TaoCommonStruct.dummy_target (0D_NOT_real - Dummy varaible
      .def_property(
          "dummy_target",
          &TaoCommonStruct::dummy_target,
          &TaoCommonStruct::set_dummy_target
      )
      // TaoCommonStruct.n_alias (0D_NOT_integer -
      .def_property("n_alias", &TaoCommonStruct::n_alias, &TaoCommonStruct::set_n_alias)
      // TaoCommonStruct.cmd_file_level (0D_NOT_integer - For nested command files. 0 -> no command
      // file.
      .def_property(
          "cmd_file_level",
          &TaoCommonStruct::cmd_file_level,
          &TaoCommonStruct::set_cmd_file_level
      )
      // TaoCommonStruct.ix_key_bank (0D_NOT_integer - For single mode.
      .def_property("ix_key_bank", &TaoCommonStruct::ix_key_bank, &TaoCommonStruct::set_ix_key_bank)
      // TaoCommonStruct.ix_history (0D_NOT_integer - Index to latest command in the history
      // circular buffer.
      .def_property("ix_history", &TaoCommonStruct::ix_history, &TaoCommonStruct::set_ix_history)
      // TaoCommonStruct.n_history (0D_NOT_integer - Number of commands issued from beginning of
      // starting Tao.
      .def_property("n_history", &TaoCommonStruct::n_history, &TaoCommonStruct::set_n_history)
      // TaoCommonStruct.lev_loop (0D_NOT_integer - in do loop nest level
      .def_property("lev_loop", &TaoCommonStruct::lev_loop, &TaoCommonStruct::set_lev_loop)
      // TaoCommonStruct.n_err_messages_printed (0D_NOT_integer - Used by tao_set_invalid to limit
      // number of messages.
      .def_property(
          "n_err_messages_printed",
          &TaoCommonStruct::n_err_messages_printed,
          &TaoCommonStruct::set_n_err_messages_printed
      )
      // TaoCommonStruct.n_universes (0D_NOT_integer -
      .def_property("n_universes", &TaoCommonStruct::n_universes, &TaoCommonStruct::set_n_universes)
      // TaoCommonStruct.ix_beam_track_active_element (0D_NOT_integer - Element being tracked
      // through `tao_beam_track`.
      .def_property(
          "ix_beam_track_active_element",
          &TaoCommonStruct::ix_beam_track_active_element,
          &TaoCommonStruct::set_ix_beam_track_active_element
      )
      // TaoCommonStruct.cmd_file_paused (0D_NOT_logical -
      .def_property(
          "cmd_file_paused",
          &TaoCommonStruct::cmd_file_paused,
          &TaoCommonStruct::set_cmd_file_paused
      )
      // TaoCommonStruct.use_cmd_here (0D_NOT_logical - Used for commands recalled from the cmd
      // history stack
      .def_property(
          "use_cmd_here",
          &TaoCommonStruct::use_cmd_here,
          &TaoCommonStruct::set_use_cmd_here
      )
      // TaoCommonStruct.cmd_from_cmd_file (0D_NOT_logical - was command from a command file?
      .def_property(
          "cmd_from_cmd_file",
          &TaoCommonStruct::cmd_from_cmd_file,
          &TaoCommonStruct::set_cmd_from_cmd_file
      )
      // TaoCommonStruct.use_saved_beam_in_tracking (0D_NOT_logical -
      .def_property(
          "use_saved_beam_in_tracking",
          &TaoCommonStruct::use_saved_beam_in_tracking,
          &TaoCommonStruct::set_use_saved_beam_in_tracking
      )
      // TaoCommonStruct.single_mode (0D_NOT_logical -
      .def_property("single_mode", &TaoCommonStruct::single_mode, &TaoCommonStruct::set_single_mode)
      // TaoCommonStruct.combine_consecutive_elements_of_like_name (0D_NOT_logical -
      .def_property(
          "combine_consecutive_elements_of_like_name",
          &TaoCommonStruct::combine_consecutive_elements_of_like_name,
          &TaoCommonStruct::set_combine_consecutive_elements_of_like_name
      )
      // TaoCommonStruct.have_tracked_beam (0D_NOT_logical - Used to catch error when beam plotting
      // without having tracked a beam.
      .def_property(
          "have_tracked_beam",
          &TaoCommonStruct::have_tracked_beam,
          &TaoCommonStruct::set_have_tracked_beam
      )
      // TaoCommonStruct.init_plot_needed (0D_NOT_logical - reinitialize plotting?
      .def_property(
          "init_plot_needed",
          &TaoCommonStruct::init_plot_needed,
          &TaoCommonStruct::set_init_plot_needed
      )
      // TaoCommonStruct.init_beam (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property("init_beam", &TaoCommonStruct::init_beam, &TaoCommonStruct::set_init_beam)
      // TaoCommonStruct.init_var (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property("init_var", &TaoCommonStruct::init_var, &TaoCommonStruct::set_init_var)
      // TaoCommonStruct.init_read_lat_info (0D_NOT_logical - Used by custom programs to control Tao
      // init
      .def_property(
          "init_read_lat_info",
          &TaoCommonStruct::init_read_lat_info,
          &TaoCommonStruct::set_init_read_lat_info
      )
      // TaoCommonStruct.optimizer_running (0D_NOT_logical -
      .def_property(
          "optimizer_running",
          &TaoCommonStruct::optimizer_running,
          &TaoCommonStruct::set_optimizer_running
      )
      // TaoCommonStruct.have_datums_using_expressions (0D_NOT_logical -
      .def_property(
          "have_datums_using_expressions",
          &TaoCommonStruct::have_datums_using_expressions,
          &TaoCommonStruct::set_have_datums_using_expressions
      )
      // TaoCommonStruct.print_to_terminal (0D_NOT_logical - Print command prompt to the terminal?
      // For use with GUIs.
      .def_property(
          "print_to_terminal",
          &TaoCommonStruct::print_to_terminal,
          &TaoCommonStruct::set_print_to_terminal
      )
      // TaoCommonStruct.lattice_calc_done (0D_NOT_logical - Used by GUI for deciding when to
      // refresh.
      .def_property(
          "lattice_calc_done",
          &TaoCommonStruct::lattice_calc_done,
          &TaoCommonStruct::set_lattice_calc_done
      )
      // TaoCommonStruct.add_measurement_noise (0D_NOT_logical - Turn off to take data derivatives.
      .def_property(
          "add_measurement_noise",
          &TaoCommonStruct::add_measurement_noise,
          &TaoCommonStruct::set_add_measurement_noise
      )
      // TaoCommonStruct.is_err_message_printed (1D_NOT_logical - Used by tao_set_invalid
      .def_property(
          "is_err_message_printed",
          &TaoCommonStruct::is_err_message_printed,
          &TaoCommonStruct::set_is_err_message_printed
      )
      // TaoCommonStruct.command_arg_has_been_executed (0D_NOT_logical - Has the -command command
      // line argument been executed?
      .def_property(
          "command_arg_has_been_executed",
          &TaoCommonStruct::command_arg_has_been_executed,
          &TaoCommonStruct::set_command_arg_has_been_executed
      )
      // TaoCommonStruct.all_merit_weights_positive (0D_NOT_logical -
      .def_property(
          "all_merit_weights_positive",
          &TaoCommonStruct::all_merit_weights_positive,
          &TaoCommonStruct::set_all_merit_weights_positive
      )
      // TaoCommonStruct.multi_turn_orbit_is_plotted (0D_NOT_logical - Is a multi_turn_orbit being
      // plotted?
      .def_property(
          "multi_turn_orbit_is_plotted",
          &TaoCommonStruct::multi_turn_orbit_is_plotted,
          &TaoCommonStruct::set_multi_turn_orbit_is_plotted
      )
      // TaoCommonStruct.force_chrom_calc (0D_NOT_logical - Used by a routine to force a single
      // chromaticity calculation.
      .def_property(
          "force_chrom_calc",
          &TaoCommonStruct::force_chrom_calc,
          &TaoCommonStruct::set_force_chrom_calc
      )
      // TaoCommonStruct.force_rad_int_calc (0D_NOT_logical - Used by a routine to force a single
      // radiation integrals calculation
      .def_property(
          "force_rad_int_calc",
          &TaoCommonStruct::force_rad_int_calc,
          &TaoCommonStruct::set_force_rad_int_calc
      )
      // TaoCommonStruct.rad_int_ri_calc_on (0D_NOT_logical - 'Classical' radiation integrals
      // calculation on/off.
      .def_property(
          "rad_int_ri_calc_on",
          &TaoCommonStruct::rad_int_ri_calc_on,
          &TaoCommonStruct::set_rad_int_ri_calc_on
      )
      // TaoCommonStruct.rad_int_6d_calc_on (0D_NOT_logical - 6D Radiation integrals calculation
      // on/off.
      .def_property(
          "rad_int_6d_calc_on",
          &TaoCommonStruct::rad_int_6d_calc_on,
          &TaoCommonStruct::set_rad_int_6d_calc_on
      )
      // TaoCommonStruct.valid_plot_who (1D_NOT_character - model, base, ref etc...
      .def_property_readonly("valid_plot_who", &TaoCommonStruct::valid_plot_who)
      // TaoCommonStruct.single_mode_buffer (0D_NOT_character -
      .def_property(
          "single_mode_buffer",
          &TaoCommonStruct::single_mode_buffer,
          &TaoCommonStruct::set_single_mode_buffer
      )
      // TaoCommonStruct.cmd (0D_NOT_character - Used for the cmd history
      .def_property("cmd", &TaoCommonStruct::cmd, &TaoCommonStruct::set_cmd)
      // TaoCommonStruct.saved_cmd_line (0D_NOT_character - Saved part of command line when there
      // are mulitple commands on a line
      .def_property(
          "saved_cmd_line",
          &TaoCommonStruct::saved_cmd_line,
          &TaoCommonStruct::set_saved_cmd_line
      )

      .def("__repr__", [](const TaoCommonStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCommonStruct &self) {
            return TaoCommonStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCommonStruct &self, py::dict &memo) { return TaoCommonStruct(self); }
      )

      ;

  // 1D TaoCommonStruct arrays are not used in structs/routines
  // 2D TaoCommonStruct arrays are not used in structs/routines
  // 3D TaoCommonStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_color_struct
void init_tao_curve_color_struct(py::module &m, py::class_<TaoCurveColorStruct> &cls) {
  cls.def(py::init<>())
      // TaoCurveColorStruct.data_type (0D_NOT_character - Datum type to use for z-axis.
      .def_property(
          "data_type",
          &TaoCurveColorStruct::data_type,
          &TaoCurveColorStruct::set_data_type
      )
      // TaoCurveColorStruct.is_on (0D_NOT_logical - On/Off
      .def_property("is_on", &TaoCurveColorStruct::is_on, &TaoCurveColorStruct::set_is_on)
      // TaoCurveColorStruct.min (0D_NOT_real - Min and max values for mapping z-axis to color.
      .def_property("min", &TaoCurveColorStruct::min, &TaoCurveColorStruct::set_min)
      // TaoCurveColorStruct.max (0D_NOT_real - Min and max values for mapping z-axis to color.
      .def_property("max", &TaoCurveColorStruct::max, &TaoCurveColorStruct::set_max)
      // TaoCurveColorStruct.autoscale (0D_NOT_logical - Set %min, %max automatically to the limits
      // of %data_type
      .def_property(
          "autoscale",
          &TaoCurveColorStruct::autoscale,
          &TaoCurveColorStruct::set_autoscale
      )

      .def("__repr__", [](const TaoCurveColorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveColorStruct &self) {
            return TaoCurveColorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCurveColorStruct &self, py::dict &memo) { return TaoCurveColorStruct(self); }
      )

      ;

  // 1D TaoCurveColorStruct arrays are not used in structs/routines
  // 2D TaoCurveColorStruct arrays are not used in structs/routines
  // 3D TaoCurveColorStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_orbit_struct
void init_tao_curve_orbit_struct(py::module &m, py::class_<TaoCurveOrbitStruct> &cls) {
  cls.def(py::init<>())
      // TaoCurveOrbitStruct.x (0D_NOT_real - Transverse offset
      .def_property("x", &TaoCurveOrbitStruct::x, &TaoCurveOrbitStruct::set_x)
      // TaoCurveOrbitStruct.y (0D_NOT_real - Transverse offset
      .def_property("y", &TaoCurveOrbitStruct::y, &TaoCurveOrbitStruct::set_y)
      // TaoCurveOrbitStruct.t (0D_NOT_real - Time
      .def_property("t", &TaoCurveOrbitStruct::t, &TaoCurveOrbitStruct::set_t)

      .def("__repr__", [](const TaoCurveOrbitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveOrbitStruct &self) {
            return TaoCurveOrbitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCurveOrbitStruct &self, py::dict &memo) { return TaoCurveOrbitStruct(self); }
      )

      ;

  // 1D TaoCurveOrbitStruct arrays are not used in structs/routines
  // 2D TaoCurveOrbitStruct arrays are not used in structs/routines
  // 3D TaoCurveOrbitStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_struct
void init_tao_curve_struct(py::module &m, py::class_<TaoCurveStruct> &cls) {
  cls.def(py::init<>())
      // TaoCurveStruct.name (0D_NOT_character - Name identifying the curve.
      .def_property("name", &TaoCurveStruct::name, &TaoCurveStruct::set_name)
      // TaoCurveStruct.data_source (0D_NOT_character - 'lat', 'beam', 'data' (deprecated: 'dat'),
      // 'var', 'multi_turn_orbit'
      .def_property("data_source", &TaoCurveStruct::data_source, &TaoCurveStruct::set_data_source)
      // TaoCurveStruct.data_index (0D_NOT_character - Used for calculating %ix_symb(:).
      .def_property("data_index", &TaoCurveStruct::data_index, &TaoCurveStruct::set_data_index)
      // TaoCurveStruct.data_type_x (0D_NOT_character - Used for data slices and phase space plots.
      .def_property("data_type_x", &TaoCurveStruct::data_type_x, &TaoCurveStruct::set_data_type_x)
      // TaoCurveStruct.data_type (0D_ALLOC_character - 'orbit.x', etc.
      .def_property("data_type", &TaoCurveStruct::data_type, &TaoCurveStruct::set_data_type)
      // TaoCurveStruct.ele_ref_name (0D_NOT_character - Reference element.
      .def_property(
          "ele_ref_name",
          &TaoCurveStruct::ele_ref_name,
          &TaoCurveStruct::set_ele_ref_name
      )
      // TaoCurveStruct.legend_text (0D_NOT_character - String to draw in a curve legend.
      .def_property("legend_text", &TaoCurveStruct::legend_text, &TaoCurveStruct::set_legend_text)
      // TaoCurveStruct.message_text (0D_NOT_character - Informational message to draw with graph.
      .def_property(
          "message_text",
          &TaoCurveStruct::message_text,
          &TaoCurveStruct::set_message_text
      )
      // TaoCurveStruct.component (0D_NOT_character - Who to plot. Eg: 'meas - design'
      .def_property("component", &TaoCurveStruct::component, &TaoCurveStruct::set_component)
      // TaoCurveStruct.why_invalid (0D_NOT_character - Informative string to print.
      .def_property("why_invalid", &TaoCurveStruct::why_invalid, &TaoCurveStruct::set_why_invalid)
      // TaoCurveStruct.g (0D_PTR_type - pointer to parent graph
      .def_property("g", &TaoCurveStruct::g, &TaoCurveStruct::set_g)
      // TaoCurveStruct.hist (0D_NOT_type -
      .def_property("hist", &TaoCurveStruct::hist, &TaoCurveStruct::set_hist)
      // TaoCurveStruct.z_color (0D_NOT_type -
      .def_property("z_color", &TaoCurveStruct::z_color, &TaoCurveStruct::set_z_color)
      // TaoCurveStruct.x_line (1D_ALLOC_real - Coords for drawing a curve
      .def_property("x_line", &TaoCurveStruct::x_line, &TaoCurveStruct::set_x_line)
      // TaoCurveStruct.y_line (1D_ALLOC_real -
      .def_property("y_line", &TaoCurveStruct::y_line, &TaoCurveStruct::set_y_line)
      // TaoCurveStruct.y2_line (1D_ALLOC_real - Second array needed for beam chamber curve.
      .def_property("y2_line", &TaoCurveStruct::y2_line, &TaoCurveStruct::set_y2_line)
      // TaoCurveStruct.ix_line (1D_ALLOC_integer - Used by wave and aperture curves.
      .def_property("ix_line", &TaoCurveStruct::ix_line, &TaoCurveStruct::set_ix_line)
      // TaoCurveStruct.x_symb (1D_ALLOC_real - Coords for drawing the symbols
      .def_property("x_symb", &TaoCurveStruct::x_symb, &TaoCurveStruct::set_x_symb)
      // TaoCurveStruct.y_symb (1D_ALLOC_real -
      .def_property("y_symb", &TaoCurveStruct::y_symb, &TaoCurveStruct::set_y_symb)
      // TaoCurveStruct.z_symb (1D_ALLOC_real - Symbol color
      .def_property("z_symb", &TaoCurveStruct::z_symb, &TaoCurveStruct::set_z_symb)
      // TaoCurveStruct.err_symb (1D_ALLOC_real - Error bars
      .def_property("err_symb", &TaoCurveStruct::err_symb, &TaoCurveStruct::set_err_symb)
      // TaoCurveStruct.symb_size (1D_ALLOC_real - Symbol size. Used with symbol_size_scale.
      .def_property("symb_size", &TaoCurveStruct::symb_size, &TaoCurveStruct::set_symb_size)
      // TaoCurveStruct.ix_symb (1D_ALLOC_integer - Corresponding index in d1_data%d(:) array.
      .def_property("ix_symb", &TaoCurveStruct::ix_symb, &TaoCurveStruct::set_ix_symb)
      // TaoCurveStruct.y_axis_scale_factor (0D_NOT_real - y-axis conversion from internal to
      // plotting units.
      .def_property(
          "y_axis_scale_factor",
          &TaoCurveStruct::y_axis_scale_factor,
          &TaoCurveStruct::set_y_axis_scale_factor
      )
      // TaoCurveStruct.line (0D_NOT_type - Line attributes
      .def_property("line", &TaoCurveStruct::line, &TaoCurveStruct::set_line)
      // TaoCurveStruct.symbol (0D_NOT_type - Symbol attributes
      .def_property("symbol", &TaoCurveStruct::symbol, &TaoCurveStruct::set_symbol)
      // TaoCurveStruct.orbit (0D_NOT_type - Used for E/B field plotting.
      .def_property("orbit", &TaoCurveStruct::orbit, &TaoCurveStruct::set_orbit)
      // TaoCurveStruct.ix_universe (0D_NOT_integer - Universe where data is. -1 => use
      // s%global%default_universe
      .def_property("ix_universe", &TaoCurveStruct::ix_universe, &TaoCurveStruct::set_ix_universe)
      // TaoCurveStruct.symbol_every (0D_NOT_integer - Symbol every how many points.
      .def_property(
          "symbol_every",
          &TaoCurveStruct::symbol_every,
          &TaoCurveStruct::set_symbol_every
      )
      // TaoCurveStruct.ix_branch (0D_NOT_integer -
      .def_property("ix_branch", &TaoCurveStruct::ix_branch, &TaoCurveStruct::set_ix_branch)
      // TaoCurveStruct.ix_bunch (0D_NOT_integer - Bunch to plot.
      .def_property("ix_bunch", &TaoCurveStruct::ix_bunch, &TaoCurveStruct::set_ix_bunch)
      // TaoCurveStruct.n_turn (0D_NOT_integer - Used for multi_turn_orbit plotting
      .def_property("n_turn", &TaoCurveStruct::n_turn, &TaoCurveStruct::set_n_turn)
      // TaoCurveStruct.use_y2 (0D_NOT_logical - Use y2 axis?
      .def_property("use_y2", &TaoCurveStruct::use_y2, &TaoCurveStruct::set_use_y2)
      // TaoCurveStruct.draw_line (0D_NOT_logical - Draw a line through the data points?
      .def_property("draw_line", &TaoCurveStruct::draw_line, &TaoCurveStruct::set_draw_line)
      // TaoCurveStruct.draw_symbols (0D_NOT_logical - Draw a symbol at the data points?
      .def_property(
          "draw_symbols",
          &TaoCurveStruct::draw_symbols,
          &TaoCurveStruct::set_draw_symbols
      )
      // TaoCurveStruct.draw_symbol_index (0D_NOT_logical - Draw the symbol index number
      // curve%ix_symb?
      .def_property(
          "draw_symbol_index",
          &TaoCurveStruct::draw_symbol_index,
          &TaoCurveStruct::set_draw_symbol_index
      )
      // TaoCurveStruct.draw_error_bars (0D_NOT_logical - Draw error bars based upon data%error_rms
      // if drawing data? !! logical :: draw_rms = .false.          ! Show mean and RMS values with
      // legend?
      .def_property(
          "draw_error_bars",
          &TaoCurveStruct::draw_error_bars,
          &TaoCurveStruct::set_draw_error_bars
      )
      // TaoCurveStruct.smooth_line_calc (0D_NOT_logical - Calculate data between element edge
      // points?
      .def_property(
          "smooth_line_calc",
          &TaoCurveStruct::smooth_line_calc,
          &TaoCurveStruct::set_smooth_line_calc
      )
      // TaoCurveStruct.valid (0D_NOT_logical - valid data?
      .def_property("valid", &TaoCurveStruct::valid, &TaoCurveStruct::set_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoCurveStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoCurveStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveStruct &self) {
            return TaoCurveStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCurveStruct &self, py::dict &memo) { return TaoCurveStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoCurveStructArray1D>(m, "TaoCurveStructArray1D");
  bind_FTypeAlloc1D<TaoCurveStructAlloc1D>(m, "TaoCurveStructAlloc1D");
  // 2D TaoCurveStruct arrays are not used in structs/routines
  // 3D TaoCurveStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d1_data_struct
void init_tao_d1_data_struct(py::module &m, py::class_<TaoD1DataStruct> &cls) {
  cls.def(py::init<>())
      // TaoD1DataStruct.name (0D_NOT_character - Eg: 'x', etc.
      .def_property("name", &TaoD1DataStruct::name, &TaoD1DataStruct::set_name)
      // TaoD1DataStruct.d2 (0D_PTR_type - ptr to parent d2_data
      .def_property("d2", &TaoD1DataStruct::d2, &TaoD1DataStruct::set_d2)
      // TaoD1DataStruct.d (1D_PTR_type - Pointer to the appropriate section in u%data
      .def_property_readonly("d", &TaoD1DataStruct::d)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoD1DataStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoD1DataStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD1DataStruct &self) {
            return TaoD1DataStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoD1DataStruct &self, py::dict &memo) { return TaoD1DataStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoD1DataStructArray1D>(m, "TaoD1DataStructArray1D");
  bind_FTypeAlloc1D<TaoD1DataStructAlloc1D>(m, "TaoD1DataStructAlloc1D");
  // 2D TaoD1DataStruct arrays are not used in structs/routines
  // 3D TaoD1DataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d2_data_struct
void init_tao_d2_data_struct(py::module &m, py::class_<TaoD2DataStruct> &cls) {
  cls.def(py::init<>())
      // TaoD2DataStruct.name (0D_NOT_character - Name to be used with commands.
      .def_property("name", &TaoD2DataStruct::name, &TaoD2DataStruct::set_name)
      // TaoD2DataStruct.data_file_name (0D_NOT_character - Data file name .
      .def_property(
          "data_file_name",
          &TaoD2DataStruct::data_file_name,
          &TaoD2DataStruct::set_data_file_name
      )
      // TaoD2DataStruct.ref_file_name (0D_NOT_character - Reference file name.
      .def_property(
          "ref_file_name",
          &TaoD2DataStruct::ref_file_name,
          &TaoD2DataStruct::set_ref_file_name
      )
      // TaoD2DataStruct.data_date (0D_NOT_character - Data measurement date.
      .def_property("data_date", &TaoD2DataStruct::data_date, &TaoD2DataStruct::set_data_date)
      // TaoD2DataStruct.ref_date (0D_NOT_character - Reference data measurement date.
      .def_property("ref_date", &TaoD2DataStruct::ref_date, &TaoD2DataStruct::set_ref_date)
      // TaoD2DataStruct.descrip (1D_NOT_character - Array for descriptive information.
      .def_property_readonly("descrip", &TaoD2DataStruct::descrip)
      // TaoD2DataStruct.d1 (1D_ALLOC_type - Points to children
      .def_property_readonly("d1", &TaoD2DataStruct::d1)
      // TaoD2DataStruct.ix_universe (0D_NOT_integer - Index of universe this is in.
      .def_property("ix_universe", &TaoD2DataStruct::ix_universe, &TaoD2DataStruct::set_ix_universe)
      // TaoD2DataStruct.ix_d2_data (0D_NOT_integer - Index in u%d2_data(:) array.
      .def_property("ix_d2_data", &TaoD2DataStruct::ix_d2_data, &TaoD2DataStruct::set_ix_d2_data)
      // TaoD2DataStruct.ix_ref (0D_NOT_integer - Index of the reference data set.
      .def_property("ix_ref", &TaoD2DataStruct::ix_ref, &TaoD2DataStruct::set_ix_ref)
      // TaoD2DataStruct.data_read_in (0D_NOT_logical - A data set has been read in?
      .def_property(
          "data_read_in",
          &TaoD2DataStruct::data_read_in,
          &TaoD2DataStruct::set_data_read_in
      )
      // TaoD2DataStruct.ref_read_in (0D_NOT_logical - A reference data set has been read in?
      .def_property("ref_read_in", &TaoD2DataStruct::ref_read_in, &TaoD2DataStruct::set_ref_read_in)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoD2DataStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoD2DataStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD2DataStruct &self) {
            return TaoD2DataStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoD2DataStruct &self, py::dict &memo) { return TaoD2DataStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoD2DataStructArray1D>(m, "TaoD2DataStructArray1D");
  bind_FTypeAlloc1D<TaoD2DataStructAlloc1D>(m, "TaoD2DataStructAlloc1D");
  // 2D TaoD2DataStruct arrays are not used in structs/routines
  // 3D TaoD2DataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_data_struct
void init_tao_data_struct(py::module &m, py::class_<TaoDataStruct> &cls) {
  cls.def(py::init<>())
      // TaoDataStruct.ele_name (0D_NOT_character - Name of the lattice element where datum is
      // evaluated.
      .def_property("ele_name", &TaoDataStruct::ele_name, &TaoDataStruct::set_ele_name)
      // TaoDataStruct.ele_start_name (0D_NOT_character - Name of starting lattice element when
      // there is a range
      .def_property(
          "ele_start_name",
          &TaoDataStruct::ele_start_name,
          &TaoDataStruct::set_ele_start_name
      )
      // TaoDataStruct.ele_ref_name (0D_NOT_character - Name of reference lattice element
      .def_property("ele_ref_name", &TaoDataStruct::ele_ref_name, &TaoDataStruct::set_ele_ref_name)
      // TaoDataStruct.data_type (0D_ALLOC_character - Type of data: 'orbit.x', etc.
      .def_property("data_type", &TaoDataStruct::data_type, &TaoDataStruct::set_data_type)
      // TaoDataStruct.merit_type (0D_NOT_character - Type of constraint: 'target', 'max', 'min',
      // etc.
      .def_property("merit_type", &TaoDataStruct::merit_type, &TaoDataStruct::set_merit_type)
      // TaoDataStruct.id (0D_NOT_character - Used by Tao extension code. Not used by Tao directly.
      .def_property("id", &TaoDataStruct::id, &TaoDataStruct::set_id)
      // TaoDataStruct.data_source (0D_NOT_character - 'lat', 'beam', 'data' or 'var'. Last two used
      // for expressions.
      .def_property("data_source", &TaoDataStruct::data_source, &TaoDataStruct::set_data_source)
      // TaoDataStruct.why_invalid (0D_NOT_character - Informational string if there is a problem.
      .def_property("why_invalid", &TaoDataStruct::why_invalid, &TaoDataStruct::set_why_invalid)
      // TaoDataStruct.ix_uni (0D_NOT_integer - Universe index of datum.
      .def_property("ix_uni", &TaoDataStruct::ix_uni, &TaoDataStruct::set_ix_uni)
      // TaoDataStruct.ix_bunch (0D_NOT_integer - Bunch number to get the data from.
      .def_property("ix_bunch", &TaoDataStruct::ix_bunch, &TaoDataStruct::set_ix_bunch)
      // TaoDataStruct.ix_branch (0D_NOT_integer - Index of the associated lattice branch.
      .def_property("ix_branch", &TaoDataStruct::ix_branch, &TaoDataStruct::set_ix_branch)
      // TaoDataStruct.ix_ele (0D_NOT_integer - Index of the lattice element corresponding to
      // ele_name
      .def_property("ix_ele", &TaoDataStruct::ix_ele, &TaoDataStruct::set_ix_ele)
      // TaoDataStruct.ix_ele_start (0D_NOT_integer - Index of lattice elment when there is a range
      .def_property("ix_ele_start", &TaoDataStruct::ix_ele_start, &TaoDataStruct::set_ix_ele_start)
      // TaoDataStruct.ix_ele_ref (0D_NOT_integer - Index of lattice elment when there is a
      // reference.
      .def_property("ix_ele_ref", &TaoDataStruct::ix_ele_ref, &TaoDataStruct::set_ix_ele_ref)
      // TaoDataStruct.ix_ele_merit (0D_NOT_integer - Index of lattice elment where merit is
      // evaluated.
      .def_property("ix_ele_merit", &TaoDataStruct::ix_ele_merit, &TaoDataStruct::set_ix_ele_merit)
      // TaoDataStruct.ix_d1 (0D_NOT_integer - Index number in u%d2_data(i)%d1_data(j)%d(:) array.
      .def_property("ix_d1", &TaoDataStruct::ix_d1, &TaoDataStruct::set_ix_d1)
      // TaoDataStruct.ix_data (0D_NOT_integer - Index of this datum in the u%data(:) array of
      // data_structs.
      .def_property("ix_data", &TaoDataStruct::ix_data, &TaoDataStruct::set_ix_data)
      // TaoDataStruct.ix_dModel (0D_NOT_integer - Row number in the dModel_dVar derivative matrix.
      .def_property("ix_dModel", &TaoDataStruct::ix_dModel, &TaoDataStruct::set_ix_dModel)
      // TaoDataStruct.eval_point (0D_NOT_integer - or anchor_center$, anchor_beginning$. Where to
      // evaluate data relative to the element.
      .def_property("eval_point", &TaoDataStruct::eval_point, &TaoDataStruct::set_eval_point)
      // TaoDataStruct.meas_value (0D_NOT_real - Measured datum value.
      .def_property("meas_value", &TaoDataStruct::meas_value, &TaoDataStruct::set_meas_value)
      // TaoDataStruct.ref_value (0D_NOT_real - Measured datum value from the reference data set.
      .def_property("ref_value", &TaoDataStruct::ref_value, &TaoDataStruct::set_ref_value)
      // TaoDataStruct.model_value (0D_NOT_real - Datum value as calculated from the model.
      .def_property("model_value", &TaoDataStruct::model_value, &TaoDataStruct::set_model_value)
      // TaoDataStruct.design_value (0D_NOT_real - What the datum value is in the design lattice.
      .def_property("design_value", &TaoDataStruct::design_value, &TaoDataStruct::set_design_value)
      // TaoDataStruct.old_value (0D_NOT_real - The model_value at some previous time.
      .def_property("old_value", &TaoDataStruct::old_value, &TaoDataStruct::set_old_value)
      // TaoDataStruct.base_value (0D_NOT_real - The value as calculated from the base model.
      .def_property("base_value", &TaoDataStruct::base_value, &TaoDataStruct::set_base_value)
      // TaoDataStruct.error_rms (0D_NOT_real - Measurement error RMS. Used in plotting.
      .def_property("error_rms", &TaoDataStruct::error_rms, &TaoDataStruct::set_error_rms)
      // TaoDataStruct.delta_merit (0D_NOT_real - Diff used to calculate the merit function term.
      .def_property("delta_merit", &TaoDataStruct::delta_merit, &TaoDataStruct::set_delta_merit)
      // TaoDataStruct.weight (0D_NOT_real - Weight for the merit function term.
      .def_property("weight", &TaoDataStruct::weight, &TaoDataStruct::set_weight)
      // TaoDataStruct.invalid_value (0D_NOT_real - Value used in merit calc if good_model = F (or
      // possibly good_design & good_base).
      .def_property(
          "invalid_value",
          &TaoDataStruct::invalid_value,
          &TaoDataStruct::set_invalid_value
      )
      // TaoDataStruct.merit (0D_NOT_real - Merit function term value: weight * delta_merit^2
      .def_property("merit", &TaoDataStruct::merit, &TaoDataStruct::set_merit)
      // TaoDataStruct.s (0D_NOT_real - longitudinal position of ele.
      .def_property("s", &TaoDataStruct::s, &TaoDataStruct::set_s)
      // TaoDataStruct.s_offset (0D_NOT_real - Offset of the evaluation point.
      .def_property("s_offset", &TaoDataStruct::s_offset, &TaoDataStruct::set_s_offset)
      // TaoDataStruct.ref_s_offset (0D_NOT_real - Offset of the reference point. In development.
      .def_property("ref_s_offset", &TaoDataStruct::ref_s_offset, &TaoDataStruct::set_ref_s_offset)
      // TaoDataStruct.err_message_printed (0D_NOT_logical - Used to prevent zillions of error
      // messages being generated
      .def_property(
          "err_message_printed",
          &TaoDataStruct::err_message_printed,
          &TaoDataStruct::set_err_message_printed
      )
      // TaoDataStruct.exists (0D_NOT_logical - See above
      .def_property("exists", &TaoDataStruct::exists, &TaoDataStruct::set_exists)
      // TaoDataStruct.good_model (0D_NOT_logical - See above
      .def_property("good_model", &TaoDataStruct::good_model, &TaoDataStruct::set_good_model)
      // TaoDataStruct.good_base (0D_NOT_logical - See above
      .def_property("good_base", &TaoDataStruct::good_base, &TaoDataStruct::set_good_base)
      // TaoDataStruct.good_design (0D_NOT_logical - See above
      .def_property("good_design", &TaoDataStruct::good_design, &TaoDataStruct::set_good_design)
      // TaoDataStruct.good_meas (0D_NOT_logical - See above
      .def_property("good_meas", &TaoDataStruct::good_meas, &TaoDataStruct::set_good_meas)
      // TaoDataStruct.good_ref (0D_NOT_logical - See above
      .def_property("good_ref", &TaoDataStruct::good_ref, &TaoDataStruct::set_good_ref)
      // TaoDataStruct.good_user (0D_NOT_logical - See above
      .def_property("good_user", &TaoDataStruct::good_user, &TaoDataStruct::set_good_user)
      // TaoDataStruct.good_opt (0D_NOT_logical - See above
      .def_property("good_opt", &TaoDataStruct::good_opt, &TaoDataStruct::set_good_opt)
      // TaoDataStruct.good_plot (0D_NOT_logical - See above
      .def_property("good_plot", &TaoDataStruct::good_plot, &TaoDataStruct::set_good_plot)
      // TaoDataStruct.useit_plot (0D_NOT_logical - See above
      .def_property("useit_plot", &TaoDataStruct::useit_plot, &TaoDataStruct::set_useit_plot)
      // TaoDataStruct.useit_opt (0D_NOT_logical - See above
      .def_property("useit_opt", &TaoDataStruct::useit_opt, &TaoDataStruct::set_useit_opt)
      // TaoDataStruct.spin_map (0D_NOT_type -
      .def_property("spin_map", &TaoDataStruct::spin_map, &TaoDataStruct::set_spin_map)
      // TaoDataStruct.d1 (0D_PTR_type - Pointer to the parent d1_data_struct
      .def_property("d1", &TaoDataStruct::d1, &TaoDataStruct::set_d1)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoDataStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoDataStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDataStruct &self) {
            return TaoDataStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDataStruct &self, py::dict &memo) { return TaoDataStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoDataStructArray1D>(m, "TaoDataStructArray1D");
  bind_FTypeAlloc1D<TaoDataStructAlloc1D>(m, "TaoDataStructAlloc1D");
  // 2D TaoDataStruct arrays are not used in structs/routines
  // 3D TaoDataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_data_var_component_struct
void init_tao_data_var_component_struct(py::module &m, py::class_<TaoDataVarComponentStruct> &cls) {
  cls.def(py::init<>())
      // TaoDataVarComponentStruct.name (0D_NOT_character - Eg: 'meas', 'ref', 'model', etc.
      .def_property("name", &TaoDataVarComponentStruct::name, &TaoDataVarComponentStruct::set_name)
      // TaoDataVarComponentStruct.sign (0D_NOT_real - +1 or -1
      .def_property("sign", &TaoDataVarComponentStruct::sign, &TaoDataVarComponentStruct::set_sign)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoDataVarComponentStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoDataVarComponentStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDataVarComponentStruct &self) {
            return TaoDataVarComponentStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDataVarComponentStruct &self, py::dict &memo) {
            return TaoDataVarComponentStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoDataVarComponentStructArray1D>(m, "TaoDataVarComponentStructArray1D");
  bind_FTypeAlloc1D<TaoDataVarComponentStructAlloc1D>(m, "TaoDataVarComponentStructAlloc1D");
  // 2D TaoDataVarComponentStruct arrays are not used in structs/routines
  // 3D TaoDataVarComponentStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_drawing_struct
void init_tao_drawing_struct(py::module &m, py::class_<TaoDrawingStruct> &cls) {
  cls.def(py::init<>())
      // TaoDrawingStruct.ele_shape (1D_ALLOC_type -
      .def_property_readonly("ele_shape", &TaoDrawingStruct::ele_shape)

      .def("__repr__", [](const TaoDrawingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDrawingStruct &self) {
            return TaoDrawingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDrawingStruct &self, py::dict &memo) { return TaoDrawingStruct(self); }
      )

      ;

  // 1D TaoDrawingStruct arrays are not used in structs/routines
  // 2D TaoDrawingStruct arrays are not used in structs/routines
  // 3D TaoDrawingStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_dynamic_aperture_struct
void init_tao_dynamic_aperture_struct(py::module &m, py::class_<TaoDynamicApertureStruct> &cls) {
  cls.def(py::init<>())
      // TaoDynamicApertureStruct.param (0D_NOT_type -
      .def_property("param", &TaoDynamicApertureStruct::param, &TaoDynamicApertureStruct::set_param)
      // TaoDynamicApertureStruct.scan (1D_ALLOC_type - One scan for each pz.
      .def_property_readonly("scan", &TaoDynamicApertureStruct::scan)
      // TaoDynamicApertureStruct.pz (1D_ALLOC_real -
      .def_property("pz", &TaoDynamicApertureStruct::pz, &TaoDynamicApertureStruct::set_pz)
      // TaoDynamicApertureStruct.ellipse_scale (0D_NOT_real -
      .def_property(
          "ellipse_scale",
          &TaoDynamicApertureStruct::ellipse_scale,
          &TaoDynamicApertureStruct::set_ellipse_scale
      )
      // TaoDynamicApertureStruct.a_emit (0D_NOT_real -
      .def_property(
          "a_emit",
          &TaoDynamicApertureStruct::a_emit,
          &TaoDynamicApertureStruct::set_a_emit
      )
      // TaoDynamicApertureStruct.b_emit (0D_NOT_real -
      .def_property(
          "b_emit",
          &TaoDynamicApertureStruct::b_emit,
          &TaoDynamicApertureStruct::set_b_emit
      )

      .def("__repr__", [](const TaoDynamicApertureStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDynamicApertureStruct &self) {
            return TaoDynamicApertureStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDynamicApertureStruct &self, py::dict &memo) {
            return TaoDynamicApertureStruct(self);
          }
      )

      ;

  // 1D TaoDynamicApertureStruct arrays are not used in structs/routines
  // 2D TaoDynamicApertureStruct arrays are not used in structs/routines
  // 3D TaoDynamicApertureStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_pointer_struct
void init_tao_ele_pointer_struct(py::module &m, py::class_<TaoElePointerStruct> &cls) {
  cls.def(py::init<>())
      // TaoElePointerStruct.eles (1D_ALLOC_type -
      .def_property_readonly("eles", &TaoElePointerStruct::eles)
      // TaoElePointerStruct.n_loc (0D_NOT_integer -
      .def_property("n_loc", &TaoElePointerStruct::n_loc, &TaoElePointerStruct::set_n_loc)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoElePointerStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoElePointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoElePointerStruct &self) {
            return TaoElePointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoElePointerStruct &self, py::dict &memo) { return TaoElePointerStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoElePointerStructArray1D>(m, "TaoElePointerStructArray1D");
  bind_FTypeAlloc1D<TaoElePointerStructAlloc1D>(m, "TaoElePointerStructAlloc1D");
  // 2D TaoElePointerStruct arrays are not used in structs/routines
  // 3D TaoElePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_shape_struct
void init_tao_ele_shape_struct(py::module &m, py::class_<TaoEleShapeStruct> &cls) {
  cls.def(py::init<>())
      // TaoEleShapeStruct.ele_id (0D_NOT_character - element 'key::name' to match to.
      .def_property("ele_id", &TaoEleShapeStruct::ele_id, &TaoEleShapeStruct::set_ele_id)
      // TaoEleShapeStruct.shape (0D_NOT_character - Shape to draw
      .def_property("shape", &TaoEleShapeStruct::shape, &TaoEleShapeStruct::set_shape)
      // TaoEleShapeStruct.color (0D_NOT_character - Color of shape
      .def_property("color", &TaoEleShapeStruct::color, &TaoEleShapeStruct::set_color)
      // TaoEleShapeStruct.size (0D_NOT_real - plot vertical height
      .def_property("size", &TaoEleShapeStruct::size, &TaoEleShapeStruct::set_size)
      // TaoEleShapeStruct.label (0D_NOT_character - Can be: 'name', 's', 'none'
      .def_property("label", &TaoEleShapeStruct::label, &TaoEleShapeStruct::set_label)
      // TaoEleShapeStruct.draw (0D_NOT_logical - Draw the shape?
      .def_property("draw", &TaoEleShapeStruct::draw, &TaoEleShapeStruct::set_draw)
      // TaoEleShapeStruct.multi (0D_NOT_logical - Can be part of a multi-shape.
      .def_property("multi", &TaoEleShapeStruct::multi, &TaoEleShapeStruct::set_multi)
      // TaoEleShapeStruct.line_width (0D_NOT_integer - Width of lines used to draw the shape.
      .def_property(
          "line_width",
          &TaoEleShapeStruct::line_width,
          &TaoEleShapeStruct::set_line_width
      )
      // TaoEleShapeStruct.offset (0D_NOT_real - Vertical offset.
      .def_property("offset", &TaoEleShapeStruct::offset, &TaoEleShapeStruct::set_offset)
      // TaoEleShapeStruct.ix_key (0D_NOT_integer - Extracted from ele_id. 0 => all classes
      // (quadrupole, etc.)
      .def_property("ix_key", &TaoEleShapeStruct::ix_key, &TaoEleShapeStruct::set_ix_key)
      // TaoEleShapeStruct.name_ele (0D_NOT_character - Name of element.
      .def_property("name_ele", &TaoEleShapeStruct::name_ele, &TaoEleShapeStruct::set_name_ele)
      // TaoEleShapeStruct.uni (1D_ALLOC_type -
      .def_property_readonly("uni", &TaoEleShapeStruct::uni)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoEleShapeStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoEleShapeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoEleShapeStruct &self) {
            return TaoEleShapeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoEleShapeStruct &self, py::dict &memo) { return TaoEleShapeStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoEleShapeStructArray1D>(m, "TaoEleShapeStructArray1D");
  bind_FTypeAlloc1D<TaoEleShapeStructAlloc1D>(m, "TaoEleShapeStructAlloc1D");
  // 2D TaoEleShapeStruct arrays are not used in structs/routines
  // 3D TaoEleShapeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_eval_node_struct
void init_tao_eval_node_struct(py::module &m, py::class_<TaoEvalNodeStruct> &cls) {
  cls.def(py::init<>())
      // TaoEvalNodeStruct.type (0D_NOT_integer -
      .def_property("type", &TaoEvalNodeStruct::type, &TaoEvalNodeStruct::set_type)
      // TaoEvalNodeStruct.name (0D_NOT_character -
      .def_property("name", &TaoEvalNodeStruct::name, &TaoEvalNodeStruct::set_name)
      // TaoEvalNodeStruct.scale (0D_NOT_real - Scale factor for ping data
      .def_property("scale", &TaoEvalNodeStruct::scale, &TaoEvalNodeStruct::set_scale)
      // TaoEvalNodeStruct.value (1D_ALLOC_real -
      .def_property("value", &TaoEvalNodeStruct::value, &TaoEvalNodeStruct::set_value)
      // TaoEvalNodeStruct.info (1D_ALLOC_type -
      .def_property_readonly("info", &TaoEvalNodeStruct::info)
      // TaoEvalNodeStruct.node (1D_PTR_type - Child nodes for tree construction.
      .def_property_readonly("node", &TaoEvalNodeStruct::node)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoEvalNodeStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoEvalNodeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoEvalNodeStruct &self) {
            return TaoEvalNodeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoEvalNodeStruct &self, py::dict &memo) { return TaoEvalNodeStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoEvalNodeStructArray1D>(m, "TaoEvalNodeStructArray1D");
  bind_FTypeAlloc1D<TaoEvalNodeStructAlloc1D>(m, "TaoEvalNodeStructAlloc1D");
  // 2D TaoEvalNodeStruct arrays are not used in structs/routines
  // 3D TaoEvalNodeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_expression_info_struct
void init_tao_expression_info_struct(py::module &m, py::class_<TaoExpressionInfoStruct> &cls) {
  cls.def(py::init<>())
      // TaoExpressionInfoStruct.good (0D_NOT_logical - Expression is valid.
      .def_property("good", &TaoExpressionInfoStruct::good, &TaoExpressionInfoStruct::set_good)
      // TaoExpressionInfoStruct.ele (0D_PTR_type - Associated ele if it exists
      .def_property("ele", &TaoExpressionInfoStruct::ele, &TaoExpressionInfoStruct::set_ele)
      // TaoExpressionInfoStruct.s (0D_NOT_real - Longitudinal position of expression.
      .def_property("s", &TaoExpressionInfoStruct::s, &TaoExpressionInfoStruct::set_s)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoExpressionInfoStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoExpressionInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoExpressionInfoStruct &self) {
            return TaoExpressionInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoExpressionInfoStruct &self, py::dict &memo) {
            return TaoExpressionInfoStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoExpressionInfoStructArray1D>(m, "TaoExpressionInfoStructArray1D");
  bind_FTypeAlloc1D<TaoExpressionInfoStructAlloc1D>(m, "TaoExpressionInfoStructAlloc1D");
  // 2D TaoExpressionInfoStruct arrays are not used in structs/routines
  // 3D TaoExpressionInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_floor_plan_struct
void init_tao_floor_plan_struct(py::module &m, py::class_<TaoFloorPlanStruct> &cls) {
  cls.def(py::init<>())
      // TaoFloorPlanStruct.view (0D_NOT_character - or 'xz'.
      .def_property("view", &TaoFloorPlanStruct::view, &TaoFloorPlanStruct::set_view)
      // TaoFloorPlanStruct.rotation (0D_NOT_real - Rotation of floor plan plot: 1.0 -> 360^deg
      .def_property("rotation", &TaoFloorPlanStruct::rotation, &TaoFloorPlanStruct::set_rotation)
      // TaoFloorPlanStruct.correct_distortion (0D_NOT_logical - T -> Shrink one axis so x-scale =
      // y-scale.
      .def_property(
          "correct_distortion",
          &TaoFloorPlanStruct::correct_distortion,
          &TaoFloorPlanStruct::set_correct_distortion
      )
      // TaoFloorPlanStruct.flip_label_side (0D_NOT_logical - Draw element label on other side of
      // element?
      .def_property(
          "flip_label_side",
          &TaoFloorPlanStruct::flip_label_side,
          &TaoFloorPlanStruct::set_flip_label_side
      )
      // TaoFloorPlanStruct.size_is_absolute (0D_NOT_logical - Are shape sizes in meters or window
      // pixels?
      .def_property(
          "size_is_absolute",
          &TaoFloorPlanStruct::size_is_absolute,
          &TaoFloorPlanStruct::set_size_is_absolute
      )
      // TaoFloorPlanStruct.draw_only_first_pass (0D_NOT_logical - Draw only first pass with
      // multipass elements?
      .def_property(
          "draw_only_first_pass",
          &TaoFloorPlanStruct::draw_only_first_pass,
          &TaoFloorPlanStruct::set_draw_only_first_pass
      )
      // TaoFloorPlanStruct.draw_building_wall (0D_NOT_logical - Draw the building wall?
      .def_property(
          "draw_building_wall",
          &TaoFloorPlanStruct::draw_building_wall,
          &TaoFloorPlanStruct::set_draw_building_wall
      )
      // TaoFloorPlanStruct.orbit_scale (0D_NOT_real - Scale factor for drawing orbits. 0 -> Do not
      // draw.
      .def_property(
          "orbit_scale",
          &TaoFloorPlanStruct::orbit_scale,
          &TaoFloorPlanStruct::set_orbit_scale
      )
      // TaoFloorPlanStruct.orbit_color (0D_NOT_character -
      .def_property(
          "orbit_color",
          &TaoFloorPlanStruct::orbit_color,
          &TaoFloorPlanStruct::set_orbit_color
      )
      // TaoFloorPlanStruct.orbit_pattern (0D_NOT_character -
      .def_property(
          "orbit_pattern",
          &TaoFloorPlanStruct::orbit_pattern,
          &TaoFloorPlanStruct::set_orbit_pattern
      )
      // TaoFloorPlanStruct.orbit_lattice (0D_NOT_character - Or 'design' or 'base'
      .def_property(
          "orbit_lattice",
          &TaoFloorPlanStruct::orbit_lattice,
          &TaoFloorPlanStruct::set_orbit_lattice
      )
      // TaoFloorPlanStruct.orbit_width (0D_NOT_integer -
      .def_property(
          "orbit_width",
          &TaoFloorPlanStruct::orbit_width,
          &TaoFloorPlanStruct::set_orbit_width
      )

      .def("__repr__", [](const TaoFloorPlanStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoFloorPlanStruct &self) {
            return TaoFloorPlanStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoFloorPlanStruct &self, py::dict &memo) { return TaoFloorPlanStruct(self); }
      )

      ;

  // 1D TaoFloorPlanStruct arrays are not used in structs/routines
  // 2D TaoFloorPlanStruct arrays are not used in structs/routines
  // 3D TaoFloorPlanStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_global_struct
void init_tao_global_struct(py::module &m, py::class_<TaoGlobalStruct> &cls) {
  cls.def(py::init<>())
      // TaoGlobalStruct.beam_dead_cutoff (0D_NOT_real - Percentage of dead particles at which beam
      // tracking is stopped.
      .def_property(
          "beam_dead_cutoff",
          &TaoGlobalStruct::beam_dead_cutoff,
          &TaoGlobalStruct::set_beam_dead_cutoff
      )
      // TaoGlobalStruct.lm_opt_deriv_reinit (0D_NOT_real - Reinit derivative matrix cutoff
      .def_property(
          "lm_opt_deriv_reinit",
          &TaoGlobalStruct::lm_opt_deriv_reinit,
          &TaoGlobalStruct::set_lm_opt_deriv_reinit
      )
      // TaoGlobalStruct.de_lm_step_ratio (0D_NOT_real - Scaling for step sizes between DE and LM
      // optimizers.
      .def_property(
          "de_lm_step_ratio",
          &TaoGlobalStruct::de_lm_step_ratio,
          &TaoGlobalStruct::set_de_lm_step_ratio
      )
      // TaoGlobalStruct.de_var_to_population_factor (0D_NOT_real - DE population =
      // max(n_var*factor, 20)
      .def_property(
          "de_var_to_population_factor",
          &TaoGlobalStruct::de_var_to_population_factor,
          &TaoGlobalStruct::set_de_var_to_population_factor
      )
      // TaoGlobalStruct.lmdif_eps (0D_NOT_real - Tollerance for lmdif optimizer.
      .def_property("lmdif_eps", &TaoGlobalStruct::lmdif_eps, &TaoGlobalStruct::set_lmdif_eps)
      // TaoGlobalStruct.lmdif_negligible_merit (0D_NOT_real -
      .def_property(
          "lmdif_negligible_merit",
          &TaoGlobalStruct::lmdif_negligible_merit,
          &TaoGlobalStruct::set_lmdif_negligible_merit
      )
      // TaoGlobalStruct.svd_cutoff (0D_NOT_real - SVD singular value cutoff.
      .def_property("svd_cutoff", &TaoGlobalStruct::svd_cutoff, &TaoGlobalStruct::set_svd_cutoff)
      // TaoGlobalStruct.unstable_penalty (0D_NOT_real - Used in unstable_ring datum merit
      // calculation.
      .def_property(
          "unstable_penalty",
          &TaoGlobalStruct::unstable_penalty,
          &TaoGlobalStruct::set_unstable_penalty
      )
      // TaoGlobalStruct.merit_stop_value (0D_NOT_real - Merit value below which an optimizer will
      // stop.
      .def_property(
          "merit_stop_value",
          &TaoGlobalStruct::merit_stop_value,
          &TaoGlobalStruct::set_merit_stop_value
      )
      // TaoGlobalStruct.dmerit_stop_value (0D_NOT_real - Fractional Merit change below which an
      // optimizer will stop.
      .def_property(
          "dmerit_stop_value",
          &TaoGlobalStruct::dmerit_stop_value,
          &TaoGlobalStruct::set_dmerit_stop_value
      )
      // TaoGlobalStruct.random_sigma_cutoff (0D_NOT_real - Cut-off in sigmas.
      .def_property(
          "random_sigma_cutoff",
          &TaoGlobalStruct::random_sigma_cutoff,
          &TaoGlobalStruct::set_random_sigma_cutoff
      )
      // TaoGlobalStruct.delta_e_chrom (0D_NOT_real - Delta E used from chrom calc.
      .def_property(
          "delta_e_chrom",
          &TaoGlobalStruct::delta_e_chrom,
          &TaoGlobalStruct::set_delta_e_chrom
      )
      // TaoGlobalStruct.max_plot_time (0D_NOT_real - If plotting time (seconds) exceeds this than a
      // message is generated.
      .def_property(
          "max_plot_time",
          &TaoGlobalStruct::max_plot_time,
          &TaoGlobalStruct::set_max_plot_time
      )
      // TaoGlobalStruct.default_universe (0D_NOT_integer - Default universe to work with.
      .def_property(
          "default_universe",
          &TaoGlobalStruct::default_universe,
          &TaoGlobalStruct::set_default_universe
      )
      // TaoGlobalStruct.default_branch (0D_NOT_integer - Default lattice branch to work with.
      .def_property(
          "default_branch",
          &TaoGlobalStruct::default_branch,
          &TaoGlobalStruct::set_default_branch
      )
      // TaoGlobalStruct.n_opti_cycles (0D_NOT_integer - Number of optimization cycles
      .def_property(
          "n_opti_cycles",
          &TaoGlobalStruct::n_opti_cycles,
          &TaoGlobalStruct::set_n_opti_cycles
      )
      // TaoGlobalStruct.n_opti_loops (0D_NOT_integer - Number of optimization loops
      .def_property(
          "n_opti_loops",
          &TaoGlobalStruct::n_opti_loops,
          &TaoGlobalStruct::set_n_opti_loops
      )
      // TaoGlobalStruct.n_threads (0D_NOT_integer - Number of OpenMP threads for parallel
      // calculations.
      .def_property("n_threads", &TaoGlobalStruct::n_threads, &TaoGlobalStruct::set_n_threads)
      // TaoGlobalStruct.phase_units (0D_NOT_integer - Phase units on output.
      .def_property("phase_units", &TaoGlobalStruct::phase_units, &TaoGlobalStruct::set_phase_units)
      // TaoGlobalStruct.bunch_to_plot (0D_NOT_integer - Which bunch to plot
      .def_property(
          "bunch_to_plot",
          &TaoGlobalStruct::bunch_to_plot,
          &TaoGlobalStruct::set_bunch_to_plot
      )
      // TaoGlobalStruct.random_seed (0D_NOT_integer - Use system clock by default
      .def_property("random_seed", &TaoGlobalStruct::random_seed, &TaoGlobalStruct::set_random_seed)
      // TaoGlobalStruct.n_top10_merit (0D_NOT_integer - Number of top merit constraints to print.
      .def_property(
          "n_top10_merit",
          &TaoGlobalStruct::n_top10_merit,
          &TaoGlobalStruct::set_n_top10_merit
      )
      // TaoGlobalStruct.srdt_gen_n_slices (0D_NOT_integer - Number times to slice elements for
      // summation RDT calculation
      .def_property(
          "srdt_gen_n_slices",
          &TaoGlobalStruct::srdt_gen_n_slices,
          &TaoGlobalStruct::set_srdt_gen_n_slices
      )
      // TaoGlobalStruct.datum_err_messages_max (0D_NOT_integer - Maximum number of error messages
      // per call to lattice_calc.
      .def_property(
          "datum_err_messages_max",
          &TaoGlobalStruct::datum_err_messages_max,
          &TaoGlobalStruct::set_datum_err_messages_max
      )
      // TaoGlobalStruct.srdt_sxt_n_slices (0D_NOT_integer - Number times to slice sextupoles for
      // summation RDT calculation
      .def_property(
          "srdt_sxt_n_slices",
          &TaoGlobalStruct::srdt_sxt_n_slices,
          &TaoGlobalStruct::set_srdt_sxt_n_slices
      )
      // TaoGlobalStruct.srdt_use_cache (0D_NOT_logical - Create cache for SRDT calculations.  Can
      // use lots of memory if srdt_*_n_slices large.
      .def_property(
          "srdt_use_cache",
          &TaoGlobalStruct::srdt_use_cache,
          &TaoGlobalStruct::set_srdt_use_cache
      )
      // TaoGlobalStruct.quiet (0D_NOT_character - Print I/O when running a command file?
      .def_property("quiet", &TaoGlobalStruct::quiet, &TaoGlobalStruct::set_quiet)
      // TaoGlobalStruct.random_engine (0D_NOT_character - Non-beam random number engine
      .def_property(
          "random_engine",
          &TaoGlobalStruct::random_engine,
          &TaoGlobalStruct::set_random_engine
      )
      // TaoGlobalStruct.random_gauss_converter (0D_NOT_character - Non-beam
      .def_property(
          "random_gauss_converter",
          &TaoGlobalStruct::random_gauss_converter,
          &TaoGlobalStruct::set_random_gauss_converter
      )
      // TaoGlobalStruct.track_type (0D_NOT_character - or 'beam'
      .def_property("track_type", &TaoGlobalStruct::track_type, &TaoGlobalStruct::set_track_type)
      // TaoGlobalStruct.lat_sigma_calc_uses_emit_from (0D_NOT_character - Lattice derived sigma
      // matrix uses emit values from where? Other possibilities: 'beam', 'beam_init'.
      .def_property(
          "lat_sigma_calc_uses_emit_from",
          &TaoGlobalStruct::lat_sigma_calc_uses_emit_from,
          &TaoGlobalStruct::set_lat_sigma_calc_uses_emit_from
      )
      // TaoGlobalStruct.prompt_string (0D_NOT_character -
      .def_property(
          "prompt_string",
          &TaoGlobalStruct::prompt_string,
          &TaoGlobalStruct::set_prompt_string
      )
      // TaoGlobalStruct.prompt_color (0D_NOT_character - See read_a_line routine for possible
      // settings.
      .def_property(
          "prompt_color",
          &TaoGlobalStruct::prompt_color,
          &TaoGlobalStruct::set_prompt_color
      )
      // TaoGlobalStruct.optimizer (0D_NOT_character - optimizer to use.
      .def_property("optimizer", &TaoGlobalStruct::optimizer, &TaoGlobalStruct::set_optimizer)
      // TaoGlobalStruct.print_command (0D_NOT_character -
      .def_property(
          "print_command",
          &TaoGlobalStruct::print_command,
          &TaoGlobalStruct::set_print_command
      )
      // TaoGlobalStruct.var_out_file (0D_NOT_character -
      .def_property(
          "var_out_file",
          &TaoGlobalStruct::var_out_file,
          &TaoGlobalStruct::set_var_out_file
      )
      // TaoGlobalStruct.history_file (0D_NOT_character -
      .def_property(
          "history_file",
          &TaoGlobalStruct::history_file,
          &TaoGlobalStruct::set_history_file
      )
      // TaoGlobalStruct.beam_timer_on (0D_NOT_logical - For timing the beam tracking calculation.
      .def_property(
          "beam_timer_on",
          &TaoGlobalStruct::beam_timer_on,
          &TaoGlobalStruct::set_beam_timer_on
      )
      // TaoGlobalStruct.box_plots (0D_NOT_logical - For debugging plot layout issues.
      .def_property("box_plots", &TaoGlobalStruct::box_plots, &TaoGlobalStruct::set_box_plots)
      // TaoGlobalStruct.blank_line_between_commands (0D_NOT_logical - Add a blank line between
      // command output?
      .def_property(
          "blank_line_between_commands",
          &TaoGlobalStruct::blank_line_between_commands,
          &TaoGlobalStruct::set_blank_line_between_commands
      )
      // TaoGlobalStruct.cmd_file_abort_on_error (0D_NOT_logical - Abort open command files if there
      // is an error?
      .def_property(
          "cmd_file_abort_on_error",
          &TaoGlobalStruct::cmd_file_abort_on_error,
          &TaoGlobalStruct::set_cmd_file_abort_on_error
      )
      // TaoGlobalStruct.concatenate_maps (0D_NOT_logical - False => tracking using DA.
      .def_property(
          "concatenate_maps",
          &TaoGlobalStruct::concatenate_maps,
          &TaoGlobalStruct::set_concatenate_maps
      )
      // TaoGlobalStruct.derivative_recalc (0D_NOT_logical - Recalc before each optimizer run?
      .def_property(
          "derivative_recalc",
          &TaoGlobalStruct::derivative_recalc,
          &TaoGlobalStruct::set_derivative_recalc
      )
      // TaoGlobalStruct.derivative_uses_design (0D_NOT_logical - Derivative calc uses design
      // lattice instead of model?
      .def_property(
          "derivative_uses_design",
          &TaoGlobalStruct::derivative_uses_design,
          &TaoGlobalStruct::set_derivative_uses_design
      )
      // TaoGlobalStruct.disable_smooth_line_calc (0D_NOT_logical - Global disable of the smooth
      // line calculation.
      .def_property(
          "disable_smooth_line_calc",
          &TaoGlobalStruct::disable_smooth_line_calc,
          &TaoGlobalStruct::set_disable_smooth_line_calc
      )
      // TaoGlobalStruct.draw_curve_off_scale_warn (0D_NOT_logical - Display warning on graphs?
      .def_property(
          "draw_curve_off_scale_warn",
          &TaoGlobalStruct::draw_curve_off_scale_warn,
          &TaoGlobalStruct::set_draw_curve_off_scale_warn
      )
      // TaoGlobalStruct.external_plotting (0D_NOT_logical - Used with matplotlib and gui.
      .def_property(
          "external_plotting",
          &TaoGlobalStruct::external_plotting,
          &TaoGlobalStruct::set_external_plotting
      )
      // TaoGlobalStruct.label_lattice_elements (0D_NOT_logical - For lat_layout plots
      .def_property(
          "label_lattice_elements",
          &TaoGlobalStruct::label_lattice_elements,
          &TaoGlobalStruct::set_label_lattice_elements
      )
      // TaoGlobalStruct.label_keys (0D_NOT_logical - For lat_layout plots
      .def_property("label_keys", &TaoGlobalStruct::label_keys, &TaoGlobalStruct::set_label_keys)
      // TaoGlobalStruct.lattice_calc_on (0D_NOT_logical - Turn on/off beam and single particle
      // calculations.
      .def_property(
          "lattice_calc_on",
          &TaoGlobalStruct::lattice_calc_on,
          &TaoGlobalStruct::set_lattice_calc_on
      )
      // TaoGlobalStruct.only_limit_opt_vars (0D_NOT_logical - Only apply limits to variables used
      // in optimization.
      .def_property(
          "only_limit_opt_vars",
          &TaoGlobalStruct::only_limit_opt_vars,
          &TaoGlobalStruct::set_only_limit_opt_vars
      )
      // TaoGlobalStruct.opt_with_ref (0D_NOT_logical - Use reference data in optimization?
      .def_property(
          "opt_with_ref",
          &TaoGlobalStruct::opt_with_ref,
          &TaoGlobalStruct::set_opt_with_ref
      )
      // TaoGlobalStruct.opt_with_base (0D_NOT_logical - Use base data in optimization?
      .def_property(
          "opt_with_base",
          &TaoGlobalStruct::opt_with_base,
          &TaoGlobalStruct::set_opt_with_base
      )
      // TaoGlobalStruct.opt_match_auto_recalc (0D_NOT_logical - Set recalc = True for match
      // elements before each cycle?
      .def_property(
          "opt_match_auto_recalc",
          &TaoGlobalStruct::opt_match_auto_recalc,
          &TaoGlobalStruct::set_opt_match_auto_recalc
      )
      // TaoGlobalStruct.opti_write_var_file (0D_NOT_logical - 'run' command writes var_out_file
      .def_property(
          "opti_write_var_file",
          &TaoGlobalStruct::opti_write_var_file,
          &TaoGlobalStruct::set_opti_write_var_file
      )
      // TaoGlobalStruct.optimizer_allow_user_abort (0D_NOT_logical - See Tao manual for more
      // details.
      .def_property(
          "optimizer_allow_user_abort",
          &TaoGlobalStruct::optimizer_allow_user_abort,
          &TaoGlobalStruct::set_optimizer_allow_user_abort
      )
      // TaoGlobalStruct.optimizer_var_limit_warn (0D_NOT_logical - Warn when vars reach a limit
      // with optimization.
      .def_property(
          "optimizer_var_limit_warn",
          &TaoGlobalStruct::optimizer_var_limit_warn,
          &TaoGlobalStruct::set_optimizer_var_limit_warn
      )
      // TaoGlobalStruct.plot_on (0D_NOT_logical - Do plotting?
      .def_property("plot_on", &TaoGlobalStruct::plot_on, &TaoGlobalStruct::set_plot_on)
      // TaoGlobalStruct.rad_int_user_calc_on (0D_NOT_logical - User set radiation integrals
      // calculation on/off.
      .def_property(
          "rad_int_user_calc_on",
          &TaoGlobalStruct::rad_int_user_calc_on,
          &TaoGlobalStruct::set_rad_int_user_calc_on
      )
      // TaoGlobalStruct.rf_on (0D_NOT_logical - RFcavities on or off? Does not affect lcavities.
      .def_property("rf_on", &TaoGlobalStruct::rf_on, &TaoGlobalStruct::set_rf_on)
      // TaoGlobalStruct.single_step (0D_NOT_logical - For debugging and demonstrations: Single step
      // through a command file?
      .def_property("single_step", &TaoGlobalStruct::single_step, &TaoGlobalStruct::set_single_step)
      // TaoGlobalStruct.stop_on_error (0D_NOT_logical - For debugging: False prevents tao from
      // exiting on an error.
      .def_property(
          "stop_on_error",
          &TaoGlobalStruct::stop_on_error,
          &TaoGlobalStruct::set_stop_on_error
      )
      // TaoGlobalStruct.svd_retreat_on_merit_increase (0D_NOT_logical -
      .def_property(
          "svd_retreat_on_merit_increase",
          &TaoGlobalStruct::svd_retreat_on_merit_increase,
          &TaoGlobalStruct::set_svd_retreat_on_merit_increase
      )
      // TaoGlobalStruct.var_limits_on (0D_NOT_logical - Respect the variable limits?
      .def_property(
          "var_limits_on",
          &TaoGlobalStruct::var_limits_on,
          &TaoGlobalStruct::set_var_limits_on
      )
      // TaoGlobalStruct.wait_for_CR_in_single_mode (0D_NOT_logical - For use with a python GUI.
      .def_property(
          "wait_for_CR_in_single_mode",
          &TaoGlobalStruct::wait_for_CR_in_single_mode,
          &TaoGlobalStruct::set_wait_for_CR_in_single_mode
      )
      // TaoGlobalStruct.symbol_import (0D_NOT_logical - Import symbols from lattice file(s)?
      // Internal stuff
      .def_property(
          "symbol_import",
          &TaoGlobalStruct::symbol_import,
          &TaoGlobalStruct::set_symbol_import
      )
      // TaoGlobalStruct.debug_on (0D_NOT_logical - For debugging.
      .def_property("debug_on", &TaoGlobalStruct::debug_on, &TaoGlobalStruct::set_debug_on)
      // TaoGlobalStruct.expression_tree_on (0D_NOT_logical - Use an expression tree instead of a
      // stack?
      .def_property(
          "expression_tree_on",
          &TaoGlobalStruct::expression_tree_on,
          &TaoGlobalStruct::set_expression_tree_on
      )
      // TaoGlobalStruct.verbose_on (0D_NOT_logical - For verbose output. Used with debugging.
      .def_property("verbose_on", &TaoGlobalStruct::verbose_on, &TaoGlobalStruct::set_verbose_on)

      .def("__repr__", [](const TaoGlobalStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoGlobalStruct &self) {
            return TaoGlobalStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoGlobalStruct &self, py::dict &memo) { return TaoGlobalStruct(self); }
      )

      ;

  // 1D TaoGlobalStruct arrays are not used in structs/routines
  // 2D TaoGlobalStruct arrays are not used in structs/routines
  // 3D TaoGlobalStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_graph_struct
void init_tao_graph_struct(py::module &m, py::class_<TaoGraphStruct> &cls) {
  cls.def(py::init<>())
      // TaoGraphStruct.name (0D_NOT_character - Name identifying the graph
      .def_property("name", &TaoGraphStruct::name, &TaoGraphStruct::set_name)
      // TaoGraphStruct.type (0D_NOT_character - 'data', 'lat_layout', 'phase_space', 'histogram',
      // 'dynamic_aperture'
      .def_property("type", &TaoGraphStruct::type, &TaoGraphStruct::set_type)
      // TaoGraphStruct.title (0D_NOT_character -
      .def_property("title", &TaoGraphStruct::title, &TaoGraphStruct::set_title)
      // TaoGraphStruct.title_suffix (0D_NOT_character -
      .def_property(
          "title_suffix",
          &TaoGraphStruct::title_suffix,
          &TaoGraphStruct::set_title_suffix
      )
      // TaoGraphStruct.text_legend (1D_NOT_character - Array for holding descriptive info.
      .def_property_readonly("text_legend", &TaoGraphStruct::text_legend)
      // TaoGraphStruct.text_legend_out (1D_NOT_character - Array for holding descriptive info.
      .def_property_readonly("text_legend_out", &TaoGraphStruct::text_legend_out)
      // TaoGraphStruct.why_invalid (0D_NOT_character - Informative string to print.
      .def_property("why_invalid", &TaoGraphStruct::why_invalid, &TaoGraphStruct::set_why_invalid)
      // TaoGraphStruct.curve (1D_ALLOC_type -
      .def_property_readonly("curve", &TaoGraphStruct::curve)
      // TaoGraphStruct.p (0D_PTR_type - pointer to parent plot
      .def_property("p", &TaoGraphStruct::p, &TaoGraphStruct::set_p)
      // TaoGraphStruct.floor_plan (0D_NOT_type -
      .def_property("floor_plan", &TaoGraphStruct::floor_plan, &TaoGraphStruct::set_floor_plan)
      // TaoGraphStruct.text_legend_origin (0D_NOT_type -
      .def_property(
          "text_legend_origin",
          &TaoGraphStruct::text_legend_origin,
          &TaoGraphStruct::set_text_legend_origin
      )
      // TaoGraphStruct.curve_legend_origin (0D_NOT_type -
      .def_property(
          "curve_legend_origin",
          &TaoGraphStruct::curve_legend_origin,
          &TaoGraphStruct::set_curve_legend_origin
      )
      // TaoGraphStruct.curve_legend (0D_NOT_type -
      .def_property(
          "curve_legend",
          &TaoGraphStruct::curve_legend,
          &TaoGraphStruct::set_curve_legend
      )
      // TaoGraphStruct.x (0D_NOT_type - X-axis parameters.
      .def_property("x", &TaoGraphStruct::x, &TaoGraphStruct::set_x)
      // TaoGraphStruct.y (0D_NOT_type - Y-axis attributes.
      .def_property("y", &TaoGraphStruct::y, &TaoGraphStruct::set_y)
      // TaoGraphStruct.x2 (0D_NOT_type - X2-axis attributes (Not currently used).
      .def_property("x2", &TaoGraphStruct::x2, &TaoGraphStruct::set_x2)
      // TaoGraphStruct.y2 (0D_NOT_type - Y2-axis attributes.
      .def_property("y2", &TaoGraphStruct::y2, &TaoGraphStruct::set_y2)
      // TaoGraphStruct.margin (0D_NOT_type - Margin around the graph.
      .def_property("margin", &TaoGraphStruct::margin, &TaoGraphStruct::set_margin)
      // TaoGraphStruct.scale_margin (0D_NOT_type - Margin for scaling
      .def_property(
          "scale_margin",
          &TaoGraphStruct::scale_margin,
          &TaoGraphStruct::set_scale_margin
      )
      // TaoGraphStruct.x_axis_scale_factor (0D_NOT_real - x-axis conversion from internal to
      // plotting units.
      .def_property(
          "x_axis_scale_factor",
          &TaoGraphStruct::x_axis_scale_factor,
          &TaoGraphStruct::set_x_axis_scale_factor
      )
      // TaoGraphStruct.symbol_size_scale (0D_NOT_real - Symbol size scale factor for phase_space
      // plots.
      .def_property(
          "symbol_size_scale",
          &TaoGraphStruct::symbol_size_scale,
          &TaoGraphStruct::set_symbol_size_scale
      )
      // TaoGraphStruct.box (1D_NOT_integer - Defines which box the plot is put in.
      .def_property("box", &TaoGraphStruct::box, &TaoGraphStruct::set_box)
      // TaoGraphStruct.ix_branch (0D_NOT_integer - Branch in lattice. Used when there are no
      // associated curves.
      .def_property("ix_branch", &TaoGraphStruct::ix_branch, &TaoGraphStruct::set_ix_branch)
      // TaoGraphStruct.ix_universe (0D_NOT_integer - Used for lat_layout plots.
      .def_property("ix_universe", &TaoGraphStruct::ix_universe, &TaoGraphStruct::set_ix_universe)
      // TaoGraphStruct.clip (0D_NOT_logical - Clip plot at graph boundary.
      .def_property("clip", &TaoGraphStruct::clip, &TaoGraphStruct::set_clip)
      // TaoGraphStruct.y2_mirrors_y (0D_NOT_logical - Y2-axis same as Y-axis?
      .def_property(
          "y2_mirrors_y",
          &TaoGraphStruct::y2_mirrors_y,
          &TaoGraphStruct::set_y2_mirrors_y
      )
      // TaoGraphStruct.limited (0D_NOT_logical - True if at least one data point past graph bounds.
      .def_property("limited", &TaoGraphStruct::limited, &TaoGraphStruct::set_limited)
      // TaoGraphStruct.draw_axes (0D_NOT_logical - Draw axes, labels, etc?
      .def_property("draw_axes", &TaoGraphStruct::draw_axes, &TaoGraphStruct::set_draw_axes)
      // TaoGraphStruct.draw_curve_legend (0D_NOT_logical - Legend for displaying curve info.
      .def_property(
          "draw_curve_legend",
          &TaoGraphStruct::draw_curve_legend,
          &TaoGraphStruct::set_draw_curve_legend
      )
      // TaoGraphStruct.draw_grid (0D_NOT_logical - Draw a grid?
      .def_property("draw_grid", &TaoGraphStruct::draw_grid, &TaoGraphStruct::set_draw_grid)
      // TaoGraphStruct.draw_title (0D_NOT_logical -
      .def_property("draw_title", &TaoGraphStruct::draw_title, &TaoGraphStruct::set_draw_title)
      // TaoGraphStruct.draw_only_good_user_data_or_vars (0D_NOT_logical -
      .def_property(
          "draw_only_good_user_data_or_vars",
          &TaoGraphStruct::draw_only_good_user_data_or_vars,
          &TaoGraphStruct::set_draw_only_good_user_data_or_vars
      )
      // TaoGraphStruct.allow_wrap_around (0D_NOT_logical - 'Wrap' curves to extend past lattice
      // boundaries?
      .def_property(
          "allow_wrap_around",
          &TaoGraphStruct::allow_wrap_around,
          &TaoGraphStruct::set_allow_wrap_around
      )
      // TaoGraphStruct.is_valid (0D_NOT_logical - EG: Bad x_axis_type.
      .def_property("is_valid", &TaoGraphStruct::is_valid, &TaoGraphStruct::set_is_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoGraphStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoGraphStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoGraphStruct &self) {
            return TaoGraphStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoGraphStruct &self, py::dict &memo) { return TaoGraphStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoGraphStructArray1D>(m, "TaoGraphStructArray1D");
  bind_FTypeAlloc1D<TaoGraphStructAlloc1D>(m, "TaoGraphStructAlloc1D");
  // 2D TaoGraphStruct arrays are not used in structs/routines
  // 3D TaoGraphStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_histogram_struct
void init_tao_histogram_struct(py::module &m, py::class_<TaoHistogramStruct> &cls) {
  cls.def(py::init<>())
      // TaoHistogramStruct.density_normalized (0D_NOT_logical -
      .def_property(
          "density_normalized",
          &TaoHistogramStruct::density_normalized,
          &TaoHistogramStruct::set_density_normalized
      )
      // TaoHistogramStruct.weight_by_charge (0D_NOT_logical -
      .def_property(
          "weight_by_charge",
          &TaoHistogramStruct::weight_by_charge,
          &TaoHistogramStruct::set_weight_by_charge
      )
      // TaoHistogramStruct.minimum (0D_NOT_real - Computed by Tao. Not User settable.
      .def_property("minimum", &TaoHistogramStruct::minimum, &TaoHistogramStruct::set_minimum)
      // TaoHistogramStruct.maximum (0D_NOT_real - Computed by Tao. Not User settable.
      .def_property("maximum", &TaoHistogramStruct::maximum, &TaoHistogramStruct::set_maximum)
      // TaoHistogramStruct.width (0D_NOT_real -
      .def_property("width", &TaoHistogramStruct::width, &TaoHistogramStruct::set_width)
      // TaoHistogramStruct.center (0D_NOT_real -
      .def_property("center", &TaoHistogramStruct::center, &TaoHistogramStruct::set_center)
      // TaoHistogramStruct.number (0D_NOT_integer -
      .def_property("number", &TaoHistogramStruct::number, &TaoHistogramStruct::set_number)

      .def("__repr__", [](const TaoHistogramStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoHistogramStruct &self) {
            return TaoHistogramStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoHistogramStruct &self, py::dict &memo) { return TaoHistogramStruct(self); }
      )

      ;

  // 1D TaoHistogramStruct arrays are not used in structs/routines
  // 2D TaoHistogramStruct arrays are not used in structs/routines
  // 3D TaoHistogramStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_init_struct
void init_tao_init_struct(py::module &m, py::class_<TaoInitStruct> &cls) {
  cls.def(py::init<>())
      // TaoInitStruct.parse_cmd_args (0D_NOT_logical - Used by custom programs to control Tao init
      .def_property(
          "parse_cmd_args",
          &TaoInitStruct::parse_cmd_args,
          &TaoInitStruct::set_parse_cmd_args
      )
      // TaoInitStruct.debug_switch (0D_NOT_logical - Is the '-debug' switch present?
      .def_property("debug_switch", &TaoInitStruct::debug_switch, &TaoInitStruct::set_debug_switch)
      // TaoInitStruct.external_plotting_switch (0D_NOT_logical - Is '-external_plotting' switch
      // present?
      .def_property(
          "external_plotting_switch",
          &TaoInitStruct::external_plotting_switch,
          &TaoInitStruct::set_external_plotting_switch
      )
      // TaoInitStruct.init_name (0D_NOT_character - label for initialization
      .def_property("init_name", &TaoInitStruct::init_name, &TaoInitStruct::set_init_name)
      // TaoInitStruct.hook_init_file (0D_NOT_character -
      .def_property(
          "hook_init_file",
          &TaoInitStruct::hook_init_file,
          &TaoInitStruct::set_hook_init_file
      )
      // TaoInitStruct.hook_lat_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_lat_file",
          &TaoInitStruct::hook_lat_file,
          &TaoInitStruct::set_hook_lat_file
      )
      // TaoInitStruct.hook_beam_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_beam_file",
          &TaoInitStruct::hook_beam_file,
          &TaoInitStruct::set_hook_beam_file
      )
      // TaoInitStruct.hook_data_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_data_file",
          &TaoInitStruct::hook_data_file,
          &TaoInitStruct::set_hook_data_file
      )
      // TaoInitStruct.hook_plot_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_plot_file",
          &TaoInitStruct::hook_plot_file,
          &TaoInitStruct::set_hook_plot_file
      )
      // TaoInitStruct.hook_startup_file (0D_NOT_character - To be set by
      // tao_hook_parse_command_args
      .def_property(
          "hook_startup_file",
          &TaoInitStruct::hook_startup_file,
          &TaoInitStruct::set_hook_startup_file
      )
      // TaoInitStruct.hook_var_file (0D_NOT_character - To be set by tao_hook_parse_command_args
      .def_property(
          "hook_var_file",
          &TaoInitStruct::hook_var_file,
          &TaoInitStruct::set_hook_var_file
      )
      // TaoInitStruct.hook_building_wall_file (0D_NOT_character - To be set by
      // tao_hook_parse_command_args
      .def_property(
          "hook_building_wall_file",
          &TaoInitStruct::hook_building_wall_file,
          &TaoInitStruct::set_hook_building_wall_file
      )
      // TaoInitStruct.init_file_arg_path (0D_NOT_character - Path part of init_tao_file
      .def_property(
          "init_file_arg_path",
          &TaoInitStruct::init_file_arg_path,
          &TaoInitStruct::set_init_file_arg_path
      )
      // TaoInitStruct.lattice_file_arg (0D_NOT_character - -lattice_file        command line
      // argument.
      .def_property(
          "lattice_file_arg",
          &TaoInitStruct::lattice_file_arg,
          &TaoInitStruct::set_lattice_file_arg
      )
      // TaoInitStruct.hook_init_file_arg (0D_NOT_character - -hook_init_file      command line
      // argument
      .def_property(
          "hook_init_file_arg",
          &TaoInitStruct::hook_init_file_arg,
          &TaoInitStruct::set_hook_init_file_arg
      )
      // TaoInitStruct.init_file_arg (0D_NOT_character - -init_file           command line argument.
      .def_property(
          "init_file_arg",
          &TaoInitStruct::init_file_arg,
          &TaoInitStruct::set_init_file_arg
      )
      // TaoInitStruct.beam_file_arg (0D_NOT_character - -beam_file           command line argument.
      .def_property(
          "beam_file_arg",
          &TaoInitStruct::beam_file_arg,
          &TaoInitStruct::set_beam_file_arg
      )
      // TaoInitStruct.beam_init_position_file_arg (0D_NOT_character - -beam_init_position_file
      // command line argument.
      .def_property(
          "beam_init_position_file_arg",
          &TaoInitStruct::beam_init_position_file_arg,
          &TaoInitStruct::set_beam_init_position_file_arg
      )
      // TaoInitStruct.command_arg (0D_NOT_character - -command             command line argument.
      .def_property("command_arg", &TaoInitStruct::command_arg, &TaoInitStruct::set_command_arg)
      // TaoInitStruct.data_file_arg (0D_NOT_character - -data_file           command line argument.
      .def_property(
          "data_file_arg",
          &TaoInitStruct::data_file_arg,
          &TaoInitStruct::set_data_file_arg
      )
      // TaoInitStruct.plot_file_arg (0D_NOT_character - -plot_file           command line argument.
      .def_property(
          "plot_file_arg",
          &TaoInitStruct::plot_file_arg,
          &TaoInitStruct::set_plot_file_arg
      )
      // TaoInitStruct.startup_file_arg (0D_NOT_character - -startup_file        command line
      // argument.
      .def_property(
          "startup_file_arg",
          &TaoInitStruct::startup_file_arg,
          &TaoInitStruct::set_startup_file_arg
      )
      // TaoInitStruct.var_file_arg (0D_NOT_character - -var_file            command line argument.
      .def_property("var_file_arg", &TaoInitStruct::var_file_arg, &TaoInitStruct::set_var_file_arg)
      // TaoInitStruct.building_wall_file_arg (0D_NOT_character - -building_wall_file  command line
      // argument.
      .def_property(
          "building_wall_file_arg",
          &TaoInitStruct::building_wall_file_arg,
          &TaoInitStruct::set_building_wall_file_arg
      )
      // TaoInitStruct.geometry_arg (0D_NOT_character - -geometry            command line argument.
      .def_property("geometry_arg", &TaoInitStruct::geometry_arg, &TaoInitStruct::set_geometry_arg)
      // TaoInitStruct.slice_lattice_arg (0D_NOT_character - -slice_lattice       command line
      // argument.
      .def_property(
          "slice_lattice_arg",
          &TaoInitStruct::slice_lattice_arg,
          &TaoInitStruct::set_slice_lattice_arg
      )
      // TaoInitStruct.start_branch_at_arg (0D_NOT_character - -start_branch_at     command line
      // argument.
      .def_property(
          "start_branch_at_arg",
          &TaoInitStruct::start_branch_at_arg,
          &TaoInitStruct::set_start_branch_at_arg
      )
      // TaoInitStruct.log_startup_arg (0D_NOT_character - -log_startup         command line
      // argument
      .def_property(
          "log_startup_arg",
          &TaoInitStruct::log_startup_arg,
          &TaoInitStruct::set_log_startup_arg
      )
      // TaoInitStruct.no_stopping_arg (0D_NOT_character - -no_stopping         command line
      // argument
      .def_property(
          "no_stopping_arg",
          &TaoInitStruct::no_stopping_arg,
          &TaoInitStruct::set_no_stopping_arg
      )
      // TaoInitStruct.noplot_arg (0D_NOT_character - -noplot              command line argument
      .def_property("noplot_arg", &TaoInitStruct::noplot_arg, &TaoInitStruct::set_noplot_arg)
      // TaoInitStruct.no_rad_int_arg (0D_NOT_character - -no_rad_int          command line argument
      .def_property(
          "no_rad_int_arg",
          &TaoInitStruct::no_rad_int_arg,
          &TaoInitStruct::set_no_rad_int_arg
      )
      // TaoInitStruct.reverse_arg (0D_NOT_character - -reverse             command line argument
      .def_property("reverse_arg", &TaoInitStruct::reverse_arg, &TaoInitStruct::set_reverse_arg)
      // TaoInitStruct.debug_arg (0D_NOT_character - -debug               command line argument
      .def_property("debug_arg", &TaoInitStruct::debug_arg, &TaoInitStruct::set_debug_arg)
      // TaoInitStruct.disable_smooth_line_calc_arg (0D_NOT_character - -disable_smooth_line_calc
      .def_property(
          "disable_smooth_line_calc_arg",
          &TaoInitStruct::disable_smooth_line_calc_arg,
          &TaoInitStruct::set_disable_smooth_line_calc_arg
      )
      // TaoInitStruct.rf_on_arg (0D_NOT_character - -rf_on               command line argument
      .def_property("rf_on_arg", &TaoInitStruct::rf_on_arg, &TaoInitStruct::set_rf_on_arg)
      // TaoInitStruct.prompt_color_arg (0D_NOT_character - -prompt_color        command line
      // argument
      .def_property(
          "prompt_color_arg",
          &TaoInitStruct::prompt_color_arg,
          &TaoInitStruct::set_prompt_color_arg
      )
      // TaoInitStruct.quiet_arg (0D_NOT_character - -quiet               command line argument
      .def_property("quiet_arg", &TaoInitStruct::quiet_arg, &TaoInitStruct::set_quiet_arg)
      // TaoInitStruct.noinit_arg (0D_NOT_character - -noinit              command line argument
      .def_property("noinit_arg", &TaoInitStruct::noinit_arg, &TaoInitStruct::set_noinit_arg)
      // TaoInitStruct.nostartup_arg (0D_NOT_character - -nostartup           command line argument
      .def_property(
          "nostartup_arg",
          &TaoInitStruct::nostartup_arg,
          &TaoInitStruct::set_nostartup_arg
      )
      // TaoInitStruct.symbol_import_arg (0D_NOT_character - -symbol_import       command line
      // argument
      .def_property(
          "symbol_import_arg",
          &TaoInitStruct::symbol_import_arg,
          &TaoInitStruct::set_symbol_import_arg
      )
      // TaoInitStruct.unique_name_suffix (0D_NOT_character -
      .def_property(
          "unique_name_suffix",
          &TaoInitStruct::unique_name_suffix,
          &TaoInitStruct::set_unique_name_suffix
      )

      .def("__repr__", [](const TaoInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoInitStruct &self) {
            return TaoInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoInitStruct &self, py::dict &memo) { return TaoInitStruct(self); }
      )

      ;

  // 1D TaoInitStruct arrays are not used in structs/routines
  // 2D TaoInitStruct arrays are not used in structs/routines
  // 3D TaoInitStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lat_sigma_struct
void init_tao_lat_sigma_struct(py::module &m, py::class_<TaoLatSigmaStruct> &cls) {
  cls.def(py::init<>())
      // TaoLatSigmaStruct.mat (2D_NOT_real -
      .def_property("mat", &TaoLatSigmaStruct::mat, &TaoLatSigmaStruct::set_mat)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoLatSigmaStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoLatSigmaStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatSigmaStruct &self) {
            return TaoLatSigmaStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoLatSigmaStruct &self, py::dict &memo) { return TaoLatSigmaStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoLatSigmaStructArray1D>(m, "TaoLatSigmaStructArray1D");
  bind_FTypeAlloc1D<TaoLatSigmaStructAlloc1D>(m, "TaoLatSigmaStructAlloc1D");
  // 2D TaoLatSigmaStruct arrays are not used in structs/routines
  // 3D TaoLatSigmaStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_branch_struct
void init_tao_lattice_branch_struct(py::module &m, py::class_<TaoLatticeBranchStruct> &cls) {
  cls.def(py::init<>())
      // TaoLatticeBranchStruct.tao_lat (0D_PTR_type - Parent tao_lat
      .def_property(
          "tao_lat",
          &TaoLatticeBranchStruct::tao_lat,
          &TaoLatticeBranchStruct::set_tao_lat
      )
      // TaoLatticeBranchStruct.lat_sigma (1D_ALLOC_type - Sigma matrix derived from lattice (not
      // beam).
      .def_property_readonly("lat_sigma", &TaoLatticeBranchStruct::lat_sigma)
      // TaoLatticeBranchStruct.spin_ele (1D_ALLOC_type - Spin stuff
      .def_property_readonly("spin_ele", &TaoLatticeBranchStruct::spin_ele)
      // TaoLatticeBranchStruct.bunch_params (1D_ALLOC_type - Per element
      .def_property_readonly("bunch_params", &TaoLatticeBranchStruct::bunch_params)
      // TaoLatticeBranchStruct.bunch_params_comb (1D_ALLOC_type - A comb for each bunch in beam.
      .def_property_readonly("bunch_params_comb", &TaoLatticeBranchStruct::bunch_params_comb)
      // TaoLatticeBranchStruct.orbit (1D_ALLOC_type -
      .def_property_readonly("orbit", &TaoLatticeBranchStruct::orbit)
      // TaoLatticeBranchStruct.plot_cache (1D_ALLOC_type - Plotting data cache
      .def_property_readonly("plot_cache", &TaoLatticeBranchStruct::plot_cache)
      // TaoLatticeBranchStruct.spin (0D_NOT_type -
      .def_property("spin", &TaoLatticeBranchStruct::spin, &TaoLatticeBranchStruct::set_spin)
      // TaoLatticeBranchStruct.srdt (0D_NOT_type -
      .def_property("srdt", &TaoLatticeBranchStruct::srdt, &TaoLatticeBranchStruct::set_srdt)
      // TaoLatticeBranchStruct.orb0 (0D_NOT_type - For saving beginning orbit
      .def_property("orb0", &TaoLatticeBranchStruct::orb0, &TaoLatticeBranchStruct::set_orb0)
      // TaoLatticeBranchStruct.modes_ri (0D_NOT_type - Synchrotron integrals stuff
      .def_property(
          "modes_ri",
          &TaoLatticeBranchStruct::modes_ri,
          &TaoLatticeBranchStruct::set_modes_ri
      )
      // TaoLatticeBranchStruct.modes_6d (0D_NOT_type - 6D radiation matrices.
      .def_property(
          "modes_6d",
          &TaoLatticeBranchStruct::modes_6d,
          &TaoLatticeBranchStruct::set_modes_6d
      )
      // TaoLatticeBranchStruct.ptc_normal_form (0D_NOT_type - Collection of normal form structures
      // defined in PTC
      .def_property(
          "ptc_normal_form",
          &TaoLatticeBranchStruct::ptc_normal_form,
          &TaoLatticeBranchStruct::set_ptc_normal_form
      )
      // TaoLatticeBranchStruct.bmad_normal_form (0D_NOT_type - Collection of normal form structures
      // defined in Bmad
      .def_property(
          "bmad_normal_form",
          &TaoLatticeBranchStruct::bmad_normal_form,
          &TaoLatticeBranchStruct::set_bmad_normal_form
      )
      // TaoLatticeBranchStruct.high_E_orb (1D_ALLOC_type -
      .def_property_readonly("high_E_orb", &TaoLatticeBranchStruct::high_E_orb)
      // TaoLatticeBranchStruct.low_E_orb (1D_ALLOC_type -
      .def_property_readonly("low_E_orb", &TaoLatticeBranchStruct::low_E_orb)
      // TaoLatticeBranchStruct.taylor_save (1D_NOT_type - Save to reduce computation time.
      .def_property_readonly("taylor_save", &TaoLatticeBranchStruct::taylor_save)
      // TaoLatticeBranchStruct.cache_x_min (0D_NOT_real -
      .def_property(
          "cache_x_min",
          &TaoLatticeBranchStruct::cache_x_min,
          &TaoLatticeBranchStruct::set_cache_x_min
      )
      // TaoLatticeBranchStruct.cache_x_max (0D_NOT_real -
      .def_property(
          "cache_x_max",
          &TaoLatticeBranchStruct::cache_x_max,
          &TaoLatticeBranchStruct::set_cache_x_max
      )
      // TaoLatticeBranchStruct.comb_ds_save (0D_NOT_real - Master parameter for
      // %bunch_params_comb(:)%ds_save
      .def_property(
          "comb_ds_save",
          &TaoLatticeBranchStruct::comb_ds_save,
          &TaoLatticeBranchStruct::set_comb_ds_save
      )
      // TaoLatticeBranchStruct.ix_ref_taylor (0D_NOT_integer -
      .def_property(
          "ix_ref_taylor",
          &TaoLatticeBranchStruct::ix_ref_taylor,
          &TaoLatticeBranchStruct::set_ix_ref_taylor
      )
      // TaoLatticeBranchStruct.ix_ele_taylor (0D_NOT_integer -
      .def_property(
          "ix_ele_taylor",
          &TaoLatticeBranchStruct::ix_ele_taylor,
          &TaoLatticeBranchStruct::set_ix_ele_taylor
      )
      // TaoLatticeBranchStruct.track_state (0D_NOT_integer -
      .def_property(
          "track_state",
          &TaoLatticeBranchStruct::track_state,
          &TaoLatticeBranchStruct::set_track_state
      )
      // TaoLatticeBranchStruct.cache_n_pts (0D_NOT_integer -
      .def_property(
          "cache_n_pts",
          &TaoLatticeBranchStruct::cache_n_pts,
          &TaoLatticeBranchStruct::set_cache_n_pts
      )
      // TaoLatticeBranchStruct.ix_rad_int_cache (0D_NOT_integer - Radiation integrals cache index.
      .def_property(
          "ix_rad_int_cache",
          &TaoLatticeBranchStruct::ix_rad_int_cache,
          &TaoLatticeBranchStruct::set_ix_rad_int_cache
      )
      // TaoLatticeBranchStruct.has_open_match_element (0D_NOT_logical -
      .def_property(
          "has_open_match_element",
          &TaoLatticeBranchStruct::has_open_match_element,
          &TaoLatticeBranchStruct::set_has_open_match_element
      )
      // TaoLatticeBranchStruct.plot_cache_valid (0D_NOT_logical - Valid plotting data cache?
      .def_property(
          "plot_cache_valid",
          &TaoLatticeBranchStruct::plot_cache_valid,
          &TaoLatticeBranchStruct::set_plot_cache_valid
      )
      // TaoLatticeBranchStruct.spin_map_valid (0D_NOT_logical -
      .def_property(
          "spin_map_valid",
          &TaoLatticeBranchStruct::spin_map_valid,
          &TaoLatticeBranchStruct::set_spin_map_valid
      )
      // TaoLatticeBranchStruct.twiss_valid (0D_NOT_logical - Invalid EG with unstable 1-turn matrix
      // with a closed branch. With open branch: twiss_valid = T even if some Twiss (and orbit) is
      // invalid.
      .def_property(
          "twiss_valid",
          &TaoLatticeBranchStruct::twiss_valid,
          &TaoLatticeBranchStruct::set_twiss_valid
      )
      // TaoLatticeBranchStruct.mode_flip_here (0D_NOT_logical - Twiss parameter mode flip seen?
      .def_property(
          "mode_flip_here",
          &TaoLatticeBranchStruct::mode_flip_here,
          &TaoLatticeBranchStruct::set_mode_flip_here
      )
      // TaoLatticeBranchStruct.chrom_calc_ok (0D_NOT_logical -
      .def_property(
          "chrom_calc_ok",
          &TaoLatticeBranchStruct::chrom_calc_ok,
          &TaoLatticeBranchStruct::set_chrom_calc_ok
      )
      // TaoLatticeBranchStruct.rad_int_calc_ok (0D_NOT_logical -
      .def_property(
          "rad_int_calc_ok",
          &TaoLatticeBranchStruct::rad_int_calc_ok,
          &TaoLatticeBranchStruct::set_rad_int_calc_ok
      )
      // TaoLatticeBranchStruct.emit_6d_calc_ok (0D_NOT_logical -
      .def_property(
          "emit_6d_calc_ok",
          &TaoLatticeBranchStruct::emit_6d_calc_ok,
          &TaoLatticeBranchStruct::set_emit_6d_calc_ok
      )
      // TaoLatticeBranchStruct.sigma_track_ok (0D_NOT_logical -
      .def_property(
          "sigma_track_ok",
          &TaoLatticeBranchStruct::sigma_track_ok,
          &TaoLatticeBranchStruct::set_sigma_track_ok
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoLatticeBranchStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoLatticeBranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatticeBranchStruct &self) {
            return TaoLatticeBranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoLatticeBranchStruct &self, py::dict &memo) {
            return TaoLatticeBranchStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoLatticeBranchStructArray1D>(m, "TaoLatticeBranchStructArray1D");
  bind_FTypeAlloc1D<TaoLatticeBranchStructAlloc1D>(m, "TaoLatticeBranchStructAlloc1D");
  // 2D TaoLatticeBranchStruct arrays are not used in structs/routines
  // 3D TaoLatticeBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_struct
void init_tao_lattice_struct(py::module &m, py::class_<TaoLatticeStruct> &cls) {
  cls.def(py::init<>())
      // TaoLatticeStruct.name (0D_NOT_character - 'model', 'base', or 'design'.
      .def_property("name", &TaoLatticeStruct::name, &TaoLatticeStruct::set_name)
      // TaoLatticeStruct.lat (0D_NOT_type - lattice structures
      .def_property("lat", &TaoLatticeStruct::lat, &TaoLatticeStruct::set_lat)
      // TaoLatticeStruct.high_E_lat (0D_NOT_type - For chrom calc.
      .def_property("high_E_lat", &TaoLatticeStruct::high_E_lat, &TaoLatticeStruct::set_high_E_lat)
      // TaoLatticeStruct.low_E_lat (0D_NOT_type - For chrom calc.
      .def_property("low_E_lat", &TaoLatticeStruct::low_E_lat, &TaoLatticeStruct::set_low_E_lat)
      // TaoLatticeStruct.rad_int_by_ele_ri (0D_NOT_type -
      .def_property(
          "rad_int_by_ele_ri",
          &TaoLatticeStruct::rad_int_by_ele_ri,
          &TaoLatticeStruct::set_rad_int_by_ele_ri
      )
      // TaoLatticeStruct.rad_int_by_ele_6d (0D_NOT_type -
      .def_property(
          "rad_int_by_ele_6d",
          &TaoLatticeStruct::rad_int_by_ele_6d,
          &TaoLatticeStruct::set_rad_int_by_ele_6d
      )
      // TaoLatticeStruct.tao_branch (1D_ALLOC_type -
      .def_property_readonly("tao_branch", &TaoLatticeStruct::tao_branch)

      .def("__repr__", [](const TaoLatticeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatticeStruct &self) {
            return TaoLatticeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoLatticeStruct &self, py::dict &memo) { return TaoLatticeStruct(self); }
      )

      ;

  // 1D TaoLatticeStruct arrays are not used in structs/routines
  // 2D TaoLatticeStruct arrays are not used in structs/routines
  // 3D TaoLatticeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_model_branch_struct
void init_tao_model_branch_struct(py::module &m, py::class_<TaoModelBranchStruct> &cls) {
  cls.def(py::init<>())
      // TaoModelBranchStruct.ele (1D_ALLOC_type - Per element information
      .def_property_readonly("ele", &TaoModelBranchStruct::ele)
      // TaoModelBranchStruct.beam (0D_NOT_type -
      .def_property("beam", &TaoModelBranchStruct::beam, &TaoModelBranchStruct::set_beam)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoModelBranchStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoModelBranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoModelBranchStruct &self) {
            return TaoModelBranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoModelBranchStruct &self, py::dict &memo) {
            return TaoModelBranchStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoModelBranchStructArray1D>(m, "TaoModelBranchStructArray1D");
  bind_FTypeAlloc1D<TaoModelBranchStructAlloc1D>(m, "TaoModelBranchStructAlloc1D");
  // 2D TaoModelBranchStruct arrays are not used in structs/routines
  // 3D TaoModelBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_model_element_struct
void init_tao_model_element_struct(py::module &m, py::class_<TaoModelElementStruct> &cls) {
  cls.def(py::init<>())
      // TaoModelElementStruct.beam (0D_NOT_type - Beam distribution at element.
      .def_property("beam", &TaoModelElementStruct::beam, &TaoModelElementStruct::set_beam)
      // TaoModelElementStruct.save_beam_internally (0D_NOT_logical - Save beam here? Beam also
      // saved at fork elements and at track ends.
      .def_property(
          "save_beam_internally",
          &TaoModelElementStruct::save_beam_internally,
          &TaoModelElementStruct::set_save_beam_internally
      )
      // TaoModelElementStruct.save_beam_to_file (0D_NOT_logical - Save beam to a file? Beam also
      // saved at fork elements and at track ends.
      .def_property(
          "save_beam_to_file",
          &TaoModelElementStruct::save_beam_to_file,
          &TaoModelElementStruct::set_save_beam_to_file
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoModelElementStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoModelElementStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoModelElementStruct &self) {
            return TaoModelElementStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoModelElementStruct &self, py::dict &memo) {
            return TaoModelElementStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoModelElementStructArray1D>(m, "TaoModelElementStructArray1D");
  bind_FTypeAlloc1D<TaoModelElementStructAlloc1D>(m, "TaoModelElementStructAlloc1D");
  // 2D TaoModelElementStruct arrays are not used in structs/routines
  // 3D TaoModelElementStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ping_scale_struct
void init_tao_ping_scale_struct(py::module &m, py::class_<TaoPingScaleStruct> &cls) {
  cls.def(py::init<>())
      // TaoPingScaleStruct.a_mode_meas (0D_NOT_real -
      .def_property(
          "a_mode_meas",
          &TaoPingScaleStruct::a_mode_meas,
          &TaoPingScaleStruct::set_a_mode_meas
      )
      // TaoPingScaleStruct.a_mode_ref (0D_NOT_real -
      .def_property(
          "a_mode_ref",
          &TaoPingScaleStruct::a_mode_ref,
          &TaoPingScaleStruct::set_a_mode_ref
      )
      // TaoPingScaleStruct.b_mode_meas (0D_NOT_real -
      .def_property(
          "b_mode_meas",
          &TaoPingScaleStruct::b_mode_meas,
          &TaoPingScaleStruct::set_b_mode_meas
      )
      // TaoPingScaleStruct.b_mode_ref (0D_NOT_real -
      .def_property(
          "b_mode_ref",
          &TaoPingScaleStruct::b_mode_ref,
          &TaoPingScaleStruct::set_b_mode_ref
      )

      .def("__repr__", [](const TaoPingScaleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPingScaleStruct &self) {
            return TaoPingScaleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPingScaleStruct &self, py::dict &memo) { return TaoPingScaleStruct(self); }
      )

      ;

  // 1D TaoPingScaleStruct arrays are not used in structs/routines
  // 2D TaoPingScaleStruct arrays are not used in structs/routines
  // 3D TaoPingScaleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_cache_struct
void init_tao_plot_cache_struct(py::module &m, py::class_<TaoPlotCacheStruct> &cls) {
  cls.def(py::init<>())
      // TaoPlotCacheStruct.ele_to_s (0D_NOT_type - Integrated element from branch beginning. Will
      // be marked as a hybrid element.
      .def_property("ele_to_s", &TaoPlotCacheStruct::ele_to_s, &TaoPlotCacheStruct::set_ele_to_s)
      // TaoPlotCacheStruct.orbit (0D_NOT_type -
      .def_property("orbit", &TaoPlotCacheStruct::orbit, &TaoPlotCacheStruct::set_orbit)
      // TaoPlotCacheStruct.err (0D_NOT_logical -
      .def_property("err", &TaoPlotCacheStruct::err, &TaoPlotCacheStruct::set_err)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoPlotCacheStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoPlotCacheStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotCacheStruct &self) {
            return TaoPlotCacheStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotCacheStruct &self, py::dict &memo) { return TaoPlotCacheStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoPlotCacheStructArray1D>(m, "TaoPlotCacheStructArray1D");
  bind_FTypeAlloc1D<TaoPlotCacheStructAlloc1D>(m, "TaoPlotCacheStructAlloc1D");
  // 2D TaoPlotCacheStruct arrays are not used in structs/routines
  // 3D TaoPlotCacheStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_page_struct
void init_tao_plot_page_struct(py::module &m, py::class_<TaoPlotPageStruct> &cls) {
  cls.def(py::init<>())
      // TaoPlotPageStruct.title (0D_NOT_type - Title  at top of page.
      .def_property("title", &TaoPlotPageStruct::title, &TaoPlotPageStruct::set_title)
      // TaoPlotPageStruct.subtitle (0D_NOT_type - Subtitle below title at top of page.
      .def_property("subtitle", &TaoPlotPageStruct::subtitle, &TaoPlotPageStruct::set_subtitle)
      // TaoPlotPageStruct.border (0D_NOT_type - Border around plots edge of page.
      .def_property("border", &TaoPlotPageStruct::border, &TaoPlotPageStruct::set_border)
      // TaoPlotPageStruct.floor_plan (0D_NOT_type -
      .def_property(
          "floor_plan",
          &TaoPlotPageStruct::floor_plan,
          &TaoPlotPageStruct::set_floor_plan
      )
      // TaoPlotPageStruct.lat_layout (0D_NOT_type -
      .def_property(
          "lat_layout",
          &TaoPlotPageStruct::lat_layout,
          &TaoPlotPageStruct::set_lat_layout
      )
      // TaoPlotPageStruct.pattern (1D_ALLOC_type -
      .def_property_readonly("pattern", &TaoPlotPageStruct::pattern)
      // TaoPlotPageStruct.template_ (1D_ALLOC_type - Templates for the plots.
      .def_property_readonly("template_", &TaoPlotPageStruct::template_)
      // TaoPlotPageStruct.region (1D_ALLOC_type -
      .def_property_readonly("region", &TaoPlotPageStruct::region)
      // TaoPlotPageStruct.plot_display_type (0D_NOT_character - 'X' or 'TK'
      .def_property(
          "plot_display_type",
          &TaoPlotPageStruct::plot_display_type,
          &TaoPlotPageStruct::set_plot_display_type
      )
      // TaoPlotPageStruct.size (1D_NOT_real - width and height of plot window in pixels.
      .def_property("size", &TaoPlotPageStruct::size, &TaoPlotPageStruct::set_size)
      // TaoPlotPageStruct.text_height (0D_NOT_real - In points. Scales the height of all text
      .def_property(
          "text_height",
          &TaoPlotPageStruct::text_height,
          &TaoPlotPageStruct::set_text_height
      )
      // TaoPlotPageStruct.main_title_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "main_title_text_scale",
          &TaoPlotPageStruct::main_title_text_scale,
          &TaoPlotPageStruct::set_main_title_text_scale
      )
      // TaoPlotPageStruct.graph_title_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "graph_title_text_scale",
          &TaoPlotPageStruct::graph_title_text_scale,
          &TaoPlotPageStruct::set_graph_title_text_scale
      )
      // TaoPlotPageStruct.axis_number_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "axis_number_text_scale",
          &TaoPlotPageStruct::axis_number_text_scale,
          &TaoPlotPageStruct::set_axis_number_text_scale
      )
      // TaoPlotPageStruct.axis_label_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "axis_label_text_scale",
          &TaoPlotPageStruct::axis_label_text_scale,
          &TaoPlotPageStruct::set_axis_label_text_scale
      )
      // TaoPlotPageStruct.legend_text_scale (0D_NOT_real - Relative to text_height. For legends,
      // plot_page, and lat_layout
      .def_property(
          "legend_text_scale",
          &TaoPlotPageStruct::legend_text_scale,
          &TaoPlotPageStruct::set_legend_text_scale
      )
      // TaoPlotPageStruct.key_table_text_scale (0D_NOT_real - Relative to text_height
      .def_property(
          "key_table_text_scale",
          &TaoPlotPageStruct::key_table_text_scale,
          &TaoPlotPageStruct::set_key_table_text_scale
      )
      // TaoPlotPageStruct.floor_plan_shape_scale (0D_NOT_real -
      .def_property(
          "floor_plan_shape_scale",
          &TaoPlotPageStruct::floor_plan_shape_scale,
          &TaoPlotPageStruct::set_floor_plan_shape_scale
      )
      // TaoPlotPageStruct.floor_plan_text_scale (0D_NOT_real - Scale used = floor_plan_text_scale *
      // legend_text_scale
      .def_property(
          "floor_plan_text_scale",
          &TaoPlotPageStruct::floor_plan_text_scale,
          &TaoPlotPageStruct::set_floor_plan_text_scale
      )
      // TaoPlotPageStruct.lat_layout_shape_scale (0D_NOT_real -
      .def_property(
          "lat_layout_shape_scale",
          &TaoPlotPageStruct::lat_layout_shape_scale,
          &TaoPlotPageStruct::set_lat_layout_shape_scale
      )
      // TaoPlotPageStruct.lat_layout_text_scale (0D_NOT_real - Scale used = lat_layout_text_scale *
      // legend_text_scale
      .def_property(
          "lat_layout_text_scale",
          &TaoPlotPageStruct::lat_layout_text_scale,
          &TaoPlotPageStruct::set_lat_layout_text_scale
      )
      // TaoPlotPageStruct.n_curve_pts (0D_NOT_integer - Default number of points for plotting a
      // smooth curve.
      .def_property(
          "n_curve_pts",
          &TaoPlotPageStruct::n_curve_pts,
          &TaoPlotPageStruct::set_n_curve_pts
      )
      // TaoPlotPageStruct.id_window (0D_NOT_integer - X window id number.
      .def_property("id_window", &TaoPlotPageStruct::id_window, &TaoPlotPageStruct::set_id_window)
      // TaoPlotPageStruct.delete_overlapping_plots (0D_NOT_logical - Delete overlapping plots when
      // a plot is placed?
      .def_property(
          "delete_overlapping_plots",
          &TaoPlotPageStruct::delete_overlapping_plots,
          &TaoPlotPageStruct::set_delete_overlapping_plots
      )
      // TaoPlotPageStruct.draw_graph_title_suffix (0D_NOT_logical - Draw the graph title suffix?
      .def_property(
          "draw_graph_title_suffix",
          &TaoPlotPageStruct::draw_graph_title_suffix,
          &TaoPlotPageStruct::set_draw_graph_title_suffix
      )

      .def("__repr__", [](const TaoPlotPageStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotPageStruct &self) {
            return TaoPlotPageStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotPageStruct &self, py::dict &memo) { return TaoPlotPageStruct(self); }
      )

      ;

  // 1D TaoPlotPageStruct arrays are not used in structs/routines
  // 2D TaoPlotPageStruct arrays are not used in structs/routines
  // 3D TaoPlotPageStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_region_struct
void init_tao_plot_region_struct(py::module &m, py::class_<TaoPlotRegionStruct> &cls) {
  cls.def(py::init<>())
      // TaoPlotRegionStruct.name (0D_NOT_character - Region name. Eg: 'r13', etc.
      .def_property("name", &TaoPlotRegionStruct::name, &TaoPlotRegionStruct::set_name)
      // TaoPlotRegionStruct.plot (0D_NOT_type - Plot associated with this region
      .def_property("plot", &TaoPlotRegionStruct::plot, &TaoPlotRegionStruct::set_plot)
      // TaoPlotRegionStruct.location (1D_NOT_real - [x1, x2, y1, y2] location on page.
      .def_property("location", &TaoPlotRegionStruct::location, &TaoPlotRegionStruct::set_location)
      // TaoPlotRegionStruct.visible (0D_NOT_logical - To draw or not to draw.
      .def_property("visible", &TaoPlotRegionStruct::visible, &TaoPlotRegionStruct::set_visible)
      // TaoPlotRegionStruct.list_with_show_plot_command (0D_NOT_logical - False used for default
      // plots to shorten the output of 'show plot'
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotRegionStruct::list_with_show_plot_command,
          &TaoPlotRegionStruct::set_list_with_show_plot_command
      )
      // TaoPlotRegionStruct.setup_done (0D_NOT_logical - Used for plot bookkeeping.
      .def_property(
          "setup_done",
          &TaoPlotRegionStruct::setup_done,
          &TaoPlotRegionStruct::set_setup_done
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoPlotRegionStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoPlotRegionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotRegionStruct &self) {
            return TaoPlotRegionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotRegionStruct &self, py::dict &memo) { return TaoPlotRegionStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoPlotRegionStructArray1D>(m, "TaoPlotRegionStructArray1D");
  bind_FTypeAlloc1D<TaoPlotRegionStructAlloc1D>(m, "TaoPlotRegionStructAlloc1D");
  // 2D TaoPlotRegionStruct arrays are not used in structs/routines
  // 3D TaoPlotRegionStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_struct
void init_tao_plot_struct(py::module &m, py::class_<TaoPlotStruct> &cls) {
  cls.def(py::init<>())
      // TaoPlotStruct.name (0D_NOT_character - Identifying name. Rule: If name is blank, plot is
      // not valid.
      .def_property("name", &TaoPlotStruct::name, &TaoPlotStruct::set_name)
      // TaoPlotStruct.description (0D_NOT_character - Descriptive string.
      .def_property("description", &TaoPlotStruct::description, &TaoPlotStruct::set_description)
      // TaoPlotStruct.graph (1D_ALLOC_type - individual graphs of a plot
      .def_property_readonly("graph", &TaoPlotStruct::graph)
      // TaoPlotStruct.r (0D_PTR_type - pointer to parent.
      .def_property("r", &TaoPlotStruct::r, &TaoPlotStruct::set_r)
      // TaoPlotStruct.ix_plot (0D_NOT_integer - Index in s%plot_page%template(:) or %region(:)
      // arrays.
      .def_property("ix_plot", &TaoPlotStruct::ix_plot, &TaoPlotStruct::set_ix_plot)
      // TaoPlotStruct.n_curve_pts (0D_NOT_integer - Overrides s%plot_page%n_curve_pts.
      .def_property("n_curve_pts", &TaoPlotStruct::n_curve_pts, &TaoPlotStruct::set_n_curve_pts)
      // TaoPlotStruct.type (0D_NOT_character - or 'wave'
      .def_property("type", &TaoPlotStruct::type, &TaoPlotStruct::set_type)
      // TaoPlotStruct.x_axis_type (0D_NOT_character - 'index', 'ele_index', 's', 'none', 'floor',
      // 'phase_space', etc.
      .def_property("x_axis_type", &TaoPlotStruct::x_axis_type, &TaoPlotStruct::set_x_axis_type)
      // TaoPlotStruct.autoscale_x (0D_NOT_logical - Horizontal autoscale.
      .def_property("autoscale_x", &TaoPlotStruct::autoscale_x, &TaoPlotStruct::set_autoscale_x)
      // TaoPlotStruct.autoscale_y (0D_NOT_logical - Vertical autoscale.
      .def_property("autoscale_y", &TaoPlotStruct::autoscale_y, &TaoPlotStruct::set_autoscale_y)
      // TaoPlotStruct.autoscale_gang_x (0D_NOT_logical - scale cmd scales graphs together?
      .def_property(
          "autoscale_gang_x",
          &TaoPlotStruct::autoscale_gang_x,
          &TaoPlotStruct::set_autoscale_gang_x
      )
      // TaoPlotStruct.autoscale_gang_y (0D_NOT_logical - scale cmd scales graphs together?
      .def_property(
          "autoscale_gang_y",
          &TaoPlotStruct::autoscale_gang_y,
          &TaoPlotStruct::set_autoscale_gang_y
      )
      // TaoPlotStruct.list_with_show_plot_command (0D_NOT_logical - False used for default plots to
      // shorten the output of 'show plot'
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotStruct::list_with_show_plot_command,
          &TaoPlotStruct::set_list_with_show_plot_command
      )
      // TaoPlotStruct.phantom (0D_NOT_logical - Used by tao_plot_init to add info lines to 'show
      // plot -templates'
      .def_property("phantom", &TaoPlotStruct::phantom, &TaoPlotStruct::set_phantom)
      // TaoPlotStruct.default_plot (0D_NOT_logical - One of Tao's default plots?
      .def_property("default_plot", &TaoPlotStruct::default_plot, &TaoPlotStruct::set_default_plot)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoPlotStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoPlotStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotStruct &self) {
            return TaoPlotStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotStruct &self, py::dict &memo) { return TaoPlotStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoPlotStructArray1D>(m, "TaoPlotStructArray1D");
  bind_FTypeAlloc1D<TaoPlotStructAlloc1D>(m, "TaoPlotStructAlloc1D");
  // 2D TaoPlotStruct arrays are not used in structs/routines
  // 3D TaoPlotStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_point_struct
void init_tao_shape_pattern_point_struct(
    py::module &m,
    py::class_<TaoShapePatternPointStruct> &cls
) {
  cls.def(py::init<>())
      // TaoShapePatternPointStruct.s (0D_NOT_real -
      .def_property("s", &TaoShapePatternPointStruct::s, &TaoShapePatternPointStruct::set_s)
      // TaoShapePatternPointStruct.y (0D_NOT_real -
      .def_property("y", &TaoShapePatternPointStruct::y, &TaoShapePatternPointStruct::set_y)
      // TaoShapePatternPointStruct.radius (0D_NOT_real -
      .def_property(
          "radius",
          &TaoShapePatternPointStruct::radius,
          &TaoShapePatternPointStruct::set_radius
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoShapePatternPointStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoShapePatternPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoShapePatternPointStruct &self) {
            return TaoShapePatternPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoShapePatternPointStruct &self, py::dict &memo) {
            return TaoShapePatternPointStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoShapePatternPointStructArray1D>(m, "TaoShapePatternPointStructArray1D");
  bind_FTypeAlloc1D<TaoShapePatternPointStructAlloc1D>(m, "TaoShapePatternPointStructAlloc1D");
  // 2D TaoShapePatternPointStruct arrays are not used in structs/routines
  // 3D TaoShapePatternPointStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_struct
void init_tao_shape_pattern_struct(py::module &m, py::class_<TaoShapePatternStruct> &cls) {
  cls.def(py::init<>())
      // TaoShapePatternStruct.name (0D_NOT_character -
      .def_property("name", &TaoShapePatternStruct::name, &TaoShapePatternStruct::set_name)
      // TaoShapePatternStruct.line (0D_NOT_type - Line color and pattern set by shape using this
      // pattern.
      .def_property("line", &TaoShapePatternStruct::line, &TaoShapePatternStruct::set_line)
      // TaoShapePatternStruct.pt (1D_ALLOC_type -
      .def_property_readonly("pt", &TaoShapePatternStruct::pt)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoShapePatternStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoShapePatternStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoShapePatternStruct &self) {
            return TaoShapePatternStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoShapePatternStruct &self, py::dict &memo) {
            return TaoShapePatternStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoShapePatternStructArray1D>(m, "TaoShapePatternStructArray1D");
  bind_FTypeAlloc1D<TaoShapePatternStructAlloc1D>(m, "TaoShapePatternStructAlloc1D");
  // 2D TaoShapePatternStruct arrays are not used in structs/routines
  // 3D TaoShapePatternStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_dn_dpz_struct
void init_tao_spin_dn_dpz_struct(py::module &m, py::class_<TaoSpinDnDpzStruct> &cls) {
  cls.def(py::init<>())
      // TaoSpinDnDpzStruct.vec (1D_NOT_real - n0 derivative wrt pz.
      .def_property("vec", &TaoSpinDnDpzStruct::vec, &TaoSpinDnDpzStruct::set_vec)
      // TaoSpinDnDpzStruct.partial (2D_NOT_real - partial(i:) is spin n0 derivative wrt pz for i^th
      // oscillation mode (1 => a-mode, etc.)
      .def_property("partial", &TaoSpinDnDpzStruct::partial, &TaoSpinDnDpzStruct::set_partial)
      // TaoSpinDnDpzStruct.partial2 (2D_NOT_real - partial(i:) is spin n0 derivative wrt pz with
      // i^th oscillation mode missing (1 => a-mode, etc.)
      .def_property("partial2", &TaoSpinDnDpzStruct::partial2, &TaoSpinDnDpzStruct::set_partial2)

      .def("__repr__", [](const TaoSpinDnDpzStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinDnDpzStruct &self) {
            return TaoSpinDnDpzStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoSpinDnDpzStruct &self, py::dict &memo) { return TaoSpinDnDpzStruct(self); }
      )

      ;

  // 1D TaoSpinDnDpzStruct arrays are not used in structs/routines
  // 2D TaoSpinDnDpzStruct arrays are not used in structs/routines
  // 3D TaoSpinDnDpzStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_ele_struct
void init_tao_spin_ele_struct(py::module &m, py::class_<TaoSpinEleStruct> &cls) {
  cls.def(py::init<>())
      // TaoSpinEleStruct.dn_dpz (0D_NOT_type -
      .def_property("dn_dpz", &TaoSpinEleStruct::dn_dpz, &TaoSpinEleStruct::set_dn_dpz)
      // TaoSpinEleStruct.orb_eigen_val (1D_NOT_real -
      .def_property(
          "orb_eigen_val",
          &TaoSpinEleStruct::orb_eigen_val,
          &TaoSpinEleStruct::set_orb_eigen_val
      )
      // TaoSpinEleStruct.orb_eigen_vec (2D_NOT_real - (j,:) is j^th vector
      .def_property(
          "orb_eigen_vec",
          &TaoSpinEleStruct::orb_eigen_vec,
          &TaoSpinEleStruct::set_orb_eigen_vec
      )
      // TaoSpinEleStruct.spin_eigen_vec (2D_NOT_real - (j,:) is j^th vector
      .def_property(
          "spin_eigen_vec",
          &TaoSpinEleStruct::spin_eigen_vec,
          &TaoSpinEleStruct::set_spin_eigen_vec
      )
      // TaoSpinEleStruct.valid (0D_NOT_logical -
      .def_property("valid", &TaoSpinEleStruct::valid, &TaoSpinEleStruct::set_valid)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoSpinEleStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoSpinEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinEleStruct &self) {
            return TaoSpinEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoSpinEleStruct &self, py::dict &memo) { return TaoSpinEleStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoSpinEleStructArray1D>(m, "TaoSpinEleStructArray1D");
  bind_FTypeAlloc1D<TaoSpinEleStructAlloc1D>(m, "TaoSpinEleStructAlloc1D");
  // 2D TaoSpinEleStruct arrays are not used in structs/routines
  // 3D TaoSpinEleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_map_struct
void init_tao_spin_map_struct(py::module &m, py::class_<TaoSpinMapStruct> &cls) {
  cls.def(py::init<>())
      // TaoSpinMapStruct.valid (0D_NOT_logical -
      .def_property("valid", &TaoSpinMapStruct::valid, &TaoSpinMapStruct::set_valid)
      // TaoSpinMapStruct.map1 (0D_NOT_type -
      .def_property("map1", &TaoSpinMapStruct::map1, &TaoSpinMapStruct::set_map1)
      // TaoSpinMapStruct.axis_input (0D_NOT_type - Input axes.
      .def_property("axis_input", &TaoSpinMapStruct::axis_input, &TaoSpinMapStruct::set_axis_input)
      // TaoSpinMapStruct.axis0 (0D_NOT_type - Initial axes.
      .def_property("axis0", &TaoSpinMapStruct::axis0, &TaoSpinMapStruct::set_axis0)
      // TaoSpinMapStruct.axis1 (0D_NOT_type - Final axes.
      .def_property("axis1", &TaoSpinMapStruct::axis1, &TaoSpinMapStruct::set_axis1)
      // TaoSpinMapStruct.ix_ele (0D_NOT_integer -
      .def_property("ix_ele", &TaoSpinMapStruct::ix_ele, &TaoSpinMapStruct::set_ix_ele)
      // TaoSpinMapStruct.ix_ref (0D_NOT_integer -
      .def_property("ix_ref", &TaoSpinMapStruct::ix_ref, &TaoSpinMapStruct::set_ix_ref)
      // TaoSpinMapStruct.ix_uni (0D_NOT_integer -
      .def_property("ix_uni", &TaoSpinMapStruct::ix_uni, &TaoSpinMapStruct::set_ix_uni)
      // TaoSpinMapStruct.ix_branch (0D_NOT_integer -
      .def_property("ix_branch", &TaoSpinMapStruct::ix_branch, &TaoSpinMapStruct::set_ix_branch)
      // TaoSpinMapStruct.mat8 (2D_NOT_real -
      .def_property("mat8", &TaoSpinMapStruct::mat8, &TaoSpinMapStruct::set_mat8)

      .def("__repr__", [](const TaoSpinMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinMapStruct &self) {
            return TaoSpinMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoSpinMapStruct &self, py::dict &memo) { return TaoSpinMapStruct(self); }
      )

      ;

  // 1D TaoSpinMapStruct arrays are not used in structs/routines
  // 2D TaoSpinMapStruct arrays are not used in structs/routines
  // 3D TaoSpinMapStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_polarization_struct
void init_tao_spin_polarization_struct(py::module &m, py::class_<TaoSpinPolarizationStruct> &cls) {
  cls.def(py::init<>())
      // TaoSpinPolarizationStruct.tune (0D_NOT_real -
      .def_property("tune", &TaoSpinPolarizationStruct::tune, &TaoSpinPolarizationStruct::set_tune)
      // TaoSpinPolarizationStruct.pol_limit_st (0D_NOT_real - Polarization calculated using
      // Sokolov-Ternov formula.
      .def_property(
          "pol_limit_st",
          &TaoSpinPolarizationStruct::pol_limit_st,
          &TaoSpinPolarizationStruct::set_pol_limit_st
      )
      // TaoSpinPolarizationStruct.pol_limit_dk (0D_NOT_real - Equalibrium Polarization calculated
      // via the Derbenev-Kondratenko-Mane formula.
      .def_property(
          "pol_limit_dk",
          &TaoSpinPolarizationStruct::pol_limit_dk,
          &TaoSpinPolarizationStruct::set_pol_limit_dk
      )
      // TaoSpinPolarizationStruct.pol_limit_dk_partial (1D_NOT_real - Limit using only single mode
      // to calc dn_dpz
      .def_property(
          "pol_limit_dk_partial",
          &TaoSpinPolarizationStruct::pol_limit_dk_partial,
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial
      )
      // TaoSpinPolarizationStruct.pol_limit_dk_partial2 (1D_NOT_real - Limit using only single mode
      // to calc dn_dpz
      .def_property(
          "pol_limit_dk_partial2",
          &TaoSpinPolarizationStruct::pol_limit_dk_partial2,
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial2
      )
      // TaoSpinPolarizationStruct.pol_rate_bks (0D_NOT_real - BKS Polarization rate (1/sec).
      .def_property(
          "pol_rate_bks",
          &TaoSpinPolarizationStruct::pol_rate_bks,
          &TaoSpinPolarizationStruct::set_pol_rate_bks
      )
      // TaoSpinPolarizationStruct.depol_rate (0D_NOT_real - Depolarization rate (1/sec).
      .def_property(
          "depol_rate",
          &TaoSpinPolarizationStruct::depol_rate,
          &TaoSpinPolarizationStruct::set_depol_rate
      )
      // TaoSpinPolarizationStruct.depol_rate_partial (1D_NOT_real - Depolarization rate (1/sec)
      // using only single mode to calc dn_dpz.
      .def_property(
          "depol_rate_partial",
          &TaoSpinPolarizationStruct::depol_rate_partial,
          &TaoSpinPolarizationStruct::set_depol_rate_partial
      )
      // TaoSpinPolarizationStruct.depol_rate_partial2 (1D_NOT_real - Depolarization rate (1/sec)
      // using only two modes to calc dn_dpz.
      .def_property(
          "depol_rate_partial2",
          &TaoSpinPolarizationStruct::depol_rate_partial2,
          &TaoSpinPolarizationStruct::set_depol_rate_partial2
      )
      // TaoSpinPolarizationStruct.integral_bn (0D_NOT_real - Integral of g^3 * b_hat * n_0
      .def_property(
          "integral_bn",
          &TaoSpinPolarizationStruct::integral_bn,
          &TaoSpinPolarizationStruct::set_integral_bn
      )
      // TaoSpinPolarizationStruct.integral_bdn (0D_NOT_real - Integral of g^3 * b_hat * dn/ddelta
      .def_property(
          "integral_bdn",
          &TaoSpinPolarizationStruct::integral_bdn,
          &TaoSpinPolarizationStruct::set_integral_bdn
      )
      // TaoSpinPolarizationStruct.integral_1ns (0D_NOT_real - Integral of g^3 (1 - 2(n * s_hat)/9)
      .def_property(
          "integral_1ns",
          &TaoSpinPolarizationStruct::integral_1ns,
          &TaoSpinPolarizationStruct::set_integral_1ns
      )
      // TaoSpinPolarizationStruct.integral_dn2 (0D_NOT_real - Integral of g^3 * 11 (dn/ddelta)^2 /
      // 9
      .def_property(
          "integral_dn2",
          &TaoSpinPolarizationStruct::integral_dn2,
          &TaoSpinPolarizationStruct::set_integral_dn2
      )
      // TaoSpinPolarizationStruct.valid (0D_NOT_logical -
      .def_property(
          "valid",
          &TaoSpinPolarizationStruct::valid,
          &TaoSpinPolarizationStruct::set_valid
      )
      // TaoSpinPolarizationStruct.q_1turn (0D_NOT_type - Save results from spin_concat_linear_maps
      // in tao_spin_polarization.
      .def_property(
          "q_1turn",
          &TaoSpinPolarizationStruct::q_1turn,
          &TaoSpinPolarizationStruct::set_q_1turn
      )
      // TaoSpinPolarizationStruct.q_ele (1D_ALLOC_type - Save results from spin_concat_linear_maps
      // in tao_spin_polarization.
      .def_property_readonly("q_ele", &TaoSpinPolarizationStruct::q_ele)

      .def("__repr__", [](const TaoSpinPolarizationStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSpinPolarizationStruct &self) {
            return TaoSpinPolarizationStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoSpinPolarizationStruct &self, py::dict &memo) {
            return TaoSpinPolarizationStruct(self);
          }
      )

      ;

  // 1D TaoSpinPolarizationStruct arrays are not used in structs/routines
  // 2D TaoSpinPolarizationStruct arrays are not used in structs/routines
  // 3D TaoSpinPolarizationStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_super_universe_struct
void init_tao_super_universe_struct(py::module &m, py::class_<TaoSuperUniverseStruct> &cls) {
  cls.def(py::init<>())
      // TaoSuperUniverseStruct.global (0D_NOT_type - User accessible global variables.
      .def_property("global_", &TaoSuperUniverseStruct::global, &TaoSuperUniverseStruct::set_global)
      // TaoSuperUniverseStruct.init (0D_NOT_type - Initialization parameters
      .def_property("init", &TaoSuperUniverseStruct::init, &TaoSuperUniverseStruct::set_init)
      // TaoSuperUniverseStruct.com (0D_NOT_type - Non-initialization common parameters
      .def_property("com", &TaoSuperUniverseStruct::com, &TaoSuperUniverseStruct::set_com)
      // TaoSuperUniverseStruct.plot_page (0D_NOT_type - Defines the plot window.
      .def_property(
          "plot_page",
          &TaoSuperUniverseStruct::plot_page,
          &TaoSuperUniverseStruct::set_plot_page
      )
      // TaoSuperUniverseStruct.v1_var (1D_ALLOC_type - The variable types
      .def_property_readonly("v1_var", &TaoSuperUniverseStruct::v1_var)
      // TaoSuperUniverseStruct.var (1D_ALLOC_type - array of all variables.
      .def_property_readonly("var", &TaoSuperUniverseStruct::var)
      // TaoSuperUniverseStruct.u (1D_ALLOC_type - array of universes.
      .def_property_readonly("u", &TaoSuperUniverseStruct::u)
      // TaoSuperUniverseStruct.key (1D_ALLOC_integer -
      .def_property("key", &TaoSuperUniverseStruct::key, &TaoSuperUniverseStruct::set_key)
      // TaoSuperUniverseStruct.building_wall (0D_NOT_type -
      .def_property(
          "building_wall",
          &TaoSuperUniverseStruct::building_wall,
          &TaoSuperUniverseStruct::set_building_wall
      )
      // TaoSuperUniverseStruct.wave (0D_NOT_type -
      .def_property("wave", &TaoSuperUniverseStruct::wave, &TaoSuperUniverseStruct::set_wave)
      // TaoSuperUniverseStruct.n_var_used (0D_NOT_integer -
      .def_property(
          "n_var_used",
          &TaoSuperUniverseStruct::n_var_used,
          &TaoSuperUniverseStruct::set_n_var_used
      )
      // TaoSuperUniverseStruct.n_v1_var_used (0D_NOT_integer -
      .def_property(
          "n_v1_var_used",
          &TaoSuperUniverseStruct::n_v1_var_used,
          &TaoSuperUniverseStruct::set_n_v1_var_used
      )
      // TaoSuperUniverseStruct.history (1D_NOT_type - command history
      .def_property_readonly("history", &TaoSuperUniverseStruct::history)
      // TaoSuperUniverseStruct.initialized (0D_NOT_logical - Does tao_init() need to be called?
      .def_property(
          "initialized",
          &TaoSuperUniverseStruct::initialized,
          &TaoSuperUniverseStruct::set_initialized
      )

      .def("__repr__", [](const TaoSuperUniverseStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoSuperUniverseStruct &self) {
            return TaoSuperUniverseStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoSuperUniverseStruct &self, py::dict &memo) {
            return TaoSuperUniverseStruct(self);
          }
      )

      ;

  // 1D TaoSuperUniverseStruct arrays are not used in structs/routines
  // 2D TaoSuperUniverseStruct arrays are not used in structs/routines
  // 3D TaoSuperUniverseStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_title_struct
void init_tao_title_struct(py::module &m, py::class_<TaoTitleStruct> &cls) {
  cls.def(py::init<>())
      // TaoTitleStruct.string (0D_NOT_character - title character string.
      .def_property("string", &TaoTitleStruct::string, &TaoTitleStruct::set_string)
      // TaoTitleStruct.x (0D_NOT_real - x, y rwt lower left corner
      .def_property("x", &TaoTitleStruct::x, &TaoTitleStruct::set_x)
      // TaoTitleStruct.y (0D_NOT_real - x, y rwt lower left corner
      .def_property("y", &TaoTitleStruct::y, &TaoTitleStruct::set_y)
      // TaoTitleStruct.units (0D_NOT_character - %BOX, POINTS, etc...
      .def_property("units", &TaoTitleStruct::units, &TaoTitleStruct::set_units)
      // TaoTitleStruct.justify (0D_NOT_character - Left, Center, or Right justification.
      .def_property("justify", &TaoTitleStruct::justify, &TaoTitleStruct::set_justify)
      // TaoTitleStruct.draw_it (0D_NOT_logical - draw the title?
      .def_property("draw_it", &TaoTitleStruct::draw_it, &TaoTitleStruct::set_draw_it)

      .def("__repr__", [](const TaoTitleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoTitleStruct &self) {
            return TaoTitleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoTitleStruct &self, py::dict &memo) { return TaoTitleStruct(self); }
      )

      ;

  // 1D TaoTitleStruct arrays are not used in structs/routines
  // 2D TaoTitleStruct arrays are not used in structs/routines
  // 3D TaoTitleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_calc_struct
void init_tao_universe_calc_struct(py::module &m, py::class_<TaoUniverseCalcStruct> &cls) {
  cls.def(py::init<>())
      // TaoUniverseCalcStruct.srdt_for_data (0D_NOT_integer - 0 = false, 1 = 1st order, 2 = 1st &
      // 2nd order
      .def_property(
          "srdt_for_data",
          &TaoUniverseCalcStruct::srdt_for_data,
          &TaoUniverseCalcStruct::set_srdt_for_data
      )
      // TaoUniverseCalcStruct.rad_int_for_data (0D_NOT_logical - Do the radiation integrals need to
      // be computed for
      .def_property(
          "rad_int_for_data",
          &TaoUniverseCalcStruct::rad_int_for_data,
          &TaoUniverseCalcStruct::set_rad_int_for_data
      )
      // TaoUniverseCalcStruct.rad_int_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "rad_int_for_plotting",
          &TaoUniverseCalcStruct::rad_int_for_plotting,
          &TaoUniverseCalcStruct::set_rad_int_for_plotting
      )
      // TaoUniverseCalcStruct.chrom_for_data (0D_NOT_logical - Does the chromaticity need to be
      // computed for
      .def_property(
          "chrom_for_data",
          &TaoUniverseCalcStruct::chrom_for_data,
          &TaoUniverseCalcStruct::set_chrom_for_data
      )
      // TaoUniverseCalcStruct.chrom_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "chrom_for_plotting",
          &TaoUniverseCalcStruct::chrom_for_plotting,
          &TaoUniverseCalcStruct::set_chrom_for_plotting
      )
      // TaoUniverseCalcStruct.lat_sigma_for_data (0D_NOT_logical - Do the beam sigmas need to be
      // computed for
      .def_property(
          "lat_sigma_for_data",
          &TaoUniverseCalcStruct::lat_sigma_for_data,
          &TaoUniverseCalcStruct::set_lat_sigma_for_data
      )
      // TaoUniverseCalcStruct.lat_sigma_for_plotting (0D_NOT_logical - data or plotting?
      .def_property(
          "lat_sigma_for_plotting",
          &TaoUniverseCalcStruct::lat_sigma_for_plotting,
          &TaoUniverseCalcStruct::set_lat_sigma_for_plotting
      )
      // TaoUniverseCalcStruct.dynamic_aperture (0D_NOT_logical - Do the dynamic_aperture calc?
      .def_property(
          "dynamic_aperture",
          &TaoUniverseCalcStruct::dynamic_aperture,
          &TaoUniverseCalcStruct::set_dynamic_aperture
      )
      // TaoUniverseCalcStruct.one_turn_map (0D_NOT_logical - Compute the one turn map?
      .def_property(
          "one_turn_map",
          &TaoUniverseCalcStruct::one_turn_map,
          &TaoUniverseCalcStruct::set_one_turn_map
      )
      // TaoUniverseCalcStruct.lattice (0D_NOT_logical - Used to indicate which lattices need
      // tracking done.
      .def_property("lattice", &TaoUniverseCalcStruct::lattice, &TaoUniverseCalcStruct::set_lattice)
      // TaoUniverseCalcStruct.twiss (0D_NOT_logical - calc linear transfer matrix?
      .def_property("twiss", &TaoUniverseCalcStruct::twiss, &TaoUniverseCalcStruct::set_twiss)
      // TaoUniverseCalcStruct.track (0D_NOT_logical - tracking needs to be done?
      .def_property("track", &TaoUniverseCalcStruct::track, &TaoUniverseCalcStruct::set_track)
      // TaoUniverseCalcStruct.spin_matrices (0D_NOT_logical - Calculate G and D spin matrices?
      .def_property(
          "spin_matrices",
          &TaoUniverseCalcStruct::spin_matrices,
          &TaoUniverseCalcStruct::set_spin_matrices
      )

      .def("__repr__", [](const TaoUniverseCalcStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniverseCalcStruct &self) {
            return TaoUniverseCalcStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoUniverseCalcStruct &self, py::dict &memo) {
            return TaoUniverseCalcStruct(self);
          }
      )

      ;

  // 1D TaoUniverseCalcStruct arrays are not used in structs/routines
  // 2D TaoUniverseCalcStruct arrays are not used in structs/routines
  // 3D TaoUniverseCalcStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_pointer_struct
void init_tao_universe_pointer_struct(py::module &m, py::class_<TaoUniversePointerStruct> &cls) {
  cls.def(py::init<>())
      // TaoUniversePointerStruct.u (0D_PTR_type -
      .def_property("u", &TaoUniversePointerStruct::u, &TaoUniversePointerStruct::set_u)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoUniversePointerStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoUniversePointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniversePointerStruct &self) {
            return TaoUniversePointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoUniversePointerStruct &self, py::dict &memo) {
            return TaoUniversePointerStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<TaoUniversePointerStructArray1D>(m, "TaoUniversePointerStructArray1D");
  bind_FTypeAlloc1D<TaoUniversePointerStructAlloc1D>(m, "TaoUniversePointerStructAlloc1D");
  // 2D TaoUniversePointerStruct arrays are not used in structs/routines
  // 3D TaoUniversePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_struct
void init_tao_universe_struct(py::module &m, py::class_<TaoUniverseStruct> &cls) {
  cls.def(py::init<>())
      // TaoUniverseStruct.model (0D_PTR_type -
      .def_property("model", &TaoUniverseStruct::model, &TaoUniverseStruct::set_model)
      // TaoUniverseStruct.design (0D_PTR_type -
      .def_property("design", &TaoUniverseStruct::design, &TaoUniverseStruct::set_design)
      // TaoUniverseStruct.base (0D_PTR_type -
      .def_property("base", &TaoUniverseStruct::base, &TaoUniverseStruct::set_base)
      // TaoUniverseStruct.beam (0D_NOT_type -
      .def_property("beam", &TaoUniverseStruct::beam, &TaoUniverseStruct::set_beam)
      // TaoUniverseStruct.dynamic_aperture (0D_NOT_type -
      .def_property(
          "dynamic_aperture",
          &TaoUniverseStruct::dynamic_aperture,
          &TaoUniverseStruct::set_dynamic_aperture
      )
      // TaoUniverseStruct.model_branch (1D_PTR_type - model specific information
      .def_property_readonly("model_branch", &TaoUniverseStruct::model_branch)
      // TaoUniverseStruct.d2_data (1D_ALLOC_type - The data types
      .def_property_readonly("d2_data", &TaoUniverseStruct::d2_data)
      // TaoUniverseStruct.data (1D_ALLOC_type - Array of all data.
      .def_property_readonly("data", &TaoUniverseStruct::data)
      // TaoUniverseStruct.ping_scale (0D_NOT_type -
      .def_property(
          "ping_scale",
          &TaoUniverseStruct::ping_scale,
          &TaoUniverseStruct::set_ping_scale
      )
      // TaoUniverseStruct.scratch_lat (0D_NOT_type - Scratch area.
      .def_property(
          "scratch_lat",
          &TaoUniverseStruct::scratch_lat,
          &TaoUniverseStruct::set_scratch_lat
      )
      // TaoUniverseStruct.calc (0D_NOT_type - What needs to be calculated?
      .def_property("calc", &TaoUniverseStruct::calc, &TaoUniverseStruct::set_calc)
      // TaoUniverseStruct.ele_order (0D_NOT_type - Order of elements with same name.
      .def_property("ele_order", &TaoUniverseStruct::ele_order, &TaoUniverseStruct::set_ele_order)
      // TaoUniverseStruct.spin_map (0D_NOT_type -
      .def_property("spin_map", &TaoUniverseStruct::spin_map, &TaoUniverseStruct::set_spin_map)
      // TaoUniverseStruct.dModel_dVar (2D_ALLOC_real - Derivative matrix.
      .def_property(
          "dModel_dVar",
          &TaoUniverseStruct::dModel_dVar,
          &TaoUniverseStruct::set_dModel_dVar
      )
      // TaoUniverseStruct.ix_uni (0D_NOT_integer - Universe index.
      .def_property("ix_uni", &TaoUniverseStruct::ix_uni, &TaoUniverseStruct::set_ix_uni)
      // TaoUniverseStruct.n_d2_data_used (0D_NOT_integer - Number of used %d2_data(:) components.
      .def_property(
          "n_d2_data_used",
          &TaoUniverseStruct::n_d2_data_used,
          &TaoUniverseStruct::set_n_d2_data_used
      )
      // TaoUniverseStruct.n_data_used (0D_NOT_integer - Number of used %data(:) components.
      .def_property(
          "n_data_used",
          &TaoUniverseStruct::n_data_used,
          &TaoUniverseStruct::set_n_data_used
      )
      // TaoUniverseStruct.is_on (0D_NOT_logical - universe turned on
      .def_property("is_on", &TaoUniverseStruct::is_on, &TaoUniverseStruct::set_is_on)
      // TaoUniverseStruct.design_same_as_previous (0D_NOT_logical - Design lat same as the previous
      // uni?
      .def_property(
          "design_same_as_previous",
          &TaoUniverseStruct::design_same_as_previous,
          &TaoUniverseStruct::set_design_same_as_previous
      )
      // TaoUniverseStruct.picked_uni (0D_NOT_logical - Scratch logical.
      .def_property(
          "picked_uni",
          &TaoUniverseStruct::picked_uni,
          &TaoUniverseStruct::set_picked_uni
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoUniverseStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoUniverseStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoUniverseStruct &self) {
            return TaoUniverseStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoUniverseStruct &self, py::dict &memo) { return TaoUniverseStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoUniverseStructArray1D>(m, "TaoUniverseStructArray1D");
  bind_FTypeAlloc1D<TaoUniverseStructAlloc1D>(m, "TaoUniverseStructAlloc1D");
  // 2D TaoUniverseStruct arrays are not used in structs/routines
  // 3D TaoUniverseStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_v1_var_struct
void init_tao_v1_var_struct(py::module &m, py::class_<TaoV1VarStruct> &cls) {
  cls.def(py::init<>())
      // TaoV1VarStruct.name (0D_NOT_character - V1 variable name. Eg: 'quad_k1'.
      .def_property("name", &TaoV1VarStruct::name, &TaoV1VarStruct::set_name)
      // TaoV1VarStruct.ix_v1_var (0D_NOT_integer - Index to s%v1_var(:) array
      .def_property("ix_v1_var", &TaoV1VarStruct::ix_v1_var, &TaoV1VarStruct::set_ix_v1_var)
      // TaoV1VarStruct.v (1D_PTR_type - Pointer to the appropriate section in s%var.
      .def_property_readonly("v", &TaoV1VarStruct::v)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoV1VarStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoV1VarStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoV1VarStruct &self) {
            return TaoV1VarStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoV1VarStruct &self, py::dict &memo) { return TaoV1VarStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoV1VarStructArray1D>(m, "TaoV1VarStructArray1D");
  bind_FTypeAlloc1D<TaoV1VarStructAlloc1D>(m, "TaoV1VarStructAlloc1D");
  // 2D TaoV1VarStruct arrays are not used in structs/routines
  // 3D TaoV1VarStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_var_slave_struct
void init_tao_var_slave_struct(py::module &m, py::class_<TaoVarSlaveStruct> &cls) {
  cls.def(py::init<>())
      // TaoVarSlaveStruct.ix_uni (0D_NOT_integer - universe index.
      .def_property("ix_uni", &TaoVarSlaveStruct::ix_uni, &TaoVarSlaveStruct::set_ix_uni)
      // TaoVarSlaveStruct.ix_branch (0D_NOT_integer -
      .def_property("ix_branch", &TaoVarSlaveStruct::ix_branch, &TaoVarSlaveStruct::set_ix_branch)
      // TaoVarSlaveStruct.ix_ele (0D_NOT_integer - Index of element in the u%lattice%ele(:) array.
      .def_property("ix_ele", &TaoVarSlaveStruct::ix_ele, &TaoVarSlaveStruct::set_ix_ele)
      // TaoVarSlaveStruct.model_value (0D_PTR_real - Pointer to the variable in the model lat.
      .def_property(
          "model_value",
          &TaoVarSlaveStruct::model_value,
          &TaoVarSlaveStruct::set_model_value
      )
      // TaoVarSlaveStruct.base_value (0D_PTR_real - Pointer to the variable in the base lat.
      .def_property(
          "base_value",
          &TaoVarSlaveStruct::base_value,
          &TaoVarSlaveStruct::set_base_value
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoVarSlaveStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoVarSlaveStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoVarSlaveStruct &self) {
            return TaoVarSlaveStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoVarSlaveStruct &self, py::dict &memo) { return TaoVarSlaveStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoVarSlaveStructArray1D>(m, "TaoVarSlaveStructArray1D");
  bind_FTypeAlloc1D<TaoVarSlaveStructAlloc1D>(m, "TaoVarSlaveStructAlloc1D");
  // 2D TaoVarSlaveStruct arrays are not used in structs/routines
  // 3D TaoVarSlaveStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_var_struct
void init_tao_var_struct(py::module &m, py::class_<TaoVarStruct> &cls) {
  cls.def(py::init<>())
      // TaoVarStruct.ele_name (0D_NOT_character - Associated lattice element name.
      .def_property("ele_name", &TaoVarStruct::ele_name, &TaoVarStruct::set_ele_name)
      // TaoVarStruct.attrib_name (0D_NOT_character - Name of the attribute to vary.
      .def_property("attrib_name", &TaoVarStruct::attrib_name, &TaoVarStruct::set_attrib_name)
      // TaoVarStruct.id (0D_NOT_character - Used by Tao extension code. Not used by Tao directly.
      .def_property("id", &TaoVarStruct::id, &TaoVarStruct::set_id)
      // TaoVarStruct.slave (1D_ALLOC_type -
      .def_property_readonly("slave", &TaoVarStruct::slave)
      // TaoVarStruct.ix_v1 (0D_NOT_integer - Index of this var in the s%v1_var(i)%v(:) array.
      .def_property("ix_v1", &TaoVarStruct::ix_v1, &TaoVarStruct::set_ix_v1)
      // TaoVarStruct.ix_var (0D_NOT_integer - Index number of this var in the s%var(:) array.
      .def_property("ix_var", &TaoVarStruct::ix_var, &TaoVarStruct::set_ix_var)
      // TaoVarStruct.ix_dvar (0D_NOT_integer - Column in the dData_dVar derivative matrix.
      .def_property("ix_dvar", &TaoVarStruct::ix_dvar, &TaoVarStruct::set_ix_dvar)
      // TaoVarStruct.ix_attrib (0D_NOT_integer - Index in ele%value(:) array if appropriate.
      .def_property("ix_attrib", &TaoVarStruct::ix_attrib, &TaoVarStruct::set_ix_attrib)
      // TaoVarStruct.ix_key_table (0D_NOT_integer - Has a key binding?
      .def_property("ix_key_table", &TaoVarStruct::ix_key_table, &TaoVarStruct::set_ix_key_table)
      // TaoVarStruct.model_value (0D_PTR_real - Model value.
      .def_property("model_value", &TaoVarStruct::model_value, &TaoVarStruct::set_model_value)
      // TaoVarStruct.base_value (0D_PTR_real - Base value.
      .def_property("base_value", &TaoVarStruct::base_value, &TaoVarStruct::set_base_value)
      // TaoVarStruct.design_value (0D_NOT_real - Design value from the design lattice.
      .def_property("design_value", &TaoVarStruct::design_value, &TaoVarStruct::set_design_value)
      // TaoVarStruct.scratch_value (0D_NOT_real - Scratch space used by Tao.
      .def_property("scratch_value", &TaoVarStruct::scratch_value, &TaoVarStruct::set_scratch_value)
      // TaoVarStruct.old_value (0D_NOT_real - Scratch space used by Tao.
      .def_property("old_value", &TaoVarStruct::old_value, &TaoVarStruct::set_old_value)
      // TaoVarStruct.meas_value (0D_NOT_real - The value when the data measurement was taken.
      .def_property("meas_value", &TaoVarStruct::meas_value, &TaoVarStruct::set_meas_value)
      // TaoVarStruct.ref_value (0D_NOT_real - Value when the reference measurement was taken.
      .def_property("ref_value", &TaoVarStruct::ref_value, &TaoVarStruct::set_ref_value)
      // TaoVarStruct.correction_value (0D_NOT_real - Value determined by a fit to correct the
      // lattice.
      .def_property(
          "correction_value",
          &TaoVarStruct::correction_value,
          &TaoVarStruct::set_correction_value
      )
      // TaoVarStruct.high_lim (0D_NOT_real - High limit for the model_value.
      .def_property("high_lim", &TaoVarStruct::high_lim, &TaoVarStruct::set_high_lim)
      // TaoVarStruct.low_lim (0D_NOT_real - Low limit for the model_value.
      .def_property("low_lim", &TaoVarStruct::low_lim, &TaoVarStruct::set_low_lim)
      // TaoVarStruct.step (0D_NOT_real - Sets what is a small step for varying this var.
      .def_property("step", &TaoVarStruct::step, &TaoVarStruct::set_step)
      // TaoVarStruct.weight (0D_NOT_real - Weight for the merit function term.
      .def_property("weight", &TaoVarStruct::weight, &TaoVarStruct::set_weight)
      // TaoVarStruct.delta_merit (0D_NOT_real - Diff used to calculate the merit function term.
      .def_property("delta_merit", &TaoVarStruct::delta_merit, &TaoVarStruct::set_delta_merit)
      // TaoVarStruct.merit (0D_NOT_real - merit_term = weight * delta^2.
      .def_property("merit", &TaoVarStruct::merit, &TaoVarStruct::set_merit)
      // TaoVarStruct.dMerit_dVar (0D_NOT_real - Merit derivative.
      .def_property("dMerit_dVar", &TaoVarStruct::dMerit_dVar, &TaoVarStruct::set_dMerit_dVar)
      // TaoVarStruct.key_val0 (0D_NOT_real - Key base value
      .def_property("key_val0", &TaoVarStruct::key_val0, &TaoVarStruct::set_key_val0)
      // TaoVarStruct.key_delta (0D_NOT_real - Change in value when a key is pressed.
      .def_property("key_delta", &TaoVarStruct::key_delta, &TaoVarStruct::set_key_delta)
      // TaoVarStruct.s (0D_NOT_real - longitudinal position of ele.
      .def_property("s", &TaoVarStruct::s, &TaoVarStruct::set_s)
      // TaoVarStruct.extend_val (0D_NOT_real - For extension code. Not used by Tao.
      .def_property("extend_val", &TaoVarStruct::extend_val, &TaoVarStruct::set_extend_val)
      // TaoVarStruct.merit_type (0D_NOT_character - 'target' or 'limit'
      .def_property("merit_type", &TaoVarStruct::merit_type, &TaoVarStruct::set_merit_type)
      // TaoVarStruct.exists (0D_NOT_logical - See above
      .def_property("exists", &TaoVarStruct::exists, &TaoVarStruct::set_exists)
      // TaoVarStruct.good_var (0D_NOT_logical - See above
      .def_property("good_var", &TaoVarStruct::good_var, &TaoVarStruct::set_good_var)
      // TaoVarStruct.good_user (0D_NOT_logical - See above
      .def_property("good_user", &TaoVarStruct::good_user, &TaoVarStruct::set_good_user)
      // TaoVarStruct.good_opt (0D_NOT_logical - See above
      .def_property("good_opt", &TaoVarStruct::good_opt, &TaoVarStruct::set_good_opt)
      // TaoVarStruct.good_plot (0D_NOT_logical - See above
      .def_property("good_plot", &TaoVarStruct::good_plot, &TaoVarStruct::set_good_plot)
      // TaoVarStruct.useit_opt (0D_NOT_logical - See above
      .def_property("useit_opt", &TaoVarStruct::useit_opt, &TaoVarStruct::set_useit_opt)
      // TaoVarStruct.useit_plot (0D_NOT_logical - See above
      .def_property("useit_plot", &TaoVarStruct::useit_plot, &TaoVarStruct::set_useit_plot)
      // TaoVarStruct.key_bound (0D_NOT_logical - Variable bound to keyboard key?
      .def_property("key_bound", &TaoVarStruct::key_bound, &TaoVarStruct::set_key_bound)
      // TaoVarStruct.v1 (0D_PTR_type - Pointer to the parent.
      .def_property("v1", &TaoVarStruct::v1, &TaoVarStruct::set_v1)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoVarStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoVarStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoVarStruct &self) {
            return TaoVarStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoVarStruct &self, py::dict &memo) { return TaoVarStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoVarStructArray1D>(m, "TaoVarStructArray1D");
  bind_FTypeAlloc1D<TaoVarStructAlloc1D>(m, "TaoVarStructAlloc1D");
  // 2D TaoVarStruct arrays are not used in structs/routines
  // 3D TaoVarStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_kick_pt_struct
void init_tao_wave_kick_pt_struct(py::module &m, py::class_<TaoWaveKickPtStruct> &cls) {
  cls.def(py::init<>())
      // TaoWaveKickPtStruct.phi_s (0D_NOT_real -
      .def_property("phi_s", &TaoWaveKickPtStruct::phi_s, &TaoWaveKickPtStruct::set_phi_s)
      // TaoWaveKickPtStruct.phi_r (0D_NOT_real -
      .def_property("phi_r", &TaoWaveKickPtStruct::phi_r, &TaoWaveKickPtStruct::set_phi_r)
      // TaoWaveKickPtStruct.phi (0D_NOT_real -
      .def_property("phi", &TaoWaveKickPtStruct::phi, &TaoWaveKickPtStruct::set_phi)
      // TaoWaveKickPtStruct.amp (0D_NOT_real -
      .def_property("amp", &TaoWaveKickPtStruct::amp, &TaoWaveKickPtStruct::set_amp)
      // TaoWaveKickPtStruct.s (0D_NOT_real - s-position of kick
      .def_property("s", &TaoWaveKickPtStruct::s, &TaoWaveKickPtStruct::set_s)
      // TaoWaveKickPtStruct.ix_dat_before_kick (0D_NOT_integer - Index of datum in data array just
      // before the kick.
      .def_property(
          "ix_dat_before_kick",
          &TaoWaveKickPtStruct::ix_dat_before_kick,
          &TaoWaveKickPtStruct::set_ix_dat_before_kick
      )
      // TaoWaveKickPtStruct.ele (0D_PTR_type - lattice element at position of kick.
      .def_property("ele", &TaoWaveKickPtStruct::ele, &TaoWaveKickPtStruct::set_ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TaoWaveKickPtStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TaoWaveKickPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoWaveKickPtStruct &self) {
            return TaoWaveKickPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoWaveKickPtStruct &self, py::dict &memo) { return TaoWaveKickPtStruct(self); }
      )

      ;

  bind_FTypeArrayND<TaoWaveKickPtStructArray1D>(m, "TaoWaveKickPtStructArray1D");
  bind_FTypeAlloc1D<TaoWaveKickPtStructAlloc1D>(m, "TaoWaveKickPtStructAlloc1D");
  // 2D TaoWaveKickPtStruct arrays are not used in structs/routines
  // 3D TaoWaveKickPtStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_struct
void init_tao_wave_struct(py::module &m, py::class_<TaoWaveStruct> &cls) {
  cls.def(py::init<>())
      // TaoWaveStruct.data_type (0D_NOT_character -
      .def_property("data_type", &TaoWaveStruct::data_type, &TaoWaveStruct::set_data_type)
      // TaoWaveStruct.rms_rel_a (0D_NOT_real -
      .def_property("rms_rel_a", &TaoWaveStruct::rms_rel_a, &TaoWaveStruct::set_rms_rel_a)
      // TaoWaveStruct.rms_rel_b (0D_NOT_real -
      .def_property("rms_rel_b", &TaoWaveStruct::rms_rel_b, &TaoWaveStruct::set_rms_rel_b)
      // TaoWaveStruct.rms_rel_as (0D_NOT_real -
      .def_property("rms_rel_as", &TaoWaveStruct::rms_rel_as, &TaoWaveStruct::set_rms_rel_as)
      // TaoWaveStruct.rms_rel_bs (0D_NOT_real -
      .def_property("rms_rel_bs", &TaoWaveStruct::rms_rel_bs, &TaoWaveStruct::set_rms_rel_bs)
      // TaoWaveStruct.rms_rel_ar (0D_NOT_real -
      .def_property("rms_rel_ar", &TaoWaveStruct::rms_rel_ar, &TaoWaveStruct::set_rms_rel_ar)
      // TaoWaveStruct.rms_rel_br (0D_NOT_real -
      .def_property("rms_rel_br", &TaoWaveStruct::rms_rel_br, &TaoWaveStruct::set_rms_rel_br)
      // TaoWaveStruct.rms_rel_k (0D_NOT_real -
      .def_property("rms_rel_k", &TaoWaveStruct::rms_rel_k, &TaoWaveStruct::set_rms_rel_k)
      // TaoWaveStruct.rms_rel_ks (0D_NOT_real -
      .def_property("rms_rel_ks", &TaoWaveStruct::rms_rel_ks, &TaoWaveStruct::set_rms_rel_ks)
      // TaoWaveStruct.rms_rel_kr (0D_NOT_real -
      .def_property("rms_rel_kr", &TaoWaveStruct::rms_rel_kr, &TaoWaveStruct::set_rms_rel_kr)
      // TaoWaveStruct.rms_phi (0D_NOT_real -
      .def_property("rms_phi", &TaoWaveStruct::rms_phi, &TaoWaveStruct::set_rms_phi)
      // TaoWaveStruct.rms_phi_s (0D_NOT_real -
      .def_property("rms_phi_s", &TaoWaveStruct::rms_phi_s, &TaoWaveStruct::set_rms_phi_s)
      // TaoWaveStruct.rms_phi_r (0D_NOT_real -
      .def_property("rms_phi_r", &TaoWaveStruct::rms_phi_r, &TaoWaveStruct::set_rms_phi_r)
      // TaoWaveStruct.amp_ba_s (0D_NOT_real -
      .def_property("amp_ba_s", &TaoWaveStruct::amp_ba_s, &TaoWaveStruct::set_amp_ba_s)
      // TaoWaveStruct.amp_ba_r (0D_NOT_real -
      .def_property("amp_ba_r", &TaoWaveStruct::amp_ba_r, &TaoWaveStruct::set_amp_ba_r)
      // TaoWaveStruct.chi_a (0D_NOT_real -
      .def_property("chi_a", &TaoWaveStruct::chi_a, &TaoWaveStruct::set_chi_a)
      // TaoWaveStruct.chi_c (0D_NOT_real -
      .def_property("chi_c", &TaoWaveStruct::chi_c, &TaoWaveStruct::set_chi_c)
      // TaoWaveStruct.chi_ba (0D_NOT_real -
      .def_property("chi_ba", &TaoWaveStruct::chi_ba, &TaoWaveStruct::set_chi_ba)
      // TaoWaveStruct.amp_a (1D_NOT_real -
      .def_property("amp_a", &TaoWaveStruct::amp_a, &TaoWaveStruct::set_amp_a)
      // TaoWaveStruct.amp_b (1D_NOT_real -
      .def_property("amp_b", &TaoWaveStruct::amp_b, &TaoWaveStruct::set_amp_b)
      // TaoWaveStruct.amp_ba (1D_NOT_real -
      .def_property("amp_ba", &TaoWaveStruct::amp_ba, &TaoWaveStruct::set_amp_ba)
      // TaoWaveStruct.coef_a (1D_NOT_real -
      .def_property("coef_a", &TaoWaveStruct::coef_a, &TaoWaveStruct::set_coef_a)
      // TaoWaveStruct.coef_b (1D_NOT_real -
      .def_property("coef_b", &TaoWaveStruct::coef_b, &TaoWaveStruct::set_coef_b)
      // TaoWaveStruct.coef_ba (1D_NOT_real -
      .def_property("coef_ba", &TaoWaveStruct::coef_ba, &TaoWaveStruct::set_coef_ba)
      // TaoWaveStruct.n_func (0D_NOT_integer - Number of functions used in the fit.
      .def_property("n_func", &TaoWaveStruct::n_func, &TaoWaveStruct::set_n_func)
      // TaoWaveStruct.ix_a1 (0D_NOT_integer -
      .def_property("ix_a1", &TaoWaveStruct::ix_a1, &TaoWaveStruct::set_ix_a1)
      // TaoWaveStruct.ix_a2 (0D_NOT_integer -
      .def_property("ix_a2", &TaoWaveStruct::ix_a2, &TaoWaveStruct::set_ix_a2)
      // TaoWaveStruct.ix_b1 (0D_NOT_integer -
      .def_property("ix_b1", &TaoWaveStruct::ix_b1, &TaoWaveStruct::set_ix_b1)
      // TaoWaveStruct.ix_b2 (0D_NOT_integer -
      .def_property("ix_b2", &TaoWaveStruct::ix_b2, &TaoWaveStruct::set_ix_b2)
      // TaoWaveStruct.i_a1 (0D_NOT_integer -
      .def_property("i_a1", &TaoWaveStruct::i_a1, &TaoWaveStruct::set_i_a1)
      // TaoWaveStruct.i_a2 (0D_NOT_integer -
      .def_property("i_a2", &TaoWaveStruct::i_a2, &TaoWaveStruct::set_i_a2)
      // TaoWaveStruct.i_b1 (0D_NOT_integer -
      .def_property("i_b1", &TaoWaveStruct::i_b1, &TaoWaveStruct::set_i_b1)
      // TaoWaveStruct.i_b2 (0D_NOT_integer -
      .def_property("i_b2", &TaoWaveStruct::i_b2, &TaoWaveStruct::set_i_b2)
      // TaoWaveStruct.n_a (0D_NOT_integer -
      .def_property("n_a", &TaoWaveStruct::n_a, &TaoWaveStruct::set_n_a)
      // TaoWaveStruct.n_b (0D_NOT_integer -
      .def_property("n_b", &TaoWaveStruct::n_b, &TaoWaveStruct::set_n_b)
      // TaoWaveStruct.i_curve_wrap_pt (0D_NOT_integer - Index of last point before wrap in curve
      // array.
      .def_property(
          "i_curve_wrap_pt",
          &TaoWaveStruct::i_curve_wrap_pt,
          &TaoWaveStruct::set_i_curve_wrap_pt
      )
      // TaoWaveStruct.ix_data (1D_ALLOC_integer - Translates from plot point to datum index
      .def_property("ix_data", &TaoWaveStruct::ix_data, &TaoWaveStruct::set_ix_data)
      // TaoWaveStruct.n_kick (0D_NOT_integer -
      .def_property("n_kick", &TaoWaveStruct::n_kick, &TaoWaveStruct::set_n_kick)
      // TaoWaveStruct.kick (1D_ALLOC_type -
      .def_property_readonly("kick", &TaoWaveStruct::kick)
      // TaoWaveStruct.base_graph (0D_NOT_type - Graph before curves extended to 1.5 periods.
      .def_property("base_graph", &TaoWaveStruct::base_graph, &TaoWaveStruct::set_base_graph)
      // TaoWaveStruct.region (0D_PTR_type - Where the wave plot is
      .def_property("region", &TaoWaveStruct::region, &TaoWaveStruct::set_region)
      // TaoWaveStruct.d1_dat (0D_PTR_type - D1 data for analysis
      .def_property("d1_dat", &TaoWaveStruct::d1_dat, &TaoWaveStruct::set_d1_dat)

      .def("__repr__", [](const TaoWaveStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoWaveStruct &self) {
            return TaoWaveStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoWaveStruct &self, py::dict &memo) { return TaoWaveStruct(self); }
      )

      ;

  // 1D TaoWaveStruct arrays are not used in structs/routines
  // 2D TaoWaveStruct arrays are not used in structs/routines
  // 3D TaoWaveStruct arrays are not used in structs/routines
}

// =============================================================================
// test_sub_struct
void init_test_sub_struct(py::module &m, py::class_<TestSubStruct> &cls) {
  cls.def(py::init<>())
      // TestSubStruct.sr (0D_NOT_type -
      .def_property("sr", &TestSubStruct::sr, &TestSubStruct::set_sr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return TestSubStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const TestSubStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TestSubStruct &self) {
            return TestSubStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TestSubStruct &self, py::dict &memo) { return TestSubStruct(self); }
      )

      ;

  bind_FTypeArrayND<TestSubStructArray1D>(m, "TestSubStructArray1D");
  bind_FTypeAlloc1D<TestSubStructAlloc1D>(m, "TestSubStructAlloc1D");
  bind_FTypeArrayND<TestSubStructArray2D>(m, "TestSubStructArray2D");
  bind_FTypeArrayND<TestSubStructArray3D>(m, "TestSubStructArray3D");
}

// =============================================================================
// test_sub_sub_struct
void init_test_sub_sub_struct(py::module &m, py::class_<TestSubSubStruct> &cls) {
  cls.def(py::init<>())
      // TestSubSubStruct.aaa (0D_NOT_integer8 -
      .def_property("aaa", &TestSubSubStruct::aaa, &TestSubSubStruct::set_aaa)
      // TestSubSubStruct.bbb (0D_NOT_integer -
      .def_property("bbb", &TestSubSubStruct::bbb, &TestSubSubStruct::set_bbb)
      // TestSubSubStruct.file (0D_NOT_character -
      .def_property("file", &TestSubSubStruct::file, &TestSubSubStruct::set_file)
      // TestSubSubStruct.t_ref (0D_NOT_real - time reference value for computing the wake
      // amplitude. This is used to prevent value overflow with long trains.
      .def_property("t_ref", &TestSubSubStruct::t_ref, &TestSubSubStruct::set_t_ref)
      // TestSubSubStruct.freq_spread (0D_NOT_real - Random frequency spread of long range modes.
      .def_property(
          "freq_spread",
          &TestSubSubStruct::freq_spread,
          &TestSubSubStruct::set_freq_spread
      )

      .def("__repr__", [](const TestSubSubStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TestSubSubStruct &self) {
            return TestSubSubStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TestSubSubStruct &self, py::dict &memo) { return TestSubSubStruct(self); }
      )

      ;

  // 1D TestSubSubStruct arrays are not used in structs/routines
  // 2D TestSubSubStruct arrays are not used in structs/routines
  // 3D TestSubSubStruct arrays are not used in structs/routines
}