#include "pybmad/generated/structs_t.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

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
// tao_eval_node_struct
void init_tao_eval_node_struct(
    py::module& m,
    py::class_<TaoEvalNodeProxy>& cls) {
  cls.def(py::init<>())
      // TaoEvalNodeProxy.type (0D_NOT_integer -
      .def_property(
          "type", &TaoEvalNodeProxy::type, &TaoEvalNodeProxy::set_type)
      // TaoEvalNodeProxy.name (0D_NOT_character -
      .def_property(
          "name", &TaoEvalNodeProxy::name, &TaoEvalNodeProxy::set_name)
      // TaoEvalNodeProxy.scale (0D_NOT_real - Scale factor for ping data
      .def_property(
          "scale", &TaoEvalNodeProxy::scale, &TaoEvalNodeProxy::set_scale)
      // TaoEvalNodeProxy.value (1D_ALLOC_real -
      .def_property_readonly("value", &TaoEvalNodeProxy::value)
      // TaoEvalNodeProxy.info (1D_ALLOC_type -
      .def_property_readonly("info", &TaoEvalNodeProxy::info)
      // TaoEvalNodeProxy.node (1D_PTR_type - Child nodes for tree construction.
      .def_property_readonly("node", &TaoEvalNodeProxy::node)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoEvalNodeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoEvalNodeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoEvalNodeProxy& self) {
            return TaoEvalNodeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoEvalNodeProxy& self, py::dict& memo) {
            return TaoEvalNodeProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoEvalNodeProxyArray1D>(m, "TaoEvalNodeStructArray1D");
  bind_FTypeAlloc1D<TaoEvalNodeProxyAlloc1D>(m, "TaoEvalNodeStructAlloc1D");
  // 2D TaoEvalNodeProxy arrays are not used in structs/routines
  // 3D TaoEvalNodeProxy arrays are not used in structs/routines
}

// =============================================================================
// tao_expression_info_struct
void init_tao_expression_info_struct(
    py::module& m,
    py::class_<TaoExpressionInfoProxy>& cls) {
  cls.def(py::init<>())
      // TaoExpressionInfoProxy.good (0D_NOT_logical - Expression is valid.
      .def_property(
          "good",
          &TaoExpressionInfoProxy::good,
          &TaoExpressionInfoProxy::set_good)
      // TaoExpressionInfoProxy.ele (0D_PTR_type - Associated ele if it exists
      .def_property(
          "ele", &TaoExpressionInfoProxy::ele, &TaoExpressionInfoProxy::set_ele)
      // TaoExpressionInfoProxy.s (0D_NOT_real - Longitudinal position of expression.
      .def_property(
          "s", &TaoExpressionInfoProxy::s, &TaoExpressionInfoProxy::set_s)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return TaoExpressionInfoProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const TaoExpressionInfoProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoExpressionInfoProxy& self) {
            return TaoExpressionInfoProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const TaoExpressionInfoProxy& self, py::dict& memo) {
            return TaoExpressionInfoProxy(self);
          })

      ;

  bind_FTypeArrayND<TaoExpressionInfoProxyArray1D>(
      m, "TaoExpressionInfoStructArray1D");
  bind_FTypeAlloc1D<TaoExpressionInfoProxyAlloc1D>(
      m, "TaoExpressionInfoStructAlloc1D");
  // 2D TaoExpressionInfoProxy arrays are not used in structs/routines
  // 3D TaoExpressionInfoProxy arrays are not used in structs/routines
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