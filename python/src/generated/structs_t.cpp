#include "pybmad/generated/structs_t.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// target_point_struct
void init_target_point_struct(py::module &m, py::class_<TargetPointStruct> &cls) {
  cls.def(py::init<std::optional<std::vector<double>>>(), py::arg("r") = py::none())
      .def_property(
          "r",
          py::cpp_function(&TargetPointStruct::r, py::keep_alive<0, 1>()),
          &TargetPointStruct::set_r,
          "(x, y, z)"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TargetPointStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TargetPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TargetPointStruct &self, const TargetPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TargetPointStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TargetPointStructArray1D, TargetPointStructAlloc1D>(
      m,
      "TargetPointStructArray1D",
      "TargetPointStructAlloc1D"
  );
  // 2D TargetPointStruct arrays are not used in structs/routines
  // 3D TargetPointStruct arrays are not used in structs/routines
}

// =============================================================================
// taylor_struct
void init_taylor_struct(py::module &m, py::class_<TaylorStruct> &cls) {
  cls.def(py::init<std::optional<double>>(), py::arg("ref") = py::none())
      .def_property("ref", &TaylorStruct::ref, &TaylorStruct::set_ref)
      .def_property_readonly("term", py::cpp_function(&TaylorStruct::term, py::keep_alive<0, 1>()))
      .def_static(
          "new_array1d",
          [](int sz) { return TaylorStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaylorStruct &self, const TaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaylorStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaylorStructArray1D, TaylorStructAlloc1D>(
      m,
      "TaylorStructArray1D",
      "TaylorStructAlloc1D"
  );
  // 2D TaylorStruct arrays are not used in structs/routines
  // 3D TaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// taylor_term_struct
void init_taylor_term_struct(py::module &m, py::class_<TaylorTermStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("expn") = py::none()
  )
      .def_property("coef", &TaylorTermStruct::coef, &TaylorTermStruct::set_coef)
      .def_property(
          "expn",
          py::cpp_function(&TaylorTermStruct::expn, py::keep_alive<0, 1>()),
          &TaylorTermStruct::set_expn
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaylorTermStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaylorTermStruct &self, const TaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaylorTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaylorTermStructArray1D, TaylorTermStructAlloc1D>(
      m,
      "TaylorTermStructArray1D",
      "TaylorTermStructAlloc1D"
  );
  // 2D TaylorTermStruct arrays are not used in structs/routines
  // 3D TaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// track_point_struct
void init_track_point_struct(py::module &m, py::class_<TrackPointStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             optional_ref<const CoordStruct>,
             optional_ref<const EmFieldStruct>,
             optional_ref<const StrongBeamStruct>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("s_lab") = py::none(),
         py::arg("s_body") = py::none(),
         py::arg("orb") = py::none(),
         py::arg("field") = py::none(),
         py::arg("strong_beam") = py::none(),
         py::arg("vec0") = py::none(),
         py::arg("mat6") = py::none()
  )
      .def_property(
          "s_lab",
          &TrackPointStruct::s_lab,
          &TrackPointStruct::set_s_lab,
          "Longitudinal lab coord with respect to the upstream end."
      )
      .def_property(
          "s_body",
          &TrackPointStruct::s_body,
          &TrackPointStruct::set_s_body,
          "Longitudinal body coord with respect to the entrance end."
      )
      .def_property(
          "orb",
          py::cpp_function(&TrackPointStruct::orb, py::keep_alive<0, 1>()),
          &TrackPointStruct::set_orb,
          "Particle position in lab coords."
      )
      .def_property(
          "field",
          py::cpp_function(&TrackPointStruct::field, py::keep_alive<0, 1>()),
          &TrackPointStruct::set_field,
          "E&M fields in lab coordinates."
      )
      .def_property(
          "strong_beam",
          py::cpp_function(&TrackPointStruct::strong_beam, py::keep_alive<0, 1>()),
          &TrackPointStruct::set_strong_beam,
          "Strong beam info for beambeam element."
      )
      .def_property(
          "vec0",
          py::cpp_function(&TrackPointStruct::vec0, py::keep_alive<0, 1>()),
          &TrackPointStruct::set_vec0,
          "0th order part of xfer map from the beginning."
      )
      .def_property(
          "mat6",
          py::cpp_function(&TrackPointStruct::mat6, py::keep_alive<0, 1>()),
          &TrackPointStruct::set_mat6,
          "1st order part of xfer map (transfer matrix)."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TrackPointStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TrackPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TrackPointStruct &self, const TrackPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TrackPointStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TrackPointStructArray1D, TrackPointStructAlloc1D>(
      m,
      "TrackPointStructArray1D",
      "TrackPointStructAlloc1D"
  );
  // 2D TrackPointStruct arrays are not used in structs/routines
  // 3D TrackPointStruct arrays are not used in structs/routines
}

// =============================================================================
// track_struct
void init_track_struct(py::module &m, py::class_<TrackStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("ds_save") = py::none(),
         py::arg("n_pt") = py::none(),
         py::arg("n_bad") = py::none(),
         py::arg("n_ok") = py::none()
  )
      .def_property_readonly(
          "pt",
          py::cpp_function(&TrackStruct::pt, py::keep_alive<0, 1>()),
          "Array of track points indexed from 0."
      )
      .def_property(
          "ds_save",
          &TrackStruct::ds_save,
          &TrackStruct::set_ds_save,
          "Min distance between points. Not positive => Save at all points."
      )
      .def_property(
          "n_pt",
          &TrackStruct::n_pt,
          &TrackStruct::set_n_pt,
          "Track upper bound for %pt(0:) array. n_bad and n_ok are used by adaptive trackers to "
          "record the number of times the step length had to be shortened."
      )
      .def_property(
          "n_bad",
          &TrackStruct::n_bad,
          &TrackStruct::set_n_bad,
          "Number of 'bad' steps where the step length was shortened."
      )
      .def_property(
          "n_ok",
          &TrackStruct::n_ok,
          &TrackStruct::set_n_ok,
          "Number of 'good' steps where the step length was not shortened."
      )

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
      .def(
          "__eq__",
          [](const TrackStruct &self, const TrackStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TrackStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TrackStruct arrays are not used in structs/routines
  // 2D TrackStruct arrays are not used in structs/routines
  // 3D TrackStruct arrays are not used in structs/routines
}

// =============================================================================
// twiss_struct
void init_twiss_struct(py::module &m, py::class_<TwissStruct> &cls) {
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
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("beta") = py::none(),
         py::arg("alpha") = py::none(),
         py::arg("gamma") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("eta") = py::none(),
         py::arg("etap") = py::none(),
         py::arg("deta_ds") = py::none(),
         py::arg("sigma") = py::none(),
         py::arg("sigma_p") = py::none(),
         py::arg("emit") = py::none(),
         py::arg("norm_emit") = py::none(),
         py::arg("chrom") = py::none(),
         py::arg("dbeta_dpz") = py::none(),
         py::arg("dalpha_dpz") = py::none(),
         py::arg("deta_dpz") = py::none(),
         py::arg("detap_dpz") = py::none()
  )
      .def_property("beta", &TwissStruct::beta, &TwissStruct::set_beta)
      .def_property("alpha", &TwissStruct::alpha, &TwissStruct::set_alpha)
      .def_property("gamma", &TwissStruct::gamma, &TwissStruct::set_gamma)
      .def_property("phi", &TwissStruct::phi, &TwissStruct::set_phi)
      .def_property("eta", &TwissStruct::eta, &TwissStruct::set_eta)
      .def_property("etap", &TwissStruct::etap, &TwissStruct::set_etap)
      .def_property("deta_ds", &TwissStruct::deta_ds, &TwissStruct::set_deta_ds)
      .def_property("sigma", &TwissStruct::sigma, &TwissStruct::set_sigma)
      .def_property("sigma_p", &TwissStruct::sigma_p, &TwissStruct::set_sigma_p)
      .def_property("emit", &TwissStruct::emit, &TwissStruct::set_emit)
      .def_property("norm_emit", &TwissStruct::norm_emit, &TwissStruct::set_norm_emit)
      .def_property("chrom", &TwissStruct::chrom, &TwissStruct::set_chrom)
      .def_property("dbeta_dpz", &TwissStruct::dbeta_dpz, &TwissStruct::set_dbeta_dpz)
      .def_property("dalpha_dpz", &TwissStruct::dalpha_dpz, &TwissStruct::set_dalpha_dpz)
      .def_property("deta_dpz", &TwissStruct::deta_dpz, &TwissStruct::set_deta_dpz)
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
      .def(
          "__eq__",
          [](const TwissStruct &self, const TwissStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TwissStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TwissStruct arrays are not used in structs/routines
  // 2D TwissStruct arrays are not used in structs/routines
  // 3D TwissStruct arrays are not used in structs/routines
}

// =============================================================================
// tricubic_cmplx_coef_struct
void init_tricubic_cmplx_coef_struct(py::module &m, py::class_<TricubicCmplxCoefStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>,
             std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("i_box") = py::none()
  )
      .def_property(
          "coef",
          py::cpp_function(&TricubicCmplxCoefStruct::coef, py::keep_alive<0, 1>()),
          &TricubicCmplxCoefStruct::set_coef,
          "Coefs"
      )
      .def_property(
          "i_box",
          py::cpp_function(&TricubicCmplxCoefStruct::i_box, py::keep_alive<0, 1>()),
          &TricubicCmplxCoefStruct::set_i_box,
          "index at lower box corner."
      )

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
      .def(
          "__eq__",
          [](const TricubicCmplxCoefStruct &self, const TricubicCmplxCoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TricubicCmplxCoefStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             optional_ref<const BeamStruct>,
             optional_ref<const BeamInitStruct>,
             optional_ref<const BeamInitStruct>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("beam_at_start") = py::none(),
         py::arg("beam_init") = py::none(),
         py::arg("beam_init_used") = py::none(),
         py::arg("init_starting_distribution") = py::none(),
         py::arg("track_start") = py::none(),
         py::arg("track_end") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_track_start") = py::none(),
         py::arg("ix_track_end") = py::none()
  )
      .def_property(
          "beam_at_start",
          py::cpp_function(&TaoBeamBranchStruct::beam_at_start, py::keep_alive<0, 1>()),
          &TaoBeamBranchStruct::set_beam_at_start,
          "Initial beam"
      )
      .def_property(
          "beam_init",
          py::cpp_function(&TaoBeamBranchStruct::beam_init, py::keep_alive<0, 1>()),
          &TaoBeamBranchStruct::set_beam_init,
          "User set beam distrubution at track start."
      )
      .def_property(
          "beam_init_used",
          py::cpp_function(&TaoBeamBranchStruct::beam_init_used, py::keep_alive<0, 1>()),
          &TaoBeamBranchStruct::set_beam_init_used,
          "beam distribution with emit values set."
      )
      .def_property(
          "init_starting_distribution",
          &TaoBeamBranchStruct::init_starting_distribution,
          &TaoBeamBranchStruct::set_init_starting_distribution,
          "Init beam"
      )
      .def_property(
          "track_start",
          &TaoBeamBranchStruct::track_start,
          &TaoBeamBranchStruct::set_track_start,
          "Tracking start element."
      )
      .def_property(
          "track_end",
          &TaoBeamBranchStruct::track_end,
          &TaoBeamBranchStruct::set_track_end
      )
      .def_property(
          "ix_branch",
          &TaoBeamBranchStruct::ix_branch,
          &TaoBeamBranchStruct::set_ix_branch,
          "Branch tracked. If track_start or track_end is a lord, ix_track_start/end index will be "
          "a index of slave."
      )
      .def_property(
          "ix_track_start",
          &TaoBeamBranchStruct::ix_track_start,
          &TaoBeamBranchStruct::set_ix_track_start,
          "Element track start index."
      )
      .def_property(
          "ix_track_end",
          &TaoBeamBranchStruct::ix_track_end,
          &TaoBeamBranchStruct::set_ix_track_end,
          "Element track end index"
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
      .def(
          "__eq__",
          [](const TaoBeamBranchStruct &self, const TaoBeamBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBeamBranchStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoBeamBranchStruct arrays are not used in structs/routines
  // 2D TaoBeamBranchStruct arrays are not used in structs/routines
  // 3D TaoBeamBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_beam_uni_struct
void init_tao_beam_uni_struct(py::module &m, py::class_<TaoBeamUniStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("saved_at") = py::none(),
         py::arg("dump_file") = py::none(),
         py::arg("dump_at") = py::none(),
         py::arg("track_beam_in_universe") = py::none(),
         py::arg("always_reinit") = py::none()
  )
      .def_property("saved_at", &TaoBeamUniStruct::saved_at, &TaoBeamUniStruct::set_saved_at)
      .def_property("dump_file", &TaoBeamUniStruct::dump_file, &TaoBeamUniStruct::set_dump_file)
      .def_property("dump_at", &TaoBeamUniStruct::dump_at, &TaoBeamUniStruct::set_dump_at)
      .def_property(
          "track_beam_in_universe",
          &TaoBeamUniStruct::track_beam_in_universe,
          &TaoBeamUniStruct::set_track_beam_in_universe,
          "Beam tracking enabled in this universe?"
      )
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
      .def(
          "__eq__",
          [](const TaoBeamUniStruct &self, const TaoBeamUniStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBeamUniStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         py::arg("theta") = py::none(),
         py::arg("x_offset") = py::none(),
         py::arg("z_offset") = py::none()
  )
      .def_property(
          "theta",
          &TaoBuildingWallOrientationStruct::theta,
          &TaoBuildingWallOrientationStruct::set_theta
      )
      .def_property(
          "x_offset",
          &TaoBuildingWallOrientationStruct::x_offset,
          &TaoBuildingWallOrientationStruct::set_x_offset
      )
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
      .def(
          "__eq__",
          [](const TaoBuildingWallOrientationStruct &self,
             const TaoBuildingWallOrientationStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBuildingWallOrientationStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("z") = py::none(),
         py::arg("x") = py::none(),
         py::arg("radius") = py::none(),
         py::arg("z_center") = py::none(),
         py::arg("x_center") = py::none()
  )
      .def_property(
          "z",
          &TaoBuildingWallPointStruct::z,
          &TaoBuildingWallPointStruct::set_z,
          "Global floor position"
      )
      .def_property(
          "x",
          &TaoBuildingWallPointStruct::x,
          &TaoBuildingWallPointStruct::set_x,
          "Global floor position"
      )
      .def_property(
          "radius",
          &TaoBuildingWallPointStruct::radius,
          &TaoBuildingWallPointStruct::set_radius,
          "Arc radius. +r -> CW rotation, same as bends."
      )
      .def_property(
          "z_center",
          &TaoBuildingWallPointStruct::z_center,
          &TaoBuildingWallPointStruct::set_z_center,
          "Arc center."
      )
      .def_property(
          "x_center",
          &TaoBuildingWallPointStruct::x_center,
          &TaoBuildingWallPointStruct::set_x_center,
          "Arc center."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoBuildingWallPointStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoBuildingWallPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoBuildingWallPointStruct &self, const TaoBuildingWallPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBuildingWallPointStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoBuildingWallPointStructArray1D, TaoBuildingWallPointStructAlloc1D>(
      m,
      "TaoBuildingWallPointStructArray1D",
      "TaoBuildingWallPointStructAlloc1D"
  );
  // 2D TaoBuildingWallPointStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallPointStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_section_struct
void init_tao_building_wall_section_struct(
    py::module &m,
    py::class_<TaoBuildingWallSectionStruct> &cls
) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<std::string>>(),
         py::arg("name") = py::none(),
         py::arg("constraint") = py::none()
  )
      .def_property(
          "name",
          &TaoBuildingWallSectionStruct::name,
          &TaoBuildingWallSectionStruct::set_name
      )
      .def_property(
          "constraint",
          &TaoBuildingWallSectionStruct::constraint,
          &TaoBuildingWallSectionStruct::set_constraint,
          "'left_side' or 'right_side' constraint."
      )
      .def_property_readonly(
          "point",
          py::cpp_function(&TaoBuildingWallSectionStruct::point, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoBuildingWallSectionStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoBuildingWallSectionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoBuildingWallSectionStruct &self, const TaoBuildingWallSectionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBuildingWallSectionStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoBuildingWallSectionStructArray1D, TaoBuildingWallSectionStructAlloc1D>(
      m,
      "TaoBuildingWallSectionStructArray1D",
      "TaoBuildingWallSectionStructAlloc1D"
  );
  // 2D TaoBuildingWallSectionStruct arrays are not used in structs/routines
  // 3D TaoBuildingWallSectionStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_building_wall_struct
void init_tao_building_wall_struct(py::module &m, py::class_<TaoBuildingWallStruct> &cls) {
  cls.def(
         py::init<optional_ref<const TaoBuildingWallOrientationStruct>>(),
         py::arg("orientation") = py::none()
  )
      .def_property(
          "orientation",
          py::cpp_function(&TaoBuildingWallStruct::orientation, py::keep_alive<0, 1>()),
          &TaoBuildingWallStruct::set_orientation
      )
      .def_property_readonly(
          "section",
          py::cpp_function(&TaoBuildingWallStruct::section, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const TaoBuildingWallStruct &self, const TaoBuildingWallStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoBuildingWallStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>>(),
         py::arg("cmd") = py::none(),
         py::arg("ix") = py::none()
  )
      .def_property("cmd", &TaoCmdHistoryStruct::cmd, &TaoCmdHistoryStruct::set_cmd, "The command")
      .def_property(
          "ix",
          &TaoCmdHistoryStruct::ix,
          &TaoCmdHistoryStruct::set_ix,
          "Command index (1st command has ix = 1, etc.) Note: Commands from command files will be "
          "assigned an index."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCmdHistoryStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCmdHistoryStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoCmdHistoryStruct &self, const TaoCmdHistoryStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCmdHistoryStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoCmdHistoryStructArray1D, TaoCmdHistoryStructAlloc1D>(
      m,
      "TaoCmdHistoryStructArray1D",
      "TaoCmdHistoryStructAlloc1D"
  );
  // 2D TaoCmdHistoryStruct arrays are not used in structs/routines
  // 3D TaoCmdHistoryStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_common_struct
void init_tao_common_struct(py::module &m, py::class_<TaoCommonStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
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
             std::optional<std::vector<bool>>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>>(),
         py::arg("covar") = py::none(),
         py::arg("alpha") = py::none(),
         py::arg("dummy_target") = py::none(),
         py::arg("n_alias") = py::none(),
         py::arg("cmd_file_level") = py::none(),
         py::arg("ix_key_bank") = py::none(),
         py::arg("ix_history") = py::none(),
         py::arg("n_history") = py::none(),
         py::arg("lev_loop") = py::none(),
         py::arg("n_err_messages_printed") = py::none(),
         py::arg("n_universes") = py::none(),
         py::arg("ix_beam_track_active_element") = py::none(),
         py::arg("cmd_file_paused") = py::none(),
         py::arg("use_cmd_here") = py::none(),
         py::arg("cmd_from_cmd_file") = py::none(),
         py::arg("use_saved_beam_in_tracking") = py::none(),
         py::arg("single_mode") = py::none(),
         py::arg("combine_consecutive_elements_of_like_name") = py::none(),
         py::arg("have_tracked_beam") = py::none(),
         py::arg("init_plot_needed") = py::none(),
         py::arg("init_beam") = py::none(),
         py::arg("init_var") = py::none(),
         py::arg("init_read_lat_info") = py::none(),
         py::arg("optimizer_running") = py::none(),
         py::arg("have_datums_using_expressions") = py::none(),
         py::arg("print_to_terminal") = py::none(),
         py::arg("lattice_calc_done") = py::none(),
         py::arg("add_measurement_noise") = py::none(),
         py::arg("is_err_message_printed") = py::none(),
         py::arg("command_arg_has_been_executed") = py::none(),
         py::arg("all_merit_weights_positive") = py::none(),
         py::arg("multi_turn_orbit_is_plotted") = py::none(),
         py::arg("force_chrom_calc") = py::none(),
         py::arg("force_rad_int_calc") = py::none(),
         py::arg("rad_int_ri_calc_on") = py::none(),
         py::arg("rad_int_6d_calc_on") = py::none(),
         py::arg("single_mode_buffer") = py::none(),
         py::arg("cmd") = py::none()
  )
      .def_property_readonly(
          "plot_place_buffer",
          py::cpp_function(&TaoCommonStruct::plot_place_buffer, py::keep_alive<0, 1>()),
          "Used when %external_plotting is on."
      )
      .def_property(
          "covar",
          py::cpp_function(&TaoCommonStruct::covar, py::keep_alive<0, 1>()),
          &TaoCommonStruct::set_covar
      )
      .def_property(
          "alpha",
          py::cpp_function(&TaoCommonStruct::alpha, py::keep_alive<0, 1>()),
          &TaoCommonStruct::set_alpha
      )
      .def_property(
          "dummy_target",
          &TaoCommonStruct::dummy_target,
          &TaoCommonStruct::set_dummy_target,
          "Dummy varaible"
      )
      .def_property("n_alias", &TaoCommonStruct::n_alias, &TaoCommonStruct::set_n_alias)
      .def_property(
          "cmd_file_level",
          &TaoCommonStruct::cmd_file_level,
          &TaoCommonStruct::set_cmd_file_level,
          "For nested command files. 0 -> no command file."
      )
      .def_property(
          "ix_key_bank",
          &TaoCommonStruct::ix_key_bank,
          &TaoCommonStruct::set_ix_key_bank,
          "For single mode."
      )
      .def_property(
          "ix_history",
          &TaoCommonStruct::ix_history,
          &TaoCommonStruct::set_ix_history,
          "Index to latest command in the history circular buffer."
      )
      .def_property(
          "n_history",
          &TaoCommonStruct::n_history,
          &TaoCommonStruct::set_n_history,
          "Number of commands issued from beginning of starting Tao."
      )
      .def_property(
          "lev_loop",
          &TaoCommonStruct::lev_loop,
          &TaoCommonStruct::set_lev_loop,
          "in do loop nest level"
      )
      .def_property(
          "n_err_messages_printed",
          &TaoCommonStruct::n_err_messages_printed,
          &TaoCommonStruct::set_n_err_messages_printed,
          "Used by tao_set_invalid to limit number of messages."
      )
      .def_property("n_universes", &TaoCommonStruct::n_universes, &TaoCommonStruct::set_n_universes)
      .def_property(
          "ix_beam_track_active_element",
          &TaoCommonStruct::ix_beam_track_active_element,
          &TaoCommonStruct::set_ix_beam_track_active_element,
          "Element being tracked through `tao_beam_track`."
      )
      .def_property(
          "cmd_file_paused",
          &TaoCommonStruct::cmd_file_paused,
          &TaoCommonStruct::set_cmd_file_paused
      )
      .def_property(
          "use_cmd_here",
          &TaoCommonStruct::use_cmd_here,
          &TaoCommonStruct::set_use_cmd_here,
          "Used for commands recalled from the cmd history stack"
      )
      .def_property(
          "cmd_from_cmd_file",
          &TaoCommonStruct::cmd_from_cmd_file,
          &TaoCommonStruct::set_cmd_from_cmd_file,
          "was command from a command file?"
      )
      .def_property(
          "use_saved_beam_in_tracking",
          &TaoCommonStruct::use_saved_beam_in_tracking,
          &TaoCommonStruct::set_use_saved_beam_in_tracking
      )
      .def_property("single_mode", &TaoCommonStruct::single_mode, &TaoCommonStruct::set_single_mode)
      .def_property(
          "combine_consecutive_elements_of_like_name",
          &TaoCommonStruct::combine_consecutive_elements_of_like_name,
          &TaoCommonStruct::set_combine_consecutive_elements_of_like_name
      )
      .def_property(
          "have_tracked_beam",
          &TaoCommonStruct::have_tracked_beam,
          &TaoCommonStruct::set_have_tracked_beam,
          "Used to catch error when beam plotting without having tracked a beam."
      )
      .def_property(
          "init_plot_needed",
          &TaoCommonStruct::init_plot_needed,
          &TaoCommonStruct::set_init_plot_needed,
          "reinitialize plotting?"
      )
      .def_property(
          "init_beam",
          &TaoCommonStruct::init_beam,
          &TaoCommonStruct::set_init_beam,
          "Used by custom programs to control Tao init"
      )
      .def_property(
          "init_var",
          &TaoCommonStruct::init_var,
          &TaoCommonStruct::set_init_var,
          "Used by custom programs to control Tao init"
      )
      .def_property(
          "init_read_lat_info",
          &TaoCommonStruct::init_read_lat_info,
          &TaoCommonStruct::set_init_read_lat_info,
          "Used by custom programs to control Tao init"
      )
      .def_property(
          "optimizer_running",
          &TaoCommonStruct::optimizer_running,
          &TaoCommonStruct::set_optimizer_running
      )
      .def_property(
          "have_datums_using_expressions",
          &TaoCommonStruct::have_datums_using_expressions,
          &TaoCommonStruct::set_have_datums_using_expressions
      )
      .def_property(
          "print_to_terminal",
          &TaoCommonStruct::print_to_terminal,
          &TaoCommonStruct::set_print_to_terminal,
          "Print command prompt to the terminal? For use with GUIs."
      )
      .def_property(
          "lattice_calc_done",
          &TaoCommonStruct::lattice_calc_done,
          &TaoCommonStruct::set_lattice_calc_done,
          "Used by GUI for deciding when to refresh."
      )
      .def_property(
          "add_measurement_noise",
          &TaoCommonStruct::add_measurement_noise,
          &TaoCommonStruct::set_add_measurement_noise,
          "Turn off to take data derivatives."
      )
      .def_property(
          "is_err_message_printed",
          py::cpp_function(&TaoCommonStruct::is_err_message_printed, py::keep_alive<0, 1>()),
          &TaoCommonStruct::set_is_err_message_printed,
          "Used by tao_set_invalid"
      )
      .def_property(
          "command_arg_has_been_executed",
          &TaoCommonStruct::command_arg_has_been_executed,
          &TaoCommonStruct::set_command_arg_has_been_executed,
          "Has the -command command line argument been executed?"
      )
      .def_property(
          "all_merit_weights_positive",
          &TaoCommonStruct::all_merit_weights_positive,
          &TaoCommonStruct::set_all_merit_weights_positive
      )
      .def_property(
          "multi_turn_orbit_is_plotted",
          &TaoCommonStruct::multi_turn_orbit_is_plotted,
          &TaoCommonStruct::set_multi_turn_orbit_is_plotted,
          "Is a multi_turn_orbit being plotted?"
      )
      .def_property(
          "force_chrom_calc",
          &TaoCommonStruct::force_chrom_calc,
          &TaoCommonStruct::set_force_chrom_calc,
          "Used by a routine to force a single chromaticity calculation."
      )
      .def_property(
          "force_rad_int_calc",
          &TaoCommonStruct::force_rad_int_calc,
          &TaoCommonStruct::set_force_rad_int_calc,
          "Used by a routine to force a single radiation integrals calculation"
      )
      .def_property(
          "rad_int_ri_calc_on",
          &TaoCommonStruct::rad_int_ri_calc_on,
          &TaoCommonStruct::set_rad_int_ri_calc_on,
          "'Classical' radiation integrals calculation on/off."
      )
      .def_property(
          "rad_int_6d_calc_on",
          &TaoCommonStruct::rad_int_6d_calc_on,
          &TaoCommonStruct::set_rad_int_6d_calc_on,
          "6D Radiation integrals calculation on/off."
      )
      .def_property_readonly(
          "valid_plot_who",
          py::cpp_function(&TaoCommonStruct::valid_plot_who, py::keep_alive<0, 1>()),
          "model, base, ref etc..."
      )
      .def_property(
          "single_mode_buffer",
          &TaoCommonStruct::single_mode_buffer,
          &TaoCommonStruct::set_single_mode_buffer
      )
      .def_property(
          "cmd",
          &TaoCommonStruct::cmd,
          &TaoCommonStruct::set_cmd,
          "Used for the cmd history"
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
      .def(
          "__eq__",
          [](const TaoCommonStruct &self, const TaoCommonStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCommonStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoCommonStruct arrays are not used in structs/routines
  // 2D TaoCommonStruct arrays are not used in structs/routines
  // 3D TaoCommonStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_color_struct
void init_tao_curve_color_struct(py::module &m, py::class_<TaoCurveColorStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>>(),
         py::arg("data_type") = py::none(),
         py::arg("is_on") = py::none(),
         py::arg("min") = py::none(),
         py::arg("max") = py::none(),
         py::arg("autoscale") = py::none()
  )
      .def_property(
          "data_type",
          &TaoCurveColorStruct::data_type,
          &TaoCurveColorStruct::set_data_type,
          "Datum type to use for z-axis."
      )
      .def_property("is_on", &TaoCurveColorStruct::is_on, &TaoCurveColorStruct::set_is_on, "On/Off")
      .def_property(
          "min",
          &TaoCurveColorStruct::min,
          &TaoCurveColorStruct::set_min,
          "Min and max values for mapping z-axis to color."
      )
      .def_property(
          "max",
          &TaoCurveColorStruct::max,
          &TaoCurveColorStruct::set_max,
          "Min and max values for mapping z-axis to color."
      )
      .def_property(
          "autoscale",
          &TaoCurveColorStruct::autoscale,
          &TaoCurveColorStruct::set_autoscale,
          "Set %min, %max automatically to the limits of %data_type"
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
      .def(
          "__eq__",
          [](const TaoCurveColorStruct &self, const TaoCurveColorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCurveColorStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoCurveColorStruct arrays are not used in structs/routines
  // 2D TaoCurveColorStruct arrays are not used in structs/routines
  // 3D TaoCurveColorStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_orbit_struct
void init_tao_curve_orbit_struct(py::module &m, py::class_<TaoCurveOrbitStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("t") = py::none()
  )
      .def_property("x", &TaoCurveOrbitStruct::x, &TaoCurveOrbitStruct::set_x, "Transverse offset")
      .def_property("y", &TaoCurveOrbitStruct::y, &TaoCurveOrbitStruct::set_y, "Transverse offset")
      .def_property("t", &TaoCurveOrbitStruct::t, &TaoCurveOrbitStruct::set_t, "Time")

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
      .def(
          "__eq__",
          [](const TaoCurveOrbitStruct &self, const TaoCurveOrbitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCurveOrbitStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoCurveOrbitStruct arrays are not used in structs/routines
  // 2D TaoCurveOrbitStruct arrays are not used in structs/routines
  // 3D TaoCurveOrbitStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_struct
void init_tao_curve_struct(py::module &m, py::class_<TaoCurveStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             optional_ref<const TaoGraphStruct>,
             optional_ref<const TaoHistogramStruct>,
             optional_ref<const TaoCurveColorStruct>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<int>>,
             std::optional<double>,
             optional_ref<const QpLineStruct>,
             optional_ref<const QpSymbolStruct>,
             optional_ref<const TaoCurveOrbitStruct>,
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
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("data_source") = py::none(),
         py::arg("data_index") = py::none(),
         py::arg("data_type_x") = py::none(),
         py::arg("data_type") = py::none(),
         py::arg("ele_ref_name") = py::none(),
         py::arg("legend_text") = py::none(),
         py::arg("message_text") = py::none(),
         py::arg("component") = py::none(),
         py::arg("why_invalid") = py::none(),
         py::arg("g") = py::none(),
         py::arg("hist") = py::none(),
         py::arg("z_color") = py::none(),
         py::arg("x_line") = py::none(),
         py::arg("y_line") = py::none(),
         py::arg("y2_line") = py::none(),
         py::arg("ix_line") = py::none(),
         py::arg("x_symb") = py::none(),
         py::arg("y_symb") = py::none(),
         py::arg("z_symb") = py::none(),
         py::arg("err_symb") = py::none(),
         py::arg("symb_size") = py::none(),
         py::arg("ix_symb") = py::none(),
         py::arg("y_axis_scale_factor") = py::none(),
         py::arg("line") = py::none(),
         py::arg("symbol") = py::none(),
         py::arg("orbit") = py::none(),
         py::arg("ix_universe") = py::none(),
         py::arg("symbol_every") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_bunch") = py::none(),
         py::arg("n_turn") = py::none(),
         py::arg("use_y2") = py::none(),
         py::arg("draw_line") = py::none(),
         py::arg("draw_symbols") = py::none(),
         py::arg("draw_symbol_index") = py::none(),
         py::arg("draw_error_bars") = py::none(),
         py::arg("smooth_line_calc") = py::none(),
         py::arg("valid") = py::none()
  )
      .def_property(
          "name",
          &TaoCurveStruct::name,
          &TaoCurveStruct::set_name,
          "Name identifying the curve."
      )
      .def_property(
          "data_source",
          &TaoCurveStruct::data_source,
          &TaoCurveStruct::set_data_source,
          "'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit'"
      )
      .def_property(
          "data_index",
          &TaoCurveStruct::data_index,
          &TaoCurveStruct::set_data_index,
          "Used for calculating %ix_symb(:)."
      )
      .def_property(
          "data_type_x",
          &TaoCurveStruct::data_type_x,
          &TaoCurveStruct::set_data_type_x,
          "Used for data slices and phase space plots."
      )
      .def_property(
          "data_type",
          &TaoCurveStruct::data_type,
          &TaoCurveStruct::set_data_type,
          "'orbit.x', etc."
      )
      .def_property(
          "ele_ref_name",
          &TaoCurveStruct::ele_ref_name,
          &TaoCurveStruct::set_ele_ref_name,
          "Reference element."
      )
      .def_property(
          "legend_text",
          &TaoCurveStruct::legend_text,
          &TaoCurveStruct::set_legend_text,
          "String to draw in a curve legend."
      )
      .def_property(
          "message_text",
          &TaoCurveStruct::message_text,
          &TaoCurveStruct::set_message_text,
          "Informational message to draw with graph."
      )
      .def_property(
          "component",
          &TaoCurveStruct::component,
          &TaoCurveStruct::set_component,
          "Who to plot. Eg: 'meas - design'"
      )
      .def_property(
          "why_invalid",
          &TaoCurveStruct::why_invalid,
          &TaoCurveStruct::set_why_invalid,
          "Informative string to print."
      )
      .def_property(
          "g",
          py::cpp_function(&TaoCurveStruct::g, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_g,
          "pointer to parent graph"
      )
      .def_property(
          "hist",
          py::cpp_function(&TaoCurveStruct::hist, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_hist
      )
      .def_property(
          "z_color",
          py::cpp_function(&TaoCurveStruct::z_color, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_z_color
      )
      .def_property(
          "x_line",
          py::cpp_function(&TaoCurveStruct::x_line, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_x_line,
          "Coords for drawing a curve"
      )
      .def_property(
          "y_line",
          py::cpp_function(&TaoCurveStruct::y_line, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_y_line
      )
      .def_property(
          "y2_line",
          py::cpp_function(&TaoCurveStruct::y2_line, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_y2_line,
          "Second array needed for beam chamber curve."
      )
      .def_property(
          "ix_line",
          py::cpp_function(&TaoCurveStruct::ix_line, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_ix_line,
          "Used by wave and aperture curves."
      )
      .def_property(
          "x_symb",
          py::cpp_function(&TaoCurveStruct::x_symb, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_x_symb,
          "Coords for drawing the symbols"
      )
      .def_property(
          "y_symb",
          py::cpp_function(&TaoCurveStruct::y_symb, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_y_symb
      )
      .def_property(
          "z_symb",
          py::cpp_function(&TaoCurveStruct::z_symb, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_z_symb,
          "Symbol color"
      )
      .def_property(
          "err_symb",
          py::cpp_function(&TaoCurveStruct::err_symb, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_err_symb,
          "Error bars"
      )
      .def_property(
          "symb_size",
          py::cpp_function(&TaoCurveStruct::symb_size, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_symb_size,
          "Symbol size. Used with symbol_size_scale."
      )
      .def_property(
          "ix_symb",
          py::cpp_function(&TaoCurveStruct::ix_symb, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_ix_symb,
          "Corresponding index in d1_data%d(:) array."
      )
      .def_property(
          "y_axis_scale_factor",
          &TaoCurveStruct::y_axis_scale_factor,
          &TaoCurveStruct::set_y_axis_scale_factor,
          "y-axis conversion from internal to plotting units."
      )
      .def_property(
          "line",
          py::cpp_function(&TaoCurveStruct::line, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_line,
          "Line attributes"
      )
      .def_property(
          "symbol",
          py::cpp_function(&TaoCurveStruct::symbol, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_symbol,
          "Symbol attributes"
      )
      .def_property(
          "orbit",
          py::cpp_function(&TaoCurveStruct::orbit, py::keep_alive<0, 1>()),
          &TaoCurveStruct::set_orbit,
          "Used for E/B field plotting."
      )
      .def_property(
          "ix_universe",
          &TaoCurveStruct::ix_universe,
          &TaoCurveStruct::set_ix_universe,
          "Universe where data is. -1 => use s%global%default_universe"
      )
      .def_property(
          "symbol_every",
          &TaoCurveStruct::symbol_every,
          &TaoCurveStruct::set_symbol_every,
          "Symbol every how many points."
      )
      .def_property("ix_branch", &TaoCurveStruct::ix_branch, &TaoCurveStruct::set_ix_branch)
      .def_property(
          "ix_bunch",
          &TaoCurveStruct::ix_bunch,
          &TaoCurveStruct::set_ix_bunch,
          "Bunch to plot."
      )
      .def_property(
          "n_turn",
          &TaoCurveStruct::n_turn,
          &TaoCurveStruct::set_n_turn,
          "Used for multi_turn_orbit plotting"
      )
      .def_property("use_y2", &TaoCurveStruct::use_y2, &TaoCurveStruct::set_use_y2, "Use y2 axis?")
      .def_property(
          "draw_line",
          &TaoCurveStruct::draw_line,
          &TaoCurveStruct::set_draw_line,
          "Draw a line through the data points?"
      )
      .def_property(
          "draw_symbols",
          &TaoCurveStruct::draw_symbols,
          &TaoCurveStruct::set_draw_symbols,
          "Draw a symbol at the data points?"
      )
      .def_property(
          "draw_symbol_index",
          &TaoCurveStruct::draw_symbol_index,
          &TaoCurveStruct::set_draw_symbol_index,
          "Draw the symbol index number curve%ix_symb?"
      )
      .def_property(
          "draw_error_bars",
          &TaoCurveStruct::draw_error_bars,
          &TaoCurveStruct::set_draw_error_bars,
          "Draw error bars based upon data%error_rms if drawing data? !! logical :: draw_rms = "
          ".false.          ! Show mean and RMS values with legend?"
      )
      .def_property(
          "smooth_line_calc",
          &TaoCurveStruct::smooth_line_calc,
          &TaoCurveStruct::set_smooth_line_calc,
          "Calculate data between element edge points?"
      )
      .def_property("valid", &TaoCurveStruct::valid, &TaoCurveStruct::set_valid, "valid data?")
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCurveStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCurveStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoCurveStruct &self, const TaoCurveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCurveStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoCurveStructArray1D, TaoCurveStructAlloc1D>(
      m,
      "TaoCurveStructArray1D",
      "TaoCurveStructAlloc1D"
  );
  // 2D TaoCurveStruct arrays are not used in structs/routines
  // 3D TaoCurveStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d1_data_struct
void init_tao_d1_data_struct(py::module &m, py::class_<TaoD1DataStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, optional_ref<const TaoD2DataStruct>>(),
         py::arg("name") = py::none(),
         py::arg("d2") = py::none()
  )
      .def_property("name", &TaoD1DataStruct::name, &TaoD1DataStruct::set_name, "Eg: 'x', etc.")
      .def_property(
          "d2",
          py::cpp_function(&TaoD1DataStruct::d2, py::keep_alive<0, 1>()),
          &TaoD1DataStruct::set_d2,
          "ptr to parent d2_data"
      )
      .def_property_readonly(
          "d",
          py::cpp_function(&TaoD1DataStruct::d, py::keep_alive<0, 1>()),
          "Pointer to the appropriate section in u%data"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD1DataStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD1DataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoD1DataStruct &self, const TaoD1DataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoD1DataStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoD1DataStructArray1D, TaoD1DataStructAlloc1D>(
      m,
      "TaoD1DataStructArray1D",
      "TaoD1DataStructAlloc1D"
  );
  // 2D TaoD1DataStruct arrays are not used in structs/routines
  // 3D TaoD1DataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d2_data_struct
void init_tao_d2_data_struct(py::module &m, py::class_<TaoD2DataStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("data_file_name") = py::none(),
         py::arg("ref_file_name") = py::none(),
         py::arg("data_date") = py::none(),
         py::arg("ref_date") = py::none(),
         py::arg("ix_universe") = py::none(),
         py::arg("ix_d2_data") = py::none(),
         py::arg("ix_ref") = py::none(),
         py::arg("data_read_in") = py::none(),
         py::arg("ref_read_in") = py::none()
  )
      .def_property(
          "name",
          &TaoD2DataStruct::name,
          &TaoD2DataStruct::set_name,
          "Name to be used with commands."
      )
      .def_property(
          "data_file_name",
          &TaoD2DataStruct::data_file_name,
          &TaoD2DataStruct::set_data_file_name,
          "Data file name ."
      )
      .def_property(
          "ref_file_name",
          &TaoD2DataStruct::ref_file_name,
          &TaoD2DataStruct::set_ref_file_name,
          "Reference file name."
      )
      .def_property(
          "data_date",
          &TaoD2DataStruct::data_date,
          &TaoD2DataStruct::set_data_date,
          "Data measurement date."
      )
      .def_property(
          "ref_date",
          &TaoD2DataStruct::ref_date,
          &TaoD2DataStruct::set_ref_date,
          "Reference data measurement date."
      )
      .def_property_readonly(
          "descrip",
          py::cpp_function(&TaoD2DataStruct::descrip, py::keep_alive<0, 1>()),
          "Array for descriptive information."
      )
      .def_property_readonly(
          "d1",
          py::cpp_function(&TaoD2DataStruct::d1, py::keep_alive<0, 1>()),
          "Points to children"
      )
      .def_property(
          "ix_universe",
          &TaoD2DataStruct::ix_universe,
          &TaoD2DataStruct::set_ix_universe,
          "Index of universe this is in."
      )
      .def_property(
          "ix_d2_data",
          &TaoD2DataStruct::ix_d2_data,
          &TaoD2DataStruct::set_ix_d2_data,
          "Index in u%d2_data(:) array."
      )
      .def_property(
          "ix_ref",
          &TaoD2DataStruct::ix_ref,
          &TaoD2DataStruct::set_ix_ref,
          "Index of the reference data set."
      )
      .def_property(
          "data_read_in",
          &TaoD2DataStruct::data_read_in,
          &TaoD2DataStruct::set_data_read_in,
          "A data set has been read in?"
      )
      .def_property(
          "ref_read_in",
          &TaoD2DataStruct::ref_read_in,
          &TaoD2DataStruct::set_ref_read_in,
          "A reference data set has been read in?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD2DataStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD2DataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoD2DataStruct &self, const TaoD2DataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoD2DataStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoD2DataStructArray1D, TaoD2DataStructAlloc1D>(
      m,
      "TaoD2DataStructArray1D",
      "TaoD2DataStructAlloc1D"
  );
  // 2D TaoD2DataStruct arrays are not used in structs/routines
  // 3D TaoD2DataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_data_struct
void init_tao_data_struct(py::module &m, py::class_<TaoDataStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
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
             optional_ref<const TaoSpinMapStruct>,
             optional_ref<const TaoD1DataStruct>>(),
         py::arg("ele_name") = py::none(),
         py::arg("ele_start_name") = py::none(),
         py::arg("ele_ref_name") = py::none(),
         py::arg("data_type") = py::none(),
         py::arg("merit_type") = py::none(),
         py::arg("id") = py::none(),
         py::arg("data_source") = py::none(),
         py::arg("why_invalid") = py::none(),
         py::arg("ix_uni") = py::none(),
         py::arg("ix_bunch") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_ele_start") = py::none(),
         py::arg("ix_ele_ref") = py::none(),
         py::arg("ix_ele_merit") = py::none(),
         py::arg("ix_d1") = py::none(),
         py::arg("ix_data") = py::none(),
         py::arg("ix_dModel") = py::none(),
         py::arg("eval_point") = py::none(),
         py::arg("meas_value") = py::none(),
         py::arg("ref_value") = py::none(),
         py::arg("model_value") = py::none(),
         py::arg("design_value") = py::none(),
         py::arg("old_value") = py::none(),
         py::arg("base_value") = py::none(),
         py::arg("error_rms") = py::none(),
         py::arg("delta_merit") = py::none(),
         py::arg("weight") = py::none(),
         py::arg("invalid_value") = py::none(),
         py::arg("merit") = py::none(),
         py::arg("s") = py::none(),
         py::arg("s_offset") = py::none(),
         py::arg("ref_s_offset") = py::none(),
         py::arg("err_message_printed") = py::none(),
         py::arg("exists") = py::none(),
         py::arg("good_model") = py::none(),
         py::arg("good_base") = py::none(),
         py::arg("good_design") = py::none(),
         py::arg("good_meas") = py::none(),
         py::arg("good_ref") = py::none(),
         py::arg("good_user") = py::none(),
         py::arg("good_opt") = py::none(),
         py::arg("good_plot") = py::none(),
         py::arg("useit_plot") = py::none(),
         py::arg("useit_opt") = py::none(),
         py::arg("spin_map") = py::none(),
         py::arg("d1") = py::none()
  )
      .def_property(
          "ele_name",
          &TaoDataStruct::ele_name,
          &TaoDataStruct::set_ele_name,
          "Name of the lattice element where datum is evaluated."
      )
      .def_property(
          "ele_start_name",
          &TaoDataStruct::ele_start_name,
          &TaoDataStruct::set_ele_start_name,
          "Name of starting lattice element when there is a range"
      )
      .def_property(
          "ele_ref_name",
          &TaoDataStruct::ele_ref_name,
          &TaoDataStruct::set_ele_ref_name,
          "Name of reference lattice element"
      )
      .def_property(
          "data_type",
          &TaoDataStruct::data_type,
          &TaoDataStruct::set_data_type,
          "Type of data: 'orbit.x', etc."
      )
      .def_property(
          "merit_type",
          &TaoDataStruct::merit_type,
          &TaoDataStruct::set_merit_type,
          "Type of constraint: 'target', 'max', 'min', etc."
      )
      .def_property(
          "id",
          &TaoDataStruct::id,
          &TaoDataStruct::set_id,
          "Used by Tao extension code. Not used by Tao directly."
      )
      .def_property(
          "data_source",
          &TaoDataStruct::data_source,
          &TaoDataStruct::set_data_source,
          "'lat', 'beam', 'data' or 'var'. Last two used for expressions."
      )
      .def_property(
          "why_invalid",
          &TaoDataStruct::why_invalid,
          &TaoDataStruct::set_why_invalid,
          "Informational string if there is a problem."
      )
      .def_property(
          "ix_uni",
          &TaoDataStruct::ix_uni,
          &TaoDataStruct::set_ix_uni,
          "Universe index of datum."
      )
      .def_property(
          "ix_bunch",
          &TaoDataStruct::ix_bunch,
          &TaoDataStruct::set_ix_bunch,
          "Bunch number to get the data from."
      )
      .def_property(
          "ix_branch",
          &TaoDataStruct::ix_branch,
          &TaoDataStruct::set_ix_branch,
          "Index of the associated lattice branch."
      )
      .def_property(
          "ix_ele",
          &TaoDataStruct::ix_ele,
          &TaoDataStruct::set_ix_ele,
          "Index of the lattice element corresponding to ele_name"
      )
      .def_property(
          "ix_ele_start",
          &TaoDataStruct::ix_ele_start,
          &TaoDataStruct::set_ix_ele_start,
          "Index of lattice elment when there is a range"
      )
      .def_property(
          "ix_ele_ref",
          &TaoDataStruct::ix_ele_ref,
          &TaoDataStruct::set_ix_ele_ref,
          "Index of lattice elment when there is a reference."
      )
      .def_property(
          "ix_ele_merit",
          &TaoDataStruct::ix_ele_merit,
          &TaoDataStruct::set_ix_ele_merit,
          "Index of lattice elment where merit is evaluated."
      )
      .def_property(
          "ix_d1",
          &TaoDataStruct::ix_d1,
          &TaoDataStruct::set_ix_d1,
          "Index number in u%d2_data(i)%d1_data(j)%d(:) array."
      )
      .def_property(
          "ix_data",
          &TaoDataStruct::ix_data,
          &TaoDataStruct::set_ix_data,
          "Index of this datum in the u%data(:) array of data_structs."
      )
      .def_property(
          "ix_dModel",
          &TaoDataStruct::ix_dModel,
          &TaoDataStruct::set_ix_dModel,
          "Row number in the dModel_dVar derivative matrix."
      )
      .def_property(
          "eval_point",
          &TaoDataStruct::eval_point,
          &TaoDataStruct::set_eval_point,
          "or anchor_center$, anchor_beginning$. Where to evaluate data relative to the element."
      )
      .def_property(
          "meas_value",
          &TaoDataStruct::meas_value,
          &TaoDataStruct::set_meas_value,
          "Measured datum value."
      )
      .def_property(
          "ref_value",
          &TaoDataStruct::ref_value,
          &TaoDataStruct::set_ref_value,
          "Measured datum value from the reference data set."
      )
      .def_property(
          "model_value",
          &TaoDataStruct::model_value,
          &TaoDataStruct::set_model_value,
          "Datum value as calculated from the model."
      )
      .def_property(
          "design_value",
          &TaoDataStruct::design_value,
          &TaoDataStruct::set_design_value,
          "What the datum value is in the design lattice."
      )
      .def_property(
          "old_value",
          &TaoDataStruct::old_value,
          &TaoDataStruct::set_old_value,
          "The model_value at some previous time."
      )
      .def_property(
          "base_value",
          &TaoDataStruct::base_value,
          &TaoDataStruct::set_base_value,
          "The value as calculated from the base model."
      )
      .def_property(
          "error_rms",
          &TaoDataStruct::error_rms,
          &TaoDataStruct::set_error_rms,
          "Measurement error RMS. Used in plotting."
      )
      .def_property(
          "delta_merit",
          &TaoDataStruct::delta_merit,
          &TaoDataStruct::set_delta_merit,
          "Diff used to calculate the merit function term."
      )
      .def_property(
          "weight",
          &TaoDataStruct::weight,
          &TaoDataStruct::set_weight,
          "Weight for the merit function term."
      )
      .def_property(
          "invalid_value",
          &TaoDataStruct::invalid_value,
          &TaoDataStruct::set_invalid_value,
          "Value used in merit calc if good_model = F (or possibly good_design & good_base)."
      )
      .def_property(
          "merit",
          &TaoDataStruct::merit,
          &TaoDataStruct::set_merit,
          "Merit function term value: weight * delta_merit^2"
      )
      .def_property("s", &TaoDataStruct::s, &TaoDataStruct::set_s, "longitudinal position of ele.")
      .def_property(
          "s_offset",
          &TaoDataStruct::s_offset,
          &TaoDataStruct::set_s_offset,
          "Offset of the evaluation point."
      )
      .def_property(
          "ref_s_offset",
          &TaoDataStruct::ref_s_offset,
          &TaoDataStruct::set_ref_s_offset,
          "Offset of the reference point. In development."
      )
      .def_property(
          "err_message_printed",
          &TaoDataStruct::err_message_printed,
          &TaoDataStruct::set_err_message_printed,
          "Used to prevent zillions of error messages being generated"
      )
      .def_property("exists", &TaoDataStruct::exists, &TaoDataStruct::set_exists, "See above")
      .def_property(
          "good_model",
          &TaoDataStruct::good_model,
          &TaoDataStruct::set_good_model,
          "See above"
      )
      .def_property(
          "good_base",
          &TaoDataStruct::good_base,
          &TaoDataStruct::set_good_base,
          "See above"
      )
      .def_property(
          "good_design",
          &TaoDataStruct::good_design,
          &TaoDataStruct::set_good_design,
          "See above"
      )
      .def_property(
          "good_meas",
          &TaoDataStruct::good_meas,
          &TaoDataStruct::set_good_meas,
          "See above"
      )
      .def_property("good_ref", &TaoDataStruct::good_ref, &TaoDataStruct::set_good_ref, "See above")
      .def_property(
          "good_user",
          &TaoDataStruct::good_user,
          &TaoDataStruct::set_good_user,
          "See above"
      )
      .def_property("good_opt", &TaoDataStruct::good_opt, &TaoDataStruct::set_good_opt, "See above")
      .def_property(
          "good_plot",
          &TaoDataStruct::good_plot,
          &TaoDataStruct::set_good_plot,
          "See above"
      )
      .def_property(
          "useit_plot",
          &TaoDataStruct::useit_plot,
          &TaoDataStruct::set_useit_plot,
          "See above"
      )
      .def_property(
          "useit_opt",
          &TaoDataStruct::useit_opt,
          &TaoDataStruct::set_useit_opt,
          "See above"
      )
      .def_property(
          "spin_map",
          py::cpp_function(&TaoDataStruct::spin_map, py::keep_alive<0, 1>()),
          &TaoDataStruct::set_spin_map
      )
      .def_property(
          "d1",
          py::cpp_function(&TaoDataStruct::d1, py::keep_alive<0, 1>()),
          &TaoDataStruct::set_d1,
          "Pointer to the parent d1_data_struct"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoDataStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoDataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoDataStruct &self, const TaoDataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoDataStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoDataStructArray1D, TaoDataStructAlloc1D>(
      m,
      "TaoDataStructArray1D",
      "TaoDataStructAlloc1D"
  );
  // 2D TaoDataStruct arrays are not used in structs/routines
  // 3D TaoDataStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_data_var_component_struct
void init_tao_data_var_component_struct(py::module &m, py::class_<TaoDataVarComponentStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("sign") = py::none()
  )
      .def_property(
          "name",
          &TaoDataVarComponentStruct::name,
          &TaoDataVarComponentStruct::set_name,
          "Eg: 'meas', 'ref', 'model', etc."
      )
      .def_property(
          "sign",
          &TaoDataVarComponentStruct::sign,
          &TaoDataVarComponentStruct::set_sign,
          "+1 or -1"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoDataVarComponentStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoDataVarComponentStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoDataVarComponentStruct &self, const TaoDataVarComponentStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoDataVarComponentStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoDataVarComponentStructArray1D, TaoDataVarComponentStructAlloc1D>(
      m,
      "TaoDataVarComponentStructArray1D",
      "TaoDataVarComponentStructAlloc1D"
  );
  // 2D TaoDataVarComponentStruct arrays are not used in structs/routines
  // 3D TaoDataVarComponentStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_drawing_struct
void init_tao_drawing_struct(py::module &m, py::class_<TaoDrawingStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly(
          "ele_shape",
          py::cpp_function(&TaoDrawingStruct::ele_shape, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const TaoDrawingStruct &self, const TaoDrawingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoDrawingStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoDrawingStruct arrays are not used in structs/routines
  // 2D TaoDrawingStruct arrays are not used in structs/routines
  // 3D TaoDrawingStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_dynamic_aperture_struct
void init_tao_dynamic_aperture_struct(py::module &m, py::class_<TaoDynamicApertureStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const ApertureParamStruct>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("param") = py::none(),
         py::arg("pz") = py::none(),
         py::arg("ellipse_scale") = py::none(),
         py::arg("a_emit") = py::none(),
         py::arg("b_emit") = py::none()
  )
      .def_property(
          "param",
          py::cpp_function(&TaoDynamicApertureStruct::param, py::keep_alive<0, 1>()),
          &TaoDynamicApertureStruct::set_param
      )
      .def_property_readonly(
          "scan",
          py::cpp_function(&TaoDynamicApertureStruct::scan, py::keep_alive<0, 1>()),
          "One scan for each pz."
      )
      .def_property(
          "pz",
          py::cpp_function(&TaoDynamicApertureStruct::pz, py::keep_alive<0, 1>()),
          &TaoDynamicApertureStruct::set_pz
      )
      .def_property(
          "ellipse_scale",
          &TaoDynamicApertureStruct::ellipse_scale,
          &TaoDynamicApertureStruct::set_ellipse_scale
      )
      .def_property(
          "a_emit",
          &TaoDynamicApertureStruct::a_emit,
          &TaoDynamicApertureStruct::set_a_emit
      )
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
      .def(
          "__eq__",
          [](const TaoDynamicApertureStruct &self, const TaoDynamicApertureStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoDynamicApertureStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(py::init<std::optional<int>>(), py::arg("n_loc") = py::none())
      .def_property_readonly(
          "eles",
          py::cpp_function(&TaoElePointerStruct::eles, py::keep_alive<0, 1>())
      )
      .def_property("n_loc", &TaoElePointerStruct::n_loc, &TaoElePointerStruct::set_n_loc)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoElePointerStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoElePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoElePointerStruct &self, const TaoElePointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoElePointerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoElePointerStructArray1D, TaoElePointerStructAlloc1D>(
      m,
      "TaoElePointerStructArray1D",
      "TaoElePointerStructAlloc1D"
  );
  // 2D TaoElePointerStruct arrays are not used in structs/routines
  // 3D TaoElePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_shape_struct
void init_tao_ele_shape_struct(py::module &m, py::class_<TaoEleShapeStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<double>,
             std::optional<int>,
             std::optional<std::string>>(),
         py::arg("ele_id") = py::none(),
         py::arg("shape") = py::none(),
         py::arg("color") = py::none(),
         py::arg("size") = py::none(),
         py::arg("label") = py::none(),
         py::arg("draw") = py::none(),
         py::arg("multi") = py::none(),
         py::arg("line_width") = py::none(),
         py::arg("offset") = py::none(),
         py::arg("ix_key") = py::none(),
         py::arg("name_ele") = py::none()
  )
      .def_property(
          "ele_id",
          &TaoEleShapeStruct::ele_id,
          &TaoEleShapeStruct::set_ele_id,
          "element 'key::name' to match to."
      )
      .def_property(
          "shape",
          &TaoEleShapeStruct::shape,
          &TaoEleShapeStruct::set_shape,
          "Shape to draw"
      )
      .def_property(
          "color",
          &TaoEleShapeStruct::color,
          &TaoEleShapeStruct::set_color,
          "Color of shape"
      )
      .def_property(
          "size",
          &TaoEleShapeStruct::size,
          &TaoEleShapeStruct::set_size,
          "plot vertical height"
      )
      .def_property(
          "label",
          &TaoEleShapeStruct::label,
          &TaoEleShapeStruct::set_label,
          "Can be: 'name', 's', 'none'"
      )
      .def_property(
          "draw",
          &TaoEleShapeStruct::draw,
          &TaoEleShapeStruct::set_draw,
          "Draw the shape?"
      )
      .def_property(
          "multi",
          &TaoEleShapeStruct::multi,
          &TaoEleShapeStruct::set_multi,
          "Can be part of a multi-shape."
      )
      .def_property(
          "line_width",
          &TaoEleShapeStruct::line_width,
          &TaoEleShapeStruct::set_line_width,
          "Width of lines used to draw the shape."
      )
      .def_property(
          "offset",
          &TaoEleShapeStruct::offset,
          &TaoEleShapeStruct::set_offset,
          "Vertical offset."
      )
      .def_property(
          "ix_key",
          &TaoEleShapeStruct::ix_key,
          &TaoEleShapeStruct::set_ix_key,
          "Extracted from ele_id. 0 => all classes (quadrupole, etc.)"
      )
      .def_property(
          "name_ele",
          &TaoEleShapeStruct::name_ele,
          &TaoEleShapeStruct::set_name_ele,
          "Name of element."
      )
      .def_property_readonly(
          "uni",
          py::cpp_function(&TaoEleShapeStruct::uni, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoEleShapeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoEleShapeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoEleShapeStruct &self, const TaoEleShapeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoEleShapeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoEleShapeStructArray1D, TaoEleShapeStructAlloc1D>(
      m,
      "TaoEleShapeStructArray1D",
      "TaoEleShapeStructAlloc1D"
  );
  // 2D TaoEleShapeStruct arrays are not used in structs/routines
  // 3D TaoEleShapeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_eval_node_struct
void init_tao_eval_node_struct(py::module &m, py::class_<TaoEvalNodeStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         py::arg("type") = py::none(),
         py::arg("name") = py::none(),
         py::arg("scale") = py::none(),
         py::arg("value") = py::none()
  )
      .def_property("type", &TaoEvalNodeStruct::type, &TaoEvalNodeStruct::set_type)
      .def_property("name", &TaoEvalNodeStruct::name, &TaoEvalNodeStruct::set_name)
      .def_property(
          "scale",
          &TaoEvalNodeStruct::scale,
          &TaoEvalNodeStruct::set_scale,
          "Scale factor for ping data"
      )
      .def_property(
          "value",
          py::cpp_function(&TaoEvalNodeStruct::value, py::keep_alive<0, 1>()),
          &TaoEvalNodeStruct::set_value
      )
      .def_property_readonly(
          "info",
          py::cpp_function(&TaoEvalNodeStruct::info, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "node",
          py::cpp_function(&TaoEvalNodeStruct::node, py::keep_alive<0, 1>()),
          "Child nodes for tree construction."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoEvalNodeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoEvalNodeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoEvalNodeStruct &self, const TaoEvalNodeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoEvalNodeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoEvalNodeStructArray1D, TaoEvalNodeStructAlloc1D>(
      m,
      "TaoEvalNodeStructArray1D",
      "TaoEvalNodeStructAlloc1D"
  );
  // 2D TaoEvalNodeStruct arrays are not used in structs/routines
  // 3D TaoEvalNodeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_expression_info_struct
void init_tao_expression_info_struct(py::module &m, py::class_<TaoExpressionInfoStruct> &cls) {
  cls.def(
         py::init<std::optional<bool>, optional_ref<const EleStruct>, std::optional<double>>(),
         py::arg("good") = py::none(),
         py::arg("ele") = py::none(),
         py::arg("s") = py::none()
  )
      .def_property(
          "good",
          &TaoExpressionInfoStruct::good,
          &TaoExpressionInfoStruct::set_good,
          "Expression is valid."
      )
      .def_property(
          "ele",
          py::cpp_function(&TaoExpressionInfoStruct::ele, py::keep_alive<0, 1>()),
          &TaoExpressionInfoStruct::set_ele,
          "Associated ele if it exists"
      )
      .def_property(
          "s",
          &TaoExpressionInfoStruct::s,
          &TaoExpressionInfoStruct::set_s,
          "Longitudinal position of expression."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoExpressionInfoStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoExpressionInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoExpressionInfoStruct &self, const TaoExpressionInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoExpressionInfoStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoExpressionInfoStructArray1D, TaoExpressionInfoStructAlloc1D>(
      m,
      "TaoExpressionInfoStructArray1D",
      "TaoExpressionInfoStructAlloc1D"
  );
  // 2D TaoExpressionInfoStruct arrays are not used in structs/routines
  // 3D TaoExpressionInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_floor_plan_struct
void init_tao_floor_plan_struct(py::module &m, py::class_<TaoFloorPlanStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>>(),
         py::arg("view") = py::none(),
         py::arg("rotation") = py::none(),
         py::arg("correct_distortion") = py::none(),
         py::arg("flip_label_side") = py::none(),
         py::arg("size_is_absolute") = py::none(),
         py::arg("draw_only_first_pass") = py::none(),
         py::arg("draw_building_wall") = py::none(),
         py::arg("orbit_scale") = py::none(),
         py::arg("orbit_color") = py::none(),
         py::arg("orbit_pattern") = py::none(),
         py::arg("orbit_lattice") = py::none(),
         py::arg("orbit_width") = py::none()
  )
      .def_property("view", &TaoFloorPlanStruct::view, &TaoFloorPlanStruct::set_view, "or 'xz'.")
      .def_property(
          "rotation",
          &TaoFloorPlanStruct::rotation,
          &TaoFloorPlanStruct::set_rotation,
          "Rotation of floor plan plot: 1.0 -> 360^deg"
      )
      .def_property(
          "correct_distortion",
          &TaoFloorPlanStruct::correct_distortion,
          &TaoFloorPlanStruct::set_correct_distortion,
          "T -> Shrink one axis so x-scale = y-scale."
      )
      .def_property(
          "flip_label_side",
          &TaoFloorPlanStruct::flip_label_side,
          &TaoFloorPlanStruct::set_flip_label_side,
          "Draw element label on other side of element?"
      )
      .def_property(
          "size_is_absolute",
          &TaoFloorPlanStruct::size_is_absolute,
          &TaoFloorPlanStruct::set_size_is_absolute,
          "Are shape sizes in meters or window pixels?"
      )
      .def_property(
          "draw_only_first_pass",
          &TaoFloorPlanStruct::draw_only_first_pass,
          &TaoFloorPlanStruct::set_draw_only_first_pass,
          "Draw only first pass with multipass elements?"
      )
      .def_property(
          "draw_building_wall",
          &TaoFloorPlanStruct::draw_building_wall,
          &TaoFloorPlanStruct::set_draw_building_wall,
          "Draw the building wall?"
      )
      .def_property(
          "orbit_scale",
          &TaoFloorPlanStruct::orbit_scale,
          &TaoFloorPlanStruct::set_orbit_scale,
          "Scale factor for drawing orbits. 0 -> Do not draw."
      )
      .def_property(
          "orbit_color",
          &TaoFloorPlanStruct::orbit_color,
          &TaoFloorPlanStruct::set_orbit_color
      )
      .def_property(
          "orbit_pattern",
          &TaoFloorPlanStruct::orbit_pattern,
          &TaoFloorPlanStruct::set_orbit_pattern
      )
      .def_property(
          "orbit_lattice",
          &TaoFloorPlanStruct::orbit_lattice,
          &TaoFloorPlanStruct::set_orbit_lattice,
          "Or 'design' or 'base'"
      )
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
      .def(
          "__eq__",
          [](const TaoFloorPlanStruct &self, const TaoFloorPlanStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoFloorPlanStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoFloorPlanStruct arrays are not used in structs/routines
  // 2D TaoFloorPlanStruct arrays are not used in structs/routines
  // 3D TaoFloorPlanStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_global_struct
void init_tao_global_struct(py::module &m, py::class_<TaoGlobalStruct> &cls) {
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
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
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
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("beam_dead_cutoff") = py::none(),
         py::arg("lm_opt_deriv_reinit") = py::none(),
         py::arg("de_lm_step_ratio") = py::none(),
         py::arg("de_var_to_population_factor") = py::none(),
         py::arg("lmdif_eps") = py::none(),
         py::arg("lmdif_negligible_merit") = py::none(),
         py::arg("svd_cutoff") = py::none(),
         py::arg("unstable_penalty") = py::none(),
         py::arg("merit_stop_value") = py::none(),
         py::arg("dmerit_stop_value") = py::none(),
         py::arg("random_sigma_cutoff") = py::none(),
         py::arg("delta_e_chrom") = py::none(),
         py::arg("max_plot_time") = py::none(),
         py::arg("default_universe") = py::none(),
         py::arg("default_branch") = py::none(),
         py::arg("n_opti_cycles") = py::none(),
         py::arg("n_opti_loops") = py::none(),
         py::arg("n_threads") = py::none(),
         py::arg("phase_units") = py::none(),
         py::arg("bunch_to_plot") = py::none(),
         py::arg("random_seed") = py::none(),
         py::arg("n_top10_merit") = py::none(),
         py::arg("srdt_gen_n_slices") = py::none(),
         py::arg("datum_err_messages_max") = py::none(),
         py::arg("srdt_sxt_n_slices") = py::none(),
         py::arg("srdt_use_cache") = py::none(),
         py::arg("quiet") = py::none(),
         py::arg("random_engine") = py::none(),
         py::arg("random_gauss_converter") = py::none(),
         py::arg("track_type") = py::none(),
         py::arg("lat_sigma_calc_uses_emit_from") = py::none(),
         py::arg("prompt_string") = py::none(),
         py::arg("prompt_color") = py::none(),
         py::arg("optimizer") = py::none(),
         py::arg("print_command") = py::none(),
         py::arg("var_out_file") = py::none(),
         py::arg("history_file") = py::none(),
         py::arg("beam_timer_on") = py::none(),
         py::arg("box_plots") = py::none(),
         py::arg("blank_line_between_commands") = py::none(),
         py::arg("cmd_file_abort_on_error") = py::none(),
         py::arg("concatenate_maps") = py::none(),
         py::arg("derivative_recalc") = py::none(),
         py::arg("derivative_uses_design") = py::none(),
         py::arg("disable_smooth_line_calc") = py::none(),
         py::arg("draw_curve_off_scale_warn") = py::none(),
         py::arg("external_plotting") = py::none(),
         py::arg("label_lattice_elements") = py::none(),
         py::arg("label_keys") = py::none(),
         py::arg("lattice_calc_on") = py::none(),
         py::arg("only_limit_opt_vars") = py::none(),
         py::arg("opt_with_ref") = py::none(),
         py::arg("opt_with_base") = py::none(),
         py::arg("opt_match_auto_recalc") = py::none(),
         py::arg("opti_write_var_file") = py::none(),
         py::arg("optimizer_allow_user_abort") = py::none(),
         py::arg("optimizer_var_limit_warn") = py::none(),
         py::arg("plot_on") = py::none(),
         py::arg("rad_int_user_calc_on") = py::none(),
         py::arg("rf_on") = py::none(),
         py::arg("single_step") = py::none(),
         py::arg("stop_on_error") = py::none(),
         py::arg("svd_retreat_on_merit_increase") = py::none(),
         py::arg("var_limits_on") = py::none(),
         py::arg("wait_for_CR_in_single_mode") = py::none(),
         py::arg("symbol_import") = py::none(),
         py::arg("debug_on") = py::none(),
         py::arg("expression_tree_on") = py::none(),
         py::arg("verbose_on") = py::none()
  )
      .def_property(
          "beam_dead_cutoff",
          &TaoGlobalStruct::beam_dead_cutoff,
          &TaoGlobalStruct::set_beam_dead_cutoff,
          "Percentage of dead particles at which beam tracking is stopped."
      )
      .def_property(
          "lm_opt_deriv_reinit",
          &TaoGlobalStruct::lm_opt_deriv_reinit,
          &TaoGlobalStruct::set_lm_opt_deriv_reinit,
          "Reinit derivative matrix cutoff"
      )
      .def_property(
          "de_lm_step_ratio",
          &TaoGlobalStruct::de_lm_step_ratio,
          &TaoGlobalStruct::set_de_lm_step_ratio,
          "Scaling for step sizes between DE and LM optimizers."
      )
      .def_property(
          "de_var_to_population_factor",
          &TaoGlobalStruct::de_var_to_population_factor,
          &TaoGlobalStruct::set_de_var_to_population_factor,
          "DE population = max(n_var*factor, 20)"
      )
      .def_property(
          "lmdif_eps",
          &TaoGlobalStruct::lmdif_eps,
          &TaoGlobalStruct::set_lmdif_eps,
          "Tollerance for lmdif optimizer."
      )
      .def_property(
          "lmdif_negligible_merit",
          &TaoGlobalStruct::lmdif_negligible_merit,
          &TaoGlobalStruct::set_lmdif_negligible_merit
      )
      .def_property(
          "svd_cutoff",
          &TaoGlobalStruct::svd_cutoff,
          &TaoGlobalStruct::set_svd_cutoff,
          "SVD singular value cutoff."
      )
      .def_property(
          "unstable_penalty",
          &TaoGlobalStruct::unstable_penalty,
          &TaoGlobalStruct::set_unstable_penalty,
          "Used in unstable_ring datum merit calculation."
      )
      .def_property(
          "merit_stop_value",
          &TaoGlobalStruct::merit_stop_value,
          &TaoGlobalStruct::set_merit_stop_value,
          "Merit value below which an optimizer will stop."
      )
      .def_property(
          "dmerit_stop_value",
          &TaoGlobalStruct::dmerit_stop_value,
          &TaoGlobalStruct::set_dmerit_stop_value,
          "Fractional Merit change below which an optimizer will stop."
      )
      .def_property(
          "random_sigma_cutoff",
          &TaoGlobalStruct::random_sigma_cutoff,
          &TaoGlobalStruct::set_random_sigma_cutoff,
          "Cut-off in sigmas."
      )
      .def_property(
          "delta_e_chrom",
          &TaoGlobalStruct::delta_e_chrom,
          &TaoGlobalStruct::set_delta_e_chrom,
          "Delta E used from chrom calc."
      )
      .def_property(
          "max_plot_time",
          &TaoGlobalStruct::max_plot_time,
          &TaoGlobalStruct::set_max_plot_time,
          "If plotting time (seconds) exceeds this than a message is generated."
      )
      .def_property(
          "default_universe",
          &TaoGlobalStruct::default_universe,
          &TaoGlobalStruct::set_default_universe,
          "Default universe to work with."
      )
      .def_property(
          "default_branch",
          &TaoGlobalStruct::default_branch,
          &TaoGlobalStruct::set_default_branch,
          "Default lattice branch to work with."
      )
      .def_property(
          "n_opti_cycles",
          &TaoGlobalStruct::n_opti_cycles,
          &TaoGlobalStruct::set_n_opti_cycles,
          "Number of optimization cycles"
      )
      .def_property(
          "n_opti_loops",
          &TaoGlobalStruct::n_opti_loops,
          &TaoGlobalStruct::set_n_opti_loops,
          "Number of optimization loops"
      )
      .def_property(
          "n_threads",
          &TaoGlobalStruct::n_threads,
          &TaoGlobalStruct::set_n_threads,
          "Number of OpenMP threads for parallel calculations."
      )
      .def_property(
          "phase_units",
          &TaoGlobalStruct::phase_units,
          &TaoGlobalStruct::set_phase_units,
          "Phase units on output."
      )
      .def_property(
          "bunch_to_plot",
          &TaoGlobalStruct::bunch_to_plot,
          &TaoGlobalStruct::set_bunch_to_plot,
          "Which bunch to plot"
      )
      .def_property(
          "random_seed",
          &TaoGlobalStruct::random_seed,
          &TaoGlobalStruct::set_random_seed,
          "Use system clock by default"
      )
      .def_property(
          "n_top10_merit",
          &TaoGlobalStruct::n_top10_merit,
          &TaoGlobalStruct::set_n_top10_merit,
          "Number of top merit constraints to print."
      )
      .def_property(
          "srdt_gen_n_slices",
          &TaoGlobalStruct::srdt_gen_n_slices,
          &TaoGlobalStruct::set_srdt_gen_n_slices,
          "Number times to slice elements for summation RDT calculation"
      )
      .def_property(
          "datum_err_messages_max",
          &TaoGlobalStruct::datum_err_messages_max,
          &TaoGlobalStruct::set_datum_err_messages_max,
          "Maximum number of error messages per call to lattice_calc."
      )
      .def_property(
          "srdt_sxt_n_slices",
          &TaoGlobalStruct::srdt_sxt_n_slices,
          &TaoGlobalStruct::set_srdt_sxt_n_slices,
          "Number times to slice sextupoles for summation RDT calculation"
      )
      .def_property(
          "srdt_use_cache",
          &TaoGlobalStruct::srdt_use_cache,
          &TaoGlobalStruct::set_srdt_use_cache,
          "Create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large."
      )
      .def_property(
          "quiet",
          &TaoGlobalStruct::quiet,
          &TaoGlobalStruct::set_quiet,
          "Print I/O when running a command file?"
      )
      .def_property(
          "random_engine",
          &TaoGlobalStruct::random_engine,
          &TaoGlobalStruct::set_random_engine,
          "Non-beam random number engine"
      )
      .def_property(
          "random_gauss_converter",
          &TaoGlobalStruct::random_gauss_converter,
          &TaoGlobalStruct::set_random_gauss_converter,
          "Non-beam"
      )
      .def_property(
          "track_type",
          &TaoGlobalStruct::track_type,
          &TaoGlobalStruct::set_track_type,
          "or 'beam'"
      )
      .def_property(
          "lat_sigma_calc_uses_emit_from",
          &TaoGlobalStruct::lat_sigma_calc_uses_emit_from,
          &TaoGlobalStruct::set_lat_sigma_calc_uses_emit_from,
          "Lattice derived sigma matrix uses emit values from where? Other possibilities: 'beam', "
          "'beam_init'."
      )
      .def_property(
          "prompt_string",
          &TaoGlobalStruct::prompt_string,
          &TaoGlobalStruct::set_prompt_string
      )
      .def_property(
          "prompt_color",
          &TaoGlobalStruct::prompt_color,
          &TaoGlobalStruct::set_prompt_color,
          "See read_a_line routine for possible settings."
      )
      .def_property(
          "optimizer",
          &TaoGlobalStruct::optimizer,
          &TaoGlobalStruct::set_optimizer,
          "optimizer to use."
      )
      .def_property(
          "print_command",
          &TaoGlobalStruct::print_command,
          &TaoGlobalStruct::set_print_command
      )
      .def_property(
          "var_out_file",
          &TaoGlobalStruct::var_out_file,
          &TaoGlobalStruct::set_var_out_file
      )
      .def_property(
          "history_file",
          &TaoGlobalStruct::history_file,
          &TaoGlobalStruct::set_history_file
      )
      .def_property(
          "beam_timer_on",
          &TaoGlobalStruct::beam_timer_on,
          &TaoGlobalStruct::set_beam_timer_on,
          "For timing the beam tracking calculation."
      )
      .def_property(
          "box_plots",
          &TaoGlobalStruct::box_plots,
          &TaoGlobalStruct::set_box_plots,
          "For debugging plot layout issues."
      )
      .def_property(
          "blank_line_between_commands",
          &TaoGlobalStruct::blank_line_between_commands,
          &TaoGlobalStruct::set_blank_line_between_commands,
          "Add a blank line between command output?"
      )
      .def_property(
          "cmd_file_abort_on_error",
          &TaoGlobalStruct::cmd_file_abort_on_error,
          &TaoGlobalStruct::set_cmd_file_abort_on_error,
          "Abort open command files if there is an error?"
      )
      .def_property(
          "concatenate_maps",
          &TaoGlobalStruct::concatenate_maps,
          &TaoGlobalStruct::set_concatenate_maps,
          "False => tracking using DA."
      )
      .def_property(
          "derivative_recalc",
          &TaoGlobalStruct::derivative_recalc,
          &TaoGlobalStruct::set_derivative_recalc,
          "Recalc before each optimizer run?"
      )
      .def_property(
          "derivative_uses_design",
          &TaoGlobalStruct::derivative_uses_design,
          &TaoGlobalStruct::set_derivative_uses_design,
          "Derivative calc uses design lattice instead of model?"
      )
      .def_property(
          "disable_smooth_line_calc",
          &TaoGlobalStruct::disable_smooth_line_calc,
          &TaoGlobalStruct::set_disable_smooth_line_calc,
          "Global disable of the smooth line calculation."
      )
      .def_property(
          "draw_curve_off_scale_warn",
          &TaoGlobalStruct::draw_curve_off_scale_warn,
          &TaoGlobalStruct::set_draw_curve_off_scale_warn,
          "Display warning on graphs?"
      )
      .def_property(
          "external_plotting",
          &TaoGlobalStruct::external_plotting,
          &TaoGlobalStruct::set_external_plotting,
          "Used with matplotlib and gui."
      )
      .def_property(
          "label_lattice_elements",
          &TaoGlobalStruct::label_lattice_elements,
          &TaoGlobalStruct::set_label_lattice_elements,
          "For lat_layout plots"
      )
      .def_property(
          "label_keys",
          &TaoGlobalStruct::label_keys,
          &TaoGlobalStruct::set_label_keys,
          "For lat_layout plots"
      )
      .def_property(
          "lattice_calc_on",
          &TaoGlobalStruct::lattice_calc_on,
          &TaoGlobalStruct::set_lattice_calc_on,
          "Turn on/off beam and single particle calculations."
      )
      .def_property(
          "only_limit_opt_vars",
          &TaoGlobalStruct::only_limit_opt_vars,
          &TaoGlobalStruct::set_only_limit_opt_vars,
          "Only apply limits to variables used in optimization."
      )
      .def_property(
          "opt_with_ref",
          &TaoGlobalStruct::opt_with_ref,
          &TaoGlobalStruct::set_opt_with_ref,
          "Use reference data in optimization?"
      )
      .def_property(
          "opt_with_base",
          &TaoGlobalStruct::opt_with_base,
          &TaoGlobalStruct::set_opt_with_base,
          "Use base data in optimization?"
      )
      .def_property(
          "opt_match_auto_recalc",
          &TaoGlobalStruct::opt_match_auto_recalc,
          &TaoGlobalStruct::set_opt_match_auto_recalc,
          "Set recalc = True for match elements before each cycle?"
      )
      .def_property(
          "opti_write_var_file",
          &TaoGlobalStruct::opti_write_var_file,
          &TaoGlobalStruct::set_opti_write_var_file,
          "'run' command writes var_out_file"
      )
      .def_property(
          "optimizer_allow_user_abort",
          &TaoGlobalStruct::optimizer_allow_user_abort,
          &TaoGlobalStruct::set_optimizer_allow_user_abort,
          "See Tao manual for more details."
      )
      .def_property(
          "optimizer_var_limit_warn",
          &TaoGlobalStruct::optimizer_var_limit_warn,
          &TaoGlobalStruct::set_optimizer_var_limit_warn,
          "Warn when vars reach a limit with optimization."
      )
      .def_property(
          "plot_on",
          &TaoGlobalStruct::plot_on,
          &TaoGlobalStruct::set_plot_on,
          "Do plotting?"
      )
      .def_property(
          "rad_int_user_calc_on",
          &TaoGlobalStruct::rad_int_user_calc_on,
          &TaoGlobalStruct::set_rad_int_user_calc_on,
          "User set radiation integrals calculation on/off."
      )
      .def_property(
          "rf_on",
          &TaoGlobalStruct::rf_on,
          &TaoGlobalStruct::set_rf_on,
          "RFcavities on or off? Does not affect lcavities."
      )
      .def_property(
          "single_step",
          &TaoGlobalStruct::single_step,
          &TaoGlobalStruct::set_single_step,
          "For debugging and demonstrations: Single step through a command file?"
      )
      .def_property(
          "stop_on_error",
          &TaoGlobalStruct::stop_on_error,
          &TaoGlobalStruct::set_stop_on_error,
          "For debugging: False prevents tao from exiting on an error."
      )
      .def_property(
          "svd_retreat_on_merit_increase",
          &TaoGlobalStruct::svd_retreat_on_merit_increase,
          &TaoGlobalStruct::set_svd_retreat_on_merit_increase
      )
      .def_property(
          "var_limits_on",
          &TaoGlobalStruct::var_limits_on,
          &TaoGlobalStruct::set_var_limits_on,
          "Respect the variable limits?"
      )
      .def_property(
          "wait_for_CR_in_single_mode",
          &TaoGlobalStruct::wait_for_CR_in_single_mode,
          &TaoGlobalStruct::set_wait_for_CR_in_single_mode,
          "For use with a python GUI."
      )
      .def_property(
          "symbol_import",
          &TaoGlobalStruct::symbol_import,
          &TaoGlobalStruct::set_symbol_import,
          "Import symbols from lattice file(s)? Internal stuff"
      )
      .def_property(
          "debug_on",
          &TaoGlobalStruct::debug_on,
          &TaoGlobalStruct::set_debug_on,
          "For debugging."
      )
      .def_property(
          "expression_tree_on",
          &TaoGlobalStruct::expression_tree_on,
          &TaoGlobalStruct::set_expression_tree_on,
          "Use an expression tree instead of a stack?"
      )
      .def_property(
          "verbose_on",
          &TaoGlobalStruct::verbose_on,
          &TaoGlobalStruct::set_verbose_on,
          "For verbose output. Used with debugging."
      )

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
      .def(
          "__eq__",
          [](const TaoGlobalStruct &self, const TaoGlobalStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoGlobalStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoGlobalStruct arrays are not used in structs/routines
  // 2D TaoGlobalStruct arrays are not used in structs/routines
  // 3D TaoGlobalStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_graph_struct
void init_tao_graph_struct(py::module &m, py::class_<TaoGraphStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             optional_ref<const TaoPlotStruct>,
             optional_ref<const TaoFloorPlanStruct>,
             optional_ref<const QpPointStruct>,
             optional_ref<const QpPointStruct>,
             optional_ref<const QpLegendStruct>,
             optional_ref<const QpAxisStruct>,
             optional_ref<const QpAxisStruct>,
             optional_ref<const QpAxisStruct>,
             optional_ref<const QpAxisStruct>,
             optional_ref<const QpRectStruct>,
             optional_ref<const QpRectStruct>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<int>>,
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
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("title") = py::none(),
         py::arg("title_suffix") = py::none(),
         py::arg("why_invalid") = py::none(),
         py::arg("p") = py::none(),
         py::arg("floor_plan") = py::none(),
         py::arg("text_legend_origin") = py::none(),
         py::arg("curve_legend_origin") = py::none(),
         py::arg("curve_legend") = py::none(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("x2") = py::none(),
         py::arg("y2") = py::none(),
         py::arg("margin") = py::none(),
         py::arg("scale_margin") = py::none(),
         py::arg("x_axis_scale_factor") = py::none(),
         py::arg("symbol_size_scale") = py::none(),
         py::arg("box") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_universe") = py::none(),
         py::arg("clip") = py::none(),
         py::arg("y2_mirrors_y") = py::none(),
         py::arg("limited") = py::none(),
         py::arg("draw_axes") = py::none(),
         py::arg("draw_curve_legend") = py::none(),
         py::arg("draw_grid") = py::none(),
         py::arg("draw_title") = py::none(),
         py::arg("draw_only_good_user_data_or_vars") = py::none(),
         py::arg("allow_wrap_around") = py::none(),
         py::arg("is_valid") = py::none()
  )
      .def_property(
          "name",
          &TaoGraphStruct::name,
          &TaoGraphStruct::set_name,
          "Name identifying the graph"
      )
      .def_property(
          "type",
          &TaoGraphStruct::type,
          &TaoGraphStruct::set_type,
          "'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture'"
      )
      .def_property("title", &TaoGraphStruct::title, &TaoGraphStruct::set_title)
      .def_property(
          "title_suffix",
          &TaoGraphStruct::title_suffix,
          &TaoGraphStruct::set_title_suffix
      )
      .def_property_readonly(
          "text_legend",
          py::cpp_function(&TaoGraphStruct::text_legend, py::keep_alive<0, 1>()),
          "Array for holding descriptive info."
      )
      .def_property_readonly(
          "text_legend_out",
          py::cpp_function(&TaoGraphStruct::text_legend_out, py::keep_alive<0, 1>()),
          "Array for holding descriptive info."
      )
      .def_property(
          "why_invalid",
          &TaoGraphStruct::why_invalid,
          &TaoGraphStruct::set_why_invalid,
          "Informative string to print."
      )
      .def_property_readonly(
          "curve",
          py::cpp_function(&TaoGraphStruct::curve, py::keep_alive<0, 1>())
      )
      .def_property(
          "p",
          py::cpp_function(&TaoGraphStruct::p, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_p,
          "pointer to parent plot"
      )
      .def_property(
          "floor_plan",
          py::cpp_function(&TaoGraphStruct::floor_plan, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_floor_plan
      )
      .def_property(
          "text_legend_origin",
          py::cpp_function(&TaoGraphStruct::text_legend_origin, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_text_legend_origin
      )
      .def_property(
          "curve_legend_origin",
          py::cpp_function(&TaoGraphStruct::curve_legend_origin, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_curve_legend_origin
      )
      .def_property(
          "curve_legend",
          py::cpp_function(&TaoGraphStruct::curve_legend, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_curve_legend
      )
      .def_property(
          "x",
          py::cpp_function(&TaoGraphStruct::x, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_x,
          "X-axis parameters."
      )
      .def_property(
          "y",
          py::cpp_function(&TaoGraphStruct::y, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_y,
          "Y-axis attributes."
      )
      .def_property(
          "x2",
          py::cpp_function(&TaoGraphStruct::x2, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_x2,
          "X2-axis attributes (Not currently used)."
      )
      .def_property(
          "y2",
          py::cpp_function(&TaoGraphStruct::y2, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_y2,
          "Y2-axis attributes."
      )
      .def_property(
          "margin",
          py::cpp_function(&TaoGraphStruct::margin, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_margin,
          "Margin around the graph."
      )
      .def_property(
          "scale_margin",
          py::cpp_function(&TaoGraphStruct::scale_margin, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_scale_margin,
          "Margin for scaling"
      )
      .def_property(
          "x_axis_scale_factor",
          &TaoGraphStruct::x_axis_scale_factor,
          &TaoGraphStruct::set_x_axis_scale_factor,
          "x-axis conversion from internal to plotting units."
      )
      .def_property(
          "symbol_size_scale",
          &TaoGraphStruct::symbol_size_scale,
          &TaoGraphStruct::set_symbol_size_scale,
          "Symbol size scale factor for phase_space plots."
      )
      .def_property(
          "box",
          py::cpp_function(&TaoGraphStruct::box, py::keep_alive<0, 1>()),
          &TaoGraphStruct::set_box,
          "Defines which box the plot is put in."
      )
      .def_property(
          "ix_branch",
          &TaoGraphStruct::ix_branch,
          &TaoGraphStruct::set_ix_branch,
          "Branch in lattice. Used when there are no associated curves."
      )
      .def_property(
          "ix_universe",
          &TaoGraphStruct::ix_universe,
          &TaoGraphStruct::set_ix_universe,
          "Used for lat_layout plots."
      )
      .def_property(
          "clip",
          &TaoGraphStruct::clip,
          &TaoGraphStruct::set_clip,
          "Clip plot at graph boundary."
      )
      .def_property(
          "y2_mirrors_y",
          &TaoGraphStruct::y2_mirrors_y,
          &TaoGraphStruct::set_y2_mirrors_y,
          "Y2-axis same as Y-axis?"
      )
      .def_property(
          "limited",
          &TaoGraphStruct::limited,
          &TaoGraphStruct::set_limited,
          "True if at least one data point past graph bounds."
      )
      .def_property(
          "draw_axes",
          &TaoGraphStruct::draw_axes,
          &TaoGraphStruct::set_draw_axes,
          "Draw axes, labels, etc?"
      )
      .def_property(
          "draw_curve_legend",
          &TaoGraphStruct::draw_curve_legend,
          &TaoGraphStruct::set_draw_curve_legend,
          "Legend for displaying curve info."
      )
      .def_property(
          "draw_grid",
          &TaoGraphStruct::draw_grid,
          &TaoGraphStruct::set_draw_grid,
          "Draw a grid?"
      )
      .def_property("draw_title", &TaoGraphStruct::draw_title, &TaoGraphStruct::set_draw_title)
      .def_property(
          "draw_only_good_user_data_or_vars",
          &TaoGraphStruct::draw_only_good_user_data_or_vars,
          &TaoGraphStruct::set_draw_only_good_user_data_or_vars
      )
      .def_property(
          "allow_wrap_around",
          &TaoGraphStruct::allow_wrap_around,
          &TaoGraphStruct::set_allow_wrap_around,
          "'Wrap' curves to extend past lattice boundaries?"
      )
      .def_property(
          "is_valid",
          &TaoGraphStruct::is_valid,
          &TaoGraphStruct::set_is_valid,
          "EG: Bad x_axis_type."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoGraphStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoGraphStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoGraphStruct &self, const TaoGraphStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoGraphStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoGraphStructArray1D, TaoGraphStructAlloc1D>(
      m,
      "TaoGraphStructArray1D",
      "TaoGraphStructAlloc1D"
  );
  // 2D TaoGraphStruct arrays are not used in structs/routines
  // 3D TaoGraphStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_histogram_struct
void init_tao_histogram_struct(py::module &m, py::class_<TaoHistogramStruct> &cls) {
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         py::arg("density_normalized") = py::none(),
         py::arg("weight_by_charge") = py::none(),
         py::arg("minimum") = py::none(),
         py::arg("maximum") = py::none(),
         py::arg("width") = py::none(),
         py::arg("center") = py::none(),
         py::arg("number") = py::none()
  )
      .def_property(
          "density_normalized",
          &TaoHistogramStruct::density_normalized,
          &TaoHistogramStruct::set_density_normalized
      )
      .def_property(
          "weight_by_charge",
          &TaoHistogramStruct::weight_by_charge,
          &TaoHistogramStruct::set_weight_by_charge
      )
      .def_property(
          "minimum",
          &TaoHistogramStruct::minimum,
          &TaoHistogramStruct::set_minimum,
          "Computed by Tao. Not User settable."
      )
      .def_property(
          "maximum",
          &TaoHistogramStruct::maximum,
          &TaoHistogramStruct::set_maximum,
          "Computed by Tao. Not User settable."
      )
      .def_property("width", &TaoHistogramStruct::width, &TaoHistogramStruct::set_width)
      .def_property("center", &TaoHistogramStruct::center, &TaoHistogramStruct::set_center)
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
      .def(
          "__eq__",
          [](const TaoHistogramStruct &self, const TaoHistogramStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoHistogramStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoHistogramStruct arrays are not used in structs/routines
  // 2D TaoHistogramStruct arrays are not used in structs/routines
  // 3D TaoHistogramStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_init_struct
void init_tao_init_struct(py::module &m, py::class_<TaoInitStruct> &cls) {
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>>(),
         py::arg("parse_cmd_args") = py::none(),
         py::arg("debug_switch") = py::none(),
         py::arg("external_plotting_switch") = py::none(),
         py::arg("init_name") = py::none(),
         py::arg("hook_init_file") = py::none(),
         py::arg("hook_lat_file") = py::none(),
         py::arg("hook_beam_file") = py::none(),
         py::arg("hook_data_file") = py::none(),
         py::arg("hook_plot_file") = py::none(),
         py::arg("hook_startup_file") = py::none(),
         py::arg("hook_var_file") = py::none(),
         py::arg("hook_building_wall_file") = py::none(),
         py::arg("init_file_arg_path") = py::none(),
         py::arg("lattice_file_arg") = py::none(),
         py::arg("hook_init_file_arg") = py::none(),
         py::arg("init_file_arg") = py::none(),
         py::arg("beam_file_arg") = py::none(),
         py::arg("beam_init_position_file_arg") = py::none(),
         py::arg("command_arg") = py::none(),
         py::arg("data_file_arg") = py::none(),
         py::arg("plot_file_arg") = py::none(),
         py::arg("startup_file_arg") = py::none(),
         py::arg("var_file_arg") = py::none(),
         py::arg("building_wall_file_arg") = py::none(),
         py::arg("geometry_arg") = py::none(),
         py::arg("slice_lattice_arg") = py::none(),
         py::arg("start_branch_at_arg") = py::none(),
         py::arg("log_startup_arg") = py::none(),
         py::arg("no_stopping_arg") = py::none(),
         py::arg("noplot_arg") = py::none(),
         py::arg("no_rad_int_arg") = py::none(),
         py::arg("reverse_arg") = py::none(),
         py::arg("debug_arg") = py::none(),
         py::arg("disable_smooth_line_calc_arg") = py::none(),
         py::arg("rf_on_arg") = py::none(),
         py::arg("prompt_color_arg") = py::none(),
         py::arg("quiet_arg") = py::none(),
         py::arg("noinit_arg") = py::none(),
         py::arg("nostartup_arg") = py::none(),
         py::arg("symbol_import_arg") = py::none(),
         py::arg("unique_name_suffix") = py::none()
  )
      .def_property(
          "parse_cmd_args",
          &TaoInitStruct::parse_cmd_args,
          &TaoInitStruct::set_parse_cmd_args,
          "Used by custom programs to control Tao init"
      )
      .def_property(
          "debug_switch",
          &TaoInitStruct::debug_switch,
          &TaoInitStruct::set_debug_switch,
          "Is the '-debug' switch present?"
      )
      .def_property(
          "external_plotting_switch",
          &TaoInitStruct::external_plotting_switch,
          &TaoInitStruct::set_external_plotting_switch,
          "Is '-external_plotting' switch present?"
      )
      .def_property(
          "init_name",
          &TaoInitStruct::init_name,
          &TaoInitStruct::set_init_name,
          "label for initialization"
      )
      .def_property(
          "hook_init_file",
          &TaoInitStruct::hook_init_file,
          &TaoInitStruct::set_hook_init_file
      )
      .def_property(
          "hook_lat_file",
          &TaoInitStruct::hook_lat_file,
          &TaoInitStruct::set_hook_lat_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_beam_file",
          &TaoInitStruct::hook_beam_file,
          &TaoInitStruct::set_hook_beam_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_data_file",
          &TaoInitStruct::hook_data_file,
          &TaoInitStruct::set_hook_data_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_plot_file",
          &TaoInitStruct::hook_plot_file,
          &TaoInitStruct::set_hook_plot_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_startup_file",
          &TaoInitStruct::hook_startup_file,
          &TaoInitStruct::set_hook_startup_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_var_file",
          &TaoInitStruct::hook_var_file,
          &TaoInitStruct::set_hook_var_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "hook_building_wall_file",
          &TaoInitStruct::hook_building_wall_file,
          &TaoInitStruct::set_hook_building_wall_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_property(
          "init_file_arg_path",
          &TaoInitStruct::init_file_arg_path,
          &TaoInitStruct::set_init_file_arg_path,
          "Path part of init_tao_file"
      )
      .def_property(
          "lattice_file_arg",
          &TaoInitStruct::lattice_file_arg,
          &TaoInitStruct::set_lattice_file_arg,
          "-lattice_file        command line argument."
      )
      .def_property(
          "hook_init_file_arg",
          &TaoInitStruct::hook_init_file_arg,
          &TaoInitStruct::set_hook_init_file_arg,
          "-hook_init_file      command line argument"
      )
      .def_property(
          "init_file_arg",
          &TaoInitStruct::init_file_arg,
          &TaoInitStruct::set_init_file_arg,
          "-init_file           command line argument."
      )
      .def_property(
          "beam_file_arg",
          &TaoInitStruct::beam_file_arg,
          &TaoInitStruct::set_beam_file_arg,
          "-beam_file           command line argument."
      )
      .def_property(
          "beam_init_position_file_arg",
          &TaoInitStruct::beam_init_position_file_arg,
          &TaoInitStruct::set_beam_init_position_file_arg,
          "-beam_init_position_file command line argument."
      )
      .def_property(
          "command_arg",
          &TaoInitStruct::command_arg,
          &TaoInitStruct::set_command_arg,
          "-command             command line argument."
      )
      .def_property(
          "data_file_arg",
          &TaoInitStruct::data_file_arg,
          &TaoInitStruct::set_data_file_arg,
          "-data_file           command line argument."
      )
      .def_property(
          "plot_file_arg",
          &TaoInitStruct::plot_file_arg,
          &TaoInitStruct::set_plot_file_arg,
          "-plot_file           command line argument."
      )
      .def_property(
          "startup_file_arg",
          &TaoInitStruct::startup_file_arg,
          &TaoInitStruct::set_startup_file_arg,
          "-startup_file        command line argument."
      )
      .def_property(
          "var_file_arg",
          &TaoInitStruct::var_file_arg,
          &TaoInitStruct::set_var_file_arg,
          "-var_file            command line argument."
      )
      .def_property(
          "building_wall_file_arg",
          &TaoInitStruct::building_wall_file_arg,
          &TaoInitStruct::set_building_wall_file_arg,
          "-building_wall_file  command line argument."
      )
      .def_property(
          "geometry_arg",
          &TaoInitStruct::geometry_arg,
          &TaoInitStruct::set_geometry_arg,
          "-geometry            command line argument."
      )
      .def_property(
          "slice_lattice_arg",
          &TaoInitStruct::slice_lattice_arg,
          &TaoInitStruct::set_slice_lattice_arg,
          "-slice_lattice       command line argument."
      )
      .def_property(
          "start_branch_at_arg",
          &TaoInitStruct::start_branch_at_arg,
          &TaoInitStruct::set_start_branch_at_arg,
          "-start_branch_at     command line argument."
      )
      .def_property(
          "log_startup_arg",
          &TaoInitStruct::log_startup_arg,
          &TaoInitStruct::set_log_startup_arg,
          "-log_startup         command line argument"
      )
      .def_property(
          "no_stopping_arg",
          &TaoInitStruct::no_stopping_arg,
          &TaoInitStruct::set_no_stopping_arg,
          "-no_stopping         command line argument"
      )
      .def_property(
          "noplot_arg",
          &TaoInitStruct::noplot_arg,
          &TaoInitStruct::set_noplot_arg,
          "-noplot              command line argument"
      )
      .def_property(
          "no_rad_int_arg",
          &TaoInitStruct::no_rad_int_arg,
          &TaoInitStruct::set_no_rad_int_arg,
          "-no_rad_int          command line argument"
      )
      .def_property(
          "reverse_arg",
          &TaoInitStruct::reverse_arg,
          &TaoInitStruct::set_reverse_arg,
          "-reverse             command line argument"
      )
      .def_property(
          "debug_arg",
          &TaoInitStruct::debug_arg,
          &TaoInitStruct::set_debug_arg,
          "-debug               command line argument"
      )
      .def_property(
          "disable_smooth_line_calc_arg",
          &TaoInitStruct::disable_smooth_line_calc_arg,
          &TaoInitStruct::set_disable_smooth_line_calc_arg,
          "-disable_smooth_line_calc"
      )
      .def_property(
          "rf_on_arg",
          &TaoInitStruct::rf_on_arg,
          &TaoInitStruct::set_rf_on_arg,
          "-rf_on               command line argument"
      )
      .def_property(
          "prompt_color_arg",
          &TaoInitStruct::prompt_color_arg,
          &TaoInitStruct::set_prompt_color_arg,
          "-prompt_color        command line argument"
      )
      .def_property(
          "quiet_arg",
          &TaoInitStruct::quiet_arg,
          &TaoInitStruct::set_quiet_arg,
          "-quiet               command line argument"
      )
      .def_property(
          "noinit_arg",
          &TaoInitStruct::noinit_arg,
          &TaoInitStruct::set_noinit_arg,
          "-noinit              command line argument"
      )
      .def_property(
          "nostartup_arg",
          &TaoInitStruct::nostartup_arg,
          &TaoInitStruct::set_nostartup_arg,
          "-nostartup           command line argument"
      )
      .def_property(
          "symbol_import_arg",
          &TaoInitStruct::symbol_import_arg,
          &TaoInitStruct::set_symbol_import_arg,
          "-symbol_import       command line argument"
      )
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
      .def(
          "__eq__",
          [](const TaoInitStruct &self, const TaoInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoInitStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoInitStruct arrays are not used in structs/routines
  // 2D TaoInitStruct arrays are not used in structs/routines
  // 3D TaoInitStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lat_sigma_struct
void init_tao_lat_sigma_struct(py::module &m, py::class_<TaoLatSigmaStruct> &cls) {
  cls.def(py::init<std::optional<std::vector<std::vector<double>>>>(), py::arg("mat") = py::none())
      .def_property(
          "mat",
          py::cpp_function(&TaoLatSigmaStruct::mat, py::keep_alive<0, 1>()),
          &TaoLatSigmaStruct::set_mat
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoLatSigmaStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoLatSigmaStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoLatSigmaStruct &self, const TaoLatSigmaStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoLatSigmaStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoLatSigmaStructArray1D, TaoLatSigmaStructAlloc1D>(
      m,
      "TaoLatSigmaStructArray1D",
      "TaoLatSigmaStructAlloc1D"
  );
  // 2D TaoLatSigmaStruct arrays are not used in structs/routines
  // 3D TaoLatSigmaStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_branch_struct
void init_tao_lattice_branch_struct(py::module &m, py::class_<TaoLatticeBranchStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const TaoLatticeStruct>,
             optional_ref<const TaoSpinPolarizationStruct>,
             optional_ref<const SummationRdtStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const NormalModesStruct>,
             optional_ref<const NormalModesStruct>,
             optional_ref<const PtcNormalFormStruct>,
             optional_ref<const BmadNormalFormStruct>,
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
             std::optional<bool>>(),
         py::arg("tao_lat") = py::none(),
         py::arg("spin") = py::none(),
         py::arg("srdt") = py::none(),
         py::arg("orb0") = py::none(),
         py::arg("modes_ri") = py::none(),
         py::arg("modes_6d") = py::none(),
         py::arg("ptc_normal_form") = py::none(),
         py::arg("bmad_normal_form") = py::none(),
         py::arg("cache_x_min") = py::none(),
         py::arg("cache_x_max") = py::none(),
         py::arg("comb_ds_save") = py::none(),
         py::arg("ix_ref_taylor") = py::none(),
         py::arg("ix_ele_taylor") = py::none(),
         py::arg("track_state") = py::none(),
         py::arg("cache_n_pts") = py::none(),
         py::arg("ix_rad_int_cache") = py::none(),
         py::arg("has_open_match_element") = py::none(),
         py::arg("plot_cache_valid") = py::none(),
         py::arg("spin_map_valid") = py::none(),
         py::arg("twiss_valid") = py::none(),
         py::arg("mode_flip_here") = py::none(),
         py::arg("chrom_calc_ok") = py::none(),
         py::arg("rad_int_calc_ok") = py::none(),
         py::arg("emit_6d_calc_ok") = py::none(),
         py::arg("sigma_track_ok") = py::none()
  )
      .def_property(
          "tao_lat",
          py::cpp_function(&TaoLatticeBranchStruct::tao_lat, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_tao_lat,
          "Parent tao_lat"
      )
      .def_property_readonly(
          "lat_sigma",
          py::cpp_function(&TaoLatticeBranchStruct::lat_sigma, py::keep_alive<0, 1>()),
          "Sigma matrix derived from lattice (not beam)."
      )
      .def_property_readonly(
          "spin_ele",
          py::cpp_function(&TaoLatticeBranchStruct::spin_ele, py::keep_alive<0, 1>()),
          "Spin stuff"
      )
      .def_property_readonly(
          "bunch_params",
          py::cpp_function(&TaoLatticeBranchStruct::bunch_params, py::keep_alive<0, 1>()),
          "Per element"
      )
      .def_property_readonly(
          "bunch_params_comb",
          py::cpp_function(&TaoLatticeBranchStruct::bunch_params_comb, py::keep_alive<0, 1>()),
          "A comb for each bunch in beam."
      )
      .def_property_readonly(
          "orbit",
          py::cpp_function(&TaoLatticeBranchStruct::orbit, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "plot_cache",
          py::cpp_function(&TaoLatticeBranchStruct::plot_cache, py::keep_alive<0, 1>()),
          "Plotting data cache"
      )
      .def_property(
          "spin",
          py::cpp_function(&TaoLatticeBranchStruct::spin, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_spin
      )
      .def_property(
          "srdt",
          py::cpp_function(&TaoLatticeBranchStruct::srdt, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_srdt
      )
      .def_property(
          "orb0",
          py::cpp_function(&TaoLatticeBranchStruct::orb0, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_orb0,
          "For saving beginning orbit in closed geometry branches. orb0 can then be used as an "
          "initial guess when closed_orbit is called again."
      )
      .def_property(
          "modes_ri",
          py::cpp_function(&TaoLatticeBranchStruct::modes_ri, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_modes_ri,
          "Synchrotron integrals stuff"
      )
      .def_property(
          "modes_6d",
          py::cpp_function(&TaoLatticeBranchStruct::modes_6d, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_modes_6d,
          "6D radiation matrices."
      )
      .def_property(
          "ptc_normal_form",
          py::cpp_function(&TaoLatticeBranchStruct::ptc_normal_form, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_ptc_normal_form,
          "Collection of normal form structures defined in PTC"
      )
      .def_property(
          "bmad_normal_form",
          py::cpp_function(&TaoLatticeBranchStruct::bmad_normal_form, py::keep_alive<0, 1>()),
          &TaoLatticeBranchStruct::set_bmad_normal_form,
          "Collection of normal form structures defined in Bmad"
      )
      .def_property_readonly(
          "high_E_orb",
          py::cpp_function(&TaoLatticeBranchStruct::high_E_orb, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "low_E_orb",
          py::cpp_function(&TaoLatticeBranchStruct::low_E_orb, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "taylor_save",
          py::cpp_function(&TaoLatticeBranchStruct::taylor_save, py::keep_alive<0, 1>()),
          "Save to reduce computation time."
      )
      .def_property(
          "cache_x_min",
          &TaoLatticeBranchStruct::cache_x_min,
          &TaoLatticeBranchStruct::set_cache_x_min
      )
      .def_property(
          "cache_x_max",
          &TaoLatticeBranchStruct::cache_x_max,
          &TaoLatticeBranchStruct::set_cache_x_max
      )
      .def_property(
          "comb_ds_save",
          &TaoLatticeBranchStruct::comb_ds_save,
          &TaoLatticeBranchStruct::set_comb_ds_save,
          "Master parameter for %bunch_params_comb(:)%ds_save"
      )
      .def_property(
          "ix_ref_taylor",
          &TaoLatticeBranchStruct::ix_ref_taylor,
          &TaoLatticeBranchStruct::set_ix_ref_taylor
      )
      .def_property(
          "ix_ele_taylor",
          &TaoLatticeBranchStruct::ix_ele_taylor,
          &TaoLatticeBranchStruct::set_ix_ele_taylor
      )
      .def_property(
          "track_state",
          &TaoLatticeBranchStruct::track_state,
          &TaoLatticeBranchStruct::set_track_state
      )
      .def_property(
          "cache_n_pts",
          &TaoLatticeBranchStruct::cache_n_pts,
          &TaoLatticeBranchStruct::set_cache_n_pts
      )
      .def_property(
          "ix_rad_int_cache",
          &TaoLatticeBranchStruct::ix_rad_int_cache,
          &TaoLatticeBranchStruct::set_ix_rad_int_cache,
          "Radiation integrals cache index."
      )
      .def_property(
          "has_open_match_element",
          &TaoLatticeBranchStruct::has_open_match_element,
          &TaoLatticeBranchStruct::set_has_open_match_element
      )
      .def_property(
          "plot_cache_valid",
          &TaoLatticeBranchStruct::plot_cache_valid,
          &TaoLatticeBranchStruct::set_plot_cache_valid,
          "Valid plotting data cache?"
      )
      .def_property(
          "spin_map_valid",
          &TaoLatticeBranchStruct::spin_map_valid,
          &TaoLatticeBranchStruct::set_spin_map_valid
      )
      .def_property(
          "twiss_valid",
          &TaoLatticeBranchStruct::twiss_valid,
          &TaoLatticeBranchStruct::set_twiss_valid,
          "Invalid EG with unstable 1-turn matrix with a closed branch. With open branch: "
          "twiss_valid = T even if some Twiss (and orbit) is invalid."
      )
      .def_property(
          "mode_flip_here",
          &TaoLatticeBranchStruct::mode_flip_here,
          &TaoLatticeBranchStruct::set_mode_flip_here,
          "Twiss parameter mode flip seen?"
      )
      .def_property(
          "chrom_calc_ok",
          &TaoLatticeBranchStruct::chrom_calc_ok,
          &TaoLatticeBranchStruct::set_chrom_calc_ok
      )
      .def_property(
          "rad_int_calc_ok",
          &TaoLatticeBranchStruct::rad_int_calc_ok,
          &TaoLatticeBranchStruct::set_rad_int_calc_ok
      )
      .def_property(
          "emit_6d_calc_ok",
          &TaoLatticeBranchStruct::emit_6d_calc_ok,
          &TaoLatticeBranchStruct::set_emit_6d_calc_ok
      )
      .def_property(
          "sigma_track_ok",
          &TaoLatticeBranchStruct::sigma_track_ok,
          &TaoLatticeBranchStruct::set_sigma_track_ok
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoLatticeBranchStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoLatticeBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoLatticeBranchStruct &self, const TaoLatticeBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoLatticeBranchStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoLatticeBranchStructArray1D, TaoLatticeBranchStructAlloc1D>(
      m,
      "TaoLatticeBranchStructArray1D",
      "TaoLatticeBranchStructAlloc1D"
  );
  // 2D TaoLatticeBranchStruct arrays are not used in structs/routines
  // 3D TaoLatticeBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lattice_struct
void init_tao_lattice_struct(py::module &m, py::class_<TaoLatticeStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             optional_ref<const LatStruct>,
             optional_ref<const LatStruct>,
             optional_ref<const LatStruct>,
             optional_ref<const RadIntAllEleStruct>,
             optional_ref<const RadIntAllEleStruct>>(),
         py::arg("name") = py::none(),
         py::arg("lat") = py::none(),
         py::arg("high_E_lat") = py::none(),
         py::arg("low_E_lat") = py::none(),
         py::arg("rad_int_by_ele_ri") = py::none(),
         py::arg("rad_int_by_ele_6d") = py::none()
  )
      .def_property(
          "name",
          &TaoLatticeStruct::name,
          &TaoLatticeStruct::set_name,
          "'model', 'base', or 'design'."
      )
      .def_property(
          "lat",
          py::cpp_function(&TaoLatticeStruct::lat, py::keep_alive<0, 1>()),
          &TaoLatticeStruct::set_lat,
          "lattice structures"
      )
      .def_property(
          "high_E_lat",
          py::cpp_function(&TaoLatticeStruct::high_E_lat, py::keep_alive<0, 1>()),
          &TaoLatticeStruct::set_high_E_lat,
          "For chrom calc."
      )
      .def_property(
          "low_E_lat",
          py::cpp_function(&TaoLatticeStruct::low_E_lat, py::keep_alive<0, 1>()),
          &TaoLatticeStruct::set_low_E_lat,
          "For chrom calc."
      )
      .def_property(
          "rad_int_by_ele_ri",
          py::cpp_function(&TaoLatticeStruct::rad_int_by_ele_ri, py::keep_alive<0, 1>()),
          &TaoLatticeStruct::set_rad_int_by_ele_ri
      )
      .def_property(
          "rad_int_by_ele_6d",
          py::cpp_function(&TaoLatticeStruct::rad_int_by_ele_6d, py::keep_alive<0, 1>()),
          &TaoLatticeStruct::set_rad_int_by_ele_6d
      )
      .def_property_readonly(
          "tao_branch",
          py::cpp_function(&TaoLatticeStruct::tao_branch, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const TaoLatticeStruct &self, const TaoLatticeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoLatticeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoLatticeStruct arrays are not used in structs/routines
  // 2D TaoLatticeStruct arrays are not used in structs/routines
  // 3D TaoLatticeStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_model_branch_struct
void init_tao_model_branch_struct(py::module &m, py::class_<TaoModelBranchStruct> &cls) {
  cls.def(py::init<optional_ref<const TaoBeamBranchStruct>>(), py::arg("beam") = py::none())
      .def_property_readonly(
          "ele",
          py::cpp_function(&TaoModelBranchStruct::ele, py::keep_alive<0, 1>()),
          "Per element information"
      )
      .def_property(
          "beam",
          py::cpp_function(&TaoModelBranchStruct::beam, py::keep_alive<0, 1>()),
          &TaoModelBranchStruct::set_beam
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoModelBranchStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoModelBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoModelBranchStruct &self, const TaoModelBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoModelBranchStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoModelBranchStructArray1D, TaoModelBranchStructAlloc1D>(
      m,
      "TaoModelBranchStructArray1D",
      "TaoModelBranchStructAlloc1D"
  );
  // 2D TaoModelBranchStruct arrays are not used in structs/routines
  // 3D TaoModelBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_model_element_struct
void init_tao_model_element_struct(py::module &m, py::class_<TaoModelElementStruct> &cls) {
  cls.def(
         py::init<optional_ref<const BeamStruct>, std::optional<bool>, std::optional<bool>>(),
         py::arg("beam") = py::none(),
         py::arg("save_beam_internally") = py::none(),
         py::arg("save_beam_to_file") = py::none()
  )
      .def_property(
          "beam",
          py::cpp_function(&TaoModelElementStruct::beam, py::keep_alive<0, 1>()),
          &TaoModelElementStruct::set_beam,
          "Beam distribution at element."
      )
      .def_property(
          "save_beam_internally",
          &TaoModelElementStruct::save_beam_internally,
          &TaoModelElementStruct::set_save_beam_internally,
          "Save beam here? Beam also saved at fork elements and at track ends."
      )
      .def_property(
          "save_beam_to_file",
          &TaoModelElementStruct::save_beam_to_file,
          &TaoModelElementStruct::set_save_beam_to_file,
          "Save beam to a file? Beam also saved at fork elements and at track ends."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoModelElementStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoModelElementStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoModelElementStruct &self, const TaoModelElementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoModelElementStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoModelElementStructArray1D, TaoModelElementStructAlloc1D>(
      m,
      "TaoModelElementStructArray1D",
      "TaoModelElementStructAlloc1D"
  );
  // 2D TaoModelElementStruct arrays are not used in structs/routines
  // 3D TaoModelElementStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ping_scale_struct
void init_tao_ping_scale_struct(py::module &m, py::class_<TaoPingScaleStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("a_mode_meas") = py::none(),
         py::arg("a_mode_ref") = py::none(),
         py::arg("b_mode_meas") = py::none(),
         py::arg("b_mode_ref") = py::none()
  )
      .def_property(
          "a_mode_meas",
          &TaoPingScaleStruct::a_mode_meas,
          &TaoPingScaleStruct::set_a_mode_meas
      )
      .def_property(
          "a_mode_ref",
          &TaoPingScaleStruct::a_mode_ref,
          &TaoPingScaleStruct::set_a_mode_ref
      )
      .def_property(
          "b_mode_meas",
          &TaoPingScaleStruct::b_mode_meas,
          &TaoPingScaleStruct::set_b_mode_meas
      )
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
      .def(
          "__eq__",
          [](const TaoPingScaleStruct &self, const TaoPingScaleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPingScaleStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoPingScaleStruct arrays are not used in structs/routines
  // 2D TaoPingScaleStruct arrays are not used in structs/routines
  // 3D TaoPingScaleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_cache_struct
void init_tao_plot_cache_struct(py::module &m, py::class_<TaoPlotCacheStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const EleStruct>,
             optional_ref<const CoordStruct>,
             std::optional<bool>>(),
         py::arg("ele_to_s") = py::none(),
         py::arg("orbit") = py::none(),
         py::arg("err") = py::none()
  )
      .def_property(
          "ele_to_s",
          py::cpp_function(&TaoPlotCacheStruct::ele_to_s, py::keep_alive<0, 1>()),
          &TaoPlotCacheStruct::set_ele_to_s,
          "Integrated element from branch beginning. Will be marked as a hybrid element."
      )
      .def_property(
          "orbit",
          py::cpp_function(&TaoPlotCacheStruct::orbit, py::keep_alive<0, 1>()),
          &TaoPlotCacheStruct::set_orbit
      )
      .def_property("err", &TaoPlotCacheStruct::err, &TaoPlotCacheStruct::set_err)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotCacheStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotCacheStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoPlotCacheStruct &self, const TaoPlotCacheStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotCacheStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoPlotCacheStructArray1D, TaoPlotCacheStructAlloc1D>(
      m,
      "TaoPlotCacheStructArray1D",
      "TaoPlotCacheStructAlloc1D"
  );
  // 2D TaoPlotCacheStruct arrays are not used in structs/routines
  // 3D TaoPlotCacheStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_page_struct
void init_tao_plot_page_struct(py::module &m, py::class_<TaoPlotPageStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const TaoTitleStruct>,
             optional_ref<const TaoTitleStruct>,
             optional_ref<const QpRectStruct>,
             optional_ref<const TaoDrawingStruct>,
             optional_ref<const TaoDrawingStruct>,
             std::optional<std::string>,
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
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("title") = py::none(),
         py::arg("subtitle") = py::none(),
         py::arg("border") = py::none(),
         py::arg("floor_plan") = py::none(),
         py::arg("lat_layout") = py::none(),
         py::arg("plot_display_type") = py::none(),
         py::arg("size") = py::none(),
         py::arg("text_height") = py::none(),
         py::arg("main_title_text_scale") = py::none(),
         py::arg("graph_title_text_scale") = py::none(),
         py::arg("axis_number_text_scale") = py::none(),
         py::arg("axis_label_text_scale") = py::none(),
         py::arg("legend_text_scale") = py::none(),
         py::arg("key_table_text_scale") = py::none(),
         py::arg("floor_plan_shape_scale") = py::none(),
         py::arg("floor_plan_text_scale") = py::none(),
         py::arg("lat_layout_shape_scale") = py::none(),
         py::arg("lat_layout_text_scale") = py::none(),
         py::arg("n_curve_pts") = py::none(),
         py::arg("id_window") = py::none(),
         py::arg("delete_overlapping_plots") = py::none(),
         py::arg("draw_graph_title_suffix") = py::none()
  )
      .def_property(
          "title",
          py::cpp_function(&TaoPlotPageStruct::title, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_title,
          "Title  at top of page."
      )
      .def_property(
          "subtitle",
          py::cpp_function(&TaoPlotPageStruct::subtitle, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_subtitle,
          "Subtitle below title at top of page."
      )
      .def_property(
          "border",
          py::cpp_function(&TaoPlotPageStruct::border, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_border,
          "Border around plots edge of page."
      )
      .def_property(
          "floor_plan",
          py::cpp_function(&TaoPlotPageStruct::floor_plan, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_floor_plan
      )
      .def_property(
          "lat_layout",
          py::cpp_function(&TaoPlotPageStruct::lat_layout, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_lat_layout
      )
      .def_property_readonly(
          "pattern",
          py::cpp_function(&TaoPlotPageStruct::pattern, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "template_",
          py::cpp_function(&TaoPlotPageStruct::template_, py::keep_alive<0, 1>()),
          "Templates for the plots."
      )
      .def_property_readonly(
          "region",
          py::cpp_function(&TaoPlotPageStruct::region, py::keep_alive<0, 1>())
      )
      .def_property(
          "plot_display_type",
          &TaoPlotPageStruct::plot_display_type,
          &TaoPlotPageStruct::set_plot_display_type,
          "'X' or 'TK'"
      )
      .def_property(
          "size",
          py::cpp_function(&TaoPlotPageStruct::size, py::keep_alive<0, 1>()),
          &TaoPlotPageStruct::set_size,
          "width and height of plot window in pixels."
      )
      .def_property(
          "text_height",
          &TaoPlotPageStruct::text_height,
          &TaoPlotPageStruct::set_text_height,
          "In points. Scales the height of all text"
      )
      .def_property(
          "main_title_text_scale",
          &TaoPlotPageStruct::main_title_text_scale,
          &TaoPlotPageStruct::set_main_title_text_scale,
          "Relative to text_height"
      )
      .def_property(
          "graph_title_text_scale",
          &TaoPlotPageStruct::graph_title_text_scale,
          &TaoPlotPageStruct::set_graph_title_text_scale,
          "Relative to text_height"
      )
      .def_property(
          "axis_number_text_scale",
          &TaoPlotPageStruct::axis_number_text_scale,
          &TaoPlotPageStruct::set_axis_number_text_scale,
          "Relative to text_height"
      )
      .def_property(
          "axis_label_text_scale",
          &TaoPlotPageStruct::axis_label_text_scale,
          &TaoPlotPageStruct::set_axis_label_text_scale,
          "Relative to text_height"
      )
      .def_property(
          "legend_text_scale",
          &TaoPlotPageStruct::legend_text_scale,
          &TaoPlotPageStruct::set_legend_text_scale,
          "Relative to text_height. For legends, plot_page, and lat_layout"
      )
      .def_property(
          "key_table_text_scale",
          &TaoPlotPageStruct::key_table_text_scale,
          &TaoPlotPageStruct::set_key_table_text_scale,
          "Relative to text_height"
      )
      .def_property(
          "floor_plan_shape_scale",
          &TaoPlotPageStruct::floor_plan_shape_scale,
          &TaoPlotPageStruct::set_floor_plan_shape_scale
      )
      .def_property(
          "floor_plan_text_scale",
          &TaoPlotPageStruct::floor_plan_text_scale,
          &TaoPlotPageStruct::set_floor_plan_text_scale,
          "Scale used = floor_plan_text_scale * legend_text_scale"
      )
      .def_property(
          "lat_layout_shape_scale",
          &TaoPlotPageStruct::lat_layout_shape_scale,
          &TaoPlotPageStruct::set_lat_layout_shape_scale
      )
      .def_property(
          "lat_layout_text_scale",
          &TaoPlotPageStruct::lat_layout_text_scale,
          &TaoPlotPageStruct::set_lat_layout_text_scale,
          "Scale used = lat_layout_text_scale * legend_text_scale"
      )
      .def_property(
          "n_curve_pts",
          &TaoPlotPageStruct::n_curve_pts,
          &TaoPlotPageStruct::set_n_curve_pts,
          "Default number of points for plotting a smooth curve."
      )
      .def_property(
          "id_window",
          &TaoPlotPageStruct::id_window,
          &TaoPlotPageStruct::set_id_window,
          "X window id number."
      )
      .def_property(
          "delete_overlapping_plots",
          &TaoPlotPageStruct::delete_overlapping_plots,
          &TaoPlotPageStruct::set_delete_overlapping_plots,
          "Delete overlapping plots when a plot is placed?"
      )
      .def_property(
          "draw_graph_title_suffix",
          &TaoPlotPageStruct::draw_graph_title_suffix,
          &TaoPlotPageStruct::set_draw_graph_title_suffix,
          "Draw the graph title suffix?"
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
      .def(
          "__eq__",
          [](const TaoPlotPageStruct &self, const TaoPlotPageStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotPageStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoPlotPageStruct arrays are not used in structs/routines
  // 2D TaoPlotPageStruct arrays are not used in structs/routines
  // 3D TaoPlotPageStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_region_struct
void init_tao_plot_region_struct(py::module &m, py::class_<TaoPlotRegionStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             optional_ref<const TaoPlotStruct>,
             std::optional<std::vector<double>>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("plot") = py::none(),
         py::arg("location") = py::none(),
         py::arg("visible") = py::none(),
         py::arg("list_with_show_plot_command") = py::none(),
         py::arg("setup_done") = py::none()
  )
      .def_property(
          "name",
          &TaoPlotRegionStruct::name,
          &TaoPlotRegionStruct::set_name,
          "Region name. Eg: 'r13', etc."
      )
      .def_property(
          "plot",
          py::cpp_function(&TaoPlotRegionStruct::plot, py::keep_alive<0, 1>()),
          &TaoPlotRegionStruct::set_plot,
          "Plot associated with this region"
      )
      .def_property(
          "location",
          py::cpp_function(&TaoPlotRegionStruct::location, py::keep_alive<0, 1>()),
          &TaoPlotRegionStruct::set_location,
          "[x1, x2, y1, y2] location on page."
      )
      .def_property(
          "visible",
          &TaoPlotRegionStruct::visible,
          &TaoPlotRegionStruct::set_visible,
          "To draw or not to draw."
      )
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotRegionStruct::list_with_show_plot_command,
          &TaoPlotRegionStruct::set_list_with_show_plot_command,
          "False used for default plots to shorten the output of 'show plot'"
      )
      .def_property(
          "setup_done",
          &TaoPlotRegionStruct::setup_done,
          &TaoPlotRegionStruct::set_setup_done,
          "Used for plot bookkeeping."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotRegionStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotRegionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoPlotRegionStruct &self, const TaoPlotRegionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotRegionStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoPlotRegionStructArray1D, TaoPlotRegionStructAlloc1D>(
      m,
      "TaoPlotRegionStructArray1D",
      "TaoPlotRegionStructAlloc1D"
  );
  // 2D TaoPlotRegionStruct arrays are not used in structs/routines
  // 3D TaoPlotRegionStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_struct
void init_tao_plot_struct(py::module &m, py::class_<TaoPlotStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             optional_ref<const TaoPlotRegionStruct>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("description") = py::none(),
         py::arg("r") = py::none(),
         py::arg("ix_plot") = py::none(),
         py::arg("n_curve_pts") = py::none(),
         py::arg("type") = py::none(),
         py::arg("x_axis_type") = py::none(),
         py::arg("autoscale_x") = py::none(),
         py::arg("autoscale_y") = py::none(),
         py::arg("autoscale_gang_x") = py::none(),
         py::arg("autoscale_gang_y") = py::none(),
         py::arg("list_with_show_plot_command") = py::none(),
         py::arg("phantom") = py::none(),
         py::arg("default_plot") = py::none()
  )
      .def_property(
          "name",
          &TaoPlotStruct::name,
          &TaoPlotStruct::set_name,
          "Identifying name. Rule: If name is blank, plot is not valid."
      )
      .def_property(
          "description",
          &TaoPlotStruct::description,
          &TaoPlotStruct::set_description,
          "Descriptive string."
      )
      .def_property_readonly(
          "graph",
          py::cpp_function(&TaoPlotStruct::graph, py::keep_alive<0, 1>()),
          "individual graphs of a plot"
      )
      .def_property(
          "r",
          py::cpp_function(&TaoPlotStruct::r, py::keep_alive<0, 1>()),
          &TaoPlotStruct::set_r,
          "pointer to parent."
      )
      .def_property(
          "ix_plot",
          &TaoPlotStruct::ix_plot,
          &TaoPlotStruct::set_ix_plot,
          "Index in s%plot_page%template(:) or %region(:) arrays."
      )
      .def_property(
          "n_curve_pts",
          &TaoPlotStruct::n_curve_pts,
          &TaoPlotStruct::set_n_curve_pts,
          "Overrides s%plot_page%n_curve_pts."
      )
      .def_property("type", &TaoPlotStruct::type, &TaoPlotStruct::set_type, "or 'wave'")
      .def_property(
          "x_axis_type",
          &TaoPlotStruct::x_axis_type,
          &TaoPlotStruct::set_x_axis_type,
          "'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc."
      )
      .def_property(
          "autoscale_x",
          &TaoPlotStruct::autoscale_x,
          &TaoPlotStruct::set_autoscale_x,
          "Horizontal autoscale."
      )
      .def_property(
          "autoscale_y",
          &TaoPlotStruct::autoscale_y,
          &TaoPlotStruct::set_autoscale_y,
          "Vertical autoscale."
      )
      .def_property(
          "autoscale_gang_x",
          &TaoPlotStruct::autoscale_gang_x,
          &TaoPlotStruct::set_autoscale_gang_x,
          "scale cmd scales graphs together?"
      )
      .def_property(
          "autoscale_gang_y",
          &TaoPlotStruct::autoscale_gang_y,
          &TaoPlotStruct::set_autoscale_gang_y,
          "scale cmd scales graphs together?"
      )
      .def_property(
          "list_with_show_plot_command",
          &TaoPlotStruct::list_with_show_plot_command,
          &TaoPlotStruct::set_list_with_show_plot_command,
          "False used for default plots to shorten the output of 'show plot'"
      )
      .def_property(
          "phantom",
          &TaoPlotStruct::phantom,
          &TaoPlotStruct::set_phantom,
          "Used by tao_plot_init to add info lines to 'show plot -templates'"
      )
      .def_property(
          "default_plot",
          &TaoPlotStruct::default_plot,
          &TaoPlotStruct::set_default_plot,
          "One of Tao's default plots?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoPlotStruct &self, const TaoPlotStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoPlotStructArray1D, TaoPlotStructAlloc1D>(
      m,
      "TaoPlotStructArray1D",
      "TaoPlotStructAlloc1D"
  );
  // 2D TaoPlotStruct arrays are not used in structs/routines
  // 3D TaoPlotStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_point_struct
void init_tao_shape_pattern_point_struct(
    py::module &m,
    py::class_<TaoShapePatternPointStruct> &cls
) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         py::arg("s") = py::none(),
         py::arg("y") = py::none(),
         py::arg("radius") = py::none()
  )
      .def_property("s", &TaoShapePatternPointStruct::s, &TaoShapePatternPointStruct::set_s)
      .def_property("y", &TaoShapePatternPointStruct::y, &TaoShapePatternPointStruct::set_y)
      .def_property(
          "radius",
          &TaoShapePatternPointStruct::radius,
          &TaoShapePatternPointStruct::set_radius
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoShapePatternPointStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoShapePatternPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoShapePatternPointStruct &self, const TaoShapePatternPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoShapePatternPointStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoShapePatternPointStructArray1D, TaoShapePatternPointStructAlloc1D>(
      m,
      "TaoShapePatternPointStructArray1D",
      "TaoShapePatternPointStructAlloc1D"
  );
  // 2D TaoShapePatternPointStruct arrays are not used in structs/routines
  // 3D TaoShapePatternPointStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_struct
void init_tao_shape_pattern_struct(py::module &m, py::class_<TaoShapePatternStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, optional_ref<const QpLineStruct>>(),
         py::arg("name") = py::none(),
         py::arg("line") = py::none()
  )
      .def_property("name", &TaoShapePatternStruct::name, &TaoShapePatternStruct::set_name)
      .def_property(
          "line",
          py::cpp_function(&TaoShapePatternStruct::line, py::keep_alive<0, 1>()),
          &TaoShapePatternStruct::set_line,
          "Line color and pattern set by shape using this pattern."
      )
      .def_property_readonly(
          "pt",
          py::cpp_function(&TaoShapePatternStruct::pt, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoShapePatternStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoShapePatternStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoShapePatternStruct &self, const TaoShapePatternStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoShapePatternStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoShapePatternStructArray1D, TaoShapePatternStructAlloc1D>(
      m,
      "TaoShapePatternStructArray1D",
      "TaoShapePatternStructAlloc1D"
  );
  // 2D TaoShapePatternStruct arrays are not used in structs/routines
  // 3D TaoShapePatternStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_dn_dpz_struct
void init_tao_spin_dn_dpz_struct(py::module &m, py::class_<TaoSpinDnDpzStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("vec") = py::none(),
         py::arg("partial") = py::none(),
         py::arg("partial2") = py::none()
  )
      .def_property(
          "vec",
          py::cpp_function(&TaoSpinDnDpzStruct::vec, py::keep_alive<0, 1>()),
          &TaoSpinDnDpzStruct::set_vec,
          "n0 derivative wrt pz."
      )
      .def_property(
          "partial",
          py::cpp_function(&TaoSpinDnDpzStruct::partial, py::keep_alive<0, 1>()),
          &TaoSpinDnDpzStruct::set_partial,
          "partial(i:) is spin n0 derivative wrt pz for i^th oscillation mode (1 => a-mode, etc.)"
      )
      .def_property(
          "partial2",
          py::cpp_function(&TaoSpinDnDpzStruct::partial2, py::keep_alive<0, 1>()),
          &TaoSpinDnDpzStruct::set_partial2,
          "partial(i:) is spin n0 derivative wrt pz with i^th oscillation mode missing (1 => "
          "a-mode, etc.)"
      )

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
      .def(
          "__eq__",
          [](const TaoSpinDnDpzStruct &self, const TaoSpinDnDpzStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoSpinDnDpzStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoSpinDnDpzStruct arrays are not used in structs/routines
  // 2D TaoSpinDnDpzStruct arrays are not used in structs/routines
  // 3D TaoSpinDnDpzStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_ele_struct
void init_tao_spin_ele_struct(py::module &m, py::class_<TaoSpinEleStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const TaoSpinDnDpzStruct>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<bool>>(),
         py::arg("dn_dpz") = py::none(),
         py::arg("orb_eigen_val") = py::none(),
         py::arg("orb_eigen_vec") = py::none(),
         py::arg("spin_eigen_vec") = py::none(),
         py::arg("valid") = py::none()
  )
      .def_property(
          "dn_dpz",
          py::cpp_function(&TaoSpinEleStruct::dn_dpz, py::keep_alive<0, 1>()),
          &TaoSpinEleStruct::set_dn_dpz
      )
      .def_property(
          "orb_eigen_val",
          py::cpp_function(&TaoSpinEleStruct::orb_eigen_val, py::keep_alive<0, 1>()),
          &TaoSpinEleStruct::set_orb_eigen_val
      )
      .def_property(
          "orb_eigen_vec",
          py::cpp_function(&TaoSpinEleStruct::orb_eigen_vec, py::keep_alive<0, 1>()),
          &TaoSpinEleStruct::set_orb_eigen_vec,
          "(j,:) is j^th vector"
      )
      .def_property(
          "spin_eigen_vec",
          py::cpp_function(&TaoSpinEleStruct::spin_eigen_vec, py::keep_alive<0, 1>()),
          &TaoSpinEleStruct::set_spin_eigen_vec,
          "(j,:) is j^th vector"
      )
      .def_property("valid", &TaoSpinEleStruct::valid, &TaoSpinEleStruct::set_valid)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoSpinEleStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoSpinEleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoSpinEleStruct &self, const TaoSpinEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoSpinEleStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoSpinEleStructArray1D, TaoSpinEleStructAlloc1D>(
      m,
      "TaoSpinEleStructArray1D",
      "TaoSpinEleStructAlloc1D"
  );
  // 2D TaoSpinEleStruct arrays are not used in structs/routines
  // 3D TaoSpinEleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_map_struct
void init_tao_spin_map_struct(py::module &m, py::class_<TaoSpinMapStruct> &cls) {
  cls.def(
         py::init<
             std::optional<bool>,
             optional_ref<const SpinOrbitMap1Struct>,
             optional_ref<const SpinAxisStruct>,
             optional_ref<const SpinAxisStruct>,
             optional_ref<const SpinAxisStruct>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("valid") = py::none(),
         py::arg("map1") = py::none(),
         py::arg("axis_input") = py::none(),
         py::arg("axis0") = py::none(),
         py::arg("axis1") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_ref") = py::none(),
         py::arg("ix_uni") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("mat8") = py::none()
  )
      .def_property("valid", &TaoSpinMapStruct::valid, &TaoSpinMapStruct::set_valid)
      .def_property(
          "map1",
          py::cpp_function(&TaoSpinMapStruct::map1, py::keep_alive<0, 1>()),
          &TaoSpinMapStruct::set_map1
      )
      .def_property(
          "axis_input",
          py::cpp_function(&TaoSpinMapStruct::axis_input, py::keep_alive<0, 1>()),
          &TaoSpinMapStruct::set_axis_input,
          "Input axes."
      )
      .def_property(
          "axis0",
          py::cpp_function(&TaoSpinMapStruct::axis0, py::keep_alive<0, 1>()),
          &TaoSpinMapStruct::set_axis0,
          "Initial axes."
      )
      .def_property(
          "axis1",
          py::cpp_function(&TaoSpinMapStruct::axis1, py::keep_alive<0, 1>()),
          &TaoSpinMapStruct::set_axis1,
          "Final axes."
      )
      .def_property("ix_ele", &TaoSpinMapStruct::ix_ele, &TaoSpinMapStruct::set_ix_ele)
      .def_property("ix_ref", &TaoSpinMapStruct::ix_ref, &TaoSpinMapStruct::set_ix_ref)
      .def_property("ix_uni", &TaoSpinMapStruct::ix_uni, &TaoSpinMapStruct::set_ix_uni)
      .def_property("ix_branch", &TaoSpinMapStruct::ix_branch, &TaoSpinMapStruct::set_ix_branch)
      .def_property(
          "mat8",
          py::cpp_function(&TaoSpinMapStruct::mat8, py::keep_alive<0, 1>()),
          &TaoSpinMapStruct::set_mat8
      )

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
      .def(
          "__eq__",
          [](const TaoSpinMapStruct &self, const TaoSpinMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoSpinMapStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoSpinMapStruct arrays are not used in structs/routines
  // 2D TaoSpinMapStruct arrays are not used in structs/routines
  // 3D TaoSpinMapStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_spin_polarization_struct
void init_tao_spin_polarization_struct(py::module &m, py::class_<TaoSpinPolarizationStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             optional_ref<const SpinOrbitMap1Struct>>(),
         py::arg("tune") = py::none(),
         py::arg("pol_limit_st") = py::none(),
         py::arg("pol_limit_dk") = py::none(),
         py::arg("pol_limit_dk_partial") = py::none(),
         py::arg("pol_limit_dk_partial2") = py::none(),
         py::arg("pol_rate_bks") = py::none(),
         py::arg("depol_rate") = py::none(),
         py::arg("depol_rate_partial") = py::none(),
         py::arg("depol_rate_partial2") = py::none(),
         py::arg("integral_bn") = py::none(),
         py::arg("integral_bdn") = py::none(),
         py::arg("integral_1ns") = py::none(),
         py::arg("integral_dn2") = py::none(),
         py::arg("valid") = py::none(),
         py::arg("q_1turn") = py::none()
  )
      .def_property("tune", &TaoSpinPolarizationStruct::tune, &TaoSpinPolarizationStruct::set_tune)
      .def_property(
          "pol_limit_st",
          &TaoSpinPolarizationStruct::pol_limit_st,
          &TaoSpinPolarizationStruct::set_pol_limit_st,
          "Polarization calculated using Sokolov-Ternov formula."
      )
      .def_property(
          "pol_limit_dk",
          &TaoSpinPolarizationStruct::pol_limit_dk,
          &TaoSpinPolarizationStruct::set_pol_limit_dk,
          "Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula."
      )
      .def_property(
          "pol_limit_dk_partial",
          py::cpp_function(
              &TaoSpinPolarizationStruct::pol_limit_dk_partial,
              py::keep_alive<0, 1>()
          ),
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial,
          "Limit using only single mode to calc dn_dpz"
      )
      .def_property(
          "pol_limit_dk_partial2",
          py::cpp_function(
              &TaoSpinPolarizationStruct::pol_limit_dk_partial2,
              py::keep_alive<0, 1>()
          ),
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial2,
          "Limit using only single mode to calc dn_dpz"
      )
      .def_property(
          "pol_rate_bks",
          &TaoSpinPolarizationStruct::pol_rate_bks,
          &TaoSpinPolarizationStruct::set_pol_rate_bks,
          "BKS Polarization rate (1/sec)."
      )
      .def_property(
          "depol_rate",
          &TaoSpinPolarizationStruct::depol_rate,
          &TaoSpinPolarizationStruct::set_depol_rate,
          "Depolarization rate (1/sec)."
      )
      .def_property(
          "depol_rate_partial",
          py::cpp_function(&TaoSpinPolarizationStruct::depol_rate_partial, py::keep_alive<0, 1>()),
          &TaoSpinPolarizationStruct::set_depol_rate_partial,
          "Depolarization rate (1/sec) using only single mode to calc dn_dpz."
      )
      .def_property(
          "depol_rate_partial2",
          py::cpp_function(&TaoSpinPolarizationStruct::depol_rate_partial2, py::keep_alive<0, 1>()),
          &TaoSpinPolarizationStruct::set_depol_rate_partial2,
          "Depolarization rate (1/sec) using only two modes to calc dn_dpz."
      )
      .def_property(
          "integral_bn",
          &TaoSpinPolarizationStruct::integral_bn,
          &TaoSpinPolarizationStruct::set_integral_bn,
          "Integral of g^3 * b_hat * n_0"
      )
      .def_property(
          "integral_bdn",
          &TaoSpinPolarizationStruct::integral_bdn,
          &TaoSpinPolarizationStruct::set_integral_bdn,
          "Integral of g^3 * b_hat * dn/ddelta"
      )
      .def_property(
          "integral_1ns",
          &TaoSpinPolarizationStruct::integral_1ns,
          &TaoSpinPolarizationStruct::set_integral_1ns,
          "Integral of g^3 (1 - 2(n * s_hat)/9)"
      )
      .def_property(
          "integral_dn2",
          &TaoSpinPolarizationStruct::integral_dn2,
          &TaoSpinPolarizationStruct::set_integral_dn2,
          "Integral of g^3 * 11 (dn/ddelta)^2 / 9"
      )
      .def_property(
          "valid",
          &TaoSpinPolarizationStruct::valid,
          &TaoSpinPolarizationStruct::set_valid
      )
      .def_property(
          "q_1turn",
          py::cpp_function(&TaoSpinPolarizationStruct::q_1turn, py::keep_alive<0, 1>()),
          &TaoSpinPolarizationStruct::set_q_1turn,
          "Save results from spin_concat_linear_maps in tao_spin_polarization."
      )
      .def_property_readonly(
          "q_ele",
          py::cpp_function(&TaoSpinPolarizationStruct::q_ele, py::keep_alive<0, 1>()),
          "Save results from spin_concat_linear_maps in tao_spin_polarization."
      )

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
      .def(
          "__eq__",
          [](const TaoSpinPolarizationStruct &self, const TaoSpinPolarizationStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoSpinPolarizationStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             optional_ref<const TaoGlobalStruct>,
             optional_ref<const TaoInitStruct>,
             optional_ref<const TaoCommonStruct>,
             optional_ref<const TaoPlotPageStruct>,
             std::optional<std::vector<int>>,
             optional_ref<const TaoBuildingWallStruct>,
             optional_ref<const TaoWaveStruct>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("global_") = py::none(),
         py::arg("init") = py::none(),
         py::arg("com") = py::none(),
         py::arg("plot_page") = py::none(),
         py::arg("key") = py::none(),
         py::arg("building_wall") = py::none(),
         py::arg("wave") = py::none(),
         py::arg("n_var_used") = py::none(),
         py::arg("n_v1_var_used") = py::none(),
         py::arg("initialized") = py::none()
  )
      .def_property(
          "global_",
          py::cpp_function(&TaoSuperUniverseStruct::global, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_global,
          "User accessible global variables."
      )
      .def_property(
          "init",
          py::cpp_function(&TaoSuperUniverseStruct::init, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_init,
          "Initialization parameters"
      )
      .def_property(
          "com",
          py::cpp_function(&TaoSuperUniverseStruct::com, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_com,
          "Non-initialization common parameters"
      )
      .def_property(
          "plot_page",
          py::cpp_function(&TaoSuperUniverseStruct::plot_page, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_plot_page,
          "Defines the plot window."
      )
      .def_property_readonly(
          "v1_var",
          py::cpp_function(&TaoSuperUniverseStruct::v1_var, py::keep_alive<0, 1>()),
          "The variable types"
      )
      .def_property_readonly(
          "var",
          py::cpp_function(&TaoSuperUniverseStruct::var, py::keep_alive<0, 1>()),
          "array of all variables."
      )
      .def_property_readonly(
          "u",
          py::cpp_function(&TaoSuperUniverseStruct::u, py::keep_alive<0, 1>()),
          "array of universes."
      )
      .def_property(
          "key",
          py::cpp_function(&TaoSuperUniverseStruct::key, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_key
      )
      .def_property(
          "building_wall",
          py::cpp_function(&TaoSuperUniverseStruct::building_wall, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_building_wall
      )
      .def_property(
          "wave",
          py::cpp_function(&TaoSuperUniverseStruct::wave, py::keep_alive<0, 1>()),
          &TaoSuperUniverseStruct::set_wave
      )
      .def_property(
          "n_var_used",
          &TaoSuperUniverseStruct::n_var_used,
          &TaoSuperUniverseStruct::set_n_var_used
      )
      .def_property(
          "n_v1_var_used",
          &TaoSuperUniverseStruct::n_v1_var_used,
          &TaoSuperUniverseStruct::set_n_v1_var_used
      )
      .def_property_readonly(
          "history",
          py::cpp_function(&TaoSuperUniverseStruct::history, py::keep_alive<0, 1>()),
          "command history"
      )
      .def_property(
          "initialized",
          &TaoSuperUniverseStruct::initialized,
          &TaoSuperUniverseStruct::set_initialized,
          "Does tao_init() need to be called?"
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
      .def(
          "__eq__",
          [](const TaoSuperUniverseStruct &self, const TaoSuperUniverseStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoSuperUniverseStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>>(),
         py::arg("string") = py::none(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("units") = py::none(),
         py::arg("justify") = py::none(),
         py::arg("draw_it") = py::none()
  )
      .def_property(
          "string",
          &TaoTitleStruct::string,
          &TaoTitleStruct::set_string,
          "title character string."
      )
      .def_property("x", &TaoTitleStruct::x, &TaoTitleStruct::set_x, "x, y rwt lower left corner")
      .def_property("y", &TaoTitleStruct::y, &TaoTitleStruct::set_y, "x, y rwt lower left corner")
      .def_property(
          "units",
          &TaoTitleStruct::units,
          &TaoTitleStruct::set_units,
          "%BOX, POINTS, etc..."
      )
      .def_property(
          "justify",
          &TaoTitleStruct::justify,
          &TaoTitleStruct::set_justify,
          "Left, Center, or Right justification."
      )
      .def_property(
          "draw_it",
          &TaoTitleStruct::draw_it,
          &TaoTitleStruct::set_draw_it,
          "draw the title?"
      )

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
      .def(
          "__eq__",
          [](const TaoTitleStruct &self, const TaoTitleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoTitleStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoTitleStruct arrays are not used in structs/routines
  // 2D TaoTitleStruct arrays are not used in structs/routines
  // 3D TaoTitleStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_calc_struct
void init_tao_universe_calc_struct(py::module &m, py::class_<TaoUniverseCalcStruct> &cls) {
  cls.def(
         py::init<
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
             std::optional<bool>>(),
         py::arg("srdt_for_data") = py::none(),
         py::arg("rad_int_for_data") = py::none(),
         py::arg("rad_int_for_plotting") = py::none(),
         py::arg("chrom_for_data") = py::none(),
         py::arg("chrom_for_plotting") = py::none(),
         py::arg("lat_sigma_for_data") = py::none(),
         py::arg("lat_sigma_for_plotting") = py::none(),
         py::arg("dynamic_aperture") = py::none(),
         py::arg("one_turn_map") = py::none(),
         py::arg("lattice") = py::none(),
         py::arg("twiss") = py::none(),
         py::arg("track") = py::none(),
         py::arg("spin_matrices") = py::none()
  )
      .def_property(
          "srdt_for_data",
          &TaoUniverseCalcStruct::srdt_for_data,
          &TaoUniverseCalcStruct::set_srdt_for_data,
          "0 = false, 1 = 1st order, 2 = 1st & 2nd order"
      )
      .def_property(
          "rad_int_for_data",
          &TaoUniverseCalcStruct::rad_int_for_data,
          &TaoUniverseCalcStruct::set_rad_int_for_data,
          "Do the radiation integrals need to be computed for"
      )
      .def_property(
          "rad_int_for_plotting",
          &TaoUniverseCalcStruct::rad_int_for_plotting,
          &TaoUniverseCalcStruct::set_rad_int_for_plotting,
          "data or plotting?"
      )
      .def_property(
          "chrom_for_data",
          &TaoUniverseCalcStruct::chrom_for_data,
          &TaoUniverseCalcStruct::set_chrom_for_data,
          "Does the chromaticity need to be computed for"
      )
      .def_property(
          "chrom_for_plotting",
          &TaoUniverseCalcStruct::chrom_for_plotting,
          &TaoUniverseCalcStruct::set_chrom_for_plotting,
          "data or plotting?"
      )
      .def_property(
          "lat_sigma_for_data",
          &TaoUniverseCalcStruct::lat_sigma_for_data,
          &TaoUniverseCalcStruct::set_lat_sigma_for_data,
          "Do the beam sigmas need to be computed for"
      )
      .def_property(
          "lat_sigma_for_plotting",
          &TaoUniverseCalcStruct::lat_sigma_for_plotting,
          &TaoUniverseCalcStruct::set_lat_sigma_for_plotting,
          "data or plotting?"
      )
      .def_property(
          "dynamic_aperture",
          &TaoUniverseCalcStruct::dynamic_aperture,
          &TaoUniverseCalcStruct::set_dynamic_aperture,
          "Do the dynamic_aperture calc?"
      )
      .def_property(
          "one_turn_map",
          &TaoUniverseCalcStruct::one_turn_map,
          &TaoUniverseCalcStruct::set_one_turn_map,
          "Compute the one turn map?"
      )
      .def_property(
          "lattice",
          &TaoUniverseCalcStruct::lattice,
          &TaoUniverseCalcStruct::set_lattice,
          "Used to indicate which lattices need tracking done."
      )
      .def_property(
          "twiss",
          &TaoUniverseCalcStruct::twiss,
          &TaoUniverseCalcStruct::set_twiss,
          "calc linear transfer matrix?"
      )
      .def_property(
          "track",
          &TaoUniverseCalcStruct::track,
          &TaoUniverseCalcStruct::set_track,
          "tracking needs to be done?"
      )
      .def_property(
          "spin_matrices",
          &TaoUniverseCalcStruct::spin_matrices,
          &TaoUniverseCalcStruct::set_spin_matrices,
          "Calculate G and D spin matrices?"
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
      .def(
          "__eq__",
          [](const TaoUniverseCalcStruct &self, const TaoUniverseCalcStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoUniverseCalcStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(py::init<optional_ref<const TaoUniverseStruct>>(), py::arg("u") = py::none())
      .def_property(
          "u",
          py::cpp_function(&TaoUniversePointerStruct::u, py::keep_alive<0, 1>()),
          &TaoUniversePointerStruct::set_u
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoUniversePointerStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoUniversePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoUniversePointerStruct &self, const TaoUniversePointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoUniversePointerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoUniversePointerStructArray1D, TaoUniversePointerStructAlloc1D>(
      m,
      "TaoUniversePointerStructArray1D",
      "TaoUniversePointerStructAlloc1D"
  );
  // 2D TaoUniversePointerStruct arrays are not used in structs/routines
  // 3D TaoUniversePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_universe_struct
void init_tao_universe_struct(py::module &m, py::class_<TaoUniverseStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const TaoLatticeStruct>,
             optional_ref<const TaoLatticeStruct>,
             optional_ref<const TaoLatticeStruct>,
             optional_ref<const TaoBeamUniStruct>,
             optional_ref<const TaoDynamicApertureStruct>,
             optional_ref<const TaoPingScaleStruct>,
             optional_ref<const LatStruct>,
             optional_ref<const TaoUniverseCalcStruct>,
             optional_ref<const LatEleOrderStruct>,
             optional_ref<const TaoSpinMapStruct>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("model") = py::none(),
         py::arg("design") = py::none(),
         py::arg("base") = py::none(),
         py::arg("beam") = py::none(),
         py::arg("dynamic_aperture") = py::none(),
         py::arg("ping_scale") = py::none(),
         py::arg("scratch_lat") = py::none(),
         py::arg("calc") = py::none(),
         py::arg("ele_order") = py::none(),
         py::arg("spin_map") = py::none(),
         py::arg("dModel_dVar") = py::none(),
         py::arg("ix_uni") = py::none(),
         py::arg("n_d2_data_used") = py::none(),
         py::arg("n_data_used") = py::none(),
         py::arg("is_on") = py::none(),
         py::arg("design_same_as_previous") = py::none(),
         py::arg("picked_uni") = py::none()
  )
      .def_property(
          "model",
          py::cpp_function(&TaoUniverseStruct::model, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_model
      )
      .def_property(
          "design",
          py::cpp_function(&TaoUniverseStruct::design, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_design
      )
      .def_property(
          "base",
          py::cpp_function(&TaoUniverseStruct::base, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_base
      )
      .def_property(
          "beam",
          py::cpp_function(&TaoUniverseStruct::beam, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_beam
      )
      .def_property(
          "dynamic_aperture",
          py::cpp_function(&TaoUniverseStruct::dynamic_aperture, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_dynamic_aperture
      )
      .def_property_readonly(
          "model_branch",
          py::cpp_function(&TaoUniverseStruct::model_branch, py::keep_alive<0, 1>()),
          "model specific information"
      )
      .def_property_readonly(
          "d2_data",
          py::cpp_function(&TaoUniverseStruct::d2_data, py::keep_alive<0, 1>()),
          "The data types"
      )
      .def_property_readonly(
          "data",
          py::cpp_function(&TaoUniverseStruct::data, py::keep_alive<0, 1>()),
          "Array of all data."
      )
      .def_property(
          "ping_scale",
          py::cpp_function(&TaoUniverseStruct::ping_scale, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_ping_scale
      )
      .def_property(
          "scratch_lat",
          py::cpp_function(&TaoUniverseStruct::scratch_lat, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_scratch_lat,
          "Scratch area."
      )
      .def_property(
          "calc",
          py::cpp_function(&TaoUniverseStruct::calc, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_calc,
          "What needs to be calculated?"
      )
      .def_property(
          "ele_order",
          py::cpp_function(&TaoUniverseStruct::ele_order, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_ele_order,
          "Order of elements with same name."
      )
      .def_property(
          "spin_map",
          py::cpp_function(&TaoUniverseStruct::spin_map, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_spin_map
      )
      .def_property(
          "dModel_dVar",
          py::cpp_function(&TaoUniverseStruct::dModel_dVar, py::keep_alive<0, 1>()),
          &TaoUniverseStruct::set_dModel_dVar,
          "Derivative matrix."
      )
      .def_property(
          "ix_uni",
          &TaoUniverseStruct::ix_uni,
          &TaoUniverseStruct::set_ix_uni,
          "Universe index."
      )
      .def_property(
          "n_d2_data_used",
          &TaoUniverseStruct::n_d2_data_used,
          &TaoUniverseStruct::set_n_d2_data_used,
          "Number of used %d2_data(:) components."
      )
      .def_property(
          "n_data_used",
          &TaoUniverseStruct::n_data_used,
          &TaoUniverseStruct::set_n_data_used,
          "Number of used %data(:) components."
      )
      .def_property(
          "is_on",
          &TaoUniverseStruct::is_on,
          &TaoUniverseStruct::set_is_on,
          "universe turned on"
      )
      .def_property(
          "design_same_as_previous",
          &TaoUniverseStruct::design_same_as_previous,
          &TaoUniverseStruct::set_design_same_as_previous,
          "Design lat same as the previous uni?"
      )
      .def_property(
          "picked_uni",
          &TaoUniverseStruct::picked_uni,
          &TaoUniverseStruct::set_picked_uni,
          "Scratch logical."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoUniverseStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoUniverseStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoUniverseStruct &self, const TaoUniverseStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoUniverseStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoUniverseStructArray1D, TaoUniverseStructAlloc1D>(
      m,
      "TaoUniverseStructArray1D",
      "TaoUniverseStructAlloc1D"
  );
  // 2D TaoUniverseStruct arrays are not used in structs/routines
  // 3D TaoUniverseStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_v1_var_struct
void init_tao_v1_var_struct(py::module &m, py::class_<TaoV1VarStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>>(),
         py::arg("name") = py::none(),
         py::arg("ix_v1_var") = py::none()
  )
      .def_property(
          "name",
          &TaoV1VarStruct::name,
          &TaoV1VarStruct::set_name,
          "V1 variable name. Eg: 'quad_k1'."
      )
      .def_property(
          "ix_v1_var",
          &TaoV1VarStruct::ix_v1_var,
          &TaoV1VarStruct::set_ix_v1_var,
          "Index to s%v1_var(:) array"
      )
      .def_property_readonly(
          "v",
          py::cpp_function(&TaoV1VarStruct::v, py::keep_alive<0, 1>()),
          "Pointer to the appropriate section in s%var."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoV1VarStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoV1VarStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoV1VarStruct &self, const TaoV1VarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoV1VarStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoV1VarStructArray1D, TaoV1VarStructAlloc1D>(
      m,
      "TaoV1VarStructArray1D",
      "TaoV1VarStructAlloc1D"
  );
  // 2D TaoV1VarStruct arrays are not used in structs/routines
  // 3D TaoV1VarStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_var_slave_struct
void init_tao_var_slave_struct(py::module &m, py::class_<TaoVarSlaveStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("ix_uni") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("model_value") = py::none(),
         py::arg("base_value") = py::none()
  )
      .def_property(
          "ix_uni",
          &TaoVarSlaveStruct::ix_uni,
          &TaoVarSlaveStruct::set_ix_uni,
          "universe index."
      )
      .def_property("ix_branch", &TaoVarSlaveStruct::ix_branch, &TaoVarSlaveStruct::set_ix_branch)
      .def_property(
          "ix_ele",
          &TaoVarSlaveStruct::ix_ele,
          &TaoVarSlaveStruct::set_ix_ele,
          "Index of element in the u%lattice%ele(:) array."
      )
      .def_property(
          "model_value",
          &TaoVarSlaveStruct::model_value,
          &TaoVarSlaveStruct::set_model_value,
          "Pointer to the variable in the model lat."
      )
      .def_property(
          "base_value",
          &TaoVarSlaveStruct::base_value,
          &TaoVarSlaveStruct::set_base_value,
          "Pointer to the variable in the base lat."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoVarSlaveStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoVarSlaveStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoVarSlaveStruct &self, const TaoVarSlaveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoVarSlaveStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoVarSlaveStructArray1D, TaoVarSlaveStructAlloc1D>(
      m,
      "TaoVarSlaveStructArray1D",
      "TaoVarSlaveStructAlloc1D"
  );
  // 2D TaoVarSlaveStruct arrays are not used in structs/routines
  // 3D TaoVarSlaveStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_var_struct
void init_tao_var_struct(py::module &m, py::class_<TaoVarStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
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
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             optional_ref<const TaoV1VarStruct>>(),
         py::arg("ele_name") = py::none(),
         py::arg("attrib_name") = py::none(),
         py::arg("id") = py::none(),
         py::arg("ix_v1") = py::none(),
         py::arg("ix_var") = py::none(),
         py::arg("ix_dvar") = py::none(),
         py::arg("ix_attrib") = py::none(),
         py::arg("ix_key_table") = py::none(),
         py::arg("model_value") = py::none(),
         py::arg("base_value") = py::none(),
         py::arg("design_value") = py::none(),
         py::arg("scratch_value") = py::none(),
         py::arg("old_value") = py::none(),
         py::arg("meas_value") = py::none(),
         py::arg("ref_value") = py::none(),
         py::arg("correction_value") = py::none(),
         py::arg("high_lim") = py::none(),
         py::arg("low_lim") = py::none(),
         py::arg("step") = py::none(),
         py::arg("weight") = py::none(),
         py::arg("delta_merit") = py::none(),
         py::arg("merit") = py::none(),
         py::arg("dMerit_dVar") = py::none(),
         py::arg("key_val0") = py::none(),
         py::arg("key_delta") = py::none(),
         py::arg("s") = py::none(),
         py::arg("extend_val") = py::none(),
         py::arg("merit_type") = py::none(),
         py::arg("exists") = py::none(),
         py::arg("good_var") = py::none(),
         py::arg("good_user") = py::none(),
         py::arg("good_opt") = py::none(),
         py::arg("good_plot") = py::none(),
         py::arg("useit_opt") = py::none(),
         py::arg("useit_plot") = py::none(),
         py::arg("key_bound") = py::none(),
         py::arg("v1") = py::none()
  )
      .def_property(
          "ele_name",
          &TaoVarStruct::ele_name,
          &TaoVarStruct::set_ele_name,
          "Associated lattice element name."
      )
      .def_property(
          "attrib_name",
          &TaoVarStruct::attrib_name,
          &TaoVarStruct::set_attrib_name,
          "Name of the attribute to vary."
      )
      .def_property(
          "id",
          &TaoVarStruct::id,
          &TaoVarStruct::set_id,
          "Used by Tao extension code. Not used by Tao directly."
      )
      .def_property_readonly(
          "slave",
          py::cpp_function(&TaoVarStruct::slave, py::keep_alive<0, 1>())
      )
      .def_property(
          "ix_v1",
          &TaoVarStruct::ix_v1,
          &TaoVarStruct::set_ix_v1,
          "Index of this var in the s%v1_var(i)%v(:) array."
      )
      .def_property(
          "ix_var",
          &TaoVarStruct::ix_var,
          &TaoVarStruct::set_ix_var,
          "Index number of this var in the s%var(:) array."
      )
      .def_property(
          "ix_dvar",
          &TaoVarStruct::ix_dvar,
          &TaoVarStruct::set_ix_dvar,
          "Column in the dData_dVar derivative matrix."
      )
      .def_property(
          "ix_attrib",
          &TaoVarStruct::ix_attrib,
          &TaoVarStruct::set_ix_attrib,
          "Index in ele%value(:) array if appropriate."
      )
      .def_property(
          "ix_key_table",
          &TaoVarStruct::ix_key_table,
          &TaoVarStruct::set_ix_key_table,
          "Has a key binding?"
      )
      .def_property(
          "model_value",
          &TaoVarStruct::model_value,
          &TaoVarStruct::set_model_value,
          "Model value."
      )
      .def_property(
          "base_value",
          &TaoVarStruct::base_value,
          &TaoVarStruct::set_base_value,
          "Base value."
      )
      .def_property(
          "design_value",
          &TaoVarStruct::design_value,
          &TaoVarStruct::set_design_value,
          "Design value from the design lattice."
      )
      .def_property(
          "scratch_value",
          &TaoVarStruct::scratch_value,
          &TaoVarStruct::set_scratch_value,
          "Scratch space used by Tao."
      )
      .def_property(
          "old_value",
          &TaoVarStruct::old_value,
          &TaoVarStruct::set_old_value,
          "Scratch space used by Tao."
      )
      .def_property(
          "meas_value",
          &TaoVarStruct::meas_value,
          &TaoVarStruct::set_meas_value,
          "The value when the data measurement was taken."
      )
      .def_property(
          "ref_value",
          &TaoVarStruct::ref_value,
          &TaoVarStruct::set_ref_value,
          "Value when the reference measurement was taken."
      )
      .def_property(
          "correction_value",
          &TaoVarStruct::correction_value,
          &TaoVarStruct::set_correction_value,
          "Value determined by a fit to correct the lattice."
      )
      .def_property(
          "high_lim",
          &TaoVarStruct::high_lim,
          &TaoVarStruct::set_high_lim,
          "High limit for the model_value."
      )
      .def_property(
          "low_lim",
          &TaoVarStruct::low_lim,
          &TaoVarStruct::set_low_lim,
          "Low limit for the model_value."
      )
      .def_property(
          "step",
          &TaoVarStruct::step,
          &TaoVarStruct::set_step,
          "Sets what is a small step for varying this var."
      )
      .def_property(
          "weight",
          &TaoVarStruct::weight,
          &TaoVarStruct::set_weight,
          "Weight for the merit function term."
      )
      .def_property(
          "delta_merit",
          &TaoVarStruct::delta_merit,
          &TaoVarStruct::set_delta_merit,
          "Diff used to calculate the merit function term."
      )
      .def_property(
          "merit",
          &TaoVarStruct::merit,
          &TaoVarStruct::set_merit,
          "merit_term = weight * delta^2."
      )
      .def_property(
          "dMerit_dVar",
          &TaoVarStruct::dMerit_dVar,
          &TaoVarStruct::set_dMerit_dVar,
          "Merit derivative."
      )
      .def_property(
          "key_val0",
          &TaoVarStruct::key_val0,
          &TaoVarStruct::set_key_val0,
          "Key base value"
      )
      .def_property(
          "key_delta",
          &TaoVarStruct::key_delta,
          &TaoVarStruct::set_key_delta,
          "Change in value when a key is pressed."
      )
      .def_property("s", &TaoVarStruct::s, &TaoVarStruct::set_s, "longitudinal position of ele.")
      .def_property(
          "extend_val",
          &TaoVarStruct::extend_val,
          &TaoVarStruct::set_extend_val,
          "For extension code. Not used by Tao."
      )
      .def_property(
          "merit_type",
          &TaoVarStruct::merit_type,
          &TaoVarStruct::set_merit_type,
          "'target' or 'limit'"
      )
      .def_property("exists", &TaoVarStruct::exists, &TaoVarStruct::set_exists, "See above")
      .def_property("good_var", &TaoVarStruct::good_var, &TaoVarStruct::set_good_var, "See above")
      .def_property(
          "good_user",
          &TaoVarStruct::good_user,
          &TaoVarStruct::set_good_user,
          "See above"
      )
      .def_property("good_opt", &TaoVarStruct::good_opt, &TaoVarStruct::set_good_opt, "See above")
      .def_property(
          "good_plot",
          &TaoVarStruct::good_plot,
          &TaoVarStruct::set_good_plot,
          "See above"
      )
      .def_property(
          "useit_opt",
          &TaoVarStruct::useit_opt,
          &TaoVarStruct::set_useit_opt,
          "See above"
      )
      .def_property(
          "useit_plot",
          &TaoVarStruct::useit_plot,
          &TaoVarStruct::set_useit_plot,
          "See above"
      )
      .def_property(
          "key_bound",
          &TaoVarStruct::key_bound,
          &TaoVarStruct::set_key_bound,
          "Variable bound to keyboard key?"
      )
      .def_property(
          "v1",
          py::cpp_function(&TaoVarStruct::v1, py::keep_alive<0, 1>()),
          &TaoVarStruct::set_v1,
          "Pointer to the parent."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoVarStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoVarStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoVarStruct &self, const TaoVarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoVarStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoVarStructArray1D, TaoVarStructAlloc1D>(
      m,
      "TaoVarStructArray1D",
      "TaoVarStructAlloc1D"
  );
  // 2D TaoVarStruct arrays are not used in structs/routines
  // 3D TaoVarStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_kick_pt_struct
void init_tao_wave_kick_pt_struct(py::module &m, py::class_<TaoWaveKickPtStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             optional_ref<const EleStruct>>(),
         py::arg("phi_s") = py::none(),
         py::arg("phi_r") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("amp") = py::none(),
         py::arg("s") = py::none(),
         py::arg("ix_dat_before_kick") = py::none(),
         py::arg("ele") = py::none()
  )
      .def_property("phi_s", &TaoWaveKickPtStruct::phi_s, &TaoWaveKickPtStruct::set_phi_s)
      .def_property("phi_r", &TaoWaveKickPtStruct::phi_r, &TaoWaveKickPtStruct::set_phi_r)
      .def_property("phi", &TaoWaveKickPtStruct::phi, &TaoWaveKickPtStruct::set_phi)
      .def_property("amp", &TaoWaveKickPtStruct::amp, &TaoWaveKickPtStruct::set_amp)
      .def_property("s", &TaoWaveKickPtStruct::s, &TaoWaveKickPtStruct::set_s, "s-position of kick")
      .def_property(
          "ix_dat_before_kick",
          &TaoWaveKickPtStruct::ix_dat_before_kick,
          &TaoWaveKickPtStruct::set_ix_dat_before_kick,
          "Index of datum in data array just before the kick."
      )
      .def_property(
          "ele",
          py::cpp_function(&TaoWaveKickPtStruct::ele, py::keep_alive<0, 1>()),
          &TaoWaveKickPtStruct::set_ele,
          "lattice element at position of kick."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoWaveKickPtStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoWaveKickPtStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TaoWaveKickPtStruct &self, const TaoWaveKickPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoWaveKickPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoWaveKickPtStructArray1D, TaoWaveKickPtStructAlloc1D>(
      m,
      "TaoWaveKickPtStructArray1D",
      "TaoWaveKickPtStructAlloc1D"
  );
  // 2D TaoWaveKickPtStruct arrays are not used in structs/routines
  // 3D TaoWaveKickPtStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_wave_struct
void init_tao_wave_struct(py::module &m, py::class_<TaoWaveStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
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
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<int>>,
             std::optional<int>,
             optional_ref<const TaoGraphStruct>,
             optional_ref<const TaoPlotRegionStruct>,
             optional_ref<const TaoD1DataStruct>>(),
         py::arg("data_type") = py::none(),
         py::arg("rms_rel_a") = py::none(),
         py::arg("rms_rel_b") = py::none(),
         py::arg("rms_rel_as") = py::none(),
         py::arg("rms_rel_bs") = py::none(),
         py::arg("rms_rel_ar") = py::none(),
         py::arg("rms_rel_br") = py::none(),
         py::arg("rms_rel_k") = py::none(),
         py::arg("rms_rel_ks") = py::none(),
         py::arg("rms_rel_kr") = py::none(),
         py::arg("rms_phi") = py::none(),
         py::arg("rms_phi_s") = py::none(),
         py::arg("rms_phi_r") = py::none(),
         py::arg("amp_ba_s") = py::none(),
         py::arg("amp_ba_r") = py::none(),
         py::arg("chi_a") = py::none(),
         py::arg("chi_c") = py::none(),
         py::arg("chi_ba") = py::none(),
         py::arg("amp_a") = py::none(),
         py::arg("amp_b") = py::none(),
         py::arg("amp_ba") = py::none(),
         py::arg("coef_a") = py::none(),
         py::arg("coef_b") = py::none(),
         py::arg("coef_ba") = py::none(),
         py::arg("n_func") = py::none(),
         py::arg("ix_a1") = py::none(),
         py::arg("ix_a2") = py::none(),
         py::arg("ix_b1") = py::none(),
         py::arg("ix_b2") = py::none(),
         py::arg("i_a1") = py::none(),
         py::arg("i_a2") = py::none(),
         py::arg("i_b1") = py::none(),
         py::arg("i_b2") = py::none(),
         py::arg("n_a") = py::none(),
         py::arg("n_b") = py::none(),
         py::arg("i_curve_wrap_pt") = py::none(),
         py::arg("ix_data") = py::none(),
         py::arg("n_kick") = py::none(),
         py::arg("base_graph") = py::none(),
         py::arg("region") = py::none(),
         py::arg("d1_dat") = py::none()
  )
      .def_property("data_type", &TaoWaveStruct::data_type, &TaoWaveStruct::set_data_type)
      .def_property("rms_rel_a", &TaoWaveStruct::rms_rel_a, &TaoWaveStruct::set_rms_rel_a)
      .def_property("rms_rel_b", &TaoWaveStruct::rms_rel_b, &TaoWaveStruct::set_rms_rel_b)
      .def_property("rms_rel_as", &TaoWaveStruct::rms_rel_as, &TaoWaveStruct::set_rms_rel_as)
      .def_property("rms_rel_bs", &TaoWaveStruct::rms_rel_bs, &TaoWaveStruct::set_rms_rel_bs)
      .def_property("rms_rel_ar", &TaoWaveStruct::rms_rel_ar, &TaoWaveStruct::set_rms_rel_ar)
      .def_property("rms_rel_br", &TaoWaveStruct::rms_rel_br, &TaoWaveStruct::set_rms_rel_br)
      .def_property("rms_rel_k", &TaoWaveStruct::rms_rel_k, &TaoWaveStruct::set_rms_rel_k)
      .def_property("rms_rel_ks", &TaoWaveStruct::rms_rel_ks, &TaoWaveStruct::set_rms_rel_ks)
      .def_property("rms_rel_kr", &TaoWaveStruct::rms_rel_kr, &TaoWaveStruct::set_rms_rel_kr)
      .def_property("rms_phi", &TaoWaveStruct::rms_phi, &TaoWaveStruct::set_rms_phi)
      .def_property("rms_phi_s", &TaoWaveStruct::rms_phi_s, &TaoWaveStruct::set_rms_phi_s)
      .def_property("rms_phi_r", &TaoWaveStruct::rms_phi_r, &TaoWaveStruct::set_rms_phi_r)
      .def_property("amp_ba_s", &TaoWaveStruct::amp_ba_s, &TaoWaveStruct::set_amp_ba_s)
      .def_property("amp_ba_r", &TaoWaveStruct::amp_ba_r, &TaoWaveStruct::set_amp_ba_r)
      .def_property("chi_a", &TaoWaveStruct::chi_a, &TaoWaveStruct::set_chi_a)
      .def_property("chi_c", &TaoWaveStruct::chi_c, &TaoWaveStruct::set_chi_c)
      .def_property("chi_ba", &TaoWaveStruct::chi_ba, &TaoWaveStruct::set_chi_ba)
      .def_property(
          "amp_a",
          py::cpp_function(&TaoWaveStruct::amp_a, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_amp_a
      )
      .def_property(
          "amp_b",
          py::cpp_function(&TaoWaveStruct::amp_b, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_amp_b
      )
      .def_property(
          "amp_ba",
          py::cpp_function(&TaoWaveStruct::amp_ba, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_amp_ba
      )
      .def_property(
          "coef_a",
          py::cpp_function(&TaoWaveStruct::coef_a, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_coef_a
      )
      .def_property(
          "coef_b",
          py::cpp_function(&TaoWaveStruct::coef_b, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_coef_b
      )
      .def_property(
          "coef_ba",
          py::cpp_function(&TaoWaveStruct::coef_ba, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_coef_ba
      )
      .def_property(
          "n_func",
          &TaoWaveStruct::n_func,
          &TaoWaveStruct::set_n_func,
          "Number of functions used in the fit."
      )
      .def_property("ix_a1", &TaoWaveStruct::ix_a1, &TaoWaveStruct::set_ix_a1)
      .def_property("ix_a2", &TaoWaveStruct::ix_a2, &TaoWaveStruct::set_ix_a2)
      .def_property("ix_b1", &TaoWaveStruct::ix_b1, &TaoWaveStruct::set_ix_b1)
      .def_property("ix_b2", &TaoWaveStruct::ix_b2, &TaoWaveStruct::set_ix_b2)
      .def_property("i_a1", &TaoWaveStruct::i_a1, &TaoWaveStruct::set_i_a1)
      .def_property("i_a2", &TaoWaveStruct::i_a2, &TaoWaveStruct::set_i_a2)
      .def_property("i_b1", &TaoWaveStruct::i_b1, &TaoWaveStruct::set_i_b1)
      .def_property("i_b2", &TaoWaveStruct::i_b2, &TaoWaveStruct::set_i_b2)
      .def_property("n_a", &TaoWaveStruct::n_a, &TaoWaveStruct::set_n_a)
      .def_property("n_b", &TaoWaveStruct::n_b, &TaoWaveStruct::set_n_b)
      .def_property(
          "i_curve_wrap_pt",
          &TaoWaveStruct::i_curve_wrap_pt,
          &TaoWaveStruct::set_i_curve_wrap_pt,
          "Index of last point before wrap in curve array."
      )
      .def_property(
          "ix_data",
          py::cpp_function(&TaoWaveStruct::ix_data, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_ix_data,
          "Translates from plot point to datum index"
      )
      .def_property("n_kick", &TaoWaveStruct::n_kick, &TaoWaveStruct::set_n_kick)
      .def_property_readonly("kick", py::cpp_function(&TaoWaveStruct::kick, py::keep_alive<0, 1>()))
      .def_property(
          "base_graph",
          py::cpp_function(&TaoWaveStruct::base_graph, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_base_graph,
          "Graph before curves extended to 1.5 periods."
      )
      .def_property(
          "region",
          py::cpp_function(&TaoWaveStruct::region, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_region,
          "Where the wave plot is"
      )
      .def_property(
          "d1_dat",
          py::cpp_function(&TaoWaveStruct::d1_dat, py::keep_alive<0, 1>()),
          &TaoWaveStruct::set_d1_dat,
          "D1 data for analysis"
      )

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
      .def(
          "__eq__",
          [](const TaoWaveStruct &self, const TaoWaveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoWaveStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoWaveStruct arrays are not used in structs/routines
  // 2D TaoWaveStruct arrays are not used in structs/routines
  // 3D TaoWaveStruct arrays are not used in structs/routines
}

// =============================================================================
// test_sub_struct
void init_test_sub_struct(py::module &m, py::class_<TestSubStruct> &cls) {
  cls.def(py::init<optional_ref<const TestSubSubStruct>>(), py::arg("sr") = py::none())
      .def_property(
          "sr",
          py::cpp_function(&TestSubStruct::sr, py::keep_alive<0, 1>()),
          &TestSubStruct::set_sr
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TestSubStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TestSubStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const TestSubStruct &self, const TestSubStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TestSubStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TestSubStructArray1D, TestSubStructAlloc1D>(
      m,
      "TestSubStructArray1D",
      "TestSubStructAlloc1D"
  );
  bind_FTypeArrayND<TestSubStructArray2D>(m, "TestSubStructArray2D");
  bind_FTypeArrayND<TestSubStructArray3D>(m, "TestSubStructArray3D");
}

// =============================================================================
// test_sub_sub_struct
void init_test_sub_sub_struct(py::module &m, py::class_<TestSubSubStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int64_t>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("aaa") = py::none(),
         py::arg("bbb") = py::none(),
         py::arg("file") = py::none(),
         py::arg("t_ref") = py::none(),
         py::arg("freq_spread") = py::none()
  )
      .def_property("aaa", &TestSubSubStruct::aaa, &TestSubSubStruct::set_aaa)
      .def_property("bbb", &TestSubSubStruct::bbb, &TestSubSubStruct::set_bbb)
      .def_property("file", &TestSubSubStruct::file, &TestSubSubStruct::set_file)
      .def_property(
          "t_ref",
          &TestSubSubStruct::t_ref,
          &TestSubSubStruct::set_t_ref,
          "time reference value for computing the wake amplitude. This is used to prevent value "
          "overflow with long trains."
      )
      .def_property(
          "freq_spread",
          &TestSubSubStruct::freq_spread,
          &TestSubSubStruct::set_freq_spread,
          "Random frequency spread of long range modes."
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
      .def(
          "__eq__",
          [](const TestSubSubStruct &self, const TestSubSubStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const TestSubSubStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TestSubSubStruct arrays are not used in structs/routines
  // 2D TestSubSubStruct arrays are not used in structs/routines
  // 3D TestSubSubStruct arrays are not used in structs/routines
}