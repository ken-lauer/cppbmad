#include "pybmad/generated/structs_t.hpp"

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
// target_point_struct
void init_target_point_struct(nb::module_ &m, nb::class_<TargetPointStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<double>>>(), nb::arg("r") = nb::none())
      .def_prop_rw(
          "r",
          &TargetPointStruct::r,
          &TargetPointStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(x, y, z)"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TargetPointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TargetPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TargetPointStruct &self, nb::dict &memo) { return TargetPointStruct(self); }
      )
      .def(
          "__eq__",
          [](const TargetPointStruct &self, const TargetPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_taylor_struct(nb::module_ &m, nb::class_<TaylorStruct> &cls) {
  cls.def(nb::init<std::optional<double>>(), nb::arg("ref") = nb::none())
      .def_prop_rw("ref", &TaylorStruct::ref, &TaylorStruct::set_ref)
      .def_prop_ro("term", &TaylorStruct::term, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return TaylorStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaylorStruct &self, nb::dict &memo) { return TaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaylorStruct &self, const TaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_taylor_term_struct(nb::module_ &m, nb::class_<TaylorTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("expn") = nb::none()
  )
      .def_prop_rw("coef", &TaylorTermStruct::coef, &TaylorTermStruct::set_coef)
      .def_prop_rw(
          "expn",
          &TaylorTermStruct::expn,
          &TaylorTermStruct::set_expn,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaylorTermStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaylorTermStruct &self, nb::dict &memo) { return TaylorTermStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaylorTermStruct &self, const TaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_track_point_struct(nb::module_ &m, nb::class_<TrackPointStruct> &cls) {
  cls.def(
         "__init__",
         [](TrackPointStruct *self,
            std::optional<double> s_lab,
            std::optional<double> s_body,
            const CoordStruct *orb,
            const EmFieldStruct *field,
            const StrongBeamStruct *strong_beam,
            std::optional<std::vector<double>> vec0,
            std::optional<std::vector<std::vector<double>>> mat6) {
           new (self) TrackPointStruct(
               s_lab,
               s_body,
               ptr_to_opt_ref(orb),
               ptr_to_opt_ref(field),
               ptr_to_opt_ref(strong_beam),
               vec0,
               mat6
           );
         },
         nb::arg("s_lab") = nb::none(),
         nb::arg("s_body") = nb::none(),
         nb::arg("orb") = nb::none(),
         nb::arg("field") = nb::none(),
         nb::arg("strong_beam") = nb::none(),
         nb::arg("vec0") = nb::none(),
         nb::arg("mat6") = nb::none()
  )
      .def_prop_rw(
          "s_lab",
          &TrackPointStruct::s_lab,
          &TrackPointStruct::set_s_lab,
          "Longitudinal lab coord with respect to the upstream end."
      )
      .def_prop_rw(
          "s_body",
          &TrackPointStruct::s_body,
          &TrackPointStruct::set_s_body,
          "Longitudinal body coord with respect to the entrance end."
      )
      .def_prop_rw(
          "orb",
          &TrackPointStruct::orb,
          &TrackPointStruct::set_orb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Particle position in lab coords."
      )
      .def_prop_rw(
          "field",
          &TrackPointStruct::field,
          &TrackPointStruct::set_field,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "E&M fields in lab coordinates."
      )
      .def_prop_rw(
          "strong_beam",
          &TrackPointStruct::strong_beam,
          &TrackPointStruct::set_strong_beam,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Strong beam info for beambeam element."
      )
      .def_prop_rw(
          "vec0",
          &TrackPointStruct::vec0,
          &TrackPointStruct::set_vec0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th order part of xfer map from the beginning."
      )
      .def_prop_rw(
          "mat6",
          &TrackPointStruct::mat6,
          &TrackPointStruct::set_mat6,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1st order part of xfer map (transfer matrix)."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TrackPointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TrackPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TrackPointStruct &self, nb::dict &memo) { return TrackPointStruct(self); }
      )
      .def(
          "__eq__",
          [](const TrackPointStruct &self, const TrackPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_track_struct(nb::module_ &m, nb::class_<TrackStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("ds_save") = nb::none(),
         nb::arg("n_pt") = nb::none(),
         nb::arg("n_bad") = nb::none(),
         nb::arg("n_ok") = nb::none()
  )
      .def_prop_ro(
          "pt",
          &TrackStruct::pt,
          nb::keep_alive<0, 1>(),
          "Array of track points indexed from 0."
      )
      .def_prop_rw(
          "ds_save",
          &TrackStruct::ds_save,
          &TrackStruct::set_ds_save,
          "Min distance between points. Not positive => Save at all points."
      )
      .def_prop_rw(
          "n_pt",
          &TrackStruct::n_pt,
          &TrackStruct::set_n_pt,
          "Track upper bound for %pt(0:) array. n_bad and n_ok are used by adaptive trackers to "
          "record the number of times the step length had to be shortened."
      )
      .def_prop_rw(
          "n_bad",
          &TrackStruct::n_bad,
          &TrackStruct::set_n_bad,
          "Number of 'bad' steps where the step length was shortened."
      )
      .def_prop_rw(
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
          [](const TrackStruct &self, nb::dict &memo) { return TrackStruct(self); }
      )
      .def(
          "__eq__",
          [](const TrackStruct &self, const TrackStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_twiss_struct(nb::module_ &m, nb::class_<TwissStruct> &cls) {
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
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("beta") = nb::none(),
         nb::arg("alpha") = nb::none(),
         nb::arg("gamma") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("eta") = nb::none(),
         nb::arg("etap") = nb::none(),
         nb::arg("deta_ds") = nb::none(),
         nb::arg("sigma") = nb::none(),
         nb::arg("sigma_p") = nb::none(),
         nb::arg("emit") = nb::none(),
         nb::arg("norm_emit") = nb::none(),
         nb::arg("chrom") = nb::none(),
         nb::arg("dbeta_dpz") = nb::none(),
         nb::arg("dalpha_dpz") = nb::none(),
         nb::arg("deta_dpz") = nb::none(),
         nb::arg("detap_dpz") = nb::none()
  )
      .def_prop_rw("beta", &TwissStruct::beta, &TwissStruct::set_beta)
      .def_prop_rw("alpha", &TwissStruct::alpha, &TwissStruct::set_alpha)
      .def_prop_rw("gamma", &TwissStruct::gamma, &TwissStruct::set_gamma)
      .def_prop_rw("phi", &TwissStruct::phi, &TwissStruct::set_phi)
      .def_prop_rw("eta", &TwissStruct::eta, &TwissStruct::set_eta)
      .def_prop_rw("etap", &TwissStruct::etap, &TwissStruct::set_etap)
      .def_prop_rw("deta_ds", &TwissStruct::deta_ds, &TwissStruct::set_deta_ds)
      .def_prop_rw("sigma", &TwissStruct::sigma, &TwissStruct::set_sigma)
      .def_prop_rw("sigma_p", &TwissStruct::sigma_p, &TwissStruct::set_sigma_p)
      .def_prop_rw("emit", &TwissStruct::emit, &TwissStruct::set_emit)
      .def_prop_rw("norm_emit", &TwissStruct::norm_emit, &TwissStruct::set_norm_emit)
      .def_prop_rw("chrom", &TwissStruct::chrom, &TwissStruct::set_chrom)
      .def_prop_rw("dbeta_dpz", &TwissStruct::dbeta_dpz, &TwissStruct::set_dbeta_dpz)
      .def_prop_rw("dalpha_dpz", &TwissStruct::dalpha_dpz, &TwissStruct::set_dalpha_dpz)
      .def_prop_rw("deta_dpz", &TwissStruct::deta_dpz, &TwissStruct::set_deta_dpz)
      .def_prop_rw("detap_dpz", &TwissStruct::detap_dpz, &TwissStruct::set_detap_dpz)

      .def("__repr__", [](const TwissStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TwissStruct &self) {
            return TwissStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TwissStruct &self, nb::dict &memo) { return TwissStruct(self); }
      )
      .def(
          "__eq__",
          [](const TwissStruct &self, const TwissStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tricubic_cmplx_coef_struct(nb::module_ &m, nb::class_<TricubicCmplxCoefStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>,
             std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("i_box") = nb::none()
  )
      .def_prop_rw(
          "coef",
          &TricubicCmplxCoefStruct::coef,
          &TricubicCmplxCoefStruct::set_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Coefs"
      )
      .def_prop_rw(
          "i_box",
          &TricubicCmplxCoefStruct::i_box,
          &TricubicCmplxCoefStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const TricubicCmplxCoefStruct &self, nb::dict &memo) {
            return TricubicCmplxCoefStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TricubicCmplxCoefStruct &self, const TricubicCmplxCoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tricubic_coef_struct
void init_tricubic_coef_struct(nb::module_ &m, nb::class_<TricubicCoefStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("i_box") = nb::none()
  )
      .def_prop_rw(
          "coef",
          &TricubicCoefStruct::coef,
          &TricubicCoefStruct::set_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Coefs"
      )
      .def_prop_rw(
          "i_box",
          &TricubicCoefStruct::i_box,
          &TricubicCoefStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "index at lower box corner."
      )

      .def("__repr__", [](const TricubicCoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TricubicCoefStruct &self) {
            return TricubicCoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TricubicCoefStruct &self, nb::dict &memo) { return TricubicCoefStruct(self); }
      )
      .def(
          "__eq__",
          [](const TricubicCoefStruct &self, const TricubicCoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TricubicCoefStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TricubicCoefStruct arrays are not used in structs/routines
  // 2D TricubicCoefStruct arrays are not used in structs/routines
  // 3D TricubicCoefStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_ele_shape_input
void init_tao_ele_shape_input(nb::module_ &m, nb::class_<TaoEleShapeInput> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<double>>(),
         nb::arg("ele_id") = nb::none(),
         nb::arg("shape") = nb::none(),
         nb::arg("color") = nb::none(),
         nb::arg("size") = nb::none(),
         nb::arg("label") = nb::none(),
         nb::arg("draw") = nb::none(),
         nb::arg("multi") = nb::none(),
         nb::arg("line_width") = nb::none(),
         nb::arg("offset") = nb::none()
  )
      .def_prop_rw(
          "ele_id",
          &TaoEleShapeInput::ele_id,
          &TaoEleShapeInput::set_ele_id,
          "element 'key::name' to match to."
      )
      .def_prop_rw("shape", &TaoEleShapeInput::shape, &TaoEleShapeInput::set_shape, "Shape to draw")
      .def_prop_rw(
          "color",
          &TaoEleShapeInput::color,
          &TaoEleShapeInput::set_color,
          "Color of shape"
      )
      .def_prop_rw(
          "size",
          &TaoEleShapeInput::size,
          &TaoEleShapeInput::set_size,
          "plot vertical height"
      )
      .def_prop_rw(
          "label",
          &TaoEleShapeInput::label,
          &TaoEleShapeInput::set_label,
          "Can be: 'name', 's', 'none'"
      )
      .def_prop_rw("draw", &TaoEleShapeInput::draw, &TaoEleShapeInput::set_draw, "Draw the shape?")
      .def_prop_rw(
          "multi",
          &TaoEleShapeInput::multi,
          &TaoEleShapeInput::set_multi,
          "Can be part of a multi-shape."
      )
      .def_prop_rw(
          "line_width",
          &TaoEleShapeInput::line_width,
          &TaoEleShapeInput::set_line_width,
          "Width of lines used to draw the shape."
      )
      .def_prop_rw(
          "offset",
          &TaoEleShapeInput::offset,
          &TaoEleShapeInput::set_offset,
          "Vertical offset."
      )

      .def("__repr__", [](const TaoEleShapeInput &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoEleShapeInput &self) {
            return TaoEleShapeInput(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoEleShapeInput &self, nb::dict &memo) { return TaoEleShapeInput(self); }
      )
      .def(
          "__eq__",
          [](const TaoEleShapeInput &self, const TaoEleShapeInput &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoEleShapeInput &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoEleShapeInput arrays are not used in structs/routines
  // 2D TaoEleShapeInput arrays are not used in structs/routines
  // 3D TaoEleShapeInput arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_page_input
void init_tao_plot_page_input(nb::module_ &m, nb::class_<TaoPlotPageInput> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotPageInput *self,
            const TaoTitleStruct *title,
            const TaoTitleStruct *subtitle,
            const QpRectStruct *border,
            std::optional<std::string> plot_display_type,
            std::optional<std::vector<double>> size,
            std::optional<double> text_height,
            std::optional<double> main_title_text_scale,
            std::optional<double> graph_title_text_scale,
            std::optional<double> axis_number_text_scale,
            std::optional<double> axis_label_text_scale,
            std::optional<double> legend_text_scale,
            std::optional<double> key_table_text_scale,
            std::optional<double> floor_plan_shape_scale,
            std::optional<double> floor_plan_text_scale,
            std::optional<double> lat_layout_shape_scale,
            std::optional<double> lat_layout_text_scale,
            std::optional<double> curve_legend_line_len,
            std::optional<double> curve_legend_text_offset,
            std::optional<int> n_curve_pts,
            std::optional<bool> delete_overlapping_plots,
            std::optional<bool> draw_graph_title_suffix) {
           new (self) TaoPlotPageInput(
               ptr_to_opt_ref(title),
               ptr_to_opt_ref(subtitle),
               ptr_to_opt_ref(border),
               plot_display_type,
               size,
               text_height,
               main_title_text_scale,
               graph_title_text_scale,
               axis_number_text_scale,
               axis_label_text_scale,
               legend_text_scale,
               key_table_text_scale,
               floor_plan_shape_scale,
               floor_plan_text_scale,
               lat_layout_shape_scale,
               lat_layout_text_scale,
               curve_legend_line_len,
               curve_legend_text_offset,
               n_curve_pts,
               delete_overlapping_plots,
               draw_graph_title_suffix
           );
         },
         nb::arg("title") = nb::none(),
         nb::arg("subtitle") = nb::none(),
         nb::arg("border") = nb::none(),
         nb::arg("plot_display_type") = nb::none(),
         nb::arg("size") = nb::none(),
         nb::arg("text_height") = nb::none(),
         nb::arg("main_title_text_scale") = nb::none(),
         nb::arg("graph_title_text_scale") = nb::none(),
         nb::arg("axis_number_text_scale") = nb::none(),
         nb::arg("axis_label_text_scale") = nb::none(),
         nb::arg("legend_text_scale") = nb::none(),
         nb::arg("key_table_text_scale") = nb::none(),
         nb::arg("floor_plan_shape_scale") = nb::none(),
         nb::arg("floor_plan_text_scale") = nb::none(),
         nb::arg("lat_layout_shape_scale") = nb::none(),
         nb::arg("lat_layout_text_scale") = nb::none(),
         nb::arg("curve_legend_line_len") = nb::none(),
         nb::arg("curve_legend_text_offset") = nb::none(),
         nb::arg("n_curve_pts") = nb::none(),
         nb::arg("delete_overlapping_plots") = nb::none(),
         nb::arg("draw_graph_title_suffix") = nb::none()
  )
      .def_prop_rw(
          "title",
          &TaoPlotPageInput::title,
          &TaoPlotPageInput::set_title,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Title  at top of page."
      )
      .def_prop_rw(
          "subtitle",
          &TaoPlotPageInput::subtitle,
          &TaoPlotPageInput::set_subtitle,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Subtitle at top of page."
      )
      .def_prop_rw(
          "border",
          &TaoPlotPageInput::border,
          &TaoPlotPageInput::set_border,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Border around plots edge of page."
      )
      .def_prop_rw(
          "plot_display_type",
          &TaoPlotPageInput::plot_display_type,
          &TaoPlotPageInput::set_plot_display_type
      )
      .def_prop_rw(
          "size",
          &TaoPlotPageInput::size,
          &TaoPlotPageInput::set_size,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "width and height of window in pixels."
      )
      .def_prop_rw(
          "text_height",
          &TaoPlotPageInput::text_height,
          &TaoPlotPageInput::set_text_height,
          "In points. Scales the height of all text"
      )
      .def_prop_rw(
          "main_title_text_scale",
          &TaoPlotPageInput::main_title_text_scale,
          &TaoPlotPageInput::set_main_title_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "graph_title_text_scale",
          &TaoPlotPageInput::graph_title_text_scale,
          &TaoPlotPageInput::set_graph_title_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "axis_number_text_scale",
          &TaoPlotPageInput::axis_number_text_scale,
          &TaoPlotPageInput::set_axis_number_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "axis_label_text_scale",
          &TaoPlotPageInput::axis_label_text_scale,
          &TaoPlotPageInput::set_axis_label_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "legend_text_scale",
          &TaoPlotPageInput::legend_text_scale,
          &TaoPlotPageInput::set_legend_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "key_table_text_scale",
          &TaoPlotPageInput::key_table_text_scale,
          &TaoPlotPageInput::set_key_table_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "floor_plan_shape_scale",
          &TaoPlotPageInput::floor_plan_shape_scale,
          &TaoPlotPageInput::set_floor_plan_shape_scale
      )
      .def_prop_rw(
          "floor_plan_text_scale",
          &TaoPlotPageInput::floor_plan_text_scale,
          &TaoPlotPageInput::set_floor_plan_text_scale,
          "Scale used = floor_plan_text_scale * legend_text_scale"
      )
      .def_prop_rw(
          "lat_layout_shape_scale",
          &TaoPlotPageInput::lat_layout_shape_scale,
          &TaoPlotPageInput::set_lat_layout_shape_scale
      )
      .def_prop_rw(
          "lat_layout_text_scale",
          &TaoPlotPageInput::lat_layout_text_scale,
          &TaoPlotPageInput::set_lat_layout_text_scale,
          "Scale used = lat_layout_text_scale * legend_text_scale"
      )
      .def_prop_rw(
          "curve_legend_line_len",
          &TaoPlotPageInput::curve_legend_line_len,
          &TaoPlotPageInput::set_curve_legend_line_len,
          "OLD STYLE. Points."
      )
      .def_prop_rw(
          "curve_legend_text_offset",
          &TaoPlotPageInput::curve_legend_text_offset,
          &TaoPlotPageInput::set_curve_legend_text_offset,
          "OLD STYLE. Points."
      )
      .def_prop_rw(
          "n_curve_pts",
          &TaoPlotPageInput::n_curve_pts,
          &TaoPlotPageInput::set_n_curve_pts,
          "Number of points for plotting a smooth curve"
      )
      .def_prop_rw(
          "delete_overlapping_plots",
          &TaoPlotPageInput::delete_overlapping_plots,
          &TaoPlotPageInput::set_delete_overlapping_plots,
          "Delete overlapping plots when a plot is placed?"
      )
      .def_prop_rw(
          "draw_graph_title_suffix",
          &TaoPlotPageInput::draw_graph_title_suffix,
          &TaoPlotPageInput::set_draw_graph_title_suffix
      )

      .def("__repr__", [](const TaoPlotPageInput &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotPageInput &self) {
            return TaoPlotPageInput(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotPageInput &self, nb::dict &memo) { return TaoPlotPageInput(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotPageInput &self, const TaoPlotPageInput &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotPageInput &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D TaoPlotPageInput arrays are not used in structs/routines
  // 2D TaoPlotPageInput arrays are not used in structs/routines
  // 3D TaoPlotPageInput arrays are not used in structs/routines
}

// =============================================================================
// tao_alias_struct
void init_tao_alias_struct(nb::module_ &m, nb::class_<TaoAliasStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<std::string>>(),
         nb::arg("name") = nb::none(),
         nb::arg("expanded_str") = nb::none()
  )
      .def_prop_rw("name", &TaoAliasStruct::name, &TaoAliasStruct::set_name)
      .def_prop_rw("expanded_str", &TaoAliasStruct::expanded_str, &TaoAliasStruct::set_expanded_str)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoAliasStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoAliasStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoAliasStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoAliasStruct &self) {
            return TaoAliasStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoAliasStruct &self, nb::dict &memo) { return TaoAliasStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoAliasStruct &self, const TaoAliasStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoAliasStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoAliasStructArray1D, TaoAliasStructAlloc1D>(
      m,
      "TaoAliasStructArray1D",
      "TaoAliasStructAlloc1D"
  );
  // 2D TaoAliasStruct arrays are not used in structs/routines
  // 3D TaoAliasStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_beam_branch_struct
void init_tao_beam_branch_struct(nb::module_ &m, nb::class_<TaoBeamBranchStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoBeamBranchStruct *self,
            const BeamStruct *beam_at_start,
            const BeamInitStruct *beam_init,
            const BeamInitStruct *beam_init_used,
            std::optional<bool> init_starting_distribution,
            std::optional<std::string> track_start,
            std::optional<std::string> track_end,
            std::optional<int> ix_branch,
            std::optional<int> ix_track_start,
            std::optional<int> ix_track_end) {
           new (self) TaoBeamBranchStruct(
               ptr_to_opt_ref(beam_at_start),
               ptr_to_opt_ref(beam_init),
               ptr_to_opt_ref(beam_init_used),
               init_starting_distribution,
               track_start,
               track_end,
               ix_branch,
               ix_track_start,
               ix_track_end
           );
         },
         nb::arg("beam_at_start") = nb::none(),
         nb::arg("beam_init") = nb::none(),
         nb::arg("beam_init_used") = nb::none(),
         nb::arg("init_starting_distribution") = nb::none(),
         nb::arg("track_start") = nb::none(),
         nb::arg("track_end") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_track_start") = nb::none(),
         nb::arg("ix_track_end") = nb::none()
  )
      .def_prop_rw(
          "beam_at_start",
          &TaoBeamBranchStruct::beam_at_start,
          &TaoBeamBranchStruct::set_beam_at_start,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Initial beam"
      )
      .def_prop_rw(
          "beam_init",
          &TaoBeamBranchStruct::beam_init,
          &TaoBeamBranchStruct::set_beam_init,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "User set beam distrubution at track start."
      )
      .def_prop_rw(
          "beam_init_used",
          &TaoBeamBranchStruct::beam_init_used,
          &TaoBeamBranchStruct::set_beam_init_used,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "beam distribution with emit values set."
      )
      .def_prop_rw(
          "init_starting_distribution",
          &TaoBeamBranchStruct::init_starting_distribution,
          &TaoBeamBranchStruct::set_init_starting_distribution,
          "Init beam"
      )
      .def_prop_rw(
          "track_start",
          &TaoBeamBranchStruct::track_start,
          &TaoBeamBranchStruct::set_track_start,
          "Tracking start element."
      )
      .def_prop_rw(
          "track_end",
          &TaoBeamBranchStruct::track_end,
          &TaoBeamBranchStruct::set_track_end
      )
      .def_prop_rw(
          "ix_branch",
          &TaoBeamBranchStruct::ix_branch,
          &TaoBeamBranchStruct::set_ix_branch,
          "Branch tracked. If track_start or track_end is a lord, ix_track_start/end index will be "
          "a index of slave."
      )
      .def_prop_rw(
          "ix_track_start",
          &TaoBeamBranchStruct::ix_track_start,
          &TaoBeamBranchStruct::set_ix_track_start,
          "Element track start index."
      )
      .def_prop_rw(
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
          [](const TaoBeamBranchStruct &self, nb::dict &memo) { return TaoBeamBranchStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoBeamBranchStruct &self, const TaoBeamBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_beam_uni_struct(nb::module_ &m, nb::class_<TaoBeamUniStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("saved_at") = nb::none(),
         nb::arg("dump_file") = nb::none(),
         nb::arg("dump_at") = nb::none(),
         nb::arg("track_beam_in_universe") = nb::none(),
         nb::arg("always_reinit") = nb::none()
  )
      .def_prop_rw("saved_at", &TaoBeamUniStruct::saved_at, &TaoBeamUniStruct::set_saved_at)
      .def_prop_rw("dump_file", &TaoBeamUniStruct::dump_file, &TaoBeamUniStruct::set_dump_file)
      .def_prop_rw("dump_at", &TaoBeamUniStruct::dump_at, &TaoBeamUniStruct::set_dump_at)
      .def_prop_rw(
          "track_beam_in_universe",
          &TaoBeamUniStruct::track_beam_in_universe,
          &TaoBeamUniStruct::set_track_beam_in_universe,
          "Beam tracking enabled in this universe?"
      )
      .def_prop_rw(
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
          [](const TaoBeamUniStruct &self, nb::dict &memo) { return TaoBeamUniStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoBeamUniStruct &self, const TaoBeamUniStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
    nb::module_ &m,
    nb::class_<TaoBuildingWallOrientationStruct> &cls
) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("theta") = nb::none(),
         nb::arg("x_offset") = nb::none(),
         nb::arg("z_offset") = nb::none()
  )
      .def_prop_rw(
          "theta",
          &TaoBuildingWallOrientationStruct::theta,
          &TaoBuildingWallOrientationStruct::set_theta
      )
      .def_prop_rw(
          "x_offset",
          &TaoBuildingWallOrientationStruct::x_offset,
          &TaoBuildingWallOrientationStruct::set_x_offset
      )
      .def_prop_rw(
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
          [](const TaoBuildingWallOrientationStruct &self, nb::dict &memo) {
            return TaoBuildingWallOrientationStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoBuildingWallOrientationStruct &self,
             const TaoBuildingWallOrientationStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
    nb::module_ &m,
    nb::class_<TaoBuildingWallPointStruct> &cls
) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("z") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("radius") = nb::none(),
         nb::arg("z_center") = nb::none(),
         nb::arg("x_center") = nb::none()
  )
      .def_prop_rw(
          "z",
          &TaoBuildingWallPointStruct::z,
          &TaoBuildingWallPointStruct::set_z,
          "Global floor position"
      )
      .def_prop_rw(
          "x",
          &TaoBuildingWallPointStruct::x,
          &TaoBuildingWallPointStruct::set_x,
          "Global floor position"
      )
      .def_prop_rw(
          "radius",
          &TaoBuildingWallPointStruct::radius,
          &TaoBuildingWallPointStruct::set_radius,
          "Arc radius. +r -> CW rotation, same as bends."
      )
      .def_prop_rw(
          "z_center",
          &TaoBuildingWallPointStruct::z_center,
          &TaoBuildingWallPointStruct::set_z_center,
          "Arc center."
      )
      .def_prop_rw(
          "x_center",
          &TaoBuildingWallPointStruct::x_center,
          &TaoBuildingWallPointStruct::set_x_center,
          "Arc center."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoBuildingWallPointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoBuildingWallPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoBuildingWallPointStruct &self, nb::dict &memo) {
            return TaoBuildingWallPointStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoBuildingWallPointStruct &self, const TaoBuildingWallPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
    nb::module_ &m,
    nb::class_<TaoBuildingWallSectionStruct> &cls
) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<std::string>>(),
         nb::arg("name") = nb::none(),
         nb::arg("constraint") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoBuildingWallSectionStruct::name,
          &TaoBuildingWallSectionStruct::set_name
      )
      .def_prop_rw(
          "constraint",
          &TaoBuildingWallSectionStruct::constraint,
          &TaoBuildingWallSectionStruct::set_constraint,
          "'left_side' or 'right_side' constraint."
      )
      .def_prop_ro("point", &TaoBuildingWallSectionStruct::point, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return TaoBuildingWallSectionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoBuildingWallSectionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoBuildingWallSectionStruct &self, nb::dict &memo) {
            return TaoBuildingWallSectionStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoBuildingWallSectionStruct &self, const TaoBuildingWallSectionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_building_wall_struct(nb::module_ &m, nb::class_<TaoBuildingWallStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoBuildingWallStruct *self, const TaoBuildingWallOrientationStruct *orientation) {
           new (self) TaoBuildingWallStruct(ptr_to_opt_ref(orientation));
         },
         nb::arg("orientation") = nb::none()
  )
      .def_prop_rw(
          "orientation",
          &TaoBuildingWallStruct::orientation,
          &TaoBuildingWallStruct::set_orientation,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("section", &TaoBuildingWallStruct::section, nb::keep_alive<0, 1>())

      .def("__repr__", [](const TaoBuildingWallStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoBuildingWallStruct &self) {
            return TaoBuildingWallStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoBuildingWallStruct &self, nb::dict &memo) {
            return TaoBuildingWallStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoBuildingWallStruct &self, const TaoBuildingWallStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_cmd_history_struct(nb::module_ &m, nb::class_<TaoCmdHistoryStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("cmd") = nb::none(),
         nb::arg("ix") = nb::none()
  )
      .def_prop_rw("cmd", &TaoCmdHistoryStruct::cmd, &TaoCmdHistoryStruct::set_cmd, "The command")
      .def_prop_rw(
          "ix",
          &TaoCmdHistoryStruct::ix,
          &TaoCmdHistoryStruct::set_ix,
          "Command index (1st command has ix = 1, etc.) Note: Commands from command files will be "
          "assigned an index."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCmdHistoryStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCmdHistoryStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoCmdHistoryStruct &self, nb::dict &memo) { return TaoCmdHistoryStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCmdHistoryStruct &self, const TaoCmdHistoryStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_command_file_struct
void init_tao_command_file_struct(nb::module_ &m, nb::class_<TaoCommandFileStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("full_name") = nb::none(),
         nb::arg("dir") = nb::none(),
         nb::arg("ix_unit") = nb::none(),
         nb::arg("quiet") = nb::none(),
         nb::arg("paused") = nb::none(),
         nb::arg("n_line") = nb::none(),
         nb::arg("reset_at_end") = nb::none(),
         nb::arg("lattice_calc_save") = nb::none(),
         nb::arg("plot_save") = nb::none()
  )
      .def_prop_rw(
          "full_name",
          &TaoCommandFileStruct::full_name,
          &TaoCommandFileStruct::set_full_name
      )
      .def_prop_rw("dir", &TaoCommandFileStruct::dir, &TaoCommandFileStruct::set_dir)
      .def_prop_rw("ix_unit", &TaoCommandFileStruct::ix_unit, &TaoCommandFileStruct::set_ix_unit)
      .def_prop_ro(
          "cmd_arg",
          &TaoCommandFileStruct::cmd_arg,
          nb::keep_alive<0, 1>(),
          "Command file arguments."
      )
      .def_prop_rw("quiet", &TaoCommandFileStruct::quiet, &TaoCommandFileStruct::set_quiet)
      .def_prop_rw(
          "paused",
          &TaoCommandFileStruct::paused,
          &TaoCommandFileStruct::set_paused,
          "Is the command file paused?"
      )
      .def_prop_rw(
          "n_line",
          &TaoCommandFileStruct::n_line,
          &TaoCommandFileStruct::set_n_line,
          "Current line number"
      )
      .def_prop_rw(
          "reset_at_end",
          &TaoCommandFileStruct::reset_at_end,
          &TaoCommandFileStruct::set_reset_at_end,
          "Reset lattice_calc_on and plot_on at end of file?"
      )
      .def_prop_rw(
          "lattice_calc_save",
          &TaoCommandFileStruct::lattice_calc_save,
          &TaoCommandFileStruct::set_lattice_calc_save
      )
      .def_prop_rw(
          "plot_save",
          &TaoCommandFileStruct::plot_save,
          &TaoCommandFileStruct::set_plot_save
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCommandFileStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCommandFileStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoCommandFileStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCommandFileStruct &self) {
            return TaoCommandFileStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCommandFileStruct &self, nb::dict &memo) {
            return TaoCommandFileStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoCommandFileStruct &self, const TaoCommandFileStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCommandFileStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoCommandFileStructArray1D, TaoCommandFileStructAlloc1D>(
      m,
      "TaoCommandFileStructArray1D",
      "TaoCommandFileStructAlloc1D"
  );
  // 2D TaoCommandFileStruct arrays are not used in structs/routines
  // 3D TaoCommandFileStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_common_struct
void init_tao_common_struct(nb::module_ &m, nb::class_<TaoCommonStruct> &cls) {
  cls.def(
         nb::init<
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
             std::optional<std::string>,
             std::optional<std::string>>(),
         nb::arg("covar") = nb::none(),
         nb::arg("alpha") = nb::none(),
         nb::arg("dummy_target") = nb::none(),
         nb::arg("n_alias") = nb::none(),
         nb::arg("cmd_file_level") = nb::none(),
         nb::arg("ix_key_bank") = nb::none(),
         nb::arg("ix_history") = nb::none(),
         nb::arg("n_history") = nb::none(),
         nb::arg("lev_loop") = nb::none(),
         nb::arg("n_err_messages_printed") = nb::none(),
         nb::arg("n_universes") = nb::none(),
         nb::arg("ix_beam_track_active_element") = nb::none(),
         nb::arg("cmd_file_paused") = nb::none(),
         nb::arg("use_cmd_here") = nb::none(),
         nb::arg("cmd_from_cmd_file") = nb::none(),
         nb::arg("use_saved_beam_in_tracking") = nb::none(),
         nb::arg("single_mode") = nb::none(),
         nb::arg("combine_consecutive_elements_of_like_name") = nb::none(),
         nb::arg("have_tracked_beam") = nb::none(),
         nb::arg("init_plot_needed") = nb::none(),
         nb::arg("init_beam") = nb::none(),
         nb::arg("init_var") = nb::none(),
         nb::arg("init_read_lat_info") = nb::none(),
         nb::arg("optimizer_running") = nb::none(),
         nb::arg("have_datums_using_expressions") = nb::none(),
         nb::arg("print_to_terminal") = nb::none(),
         nb::arg("lattice_calc_done") = nb::none(),
         nb::arg("add_measurement_noise") = nb::none(),
         nb::arg("is_err_message_printed") = nb::none(),
         nb::arg("command_arg_has_been_executed") = nb::none(),
         nb::arg("all_merit_weights_positive") = nb::none(),
         nb::arg("multi_turn_orbit_is_plotted") = nb::none(),
         nb::arg("force_chrom_calc") = nb::none(),
         nb::arg("force_rad_int_calc") = nb::none(),
         nb::arg("rad_int_ri_calc_on") = nb::none(),
         nb::arg("rad_int_6d_calc_on") = nb::none(),
         nb::arg("single_mode_buffer") = nb::none(),
         nb::arg("cmd") = nb::none(),
         nb::arg("saved_cmd_line") = nb::none()
  )
      .def_prop_ro("alias", &TaoCommonStruct::alias, nb::keep_alive<0, 1>())
      .def_prop_ro("key", &TaoCommonStruct::key, nb::keep_alive<0, 1>())
      .def_prop_ro("cmd_file", &TaoCommonStruct::cmd_file, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "symbolic_num",
          &TaoCommonStruct::symbolic_num,
          nb::keep_alive<0, 1>(),
          "Named numbers"
      )
      .def_prop_ro(
          "plot_place_buffer",
          &TaoCommonStruct::plot_place_buffer,
          nb::keep_alive<0, 1>(),
          "Used when %external_plotting is on."
      )
      .def_prop_ro("do_loop", &TaoCommonStruct::do_loop, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "covar",
          &TaoCommonStruct::covar,
          &TaoCommonStruct::set_covar,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "alpha",
          &TaoCommonStruct::alpha,
          &TaoCommonStruct::set_alpha,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dummy_target",
          &TaoCommonStruct::dummy_target,
          &TaoCommonStruct::set_dummy_target,
          "Dummy varaible"
      )
      .def_prop_rw("n_alias", &TaoCommonStruct::n_alias, &TaoCommonStruct::set_n_alias)
      .def_prop_rw(
          "cmd_file_level",
          &TaoCommonStruct::cmd_file_level,
          &TaoCommonStruct::set_cmd_file_level,
          "For nested command files. 0 -> no command file."
      )
      .def_prop_rw(
          "ix_key_bank",
          &TaoCommonStruct::ix_key_bank,
          &TaoCommonStruct::set_ix_key_bank,
          "For single mode."
      )
      .def_prop_rw(
          "ix_history",
          &TaoCommonStruct::ix_history,
          &TaoCommonStruct::set_ix_history,
          "Index to latest command in the history circular buffer."
      )
      .def_prop_rw(
          "n_history",
          &TaoCommonStruct::n_history,
          &TaoCommonStruct::set_n_history,
          "Number of commands issued from beginning of starting Tao."
      )
      .def_prop_rw(
          "lev_loop",
          &TaoCommonStruct::lev_loop,
          &TaoCommonStruct::set_lev_loop,
          "in do loop nest level"
      )
      .def_prop_rw(
          "n_err_messages_printed",
          &TaoCommonStruct::n_err_messages_printed,
          &TaoCommonStruct::set_n_err_messages_printed,
          "Used by tao_set_invalid to limit number of messages."
      )
      .def_prop_rw("n_universes", &TaoCommonStruct::n_universes, &TaoCommonStruct::set_n_universes)
      .def_prop_rw(
          "ix_beam_track_active_element",
          &TaoCommonStruct::ix_beam_track_active_element,
          &TaoCommonStruct::set_ix_beam_track_active_element,
          "Element being tracked through `tao_beam_track`."
      )
      .def_prop_rw(
          "cmd_file_paused",
          &TaoCommonStruct::cmd_file_paused,
          &TaoCommonStruct::set_cmd_file_paused
      )
      .def_prop_rw(
          "use_cmd_here",
          &TaoCommonStruct::use_cmd_here,
          &TaoCommonStruct::set_use_cmd_here,
          "Used for commands recalled from the cmd history stack"
      )
      .def_prop_rw(
          "cmd_from_cmd_file",
          &TaoCommonStruct::cmd_from_cmd_file,
          &TaoCommonStruct::set_cmd_from_cmd_file,
          "was command from a command file?"
      )
      .def_prop_rw(
          "use_saved_beam_in_tracking",
          &TaoCommonStruct::use_saved_beam_in_tracking,
          &TaoCommonStruct::set_use_saved_beam_in_tracking
      )
      .def_prop_rw("single_mode", &TaoCommonStruct::single_mode, &TaoCommonStruct::set_single_mode)
      .def_prop_rw(
          "combine_consecutive_elements_of_like_name",
          &TaoCommonStruct::combine_consecutive_elements_of_like_name,
          &TaoCommonStruct::set_combine_consecutive_elements_of_like_name
      )
      .def_prop_rw(
          "have_tracked_beam",
          &TaoCommonStruct::have_tracked_beam,
          &TaoCommonStruct::set_have_tracked_beam,
          "Used to catch error when beam plotting without having tracked a beam."
      )
      .def_prop_rw(
          "init_plot_needed",
          &TaoCommonStruct::init_plot_needed,
          &TaoCommonStruct::set_init_plot_needed,
          "reinitialize plotting?"
      )
      .def_prop_rw(
          "init_beam",
          &TaoCommonStruct::init_beam,
          &TaoCommonStruct::set_init_beam,
          "Used by custom programs to control Tao init"
      )
      .def_prop_rw(
          "init_var",
          &TaoCommonStruct::init_var,
          &TaoCommonStruct::set_init_var,
          "Used by custom programs to control Tao init"
      )
      .def_prop_rw(
          "init_read_lat_info",
          &TaoCommonStruct::init_read_lat_info,
          &TaoCommonStruct::set_init_read_lat_info,
          "Used by custom programs to control Tao init"
      )
      .def_prop_rw(
          "optimizer_running",
          &TaoCommonStruct::optimizer_running,
          &TaoCommonStruct::set_optimizer_running
      )
      .def_prop_rw(
          "have_datums_using_expressions",
          &TaoCommonStruct::have_datums_using_expressions,
          &TaoCommonStruct::set_have_datums_using_expressions
      )
      .def_prop_rw(
          "print_to_terminal",
          &TaoCommonStruct::print_to_terminal,
          &TaoCommonStruct::set_print_to_terminal,
          "Print command prompt to the terminal? For use with GUIs."
      )
      .def_prop_rw(
          "lattice_calc_done",
          &TaoCommonStruct::lattice_calc_done,
          &TaoCommonStruct::set_lattice_calc_done,
          "Used by GUI for deciding when to refresh."
      )
      .def_prop_rw(
          "add_measurement_noise",
          &TaoCommonStruct::add_measurement_noise,
          &TaoCommonStruct::set_add_measurement_noise,
          "Turn off to take data derivatives."
      )
      .def_prop_rw(
          "is_err_message_printed",
          &TaoCommonStruct::is_err_message_printed,
          &TaoCommonStruct::set_is_err_message_printed,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Used by tao_set_invalid"
      )
      .def_prop_rw(
          "command_arg_has_been_executed",
          &TaoCommonStruct::command_arg_has_been_executed,
          &TaoCommonStruct::set_command_arg_has_been_executed,
          "Has the -command command line argument been executed?"
      )
      .def_prop_rw(
          "all_merit_weights_positive",
          &TaoCommonStruct::all_merit_weights_positive,
          &TaoCommonStruct::set_all_merit_weights_positive
      )
      .def_prop_rw(
          "multi_turn_orbit_is_plotted",
          &TaoCommonStruct::multi_turn_orbit_is_plotted,
          &TaoCommonStruct::set_multi_turn_orbit_is_plotted,
          "Is a multi_turn_orbit being plotted?"
      )
      .def_prop_rw(
          "force_chrom_calc",
          &TaoCommonStruct::force_chrom_calc,
          &TaoCommonStruct::set_force_chrom_calc,
          "Used by a routine to force a single chromaticity calculation."
      )
      .def_prop_rw(
          "force_rad_int_calc",
          &TaoCommonStruct::force_rad_int_calc,
          &TaoCommonStruct::set_force_rad_int_calc,
          "Used by a routine to force a single radiation integrals calculation"
      )
      .def_prop_rw(
          "rad_int_ri_calc_on",
          &TaoCommonStruct::rad_int_ri_calc_on,
          &TaoCommonStruct::set_rad_int_ri_calc_on,
          "'Classical' radiation integrals calculation on/off."
      )
      .def_prop_rw(
          "rad_int_6d_calc_on",
          &TaoCommonStruct::rad_int_6d_calc_on,
          &TaoCommonStruct::set_rad_int_6d_calc_on,
          "6D Radiation integrals calculation on/off."
      )
      .def_prop_ro(
          "valid_plot_who",
          &TaoCommonStruct::valid_plot_who,
          nb::keep_alive<0, 1>(),
          "model, base, ref etc..."
      )
      .def_prop_rw(
          "single_mode_buffer",
          &TaoCommonStruct::single_mode_buffer,
          &TaoCommonStruct::set_single_mode_buffer
      )
      .def_prop_rw(
          "cmd",
          &TaoCommonStruct::cmd,
          &TaoCommonStruct::set_cmd,
          "Used for the cmd history"
      )
      .def_prop_rw(
          "saved_cmd_line",
          &TaoCommonStruct::saved_cmd_line,
          &TaoCommonStruct::set_saved_cmd_line,
          "Saved part of command line when there are mulitple commands on a line"
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
          [](const TaoCommonStruct &self, nb::dict &memo) { return TaoCommonStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCommonStruct &self, const TaoCommonStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_curve_array_struct
void init_tao_curve_array_struct(nb::module_ &m, nb::class_<TaoCurveArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoCurveArrayStruct *self, const TaoCurveStruct *c) {
           new (self) TaoCurveArrayStruct(ptr_to_opt_ref(c));
         },
         nb::arg("c") = nb::none()
  )
      .def_prop_rw(
          "c",
          &TaoCurveArrayStruct::c,
          &TaoCurveArrayStruct::set_c,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCurveArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCurveArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoCurveArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveArrayStruct &self) {
            return TaoCurveArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCurveArrayStruct &self, nb::dict &memo) { return TaoCurveArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCurveArrayStruct &self, const TaoCurveArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoCurveArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoCurveArrayStructArray1D, TaoCurveArrayStructAlloc1D>(
      m,
      "TaoCurveArrayStructArray1D",
      "TaoCurveArrayStructAlloc1D"
  );
  // 2D TaoCurveArrayStruct arrays are not used in structs/routines
  // 3D TaoCurveArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_curve_color_struct
void init_tao_curve_color_struct(nb::module_ &m, nb::class_<TaoCurveColorStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>>(),
         nb::arg("data_type") = nb::none(),
         nb::arg("is_on") = nb::none(),
         nb::arg("min") = nb::none(),
         nb::arg("max") = nb::none(),
         nb::arg("autoscale") = nb::none()
  )
      .def_prop_rw(
          "data_type",
          &TaoCurveColorStruct::data_type,
          &TaoCurveColorStruct::set_data_type,
          "Datum type to use for z-axis."
      )
      .def_prop_rw("is_on", &TaoCurveColorStruct::is_on, &TaoCurveColorStruct::set_is_on, "On/Off")
      .def_prop_rw(
          "min",
          &TaoCurveColorStruct::min,
          &TaoCurveColorStruct::set_min,
          "Min and max values for mapping z-axis to color."
      )
      .def_prop_rw(
          "max",
          &TaoCurveColorStruct::max,
          &TaoCurveColorStruct::set_max,
          "Min and max values for mapping z-axis to color."
      )
      .def_prop_rw(
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
          [](const TaoCurveColorStruct &self, nb::dict &memo) { return TaoCurveColorStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCurveColorStruct &self, const TaoCurveColorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_curve_orbit_struct(nb::module_ &m, nb::class_<TaoCurveOrbitStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("t") = nb::none()
  )
      .def_prop_rw("x", &TaoCurveOrbitStruct::x, &TaoCurveOrbitStruct::set_x, "Transverse offset")
      .def_prop_rw("y", &TaoCurveOrbitStruct::y, &TaoCurveOrbitStruct::set_y, "Transverse offset")
      .def_prop_rw("t", &TaoCurveOrbitStruct::t, &TaoCurveOrbitStruct::set_t, "Time")

      .def("__repr__", [](const TaoCurveOrbitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoCurveOrbitStruct &self) {
            return TaoCurveOrbitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoCurveOrbitStruct &self, nb::dict &memo) { return TaoCurveOrbitStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCurveOrbitStruct &self, const TaoCurveOrbitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_curve_struct(nb::module_ &m, nb::class_<TaoCurveStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoCurveStruct *self,
            std::optional<std::string> name,
            std::optional<std::string> data_source,
            std::optional<std::string> data_index,
            std::optional<std::string> data_type_x,
            std::optional<std::string> data_type,
            std::optional<std::string> ele_ref_name,
            std::optional<std::string> legend_text,
            std::optional<std::string> message_text,
            std::optional<std::string> component,
            std::optional<std::string> why_invalid,
            const TaoGraphStruct *g,
            const TaoHistogramStruct *hist,
            const TaoCurveColorStruct *z_color,
            std::optional<std::vector<double>> x_line,
            std::optional<std::vector<double>> y_line,
            std::optional<std::vector<double>> y2_line,
            std::optional<std::vector<int>> ix_line,
            std::optional<std::vector<double>> x_symb,
            std::optional<std::vector<double>> y_symb,
            std::optional<std::vector<double>> z_symb,
            std::optional<std::vector<double>> err_symb,
            std::optional<std::vector<double>> symb_size,
            std::optional<std::vector<int>> ix_symb,
            std::optional<double> y_axis_scale_factor,
            const QpLineStruct *line,
            const QpSymbolStruct *symbol,
            const TaoCurveOrbitStruct *orbit,
            std::optional<int> ix_universe,
            std::optional<int> symbol_every,
            std::optional<int> ix_branch,
            std::optional<int> ix_bunch,
            std::optional<int> n_turn,
            std::optional<bool> use_y2,
            std::optional<bool> draw_line,
            std::optional<bool> draw_symbols,
            std::optional<bool> draw_symbol_index,
            std::optional<bool> draw_error_bars,
            std::optional<bool> smooth_line_calc,
            std::optional<bool> valid) {
           new (self) TaoCurveStruct(
               name,
               data_source,
               data_index,
               data_type_x,
               data_type,
               ele_ref_name,
               legend_text,
               message_text,
               component,
               why_invalid,
               ptr_to_opt_ref(g),
               ptr_to_opt_ref(hist),
               ptr_to_opt_ref(z_color),
               x_line,
               y_line,
               y2_line,
               ix_line,
               x_symb,
               y_symb,
               z_symb,
               err_symb,
               symb_size,
               ix_symb,
               y_axis_scale_factor,
               ptr_to_opt_ref(line),
               ptr_to_opt_ref(symbol),
               ptr_to_opt_ref(orbit),
               ix_universe,
               symbol_every,
               ix_branch,
               ix_bunch,
               n_turn,
               use_y2,
               draw_line,
               draw_symbols,
               draw_symbol_index,
               draw_error_bars,
               smooth_line_calc,
               valid
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("data_source") = nb::none(),
         nb::arg("data_index") = nb::none(),
         nb::arg("data_type_x") = nb::none(),
         nb::arg("data_type") = nb::none(),
         nb::arg("ele_ref_name") = nb::none(),
         nb::arg("legend_text") = nb::none(),
         nb::arg("message_text") = nb::none(),
         nb::arg("component") = nb::none(),
         nb::arg("why_invalid") = nb::none(),
         nb::arg("g") = nb::none(),
         nb::arg("hist") = nb::none(),
         nb::arg("z_color") = nb::none(),
         nb::arg("x_line") = nb::none(),
         nb::arg("y_line") = nb::none(),
         nb::arg("y2_line") = nb::none(),
         nb::arg("ix_line") = nb::none(),
         nb::arg("x_symb") = nb::none(),
         nb::arg("y_symb") = nb::none(),
         nb::arg("z_symb") = nb::none(),
         nb::arg("err_symb") = nb::none(),
         nb::arg("symb_size") = nb::none(),
         nb::arg("ix_symb") = nb::none(),
         nb::arg("y_axis_scale_factor") = nb::none(),
         nb::arg("line") = nb::none(),
         nb::arg("symbol") = nb::none(),
         nb::arg("orbit") = nb::none(),
         nb::arg("ix_universe") = nb::none(),
         nb::arg("symbol_every") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_bunch") = nb::none(),
         nb::arg("n_turn") = nb::none(),
         nb::arg("use_y2") = nb::none(),
         nb::arg("draw_line") = nb::none(),
         nb::arg("draw_symbols") = nb::none(),
         nb::arg("draw_symbol_index") = nb::none(),
         nb::arg("draw_error_bars") = nb::none(),
         nb::arg("smooth_line_calc") = nb::none(),
         nb::arg("valid") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoCurveStruct::name,
          &TaoCurveStruct::set_name,
          "Name identifying the curve."
      )
      .def_prop_rw(
          "data_source",
          &TaoCurveStruct::data_source,
          &TaoCurveStruct::set_data_source,
          "'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit'"
      )
      .def_prop_rw(
          "data_index",
          &TaoCurveStruct::data_index,
          &TaoCurveStruct::set_data_index,
          "Used for calculating %ix_symb(:)."
      )
      .def_prop_rw(
          "data_type_x",
          &TaoCurveStruct::data_type_x,
          &TaoCurveStruct::set_data_type_x,
          "Used for data slices and phase space plots."
      )
      .def_prop_rw(
          "data_type",
          &TaoCurveStruct::data_type,
          &TaoCurveStruct::set_data_type,
          "'orbit.x', etc."
      )
      .def_prop_rw(
          "ele_ref_name",
          &TaoCurveStruct::ele_ref_name,
          &TaoCurveStruct::set_ele_ref_name,
          "Reference element."
      )
      .def_prop_rw(
          "legend_text",
          &TaoCurveStruct::legend_text,
          &TaoCurveStruct::set_legend_text,
          "String to draw in a curve legend."
      )
      .def_prop_rw(
          "message_text",
          &TaoCurveStruct::message_text,
          &TaoCurveStruct::set_message_text,
          "Informational message to draw with graph."
      )
      .def_prop_rw(
          "component",
          &TaoCurveStruct::component,
          &TaoCurveStruct::set_component,
          "Who to plot. Eg: 'meas - design'"
      )
      .def_prop_rw(
          "why_invalid",
          &TaoCurveStruct::why_invalid,
          &TaoCurveStruct::set_why_invalid,
          "Informative string to print."
      )
      .def_prop_rw(
          "g",
          &TaoCurveStruct::g,
          &TaoCurveStruct::set_g,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "pointer to parent graph"
      )
      .def_prop_rw(
          "hist",
          &TaoCurveStruct::hist,
          &TaoCurveStruct::set_hist,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "z_color",
          &TaoCurveStruct::z_color,
          &TaoCurveStruct::set_z_color,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "x_line",
          &TaoCurveStruct::x_line,
          &TaoCurveStruct::set_x_line,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Coords for drawing a curve"
      )
      .def_prop_rw(
          "y_line",
          &TaoCurveStruct::y_line,
          &TaoCurveStruct::set_y_line,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "y2_line",
          &TaoCurveStruct::y2_line,
          &TaoCurveStruct::set_y2_line,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Second array needed for beam chamber curve."
      )
      .def_prop_rw(
          "ix_line",
          &TaoCurveStruct::ix_line,
          &TaoCurveStruct::set_ix_line,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Used by wave and aperture curves."
      )
      .def_prop_rw(
          "x_symb",
          &TaoCurveStruct::x_symb,
          &TaoCurveStruct::set_x_symb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Coords for drawing the symbols"
      )
      .def_prop_rw(
          "y_symb",
          &TaoCurveStruct::y_symb,
          &TaoCurveStruct::set_y_symb,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "z_symb",
          &TaoCurveStruct::z_symb,
          &TaoCurveStruct::set_z_symb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Symbol color"
      )
      .def_prop_rw(
          "err_symb",
          &TaoCurveStruct::err_symb,
          &TaoCurveStruct::set_err_symb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Error bars"
      )
      .def_prop_rw(
          "symb_size",
          &TaoCurveStruct::symb_size,
          &TaoCurveStruct::set_symb_size,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Symbol size. Used with symbol_size_scale."
      )
      .def_prop_rw(
          "ix_symb",
          &TaoCurveStruct::ix_symb,
          &TaoCurveStruct::set_ix_symb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Corresponding index in d1_data%d(:) array."
      )
      .def_prop_rw(
          "y_axis_scale_factor",
          &TaoCurveStruct::y_axis_scale_factor,
          &TaoCurveStruct::set_y_axis_scale_factor,
          "y-axis conversion from internal to plotting units."
      )
      .def_prop_rw(
          "line",
          &TaoCurveStruct::line,
          &TaoCurveStruct::set_line,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Line attributes"
      )
      .def_prop_rw(
          "symbol",
          &TaoCurveStruct::symbol,
          &TaoCurveStruct::set_symbol,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Symbol attributes"
      )
      .def_prop_rw(
          "orbit",
          &TaoCurveStruct::orbit,
          &TaoCurveStruct::set_orbit,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Used for E/B field plotting."
      )
      .def_prop_rw(
          "ix_universe",
          &TaoCurveStruct::ix_universe,
          &TaoCurveStruct::set_ix_universe,
          "Universe where data is. -1 => use s%global%default_universe"
      )
      .def_prop_rw(
          "symbol_every",
          &TaoCurveStruct::symbol_every,
          &TaoCurveStruct::set_symbol_every,
          "Symbol every how many points."
      )
      .def_prop_rw("ix_branch", &TaoCurveStruct::ix_branch, &TaoCurveStruct::set_ix_branch)
      .def_prop_rw(
          "ix_bunch",
          &TaoCurveStruct::ix_bunch,
          &TaoCurveStruct::set_ix_bunch,
          "Bunch to plot."
      )
      .def_prop_rw(
          "n_turn",
          &TaoCurveStruct::n_turn,
          &TaoCurveStruct::set_n_turn,
          "Used for multi_turn_orbit plotting"
      )
      .def_prop_rw("use_y2", &TaoCurveStruct::use_y2, &TaoCurveStruct::set_use_y2, "Use y2 axis?")
      .def_prop_rw(
          "draw_line",
          &TaoCurveStruct::draw_line,
          &TaoCurveStruct::set_draw_line,
          "Draw a line through the data points?"
      )
      .def_prop_rw(
          "draw_symbols",
          &TaoCurveStruct::draw_symbols,
          &TaoCurveStruct::set_draw_symbols,
          "Draw a symbol at the data points?"
      )
      .def_prop_rw(
          "draw_symbol_index",
          &TaoCurveStruct::draw_symbol_index,
          &TaoCurveStruct::set_draw_symbol_index,
          "Draw the symbol index number curve%ix_symb?"
      )
      .def_prop_rw(
          "draw_error_bars",
          &TaoCurveStruct::draw_error_bars,
          &TaoCurveStruct::set_draw_error_bars,
          "Draw error bars based upon data%error_rms if drawing data? !! logical :: draw_rms = "
          ".false.          ! Show mean and RMS values with legend?"
      )
      .def_prop_rw(
          "smooth_line_calc",
          &TaoCurveStruct::smooth_line_calc,
          &TaoCurveStruct::set_smooth_line_calc,
          "Calculate data between element edge points?"
      )
      .def_prop_rw("valid", &TaoCurveStruct::valid, &TaoCurveStruct::set_valid, "valid data?")
      .def_static(
          "new_array1d",
          [](int sz) { return TaoCurveStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoCurveStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoCurveStruct &self, nb::dict &memo) { return TaoCurveStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoCurveStruct &self, const TaoCurveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_d1_data_array_struct
void init_tao_d1_data_array_struct(nb::module_ &m, nb::class_<TaoD1DataArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoD1DataArrayStruct *self, const TaoD1DataStruct *d1) {
           new (self) TaoD1DataArrayStruct(ptr_to_opt_ref(d1));
         },
         nb::arg("d1") = nb::none()
  )
      .def_prop_rw(
          "d1",
          &TaoD1DataArrayStruct::d1,
          &TaoD1DataArrayStruct::set_d1,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD1DataArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD1DataArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoD1DataArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD1DataArrayStruct &self) {
            return TaoD1DataArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoD1DataArrayStruct &self, nb::dict &memo) {
            return TaoD1DataArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoD1DataArrayStruct &self, const TaoD1DataArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoD1DataArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoD1DataArrayStructArray1D, TaoD1DataArrayStructAlloc1D>(
      m,
      "TaoD1DataArrayStructArray1D",
      "TaoD1DataArrayStructAlloc1D"
  );
  // 2D TaoD1DataArrayStruct arrays are not used in structs/routines
  // 3D TaoD1DataArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d1_data_struct
void init_tao_d1_data_struct(nb::module_ &m, nb::class_<TaoD1DataStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoD1DataStruct *self, std::optional<std::string> name, const TaoD2DataStruct *d2) {
           new (self) TaoD1DataStruct(name, ptr_to_opt_ref(d2));
         },
         nb::arg("name") = nb::none(),
         nb::arg("d2") = nb::none()
  )
      .def_prop_rw("name", &TaoD1DataStruct::name, &TaoD1DataStruct::set_name, "Eg: 'x', etc.")
      .def_prop_rw(
          "d2",
          &TaoD1DataStruct::d2,
          &TaoD1DataStruct::set_d2,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "ptr to parent d2_data"
      )
      .def_prop_ro(
          "d",
          &TaoD1DataStruct::d,
          nb::keep_alive<0, 1>(),
          "Pointer to the appropriate section in u%data"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD1DataStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD1DataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoD1DataStruct &self, nb::dict &memo) { return TaoD1DataStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoD1DataStruct &self, const TaoD1DataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_d2_data_array_struct
void init_tao_d2_data_array_struct(nb::module_ &m, nb::class_<TaoD2DataArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoD2DataArrayStruct *self, const TaoD2DataStruct *d2) {
           new (self) TaoD2DataArrayStruct(ptr_to_opt_ref(d2));
         },
         nb::arg("d2") = nb::none()
  )
      .def_prop_rw(
          "d2",
          &TaoD2DataArrayStruct::d2,
          &TaoD2DataArrayStruct::set_d2,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD2DataArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD2DataArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoD2DataArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoD2DataArrayStruct &self) {
            return TaoD2DataArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoD2DataArrayStruct &self, nb::dict &memo) {
            return TaoD2DataArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoD2DataArrayStruct &self, const TaoD2DataArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoD2DataArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoD2DataArrayStructArray1D, TaoD2DataArrayStructAlloc1D>(
      m,
      "TaoD2DataArrayStructArray1D",
      "TaoD2DataArrayStructAlloc1D"
  );
  // 2D TaoD2DataArrayStruct arrays are not used in structs/routines
  // 3D TaoD2DataArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_d2_data_struct
void init_tao_d2_data_struct(nb::module_ &m, nb::class_<TaoD2DataStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("name") = nb::none(),
         nb::arg("data_file_name") = nb::none(),
         nb::arg("ref_file_name") = nb::none(),
         nb::arg("data_date") = nb::none(),
         nb::arg("ref_date") = nb::none(),
         nb::arg("ix_universe") = nb::none(),
         nb::arg("ix_d2_data") = nb::none(),
         nb::arg("ix_ref") = nb::none(),
         nb::arg("data_read_in") = nb::none(),
         nb::arg("ref_read_in") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoD2DataStruct::name,
          &TaoD2DataStruct::set_name,
          "Name to be used with commands."
      )
      .def_prop_rw(
          "data_file_name",
          &TaoD2DataStruct::data_file_name,
          &TaoD2DataStruct::set_data_file_name,
          "Data file name ."
      )
      .def_prop_rw(
          "ref_file_name",
          &TaoD2DataStruct::ref_file_name,
          &TaoD2DataStruct::set_ref_file_name,
          "Reference file name."
      )
      .def_prop_rw(
          "data_date",
          &TaoD2DataStruct::data_date,
          &TaoD2DataStruct::set_data_date,
          "Data measurement date."
      )
      .def_prop_rw(
          "ref_date",
          &TaoD2DataStruct::ref_date,
          &TaoD2DataStruct::set_ref_date,
          "Reference data measurement date."
      )
      .def_prop_ro(
          "descrip",
          &TaoD2DataStruct::descrip,
          nb::keep_alive<0, 1>(),
          "Array for descriptive information."
      )
      .def_prop_ro("d1", &TaoD2DataStruct::d1, nb::keep_alive<0, 1>(), "Points to children")
      .def_prop_rw(
          "ix_universe",
          &TaoD2DataStruct::ix_universe,
          &TaoD2DataStruct::set_ix_universe,
          "Index of universe this is in."
      )
      .def_prop_rw(
          "ix_d2_data",
          &TaoD2DataStruct::ix_d2_data,
          &TaoD2DataStruct::set_ix_d2_data,
          "Index in u%d2_data(:) array."
      )
      .def_prop_rw(
          "ix_ref",
          &TaoD2DataStruct::ix_ref,
          &TaoD2DataStruct::set_ix_ref,
          "Index of the reference data set."
      )
      .def_prop_rw(
          "data_read_in",
          &TaoD2DataStruct::data_read_in,
          &TaoD2DataStruct::set_data_read_in,
          "A data set has been read in?"
      )
      .def_prop_rw(
          "ref_read_in",
          &TaoD2DataStruct::ref_read_in,
          &TaoD2DataStruct::set_ref_read_in,
          "A reference data set has been read in?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoD2DataStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoD2DataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoD2DataStruct &self, nb::dict &memo) { return TaoD2DataStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoD2DataStruct &self, const TaoD2DataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_data_array_struct
void init_tao_data_array_struct(nb::module_ &m, nb::class_<TaoDataArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoDataArrayStruct *self, const TaoDataStruct *d) {
           new (self) TaoDataArrayStruct(ptr_to_opt_ref(d));
         },
         nb::arg("d") = nb::none()
  )
      .def_prop_rw(
          "d",
          &TaoDataArrayStruct::d,
          &TaoDataArrayStruct::set_d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoDataArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoDataArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoDataArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDataArrayStruct &self) {
            return TaoDataArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDataArrayStruct &self, nb::dict &memo) { return TaoDataArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoDataArrayStruct &self, const TaoDataArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoDataArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoDataArrayStructArray1D, TaoDataArrayStructAlloc1D>(
      m,
      "TaoDataArrayStructArray1D",
      "TaoDataArrayStructAlloc1D"
  );
  // 2D TaoDataArrayStruct arrays are not used in structs/routines
  // 3D TaoDataArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_data_struct
void init_tao_data_struct(nb::module_ &m, nb::class_<TaoDataStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoDataStruct *self,
            std::optional<std::string> ele_name,
            std::optional<std::string> ele_start_name,
            std::optional<std::string> ele_ref_name,
            std::optional<std::string> data_type,
            std::optional<std::string> merit_type,
            std::optional<std::string> id,
            std::optional<std::string> data_source,
            std::optional<std::string> why_invalid,
            std::optional<int> ix_uni,
            std::optional<int> ix_bunch,
            std::optional<int> ix_branch,
            std::optional<int> ix_ele,
            std::optional<int> ix_ele_start,
            std::optional<int> ix_ele_ref,
            std::optional<int> ix_ele_merit,
            std::optional<int> ix_d1,
            std::optional<int> ix_data,
            std::optional<int> ix_dModel,
            std::optional<int> eval_point,
            std::optional<double> meas_value,
            std::optional<double> ref_value,
            std::optional<double> model_value,
            std::optional<double> design_value,
            std::optional<double> old_value,
            std::optional<double> base_value,
            std::optional<double> error_rms,
            std::optional<double> delta_merit,
            std::optional<double> weight,
            std::optional<double> invalid_value,
            std::optional<double> merit,
            std::optional<double> s,
            std::optional<double> s_offset,
            std::optional<double> ref_s_offset,
            std::optional<bool> err_message_printed,
            std::optional<bool> exists,
            std::optional<bool> good_model,
            std::optional<bool> good_base,
            std::optional<bool> good_design,
            std::optional<bool> good_meas,
            std::optional<bool> good_ref,
            std::optional<bool> good_user,
            std::optional<bool> good_opt,
            std::optional<bool> good_plot,
            std::optional<bool> useit_plot,
            std::optional<bool> useit_opt,
            const TaoSpinMapStruct *spin_map,
            const TaoD1DataStruct *d1) {
           new (self) TaoDataStruct(
               ele_name,
               ele_start_name,
               ele_ref_name,
               data_type,
               merit_type,
               id,
               data_source,
               why_invalid,
               ix_uni,
               ix_bunch,
               ix_branch,
               ix_ele,
               ix_ele_start,
               ix_ele_ref,
               ix_ele_merit,
               ix_d1,
               ix_data,
               ix_dModel,
               eval_point,
               meas_value,
               ref_value,
               model_value,
               design_value,
               old_value,
               base_value,
               error_rms,
               delta_merit,
               weight,
               invalid_value,
               merit,
               s,
               s_offset,
               ref_s_offset,
               err_message_printed,
               exists,
               good_model,
               good_base,
               good_design,
               good_meas,
               good_ref,
               good_user,
               good_opt,
               good_plot,
               useit_plot,
               useit_opt,
               ptr_to_opt_ref(spin_map),
               ptr_to_opt_ref(d1)
           );
         },
         nb::arg("ele_name") = nb::none(),
         nb::arg("ele_start_name") = nb::none(),
         nb::arg("ele_ref_name") = nb::none(),
         nb::arg("data_type") = nb::none(),
         nb::arg("merit_type") = nb::none(),
         nb::arg("id") = nb::none(),
         nb::arg("data_source") = nb::none(),
         nb::arg("why_invalid") = nb::none(),
         nb::arg("ix_uni") = nb::none(),
         nb::arg("ix_bunch") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_ele_start") = nb::none(),
         nb::arg("ix_ele_ref") = nb::none(),
         nb::arg("ix_ele_merit") = nb::none(),
         nb::arg("ix_d1") = nb::none(),
         nb::arg("ix_data") = nb::none(),
         nb::arg("ix_dModel") = nb::none(),
         nb::arg("eval_point") = nb::none(),
         nb::arg("meas_value") = nb::none(),
         nb::arg("ref_value") = nb::none(),
         nb::arg("model_value") = nb::none(),
         nb::arg("design_value") = nb::none(),
         nb::arg("old_value") = nb::none(),
         nb::arg("base_value") = nb::none(),
         nb::arg("error_rms") = nb::none(),
         nb::arg("delta_merit") = nb::none(),
         nb::arg("weight") = nb::none(),
         nb::arg("invalid_value") = nb::none(),
         nb::arg("merit") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("s_offset") = nb::none(),
         nb::arg("ref_s_offset") = nb::none(),
         nb::arg("err_message_printed") = nb::none(),
         nb::arg("exists") = nb::none(),
         nb::arg("good_model") = nb::none(),
         nb::arg("good_base") = nb::none(),
         nb::arg("good_design") = nb::none(),
         nb::arg("good_meas") = nb::none(),
         nb::arg("good_ref") = nb::none(),
         nb::arg("good_user") = nb::none(),
         nb::arg("good_opt") = nb::none(),
         nb::arg("good_plot") = nb::none(),
         nb::arg("useit_plot") = nb::none(),
         nb::arg("useit_opt") = nb::none(),
         nb::arg("spin_map") = nb::none(),
         nb::arg("d1") = nb::none()
  )
      .def_prop_rw(
          "ele_name",
          &TaoDataStruct::ele_name,
          &TaoDataStruct::set_ele_name,
          "Name of the lattice element where datum is evaluated."
      )
      .def_prop_rw(
          "ele_start_name",
          &TaoDataStruct::ele_start_name,
          &TaoDataStruct::set_ele_start_name,
          "Name of starting lattice element when there is a range"
      )
      .def_prop_rw(
          "ele_ref_name",
          &TaoDataStruct::ele_ref_name,
          &TaoDataStruct::set_ele_ref_name,
          "Name of reference lattice element"
      )
      .def_prop_rw(
          "data_type",
          &TaoDataStruct::data_type,
          &TaoDataStruct::set_data_type,
          "Type of data: 'orbit.x', etc."
      )
      .def_prop_rw(
          "merit_type",
          &TaoDataStruct::merit_type,
          &TaoDataStruct::set_merit_type,
          "Type of constraint: 'target', 'max', 'min', etc."
      )
      .def_prop_rw(
          "id",
          &TaoDataStruct::id,
          &TaoDataStruct::set_id,
          "Used by Tao extension code. Not used by Tao directly."
      )
      .def_prop_rw(
          "data_source",
          &TaoDataStruct::data_source,
          &TaoDataStruct::set_data_source,
          "'lat', 'beam', 'data' or 'var'. Last two used for expressions."
      )
      .def_prop_rw(
          "why_invalid",
          &TaoDataStruct::why_invalid,
          &TaoDataStruct::set_why_invalid,
          "Informational string if there is a problem."
      )
      .def_prop_rw(
          "ix_uni",
          &TaoDataStruct::ix_uni,
          &TaoDataStruct::set_ix_uni,
          "Universe index of datum."
      )
      .def_prop_rw(
          "ix_bunch",
          &TaoDataStruct::ix_bunch,
          &TaoDataStruct::set_ix_bunch,
          "Bunch number to get the data from."
      )
      .def_prop_rw(
          "ix_branch",
          &TaoDataStruct::ix_branch,
          &TaoDataStruct::set_ix_branch,
          "Index of the associated lattice branch."
      )
      .def_prop_rw(
          "ix_ele",
          &TaoDataStruct::ix_ele,
          &TaoDataStruct::set_ix_ele,
          "Index of the lattice element corresponding to ele_name"
      )
      .def_prop_rw(
          "ix_ele_start",
          &TaoDataStruct::ix_ele_start,
          &TaoDataStruct::set_ix_ele_start,
          "Index of lattice elment when there is a range"
      )
      .def_prop_rw(
          "ix_ele_ref",
          &TaoDataStruct::ix_ele_ref,
          &TaoDataStruct::set_ix_ele_ref,
          "Index of lattice elment when there is a reference."
      )
      .def_prop_rw(
          "ix_ele_merit",
          &TaoDataStruct::ix_ele_merit,
          &TaoDataStruct::set_ix_ele_merit,
          "Index of lattice elment where merit is evaluated."
      )
      .def_prop_rw(
          "ix_d1",
          &TaoDataStruct::ix_d1,
          &TaoDataStruct::set_ix_d1,
          "Index number in u%d2_data(i)%d1_data(j)%d(:) array."
      )
      .def_prop_rw(
          "ix_data",
          &TaoDataStruct::ix_data,
          &TaoDataStruct::set_ix_data,
          "Index of this datum in the u%data(:) array of data_structs."
      )
      .def_prop_rw(
          "ix_dModel",
          &TaoDataStruct::ix_dModel,
          &TaoDataStruct::set_ix_dModel,
          "Row number in the dModel_dVar derivative matrix."
      )
      .def_prop_rw(
          "eval_point",
          &TaoDataStruct::eval_point,
          &TaoDataStruct::set_eval_point,
          "or anchor_center$, anchor_beginning$. Where to evaluate data relative to the element."
      )
      .def_prop_rw(
          "meas_value",
          &TaoDataStruct::meas_value,
          &TaoDataStruct::set_meas_value,
          "Measured datum value."
      )
      .def_prop_rw(
          "ref_value",
          &TaoDataStruct::ref_value,
          &TaoDataStruct::set_ref_value,
          "Measured datum value from the reference data set."
      )
      .def_prop_rw(
          "model_value",
          &TaoDataStruct::model_value,
          &TaoDataStruct::set_model_value,
          "Datum value as calculated from the model."
      )
      .def_prop_rw(
          "design_value",
          &TaoDataStruct::design_value,
          &TaoDataStruct::set_design_value,
          "What the datum value is in the design lattice."
      )
      .def_prop_rw(
          "old_value",
          &TaoDataStruct::old_value,
          &TaoDataStruct::set_old_value,
          "The model_value at some previous time."
      )
      .def_prop_rw(
          "base_value",
          &TaoDataStruct::base_value,
          &TaoDataStruct::set_base_value,
          "The value as calculated from the base model."
      )
      .def_prop_rw(
          "error_rms",
          &TaoDataStruct::error_rms,
          &TaoDataStruct::set_error_rms,
          "Measurement error RMS. Used in plotting."
      )
      .def_prop_rw(
          "delta_merit",
          &TaoDataStruct::delta_merit,
          &TaoDataStruct::set_delta_merit,
          "Diff used to calculate the merit function term."
      )
      .def_prop_rw(
          "weight",
          &TaoDataStruct::weight,
          &TaoDataStruct::set_weight,
          "Weight for the merit function term."
      )
      .def_prop_rw(
          "invalid_value",
          &TaoDataStruct::invalid_value,
          &TaoDataStruct::set_invalid_value,
          "Value used in merit calc if good_model = F (or possibly good_design & good_base)."
      )
      .def_prop_rw(
          "merit",
          &TaoDataStruct::merit,
          &TaoDataStruct::set_merit,
          "Merit function term value: weight * delta_merit^2"
      )
      .def_prop_rw("s", &TaoDataStruct::s, &TaoDataStruct::set_s, "longitudinal position of ele.")
      .def_prop_rw(
          "s_offset",
          &TaoDataStruct::s_offset,
          &TaoDataStruct::set_s_offset,
          "Offset of the evaluation point."
      )
      .def_prop_rw(
          "ref_s_offset",
          &TaoDataStruct::ref_s_offset,
          &TaoDataStruct::set_ref_s_offset,
          "Offset of the reference point. In development."
      )
      .def_prop_rw(
          "err_message_printed",
          &TaoDataStruct::err_message_printed,
          &TaoDataStruct::set_err_message_printed,
          "Used to prevent zillions of error messages being generated"
      )
      .def_prop_rw("exists", &TaoDataStruct::exists, &TaoDataStruct::set_exists, "See above")
      .def_prop_rw(
          "good_model",
          &TaoDataStruct::good_model,
          &TaoDataStruct::set_good_model,
          "See above"
      )
      .def_prop_rw(
          "good_base",
          &TaoDataStruct::good_base,
          &TaoDataStruct::set_good_base,
          "See above"
      )
      .def_prop_rw(
          "good_design",
          &TaoDataStruct::good_design,
          &TaoDataStruct::set_good_design,
          "See above"
      )
      .def_prop_rw(
          "good_meas",
          &TaoDataStruct::good_meas,
          &TaoDataStruct::set_good_meas,
          "See above"
      )
      .def_prop_rw("good_ref", &TaoDataStruct::good_ref, &TaoDataStruct::set_good_ref, "See above")
      .def_prop_rw(
          "good_user",
          &TaoDataStruct::good_user,
          &TaoDataStruct::set_good_user,
          "See above"
      )
      .def_prop_rw("good_opt", &TaoDataStruct::good_opt, &TaoDataStruct::set_good_opt, "See above")
      .def_prop_rw(
          "good_plot",
          &TaoDataStruct::good_plot,
          &TaoDataStruct::set_good_plot,
          "See above"
      )
      .def_prop_rw(
          "useit_plot",
          &TaoDataStruct::useit_plot,
          &TaoDataStruct::set_useit_plot,
          "See above"
      )
      .def_prop_rw(
          "useit_opt",
          &TaoDataStruct::useit_opt,
          &TaoDataStruct::set_useit_opt,
          "See above"
      )
      .def_prop_rw(
          "spin_map",
          &TaoDataStruct::spin_map,
          &TaoDataStruct::set_spin_map,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "d1",
          &TaoDataStruct::d1,
          &TaoDataStruct::set_d1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Pointer to the parent d1_data_struct"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoDataStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoDataStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoDataStruct &self, nb::dict &memo) { return TaoDataStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoDataStruct &self, const TaoDataStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_data_var_component_struct(
    nb::module_ &m,
    nb::class_<TaoDataVarComponentStruct> &cls
) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<double>>(),
         nb::arg("name") = nb::none(),
         nb::arg("sign") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoDataVarComponentStruct::name,
          &TaoDataVarComponentStruct::set_name,
          "Eg: 'meas', 'ref', 'model', etc."
      )
      .def_prop_rw(
          "sign",
          &TaoDataVarComponentStruct::sign,
          &TaoDataVarComponentStruct::set_sign,
          "+1 or -1"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoDataVarComponentStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoDataVarComponentStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoDataVarComponentStruct &self, nb::dict &memo) {
            return TaoDataVarComponentStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoDataVarComponentStruct &self, const TaoDataVarComponentStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_drawing_struct(nb::module_ &m, nb::class_<TaoDrawingStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("ele_shape", &TaoDrawingStruct::ele_shape, nb::keep_alive<0, 1>())

      .def("__repr__", [](const TaoDrawingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoDrawingStruct &self) {
            return TaoDrawingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoDrawingStruct &self, nb::dict &memo) { return TaoDrawingStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoDrawingStruct &self, const TaoDrawingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_dynamic_aperture_struct(nb::module_ &m, nb::class_<TaoDynamicApertureStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoDynamicApertureStruct *self,
            const ApertureParamStruct *param,
            std::optional<std::vector<double>> pz,
            std::optional<double> ellipse_scale,
            std::optional<double> a_emit,
            std::optional<double> b_emit) {
           new (self)
               TaoDynamicApertureStruct(ptr_to_opt_ref(param), pz, ellipse_scale, a_emit, b_emit);
         },
         nb::arg("param") = nb::none(),
         nb::arg("pz") = nb::none(),
         nb::arg("ellipse_scale") = nb::none(),
         nb::arg("a_emit") = nb::none(),
         nb::arg("b_emit") = nb::none()
  )
      .def_prop_rw(
          "param",
          &TaoDynamicApertureStruct::param,
          &TaoDynamicApertureStruct::set_param,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro(
          "scan",
          &TaoDynamicApertureStruct::scan,
          nb::keep_alive<0, 1>(),
          "One scan for each pz."
      )
      .def_prop_rw(
          "pz",
          &TaoDynamicApertureStruct::pz,
          &TaoDynamicApertureStruct::set_pz,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ellipse_scale",
          &TaoDynamicApertureStruct::ellipse_scale,
          &TaoDynamicApertureStruct::set_ellipse_scale
      )
      .def_prop_rw(
          "a_emit",
          &TaoDynamicApertureStruct::a_emit,
          &TaoDynamicApertureStruct::set_a_emit
      )
      .def_prop_rw(
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
          [](const TaoDynamicApertureStruct &self, nb::dict &memo) {
            return TaoDynamicApertureStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoDynamicApertureStruct &self, const TaoDynamicApertureStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_ele_pointer_struct(nb::module_ &m, nb::class_<TaoElePointerStruct> &cls) {
  cls.def(nb::init<std::optional<int>>(), nb::arg("n_loc") = nb::none())
      .def_prop_ro("eles", &TaoElePointerStruct::eles, nb::keep_alive<0, 1>())
      .def_prop_rw("n_loc", &TaoElePointerStruct::n_loc, &TaoElePointerStruct::set_n_loc)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoElePointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoElePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoElePointerStruct &self, nb::dict &memo) { return TaoElePointerStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoElePointerStruct &self, const TaoElePointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_ele_shape_struct(nb::module_ &m, nb::class_<TaoEleShapeStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("ele_id") = nb::none(),
         nb::arg("shape") = nb::none(),
         nb::arg("color") = nb::none(),
         nb::arg("size") = nb::none(),
         nb::arg("label") = nb::none(),
         nb::arg("draw") = nb::none(),
         nb::arg("multi") = nb::none(),
         nb::arg("line_width") = nb::none(),
         nb::arg("offset") = nb::none(),
         nb::arg("ix_key") = nb::none(),
         nb::arg("name_ele") = nb::none()
  )
      .def_prop_rw(
          "ele_id",
          &TaoEleShapeStruct::ele_id,
          &TaoEleShapeStruct::set_ele_id,
          "element 'key::name' to match to."
      )
      .def_prop_rw(
          "shape",
          &TaoEleShapeStruct::shape,
          &TaoEleShapeStruct::set_shape,
          "Shape to draw"
      )
      .def_prop_rw(
          "color",
          &TaoEleShapeStruct::color,
          &TaoEleShapeStruct::set_color,
          "Color of shape"
      )
      .def_prop_rw(
          "size",
          &TaoEleShapeStruct::size,
          &TaoEleShapeStruct::set_size,
          "plot vertical height"
      )
      .def_prop_rw(
          "label",
          &TaoEleShapeStruct::label,
          &TaoEleShapeStruct::set_label,
          "Can be: 'name', 's', 'none'"
      )
      .def_prop_rw(
          "draw",
          &TaoEleShapeStruct::draw,
          &TaoEleShapeStruct::set_draw,
          "Draw the shape?"
      )
      .def_prop_rw(
          "multi",
          &TaoEleShapeStruct::multi,
          &TaoEleShapeStruct::set_multi,
          "Can be part of a multi-shape."
      )
      .def_prop_rw(
          "line_width",
          &TaoEleShapeStruct::line_width,
          &TaoEleShapeStruct::set_line_width,
          "Width of lines used to draw the shape."
      )
      .def_prop_rw(
          "offset",
          &TaoEleShapeStruct::offset,
          &TaoEleShapeStruct::set_offset,
          "Vertical offset."
      )
      .def_prop_rw(
          "ix_key",
          &TaoEleShapeStruct::ix_key,
          &TaoEleShapeStruct::set_ix_key,
          "Extracted from ele_id. 0 => all classes (quadrupole, etc.)"
      )
      .def_prop_rw(
          "name_ele",
          &TaoEleShapeStruct::name_ele,
          &TaoEleShapeStruct::set_name_ele,
          "Name of element."
      )
      .def_prop_ro("uni", &TaoEleShapeStruct::uni, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return TaoEleShapeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoEleShapeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoEleShapeStruct &self, nb::dict &memo) { return TaoEleShapeStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoEleShapeStruct &self, const TaoEleShapeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_eval_node_struct(nb::module_ &m, nb::class_<TaoEvalNodeStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         nb::arg("type") = nb::none(),
         nb::arg("name") = nb::none(),
         nb::arg("scale") = nb::none(),
         nb::arg("value") = nb::none()
  )
      .def_prop_rw("type", &TaoEvalNodeStruct::type, &TaoEvalNodeStruct::set_type)
      .def_prop_rw("name", &TaoEvalNodeStruct::name, &TaoEvalNodeStruct::set_name)
      .def_prop_rw(
          "scale",
          &TaoEvalNodeStruct::scale,
          &TaoEvalNodeStruct::set_scale,
          "Scale factor for ping data"
      )
      .def_prop_rw(
          "value",
          &TaoEvalNodeStruct::value,
          &TaoEvalNodeStruct::set_value,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("info", &TaoEvalNodeStruct::info, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "value_ptr",
          &TaoEvalNodeStruct::value_ptr,
          nb::keep_alive<0, 1>(),
          "Used to point to data, lattice parameters, etc"
      )
      .def_prop_ro(
          "node",
          &TaoEvalNodeStruct::node,
          nb::keep_alive<0, 1>(),
          "Child nodes for tree construction."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoEvalNodeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoEvalNodeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoEvalNodeStruct &self, nb::dict &memo) { return TaoEvalNodeStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoEvalNodeStruct &self, const TaoEvalNodeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_expression_info_struct(nb::module_ &m, nb::class_<TaoExpressionInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoExpressionInfoStruct *self,
            std::optional<bool> good,
            const EleStruct *ele,
            std::optional<double> s) {
           new (self) TaoExpressionInfoStruct(good, ptr_to_opt_ref(ele), s);
         },
         nb::arg("good") = nb::none(),
         nb::arg("ele") = nb::none(),
         nb::arg("s") = nb::none()
  )
      .def_prop_rw(
          "good",
          &TaoExpressionInfoStruct::good,
          &TaoExpressionInfoStruct::set_good,
          "Expression is valid."
      )
      .def_prop_rw(
          "ele",
          &TaoExpressionInfoStruct::ele,
          &TaoExpressionInfoStruct::set_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Associated ele if it exists"
      )
      .def_prop_rw(
          "s",
          &TaoExpressionInfoStruct::s,
          &TaoExpressionInfoStruct::set_s,
          "Longitudinal position of expression."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoExpressionInfoStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoExpressionInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoExpressionInfoStruct &self, nb::dict &memo) {
            return TaoExpressionInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoExpressionInfoStruct &self, const TaoExpressionInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_floor_plan_struct(nb::module_ &m, nb::class_<TaoFloorPlanStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("view") = nb::none(),
         nb::arg("rotation") = nb::none(),
         nb::arg("correct_distortion") = nb::none(),
         nb::arg("flip_label_side") = nb::none(),
         nb::arg("size_is_absolute") = nb::none(),
         nb::arg("draw_only_first_pass") = nb::none(),
         nb::arg("draw_building_wall") = nb::none(),
         nb::arg("orbit_scale") = nb::none(),
         nb::arg("orbit_color") = nb::none(),
         nb::arg("orbit_pattern") = nb::none(),
         nb::arg("orbit_lattice") = nb::none(),
         nb::arg("orbit_width") = nb::none()
  )
      .def_prop_rw("view", &TaoFloorPlanStruct::view, &TaoFloorPlanStruct::set_view, "or 'xz'.")
      .def_prop_rw(
          "rotation",
          &TaoFloorPlanStruct::rotation,
          &TaoFloorPlanStruct::set_rotation,
          "Rotation of floor plan plot: 1.0 -> 360^deg"
      )
      .def_prop_rw(
          "correct_distortion",
          &TaoFloorPlanStruct::correct_distortion,
          &TaoFloorPlanStruct::set_correct_distortion,
          "T -> Shrink one axis so x-scale = y-scale."
      )
      .def_prop_rw(
          "flip_label_side",
          &TaoFloorPlanStruct::flip_label_side,
          &TaoFloorPlanStruct::set_flip_label_side,
          "Draw element label on other side of element?"
      )
      .def_prop_rw(
          "size_is_absolute",
          &TaoFloorPlanStruct::size_is_absolute,
          &TaoFloorPlanStruct::set_size_is_absolute,
          "Are shape sizes in meters or window pixels?"
      )
      .def_prop_rw(
          "draw_only_first_pass",
          &TaoFloorPlanStruct::draw_only_first_pass,
          &TaoFloorPlanStruct::set_draw_only_first_pass,
          "Draw only first pass with multipass elements?"
      )
      .def_prop_rw(
          "draw_building_wall",
          &TaoFloorPlanStruct::draw_building_wall,
          &TaoFloorPlanStruct::set_draw_building_wall,
          "Draw the building wall?"
      )
      .def_prop_rw(
          "orbit_scale",
          &TaoFloorPlanStruct::orbit_scale,
          &TaoFloorPlanStruct::set_orbit_scale,
          "Scale factor for drawing orbits. 0 -> Do not draw."
      )
      .def_prop_rw(
          "orbit_color",
          &TaoFloorPlanStruct::orbit_color,
          &TaoFloorPlanStruct::set_orbit_color
      )
      .def_prop_rw(
          "orbit_pattern",
          &TaoFloorPlanStruct::orbit_pattern,
          &TaoFloorPlanStruct::set_orbit_pattern
      )
      .def_prop_rw(
          "orbit_lattice",
          &TaoFloorPlanStruct::orbit_lattice,
          &TaoFloorPlanStruct::set_orbit_lattice,
          "Or 'design' or 'base'"
      )
      .def_prop_rw(
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
          [](const TaoFloorPlanStruct &self, nb::dict &memo) { return TaoFloorPlanStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoFloorPlanStruct &self, const TaoFloorPlanStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_global_struct(nb::module_ &m, nb::class_<TaoGlobalStruct> &cls) {
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
         nb::arg("beam_dead_cutoff") = nb::none(),
         nb::arg("lm_opt_deriv_reinit") = nb::none(),
         nb::arg("de_lm_step_ratio") = nb::none(),
         nb::arg("de_var_to_population_factor") = nb::none(),
         nb::arg("lmdif_eps") = nb::none(),
         nb::arg("lmdif_negligible_merit") = nb::none(),
         nb::arg("svd_cutoff") = nb::none(),
         nb::arg("unstable_penalty") = nb::none(),
         nb::arg("merit_stop_value") = nb::none(),
         nb::arg("dmerit_stop_value") = nb::none(),
         nb::arg("random_sigma_cutoff") = nb::none(),
         nb::arg("delta_e_chrom") = nb::none(),
         nb::arg("max_plot_time") = nb::none(),
         nb::arg("default_universe") = nb::none(),
         nb::arg("default_branch") = nb::none(),
         nb::arg("n_opti_cycles") = nb::none(),
         nb::arg("n_opti_loops") = nb::none(),
         nb::arg("n_threads") = nb::none(),
         nb::arg("phase_units") = nb::none(),
         nb::arg("bunch_to_plot") = nb::none(),
         nb::arg("random_seed") = nb::none(),
         nb::arg("n_top10_merit") = nb::none(),
         nb::arg("srdt_gen_n_slices") = nb::none(),
         nb::arg("datum_err_messages_max") = nb::none(),
         nb::arg("srdt_sxt_n_slices") = nb::none(),
         nb::arg("srdt_use_cache") = nb::none(),
         nb::arg("quiet") = nb::none(),
         nb::arg("random_engine") = nb::none(),
         nb::arg("random_gauss_converter") = nb::none(),
         nb::arg("track_type") = nb::none(),
         nb::arg("lat_sigma_calc_uses_emit_from") = nb::none(),
         nb::arg("prompt_string") = nb::none(),
         nb::arg("prompt_color") = nb::none(),
         nb::arg("optimizer") = nb::none(),
         nb::arg("print_command") = nb::none(),
         nb::arg("var_out_file") = nb::none(),
         nb::arg("history_file") = nb::none(),
         nb::arg("beam_timer_on") = nb::none(),
         nb::arg("box_plots") = nb::none(),
         nb::arg("blank_line_between_commands") = nb::none(),
         nb::arg("cmd_file_abort_on_error") = nb::none(),
         nb::arg("concatenate_maps") = nb::none(),
         nb::arg("derivative_recalc") = nb::none(),
         nb::arg("derivative_uses_design") = nb::none(),
         nb::arg("disable_smooth_line_calc") = nb::none(),
         nb::arg("draw_curve_off_scale_warn") = nb::none(),
         nb::arg("external_plotting") = nb::none(),
         nb::arg("label_lattice_elements") = nb::none(),
         nb::arg("label_keys") = nb::none(),
         nb::arg("lattice_calc_on") = nb::none(),
         nb::arg("only_limit_opt_vars") = nb::none(),
         nb::arg("opt_with_ref") = nb::none(),
         nb::arg("opt_with_base") = nb::none(),
         nb::arg("opt_match_auto_recalc") = nb::none(),
         nb::arg("opti_write_var_file") = nb::none(),
         nb::arg("optimizer_allow_user_abort") = nb::none(),
         nb::arg("optimizer_var_limit_warn") = nb::none(),
         nb::arg("plot_on") = nb::none(),
         nb::arg("rad_int_user_calc_on") = nb::none(),
         nb::arg("rf_on") = nb::none(),
         nb::arg("single_step") = nb::none(),
         nb::arg("stop_on_error") = nb::none(),
         nb::arg("svd_retreat_on_merit_increase") = nb::none(),
         nb::arg("var_limits_on") = nb::none(),
         nb::arg("wait_for_CR_in_single_mode") = nb::none(),
         nb::arg("symbol_import") = nb::none(),
         nb::arg("debug_on") = nb::none(),
         nb::arg("expression_tree_on") = nb::none(),
         nb::arg("verbose_on") = nb::none()
  )
      .def_prop_rw(
          "beam_dead_cutoff",
          &TaoGlobalStruct::beam_dead_cutoff,
          &TaoGlobalStruct::set_beam_dead_cutoff,
          "Percentage of dead particles at which beam tracking is stopped."
      )
      .def_prop_rw(
          "lm_opt_deriv_reinit",
          &TaoGlobalStruct::lm_opt_deriv_reinit,
          &TaoGlobalStruct::set_lm_opt_deriv_reinit,
          "Reinit derivative matrix cutoff"
      )
      .def_prop_rw(
          "de_lm_step_ratio",
          &TaoGlobalStruct::de_lm_step_ratio,
          &TaoGlobalStruct::set_de_lm_step_ratio,
          "Scaling for step sizes between DE and LM optimizers."
      )
      .def_prop_rw(
          "de_var_to_population_factor",
          &TaoGlobalStruct::de_var_to_population_factor,
          &TaoGlobalStruct::set_de_var_to_population_factor,
          "DE population = max(n_var*factor, 20)"
      )
      .def_prop_rw(
          "lmdif_eps",
          &TaoGlobalStruct::lmdif_eps,
          &TaoGlobalStruct::set_lmdif_eps,
          "Tollerance for lmdif optimizer."
      )
      .def_prop_rw(
          "lmdif_negligible_merit",
          &TaoGlobalStruct::lmdif_negligible_merit,
          &TaoGlobalStruct::set_lmdif_negligible_merit
      )
      .def_prop_rw(
          "svd_cutoff",
          &TaoGlobalStruct::svd_cutoff,
          &TaoGlobalStruct::set_svd_cutoff,
          "SVD singular value cutoff."
      )
      .def_prop_rw(
          "unstable_penalty",
          &TaoGlobalStruct::unstable_penalty,
          &TaoGlobalStruct::set_unstable_penalty,
          "Used in unstable_ring datum merit calculation."
      )
      .def_prop_rw(
          "merit_stop_value",
          &TaoGlobalStruct::merit_stop_value,
          &TaoGlobalStruct::set_merit_stop_value,
          "Merit value below which an optimizer will stop."
      )
      .def_prop_rw(
          "dmerit_stop_value",
          &TaoGlobalStruct::dmerit_stop_value,
          &TaoGlobalStruct::set_dmerit_stop_value,
          "Fractional Merit change below which an optimizer will stop."
      )
      .def_prop_rw(
          "random_sigma_cutoff",
          &TaoGlobalStruct::random_sigma_cutoff,
          &TaoGlobalStruct::set_random_sigma_cutoff,
          "Cut-off in sigmas."
      )
      .def_prop_rw(
          "delta_e_chrom",
          &TaoGlobalStruct::delta_e_chrom,
          &TaoGlobalStruct::set_delta_e_chrom,
          "Delta E used from chrom calc."
      )
      .def_prop_rw(
          "max_plot_time",
          &TaoGlobalStruct::max_plot_time,
          &TaoGlobalStruct::set_max_plot_time,
          "If plotting time (seconds) exceeds this than a message is generated."
      )
      .def_prop_rw(
          "default_universe",
          &TaoGlobalStruct::default_universe,
          &TaoGlobalStruct::set_default_universe,
          "Default universe to work with."
      )
      .def_prop_rw(
          "default_branch",
          &TaoGlobalStruct::default_branch,
          &TaoGlobalStruct::set_default_branch,
          "Default lattice branch to work with."
      )
      .def_prop_rw(
          "n_opti_cycles",
          &TaoGlobalStruct::n_opti_cycles,
          &TaoGlobalStruct::set_n_opti_cycles,
          "Number of optimization cycles"
      )
      .def_prop_rw(
          "n_opti_loops",
          &TaoGlobalStruct::n_opti_loops,
          &TaoGlobalStruct::set_n_opti_loops,
          "Number of optimization loops"
      )
      .def_prop_rw(
          "n_threads",
          &TaoGlobalStruct::n_threads,
          &TaoGlobalStruct::set_n_threads,
          "Number of OpenMP threads for parallel calculations."
      )
      .def_prop_rw(
          "phase_units",
          &TaoGlobalStruct::phase_units,
          &TaoGlobalStruct::set_phase_units,
          "Phase units on output."
      )
      .def_prop_rw(
          "bunch_to_plot",
          &TaoGlobalStruct::bunch_to_plot,
          &TaoGlobalStruct::set_bunch_to_plot,
          "Which bunch to plot"
      )
      .def_prop_rw(
          "random_seed",
          &TaoGlobalStruct::random_seed,
          &TaoGlobalStruct::set_random_seed,
          "Use system clock by default"
      )
      .def_prop_rw(
          "n_top10_merit",
          &TaoGlobalStruct::n_top10_merit,
          &TaoGlobalStruct::set_n_top10_merit,
          "Number of top merit constraints to print."
      )
      .def_prop_rw(
          "srdt_gen_n_slices",
          &TaoGlobalStruct::srdt_gen_n_slices,
          &TaoGlobalStruct::set_srdt_gen_n_slices,
          "Number times to slice elements for summation RDT calculation"
      )
      .def_prop_rw(
          "datum_err_messages_max",
          &TaoGlobalStruct::datum_err_messages_max,
          &TaoGlobalStruct::set_datum_err_messages_max,
          "Maximum number of error messages per call to lattice_calc."
      )
      .def_prop_rw(
          "srdt_sxt_n_slices",
          &TaoGlobalStruct::srdt_sxt_n_slices,
          &TaoGlobalStruct::set_srdt_sxt_n_slices,
          "Number times to slice sextupoles for summation RDT calculation"
      )
      .def_prop_rw(
          "srdt_use_cache",
          &TaoGlobalStruct::srdt_use_cache,
          &TaoGlobalStruct::set_srdt_use_cache,
          "Create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large."
      )
      .def_prop_rw(
          "quiet",
          &TaoGlobalStruct::quiet,
          &TaoGlobalStruct::set_quiet,
          "Print I/O when running a command file?"
      )
      .def_prop_rw(
          "random_engine",
          &TaoGlobalStruct::random_engine,
          &TaoGlobalStruct::set_random_engine,
          "Non-beam random number engine"
      )
      .def_prop_rw(
          "random_gauss_converter",
          &TaoGlobalStruct::random_gauss_converter,
          &TaoGlobalStruct::set_random_gauss_converter,
          "Non-beam"
      )
      .def_prop_rw(
          "track_type",
          &TaoGlobalStruct::track_type,
          &TaoGlobalStruct::set_track_type,
          "or 'beam'"
      )
      .def_prop_rw(
          "lat_sigma_calc_uses_emit_from",
          &TaoGlobalStruct::lat_sigma_calc_uses_emit_from,
          &TaoGlobalStruct::set_lat_sigma_calc_uses_emit_from,
          "Lattice derived sigma matrix uses emit values from where? Other possibilities: 'beam', "
          "'beam_init'."
      )
      .def_prop_rw(
          "prompt_string",
          &TaoGlobalStruct::prompt_string,
          &TaoGlobalStruct::set_prompt_string
      )
      .def_prop_rw(
          "prompt_color",
          &TaoGlobalStruct::prompt_color,
          &TaoGlobalStruct::set_prompt_color,
          "See read_a_line routine for possible settings."
      )
      .def_prop_rw(
          "optimizer",
          &TaoGlobalStruct::optimizer,
          &TaoGlobalStruct::set_optimizer,
          "optimizer to use."
      )
      .def_prop_rw(
          "print_command",
          &TaoGlobalStruct::print_command,
          &TaoGlobalStruct::set_print_command
      )
      .def_prop_rw(
          "var_out_file",
          &TaoGlobalStruct::var_out_file,
          &TaoGlobalStruct::set_var_out_file
      )
      .def_prop_rw(
          "history_file",
          &TaoGlobalStruct::history_file,
          &TaoGlobalStruct::set_history_file
      )
      .def_prop_rw(
          "beam_timer_on",
          &TaoGlobalStruct::beam_timer_on,
          &TaoGlobalStruct::set_beam_timer_on,
          "For timing the beam tracking calculation."
      )
      .def_prop_rw(
          "box_plots",
          &TaoGlobalStruct::box_plots,
          &TaoGlobalStruct::set_box_plots,
          "For debugging plot layout issues."
      )
      .def_prop_rw(
          "blank_line_between_commands",
          &TaoGlobalStruct::blank_line_between_commands,
          &TaoGlobalStruct::set_blank_line_between_commands,
          "Add a blank line between command output?"
      )
      .def_prop_rw(
          "cmd_file_abort_on_error",
          &TaoGlobalStruct::cmd_file_abort_on_error,
          &TaoGlobalStruct::set_cmd_file_abort_on_error,
          "Abort open command files if there is an error?"
      )
      .def_prop_rw(
          "concatenate_maps",
          &TaoGlobalStruct::concatenate_maps,
          &TaoGlobalStruct::set_concatenate_maps,
          "False => tracking using DA."
      )
      .def_prop_rw(
          "derivative_recalc",
          &TaoGlobalStruct::derivative_recalc,
          &TaoGlobalStruct::set_derivative_recalc,
          "Recalc before each optimizer run?"
      )
      .def_prop_rw(
          "derivative_uses_design",
          &TaoGlobalStruct::derivative_uses_design,
          &TaoGlobalStruct::set_derivative_uses_design,
          "Derivative calc uses design lattice instead of model?"
      )
      .def_prop_rw(
          "disable_smooth_line_calc",
          &TaoGlobalStruct::disable_smooth_line_calc,
          &TaoGlobalStruct::set_disable_smooth_line_calc,
          "Global disable of the smooth line calculation."
      )
      .def_prop_rw(
          "draw_curve_off_scale_warn",
          &TaoGlobalStruct::draw_curve_off_scale_warn,
          &TaoGlobalStruct::set_draw_curve_off_scale_warn,
          "Display warning on graphs?"
      )
      .def_prop_rw(
          "external_plotting",
          &TaoGlobalStruct::external_plotting,
          &TaoGlobalStruct::set_external_plotting,
          "Used with matplotlib and gui."
      )
      .def_prop_rw(
          "label_lattice_elements",
          &TaoGlobalStruct::label_lattice_elements,
          &TaoGlobalStruct::set_label_lattice_elements,
          "For lat_layout plots"
      )
      .def_prop_rw(
          "label_keys",
          &TaoGlobalStruct::label_keys,
          &TaoGlobalStruct::set_label_keys,
          "For lat_layout plots"
      )
      .def_prop_rw(
          "lattice_calc_on",
          &TaoGlobalStruct::lattice_calc_on,
          &TaoGlobalStruct::set_lattice_calc_on,
          "Turn on/off beam and single particle calculations."
      )
      .def_prop_rw(
          "only_limit_opt_vars",
          &TaoGlobalStruct::only_limit_opt_vars,
          &TaoGlobalStruct::set_only_limit_opt_vars,
          "Only apply limits to variables used in optimization."
      )
      .def_prop_rw(
          "opt_with_ref",
          &TaoGlobalStruct::opt_with_ref,
          &TaoGlobalStruct::set_opt_with_ref,
          "Use reference data in optimization?"
      )
      .def_prop_rw(
          "opt_with_base",
          &TaoGlobalStruct::opt_with_base,
          &TaoGlobalStruct::set_opt_with_base,
          "Use base data in optimization?"
      )
      .def_prop_rw(
          "opt_match_auto_recalc",
          &TaoGlobalStruct::opt_match_auto_recalc,
          &TaoGlobalStruct::set_opt_match_auto_recalc,
          "Set recalc = True for match elements before each cycle?"
      )
      .def_prop_rw(
          "opti_write_var_file",
          &TaoGlobalStruct::opti_write_var_file,
          &TaoGlobalStruct::set_opti_write_var_file,
          "'run' command writes var_out_file"
      )
      .def_prop_rw(
          "optimizer_allow_user_abort",
          &TaoGlobalStruct::optimizer_allow_user_abort,
          &TaoGlobalStruct::set_optimizer_allow_user_abort,
          "See Tao manual for more details."
      )
      .def_prop_rw(
          "optimizer_var_limit_warn",
          &TaoGlobalStruct::optimizer_var_limit_warn,
          &TaoGlobalStruct::set_optimizer_var_limit_warn,
          "Warn when vars reach a limit with optimization."
      )
      .def_prop_rw(
          "plot_on",
          &TaoGlobalStruct::plot_on,
          &TaoGlobalStruct::set_plot_on,
          "Do plotting?"
      )
      .def_prop_rw(
          "rad_int_user_calc_on",
          &TaoGlobalStruct::rad_int_user_calc_on,
          &TaoGlobalStruct::set_rad_int_user_calc_on,
          "User set radiation integrals calculation on/off."
      )
      .def_prop_rw(
          "rf_on",
          &TaoGlobalStruct::rf_on,
          &TaoGlobalStruct::set_rf_on,
          "RFcavities on or off? Does not affect lcavities."
      )
      .def_prop_rw(
          "single_step",
          &TaoGlobalStruct::single_step,
          &TaoGlobalStruct::set_single_step,
          "For debugging and demonstrations: Single step through a command file?"
      )
      .def_prop_rw(
          "stop_on_error",
          &TaoGlobalStruct::stop_on_error,
          &TaoGlobalStruct::set_stop_on_error,
          "For debugging: False prevents tao from exiting on an error."
      )
      .def_prop_rw(
          "svd_retreat_on_merit_increase",
          &TaoGlobalStruct::svd_retreat_on_merit_increase,
          &TaoGlobalStruct::set_svd_retreat_on_merit_increase
      )
      .def_prop_rw(
          "var_limits_on",
          &TaoGlobalStruct::var_limits_on,
          &TaoGlobalStruct::set_var_limits_on,
          "Respect the variable limits?"
      )
      .def_prop_rw(
          "wait_for_CR_in_single_mode",
          &TaoGlobalStruct::wait_for_CR_in_single_mode,
          &TaoGlobalStruct::set_wait_for_CR_in_single_mode,
          "For use with a python GUI."
      )
      .def_prop_rw(
          "symbol_import",
          &TaoGlobalStruct::symbol_import,
          &TaoGlobalStruct::set_symbol_import,
          "Import symbols from lattice file(s)? Internal stuff"
      )
      .def_prop_rw(
          "debug_on",
          &TaoGlobalStruct::debug_on,
          &TaoGlobalStruct::set_debug_on,
          "For debugging."
      )
      .def_prop_rw(
          "expression_tree_on",
          &TaoGlobalStruct::expression_tree_on,
          &TaoGlobalStruct::set_expression_tree_on,
          "Use an expression tree instead of a stack?"
      )
      .def_prop_rw(
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
          [](const TaoGlobalStruct &self, nb::dict &memo) { return TaoGlobalStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoGlobalStruct &self, const TaoGlobalStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_graph_array_struct
void init_tao_graph_array_struct(nb::module_ &m, nb::class_<TaoGraphArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoGraphArrayStruct *self, const TaoGraphStruct *g) {
           new (self) TaoGraphArrayStruct(ptr_to_opt_ref(g));
         },
         nb::arg("g") = nb::none()
  )
      .def_prop_rw(
          "g",
          &TaoGraphArrayStruct::g,
          &TaoGraphArrayStruct::set_g,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoGraphArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoGraphArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoGraphArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoGraphArrayStruct &self) {
            return TaoGraphArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoGraphArrayStruct &self, nb::dict &memo) { return TaoGraphArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoGraphArrayStruct &self, const TaoGraphArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoGraphArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoGraphArrayStructArray1D, TaoGraphArrayStructAlloc1D>(
      m,
      "TaoGraphArrayStructArray1D",
      "TaoGraphArrayStructAlloc1D"
  );
  // 2D TaoGraphArrayStruct arrays are not used in structs/routines
  // 3D TaoGraphArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_graph_struct
void init_tao_graph_struct(nb::module_ &m, nb::class_<TaoGraphStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoGraphStruct *self,
            std::optional<std::string> name,
            std::optional<std::string> type,
            std::optional<std::string> title,
            std::optional<std::string> title_suffix,
            std::optional<std::string> why_invalid,
            const TaoPlotStruct *p,
            const TaoFloorPlanStruct *floor_plan,
            const QpPointStruct *text_legend_origin,
            const QpPointStruct *curve_legend_origin,
            const QpLegendStruct *curve_legend,
            const QpAxisStruct *x,
            const QpAxisStruct *y,
            const QpAxisStruct *x2,
            const QpAxisStruct *y2,
            const QpRectStruct *margin,
            const QpRectStruct *scale_margin,
            std::optional<double> x_axis_scale_factor,
            std::optional<double> symbol_size_scale,
            std::optional<std::vector<int>> box,
            std::optional<int> ix_branch,
            std::optional<int> ix_universe,
            std::optional<bool> clip,
            std::optional<bool> y2_mirrors_y,
            std::optional<bool> limited,
            std::optional<bool> draw_axes,
            std::optional<bool> draw_curve_legend,
            std::optional<bool> draw_grid,
            std::optional<bool> draw_title,
            std::optional<bool> draw_only_good_user_data_or_vars,
            std::optional<bool> allow_wrap_around,
            std::optional<bool> is_valid) {
           new (self) TaoGraphStruct(
               name,
               type,
               title,
               title_suffix,
               why_invalid,
               ptr_to_opt_ref(p),
               ptr_to_opt_ref(floor_plan),
               ptr_to_opt_ref(text_legend_origin),
               ptr_to_opt_ref(curve_legend_origin),
               ptr_to_opt_ref(curve_legend),
               ptr_to_opt_ref(x),
               ptr_to_opt_ref(y),
               ptr_to_opt_ref(x2),
               ptr_to_opt_ref(y2),
               ptr_to_opt_ref(margin),
               ptr_to_opt_ref(scale_margin),
               x_axis_scale_factor,
               symbol_size_scale,
               box,
               ix_branch,
               ix_universe,
               clip,
               y2_mirrors_y,
               limited,
               draw_axes,
               draw_curve_legend,
               draw_grid,
               draw_title,
               draw_only_good_user_data_or_vars,
               allow_wrap_around,
               is_valid
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("title") = nb::none(),
         nb::arg("title_suffix") = nb::none(),
         nb::arg("why_invalid") = nb::none(),
         nb::arg("p") = nb::none(),
         nb::arg("floor_plan") = nb::none(),
         nb::arg("text_legend_origin") = nb::none(),
         nb::arg("curve_legend_origin") = nb::none(),
         nb::arg("curve_legend") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("x2") = nb::none(),
         nb::arg("y2") = nb::none(),
         nb::arg("margin") = nb::none(),
         nb::arg("scale_margin") = nb::none(),
         nb::arg("x_axis_scale_factor") = nb::none(),
         nb::arg("symbol_size_scale") = nb::none(),
         nb::arg("box") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_universe") = nb::none(),
         nb::arg("clip") = nb::none(),
         nb::arg("y2_mirrors_y") = nb::none(),
         nb::arg("limited") = nb::none(),
         nb::arg("draw_axes") = nb::none(),
         nb::arg("draw_curve_legend") = nb::none(),
         nb::arg("draw_grid") = nb::none(),
         nb::arg("draw_title") = nb::none(),
         nb::arg("draw_only_good_user_data_or_vars") = nb::none(),
         nb::arg("allow_wrap_around") = nb::none(),
         nb::arg("is_valid") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoGraphStruct::name,
          &TaoGraphStruct::set_name,
          "Name identifying the graph"
      )
      .def_prop_rw(
          "type",
          &TaoGraphStruct::type,
          &TaoGraphStruct::set_type,
          "'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture'"
      )
      .def_prop_rw("title", &TaoGraphStruct::title, &TaoGraphStruct::set_title)
      .def_prop_rw("title_suffix", &TaoGraphStruct::title_suffix, &TaoGraphStruct::set_title_suffix)
      .def_prop_ro(
          "text_legend",
          &TaoGraphStruct::text_legend,
          nb::keep_alive<0, 1>(),
          "Array for holding descriptive info."
      )
      .def_prop_ro(
          "text_legend_out",
          &TaoGraphStruct::text_legend_out,
          nb::keep_alive<0, 1>(),
          "Array for holding descriptive info."
      )
      .def_prop_rw(
          "why_invalid",
          &TaoGraphStruct::why_invalid,
          &TaoGraphStruct::set_why_invalid,
          "Informative string to print."
      )
      .def_prop_ro("curve", &TaoGraphStruct::curve, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "p",
          &TaoGraphStruct::p,
          &TaoGraphStruct::set_p,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "pointer to parent plot"
      )
      .def_prop_rw(
          "floor_plan",
          &TaoGraphStruct::floor_plan,
          &TaoGraphStruct::set_floor_plan,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "text_legend_origin",
          &TaoGraphStruct::text_legend_origin,
          &TaoGraphStruct::set_text_legend_origin,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "curve_legend_origin",
          &TaoGraphStruct::curve_legend_origin,
          &TaoGraphStruct::set_curve_legend_origin,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "curve_legend",
          &TaoGraphStruct::curve_legend,
          &TaoGraphStruct::set_curve_legend,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "x",
          &TaoGraphStruct::x,
          &TaoGraphStruct::set_x,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "X-axis parameters."
      )
      .def_prop_rw(
          "y",
          &TaoGraphStruct::y,
          &TaoGraphStruct::set_y,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Y-axis attributes."
      )
      .def_prop_rw(
          "x2",
          &TaoGraphStruct::x2,
          &TaoGraphStruct::set_x2,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "X2-axis attributes (Not currently used)."
      )
      .def_prop_rw(
          "y2",
          &TaoGraphStruct::y2,
          &TaoGraphStruct::set_y2,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Y2-axis attributes."
      )
      .def_prop_rw(
          "margin",
          &TaoGraphStruct::margin,
          &TaoGraphStruct::set_margin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Margin around the graph."
      )
      .def_prop_rw(
          "scale_margin",
          &TaoGraphStruct::scale_margin,
          &TaoGraphStruct::set_scale_margin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Margin for scaling"
      )
      .def_prop_rw(
          "x_axis_scale_factor",
          &TaoGraphStruct::x_axis_scale_factor,
          &TaoGraphStruct::set_x_axis_scale_factor,
          "x-axis conversion from internal to plotting units."
      )
      .def_prop_rw(
          "symbol_size_scale",
          &TaoGraphStruct::symbol_size_scale,
          &TaoGraphStruct::set_symbol_size_scale,
          "Symbol size scale factor for phase_space plots."
      )
      .def_prop_rw(
          "box",
          &TaoGraphStruct::box,
          &TaoGraphStruct::set_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Defines which box the plot is put in."
      )
      .def_prop_rw(
          "ix_branch",
          &TaoGraphStruct::ix_branch,
          &TaoGraphStruct::set_ix_branch,
          "Branch in lattice. Used when there are no associated curves."
      )
      .def_prop_rw(
          "ix_universe",
          &TaoGraphStruct::ix_universe,
          &TaoGraphStruct::set_ix_universe,
          "Used for lat_layout plots."
      )
      .def_prop_rw(
          "clip",
          &TaoGraphStruct::clip,
          &TaoGraphStruct::set_clip,
          "Clip plot at graph boundary."
      )
      .def_prop_rw(
          "y2_mirrors_y",
          &TaoGraphStruct::y2_mirrors_y,
          &TaoGraphStruct::set_y2_mirrors_y,
          "Y2-axis same as Y-axis?"
      )
      .def_prop_rw(
          "limited",
          &TaoGraphStruct::limited,
          &TaoGraphStruct::set_limited,
          "True if at least one data point past graph bounds."
      )
      .def_prop_rw(
          "draw_axes",
          &TaoGraphStruct::draw_axes,
          &TaoGraphStruct::set_draw_axes,
          "Draw axes, labels, etc?"
      )
      .def_prop_rw(
          "draw_curve_legend",
          &TaoGraphStruct::draw_curve_legend,
          &TaoGraphStruct::set_draw_curve_legend,
          "Legend for displaying curve info."
      )
      .def_prop_rw(
          "draw_grid",
          &TaoGraphStruct::draw_grid,
          &TaoGraphStruct::set_draw_grid,
          "Draw a grid?"
      )
      .def_prop_rw("draw_title", &TaoGraphStruct::draw_title, &TaoGraphStruct::set_draw_title)
      .def_prop_rw(
          "draw_only_good_user_data_or_vars",
          &TaoGraphStruct::draw_only_good_user_data_or_vars,
          &TaoGraphStruct::set_draw_only_good_user_data_or_vars
      )
      .def_prop_rw(
          "allow_wrap_around",
          &TaoGraphStruct::allow_wrap_around,
          &TaoGraphStruct::set_allow_wrap_around,
          "'Wrap' curves to extend past lattice boundaries?"
      )
      .def_prop_rw(
          "is_valid",
          &TaoGraphStruct::is_valid,
          &TaoGraphStruct::set_is_valid,
          "EG: Bad x_axis_type."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoGraphStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoGraphStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoGraphStruct &self, nb::dict &memo) { return TaoGraphStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoGraphStruct &self, const TaoGraphStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_histogram_struct(nb::module_ &m, nb::class_<TaoHistogramStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("density_normalized") = nb::none(),
         nb::arg("weight_by_charge") = nb::none(),
         nb::arg("minimum") = nb::none(),
         nb::arg("maximum") = nb::none(),
         nb::arg("width") = nb::none(),
         nb::arg("center") = nb::none(),
         nb::arg("number") = nb::none()
  )
      .def_prop_rw(
          "density_normalized",
          &TaoHistogramStruct::density_normalized,
          &TaoHistogramStruct::set_density_normalized
      )
      .def_prop_rw(
          "weight_by_charge",
          &TaoHistogramStruct::weight_by_charge,
          &TaoHistogramStruct::set_weight_by_charge
      )
      .def_prop_rw(
          "minimum",
          &TaoHistogramStruct::minimum,
          &TaoHistogramStruct::set_minimum,
          "Computed by Tao. Not User settable."
      )
      .def_prop_rw(
          "maximum",
          &TaoHistogramStruct::maximum,
          &TaoHistogramStruct::set_maximum,
          "Computed by Tao. Not User settable."
      )
      .def_prop_rw("width", &TaoHistogramStruct::width, &TaoHistogramStruct::set_width)
      .def_prop_rw("center", &TaoHistogramStruct::center, &TaoHistogramStruct::set_center)
      .def_prop_rw("number", &TaoHistogramStruct::number, &TaoHistogramStruct::set_number)

      .def("__repr__", [](const TaoHistogramStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoHistogramStruct &self) {
            return TaoHistogramStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoHistogramStruct &self, nb::dict &memo) { return TaoHistogramStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoHistogramStruct &self, const TaoHistogramStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_init_struct(nb::module_ &m, nb::class_<TaoInitStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("parse_cmd_args") = nb::none(),
         nb::arg("debug_switch") = nb::none(),
         nb::arg("external_plotting_switch") = nb::none(),
         nb::arg("init_name") = nb::none(),
         nb::arg("hook_init_file") = nb::none(),
         nb::arg("hook_lat_file") = nb::none(),
         nb::arg("hook_beam_file") = nb::none(),
         nb::arg("hook_data_file") = nb::none(),
         nb::arg("hook_plot_file") = nb::none(),
         nb::arg("hook_startup_file") = nb::none(),
         nb::arg("hook_var_file") = nb::none(),
         nb::arg("hook_building_wall_file") = nb::none(),
         nb::arg("init_file_arg_path") = nb::none(),
         nb::arg("lattice_file_arg") = nb::none(),
         nb::arg("hook_init_file_arg") = nb::none(),
         nb::arg("init_file_arg") = nb::none(),
         nb::arg("beam_file_arg") = nb::none(),
         nb::arg("beam_init_position_file_arg") = nb::none(),
         nb::arg("command_arg") = nb::none(),
         nb::arg("data_file_arg") = nb::none(),
         nb::arg("plot_file_arg") = nb::none(),
         nb::arg("startup_file_arg") = nb::none(),
         nb::arg("var_file_arg") = nb::none(),
         nb::arg("building_wall_file_arg") = nb::none(),
         nb::arg("geometry_arg") = nb::none(),
         nb::arg("slice_lattice_arg") = nb::none(),
         nb::arg("start_branch_at_arg") = nb::none(),
         nb::arg("log_startup_arg") = nb::none(),
         nb::arg("no_stopping_arg") = nb::none(),
         nb::arg("noplot_arg") = nb::none(),
         nb::arg("no_rad_int_arg") = nb::none(),
         nb::arg("reverse_arg") = nb::none(),
         nb::arg("debug_arg") = nb::none(),
         nb::arg("disable_smooth_line_calc_arg") = nb::none(),
         nb::arg("rf_on_arg") = nb::none(),
         nb::arg("prompt_color_arg") = nb::none(),
         nb::arg("quiet_arg") = nb::none(),
         nb::arg("noinit_arg") = nb::none(),
         nb::arg("nostartup_arg") = nb::none(),
         nb::arg("symbol_import_arg") = nb::none(),
         nb::arg("unique_name_suffix") = nb::none()
  )
      .def_prop_rw(
          "parse_cmd_args",
          &TaoInitStruct::parse_cmd_args,
          &TaoInitStruct::set_parse_cmd_args,
          "Used by custom programs to control Tao init"
      )
      .def_prop_rw(
          "debug_switch",
          &TaoInitStruct::debug_switch,
          &TaoInitStruct::set_debug_switch,
          "Is the '-debug' switch present?"
      )
      .def_prop_rw(
          "external_plotting_switch",
          &TaoInitStruct::external_plotting_switch,
          &TaoInitStruct::set_external_plotting_switch,
          "Is '-external_plotting' switch present?"
      )
      .def_prop_rw(
          "init_name",
          &TaoInitStruct::init_name,
          &TaoInitStruct::set_init_name,
          "label for initialization"
      )
      .def_prop_rw(
          "hook_init_file",
          &TaoInitStruct::hook_init_file,
          &TaoInitStruct::set_hook_init_file
      )
      .def_prop_rw(
          "hook_lat_file",
          &TaoInitStruct::hook_lat_file,
          &TaoInitStruct::set_hook_lat_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_beam_file",
          &TaoInitStruct::hook_beam_file,
          &TaoInitStruct::set_hook_beam_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_data_file",
          &TaoInitStruct::hook_data_file,
          &TaoInitStruct::set_hook_data_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_plot_file",
          &TaoInitStruct::hook_plot_file,
          &TaoInitStruct::set_hook_plot_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_startup_file",
          &TaoInitStruct::hook_startup_file,
          &TaoInitStruct::set_hook_startup_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_var_file",
          &TaoInitStruct::hook_var_file,
          &TaoInitStruct::set_hook_var_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "hook_building_wall_file",
          &TaoInitStruct::hook_building_wall_file,
          &TaoInitStruct::set_hook_building_wall_file,
          "To be set by tao_hook_parse_command_args"
      )
      .def_prop_rw(
          "init_file_arg_path",
          &TaoInitStruct::init_file_arg_path,
          &TaoInitStruct::set_init_file_arg_path,
          "Path part of init_tao_file"
      )
      .def_prop_rw(
          "lattice_file_arg",
          &TaoInitStruct::lattice_file_arg,
          &TaoInitStruct::set_lattice_file_arg,
          "-lattice_file        command line argument."
      )
      .def_prop_rw(
          "hook_init_file_arg",
          &TaoInitStruct::hook_init_file_arg,
          &TaoInitStruct::set_hook_init_file_arg,
          "-hook_init_file      command line argument"
      )
      .def_prop_rw(
          "init_file_arg",
          &TaoInitStruct::init_file_arg,
          &TaoInitStruct::set_init_file_arg,
          "-init_file           command line argument."
      )
      .def_prop_rw(
          "beam_file_arg",
          &TaoInitStruct::beam_file_arg,
          &TaoInitStruct::set_beam_file_arg,
          "-beam_file           command line argument."
      )
      .def_prop_rw(
          "beam_init_position_file_arg",
          &TaoInitStruct::beam_init_position_file_arg,
          &TaoInitStruct::set_beam_init_position_file_arg,
          "-beam_init_position_file command line argument."
      )
      .def_prop_rw(
          "command_arg",
          &TaoInitStruct::command_arg,
          &TaoInitStruct::set_command_arg,
          "-command             command line argument."
      )
      .def_prop_rw(
          "data_file_arg",
          &TaoInitStruct::data_file_arg,
          &TaoInitStruct::set_data_file_arg,
          "-data_file           command line argument."
      )
      .def_prop_rw(
          "plot_file_arg",
          &TaoInitStruct::plot_file_arg,
          &TaoInitStruct::set_plot_file_arg,
          "-plot_file           command line argument."
      )
      .def_prop_rw(
          "startup_file_arg",
          &TaoInitStruct::startup_file_arg,
          &TaoInitStruct::set_startup_file_arg,
          "-startup_file        command line argument."
      )
      .def_prop_rw(
          "var_file_arg",
          &TaoInitStruct::var_file_arg,
          &TaoInitStruct::set_var_file_arg,
          "-var_file            command line argument."
      )
      .def_prop_rw(
          "building_wall_file_arg",
          &TaoInitStruct::building_wall_file_arg,
          &TaoInitStruct::set_building_wall_file_arg,
          "-building_wall_file  command line argument."
      )
      .def_prop_rw(
          "geometry_arg",
          &TaoInitStruct::geometry_arg,
          &TaoInitStruct::set_geometry_arg,
          "-geometry            command line argument."
      )
      .def_prop_rw(
          "slice_lattice_arg",
          &TaoInitStruct::slice_lattice_arg,
          &TaoInitStruct::set_slice_lattice_arg,
          "-slice_lattice       command line argument."
      )
      .def_prop_rw(
          "start_branch_at_arg",
          &TaoInitStruct::start_branch_at_arg,
          &TaoInitStruct::set_start_branch_at_arg,
          "-start_branch_at     command line argument."
      )
      .def_prop_rw(
          "log_startup_arg",
          &TaoInitStruct::log_startup_arg,
          &TaoInitStruct::set_log_startup_arg,
          "-log_startup         command line argument"
      )
      .def_prop_rw(
          "no_stopping_arg",
          &TaoInitStruct::no_stopping_arg,
          &TaoInitStruct::set_no_stopping_arg,
          "-no_stopping         command line argument"
      )
      .def_prop_rw(
          "noplot_arg",
          &TaoInitStruct::noplot_arg,
          &TaoInitStruct::set_noplot_arg,
          "-noplot              command line argument"
      )
      .def_prop_rw(
          "no_rad_int_arg",
          &TaoInitStruct::no_rad_int_arg,
          &TaoInitStruct::set_no_rad_int_arg,
          "-no_rad_int          command line argument"
      )
      .def_prop_rw(
          "reverse_arg",
          &TaoInitStruct::reverse_arg,
          &TaoInitStruct::set_reverse_arg,
          "-reverse             command line argument"
      )
      .def_prop_rw(
          "debug_arg",
          &TaoInitStruct::debug_arg,
          &TaoInitStruct::set_debug_arg,
          "-debug               command line argument"
      )
      .def_prop_rw(
          "disable_smooth_line_calc_arg",
          &TaoInitStruct::disable_smooth_line_calc_arg,
          &TaoInitStruct::set_disable_smooth_line_calc_arg,
          "-disable_smooth_line_calc"
      )
      .def_prop_rw(
          "rf_on_arg",
          &TaoInitStruct::rf_on_arg,
          &TaoInitStruct::set_rf_on_arg,
          "-rf_on               command line argument"
      )
      .def_prop_rw(
          "prompt_color_arg",
          &TaoInitStruct::prompt_color_arg,
          &TaoInitStruct::set_prompt_color_arg,
          "-prompt_color        command line argument"
      )
      .def_prop_rw(
          "quiet_arg",
          &TaoInitStruct::quiet_arg,
          &TaoInitStruct::set_quiet_arg,
          "-quiet               command line argument"
      )
      .def_prop_rw(
          "noinit_arg",
          &TaoInitStruct::noinit_arg,
          &TaoInitStruct::set_noinit_arg,
          "-noinit              command line argument"
      )
      .def_prop_rw(
          "nostartup_arg",
          &TaoInitStruct::nostartup_arg,
          &TaoInitStruct::set_nostartup_arg,
          "-nostartup           command line argument"
      )
      .def_prop_rw(
          "symbol_import_arg",
          &TaoInitStruct::symbol_import_arg,
          &TaoInitStruct::set_symbol_import_arg,
          "-symbol_import       command line argument"
      )
      .def_prop_rw(
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
          [](const TaoInitStruct &self, nb::dict &memo) { return TaoInitStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoInitStruct &self, const TaoInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_integer_array_struct
void init_tao_integer_array_struct(nb::module_ &m, nb::class_<TaoIntegerArrayStruct> &cls) {
  cls.def(nb::init<std::optional<int>>(), nb::arg("i") = nb::none())
      .def_prop_rw("i", &TaoIntegerArrayStruct::i, &TaoIntegerArrayStruct::set_i)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoIntegerArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoIntegerArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoIntegerArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoIntegerArrayStruct &self) {
            return TaoIntegerArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoIntegerArrayStruct &self, nb::dict &memo) {
            return TaoIntegerArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoIntegerArrayStruct &self, const TaoIntegerArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoIntegerArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoIntegerArrayStructArray1D, TaoIntegerArrayStructAlloc1D>(
      m,
      "TaoIntegerArrayStructArray1D",
      "TaoIntegerArrayStructAlloc1D"
  );
  // 2D TaoIntegerArrayStruct arrays are not used in structs/routines
  // 3D TaoIntegerArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_lat_sigma_struct
void init_tao_lat_sigma_struct(nb::module_ &m, nb::class_<TaoLatSigmaStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<std::vector<double>>>>(), nb::arg("mat") = nb::none())
      .def_prop_rw(
          "mat",
          &TaoLatSigmaStruct::mat,
          &TaoLatSigmaStruct::set_mat,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoLatSigmaStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoLatSigmaStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoLatSigmaStruct &self, nb::dict &memo) { return TaoLatSigmaStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoLatSigmaStruct &self, const TaoLatSigmaStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_lattice_branch_struct(nb::module_ &m, nb::class_<TaoLatticeBranchStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoLatticeBranchStruct *self,
            const TaoLatticeStruct *tao_lat,
            const TaoSpinPolarizationStruct *spin,
            const SummationRdtStruct *srdt,
            const CoordStruct *orb0,
            const NormalModesStruct *modes_ri,
            const NormalModesStruct *modes_6d,
            const PtcNormalFormStruct *ptc_normal_form,
            const BmadNormalFormStruct *bmad_normal_form,
            std::optional<double> cache_x_min,
            std::optional<double> cache_x_max,
            std::optional<double> comb_ds_save,
            std::optional<int> ix_ref_taylor,
            std::optional<int> ix_ele_taylor,
            std::optional<int> track_state,
            std::optional<int> cache_n_pts,
            std::optional<int> ix_rad_int_cache,
            std::optional<bool> has_open_match_element,
            std::optional<bool> plot_cache_valid,
            std::optional<bool> spin_map_valid,
            std::optional<bool> twiss_valid,
            std::optional<bool> mode_flip_here,
            std::optional<bool> chrom_calc_ok,
            std::optional<bool> rad_int_calc_ok,
            std::optional<bool> emit_6d_calc_ok,
            std::optional<bool> sigma_track_ok) {
           new (self) TaoLatticeBranchStruct(
               ptr_to_opt_ref(tao_lat),
               ptr_to_opt_ref(spin),
               ptr_to_opt_ref(srdt),
               ptr_to_opt_ref(orb0),
               ptr_to_opt_ref(modes_ri),
               ptr_to_opt_ref(modes_6d),
               ptr_to_opt_ref(ptc_normal_form),
               ptr_to_opt_ref(bmad_normal_form),
               cache_x_min,
               cache_x_max,
               comb_ds_save,
               ix_ref_taylor,
               ix_ele_taylor,
               track_state,
               cache_n_pts,
               ix_rad_int_cache,
               has_open_match_element,
               plot_cache_valid,
               spin_map_valid,
               twiss_valid,
               mode_flip_here,
               chrom_calc_ok,
               rad_int_calc_ok,
               emit_6d_calc_ok,
               sigma_track_ok
           );
         },
         nb::arg("tao_lat") = nb::none(),
         nb::arg("spin") = nb::none(),
         nb::arg("srdt") = nb::none(),
         nb::arg("orb0") = nb::none(),
         nb::arg("modes_ri") = nb::none(),
         nb::arg("modes_6d") = nb::none(),
         nb::arg("ptc_normal_form") = nb::none(),
         nb::arg("bmad_normal_form") = nb::none(),
         nb::arg("cache_x_min") = nb::none(),
         nb::arg("cache_x_max") = nb::none(),
         nb::arg("comb_ds_save") = nb::none(),
         nb::arg("ix_ref_taylor") = nb::none(),
         nb::arg("ix_ele_taylor") = nb::none(),
         nb::arg("track_state") = nb::none(),
         nb::arg("cache_n_pts") = nb::none(),
         nb::arg("ix_rad_int_cache") = nb::none(),
         nb::arg("has_open_match_element") = nb::none(),
         nb::arg("plot_cache_valid") = nb::none(),
         nb::arg("spin_map_valid") = nb::none(),
         nb::arg("twiss_valid") = nb::none(),
         nb::arg("mode_flip_here") = nb::none(),
         nb::arg("chrom_calc_ok") = nb::none(),
         nb::arg("rad_int_calc_ok") = nb::none(),
         nb::arg("emit_6d_calc_ok") = nb::none(),
         nb::arg("sigma_track_ok") = nb::none()
  )
      .def_prop_rw(
          "tao_lat",
          &TaoLatticeBranchStruct::tao_lat,
          &TaoLatticeBranchStruct::set_tao_lat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Parent tao_lat"
      )
      .def_prop_ro(
          "lat_sigma",
          &TaoLatticeBranchStruct::lat_sigma,
          nb::keep_alive<0, 1>(),
          "Sigma matrix derived from lattice (not beam)."
      )
      .def_prop_ro(
          "spin_ele",
          &TaoLatticeBranchStruct::spin_ele,
          nb::keep_alive<0, 1>(),
          "Spin stuff"
      )
      .def_prop_ro(
          "bunch_params",
          &TaoLatticeBranchStruct::bunch_params,
          nb::keep_alive<0, 1>(),
          "Per element"
      )
      .def_prop_ro(
          "bunch_params_comb",
          &TaoLatticeBranchStruct::bunch_params_comb,
          nb::keep_alive<0, 1>(),
          "A comb for each bunch in beam."
      )
      .def_prop_ro("orbit", &TaoLatticeBranchStruct::orbit, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "plot_cache",
          &TaoLatticeBranchStruct::plot_cache,
          nb::keep_alive<0, 1>(),
          "Plotting data cache"
      )
      .def_prop_rw(
          "spin",
          &TaoLatticeBranchStruct::spin,
          &TaoLatticeBranchStruct::set_spin,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "srdt",
          &TaoLatticeBranchStruct::srdt,
          &TaoLatticeBranchStruct::set_srdt,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "orb0",
          &TaoLatticeBranchStruct::orb0,
          &TaoLatticeBranchStruct::set_orb0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For saving beginning orbit"
      )
      .def_prop_rw(
          "modes_ri",
          &TaoLatticeBranchStruct::modes_ri,
          &TaoLatticeBranchStruct::set_modes_ri,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Synchrotron integrals stuff"
      )
      .def_prop_rw(
          "modes_6d",
          &TaoLatticeBranchStruct::modes_6d,
          &TaoLatticeBranchStruct::set_modes_6d,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "6D radiation matrices."
      )
      .def_prop_rw(
          "ptc_normal_form",
          &TaoLatticeBranchStruct::ptc_normal_form,
          &TaoLatticeBranchStruct::set_ptc_normal_form,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Collection of normal form structures defined in PTC"
      )
      .def_prop_rw(
          "bmad_normal_form",
          &TaoLatticeBranchStruct::bmad_normal_form,
          &TaoLatticeBranchStruct::set_bmad_normal_form,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Collection of normal form structures defined in Bmad"
      )
      .def_prop_ro("high_E_orb", &TaoLatticeBranchStruct::high_E_orb, nb::keep_alive<0, 1>())
      .def_prop_ro("low_E_orb", &TaoLatticeBranchStruct::low_E_orb, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "taylor_save",
          &TaoLatticeBranchStruct::taylor_save,
          nb::keep_alive<0, 1>(),
          "Save to reduce computation time."
      )
      .def_prop_rw(
          "cache_x_min",
          &TaoLatticeBranchStruct::cache_x_min,
          &TaoLatticeBranchStruct::set_cache_x_min
      )
      .def_prop_rw(
          "cache_x_max",
          &TaoLatticeBranchStruct::cache_x_max,
          &TaoLatticeBranchStruct::set_cache_x_max
      )
      .def_prop_rw(
          "comb_ds_save",
          &TaoLatticeBranchStruct::comb_ds_save,
          &TaoLatticeBranchStruct::set_comb_ds_save,
          "Master parameter for %bunch_params_comb(:)%ds_save"
      )
      .def_prop_rw(
          "ix_ref_taylor",
          &TaoLatticeBranchStruct::ix_ref_taylor,
          &TaoLatticeBranchStruct::set_ix_ref_taylor
      )
      .def_prop_rw(
          "ix_ele_taylor",
          &TaoLatticeBranchStruct::ix_ele_taylor,
          &TaoLatticeBranchStruct::set_ix_ele_taylor
      )
      .def_prop_rw(
          "track_state",
          &TaoLatticeBranchStruct::track_state,
          &TaoLatticeBranchStruct::set_track_state
      )
      .def_prop_rw(
          "cache_n_pts",
          &TaoLatticeBranchStruct::cache_n_pts,
          &TaoLatticeBranchStruct::set_cache_n_pts
      )
      .def_prop_rw(
          "ix_rad_int_cache",
          &TaoLatticeBranchStruct::ix_rad_int_cache,
          &TaoLatticeBranchStruct::set_ix_rad_int_cache,
          "Radiation integrals cache index."
      )
      .def_prop_rw(
          "has_open_match_element",
          &TaoLatticeBranchStruct::has_open_match_element,
          &TaoLatticeBranchStruct::set_has_open_match_element
      )
      .def_prop_rw(
          "plot_cache_valid",
          &TaoLatticeBranchStruct::plot_cache_valid,
          &TaoLatticeBranchStruct::set_plot_cache_valid,
          "Valid plotting data cache?"
      )
      .def_prop_rw(
          "spin_map_valid",
          &TaoLatticeBranchStruct::spin_map_valid,
          &TaoLatticeBranchStruct::set_spin_map_valid
      )
      .def_prop_rw(
          "twiss_valid",
          &TaoLatticeBranchStruct::twiss_valid,
          &TaoLatticeBranchStruct::set_twiss_valid,
          "Invalid EG with unstable 1-turn matrix with a closed branch. With open branch: "
          "twiss_valid = T even if some Twiss (and orbit) is invalid."
      )
      .def_prop_rw(
          "mode_flip_here",
          &TaoLatticeBranchStruct::mode_flip_here,
          &TaoLatticeBranchStruct::set_mode_flip_here,
          "Twiss parameter mode flip seen?"
      )
      .def_prop_rw(
          "chrom_calc_ok",
          &TaoLatticeBranchStruct::chrom_calc_ok,
          &TaoLatticeBranchStruct::set_chrom_calc_ok
      )
      .def_prop_rw(
          "rad_int_calc_ok",
          &TaoLatticeBranchStruct::rad_int_calc_ok,
          &TaoLatticeBranchStruct::set_rad_int_calc_ok
      )
      .def_prop_rw(
          "emit_6d_calc_ok",
          &TaoLatticeBranchStruct::emit_6d_calc_ok,
          &TaoLatticeBranchStruct::set_emit_6d_calc_ok
      )
      .def_prop_rw(
          "sigma_track_ok",
          &TaoLatticeBranchStruct::sigma_track_ok,
          &TaoLatticeBranchStruct::set_sigma_track_ok
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoLatticeBranchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoLatticeBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoLatticeBranchStruct &self, nb::dict &memo) {
            return TaoLatticeBranchStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoLatticeBranchStruct &self, const TaoLatticeBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_lattice_struct(nb::module_ &m, nb::class_<TaoLatticeStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoLatticeStruct *self,
            std::optional<std::string> name,
            const LatStruct *lat,
            const LatStruct *high_E_lat,
            const LatStruct *low_E_lat,
            const TaoUniverseStruct *u,
            const RadIntAllEleStruct *rad_int_by_ele_ri,
            const RadIntAllEleStruct *rad_int_by_ele_6d) {
           new (self) TaoLatticeStruct(
               name,
               ptr_to_opt_ref(lat),
               ptr_to_opt_ref(high_E_lat),
               ptr_to_opt_ref(low_E_lat),
               ptr_to_opt_ref(u),
               ptr_to_opt_ref(rad_int_by_ele_ri),
               ptr_to_opt_ref(rad_int_by_ele_6d)
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("lat") = nb::none(),
         nb::arg("high_E_lat") = nb::none(),
         nb::arg("low_E_lat") = nb::none(),
         nb::arg("u") = nb::none(),
         nb::arg("rad_int_by_ele_ri") = nb::none(),
         nb::arg("rad_int_by_ele_6d") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoLatticeStruct::name,
          &TaoLatticeStruct::set_name,
          "'model', 'base', or 'design'."
      )
      .def_prop_rw(
          "lat",
          &TaoLatticeStruct::lat,
          &TaoLatticeStruct::set_lat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "lattice structures"
      )
      .def_prop_rw(
          "high_E_lat",
          &TaoLatticeStruct::high_E_lat,
          &TaoLatticeStruct::set_high_E_lat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For chrom calc."
      )
      .def_prop_rw(
          "low_E_lat",
          &TaoLatticeStruct::low_E_lat,
          &TaoLatticeStruct::set_low_E_lat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For chrom calc."
      )
      .def_prop_rw(
          "u",
          &TaoLatticeStruct::u,
          &TaoLatticeStruct::set_u,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Parent universe"
      )
      .def_prop_rw(
          "rad_int_by_ele_ri",
          &TaoLatticeStruct::rad_int_by_ele_ri,
          &TaoLatticeStruct::set_rad_int_by_ele_ri,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "rad_int_by_ele_6d",
          &TaoLatticeStruct::rad_int_by_ele_6d,
          &TaoLatticeStruct::set_rad_int_by_ele_6d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("tao_branch", &TaoLatticeStruct::tao_branch, nb::keep_alive<0, 1>())

      .def("__repr__", [](const TaoLatticeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLatticeStruct &self) {
            return TaoLatticeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoLatticeStruct &self, nb::dict &memo) { return TaoLatticeStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoLatticeStruct &self, const TaoLatticeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_logical_array_struct
void init_tao_logical_array_struct(nb::module_ &m, nb::class_<TaoLogicalArrayStruct> &cls) {
  cls.def(nb::init<std::optional<bool>>(), nb::arg("l") = nb::none())
      .def_prop_rw("l", &TaoLogicalArrayStruct::l, &TaoLogicalArrayStruct::set_l)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoLogicalArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoLogicalArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoLogicalArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoLogicalArrayStruct &self) {
            return TaoLogicalArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoLogicalArrayStruct &self, nb::dict &memo) {
            return TaoLogicalArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoLogicalArrayStruct &self, const TaoLogicalArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoLogicalArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoLogicalArrayStructArray1D, TaoLogicalArrayStructAlloc1D>(
      m,
      "TaoLogicalArrayStructArray1D",
      "TaoLogicalArrayStructAlloc1D"
  );
  // 2D TaoLogicalArrayStruct arrays are not used in structs/routines
  // 3D TaoLogicalArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_model_branch_struct
void init_tao_model_branch_struct(nb::module_ &m, nb::class_<TaoModelBranchStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoModelBranchStruct *self, const TaoBeamBranchStruct *beam) {
           new (self) TaoModelBranchStruct(ptr_to_opt_ref(beam));
         },
         nb::arg("beam") = nb::none()
  )
      .def_prop_ro(
          "ele",
          &TaoModelBranchStruct::ele,
          nb::keep_alive<0, 1>(),
          "Per element information"
      )
      .def_prop_rw(
          "beam",
          &TaoModelBranchStruct::beam,
          &TaoModelBranchStruct::set_beam,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoModelBranchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoModelBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoModelBranchStruct &self, nb::dict &memo) {
            return TaoModelBranchStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoModelBranchStruct &self, const TaoModelBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_model_element_struct(nb::module_ &m, nb::class_<TaoModelElementStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoModelElementStruct *self,
            const BeamStruct *beam,
            std::optional<bool> save_beam_internally,
            std::optional<bool> save_beam_to_file) {
           new (self)
               TaoModelElementStruct(ptr_to_opt_ref(beam), save_beam_internally, save_beam_to_file);
         },
         nb::arg("beam") = nb::none(),
         nb::arg("save_beam_internally") = nb::none(),
         nb::arg("save_beam_to_file") = nb::none()
  )
      .def_prop_rw(
          "beam",
          &TaoModelElementStruct::beam,
          &TaoModelElementStruct::set_beam,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Beam distribution at element."
      )
      .def_prop_rw(
          "save_beam_internally",
          &TaoModelElementStruct::save_beam_internally,
          &TaoModelElementStruct::set_save_beam_internally,
          "Save beam here? Beam also saved at fork elements and at track ends."
      )
      .def_prop_rw(
          "save_beam_to_file",
          &TaoModelElementStruct::save_beam_to_file,
          &TaoModelElementStruct::set_save_beam_to_file,
          "Save beam to a file? Beam also saved at fork elements and at track ends."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoModelElementStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoModelElementStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoModelElementStruct &self, nb::dict &memo) {
            return TaoModelElementStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoModelElementStruct &self, const TaoModelElementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_ping_scale_struct(nb::module_ &m, nb::class_<TaoPingScaleStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("a_mode_meas") = nb::none(),
         nb::arg("a_mode_ref") = nb::none(),
         nb::arg("b_mode_meas") = nb::none(),
         nb::arg("b_mode_ref") = nb::none()
  )
      .def_prop_rw(
          "a_mode_meas",
          &TaoPingScaleStruct::a_mode_meas,
          &TaoPingScaleStruct::set_a_mode_meas
      )
      .def_prop_rw(
          "a_mode_ref",
          &TaoPingScaleStruct::a_mode_ref,
          &TaoPingScaleStruct::set_a_mode_ref
      )
      .def_prop_rw(
          "b_mode_meas",
          &TaoPingScaleStruct::b_mode_meas,
          &TaoPingScaleStruct::set_b_mode_meas
      )
      .def_prop_rw(
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
          [](const TaoPingScaleStruct &self, nb::dict &memo) { return TaoPingScaleStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPingScaleStruct &self, const TaoPingScaleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_plot_array_struct
void init_tao_plot_array_struct(nb::module_ &m, nb::class_<TaoPlotArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotArrayStruct *self, const TaoPlotStruct *p) {
           new (self) TaoPlotArrayStruct(ptr_to_opt_ref(p));
         },
         nb::arg("p") = nb::none()
  )
      .def_prop_rw(
          "p",
          &TaoPlotArrayStruct::p,
          &TaoPlotArrayStruct::set_p,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoPlotArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoPlotArrayStruct &self) {
            return TaoPlotArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoPlotArrayStruct &self, nb::dict &memo) { return TaoPlotArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotArrayStruct &self, const TaoPlotArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoPlotArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoPlotArrayStructArray1D, TaoPlotArrayStructAlloc1D>(
      m,
      "TaoPlotArrayStructArray1D",
      "TaoPlotArrayStructAlloc1D"
  );
  // 2D TaoPlotArrayStruct arrays are not used in structs/routines
  // 3D TaoPlotArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_plot_cache_struct
void init_tao_plot_cache_struct(nb::module_ &m, nb::class_<TaoPlotCacheStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotCacheStruct *self,
            const EleStruct *ele_to_s,
            const CoordStruct *orbit,
            std::optional<bool> err) {
           new (self) TaoPlotCacheStruct(ptr_to_opt_ref(ele_to_s), ptr_to_opt_ref(orbit), err);
         },
         nb::arg("ele_to_s") = nb::none(),
         nb::arg("orbit") = nb::none(),
         nb::arg("err") = nb::none()
  )
      .def_prop_rw(
          "ele_to_s",
          &TaoPlotCacheStruct::ele_to_s,
          &TaoPlotCacheStruct::set_ele_to_s,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Integrated element from branch beginning. Will be marked as a hybrid element."
      )
      .def_prop_rw(
          "orbit",
          &TaoPlotCacheStruct::orbit,
          &TaoPlotCacheStruct::set_orbit,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("err", &TaoPlotCacheStruct::err, &TaoPlotCacheStruct::set_err)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotCacheStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotCacheStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoPlotCacheStruct &self, nb::dict &memo) { return TaoPlotCacheStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotCacheStruct &self, const TaoPlotCacheStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_plot_page_struct(nb::module_ &m, nb::class_<TaoPlotPageStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotPageStruct *self,
            const TaoTitleStruct *title,
            const TaoTitleStruct *subtitle,
            const QpRectStruct *border,
            const TaoDrawingStruct *floor_plan,
            const TaoDrawingStruct *lat_layout,
            std::optional<std::string> plot_display_type,
            std::optional<std::vector<double>> size,
            std::optional<double> text_height,
            std::optional<double> main_title_text_scale,
            std::optional<double> graph_title_text_scale,
            std::optional<double> axis_number_text_scale,
            std::optional<double> axis_label_text_scale,
            std::optional<double> legend_text_scale,
            std::optional<double> key_table_text_scale,
            std::optional<double> floor_plan_shape_scale,
            std::optional<double> floor_plan_text_scale,
            std::optional<double> lat_layout_shape_scale,
            std::optional<double> lat_layout_text_scale,
            std::optional<int> n_curve_pts,
            std::optional<int> id_window,
            std::optional<bool> delete_overlapping_plots,
            std::optional<bool> draw_graph_title_suffix) {
           new (self) TaoPlotPageStruct(
               ptr_to_opt_ref(title),
               ptr_to_opt_ref(subtitle),
               ptr_to_opt_ref(border),
               ptr_to_opt_ref(floor_plan),
               ptr_to_opt_ref(lat_layout),
               plot_display_type,
               size,
               text_height,
               main_title_text_scale,
               graph_title_text_scale,
               axis_number_text_scale,
               axis_label_text_scale,
               legend_text_scale,
               key_table_text_scale,
               floor_plan_shape_scale,
               floor_plan_text_scale,
               lat_layout_shape_scale,
               lat_layout_text_scale,
               n_curve_pts,
               id_window,
               delete_overlapping_plots,
               draw_graph_title_suffix
           );
         },
         nb::arg("title") = nb::none(),
         nb::arg("subtitle") = nb::none(),
         nb::arg("border") = nb::none(),
         nb::arg("floor_plan") = nb::none(),
         nb::arg("lat_layout") = nb::none(),
         nb::arg("plot_display_type") = nb::none(),
         nb::arg("size") = nb::none(),
         nb::arg("text_height") = nb::none(),
         nb::arg("main_title_text_scale") = nb::none(),
         nb::arg("graph_title_text_scale") = nb::none(),
         nb::arg("axis_number_text_scale") = nb::none(),
         nb::arg("axis_label_text_scale") = nb::none(),
         nb::arg("legend_text_scale") = nb::none(),
         nb::arg("key_table_text_scale") = nb::none(),
         nb::arg("floor_plan_shape_scale") = nb::none(),
         nb::arg("floor_plan_text_scale") = nb::none(),
         nb::arg("lat_layout_shape_scale") = nb::none(),
         nb::arg("lat_layout_text_scale") = nb::none(),
         nb::arg("n_curve_pts") = nb::none(),
         nb::arg("id_window") = nb::none(),
         nb::arg("delete_overlapping_plots") = nb::none(),
         nb::arg("draw_graph_title_suffix") = nb::none()
  )
      .def_prop_rw(
          "title",
          &TaoPlotPageStruct::title,
          &TaoPlotPageStruct::set_title,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Title  at top of page."
      )
      .def_prop_rw(
          "subtitle",
          &TaoPlotPageStruct::subtitle,
          &TaoPlotPageStruct::set_subtitle,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Subtitle below title at top of page."
      )
      .def_prop_rw(
          "border",
          &TaoPlotPageStruct::border,
          &TaoPlotPageStruct::set_border,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Border around plots edge of page."
      )
      .def_prop_rw(
          "floor_plan",
          &TaoPlotPageStruct::floor_plan,
          &TaoPlotPageStruct::set_floor_plan,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "lat_layout",
          &TaoPlotPageStruct::lat_layout,
          &TaoPlotPageStruct::set_lat_layout,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("pattern", &TaoPlotPageStruct::pattern, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "template_",
          &TaoPlotPageStruct::template_,
          nb::keep_alive<0, 1>(),
          "Templates for the plots."
      )
      .def_prop_ro("region", &TaoPlotPageStruct::region, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "plot_display_type",
          &TaoPlotPageStruct::plot_display_type,
          &TaoPlotPageStruct::set_plot_display_type,
          "'X' or 'TK'"
      )
      .def_prop_rw(
          "size",
          &TaoPlotPageStruct::size,
          &TaoPlotPageStruct::set_size,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "width and height of plot window in pixels."
      )
      .def_prop_rw(
          "text_height",
          &TaoPlotPageStruct::text_height,
          &TaoPlotPageStruct::set_text_height,
          "In points. Scales the height of all text"
      )
      .def_prop_rw(
          "main_title_text_scale",
          &TaoPlotPageStruct::main_title_text_scale,
          &TaoPlotPageStruct::set_main_title_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "graph_title_text_scale",
          &TaoPlotPageStruct::graph_title_text_scale,
          &TaoPlotPageStruct::set_graph_title_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "axis_number_text_scale",
          &TaoPlotPageStruct::axis_number_text_scale,
          &TaoPlotPageStruct::set_axis_number_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "axis_label_text_scale",
          &TaoPlotPageStruct::axis_label_text_scale,
          &TaoPlotPageStruct::set_axis_label_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "legend_text_scale",
          &TaoPlotPageStruct::legend_text_scale,
          &TaoPlotPageStruct::set_legend_text_scale,
          "Relative to text_height. For legends, plot_page, and lat_layout"
      )
      .def_prop_rw(
          "key_table_text_scale",
          &TaoPlotPageStruct::key_table_text_scale,
          &TaoPlotPageStruct::set_key_table_text_scale,
          "Relative to text_height"
      )
      .def_prop_rw(
          "floor_plan_shape_scale",
          &TaoPlotPageStruct::floor_plan_shape_scale,
          &TaoPlotPageStruct::set_floor_plan_shape_scale
      )
      .def_prop_rw(
          "floor_plan_text_scale",
          &TaoPlotPageStruct::floor_plan_text_scale,
          &TaoPlotPageStruct::set_floor_plan_text_scale,
          "Scale used = floor_plan_text_scale * legend_text_scale"
      )
      .def_prop_rw(
          "lat_layout_shape_scale",
          &TaoPlotPageStruct::lat_layout_shape_scale,
          &TaoPlotPageStruct::set_lat_layout_shape_scale
      )
      .def_prop_rw(
          "lat_layout_text_scale",
          &TaoPlotPageStruct::lat_layout_text_scale,
          &TaoPlotPageStruct::set_lat_layout_text_scale,
          "Scale used = lat_layout_text_scale * legend_text_scale"
      )
      .def_prop_rw(
          "n_curve_pts",
          &TaoPlotPageStruct::n_curve_pts,
          &TaoPlotPageStruct::set_n_curve_pts,
          "Default number of points for plotting a smooth curve."
      )
      .def_prop_rw(
          "id_window",
          &TaoPlotPageStruct::id_window,
          &TaoPlotPageStruct::set_id_window,
          "X window id number."
      )
      .def_prop_rw(
          "delete_overlapping_plots",
          &TaoPlotPageStruct::delete_overlapping_plots,
          &TaoPlotPageStruct::set_delete_overlapping_plots,
          "Delete overlapping plots when a plot is placed?"
      )
      .def_prop_rw(
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
          [](const TaoPlotPageStruct &self, nb::dict &memo) { return TaoPlotPageStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotPageStruct &self, const TaoPlotPageStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_plot_region_struct(nb::module_ &m, nb::class_<TaoPlotRegionStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotRegionStruct *self,
            std::optional<std::string> name,
            const TaoPlotStruct *plot,
            std::optional<std::vector<double>> location,
            std::optional<bool> visible,
            std::optional<bool> list_with_show_plot_command,
            std::optional<bool> setup_done) {
           new (self) TaoPlotRegionStruct(
               name,
               ptr_to_opt_ref(plot),
               location,
               visible,
               list_with_show_plot_command,
               setup_done
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("plot") = nb::none(),
         nb::arg("location") = nb::none(),
         nb::arg("visible") = nb::none(),
         nb::arg("list_with_show_plot_command") = nb::none(),
         nb::arg("setup_done") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoPlotRegionStruct::name,
          &TaoPlotRegionStruct::set_name,
          "Region name. Eg: 'r13', etc."
      )
      .def_prop_rw(
          "plot",
          &TaoPlotRegionStruct::plot,
          &TaoPlotRegionStruct::set_plot,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Plot associated with this region"
      )
      .def_prop_rw(
          "location",
          &TaoPlotRegionStruct::location,
          &TaoPlotRegionStruct::set_location,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "[x1, x2, y1, y2] location on page."
      )
      .def_prop_rw(
          "visible",
          &TaoPlotRegionStruct::visible,
          &TaoPlotRegionStruct::set_visible,
          "To draw or not to draw."
      )
      .def_prop_rw(
          "list_with_show_plot_command",
          &TaoPlotRegionStruct::list_with_show_plot_command,
          &TaoPlotRegionStruct::set_list_with_show_plot_command,
          "False used for default plots to shorten the output of 'show plot'"
      )
      .def_prop_rw(
          "setup_done",
          &TaoPlotRegionStruct::setup_done,
          &TaoPlotRegionStruct::set_setup_done,
          "Used for plot bookkeeping."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotRegionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotRegionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoPlotRegionStruct &self, nb::dict &memo) { return TaoPlotRegionStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotRegionStruct &self, const TaoPlotRegionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_plot_struct(nb::module_ &m, nb::class_<TaoPlotStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoPlotStruct *self,
            std::optional<std::string> name,
            std::optional<std::string> description,
            const TaoPlotRegionStruct *r,
            std::optional<int> ix_plot,
            std::optional<int> n_curve_pts,
            std::optional<std::string> type,
            std::optional<std::string> x_axis_type,
            std::optional<bool> autoscale_x,
            std::optional<bool> autoscale_y,
            std::optional<bool> autoscale_gang_x,
            std::optional<bool> autoscale_gang_y,
            std::optional<bool> list_with_show_plot_command,
            std::optional<bool> phantom,
            std::optional<bool> default_plot) {
           new (self) TaoPlotStruct(
               name,
               description,
               ptr_to_opt_ref(r),
               ix_plot,
               n_curve_pts,
               type,
               x_axis_type,
               autoscale_x,
               autoscale_y,
               autoscale_gang_x,
               autoscale_gang_y,
               list_with_show_plot_command,
               phantom,
               default_plot
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("description") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("ix_plot") = nb::none(),
         nb::arg("n_curve_pts") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("x_axis_type") = nb::none(),
         nb::arg("autoscale_x") = nb::none(),
         nb::arg("autoscale_y") = nb::none(),
         nb::arg("autoscale_gang_x") = nb::none(),
         nb::arg("autoscale_gang_y") = nb::none(),
         nb::arg("list_with_show_plot_command") = nb::none(),
         nb::arg("phantom") = nb::none(),
         nb::arg("default_plot") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoPlotStruct::name,
          &TaoPlotStruct::set_name,
          "Identifying name. Rule: If name is blank, plot is not valid."
      )
      .def_prop_rw(
          "description",
          &TaoPlotStruct::description,
          &TaoPlotStruct::set_description,
          "Descriptive string."
      )
      .def_prop_ro(
          "graph",
          &TaoPlotStruct::graph,
          nb::keep_alive<0, 1>(),
          "individual graphs of a plot"
      )
      .def_prop_rw(
          "r",
          &TaoPlotStruct::r,
          &TaoPlotStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "pointer to parent."
      )
      .def_prop_rw(
          "ix_plot",
          &TaoPlotStruct::ix_plot,
          &TaoPlotStruct::set_ix_plot,
          "Index in s%plot_page%template(:) or %region(:) arrays."
      )
      .def_prop_rw(
          "n_curve_pts",
          &TaoPlotStruct::n_curve_pts,
          &TaoPlotStruct::set_n_curve_pts,
          "Overrides s%plot_page%n_curve_pts."
      )
      .def_prop_rw("type", &TaoPlotStruct::type, &TaoPlotStruct::set_type, "or 'wave'")
      .def_prop_rw(
          "x_axis_type",
          &TaoPlotStruct::x_axis_type,
          &TaoPlotStruct::set_x_axis_type,
          "'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc."
      )
      .def_prop_rw(
          "autoscale_x",
          &TaoPlotStruct::autoscale_x,
          &TaoPlotStruct::set_autoscale_x,
          "Horizontal autoscale."
      )
      .def_prop_rw(
          "autoscale_y",
          &TaoPlotStruct::autoscale_y,
          &TaoPlotStruct::set_autoscale_y,
          "Vertical autoscale."
      )
      .def_prop_rw(
          "autoscale_gang_x",
          &TaoPlotStruct::autoscale_gang_x,
          &TaoPlotStruct::set_autoscale_gang_x,
          "scale cmd scales graphs together?"
      )
      .def_prop_rw(
          "autoscale_gang_y",
          &TaoPlotStruct::autoscale_gang_y,
          &TaoPlotStruct::set_autoscale_gang_y,
          "scale cmd scales graphs together?"
      )
      .def_prop_rw(
          "list_with_show_plot_command",
          &TaoPlotStruct::list_with_show_plot_command,
          &TaoPlotStruct::set_list_with_show_plot_command,
          "False used for default plots to shorten the output of 'show plot'"
      )
      .def_prop_rw(
          "phantom",
          &TaoPlotStruct::phantom,
          &TaoPlotStruct::set_phantom,
          "Used by tao_plot_init to add info lines to 'show plot -templates'"
      )
      .def_prop_rw(
          "default_plot",
          &TaoPlotStruct::default_plot,
          &TaoPlotStruct::set_default_plot,
          "One of Tao's default plots?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoPlotStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoPlotStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoPlotStruct &self, nb::dict &memo) { return TaoPlotStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoPlotStruct &self, const TaoPlotStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_real_pointer_struct
void init_tao_real_pointer_struct(nb::module_ &m, nb::class_<TaoRealPointerStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<bool>, std::optional<bool>>(),
         nb::arg("r") = nb::none(),
         nb::arg("good_value") = nb::none(),
         nb::arg("good_user") = nb::none()
  )
      .def_prop_rw("r", &TaoRealPointerStruct::r, &TaoRealPointerStruct::set_r)
      .def_prop_rw(
          "good_value",
          &TaoRealPointerStruct::good_value,
          &TaoRealPointerStruct::set_good_value
      )
      .def_prop_rw(
          "good_user",
          &TaoRealPointerStruct::good_user,
          &TaoRealPointerStruct::set_good_user
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoRealPointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoRealPointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoRealPointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoRealPointerStruct &self) {
            return TaoRealPointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoRealPointerStruct &self, nb::dict &memo) {
            return TaoRealPointerStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoRealPointerStruct &self, const TaoRealPointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoRealPointerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoRealPointerStructArray1D, TaoRealPointerStructAlloc1D>(
      m,
      "TaoRealPointerStructArray1D",
      "TaoRealPointerStructAlloc1D"
  );
  // 2D TaoRealPointerStruct arrays are not used in structs/routines
  // 3D TaoRealPointerStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_shape_pattern_point_struct
void init_tao_shape_pattern_point_struct(
    nb::module_ &m,
    nb::class_<TaoShapePatternPointStruct> &cls
) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("s") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("radius") = nb::none()
  )
      .def_prop_rw("s", &TaoShapePatternPointStruct::s, &TaoShapePatternPointStruct::set_s)
      .def_prop_rw("y", &TaoShapePatternPointStruct::y, &TaoShapePatternPointStruct::set_y)
      .def_prop_rw(
          "radius",
          &TaoShapePatternPointStruct::radius,
          &TaoShapePatternPointStruct::set_radius
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoShapePatternPointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoShapePatternPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoShapePatternPointStruct &self, nb::dict &memo) {
            return TaoShapePatternPointStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoShapePatternPointStruct &self, const TaoShapePatternPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_shape_pattern_struct(nb::module_ &m, nb::class_<TaoShapePatternStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoShapePatternStruct *self,
            std::optional<std::string> name,
            const QpLineStruct *line) {
           new (self) TaoShapePatternStruct(name, ptr_to_opt_ref(line));
         },
         nb::arg("name") = nb::none(),
         nb::arg("line") = nb::none()
  )
      .def_prop_rw("name", &TaoShapePatternStruct::name, &TaoShapePatternStruct::set_name)
      .def_prop_rw(
          "line",
          &TaoShapePatternStruct::line,
          &TaoShapePatternStruct::set_line,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Line color and pattern set by shape using this pattern."
      )
      .def_prop_ro("pt", &TaoShapePatternStruct::pt, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return TaoShapePatternStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoShapePatternStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoShapePatternStruct &self, nb::dict &memo) {
            return TaoShapePatternStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoShapePatternStruct &self, const TaoShapePatternStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_spin_dn_dpz_struct(nb::module_ &m, nb::class_<TaoSpinDnDpzStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>>(),
         nb::arg("vec") = nb::none(),
         nb::arg("partial") = nb::none(),
         nb::arg("partial2") = nb::none()
  )
      .def_prop_rw(
          "vec",
          &TaoSpinDnDpzStruct::vec,
          &TaoSpinDnDpzStruct::set_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "n0 derivative wrt pz."
      )
      .def_prop_rw(
          "partial",
          &TaoSpinDnDpzStruct::partial,
          &TaoSpinDnDpzStruct::set_partial,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "partial(i:) is spin n0 derivative wrt pz for i^th oscillation mode (1 => a-mode, etc.)"
      )
      .def_prop_rw(
          "partial2",
          &TaoSpinDnDpzStruct::partial2,
          &TaoSpinDnDpzStruct::set_partial2,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const TaoSpinDnDpzStruct &self, nb::dict &memo) { return TaoSpinDnDpzStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoSpinDnDpzStruct &self, const TaoSpinDnDpzStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_spin_ele_struct(nb::module_ &m, nb::class_<TaoSpinEleStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoSpinEleStruct *self,
            const TaoSpinDnDpzStruct *dn_dpz,
            std::optional<std::vector<double>> orb_eigen_val,
            std::optional<std::vector<std::vector<double>>> orb_eigen_vec,
            std::optional<std::vector<std::vector<double>>> spin_eigen_vec,
            std::optional<bool> valid) {
           new (self) TaoSpinEleStruct(
               ptr_to_opt_ref(dn_dpz),
               orb_eigen_val,
               orb_eigen_vec,
               spin_eigen_vec,
               valid
           );
         },
         nb::arg("dn_dpz") = nb::none(),
         nb::arg("orb_eigen_val") = nb::none(),
         nb::arg("orb_eigen_vec") = nb::none(),
         nb::arg("spin_eigen_vec") = nb::none(),
         nb::arg("valid") = nb::none()
  )
      .def_prop_rw(
          "dn_dpz",
          &TaoSpinEleStruct::dn_dpz,
          &TaoSpinEleStruct::set_dn_dpz,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "orb_eigen_val",
          &TaoSpinEleStruct::orb_eigen_val,
          &TaoSpinEleStruct::set_orb_eigen_val,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "orb_eigen_vec",
          &TaoSpinEleStruct::orb_eigen_vec,
          &TaoSpinEleStruct::set_orb_eigen_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(j,:) is j^th vector"
      )
      .def_prop_rw(
          "spin_eigen_vec",
          &TaoSpinEleStruct::spin_eigen_vec,
          &TaoSpinEleStruct::set_spin_eigen_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(j,:) is j^th vector"
      )
      .def_prop_rw("valid", &TaoSpinEleStruct::valid, &TaoSpinEleStruct::set_valid)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoSpinEleStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoSpinEleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoSpinEleStruct &self, nb::dict &memo) { return TaoSpinEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoSpinEleStruct &self, const TaoSpinEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_spin_map_struct(nb::module_ &m, nb::class_<TaoSpinMapStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoSpinMapStruct *self,
            std::optional<bool> valid,
            const SpinOrbitMap1Struct *map1,
            const SpinAxisStruct *axis_input,
            const SpinAxisStruct *axis0,
            const SpinAxisStruct *axis1,
            std::optional<int> ix_ele,
            std::optional<int> ix_ref,
            std::optional<int> ix_uni,
            std::optional<int> ix_branch,
            std::optional<std::vector<std::vector<double>>> mat8) {
           new (self) TaoSpinMapStruct(
               valid,
               ptr_to_opt_ref(map1),
               ptr_to_opt_ref(axis_input),
               ptr_to_opt_ref(axis0),
               ptr_to_opt_ref(axis1),
               ix_ele,
               ix_ref,
               ix_uni,
               ix_branch,
               mat8
           );
         },
         nb::arg("valid") = nb::none(),
         nb::arg("map1") = nb::none(),
         nb::arg("axis_input") = nb::none(),
         nb::arg("axis0") = nb::none(),
         nb::arg("axis1") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_ref") = nb::none(),
         nb::arg("ix_uni") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("mat8") = nb::none()
  )
      .def_prop_rw("valid", &TaoSpinMapStruct::valid, &TaoSpinMapStruct::set_valid)
      .def_prop_rw(
          "map1",
          &TaoSpinMapStruct::map1,
          &TaoSpinMapStruct::set_map1,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "axis_input",
          &TaoSpinMapStruct::axis_input,
          &TaoSpinMapStruct::set_axis_input,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Input axes."
      )
      .def_prop_rw(
          "axis0",
          &TaoSpinMapStruct::axis0,
          &TaoSpinMapStruct::set_axis0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Initial axes."
      )
      .def_prop_rw(
          "axis1",
          &TaoSpinMapStruct::axis1,
          &TaoSpinMapStruct::set_axis1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Final axes."
      )
      .def_prop_rw("ix_ele", &TaoSpinMapStruct::ix_ele, &TaoSpinMapStruct::set_ix_ele)
      .def_prop_rw("ix_ref", &TaoSpinMapStruct::ix_ref, &TaoSpinMapStruct::set_ix_ref)
      .def_prop_rw("ix_uni", &TaoSpinMapStruct::ix_uni, &TaoSpinMapStruct::set_ix_uni)
      .def_prop_rw("ix_branch", &TaoSpinMapStruct::ix_branch, &TaoSpinMapStruct::set_ix_branch)
      .def_prop_rw(
          "mat8",
          &TaoSpinMapStruct::mat8,
          &TaoSpinMapStruct::set_mat8,
          nb::for_getter(nb::keep_alive<0, 1>())
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
          [](const TaoSpinMapStruct &self, nb::dict &memo) { return TaoSpinMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoSpinMapStruct &self, const TaoSpinMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_spin_polarization_struct(nb::module_ &m, nb::class_<TaoSpinPolarizationStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoSpinPolarizationStruct *self,
            std::optional<double> tune,
            std::optional<double> pol_limit_st,
            std::optional<double> pol_limit_dk,
            std::optional<std::vector<double>> pol_limit_dk_partial,
            std::optional<std::vector<double>> pol_limit_dk_partial2,
            std::optional<double> pol_rate_bks,
            std::optional<double> depol_rate,
            std::optional<std::vector<double>> depol_rate_partial,
            std::optional<std::vector<double>> depol_rate_partial2,
            std::optional<double> integral_bn,
            std::optional<double> integral_bdn,
            std::optional<double> integral_1ns,
            std::optional<double> integral_dn2,
            std::optional<bool> valid,
            const SpinOrbitMap1Struct *q_1turn) {
           new (self) TaoSpinPolarizationStruct(
               tune,
               pol_limit_st,
               pol_limit_dk,
               pol_limit_dk_partial,
               pol_limit_dk_partial2,
               pol_rate_bks,
               depol_rate,
               depol_rate_partial,
               depol_rate_partial2,
               integral_bn,
               integral_bdn,
               integral_1ns,
               integral_dn2,
               valid,
               ptr_to_opt_ref(q_1turn)
           );
         },
         nb::arg("tune") = nb::none(),
         nb::arg("pol_limit_st") = nb::none(),
         nb::arg("pol_limit_dk") = nb::none(),
         nb::arg("pol_limit_dk_partial") = nb::none(),
         nb::arg("pol_limit_dk_partial2") = nb::none(),
         nb::arg("pol_rate_bks") = nb::none(),
         nb::arg("depol_rate") = nb::none(),
         nb::arg("depol_rate_partial") = nb::none(),
         nb::arg("depol_rate_partial2") = nb::none(),
         nb::arg("integral_bn") = nb::none(),
         nb::arg("integral_bdn") = nb::none(),
         nb::arg("integral_1ns") = nb::none(),
         nb::arg("integral_dn2") = nb::none(),
         nb::arg("valid") = nb::none(),
         nb::arg("q_1turn") = nb::none()
  )
      .def_prop_rw("tune", &TaoSpinPolarizationStruct::tune, &TaoSpinPolarizationStruct::set_tune)
      .def_prop_rw(
          "pol_limit_st",
          &TaoSpinPolarizationStruct::pol_limit_st,
          &TaoSpinPolarizationStruct::set_pol_limit_st,
          "Polarization calculated using Sokolov-Ternov formula."
      )
      .def_prop_rw(
          "pol_limit_dk",
          &TaoSpinPolarizationStruct::pol_limit_dk,
          &TaoSpinPolarizationStruct::set_pol_limit_dk,
          "Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula."
      )
      .def_prop_rw(
          "pol_limit_dk_partial",
          &TaoSpinPolarizationStruct::pol_limit_dk_partial,
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Limit using only single mode to calc dn_dpz"
      )
      .def_prop_rw(
          "pol_limit_dk_partial2",
          &TaoSpinPolarizationStruct::pol_limit_dk_partial2,
          &TaoSpinPolarizationStruct::set_pol_limit_dk_partial2,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Limit using only single mode to calc dn_dpz"
      )
      .def_prop_rw(
          "pol_rate_bks",
          &TaoSpinPolarizationStruct::pol_rate_bks,
          &TaoSpinPolarizationStruct::set_pol_rate_bks,
          "BKS Polarization rate (1/sec)."
      )
      .def_prop_rw(
          "depol_rate",
          &TaoSpinPolarizationStruct::depol_rate,
          &TaoSpinPolarizationStruct::set_depol_rate,
          "Depolarization rate (1/sec)."
      )
      .def_prop_rw(
          "depol_rate_partial",
          &TaoSpinPolarizationStruct::depol_rate_partial,
          &TaoSpinPolarizationStruct::set_depol_rate_partial,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Depolarization rate (1/sec) using only single mode to calc dn_dpz."
      )
      .def_prop_rw(
          "depol_rate_partial2",
          &TaoSpinPolarizationStruct::depol_rate_partial2,
          &TaoSpinPolarizationStruct::set_depol_rate_partial2,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Depolarization rate (1/sec) using only two modes to calc dn_dpz."
      )
      .def_prop_rw(
          "integral_bn",
          &TaoSpinPolarizationStruct::integral_bn,
          &TaoSpinPolarizationStruct::set_integral_bn,
          "Integral of g^3 * b_hat * n_0"
      )
      .def_prop_rw(
          "integral_bdn",
          &TaoSpinPolarizationStruct::integral_bdn,
          &TaoSpinPolarizationStruct::set_integral_bdn,
          "Integral of g^3 * b_hat * dn/ddelta"
      )
      .def_prop_rw(
          "integral_1ns",
          &TaoSpinPolarizationStruct::integral_1ns,
          &TaoSpinPolarizationStruct::set_integral_1ns,
          "Integral of g^3 (1 - 2(n * s_hat)/9)"
      )
      .def_prop_rw(
          "integral_dn2",
          &TaoSpinPolarizationStruct::integral_dn2,
          &TaoSpinPolarizationStruct::set_integral_dn2,
          "Integral of g^3 * 11 (dn/ddelta)^2 / 9"
      )
      .def_prop_rw(
          "valid",
          &TaoSpinPolarizationStruct::valid,
          &TaoSpinPolarizationStruct::set_valid
      )
      .def_prop_rw(
          "q_1turn",
          &TaoSpinPolarizationStruct::q_1turn,
          &TaoSpinPolarizationStruct::set_q_1turn,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Save results from spin_concat_linear_maps in tao_spin_polarization."
      )
      .def_prop_ro(
          "q_ele",
          &TaoSpinPolarizationStruct::q_ele,
          nb::keep_alive<0, 1>(),
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
          [](const TaoSpinPolarizationStruct &self, nb::dict &memo) {
            return TaoSpinPolarizationStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoSpinPolarizationStruct &self, const TaoSpinPolarizationStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_string_array_struct
void init_tao_string_array_struct(nb::module_ &m, nb::class_<TaoStringArrayStruct> &cls) {
  cls.def(nb::init<std::optional<std::string>>(), nb::arg("s") = nb::none())
      .def_prop_rw("s", &TaoStringArrayStruct::s, &TaoStringArrayStruct::set_s)
      .def_static(
          "new_array1d",
          [](int sz) { return TaoStringArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoStringArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoStringArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoStringArrayStruct &self) {
            return TaoStringArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoStringArrayStruct &self, nb::dict &memo) {
            return TaoStringArrayStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoStringArrayStruct &self, const TaoStringArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoStringArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoStringArrayStructArray1D, TaoStringArrayStructAlloc1D>(
      m,
      "TaoStringArrayStructArray1D",
      "TaoStringArrayStructAlloc1D"
  );
  // 2D TaoStringArrayStruct arrays are not used in structs/routines
  // 3D TaoStringArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_super_universe_struct
void init_tao_super_universe_struct(nb::module_ &m, nb::class_<TaoSuperUniverseStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoSuperUniverseStruct *self,
            const TaoGlobalStruct *global,
            const TaoInitStruct *init,
            const TaoCommonStruct *com,
            const TaoPlotPageStruct *plot_page,
            std::optional<std::vector<int>> key,
            const TaoBuildingWallStruct *building_wall,
            const TaoWaveStruct *wave,
            std::optional<int> n_var_used,
            std::optional<int> n_v1_var_used,
            std::optional<bool> initialized) {
           new (self) TaoSuperUniverseStruct(
               ptr_to_opt_ref(global),
               ptr_to_opt_ref(init),
               ptr_to_opt_ref(com),
               ptr_to_opt_ref(plot_page),
               key,
               ptr_to_opt_ref(building_wall),
               ptr_to_opt_ref(wave),
               n_var_used,
               n_v1_var_used,
               initialized
           );
         },
         nb::arg("global_") = nb::none(),
         nb::arg("init") = nb::none(),
         nb::arg("com") = nb::none(),
         nb::arg("plot_page") = nb::none(),
         nb::arg("key") = nb::none(),
         nb::arg("building_wall") = nb::none(),
         nb::arg("wave") = nb::none(),
         nb::arg("n_var_used") = nb::none(),
         nb::arg("n_v1_var_used") = nb::none(),
         nb::arg("initialized") = nb::none()
  )
      .def_prop_rw(
          "global_",
          &TaoSuperUniverseStruct::global,
          &TaoSuperUniverseStruct::set_global,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "User accessible global variables."
      )
      .def_prop_rw(
          "init",
          &TaoSuperUniverseStruct::init,
          &TaoSuperUniverseStruct::set_init,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Initialization parameters"
      )
      .def_prop_rw(
          "com",
          &TaoSuperUniverseStruct::com,
          &TaoSuperUniverseStruct::set_com,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Non-initialization common parameters"
      )
      .def_prop_rw(
          "plot_page",
          &TaoSuperUniverseStruct::plot_page,
          &TaoSuperUniverseStruct::set_plot_page,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Defines the plot window."
      )
      .def_prop_ro(
          "v1_var",
          &TaoSuperUniverseStruct::v1_var,
          nb::keep_alive<0, 1>(),
          "The variable types"
      )
      .def_prop_ro(
          "var",
          &TaoSuperUniverseStruct::var,
          nb::keep_alive<0, 1>(),
          "array of all variables."
      )
      .def_prop_ro("u", &TaoSuperUniverseStruct::u, nb::keep_alive<0, 1>(), "array of universes.")
      .def_prop_rw(
          "key",
          &TaoSuperUniverseStruct::key,
          &TaoSuperUniverseStruct::set_key,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "building_wall",
          &TaoSuperUniverseStruct::building_wall,
          &TaoSuperUniverseStruct::set_building_wall,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "wave",
          &TaoSuperUniverseStruct::wave,
          &TaoSuperUniverseStruct::set_wave,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "n_var_used",
          &TaoSuperUniverseStruct::n_var_used,
          &TaoSuperUniverseStruct::set_n_var_used
      )
      .def_prop_rw(
          "n_v1_var_used",
          &TaoSuperUniverseStruct::n_v1_var_used,
          &TaoSuperUniverseStruct::set_n_v1_var_used
      )
      .def_prop_ro(
          "history",
          &TaoSuperUniverseStruct::history,
          nb::keep_alive<0, 1>(),
          "command history"
      )
      .def_prop_rw(
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
          [](const TaoSuperUniverseStruct &self, nb::dict &memo) {
            return TaoSuperUniverseStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoSuperUniverseStruct &self, const TaoSuperUniverseStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_title_struct(nb::module_ &m, nb::class_<TaoTitleStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>>(),
         nb::arg("string") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("units") = nb::none(),
         nb::arg("justify") = nb::none(),
         nb::arg("draw_it") = nb::none()
  )
      .def_prop_rw(
          "string",
          &TaoTitleStruct::string,
          &TaoTitleStruct::set_string,
          "title character string."
      )
      .def_prop_rw("x", &TaoTitleStruct::x, &TaoTitleStruct::set_x, "x, y rwt lower left corner")
      .def_prop_rw("y", &TaoTitleStruct::y, &TaoTitleStruct::set_y, "x, y rwt lower left corner")
      .def_prop_rw(
          "units",
          &TaoTitleStruct::units,
          &TaoTitleStruct::set_units,
          "%BOX, POINTS, etc..."
      )
      .def_prop_rw(
          "justify",
          &TaoTitleStruct::justify,
          &TaoTitleStruct::set_justify,
          "Left, Center, or Right justification."
      )
      .def_prop_rw(
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
          [](const TaoTitleStruct &self, nb::dict &memo) { return TaoTitleStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoTitleStruct &self, const TaoTitleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_universe_calc_struct(nb::module_ &m, nb::class_<TaoUniverseCalcStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("srdt_for_data") = nb::none(),
         nb::arg("rad_int_for_data") = nb::none(),
         nb::arg("rad_int_for_plotting") = nb::none(),
         nb::arg("chrom_for_data") = nb::none(),
         nb::arg("chrom_for_plotting") = nb::none(),
         nb::arg("lat_sigma_for_data") = nb::none(),
         nb::arg("lat_sigma_for_plotting") = nb::none(),
         nb::arg("dynamic_aperture") = nb::none(),
         nb::arg("one_turn_map") = nb::none(),
         nb::arg("lattice") = nb::none(),
         nb::arg("twiss") = nb::none(),
         nb::arg("track") = nb::none(),
         nb::arg("spin_matrices") = nb::none()
  )
      .def_prop_rw(
          "srdt_for_data",
          &TaoUniverseCalcStruct::srdt_for_data,
          &TaoUniverseCalcStruct::set_srdt_for_data,
          "0 = false, 1 = 1st order, 2 = 1st & 2nd order"
      )
      .def_prop_rw(
          "rad_int_for_data",
          &TaoUniverseCalcStruct::rad_int_for_data,
          &TaoUniverseCalcStruct::set_rad_int_for_data,
          "Do the radiation integrals need to be computed for"
      )
      .def_prop_rw(
          "rad_int_for_plotting",
          &TaoUniverseCalcStruct::rad_int_for_plotting,
          &TaoUniverseCalcStruct::set_rad_int_for_plotting,
          "data or plotting?"
      )
      .def_prop_rw(
          "chrom_for_data",
          &TaoUniverseCalcStruct::chrom_for_data,
          &TaoUniverseCalcStruct::set_chrom_for_data,
          "Does the chromaticity need to be computed for"
      )
      .def_prop_rw(
          "chrom_for_plotting",
          &TaoUniverseCalcStruct::chrom_for_plotting,
          &TaoUniverseCalcStruct::set_chrom_for_plotting,
          "data or plotting?"
      )
      .def_prop_rw(
          "lat_sigma_for_data",
          &TaoUniverseCalcStruct::lat_sigma_for_data,
          &TaoUniverseCalcStruct::set_lat_sigma_for_data,
          "Do the beam sigmas need to be computed for"
      )
      .def_prop_rw(
          "lat_sigma_for_plotting",
          &TaoUniverseCalcStruct::lat_sigma_for_plotting,
          &TaoUniverseCalcStruct::set_lat_sigma_for_plotting,
          "data or plotting?"
      )
      .def_prop_rw(
          "dynamic_aperture",
          &TaoUniverseCalcStruct::dynamic_aperture,
          &TaoUniverseCalcStruct::set_dynamic_aperture,
          "Do the dynamic_aperture calc?"
      )
      .def_prop_rw(
          "one_turn_map",
          &TaoUniverseCalcStruct::one_turn_map,
          &TaoUniverseCalcStruct::set_one_turn_map,
          "Compute the one turn map?"
      )
      .def_prop_rw(
          "lattice",
          &TaoUniverseCalcStruct::lattice,
          &TaoUniverseCalcStruct::set_lattice,
          "Used to indicate which lattices need tracking done."
      )
      .def_prop_rw(
          "twiss",
          &TaoUniverseCalcStruct::twiss,
          &TaoUniverseCalcStruct::set_twiss,
          "calc linear transfer matrix?"
      )
      .def_prop_rw(
          "track",
          &TaoUniverseCalcStruct::track,
          &TaoUniverseCalcStruct::set_track,
          "tracking needs to be done?"
      )
      .def_prop_rw(
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
          [](const TaoUniverseCalcStruct &self, nb::dict &memo) {
            return TaoUniverseCalcStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoUniverseCalcStruct &self, const TaoUniverseCalcStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_universe_pointer_struct(nb::module_ &m, nb::class_<TaoUniversePointerStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoUniversePointerStruct *self, const TaoUniverseStruct *u) {
           new (self) TaoUniversePointerStruct(ptr_to_opt_ref(u));
         },
         nb::arg("u") = nb::none()
  )
      .def_prop_rw(
          "u",
          &TaoUniversePointerStruct::u,
          &TaoUniversePointerStruct::set_u,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoUniversePointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoUniversePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoUniversePointerStruct &self, nb::dict &memo) {
            return TaoUniversePointerStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const TaoUniversePointerStruct &self, const TaoUniversePointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_universe_struct(nb::module_ &m, nb::class_<TaoUniverseStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoUniverseStruct *self,
            const TaoLatticeStruct *model,
            const TaoLatticeStruct *design,
            const TaoLatticeStruct *base,
            const TaoBeamUniStruct *beam,
            const TaoDynamicApertureStruct *dynamic_aperture,
            const TaoPingScaleStruct *ping_scale,
            const LatStruct *scratch_lat,
            const TaoUniverseCalcStruct *calc,
            const LatEleOrderStruct *ele_order,
            const TaoSpinMapStruct *spin_map,
            std::optional<std::vector<std::vector<double>>> dModel_dVar,
            std::optional<int> ix_uni,
            std::optional<int> n_d2_data_used,
            std::optional<int> n_data_used,
            std::optional<bool> is_on,
            std::optional<bool> design_same_as_previous,
            std::optional<bool> picked_uni) {
           new (self) TaoUniverseStruct(
               ptr_to_opt_ref(model),
               ptr_to_opt_ref(design),
               ptr_to_opt_ref(base),
               ptr_to_opt_ref(beam),
               ptr_to_opt_ref(dynamic_aperture),
               ptr_to_opt_ref(ping_scale),
               ptr_to_opt_ref(scratch_lat),
               ptr_to_opt_ref(calc),
               ptr_to_opt_ref(ele_order),
               ptr_to_opt_ref(spin_map),
               dModel_dVar,
               ix_uni,
               n_d2_data_used,
               n_data_used,
               is_on,
               design_same_as_previous,
               picked_uni
           );
         },
         nb::arg("model") = nb::none(),
         nb::arg("design") = nb::none(),
         nb::arg("base") = nb::none(),
         nb::arg("beam") = nb::none(),
         nb::arg("dynamic_aperture") = nb::none(),
         nb::arg("ping_scale") = nb::none(),
         nb::arg("scratch_lat") = nb::none(),
         nb::arg("calc") = nb::none(),
         nb::arg("ele_order") = nb::none(),
         nb::arg("spin_map") = nb::none(),
         nb::arg("dModel_dVar") = nb::none(),
         nb::arg("ix_uni") = nb::none(),
         nb::arg("n_d2_data_used") = nb::none(),
         nb::arg("n_data_used") = nb::none(),
         nb::arg("is_on") = nb::none(),
         nb::arg("design_same_as_previous") = nb::none(),
         nb::arg("picked_uni") = nb::none()
  )
      .def_prop_rw(
          "model",
          &TaoUniverseStruct::model,
          &TaoUniverseStruct::set_model,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "design",
          &TaoUniverseStruct::design,
          &TaoUniverseStruct::set_design,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "base",
          &TaoUniverseStruct::base,
          &TaoUniverseStruct::set_base,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "beam",
          &TaoUniverseStruct::beam,
          &TaoUniverseStruct::set_beam,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dynamic_aperture",
          &TaoUniverseStruct::dynamic_aperture,
          &TaoUniverseStruct::set_dynamic_aperture,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro(
          "model_branch",
          &TaoUniverseStruct::model_branch,
          nb::keep_alive<0, 1>(),
          "model specific information"
      )
      .def_prop_ro("d2_data", &TaoUniverseStruct::d2_data, nb::keep_alive<0, 1>(), "The data types")
      .def_prop_ro("data", &TaoUniverseStruct::data, nb::keep_alive<0, 1>(), "Array of all data.")
      .def_prop_rw(
          "ping_scale",
          &TaoUniverseStruct::ping_scale,
          &TaoUniverseStruct::set_ping_scale,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "scratch_lat",
          &TaoUniverseStruct::scratch_lat,
          &TaoUniverseStruct::set_scratch_lat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Scratch area."
      )
      .def_prop_rw(
          "calc",
          &TaoUniverseStruct::calc,
          &TaoUniverseStruct::set_calc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "What needs to be calculated?"
      )
      .def_prop_rw(
          "ele_order",
          &TaoUniverseStruct::ele_order,
          &TaoUniverseStruct::set_ele_order,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Order of elements with same name."
      )
      .def_prop_rw(
          "spin_map",
          &TaoUniverseStruct::spin_map,
          &TaoUniverseStruct::set_spin_map,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dModel_dVar",
          &TaoUniverseStruct::dModel_dVar,
          &TaoUniverseStruct::set_dModel_dVar,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Derivative matrix."
      )
      .def_prop_rw(
          "ix_uni",
          &TaoUniverseStruct::ix_uni,
          &TaoUniverseStruct::set_ix_uni,
          "Universe index."
      )
      .def_prop_rw(
          "n_d2_data_used",
          &TaoUniverseStruct::n_d2_data_used,
          &TaoUniverseStruct::set_n_d2_data_used,
          "Number of used %d2_data(:) components."
      )
      .def_prop_rw(
          "n_data_used",
          &TaoUniverseStruct::n_data_used,
          &TaoUniverseStruct::set_n_data_used,
          "Number of used %data(:) components."
      )
      .def_prop_rw(
          "is_on",
          &TaoUniverseStruct::is_on,
          &TaoUniverseStruct::set_is_on,
          "universe turned on"
      )
      .def_prop_rw(
          "design_same_as_previous",
          &TaoUniverseStruct::design_same_as_previous,
          &TaoUniverseStruct::set_design_same_as_previous,
          "Design lat same as the previous uni?"
      )
      .def_prop_rw(
          "picked_uni",
          &TaoUniverseStruct::picked_uni,
          &TaoUniverseStruct::set_picked_uni,
          "Scratch logical."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoUniverseStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoUniverseStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoUniverseStruct &self, nb::dict &memo) { return TaoUniverseStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoUniverseStruct &self, const TaoUniverseStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_v1_var_array_struct
void init_tao_v1_var_array_struct(nb::module_ &m, nb::class_<TaoV1VarArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoV1VarArrayStruct *self, const TaoV1VarStruct *v1) {
           new (self) TaoV1VarArrayStruct(ptr_to_opt_ref(v1));
         },
         nb::arg("v1") = nb::none()
  )
      .def_prop_rw(
          "v1",
          &TaoV1VarArrayStruct::v1,
          &TaoV1VarArrayStruct::set_v1,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoV1VarArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoV1VarArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoV1VarArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoV1VarArrayStruct &self) {
            return TaoV1VarArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoV1VarArrayStruct &self, nb::dict &memo) { return TaoV1VarArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoV1VarArrayStruct &self, const TaoV1VarArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoV1VarArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoV1VarArrayStructArray1D, TaoV1VarArrayStructAlloc1D>(
      m,
      "TaoV1VarArrayStructArray1D",
      "TaoV1VarArrayStructAlloc1D"
  );
  // 2D TaoV1VarArrayStruct arrays are not used in structs/routines
  // 3D TaoV1VarArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_v1_var_struct
void init_tao_v1_var_struct(nb::module_ &m, nb::class_<TaoV1VarStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("ix_v1_var") = nb::none()
  )
      .def_prop_rw(
          "name",
          &TaoV1VarStruct::name,
          &TaoV1VarStruct::set_name,
          "V1 variable name. Eg: 'quad_k1'."
      )
      .def_prop_rw(
          "ix_v1_var",
          &TaoV1VarStruct::ix_v1_var,
          &TaoV1VarStruct::set_ix_v1_var,
          "Index to s%v1_var(:) array"
      )
      .def_prop_ro(
          "v",
          &TaoV1VarStruct::v,
          nb::keep_alive<0, 1>(),
          "Pointer to the appropriate section in s%var."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoV1VarStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoV1VarStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoV1VarStruct &self, nb::dict &memo) { return TaoV1VarStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoV1VarStruct &self, const TaoV1VarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_var_array_struct
void init_tao_var_array_struct(nb::module_ &m, nb::class_<TaoVarArrayStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoVarArrayStruct *self, const TaoVarStruct *v) {
           new (self) TaoVarArrayStruct(ptr_to_opt_ref(v));
         },
         nb::arg("v") = nb::none()
  )
      .def_prop_rw(
          "v",
          &TaoVarArrayStruct::v,
          &TaoVarArrayStruct::set_v,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoVarArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoVarArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoVarArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoVarArrayStruct &self) {
            return TaoVarArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoVarArrayStruct &self, nb::dict &memo) { return TaoVarArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoVarArrayStruct &self, const TaoVarArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoVarArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoVarArrayStructArray1D, TaoVarArrayStructAlloc1D>(
      m,
      "TaoVarArrayStructArray1D",
      "TaoVarArrayStructAlloc1D"
  );
  // 2D TaoVarArrayStruct arrays are not used in structs/routines
  // 3D TaoVarArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// tao_var_slave_struct
void init_tao_var_slave_struct(nb::module_ &m, nb::class_<TaoVarSlaveStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("ix_uni") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("model_value") = nb::none(),
         nb::arg("base_value") = nb::none()
  )
      .def_prop_rw(
          "ix_uni",
          &TaoVarSlaveStruct::ix_uni,
          &TaoVarSlaveStruct::set_ix_uni,
          "universe index."
      )
      .def_prop_rw("ix_branch", &TaoVarSlaveStruct::ix_branch, &TaoVarSlaveStruct::set_ix_branch)
      .def_prop_rw(
          "ix_ele",
          &TaoVarSlaveStruct::ix_ele,
          &TaoVarSlaveStruct::set_ix_ele,
          "Index of element in the u%lattice%ele(:) array."
      )
      .def_prop_rw(
          "model_value",
          &TaoVarSlaveStruct::model_value,
          &TaoVarSlaveStruct::set_model_value,
          "Pointer to the variable in the model lat."
      )
      .def_prop_rw(
          "base_value",
          &TaoVarSlaveStruct::base_value,
          &TaoVarSlaveStruct::set_base_value,
          "Pointer to the variable in the base lat."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoVarSlaveStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoVarSlaveStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoVarSlaveStruct &self, nb::dict &memo) { return TaoVarSlaveStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoVarSlaveStruct &self, const TaoVarSlaveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_var_struct(nb::module_ &m, nb::class_<TaoVarStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoVarStruct *self,
            std::optional<std::string> ele_name,
            std::optional<std::string> attrib_name,
            std::optional<std::string> id,
            std::optional<int> ix_v1,
            std::optional<int> ix_var,
            std::optional<int> ix_dvar,
            std::optional<int> ix_attrib,
            std::optional<int> ix_key_table,
            std::optional<double> model_value,
            std::optional<double> base_value,
            std::optional<double> design_value,
            std::optional<double> scratch_value,
            std::optional<double> old_value,
            std::optional<double> meas_value,
            std::optional<double> ref_value,
            std::optional<double> correction_value,
            std::optional<double> high_lim,
            std::optional<double> low_lim,
            std::optional<double> step,
            std::optional<double> weight,
            std::optional<double> delta_merit,
            std::optional<double> merit,
            std::optional<double> dMerit_dVar,
            std::optional<double> key_val0,
            std::optional<double> key_delta,
            std::optional<double> s,
            std::optional<double> extend_val,
            std::optional<std::string> merit_type,
            std::optional<bool> exists,
            std::optional<bool> good_var,
            std::optional<bool> good_user,
            std::optional<bool> good_opt,
            std::optional<bool> good_plot,
            std::optional<bool> useit_opt,
            std::optional<bool> useit_plot,
            std::optional<bool> key_bound,
            const TaoV1VarStruct *v1) {
           new (self) TaoVarStruct(
               ele_name,
               attrib_name,
               id,
               ix_v1,
               ix_var,
               ix_dvar,
               ix_attrib,
               ix_key_table,
               model_value,
               base_value,
               design_value,
               scratch_value,
               old_value,
               meas_value,
               ref_value,
               correction_value,
               high_lim,
               low_lim,
               step,
               weight,
               delta_merit,
               merit,
               dMerit_dVar,
               key_val0,
               key_delta,
               s,
               extend_val,
               merit_type,
               exists,
               good_var,
               good_user,
               good_opt,
               good_plot,
               useit_opt,
               useit_plot,
               key_bound,
               ptr_to_opt_ref(v1)
           );
         },
         nb::arg("ele_name") = nb::none(),
         nb::arg("attrib_name") = nb::none(),
         nb::arg("id") = nb::none(),
         nb::arg("ix_v1") = nb::none(),
         nb::arg("ix_var") = nb::none(),
         nb::arg("ix_dvar") = nb::none(),
         nb::arg("ix_attrib") = nb::none(),
         nb::arg("ix_key_table") = nb::none(),
         nb::arg("model_value") = nb::none(),
         nb::arg("base_value") = nb::none(),
         nb::arg("design_value") = nb::none(),
         nb::arg("scratch_value") = nb::none(),
         nb::arg("old_value") = nb::none(),
         nb::arg("meas_value") = nb::none(),
         nb::arg("ref_value") = nb::none(),
         nb::arg("correction_value") = nb::none(),
         nb::arg("high_lim") = nb::none(),
         nb::arg("low_lim") = nb::none(),
         nb::arg("step") = nb::none(),
         nb::arg("weight") = nb::none(),
         nb::arg("delta_merit") = nb::none(),
         nb::arg("merit") = nb::none(),
         nb::arg("dMerit_dVar") = nb::none(),
         nb::arg("key_val0") = nb::none(),
         nb::arg("key_delta") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("extend_val") = nb::none(),
         nb::arg("merit_type") = nb::none(),
         nb::arg("exists") = nb::none(),
         nb::arg("good_var") = nb::none(),
         nb::arg("good_user") = nb::none(),
         nb::arg("good_opt") = nb::none(),
         nb::arg("good_plot") = nb::none(),
         nb::arg("useit_opt") = nb::none(),
         nb::arg("useit_plot") = nb::none(),
         nb::arg("key_bound") = nb::none(),
         nb::arg("v1") = nb::none()
  )
      .def_prop_rw(
          "ele_name",
          &TaoVarStruct::ele_name,
          &TaoVarStruct::set_ele_name,
          "Associated lattice element name."
      )
      .def_prop_rw(
          "attrib_name",
          &TaoVarStruct::attrib_name,
          &TaoVarStruct::set_attrib_name,
          "Name of the attribute to vary."
      )
      .def_prop_rw(
          "id",
          &TaoVarStruct::id,
          &TaoVarStruct::set_id,
          "Used by Tao extension code. Not used by Tao directly."
      )
      .def_prop_ro("slave", &TaoVarStruct::slave, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "ix_v1",
          &TaoVarStruct::ix_v1,
          &TaoVarStruct::set_ix_v1,
          "Index of this var in the s%v1_var(i)%v(:) array."
      )
      .def_prop_rw(
          "ix_var",
          &TaoVarStruct::ix_var,
          &TaoVarStruct::set_ix_var,
          "Index number of this var in the s%var(:) array."
      )
      .def_prop_rw(
          "ix_dvar",
          &TaoVarStruct::ix_dvar,
          &TaoVarStruct::set_ix_dvar,
          "Column in the dData_dVar derivative matrix."
      )
      .def_prop_rw(
          "ix_attrib",
          &TaoVarStruct::ix_attrib,
          &TaoVarStruct::set_ix_attrib,
          "Index in ele%value(:) array if appropriate."
      )
      .def_prop_rw(
          "ix_key_table",
          &TaoVarStruct::ix_key_table,
          &TaoVarStruct::set_ix_key_table,
          "Has a key binding?"
      )
      .def_prop_rw(
          "model_value",
          &TaoVarStruct::model_value,
          &TaoVarStruct::set_model_value,
          "Model value."
      )
      .def_prop_rw(
          "base_value",
          &TaoVarStruct::base_value,
          &TaoVarStruct::set_base_value,
          "Base value."
      )
      .def_prop_rw(
          "design_value",
          &TaoVarStruct::design_value,
          &TaoVarStruct::set_design_value,
          "Design value from the design lattice."
      )
      .def_prop_rw(
          "scratch_value",
          &TaoVarStruct::scratch_value,
          &TaoVarStruct::set_scratch_value,
          "Scratch space used by Tao."
      )
      .def_prop_rw(
          "old_value",
          &TaoVarStruct::old_value,
          &TaoVarStruct::set_old_value,
          "Scratch space used by Tao."
      )
      .def_prop_rw(
          "meas_value",
          &TaoVarStruct::meas_value,
          &TaoVarStruct::set_meas_value,
          "The value when the data measurement was taken."
      )
      .def_prop_rw(
          "ref_value",
          &TaoVarStruct::ref_value,
          &TaoVarStruct::set_ref_value,
          "Value when the reference measurement was taken."
      )
      .def_prop_rw(
          "correction_value",
          &TaoVarStruct::correction_value,
          &TaoVarStruct::set_correction_value,
          "Value determined by a fit to correct the lattice."
      )
      .def_prop_rw(
          "high_lim",
          &TaoVarStruct::high_lim,
          &TaoVarStruct::set_high_lim,
          "High limit for the model_value."
      )
      .def_prop_rw(
          "low_lim",
          &TaoVarStruct::low_lim,
          &TaoVarStruct::set_low_lim,
          "Low limit for the model_value."
      )
      .def_prop_rw(
          "step",
          &TaoVarStruct::step,
          &TaoVarStruct::set_step,
          "Sets what is a small step for varying this var."
      )
      .def_prop_rw(
          "weight",
          &TaoVarStruct::weight,
          &TaoVarStruct::set_weight,
          "Weight for the merit function term."
      )
      .def_prop_rw(
          "delta_merit",
          &TaoVarStruct::delta_merit,
          &TaoVarStruct::set_delta_merit,
          "Diff used to calculate the merit function term."
      )
      .def_prop_rw(
          "merit",
          &TaoVarStruct::merit,
          &TaoVarStruct::set_merit,
          "merit_term = weight * delta^2."
      )
      .def_prop_rw(
          "dMerit_dVar",
          &TaoVarStruct::dMerit_dVar,
          &TaoVarStruct::set_dMerit_dVar,
          "Merit derivative."
      )
      .def_prop_rw(
          "key_val0",
          &TaoVarStruct::key_val0,
          &TaoVarStruct::set_key_val0,
          "Key base value"
      )
      .def_prop_rw(
          "key_delta",
          &TaoVarStruct::key_delta,
          &TaoVarStruct::set_key_delta,
          "Change in value when a key is pressed."
      )
      .def_prop_rw("s", &TaoVarStruct::s, &TaoVarStruct::set_s, "longitudinal position of ele.")
      .def_prop_rw(
          "extend_val",
          &TaoVarStruct::extend_val,
          &TaoVarStruct::set_extend_val,
          "For extension code. Not used by Tao."
      )
      .def_prop_rw(
          "merit_type",
          &TaoVarStruct::merit_type,
          &TaoVarStruct::set_merit_type,
          "'target' or 'limit'"
      )
      .def_prop_rw("exists", &TaoVarStruct::exists, &TaoVarStruct::set_exists, "See above")
      .def_prop_rw("good_var", &TaoVarStruct::good_var, &TaoVarStruct::set_good_var, "See above")
      .def_prop_rw("good_user", &TaoVarStruct::good_user, &TaoVarStruct::set_good_user, "See above")
      .def_prop_rw("good_opt", &TaoVarStruct::good_opt, &TaoVarStruct::set_good_opt, "See above")
      .def_prop_rw("good_plot", &TaoVarStruct::good_plot, &TaoVarStruct::set_good_plot, "See above")
      .def_prop_rw("useit_opt", &TaoVarStruct::useit_opt, &TaoVarStruct::set_useit_opt, "See above")
      .def_prop_rw(
          "useit_plot",
          &TaoVarStruct::useit_plot,
          &TaoVarStruct::set_useit_plot,
          "See above"
      )
      .def_prop_rw(
          "key_bound",
          &TaoVarStruct::key_bound,
          &TaoVarStruct::set_key_bound,
          "Variable bound to keyboard key?"
      )
      .def_prop_rw(
          "v1",
          &TaoVarStruct::v1,
          &TaoVarStruct::set_v1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Pointer to the parent."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoVarStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoVarStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoVarStruct &self, nb::dict &memo) { return TaoVarStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoVarStruct &self, const TaoVarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_wave_kick_pt_struct(nb::module_ &m, nb::class_<TaoWaveKickPtStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoWaveKickPtStruct *self,
            std::optional<double> phi_s,
            std::optional<double> phi_r,
            std::optional<double> phi,
            std::optional<double> amp,
            std::optional<double> s,
            std::optional<int> ix_dat_before_kick,
            const EleStruct *ele) {
           new (self) TaoWaveKickPtStruct(
               phi_s,
               phi_r,
               phi,
               amp,
               s,
               ix_dat_before_kick,
               ptr_to_opt_ref(ele)
           );
         },
         nb::arg("phi_s") = nb::none(),
         nb::arg("phi_r") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("amp") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("ix_dat_before_kick") = nb::none(),
         nb::arg("ele") = nb::none()
  )
      .def_prop_rw("phi_s", &TaoWaveKickPtStruct::phi_s, &TaoWaveKickPtStruct::set_phi_s)
      .def_prop_rw("phi_r", &TaoWaveKickPtStruct::phi_r, &TaoWaveKickPtStruct::set_phi_r)
      .def_prop_rw("phi", &TaoWaveKickPtStruct::phi, &TaoWaveKickPtStruct::set_phi)
      .def_prop_rw("amp", &TaoWaveKickPtStruct::amp, &TaoWaveKickPtStruct::set_amp)
      .def_prop_rw("s", &TaoWaveKickPtStruct::s, &TaoWaveKickPtStruct::set_s, "s-position of kick")
      .def_prop_rw(
          "ix_dat_before_kick",
          &TaoWaveKickPtStruct::ix_dat_before_kick,
          &TaoWaveKickPtStruct::set_ix_dat_before_kick,
          "Index of datum in data array just before the kick."
      )
      .def_prop_rw(
          "ele",
          &TaoWaveKickPtStruct::ele,
          &TaoWaveKickPtStruct::set_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "lattice element at position of kick."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TaoWaveKickPtStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoWaveKickPtStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TaoWaveKickPtStruct &self, nb::dict &memo) { return TaoWaveKickPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoWaveKickPtStruct &self, const TaoWaveKickPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_tao_wave_struct(nb::module_ &m, nb::class_<TaoWaveStruct> &cls) {
  cls.def(
         "__init__",
         [](TaoWaveStruct *self,
            std::optional<std::string> data_type,
            std::optional<double> rms_rel_a,
            std::optional<double> rms_rel_b,
            std::optional<double> rms_rel_as,
            std::optional<double> rms_rel_bs,
            std::optional<double> rms_rel_ar,
            std::optional<double> rms_rel_br,
            std::optional<double> rms_rel_k,
            std::optional<double> rms_rel_ks,
            std::optional<double> rms_rel_kr,
            std::optional<double> rms_phi,
            std::optional<double> rms_phi_s,
            std::optional<double> rms_phi_r,
            std::optional<double> amp_ba_s,
            std::optional<double> amp_ba_r,
            std::optional<double> chi_a,
            std::optional<double> chi_c,
            std::optional<double> chi_ba,
            std::optional<std::vector<double>> amp_a,
            std::optional<std::vector<double>> amp_b,
            std::optional<std::vector<double>> amp_ba,
            std::optional<std::vector<double>> coef_a,
            std::optional<std::vector<double>> coef_b,
            std::optional<std::vector<double>> coef_ba,
            std::optional<int> n_func,
            std::optional<int> ix_a1,
            std::optional<int> ix_a2,
            std::optional<int> ix_b1,
            std::optional<int> ix_b2,
            std::optional<int> i_a1,
            std::optional<int> i_a2,
            std::optional<int> i_b1,
            std::optional<int> i_b2,
            std::optional<int> n_a,
            std::optional<int> n_b,
            std::optional<int> i_curve_wrap_pt,
            std::optional<std::vector<int>> ix_data,
            std::optional<int> n_kick,
            const TaoGraphStruct *base_graph,
            const TaoPlotRegionStruct *region,
            const TaoD1DataStruct *d1_dat) {
           new (self) TaoWaveStruct(
               data_type,
               rms_rel_a,
               rms_rel_b,
               rms_rel_as,
               rms_rel_bs,
               rms_rel_ar,
               rms_rel_br,
               rms_rel_k,
               rms_rel_ks,
               rms_rel_kr,
               rms_phi,
               rms_phi_s,
               rms_phi_r,
               amp_ba_s,
               amp_ba_r,
               chi_a,
               chi_c,
               chi_ba,
               amp_a,
               amp_b,
               amp_ba,
               coef_a,
               coef_b,
               coef_ba,
               n_func,
               ix_a1,
               ix_a2,
               ix_b1,
               ix_b2,
               i_a1,
               i_a2,
               i_b1,
               i_b2,
               n_a,
               n_b,
               i_curve_wrap_pt,
               ix_data,
               n_kick,
               ptr_to_opt_ref(base_graph),
               ptr_to_opt_ref(region),
               ptr_to_opt_ref(d1_dat)
           );
         },
         nb::arg("data_type") = nb::none(),
         nb::arg("rms_rel_a") = nb::none(),
         nb::arg("rms_rel_b") = nb::none(),
         nb::arg("rms_rel_as") = nb::none(),
         nb::arg("rms_rel_bs") = nb::none(),
         nb::arg("rms_rel_ar") = nb::none(),
         nb::arg("rms_rel_br") = nb::none(),
         nb::arg("rms_rel_k") = nb::none(),
         nb::arg("rms_rel_ks") = nb::none(),
         nb::arg("rms_rel_kr") = nb::none(),
         nb::arg("rms_phi") = nb::none(),
         nb::arg("rms_phi_s") = nb::none(),
         nb::arg("rms_phi_r") = nb::none(),
         nb::arg("amp_ba_s") = nb::none(),
         nb::arg("amp_ba_r") = nb::none(),
         nb::arg("chi_a") = nb::none(),
         nb::arg("chi_c") = nb::none(),
         nb::arg("chi_ba") = nb::none(),
         nb::arg("amp_a") = nb::none(),
         nb::arg("amp_b") = nb::none(),
         nb::arg("amp_ba") = nb::none(),
         nb::arg("coef_a") = nb::none(),
         nb::arg("coef_b") = nb::none(),
         nb::arg("coef_ba") = nb::none(),
         nb::arg("n_func") = nb::none(),
         nb::arg("ix_a1") = nb::none(),
         nb::arg("ix_a2") = nb::none(),
         nb::arg("ix_b1") = nb::none(),
         nb::arg("ix_b2") = nb::none(),
         nb::arg("i_a1") = nb::none(),
         nb::arg("i_a2") = nb::none(),
         nb::arg("i_b1") = nb::none(),
         nb::arg("i_b2") = nb::none(),
         nb::arg("n_a") = nb::none(),
         nb::arg("n_b") = nb::none(),
         nb::arg("i_curve_wrap_pt") = nb::none(),
         nb::arg("ix_data") = nb::none(),
         nb::arg("n_kick") = nb::none(),
         nb::arg("base_graph") = nb::none(),
         nb::arg("region") = nb::none(),
         nb::arg("d1_dat") = nb::none()
  )
      .def_prop_rw("data_type", &TaoWaveStruct::data_type, &TaoWaveStruct::set_data_type)
      .def_prop_rw("rms_rel_a", &TaoWaveStruct::rms_rel_a, &TaoWaveStruct::set_rms_rel_a)
      .def_prop_rw("rms_rel_b", &TaoWaveStruct::rms_rel_b, &TaoWaveStruct::set_rms_rel_b)
      .def_prop_rw("rms_rel_as", &TaoWaveStruct::rms_rel_as, &TaoWaveStruct::set_rms_rel_as)
      .def_prop_rw("rms_rel_bs", &TaoWaveStruct::rms_rel_bs, &TaoWaveStruct::set_rms_rel_bs)
      .def_prop_rw("rms_rel_ar", &TaoWaveStruct::rms_rel_ar, &TaoWaveStruct::set_rms_rel_ar)
      .def_prop_rw("rms_rel_br", &TaoWaveStruct::rms_rel_br, &TaoWaveStruct::set_rms_rel_br)
      .def_prop_rw("rms_rel_k", &TaoWaveStruct::rms_rel_k, &TaoWaveStruct::set_rms_rel_k)
      .def_prop_rw("rms_rel_ks", &TaoWaveStruct::rms_rel_ks, &TaoWaveStruct::set_rms_rel_ks)
      .def_prop_rw("rms_rel_kr", &TaoWaveStruct::rms_rel_kr, &TaoWaveStruct::set_rms_rel_kr)
      .def_prop_rw("rms_phi", &TaoWaveStruct::rms_phi, &TaoWaveStruct::set_rms_phi)
      .def_prop_rw("rms_phi_s", &TaoWaveStruct::rms_phi_s, &TaoWaveStruct::set_rms_phi_s)
      .def_prop_rw("rms_phi_r", &TaoWaveStruct::rms_phi_r, &TaoWaveStruct::set_rms_phi_r)
      .def_prop_rw("amp_ba_s", &TaoWaveStruct::amp_ba_s, &TaoWaveStruct::set_amp_ba_s)
      .def_prop_rw("amp_ba_r", &TaoWaveStruct::amp_ba_r, &TaoWaveStruct::set_amp_ba_r)
      .def_prop_rw("chi_a", &TaoWaveStruct::chi_a, &TaoWaveStruct::set_chi_a)
      .def_prop_rw("chi_c", &TaoWaveStruct::chi_c, &TaoWaveStruct::set_chi_c)
      .def_prop_rw("chi_ba", &TaoWaveStruct::chi_ba, &TaoWaveStruct::set_chi_ba)
      .def_prop_rw(
          "amp_a",
          &TaoWaveStruct::amp_a,
          &TaoWaveStruct::set_amp_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "amp_b",
          &TaoWaveStruct::amp_b,
          &TaoWaveStruct::set_amp_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "amp_ba",
          &TaoWaveStruct::amp_ba,
          &TaoWaveStruct::set_amp_ba,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "coef_a",
          &TaoWaveStruct::coef_a,
          &TaoWaveStruct::set_coef_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "coef_b",
          &TaoWaveStruct::coef_b,
          &TaoWaveStruct::set_coef_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "coef_ba",
          &TaoWaveStruct::coef_ba,
          &TaoWaveStruct::set_coef_ba,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "n_func",
          &TaoWaveStruct::n_func,
          &TaoWaveStruct::set_n_func,
          "Number of functions used in the fit."
      )
      .def_prop_rw("ix_a1", &TaoWaveStruct::ix_a1, &TaoWaveStruct::set_ix_a1)
      .def_prop_rw("ix_a2", &TaoWaveStruct::ix_a2, &TaoWaveStruct::set_ix_a2)
      .def_prop_rw("ix_b1", &TaoWaveStruct::ix_b1, &TaoWaveStruct::set_ix_b1)
      .def_prop_rw("ix_b2", &TaoWaveStruct::ix_b2, &TaoWaveStruct::set_ix_b2)
      .def_prop_rw("i_a1", &TaoWaveStruct::i_a1, &TaoWaveStruct::set_i_a1)
      .def_prop_rw("i_a2", &TaoWaveStruct::i_a2, &TaoWaveStruct::set_i_a2)
      .def_prop_rw("i_b1", &TaoWaveStruct::i_b1, &TaoWaveStruct::set_i_b1)
      .def_prop_rw("i_b2", &TaoWaveStruct::i_b2, &TaoWaveStruct::set_i_b2)
      .def_prop_rw("n_a", &TaoWaveStruct::n_a, &TaoWaveStruct::set_n_a)
      .def_prop_rw("n_b", &TaoWaveStruct::n_b, &TaoWaveStruct::set_n_b)
      .def_prop_rw(
          "i_curve_wrap_pt",
          &TaoWaveStruct::i_curve_wrap_pt,
          &TaoWaveStruct::set_i_curve_wrap_pt,
          "Index of last point before wrap in curve array."
      )
      .def_prop_rw(
          "ix_data",
          &TaoWaveStruct::ix_data,
          &TaoWaveStruct::set_ix_data,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Translates from plot point to datum index"
      )
      .def_prop_rw("n_kick", &TaoWaveStruct::n_kick, &TaoWaveStruct::set_n_kick)
      .def_prop_ro("kick", &TaoWaveStruct::kick, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "base_graph",
          &TaoWaveStruct::base_graph,
          &TaoWaveStruct::set_base_graph,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Graph before curves extended to 1.5 periods."
      )
      .def_prop_rw(
          "region",
          &TaoWaveStruct::region,
          &TaoWaveStruct::set_region,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Where the wave plot is"
      )
      .def_prop_rw(
          "d1_dat",
          &TaoWaveStruct::d1_dat,
          &TaoWaveStruct::set_d1_dat,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const TaoWaveStruct &self, nb::dict &memo) { return TaoWaveStruct(self); }
      )
      .def(
          "__eq__",
          [](const TaoWaveStruct &self, const TaoWaveStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// tao_top10_struct
void init_tao_top10_struct(nb::module_ &m, nb::class_<TaoTop10Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("name") = nb::none(),
         nb::arg("value") = nb::none(),
         nb::arg("index") = nb::none(),
         nb::arg("valid") = nb::none()
  )
      .def_prop_rw("name", &TaoTop10Struct::name, &TaoTop10Struct::set_name, "name of contributor")
      .def_prop_rw(
          "value",
          &TaoTop10Struct::value,
          &TaoTop10Struct::set_value,
          "contribution to the merit function"
      )
      .def_prop_rw(
          "index",
          &TaoTop10Struct::index,
          &TaoTop10Struct::set_index,
          "index of contributor."
      )
      .def_prop_rw("valid", &TaoTop10Struct::valid, &TaoTop10Struct::set_valid, "valid entry?")
      .def_static(
          "new_array1d",
          [](int sz) { return TaoTop10StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TaoTop10StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const TaoTop10Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const TaoTop10Struct &self) {
            return TaoTop10Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const TaoTop10Struct &self, nb::dict &memo) { return TaoTop10Struct(self); }
      )
      .def(
          "__eq__",
          [](const TaoTop10Struct &self, const TaoTop10Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const TaoTop10Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<TaoTop10StructArray1D, TaoTop10StructAlloc1D>(
      m,
      "TaoTop10StructArray1D",
      "TaoTop10StructAlloc1D"
  );
  // 2D TaoTop10Struct arrays are not used in structs/routines
  // 3D TaoTop10Struct arrays are not used in structs/routines
}

// =============================================================================
// test_sub_struct
void init_test_sub_struct(nb::module_ &m, nb::class_<TestSubStruct> &cls) {
  cls.def(
         "__init__",
         [](TestSubStruct *self, const TestSubSubStruct *sr) {
           new (self) TestSubStruct(ptr_to_opt_ref(sr));
         },
         nb::arg("sr") = nb::none()
  )
      .def_prop_rw(
          "sr",
          &TestSubStruct::sr,
          &TestSubStruct::set_sr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return TestSubStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = TestSubStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const TestSubStruct &self, nb::dict &memo) { return TestSubStruct(self); }
      )
      .def(
          "__eq__",
          [](const TestSubStruct &self, const TestSubStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_test_sub_sub_struct(nb::module_ &m, nb::class_<TestSubSubStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int64_t>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("aaa") = nb::none(),
         nb::arg("bbb") = nb::none(),
         nb::arg("file") = nb::none(),
         nb::arg("t_ref") = nb::none(),
         nb::arg("freq_spread") = nb::none()
  )
      .def_prop_rw("aaa", &TestSubSubStruct::aaa, &TestSubSubStruct::set_aaa)
      .def_prop_rw("bbb", &TestSubSubStruct::bbb, &TestSubSubStruct::set_bbb)
      .def_prop_rw("file", &TestSubSubStruct::file, &TestSubSubStruct::set_file)
      .def_prop_rw(
          "t_ref",
          &TestSubSubStruct::t_ref,
          &TestSubSubStruct::set_t_ref,
          "time reference value for computing the wake amplitude. This is used to prevent value "
          "overflow with long trains."
      )
      .def_prop_rw(
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
          [](const TestSubSubStruct &self, nb::dict &memo) { return TestSubSubStruct(self); }
      )
      .def(
          "__eq__",
          [](const TestSubSubStruct &self, const TestSubSubStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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