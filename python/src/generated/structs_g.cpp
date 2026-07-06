#include "pybmad/generated/structs_g.hpp"

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
// general_bin_struct
void init_general_bin_struct(nb::module_ &m, nb::class_<GeneralBinStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<std::vector<int>>>(),
         nb::arg("count") = nb::none(),
         nb::arg("min") = nb::none(),
         nb::arg("max") = nb::none(),
         nb::arg("delta") = nb::none(),
         nb::arg("dim") = nb::none(),
         nb::arg("n") = nb::none()
  )
      .def_prop_rw(
          "count",
          &GeneralBinStruct::count,
          &GeneralBinStruct::set_count,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Counts (or weight) in each bin"
      )
      .def_prop_rw(
          "min",
          &GeneralBinStruct::min,
          &GeneralBinStruct::set_min,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Bounds for the bins"
      )
      .def_prop_rw(
          "max",
          &GeneralBinStruct::max,
          &GeneralBinStruct::set_max,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "delta",
          &GeneralBinStruct::delta,
          &GeneralBinStruct::set_delta,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Size of a bin"
      )
      .def_prop_rw(
          "dim",
          &GeneralBinStruct::dim,
          &GeneralBinStruct::set_dim,
          "Number of dimensions"
      )
      .def_prop_rw(
          "n",
          &GeneralBinStruct::n,
          &GeneralBinStruct::set_n,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "number of bins in each dimension"
      )

      .def("__repr__", [](const GeneralBinStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GeneralBinStruct &self) {
            return GeneralBinStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GeneralBinStruct &self, nb::dict &memo) { return GeneralBinStruct(self); }
      )
      .def(
          "__eq__",
          [](const GeneralBinStruct &self, const GeneralBinStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GeneralBinStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D GeneralBinStruct arrays are not used in structs/routines
  // 2D GeneralBinStruct arrays are not used in structs/routines
  // 3D GeneralBinStruct arrays are not used in structs/routines
}

// =============================================================================
// gen_grad1_struct
void init_gen_grad1_struct(nb::module_ &m, nb::class_<GenGrad1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<std::vector<double>>>>(),
         nb::arg("m") = nb::none(),
         nb::arg("sincos") = nb::none(),
         nb::arg("n_deriv_max") = nb::none(),
         nb::arg("deriv") = nb::none()
  )
      .def_prop_rw("m", &GenGrad1Struct::m, &GenGrad1Struct::set_m, "Azimuthal index")
      .def_prop_rw("sincos", &GenGrad1Struct::sincos, &GenGrad1Struct::set_sincos, "sin$ or cos$")
      .def_prop_rw(
          "n_deriv_max",
          &GenGrad1Struct::n_deriv_max,
          &GenGrad1Struct::set_n_deriv_max,
          "Max GG derivative The derivative matrix is extended to include the interpolating spline "
          "polynomial."
      )
      .def_prop_rw(
          "deriv",
          &GenGrad1Struct::deriv,
          &GenGrad1Struct::set_deriv,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Range: (iz0:iz1, 0:2*n_deriv_max+1)"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GenGrad1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GenGrad1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GenGrad1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GenGrad1Struct &self) {
            return GenGrad1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GenGrad1Struct &self, nb::dict &memo) { return GenGrad1Struct(self); }
      )
      .def(
          "__eq__",
          [](const GenGrad1Struct &self, const GenGrad1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GenGrad1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GenGrad1StructArray1D, GenGrad1StructAlloc1D>(
      m,
      "GenGrad1StructArray1D",
      "GenGrad1StructAlloc1D"
  );
  // 2D GenGrad1Struct arrays are not used in structs/routines
  // 3D GenGrad1Struct arrays are not used in structs/routines
}

// =============================================================================
// gen_grad_map_struct
void init_gen_grad_map_struct(nb::module_ &m, nb::class_<GenGradMapStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("file") = nb::none(),
         nb::arg("ele_anchor_pt") = nb::none(),
         nb::arg("field_type") = nb::none(),
         nb::arg("iz0") = nb::none(),
         nb::arg("iz1") = nb::none(),
         nb::arg("dz") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("field_scale") = nb::none(),
         nb::arg("master_parameter") = nb::none(),
         nb::arg("curved_ref_frame") = nb::none()
  )
      .def_prop_rw(
          "file",
          &GenGradMapStruct::file,
          &GenGradMapStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_prop_ro("gg", &GenGradMapStruct::gg, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "ele_anchor_pt",
          &GenGradMapStruct::ele_anchor_pt,
          &GenGradMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_prop_rw(
          "field_type",
          &GenGradMapStruct::field_type,
          &GenGradMapStruct::set_field_type,
          "or electric$"
      )
      .def_prop_rw(
          "iz0",
          &GenGradMapStruct::iz0,
          &GenGradMapStruct::set_iz0,
          "gg%deriv(iz0:iz1, :) lower bound."
      )
      .def_prop_rw(
          "iz1",
          &GenGradMapStruct::iz1,
          &GenGradMapStruct::set_iz1,
          "gg%deriv(iz0:iz1, :) upper bound."
      )
      .def_prop_rw("dz", &GenGradMapStruct::dz, &GenGradMapStruct::set_dz, "Point spacing.")
      .def_prop_rw(
          "r0",
          &GenGradMapStruct::r0,
          &GenGradMapStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "field origin relative to ele_anchor_pt."
      )
      .def_prop_rw(
          "field_scale",
          &GenGradMapStruct::field_scale,
          &GenGradMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_prop_rw(
          "master_parameter",
          &GenGradMapStruct::master_parameter,
          &GenGradMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_prop_rw(
          "curved_ref_frame",
          &GenGradMapStruct::curved_ref_frame,
          &GenGradMapStruct::set_curved_ref_frame
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GenGradMapStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GenGradMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GenGradMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GenGradMapStruct &self) {
            return GenGradMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GenGradMapStruct &self, nb::dict &memo) { return GenGradMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const GenGradMapStruct &self, const GenGradMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GenGradMapStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GenGradMapStructArray1D, GenGradMapStructAlloc1D>(
      m,
      "GenGradMapStructArray1D",
      "GenGradMapStructAlloc1D"
  );
  // 2D GenGradMapStruct arrays are not used in structs/routines
  // 3D GenGradMapStruct arrays are not used in structs/routines
}

// =============================================================================
// gg_taylor_struct
void init_gg_taylor_struct(nb::module_ &m, nb::class_<GgTaylorStruct> &cls) {
  cls.def(nb::init<std::optional<double>>(), nb::arg("ref") = nb::none())
      .def_prop_rw("ref", &GgTaylorStruct::ref, &GgTaylorStruct::set_ref)
      .def_prop_ro("term", &GgTaylorStruct::term, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return GgTaylorStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GgTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GgTaylorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GgTaylorStruct &self) {
            return GgTaylorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GgTaylorStruct &self, nb::dict &memo) { return GgTaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const GgTaylorStruct &self, const GgTaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GgTaylorStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GgTaylorStructArray1D, GgTaylorStructAlloc1D>(
      m,
      "GgTaylorStructArray1D",
      "GgTaylorStructAlloc1D"
  );
  // 2D GgTaylorStruct arrays are not used in structs/routines
  // 3D GgTaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// gg_taylor_term_struct
void init_gg_taylor_term_struct(nb::module_ &m, nb::class_<GgTaylorTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("expn") = nb::none()
  )
      .def_prop_rw("coef", &GgTaylorTermStruct::coef, &GgTaylorTermStruct::set_coef)
      .def_prop_rw(
          "expn",
          &GgTaylorTermStruct::expn,
          &GgTaylorTermStruct::set_expn,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GgTaylorTermStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GgTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GgTaylorTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GgTaylorTermStruct &self) {
            return GgTaylorTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GgTaylorTermStruct &self, nb::dict &memo) { return GgTaylorTermStruct(self); }
      )
      .def(
          "__eq__",
          [](const GgTaylorTermStruct &self, const GgTaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GgTaylorTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GgTaylorTermStructArray1D, GgTaylorTermStructAlloc1D>(
      m,
      "GgTaylorTermStructArray1D",
      "GgTaylorTermStructAlloc1D"
  );
  // 2D GgTaylorTermStruct arrays are not used in structs/routines
  // 3D GgTaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// grid_beam_init_struct
void init_grid_beam_init_struct(nb::module_ &m, nb::class_<GridBeamInitStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("n_x") = nb::none(),
         nb::arg("n_px") = nb::none(),
         nb::arg("x_min") = nb::none(),
         nb::arg("x_max") = nb::none(),
         nb::arg("px_min") = nb::none(),
         nb::arg("px_max") = nb::none()
  )
      .def_prop_rw(
          "n_x",
          &GridBeamInitStruct::n_x,
          &GridBeamInitStruct::set_n_x,
          "Number of columns."
      )
      .def_prop_rw(
          "n_px",
          &GridBeamInitStruct::n_px,
          &GridBeamInitStruct::set_n_px,
          "Number of rows."
      )
      .def_prop_rw(
          "x_min",
          &GridBeamInitStruct::x_min,
          &GridBeamInitStruct::set_x_min,
          "Lower x limit."
      )
      .def_prop_rw(
          "x_max",
          &GridBeamInitStruct::x_max,
          &GridBeamInitStruct::set_x_max,
          "Upper x limit."
      )
      .def_prop_rw(
          "px_min",
          &GridBeamInitStruct::px_min,
          &GridBeamInitStruct::set_px_min,
          "Lower px limit."
      )
      .def_prop_rw(
          "px_max",
          &GridBeamInitStruct::px_max,
          &GridBeamInitStruct::set_px_max,
          "Upper px limit."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GridBeamInitStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GridBeamInitStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GridBeamInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridBeamInitStruct &self) {
            return GridBeamInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GridBeamInitStruct &self, nb::dict &memo) { return GridBeamInitStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridBeamInitStruct &self, const GridBeamInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GridBeamInitStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GridBeamInitStructArray1D, GridBeamInitStructAlloc1D>(
      m,
      "GridBeamInitStructArray1D",
      "GridBeamInitStructAlloc1D"
  );
  // 2D GridBeamInitStruct arrays are not used in structs/routines
  // 3D GridBeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// grid_field_pt1_struct
void init_grid_field_pt1_struct(nb::module_ &m, nb::class_<GridFieldPt1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>>(),
         nb::arg("E") = nb::none(),
         nb::arg("B") = nb::none()
  )
      .def_prop_rw(
          "E",
          &GridFieldPt1Struct::E,
          &GridFieldPt1Struct::set_E,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "B",
          &GridFieldPt1Struct::B,
          &GridFieldPt1Struct::set_B,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const GridFieldPt1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldPt1Struct &self) {
            return GridFieldPt1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GridFieldPt1Struct &self, nb::dict &memo) { return GridFieldPt1Struct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldPt1Struct &self, const GridFieldPt1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GridFieldPt1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D GridFieldPt1Struct arrays are not used in structs/routines
  // 2D GridFieldPt1Struct arrays are not used in structs/routines
  bind_FTypeArrayND<GridFieldPt1StructArray3D>(m, "GridFieldPt1StructArray3D");
}

// =============================================================================
// grid_field_pt_struct
void init_grid_field_pt_struct(nb::module_ &m, nb::class_<GridFieldPtStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("file") = nb::none(),
         nb::arg("n_link") = nb::none()
  )
      .def_prop_rw(
          "file",
          &GridFieldPtStruct::file,
          &GridFieldPtStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_prop_rw(
          "n_link",
          &GridFieldPtStruct::n_link,
          &GridFieldPtStruct::set_n_link,
          "For memory management of this structure"
      )
      .def_prop_ro("pt", &GridFieldPtStruct::pt, nb::keep_alive<0, 1>())

      .def("__repr__", [](const GridFieldPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldPtStruct &self) {
            return GridFieldPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GridFieldPtStruct &self, nb::dict &memo) { return GridFieldPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldPtStruct &self, const GridFieldPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GridFieldPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D GridFieldPtStruct arrays are not used in structs/routines
  // 2D GridFieldPtStruct arrays are not used in structs/routines
  // 3D GridFieldPtStruct arrays are not used in structs/routines
}

// =============================================================================
// grid_field_struct
void init_grid_field_struct(nb::module_ &m, nb::class_<GridFieldStruct> &cls) {
  cls.def(
         "__init__",
         [](GridFieldStruct *self,
            std::optional<int> geometry,
            std::optional<int> harmonic,
            std::optional<double> phi0_fieldmap,
            std::optional<double> field_scale,
            std::optional<int> field_type,
            std::optional<int> master_parameter,
            std::optional<int> ele_anchor_pt,
            std::optional<int> interpolation_order,
            std::optional<std::vector<double>> dr,
            std::optional<std::vector<double>> r0,
            std::optional<bool> curved_ref_frame,
            const GridFieldPtStruct *ptr) {
           new (self) GridFieldStruct(
               geometry,
               harmonic,
               phi0_fieldmap,
               field_scale,
               field_type,
               master_parameter,
               ele_anchor_pt,
               interpolation_order,
               dr,
               r0,
               curved_ref_frame,
               ptr_to_opt_ref(ptr)
           );
         },
         nb::arg("geometry") = nb::none(),
         nb::arg("harmonic") = nb::none(),
         nb::arg("phi0_fieldmap") = nb::none(),
         nb::arg("field_scale") = nb::none(),
         nb::arg("field_type") = nb::none(),
         nb::arg("master_parameter") = nb::none(),
         nb::arg("ele_anchor_pt") = nb::none(),
         nb::arg("interpolation_order") = nb::none(),
         nb::arg("dr") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("curved_ref_frame") = nb::none(),
         nb::arg("ptr") = nb::none()
  )
      .def_prop_rw(
          "geometry",
          &GridFieldStruct::geometry,
          &GridFieldStruct::set_geometry,
          "Type of grid: xyz$, or rotationally_symmetric_rz$"
      )
      .def_prop_rw(
          "harmonic",
          &GridFieldStruct::harmonic,
          &GridFieldStruct::set_harmonic,
          "Harmonic of fundamental for AC fields."
      )
      .def_prop_rw(
          "phi0_fieldmap",
          &GridFieldStruct::phi0_fieldmap,
          &GridFieldStruct::set_phi0_fieldmap,
          "Mode oscillates as: twopi * (f * t + phi0_fieldmap)"
      )
      .def_prop_rw(
          "field_scale",
          &GridFieldStruct::field_scale,
          &GridFieldStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_prop_rw(
          "field_type",
          &GridFieldStruct::field_type,
          &GridFieldStruct::set_field_type,
          "or magnetic$ or electric$"
      )
      .def_prop_rw(
          "master_parameter",
          &GridFieldStruct::master_parameter,
          &GridFieldStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_prop_rw(
          "ele_anchor_pt",
          &GridFieldStruct::ele_anchor_pt,
          &GridFieldStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_prop_rw(
          "interpolation_order",
          &GridFieldStruct::interpolation_order,
          &GridFieldStruct::set_interpolation_order,
          "Possibilities are 1 or 3."
      )
      .def_prop_rw(
          "dr",
          &GridFieldStruct::dr,
          &GridFieldStruct::set_dr,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Grid spacing."
      )
      .def_prop_rw(
          "r0",
          &GridFieldStruct::r0,
          &GridFieldStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Field origin relative to ele_anchor_pt."
      )
      .def_prop_rw(
          "curved_ref_frame",
          &GridFieldStruct::curved_ref_frame,
          &GridFieldStruct::set_curved_ref_frame
      )
      .def_prop_rw(
          "ptr",
          &GridFieldStruct::ptr,
          &GridFieldStruct::set_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro(
          "bi_coef",
          &GridFieldStruct::bi_coef,
          nb::keep_alive<0, 1>(),
          "Save computed coefs for faster tracking"
      )
      .def_prop_ro(
          "tri_coef",
          &GridFieldStruct::tri_coef,
          nb::keep_alive<0, 1>(),
          "Save computed coefs for faster tracking"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GridFieldStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GridFieldStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const GridFieldStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldStruct &self) {
            return GridFieldStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GridFieldStruct &self, nb::dict &memo) { return GridFieldStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldStruct &self, const GridFieldStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GridFieldStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<GridFieldStructArray1D, GridFieldStructAlloc1D>(
      m,
      "GridFieldStructArray1D",
      "GridFieldStructAlloc1D"
  );
  // 2D GridFieldStruct arrays are not used in structs/routines
  // 3D GridFieldStruct arrays are not used in structs/routines
}

// =============================================================================
// gpt_lat_param_struct
void init_gpt_lat_param_struct(nb::module_ &m, nb::class_<GptLatParamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<bool>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>>(),
         nb::arg("fieldmap_dimension") = nb::none(),
         nb::arg("only_write_autophase_parameters") = nb::none(),
         nb::arg("gpt_filename") = nb::none(),
         nb::arg("header_file_name") = nb::none(),
         nb::arg("tracking_end_element") = nb::none()
  )
      .def_prop_rw(
          "fieldmap_dimension",
          &GptLatParamStruct::fieldmap_dimension,
          &GptLatParamStruct::set_fieldmap_dimension,
          "Dimensions for field map. 1 or 3"
      )
      .def_prop_rw(
          "only_write_autophase_parameters",
          &GptLatParamStruct::only_write_autophase_parameters,
          &GptLatParamStruct::set_only_write_autophase_parameters,
          "Option to only write phasing info"
      )
      .def_prop_rw(
          "gpt_filename",
          &GptLatParamStruct::gpt_filename,
          &GptLatParamStruct::set_gpt_filename,
          "Blank => Append '.gpt' to Bmad lattice file name."
      )
      .def_prop_rw(
          "header_file_name",
          &GptLatParamStruct::header_file_name,
          &GptLatParamStruct::set_header_file_name,
          "Header file to include in gpt file."
      )
      .def_prop_rw(
          "tracking_end_element",
          &GptLatParamStruct::tracking_end_element,
          &GptLatParamStruct::set_tracking_end_element,
          "Bmad lattice element name or index."
      )

      .def("__repr__", [](const GptLatParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GptLatParamStruct &self) {
            return GptLatParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GptLatParamStruct &self, nb::dict &memo) { return GptLatParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const GptLatParamStruct &self, const GptLatParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const GptLatParamStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D GptLatParamStruct arrays are not used in structs/routines
  // 2D GptLatParamStruct arrays are not used in structs/routines
  // 3D GptLatParamStruct arrays are not used in structs/routines
}