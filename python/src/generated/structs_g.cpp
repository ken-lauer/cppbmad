#include "pybmad/generated/structs_g.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// gen_grad1_struct
void init_gen_grad1_struct(py::module &m, py::class_<GenGrad1Struct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("m") = py::none(),
         py::arg("sincos") = py::none(),
         py::arg("n_deriv_max") = py::none(),
         py::arg("deriv") = py::none()
  )
      .def_property("m", &GenGrad1Struct::m, &GenGrad1Struct::set_m, "Azimuthal index")
      .def_property("sincos", &GenGrad1Struct::sincos, &GenGrad1Struct::set_sincos, "sin$ or cos$")
      .def_property(
          "n_deriv_max",
          &GenGrad1Struct::n_deriv_max,
          &GenGrad1Struct::set_n_deriv_max,
          "Max GG derivative The derivative matrix is extended to include the interpolating spline "
          "polynomial."
      )
      .def_property(
          "deriv",
          py::cpp_function(&GenGrad1Struct::deriv, py::keep_alive<0, 1>()),
          &GenGrad1Struct::set_deriv,
          "Range: (iz0:iz1, 0:2*n_deriv_max+1)"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GenGrad1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GenGrad1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GenGrad1Struct &self, py::dict &memo) { return GenGrad1Struct(self); }
      )
      .def(
          "__eq__",
          [](const GenGrad1Struct &self, const GenGrad1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_gen_grad_map_struct(py::module &m, py::class_<GenGradMapStruct> &cls) {
  cls.def(
         py::init<
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
         py::arg("file") = py::none(),
         py::arg("ele_anchor_pt") = py::none(),
         py::arg("field_type") = py::none(),
         py::arg("iz0") = py::none(),
         py::arg("iz1") = py::none(),
         py::arg("dz") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("field_scale") = py::none(),
         py::arg("master_parameter") = py::none(),
         py::arg("curved_ref_frame") = py::none()
  )
      .def_property(
          "file",
          &GenGradMapStruct::file,
          &GenGradMapStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_property_readonly("gg", py::cpp_function(&GenGradMapStruct::gg, py::keep_alive<0, 1>()))
      .def_property(
          "ele_anchor_pt",
          &GenGradMapStruct::ele_anchor_pt,
          &GenGradMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_property(
          "field_type",
          &GenGradMapStruct::field_type,
          &GenGradMapStruct::set_field_type,
          "or electric$"
      )
      .def_property(
          "iz0",
          &GenGradMapStruct::iz0,
          &GenGradMapStruct::set_iz0,
          "gg%deriv(iz0:iz1, :) lower bound."
      )
      .def_property(
          "iz1",
          &GenGradMapStruct::iz1,
          &GenGradMapStruct::set_iz1,
          "gg%deriv(iz0:iz1, :) upper bound."
      )
      .def_property("dz", &GenGradMapStruct::dz, &GenGradMapStruct::set_dz, "Point spacing.")
      .def_property(
          "r0",
          py::cpp_function(&GenGradMapStruct::r0, py::keep_alive<0, 1>()),
          &GenGradMapStruct::set_r0,
          "field origin relative to ele_anchor_pt."
      )
      .def_property(
          "field_scale",
          &GenGradMapStruct::field_scale,
          &GenGradMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_property(
          "master_parameter",
          &GenGradMapStruct::master_parameter,
          &GenGradMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_property(
          "curved_ref_frame",
          &GenGradMapStruct::curved_ref_frame,
          &GenGradMapStruct::set_curved_ref_frame
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GenGradMapStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GenGradMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GenGradMapStruct &self, py::dict &memo) { return GenGradMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const GenGradMapStruct &self, const GenGradMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_gg_taylor_struct(py::module &m, py::class_<GgTaylorStruct> &cls) {
  cls.def(py::init<std::optional<double>>(), py::arg("ref") = py::none())
      .def_property("ref", &GgTaylorStruct::ref, &GgTaylorStruct::set_ref)
      .def_property_readonly(
          "term",
          py::cpp_function(&GgTaylorStruct::term, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GgTaylorStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GgTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GgTaylorStruct &self, py::dict &memo) { return GgTaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const GgTaylorStruct &self, const GgTaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_gg_taylor_term_struct(py::module &m, py::class_<GgTaylorTermStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("expn") = py::none()
  )
      .def_property("coef", &GgTaylorTermStruct::coef, &GgTaylorTermStruct::set_coef)
      .def_property(
          "expn",
          py::cpp_function(&GgTaylorTermStruct::expn, py::keep_alive<0, 1>()),
          &GgTaylorTermStruct::set_expn
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GgTaylorTermStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GgTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GgTaylorTermStruct &self, py::dict &memo) { return GgTaylorTermStruct(self); }
      )
      .def(
          "__eq__",
          [](const GgTaylorTermStruct &self, const GgTaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_grid_beam_init_struct(py::module &m, py::class_<GridBeamInitStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("n_x") = py::none(),
         py::arg("n_px") = py::none(),
         py::arg("x_min") = py::none(),
         py::arg("x_max") = py::none(),
         py::arg("px_min") = py::none(),
         py::arg("px_max") = py::none()
  )
      .def_property(
          "n_x",
          &GridBeamInitStruct::n_x,
          &GridBeamInitStruct::set_n_x,
          "Number of columns."
      )
      .def_property(
          "n_px",
          &GridBeamInitStruct::n_px,
          &GridBeamInitStruct::set_n_px,
          "Number of rows."
      )
      .def_property(
          "x_min",
          &GridBeamInitStruct::x_min,
          &GridBeamInitStruct::set_x_min,
          "Lower x limit."
      )
      .def_property(
          "x_max",
          &GridBeamInitStruct::x_max,
          &GridBeamInitStruct::set_x_max,
          "Upper x limit."
      )
      .def_property(
          "px_min",
          &GridBeamInitStruct::px_min,
          &GridBeamInitStruct::set_px_min,
          "Lower px limit."
      )
      .def_property(
          "px_max",
          &GridBeamInitStruct::px_max,
          &GridBeamInitStruct::set_px_max,
          "Upper px limit."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GridBeamInitStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GridBeamInitStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GridBeamInitStruct &self, py::dict &memo) { return GridBeamInitStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridBeamInitStruct &self, const GridBeamInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_grid_field_pt1_struct(py::module &m, py::class_<GridFieldPt1Struct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>>(),
         py::arg("E") = py::none(),
         py::arg("B") = py::none()
  )
      .def_property(
          "E",
          py::cpp_function(&GridFieldPt1Struct::E, py::keep_alive<0, 1>()),
          &GridFieldPt1Struct::set_E
      )
      .def_property(
          "B",
          py::cpp_function(&GridFieldPt1Struct::B, py::keep_alive<0, 1>()),
          &GridFieldPt1Struct::set_B
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
          [](const GridFieldPt1Struct &self, py::dict &memo) { return GridFieldPt1Struct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldPt1Struct &self, const GridFieldPt1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_grid_field_pt_struct(py::module &m, py::class_<GridFieldPtStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>>(),
         py::arg("file") = py::none(),
         py::arg("n_link") = py::none()
  )
      .def_property(
          "file",
          &GridFieldPtStruct::file,
          &GridFieldPtStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_property(
          "n_link",
          &GridFieldPtStruct::n_link,
          &GridFieldPtStruct::set_n_link,
          "For memory management of this structure"
      )
      .def_property_readonly("pt", py::cpp_function(&GridFieldPtStruct::pt, py::keep_alive<0, 1>()))

      .def("__repr__", [](const GridFieldPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const GridFieldPtStruct &self) {
            return GridFieldPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const GridFieldPtStruct &self, py::dict &memo) { return GridFieldPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldPtStruct &self, const GridFieldPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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
void init_grid_field_struct(py::module &m, py::class_<GridFieldStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<bool>,
             optional_ref<const GridFieldPtStruct>>(),
         py::arg("geometry") = py::none(),
         py::arg("harmonic") = py::none(),
         py::arg("phi0_fieldmap") = py::none(),
         py::arg("field_scale") = py::none(),
         py::arg("field_type") = py::none(),
         py::arg("master_parameter") = py::none(),
         py::arg("ele_anchor_pt") = py::none(),
         py::arg("interpolation_order") = py::none(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("curved_ref_frame") = py::none(),
         py::arg("ptr") = py::none()
  )
      .def_property(
          "geometry",
          &GridFieldStruct::geometry,
          &GridFieldStruct::set_geometry,
          "Type of grid: xyz$, or rotationally_symmetric_rz$"
      )
      .def_property(
          "harmonic",
          &GridFieldStruct::harmonic,
          &GridFieldStruct::set_harmonic,
          "Harmonic of fundamental for AC fields."
      )
      .def_property(
          "phi0_fieldmap",
          &GridFieldStruct::phi0_fieldmap,
          &GridFieldStruct::set_phi0_fieldmap,
          "Mode oscillates as: twopi * (f * t + phi0_fieldmap)"
      )
      .def_property(
          "field_scale",
          &GridFieldStruct::field_scale,
          &GridFieldStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_property(
          "field_type",
          &GridFieldStruct::field_type,
          &GridFieldStruct::set_field_type,
          "or magnetic$ or electric$"
      )
      .def_property(
          "master_parameter",
          &GridFieldStruct::master_parameter,
          &GridFieldStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_property(
          "ele_anchor_pt",
          &GridFieldStruct::ele_anchor_pt,
          &GridFieldStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_property(
          "interpolation_order",
          &GridFieldStruct::interpolation_order,
          &GridFieldStruct::set_interpolation_order,
          "Possibilities are 1 or 3."
      )
      .def_property(
          "dr",
          py::cpp_function(&GridFieldStruct::dr, py::keep_alive<0, 1>()),
          &GridFieldStruct::set_dr,
          "Grid spacing."
      )
      .def_property(
          "r0",
          py::cpp_function(&GridFieldStruct::r0, py::keep_alive<0, 1>()),
          &GridFieldStruct::set_r0,
          "Field origin relative to ele_anchor_pt."
      )
      .def_property(
          "curved_ref_frame",
          &GridFieldStruct::curved_ref_frame,
          &GridFieldStruct::set_curved_ref_frame
      )
      .def_property(
          "ptr",
          py::cpp_function(&GridFieldStruct::ptr, py::keep_alive<0, 1>()),
          &GridFieldStruct::set_ptr
      )
      .def_property_readonly(
          "bi_coef",
          py::cpp_function(&GridFieldStruct::bi_coef, py::keep_alive<0, 1>()),
          "Save computed coefs for faster tracking"
      )
      .def_property_readonly(
          "tri_coef",
          py::cpp_function(&GridFieldStruct::tri_coef, py::keep_alive<0, 1>()),
          "Save computed coefs for faster tracking"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return GridFieldStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = GridFieldStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const GridFieldStruct &self, py::dict &memo) { return GridFieldStruct(self); }
      )
      .def(
          "__eq__",
          [](const GridFieldStruct &self, const GridFieldStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
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