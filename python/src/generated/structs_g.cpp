#include "pybmad/generated/structs_g.hpp"

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
             optional_ref<const std::vector<std::vector<double>>>>(),
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
          &GenGrad1Struct::deriv,
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

      ;

  bind_FTypeArrayND<GenGrad1StructArray1D>(m, "GenGrad1StructArray1D");
  bind_FTypeAlloc1D<GenGrad1StructAlloc1D>(m, "GenGrad1StructAlloc1D");
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
             optional_ref<const std::vector<double>>,
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
      .def_property_readonly("gg", &GenGradMapStruct::gg)
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
          &GenGradMapStruct::r0,
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

      ;

  bind_FTypeArrayND<GenGradMapStructArray1D>(m, "GenGradMapStructArray1D");
  bind_FTypeAlloc1D<GenGradMapStructAlloc1D>(m, "GenGradMapStructAlloc1D");
  // 2D GenGradMapStruct arrays are not used in structs/routines
  // 3D GenGradMapStruct arrays are not used in structs/routines
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

      ;

  bind_FTypeArrayND<GridBeamInitStructArray1D>(m, "GridBeamInitStructArray1D");
  bind_FTypeAlloc1D<GridBeamInitStructAlloc1D>(m, "GridBeamInitStructAlloc1D");
  // 2D GridBeamInitStruct arrays are not used in structs/routines
  // 3D GridBeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// grid_field_pt1_struct
void init_grid_field_pt1_struct(py::module &m, py::class_<GridFieldPt1Struct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::vector<std::complex<double>>>,
             optional_ref<const std::vector<std::complex<double>>>>(),
         py::arg("E") = py::none(),
         py::arg("B") = py::none()
  )
      .def_property("E", &GridFieldPt1Struct::E, &GridFieldPt1Struct::set_E)
      .def_property("B", &GridFieldPt1Struct::B, &GridFieldPt1Struct::set_B)

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
      .def_property_readonly("pt", &GridFieldPtStruct::pt)

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
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
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
      .def_property("dr", &GridFieldStruct::dr, &GridFieldStruct::set_dr, "Grid spacing.")
      .def_property(
          "r0",
          &GridFieldStruct::r0,
          &GridFieldStruct::set_r0,
          "Field origin relative to ele_anchor_pt."
      )
      .def_property(
          "curved_ref_frame",
          &GridFieldStruct::curved_ref_frame,
          &GridFieldStruct::set_curved_ref_frame
      )
      .def_property("ptr", &GridFieldStruct::ptr, &GridFieldStruct::set_ptr)
      .def_property_readonly(
          "bi_coef",
          &GridFieldStruct::bi_coef,
          "Save computed coefs for faster tracking"
      )
      .def_property_readonly(
          "tri_coef",
          &GridFieldStruct::tri_coef,
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

      ;

  bind_FTypeArrayND<GridFieldStructArray1D>(m, "GridFieldStructArray1D");
  bind_FTypeAlloc1D<GridFieldStructAlloc1D>(m, "GridFieldStructAlloc1D");
  // 2D GridFieldStruct arrays are not used in structs/routines
  // 3D GridFieldStruct arrays are not used in structs/routines
}