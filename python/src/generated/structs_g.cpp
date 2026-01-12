#include "pybmad/generated/structs_g.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

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