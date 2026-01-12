#include "pybmad/generated/structs_w.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

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