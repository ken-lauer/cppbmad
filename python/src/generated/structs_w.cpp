#include "pybmad/generated/structs_w.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// wake_lr_mode_struct
void init_wake_lr_mode_struct(py::module &m, py::class_<WakeLrModeStruct> &cls) {
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
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("freq") = py::none(),
         py::arg("freq_in") = py::none(),
         py::arg("R_over_Q") = py::none(),
         py::arg("Q") = py::none(),
         py::arg("damp") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("angle") = py::none(),
         py::arg("b_sin") = py::none(),
         py::arg("b_cos") = py::none(),
         py::arg("a_sin") = py::none(),
         py::arg("a_cos") = py::none(),
         py::arg("m") = py::none(),
         py::arg("polarized") = py::none()
  )
      // WakeLrModeStruct.freq (0D_NOT_real - Actual Frequency in Hz.
      .def_property("freq", &WakeLrModeStruct::freq, &WakeLrModeStruct::set_freq)
      // WakeLrModeStruct.freq_in (0D_NOT_real - Input frequency in Hz.
      .def_property("freq_in", &WakeLrModeStruct::freq_in, &WakeLrModeStruct::set_freq_in)
      // WakeLrModeStruct.R_over_Q (0D_NOT_real - Strength in V/C/m^(2*m_mode).
      .def_property("R_over_Q", &WakeLrModeStruct::R_over_Q, &WakeLrModeStruct::set_R_over_Q)
      // WakeLrModeStruct.Q (0D_NOT_real - Used for backwards compatability.
      .def_property("Q", &WakeLrModeStruct::Q, &WakeLrModeStruct::set_Q)
      // WakeLrModeStruct.damp (0D_NOT_real - Damping factor = omega / 2 * Q = pi * freq / Q
      .def_property("damp", &WakeLrModeStruct::damp, &WakeLrModeStruct::set_damp)
      // WakeLrModeStruct.phi (0D_NOT_real - Phase in radians/2pi.
      .def_property("phi", &WakeLrModeStruct::phi, &WakeLrModeStruct::set_phi)
      // WakeLrModeStruct.angle (0D_NOT_real - polarization angle (radians/2pi).
      .def_property("angle", &WakeLrModeStruct::angle, &WakeLrModeStruct::set_angle)
      // WakeLrModeStruct.b_sin (0D_NOT_real - non-skew sin-like component of the wake.
      .def_property("b_sin", &WakeLrModeStruct::b_sin, &WakeLrModeStruct::set_b_sin)
      // WakeLrModeStruct.b_cos (0D_NOT_real - non-skew cos-like component of the wake.
      .def_property("b_cos", &WakeLrModeStruct::b_cos, &WakeLrModeStruct::set_b_cos)
      // WakeLrModeStruct.a_sin (0D_NOT_real - skew sin-like component of the wake.
      .def_property("a_sin", &WakeLrModeStruct::a_sin, &WakeLrModeStruct::set_a_sin)
      // WakeLrModeStruct.a_cos (0D_NOT_real - skew cos-like component of the wake.
      .def_property("a_cos", &WakeLrModeStruct::a_cos, &WakeLrModeStruct::set_a_cos)
      // WakeLrModeStruct.m (0D_NOT_integer - Mode order (1 = dipole, 2 = quad, etc.)
      .def_property("m", &WakeLrModeStruct::m, &WakeLrModeStruct::set_m)
      // WakeLrModeStruct.polarized (0D_NOT_logical - Polaraized mode?
      .def_property("polarized", &WakeLrModeStruct::polarized, &WakeLrModeStruct::set_polarized)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return WakeLrModeStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const WakeLrModeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeLrModeStruct &self) {
            return WakeLrModeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const WakeLrModeStruct &self, py::dict &memo) { return WakeLrModeStruct(self); }
      )

      ;

  bind_FTypeArrayND<WakeLrModeStructArray1D>(m, "WakeLrModeStructArray1D");
  bind_FTypeAlloc1D<WakeLrModeStructAlloc1D>(m, "WakeLrModeStructAlloc1D");
  // 2D WakeLrModeStruct arrays are not used in structs/routines
  // 3D WakeLrModeStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_lr_struct
void init_wake_lr_struct(py::module &m, py::class_<WakeLrStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>>(),
         py::arg("file") = py::none(),
         py::arg("t_ref") = py::none(),
         py::arg("freq_spread") = py::none(),
         py::arg("amp_scale") = py::none(),
         py::arg("time_scale") = py::none(),
         py::arg("self_wake_on") = py::none()
  )
      // WakeLrStruct.file (0D_NOT_character -
      .def_property("file", &WakeLrStruct::file, &WakeLrStruct::set_file)
      // WakeLrStruct.mode (1D_ALLOC_type -
      .def_property_readonly("mode", &WakeLrStruct::mode)
      // WakeLrStruct.t_ref (0D_NOT_real - time reference value for computing the wake amplitude.
      // This is used to prevent value overflow with long trains.
      .def_property("t_ref", &WakeLrStruct::t_ref, &WakeLrStruct::set_t_ref)
      // WakeLrStruct.freq_spread (0D_NOT_real - Random frequency spread of long range modes.
      .def_property("freq_spread", &WakeLrStruct::freq_spread, &WakeLrStruct::set_freq_spread)
      // WakeLrStruct.amp_scale (0D_NOT_real - Wake amplitude scale factor.
      .def_property("amp_scale", &WakeLrStruct::amp_scale, &WakeLrStruct::set_amp_scale)
      // WakeLrStruct.time_scale (0D_NOT_real - time scale factor.
      .def_property("time_scale", &WakeLrStruct::time_scale, &WakeLrStruct::set_time_scale)
      // WakeLrStruct.self_wake_on (0D_NOT_logical - Long range self-wake used in tracking?
      .def_property("self_wake_on", &WakeLrStruct::self_wake_on, &WakeLrStruct::set_self_wake_on)

      .def("__repr__", [](const WakeLrStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeLrStruct &self) {
            return WakeLrStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const WakeLrStruct &self, py::dict &memo) { return WakeLrStruct(self); }
      )

      ;

  // 1D WakeLrStruct arrays are not used in structs/routines
  // 2D WakeLrStruct arrays are not used in structs/routines
  // 3D WakeLrStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_mode_struct
void init_wake_sr_mode_struct(py::module &m, py::class_<WakeSrModeStruct> &cls) {
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
             std::optional<int>,
             std::optional<int>>(),
         py::arg("amp") = py::none(),
         py::arg("damp") = py::none(),
         py::arg("k") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("b_sin") = py::none(),
         py::arg("b_cos") = py::none(),
         py::arg("a_sin") = py::none(),
         py::arg("a_cos") = py::none(),
         py::arg("polarization") = py::none(),
         py::arg("position_dependence") = py::none()
  )
      // WakeSrModeStruct.amp (0D_NOT_real - Amplitude
      .def_property("amp", &WakeSrModeStruct::amp, &WakeSrModeStruct::set_amp)
      // WakeSrModeStruct.damp (0D_NOT_real - Dampling factor.
      .def_property("damp", &WakeSrModeStruct::damp, &WakeSrModeStruct::set_damp)
      // WakeSrModeStruct.k (0D_NOT_real - k factor
      .def_property("k", &WakeSrModeStruct::k, &WakeSrModeStruct::set_k)
      // WakeSrModeStruct.phi (0D_NOT_real - Phase in radians/2pi
      .def_property("phi", &WakeSrModeStruct::phi, &WakeSrModeStruct::set_phi)
      // WakeSrModeStruct.b_sin (0D_NOT_real - non-skew (x) sin-like component of the wake
      .def_property("b_sin", &WakeSrModeStruct::b_sin, &WakeSrModeStruct::set_b_sin)
      // WakeSrModeStruct.b_cos (0D_NOT_real - non-skew (x) cos-like component of the wake
      .def_property("b_cos", &WakeSrModeStruct::b_cos, &WakeSrModeStruct::set_b_cos)
      // WakeSrModeStruct.a_sin (0D_NOT_real - skew (y) sin-like component of the wake
      .def_property("a_sin", &WakeSrModeStruct::a_sin, &WakeSrModeStruct::set_a_sin)
      // WakeSrModeStruct.a_cos (0D_NOT_real - skew (y) cos-like component of the wake
      .def_property("a_cos", &WakeSrModeStruct::a_cos, &WakeSrModeStruct::set_a_cos)
      // WakeSrModeStruct.polarization (0D_NOT_integer - Transverse: none$, x_axis$, y_axis$. Not
      // used for longitudinal.
      .def_property(
          "polarization",
          &WakeSrModeStruct::polarization,
          &WakeSrModeStruct::set_polarization
      )
      // WakeSrModeStruct.position_dependence (0D_NOT_integer - Transverse: leading$, trailing$,
      // none$ Longitudinal: x_leading$, ..., y_trailing$, none$
      .def_property(
          "position_dependence",
          &WakeSrModeStruct::position_dependence,
          &WakeSrModeStruct::set_position_dependence
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return WakeSrModeStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const WakeSrModeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrModeStruct &self) {
            return WakeSrModeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const WakeSrModeStruct &self, py::dict &memo) { return WakeSrModeStruct(self); }
      )

      ;

  bind_FTypeArrayND<WakeSrModeStructArray1D>(m, "WakeSrModeStructArray1D");
  bind_FTypeAlloc1D<WakeSrModeStructAlloc1D>(m, "WakeSrModeStructAlloc1D");
  // 2D WakeSrModeStruct arrays are not used in structs/routines
  // 3D WakeSrModeStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_struct
void init_wake_sr_struct(py::module &m, py::class_<WakeSrStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::string>,
             optional_ref<const WakeSrZLongStruct>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>>(),
         py::arg("file") = py::none(),
         py::arg("z_long") = py::none(),
         py::arg("z_ref_long") = py::none(),
         py::arg("z_ref_trans") = py::none(),
         py::arg("z_max") = py::none(),
         py::arg("amp_scale") = py::none(),
         py::arg("z_scale") = py::none(),
         py::arg("scale_with_length") = py::none()
  )
      // WakeSrStruct.file (0D_NOT_character -
      .def_property("file", &WakeSrStruct::file, &WakeSrStruct::set_file)
      // WakeSrStruct.z_long (0D_NOT_type -
      .def_property("z_long", &WakeSrStruct::z_long, &WakeSrStruct::set_z_long)
      // WakeSrStruct.long_wake (1D_ALLOC_type -
      .def_property_readonly("long_wake", &WakeSrStruct::long_wake)
      // WakeSrStruct.trans_wake (1D_ALLOC_type -
      .def_property_readonly("trans_wake", &WakeSrStruct::trans_wake)
      // WakeSrStruct.z_ref_long (0D_NOT_real - z reference value for computing the wake amplitude.
      .def_property("z_ref_long", &WakeSrStruct::z_ref_long, &WakeSrStruct::set_z_ref_long)
      // WakeSrStruct.z_ref_trans (0D_NOT_real - This is used to prevent value overflow with long
      // bunches.
      .def_property("z_ref_trans", &WakeSrStruct::z_ref_trans, &WakeSrStruct::set_z_ref_trans)
      // WakeSrStruct.z_max (0D_NOT_real - Max allowable z value. 0-> ignore
      .def_property("z_max", &WakeSrStruct::z_max, &WakeSrStruct::set_z_max)
      // WakeSrStruct.amp_scale (0D_NOT_real - Wake amplitude scale factor.
      .def_property("amp_scale", &WakeSrStruct::amp_scale, &WakeSrStruct::set_amp_scale)
      // WakeSrStruct.z_scale (0D_NOT_real - z-distance scale factor.
      .def_property("z_scale", &WakeSrStruct::z_scale, &WakeSrStruct::set_z_scale)
      // WakeSrStruct.scale_with_length (0D_NOT_logical - Scale wake with element length?
      .def_property(
          "scale_with_length",
          &WakeSrStruct::scale_with_length,
          &WakeSrStruct::set_scale_with_length
      )

      .def("__repr__", [](const WakeSrStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrStruct &self) {
            return WakeSrStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const WakeSrStruct &self, py::dict &memo) { return WakeSrStruct(self); }
      )

      ;

  // 1D WakeSrStruct arrays are not used in structs/routines
  // 2D WakeSrStruct arrays are not used in structs/routines
  // 3D WakeSrStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_z_long_struct
void init_wake_sr_z_long_struct(py::module &m, py::class_<WakeSrZLongStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::complex<double>>>,
             optional_ref<const std::vector<std::complex<double>>>,
             optional_ref<const std::vector<std::complex<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>>(),
         py::arg("w") = py::none(),
         py::arg("fw") = py::none(),
         py::arg("fbunch") = py::none(),
         py::arg("w_out") = py::none(),
         py::arg("dz") = py::none(),
         py::arg("z0") = py::none(),
         py::arg("smoothing_sigma") = py::none(),
         py::arg("position_dependence") = py::none(),
         py::arg("time_based") = py::none()
  )
      // WakeSrZLongStruct.w (1D_ALLOC_real - Input single particle Wake. Indexed from 1.
      .def_property("w", &WakeSrZLongStruct::w, &WakeSrZLongStruct::set_w)
      // WakeSrZLongStruct.fw (1D_ALLOC_complex - Fourier transform of w.
      .def_property("fw", &WakeSrZLongStruct::fw, &WakeSrZLongStruct::set_fw)
      // WakeSrZLongStruct.fbunch (1D_ALLOC_complex - Scratch space.
      .def_property("fbunch", &WakeSrZLongStruct::fbunch, &WakeSrZLongStruct::set_fbunch)
      // WakeSrZLongStruct.w_out (1D_ALLOC_complex - Scratch space.
      .def_property("w_out", &WakeSrZLongStruct::w_out, &WakeSrZLongStruct::set_w_out)
      // WakeSrZLongStruct.dz (0D_NOT_real - Distance between points. If zero there is no wake.
      .def_property("dz", &WakeSrZLongStruct::dz, &WakeSrZLongStruct::set_dz)
      // WakeSrZLongStruct.z0 (0D_NOT_real - Wake extent is [-z0, z0].
      .def_property("z0", &WakeSrZLongStruct::z0, &WakeSrZLongStruct::set_z0)
      // WakeSrZLongStruct.smoothing_sigma (0D_NOT_real - 0 => No smoothing.
      .def_property(
          "smoothing_sigma",
          &WakeSrZLongStruct::smoothing_sigma,
          &WakeSrZLongStruct::set_smoothing_sigma
      )
      // WakeSrZLongStruct.position_dependence (0D_NOT_integer - Transverse: leading$, trailing$,
      // none$ Longitudinal: x_leading$, ..., y_trailing$, none$
      .def_property(
          "position_dependence",
          &WakeSrZLongStruct::position_dependence,
          &WakeSrZLongStruct::set_position_dependence
      )
      // WakeSrZLongStruct.time_based (0D_NOT_logical - Was input time based?
      .def_property(
          "time_based",
          &WakeSrZLongStruct::time_based,
          &WakeSrZLongStruct::set_time_based
      )

      .def("__repr__", [](const WakeSrZLongStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeSrZLongStruct &self) {
            return WakeSrZLongStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const WakeSrZLongStruct &self, py::dict &memo) { return WakeSrZLongStruct(self); }
      )

      ;

  // 1D WakeSrZLongStruct arrays are not used in structs/routines
  // 2D WakeSrZLongStruct arrays are not used in structs/routines
  // 3D WakeSrZLongStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_struct
void init_wake_struct(py::module &m, py::class_<WakeStruct> &cls) {
  cls.def(
         py::init<optional_ref<const WakeSrStruct>, optional_ref<const WakeLrStruct>>(),
         py::arg("sr") = py::none(),
         py::arg("lr") = py::none()
  )
      // WakeStruct.sr (0D_NOT_type - Short-range wake
      .def_property("sr", &WakeStruct::sr, &WakeStruct::set_sr)
      // WakeStruct.lr (0D_NOT_type - Long-range wake
      .def_property("lr", &WakeStruct::lr, &WakeStruct::set_lr)

      .def("__repr__", [](const WakeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeStruct &self) {
            return WakeStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const WakeStruct &self, py::dict &memo) { return WakeStruct(self); })

      ;

  // 1D WakeStruct arrays are not used in structs/routines
  // 2D WakeStruct arrays are not used in structs/routines
  // 3D WakeStruct arrays are not used in structs/routines
}

// =============================================================================
// wall3d_section_struct
void init_wall3d_section_struct(py::module &m, py::class_<Wall3dSectionStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             optional_ref<const PhotonReflectSurfaceStruct>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>>(),
         py::arg("name") = py::none(),
         py::arg("material") = py::none(),
         py::arg("surface") = py::none(),
         py::arg("type") = py::none(),
         py::arg("n_vertex_input") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("vertices_state") = py::none(),
         py::arg("patch_in_region") = py::none(),
         py::arg("thickness") = py::none(),
         py::arg("s") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("dx0_ds") = py::none(),
         py::arg("dy0_ds") = py::none(),
         py::arg("x0_coef") = py::none(),
         py::arg("y0_coef") = py::none(),
         py::arg("dr_ds") = py::none(),
         py::arg("p1_coef") = py::none(),
         py::arg("p2_coef") = py::none()
  )
      // Wall3dSectionStruct.name (0D_NOT_character - Identifying name
      .def_property("name", &Wall3dSectionStruct::name, &Wall3dSectionStruct::set_name)
      // Wall3dSectionStruct.material (0D_NOT_character - Material.
      .def_property("material", &Wall3dSectionStruct::material, &Wall3dSectionStruct::set_material)
      // Wall3dSectionStruct.v (1D_ALLOC_type - Array of vertices. Always stored relative.
      .def_property_readonly("v", &Wall3dSectionStruct::v)
      // Wall3dSectionStruct.surface (0D_PTR_type - Surface reflectivity tables.
      .def_property("surface", &Wall3dSectionStruct::surface, &Wall3dSectionStruct::set_surface)
      // Wall3dSectionStruct.type (0D_NOT_integer - normal$, clear$, opaque$, wall_start$, wall_end$
      .def_property("type", &Wall3dSectionStruct::type, &Wall3dSectionStruct::set_type)
      // Wall3dSectionStruct.n_vertex_input (0D_NOT_integer - Number of vertices specified by the
      // user.
      .def_property(
          "n_vertex_input",
          &Wall3dSectionStruct::n_vertex_input,
          &Wall3dSectionStruct::set_n_vertex_input
      )
      // Wall3dSectionStruct.ix_ele (0D_NOT_integer - index of lattice element containing section
      .def_property("ix_ele", &Wall3dSectionStruct::ix_ele, &Wall3dSectionStruct::set_ix_ele)
      // Wall3dSectionStruct.ix_branch (0D_NOT_integer - Index of branch lattice element is in.
      .def_property(
          "ix_branch",
          &Wall3dSectionStruct::ix_branch,
          &Wall3dSectionStruct::set_ix_branch
      )
      // Wall3dSectionStruct.vertices_state (0D_NOT_integer - absolute$, or shifted_to_relative$. If
      // set to absolute$ on input, will be changed to shifted_to_relative$ by section initalizer.
      .def_property(
          "vertices_state",
          &Wall3dSectionStruct::vertices_state,
          &Wall3dSectionStruct::set_vertices_state
      )
      // Wall3dSectionStruct.patch_in_region (0D_NOT_logical - Patch element exists between this
      // section and previous one?
      .def_property(
          "patch_in_region",
          &Wall3dSectionStruct::patch_in_region,
          &Wall3dSectionStruct::set_patch_in_region
      )
      // Wall3dSectionStruct.thickness (0D_NOT_real - Material thickness.
      .def_property(
          "thickness",
          &Wall3dSectionStruct::thickness,
          &Wall3dSectionStruct::set_thickness
      )
      // Wall3dSectionStruct.s (0D_NOT_real - Longitudinal position
      .def_property("s", &Wall3dSectionStruct::s, &Wall3dSectionStruct::set_s)
      // Wall3dSectionStruct.r0 (1D_NOT_real - Center of section Section-to-section spline
      // interpolation of the center of the section
      .def_property("r0", &Wall3dSectionStruct::r0, &Wall3dSectionStruct::set_r0)
      // Wall3dSectionStruct.dx0_ds (0D_NOT_real - Center of wall derivative
      .def_property("dx0_ds", &Wall3dSectionStruct::dx0_ds, &Wall3dSectionStruct::set_dx0_ds)
      // Wall3dSectionStruct.dy0_ds (0D_NOT_real - Center of wall derivative
      .def_property("dy0_ds", &Wall3dSectionStruct::dy0_ds, &Wall3dSectionStruct::set_dy0_ds)
      // Wall3dSectionStruct.x0_coef (1D_NOT_real - Spline coefs for x-center
      .def_property("x0_coef", &Wall3dSectionStruct::x0_coef, &Wall3dSectionStruct::set_x0_coef)
      // Wall3dSectionStruct.y0_coef (1D_NOT_real - Spline coefs for y-center Section-to_section
      // spline interpolation of the wall.
      .def_property("y0_coef", &Wall3dSectionStruct::y0_coef, &Wall3dSectionStruct::set_y0_coef)
      // Wall3dSectionStruct.dr_ds (0D_NOT_real - derivative of wall radius
      .def_property("dr_ds", &Wall3dSectionStruct::dr_ds, &Wall3dSectionStruct::set_dr_ds)
      // Wall3dSectionStruct.p1_coef (1D_NOT_real - Spline coefs for p0 function
      .def_property("p1_coef", &Wall3dSectionStruct::p1_coef, &Wall3dSectionStruct::set_p1_coef)
      // Wall3dSectionStruct.p2_coef (1D_NOT_real - Spline coefs for p1 function
      .def_property("p2_coef", &Wall3dSectionStruct::p2_coef, &Wall3dSectionStruct::set_p2_coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return Wall3dSectionStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const Wall3dSectionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dSectionStruct &self) {
            return Wall3dSectionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Wall3dSectionStruct &self, py::dict &memo) { return Wall3dSectionStruct(self); }
      )

      ;

  bind_FTypeArrayND<Wall3dSectionStructArray1D>(m, "Wall3DSectionStructArray1D");
  bind_FTypeAlloc1D<Wall3dSectionStructAlloc1D>(m, "Wall3DSectionStructAlloc1D");
  // 2D Wall3dSectionStruct arrays are not used in structs/routines
  // 3D Wall3dSectionStruct arrays are not used in structs/routines
}

// =============================================================================
// wall3d_struct
void init_wall3d_struct(py::module &m, py::class_<Wall3dStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             std::optional<bool>,
             std::optional<int>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("ix_wall3d") = py::none(),
         py::arg("n_link") = py::none(),
         py::arg("thickness") = py::none(),
         py::arg("clear_material") = py::none(),
         py::arg("opaque_material") = py::none(),
         py::arg("superimpose") = py::none(),
         py::arg("ele_anchor_pt") = py::none()
  )
      // Wall3dStruct.name (0D_NOT_character -
      .def_property("name", &Wall3dStruct::name, &Wall3dStruct::set_name)
      // Wall3dStruct.type (0D_NOT_integer - or mask_plate$
      .def_property("type", &Wall3dStruct::type, &Wall3dStruct::set_type)
      // Wall3dStruct.ix_wall3d (0D_NOT_integer - Index in branch%wall3d(:) array.
      .def_property("ix_wall3d", &Wall3dStruct::ix_wall3d, &Wall3dStruct::set_ix_wall3d)
      // Wall3dStruct.n_link (0D_NOT_integer - For memory management of ele%wall3d
      .def_property("n_link", &Wall3dStruct::n_link, &Wall3dStruct::set_n_link)
      // Wall3dStruct.thickness (0D_NOT_real - For diffraction_plate elements
      .def_property("thickness", &Wall3dStruct::thickness, &Wall3dStruct::set_thickness)
      // Wall3dStruct.clear_material (0D_NOT_character -
      .def_property(
          "clear_material",
          &Wall3dStruct::clear_material,
          &Wall3dStruct::set_clear_material
      )
      // Wall3dStruct.opaque_material (0D_NOT_character -
      .def_property(
          "opaque_material",
          &Wall3dStruct::opaque_material,
          &Wall3dStruct::set_opaque_material
      )
      // Wall3dStruct.superimpose (0D_NOT_logical - Can overlap another wall
      .def_property("superimpose", &Wall3dStruct::superimpose, &Wall3dStruct::set_superimpose)
      // Wall3dStruct.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or
      // anchor_end$
      .def_property("ele_anchor_pt", &Wall3dStruct::ele_anchor_pt, &Wall3dStruct::set_ele_anchor_pt)
      // Wall3dStruct.section (1D_ALLOC_type - Indexed from 1.
      .def_property_readonly("section", &Wall3dStruct::section)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return Wall3dStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const Wall3dStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dStruct &self) {
            return Wall3dStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Wall3dStruct &self, py::dict &memo) { return Wall3dStruct(self); }
      )

      ;

  bind_FTypeArrayND<Wall3dStructArray1D>(m, "Wall3DStructArray1D");
  bind_FTypeAlloc1D<Wall3dStructAlloc1D>(m, "Wall3DStructAlloc1D");
  // 2D Wall3dStruct arrays are not used in structs/routines
  // 3D Wall3dStruct arrays are not used in structs/routines
}

// =============================================================================
// wall3d_vertex_struct
void init_wall3d_vertex_struct(py::module &m, py::class_<Wall3dVertexStruct> &cls) {
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
             std::optional<int>>(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("radius_x") = py::none(),
         py::arg("radius_y") = py::none(),
         py::arg("tilt") = py::none(),
         py::arg("angle") = py::none(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("type") = py::none()
  )
      // Wall3dVertexStruct.x (0D_NOT_real - Coordinates of the vertex.
      .def_property("x", &Wall3dVertexStruct::x, &Wall3dVertexStruct::set_x)
      // Wall3dVertexStruct.y (0D_NOT_real - Coordinates of the vertex.
      .def_property("y", &Wall3dVertexStruct::y, &Wall3dVertexStruct::set_y)
      // Wall3dVertexStruct.radius_x (0D_NOT_real - Radius of arc or ellipse x-axis half width. 0 =>
      // Straight line.
      .def_property("radius_x", &Wall3dVertexStruct::radius_x, &Wall3dVertexStruct::set_radius_x)
      // Wall3dVertexStruct.radius_y (0D_NOT_real - Ellipse y-axis half height.
      .def_property("radius_y", &Wall3dVertexStruct::radius_y, &Wall3dVertexStruct::set_radius_y)
      // Wall3dVertexStruct.tilt (0D_NOT_real - Tilt of ellipse
      .def_property("tilt", &Wall3dVertexStruct::tilt, &Wall3dVertexStruct::set_tilt)
      // Wall3dVertexStruct.angle (0D_NOT_real - Angle of (x, y) point.
      .def_property("angle", &Wall3dVertexStruct::angle, &Wall3dVertexStruct::set_angle)
      // Wall3dVertexStruct.x0 (0D_NOT_real - Center of ellipse
      .def_property("x0", &Wall3dVertexStruct::x0, &Wall3dVertexStruct::set_x0)
      // Wall3dVertexStruct.y0 (0D_NOT_real - Center of ellipse
      .def_property("y0", &Wall3dVertexStruct::y0, &Wall3dVertexStruct::set_y0)
      // Wall3dVertexStruct.type (0D_NOT_integer - No longer used.
      .def_property("type", &Wall3dVertexStruct::type, &Wall3dVertexStruct::set_type)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return Wall3dVertexStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const Wall3dVertexStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Wall3dVertexStruct &self) {
            return Wall3dVertexStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Wall3dVertexStruct &self, py::dict &memo) { return Wall3dVertexStruct(self); }
      )

      ;

  bind_FTypeArrayND<Wall3dVertexStructArray1D>(m, "Wall3DVertexStructArray1D");
  bind_FTypeAlloc1D<Wall3dVertexStructAlloc1D>(m, "Wall3DVertexStructAlloc1D");
  // 2D Wall3dVertexStruct arrays are not used in structs/routines
  // 3D Wall3dVertexStruct arrays are not used in structs/routines
}