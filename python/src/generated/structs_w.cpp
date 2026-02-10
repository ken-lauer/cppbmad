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
      .def_property(
          "freq",
          &WakeLrModeStruct::freq,
          &WakeLrModeStruct::set_freq,
          "Actual Frequency in Hz."
      )
      .def_property(
          "freq_in",
          &WakeLrModeStruct::freq_in,
          &WakeLrModeStruct::set_freq_in,
          "Input frequency in Hz."
      )
      .def_property(
          "R_over_Q",
          &WakeLrModeStruct::R_over_Q,
          &WakeLrModeStruct::set_R_over_Q,
          "Strength in V/C/m^(2*m_mode)."
      )
      .def_property(
          "Q",
          &WakeLrModeStruct::Q,
          &WakeLrModeStruct::set_Q,
          "Used for backwards compatability."
      )
      .def_property(
          "damp",
          &WakeLrModeStruct::damp,
          &WakeLrModeStruct::set_damp,
          "Damping factor = omega / 2 * Q = pi * freq / Q"
      )
      .def_property(
          "phi",
          &WakeLrModeStruct::phi,
          &WakeLrModeStruct::set_phi,
          "Phase in radians/2pi."
      )
      .def_property(
          "angle",
          &WakeLrModeStruct::angle,
          &WakeLrModeStruct::set_angle,
          "polarization angle (radians/2pi)."
      )
      .def_property(
          "b_sin",
          &WakeLrModeStruct::b_sin,
          &WakeLrModeStruct::set_b_sin,
          "non-skew sin-like component of the wake."
      )
      .def_property(
          "b_cos",
          &WakeLrModeStruct::b_cos,
          &WakeLrModeStruct::set_b_cos,
          "non-skew cos-like component of the wake."
      )
      .def_property(
          "a_sin",
          &WakeLrModeStruct::a_sin,
          &WakeLrModeStruct::set_a_sin,
          "skew sin-like component of the wake."
      )
      .def_property(
          "a_cos",
          &WakeLrModeStruct::a_cos,
          &WakeLrModeStruct::set_a_cos,
          "skew cos-like component of the wake."
      )
      .def_property(
          "m",
          &WakeLrModeStruct::m,
          &WakeLrModeStruct::set_m,
          "Mode order (1 = dipole, 2 = quad, etc.)"
      )
      .def_property(
          "polarized",
          &WakeLrModeStruct::polarized,
          &WakeLrModeStruct::set_polarized,
          "Polaraized mode?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return WakeLrModeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = WakeLrModeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<WakeLrModeStructArray1D, WakeLrModeStructAlloc1D>(
      m,
      "WakeLrModeStructArray1D",
      "WakeLrModeStructAlloc1D"
  );
  // 2D WakeLrModeStruct arrays are not used in structs/routines
  // 3D WakeLrModeStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_lr_struct
void init_wake_lr_struct(py::module &m, py::class_<WakeLrStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
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
      .def_property("file", &WakeLrStruct::file, &WakeLrStruct::set_file)
      .def_property_readonly("mode", &WakeLrStruct::mode, py::keep_alive<0, 1>())
      .def_property(
          "t_ref",
          &WakeLrStruct::t_ref,
          &WakeLrStruct::set_t_ref,
          "time reference value for computing the wake amplitude. This is used to prevent value "
          "overflow with long trains."
      )
      .def_property(
          "freq_spread",
          &WakeLrStruct::freq_spread,
          &WakeLrStruct::set_freq_spread,
          "Random frequency spread of long range modes."
      )
      .def_property(
          "amp_scale",
          &WakeLrStruct::amp_scale,
          &WakeLrStruct::set_amp_scale,
          "Wake amplitude scale factor."
      )
      .def_property(
          "time_scale",
          &WakeLrStruct::time_scale,
          &WakeLrStruct::set_time_scale,
          "time scale factor."
      )
      .def_property(
          "self_wake_on",
          &WakeLrStruct::self_wake_on,
          &WakeLrStruct::set_self_wake_on,
          "Long range self-wake used in tracking?"
      )

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
      .def_property("amp", &WakeSrModeStruct::amp, &WakeSrModeStruct::set_amp, "Amplitude")
      .def_property(
          "damp",
          &WakeSrModeStruct::damp,
          &WakeSrModeStruct::set_damp,
          "Dampling factor."
      )
      .def_property("k", &WakeSrModeStruct::k, &WakeSrModeStruct::set_k, "k factor")
      .def_property(
          "phi",
          &WakeSrModeStruct::phi,
          &WakeSrModeStruct::set_phi,
          "Phase in radians/2pi"
      )
      .def_property(
          "b_sin",
          &WakeSrModeStruct::b_sin,
          &WakeSrModeStruct::set_b_sin,
          "non-skew (x) sin-like component of the wake"
      )
      .def_property(
          "b_cos",
          &WakeSrModeStruct::b_cos,
          &WakeSrModeStruct::set_b_cos,
          "non-skew (x) cos-like component of the wake"
      )
      .def_property(
          "a_sin",
          &WakeSrModeStruct::a_sin,
          &WakeSrModeStruct::set_a_sin,
          "skew (y) sin-like component of the wake"
      )
      .def_property(
          "a_cos",
          &WakeSrModeStruct::a_cos,
          &WakeSrModeStruct::set_a_cos,
          "skew (y) cos-like component of the wake"
      )
      .def_property(
          "polarization",
          &WakeSrModeStruct::polarization,
          &WakeSrModeStruct::set_polarization,
          "Transverse: none$, x_axis$, y_axis$. Not used for longitudinal."
      )
      .def_property(
          "position_dependence",
          &WakeSrModeStruct::position_dependence,
          &WakeSrModeStruct::set_position_dependence,
          "Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return WakeSrModeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = WakeSrModeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<WakeSrModeStructArray1D, WakeSrModeStructAlloc1D>(
      m,
      "WakeSrModeStructArray1D",
      "WakeSrModeStructAlloc1D"
  );
  // 2D WakeSrModeStruct arrays are not used in structs/routines
  // 3D WakeSrModeStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_struct
void init_wake_sr_struct(py::module &m, py::class_<WakeSrStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
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
      .def_property("file", &WakeSrStruct::file, &WakeSrStruct::set_file)
      .def_property(
          "z_long",
          &WakeSrStruct::z_long,
          &WakeSrStruct::set_z_long,
          py::keep_alive<0, 1>()
      )
      .def_property_readonly("long_wake", &WakeSrStruct::long_wake, py::keep_alive<0, 1>())
      .def_property_readonly("trans_wake", &WakeSrStruct::trans_wake, py::keep_alive<0, 1>())
      .def_property(
          "z_ref_long",
          &WakeSrStruct::z_ref_long,
          &WakeSrStruct::set_z_ref_long,
          "z reference value for computing the wake amplitude."
      )
      .def_property(
          "z_ref_trans",
          &WakeSrStruct::z_ref_trans,
          &WakeSrStruct::set_z_ref_trans,
          "This is used to prevent value overflow with long bunches."
      )
      .def_property(
          "z_max",
          &WakeSrStruct::z_max,
          &WakeSrStruct::set_z_max,
          "Max allowable z value. 0-> ignore"
      )
      .def_property(
          "amp_scale",
          &WakeSrStruct::amp_scale,
          &WakeSrStruct::set_amp_scale,
          "Wake amplitude scale factor."
      )
      .def_property(
          "z_scale",
          &WakeSrStruct::z_scale,
          &WakeSrStruct::set_z_scale,
          "z-distance scale factor."
      )
      .def_property(
          "scale_with_length",
          &WakeSrStruct::scale_with_length,
          &WakeSrStruct::set_scale_with_length,
          "Scale wake with element length?"
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
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>,
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
      .def_property(
          "w",
          &WakeSrZLongStruct::w,
          &WakeSrZLongStruct::set_w,
          py::keep_alive<0, 1>(),
          "Input single particle Wake. Indexed from 1."
      )
      .def_property(
          "fw",
          &WakeSrZLongStruct::fw,
          &WakeSrZLongStruct::set_fw,
          py::keep_alive<0, 1>(),
          "Fourier transform of w."
      )
      .def_property(
          "fbunch",
          &WakeSrZLongStruct::fbunch,
          &WakeSrZLongStruct::set_fbunch,
          py::keep_alive<0, 1>(),
          "Scratch space."
      )
      .def_property(
          "w_out",
          &WakeSrZLongStruct::w_out,
          &WakeSrZLongStruct::set_w_out,
          py::keep_alive<0, 1>(),
          "Scratch space."
      )
      .def_property(
          "dz",
          &WakeSrZLongStruct::dz,
          &WakeSrZLongStruct::set_dz,
          "Distance between points. If zero there is no wake."
      )
      .def_property(
          "z0",
          &WakeSrZLongStruct::z0,
          &WakeSrZLongStruct::set_z0,
          "Wake extent is [-z0, z0]."
      )
      .def_property(
          "smoothing_sigma",
          &WakeSrZLongStruct::smoothing_sigma,
          &WakeSrZLongStruct::set_smoothing_sigma,
          "0 => No smoothing."
      )
      .def_property(
          "position_dependence",
          &WakeSrZLongStruct::position_dependence,
          &WakeSrZLongStruct::set_position_dependence,
          "Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$"
      )
      .def_property(
          "time_based",
          &WakeSrZLongStruct::time_based,
          &WakeSrZLongStruct::set_time_based,
          "Was input time based?"
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
      .def_property(
          "sr",
          &WakeStruct::sr,
          &WakeStruct::set_sr,
          py::keep_alive<0, 1>(),
          "Short-range wake"
      )
      .def_property(
          "lr",
          &WakeStruct::lr,
          &WakeStruct::set_lr,
          py::keep_alive<0, 1>(),
          "Long-range wake"
      )

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
             std::optional<std::string>,
             std::optional<std::string>,
             optional_ref<const PhotonReflectSurfaceStruct>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
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
      .def_property(
          "name",
          &Wall3dSectionStruct::name,
          &Wall3dSectionStruct::set_name,
          "Identifying name"
      )
      .def_property(
          "material",
          &Wall3dSectionStruct::material,
          &Wall3dSectionStruct::set_material,
          "Material."
      )
      .def_property_readonly(
          "v",
          &Wall3dSectionStruct::v,
          py::keep_alive<0, 1>(),
          "Array of vertices. Always stored relative."
      )
      .def_property(
          "surface",
          &Wall3dSectionStruct::surface,
          &Wall3dSectionStruct::set_surface,
          py::keep_alive<0, 1>(),
          "Surface reflectivity tables."
      )
      .def_property(
          "type",
          &Wall3dSectionStruct::type,
          &Wall3dSectionStruct::set_type,
          "normal$, clear$, opaque$, wall_start$, wall_end$"
      )
      .def_property(
          "n_vertex_input",
          &Wall3dSectionStruct::n_vertex_input,
          &Wall3dSectionStruct::set_n_vertex_input,
          "Number of vertices specified by the user."
      )
      .def_property(
          "ix_ele",
          &Wall3dSectionStruct::ix_ele,
          &Wall3dSectionStruct::set_ix_ele,
          "index of lattice element containing section"
      )
      .def_property(
          "ix_branch",
          &Wall3dSectionStruct::ix_branch,
          &Wall3dSectionStruct::set_ix_branch,
          "Index of branch lattice element is in."
      )
      .def_property(
          "vertices_state",
          &Wall3dSectionStruct::vertices_state,
          &Wall3dSectionStruct::set_vertices_state,
          "absolute$, or shifted_to_relative$. If set to absolute$ on input, will be changed to "
          "shifted_to_relative$ by section initalizer."
      )
      .def_property(
          "patch_in_region",
          &Wall3dSectionStruct::patch_in_region,
          &Wall3dSectionStruct::set_patch_in_region,
          "Patch element exists between this section and previous one?"
      )
      .def_property(
          "thickness",
          &Wall3dSectionStruct::thickness,
          &Wall3dSectionStruct::set_thickness,
          "Material thickness."
      )
      .def_property(
          "s",
          &Wall3dSectionStruct::s,
          &Wall3dSectionStruct::set_s,
          "Longitudinal position"
      )
      .def_property(
          "r0",
          &Wall3dSectionStruct::r0,
          &Wall3dSectionStruct::set_r0,
          py::keep_alive<0, 1>(),
          "Center of section Section-to-section spline interpolation of the center of the section"
      )
      .def_property(
          "dx0_ds",
          &Wall3dSectionStruct::dx0_ds,
          &Wall3dSectionStruct::set_dx0_ds,
          "Center of wall derivative"
      )
      .def_property(
          "dy0_ds",
          &Wall3dSectionStruct::dy0_ds,
          &Wall3dSectionStruct::set_dy0_ds,
          "Center of wall derivative"
      )
      .def_property(
          "x0_coef",
          &Wall3dSectionStruct::x0_coef,
          &Wall3dSectionStruct::set_x0_coef,
          py::keep_alive<0, 1>(),
          "Spline coefs for x-center"
      )
      .def_property(
          "y0_coef",
          &Wall3dSectionStruct::y0_coef,
          &Wall3dSectionStruct::set_y0_coef,
          py::keep_alive<0, 1>(),
          "Spline coefs for y-center Section-to_section spline interpolation of the wall."
      )
      .def_property(
          "dr_ds",
          &Wall3dSectionStruct::dr_ds,
          &Wall3dSectionStruct::set_dr_ds,
          "derivative of wall radius"
      )
      .def_property(
          "p1_coef",
          &Wall3dSectionStruct::p1_coef,
          &Wall3dSectionStruct::set_p1_coef,
          py::keep_alive<0, 1>(),
          "Spline coefs for p0 function"
      )
      .def_property(
          "p2_coef",
          &Wall3dSectionStruct::p2_coef,
          &Wall3dSectionStruct::set_p2_coef,
          py::keep_alive<0, 1>(),
          "Spline coefs for p1 function"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dSectionStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dSectionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<Wall3dSectionStructArray1D, Wall3dSectionStructAlloc1D>(
      m,
      "Wall3DSectionStructArray1D",
      "Wall3DSectionStructAlloc1D"
  );
  // 2D Wall3dSectionStruct arrays are not used in structs/routines
  // 3D Wall3dSectionStruct arrays are not used in structs/routines
}

// =============================================================================
// wall3d_struct
void init_wall3d_struct(py::module &m, py::class_<Wall3dStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
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
      .def_property("name", &Wall3dStruct::name, &Wall3dStruct::set_name)
      .def_property("type", &Wall3dStruct::type, &Wall3dStruct::set_type, "or mask_plate$")
      .def_property(
          "ix_wall3d",
          &Wall3dStruct::ix_wall3d,
          &Wall3dStruct::set_ix_wall3d,
          "Index in branch%wall3d(:) array."
      )
      .def_property(
          "n_link",
          &Wall3dStruct::n_link,
          &Wall3dStruct::set_n_link,
          "For memory management of ele%wall3d"
      )
      .def_property(
          "thickness",
          &Wall3dStruct::thickness,
          &Wall3dStruct::set_thickness,
          "For diffraction_plate elements"
      )
      .def_property(
          "clear_material",
          &Wall3dStruct::clear_material,
          &Wall3dStruct::set_clear_material
      )
      .def_property(
          "opaque_material",
          &Wall3dStruct::opaque_material,
          &Wall3dStruct::set_opaque_material
      )
      .def_property(
          "superimpose",
          &Wall3dStruct::superimpose,
          &Wall3dStruct::set_superimpose,
          "Can overlap another wall"
      )
      .def_property(
          "ele_anchor_pt",
          &Wall3dStruct::ele_anchor_pt,
          &Wall3dStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_property_readonly(
          "section",
          &Wall3dStruct::section,
          py::keep_alive<0, 1>(),
          "Indexed from 1."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<Wall3dStructArray1D, Wall3dStructAlloc1D>(
      m,
      "Wall3DStructArray1D",
      "Wall3DStructAlloc1D"
  );
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
      .def_property(
          "x",
          &Wall3dVertexStruct::x,
          &Wall3dVertexStruct::set_x,
          "Coordinates of the vertex."
      )
      .def_property(
          "y",
          &Wall3dVertexStruct::y,
          &Wall3dVertexStruct::set_y,
          "Coordinates of the vertex."
      )
      .def_property(
          "radius_x",
          &Wall3dVertexStruct::radius_x,
          &Wall3dVertexStruct::set_radius_x,
          "Radius of arc or ellipse x-axis half width. 0 => Straight line."
      )
      .def_property(
          "radius_y",
          &Wall3dVertexStruct::radius_y,
          &Wall3dVertexStruct::set_radius_y,
          "Ellipse y-axis half height."
      )
      .def_property(
          "tilt",
          &Wall3dVertexStruct::tilt,
          &Wall3dVertexStruct::set_tilt,
          "Tilt of ellipse"
      )
      .def_property(
          "angle",
          &Wall3dVertexStruct::angle,
          &Wall3dVertexStruct::set_angle,
          "Angle of (x, y) point."
      )
      .def_property("x0", &Wall3dVertexStruct::x0, &Wall3dVertexStruct::set_x0, "Center of ellipse")
      .def_property("y0", &Wall3dVertexStruct::y0, &Wall3dVertexStruct::set_y0, "Center of ellipse")
      .def_property(
          "type",
          &Wall3dVertexStruct::type,
          &Wall3dVertexStruct::set_type,
          "No longer used."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dVertexStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dVertexStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<Wall3dVertexStructArray1D, Wall3dVertexStructAlloc1D>(
      m,
      "Wall3DVertexStructArray1D",
      "Wall3DVertexStructAlloc1D"
  );
  // 2D Wall3dVertexStruct arrays are not used in structs/routines
  // 3D Wall3dVertexStruct arrays are not used in structs/routines
}