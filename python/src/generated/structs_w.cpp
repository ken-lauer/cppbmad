#include "pybmad/generated/structs_w.hpp"

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
// wake_lr_mode_struct
void init_wake_lr_mode_struct(nb::module_ &m, nb::class_<WakeLrModeStruct> &cls) {
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
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("freq") = nb::none(),
         nb::arg("freq_in") = nb::none(),
         nb::arg("R_over_Q") = nb::none(),
         nb::arg("Q") = nb::none(),
         nb::arg("damp") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("angle") = nb::none(),
         nb::arg("b_sin") = nb::none(),
         nb::arg("b_cos") = nb::none(),
         nb::arg("a_sin") = nb::none(),
         nb::arg("a_cos") = nb::none(),
         nb::arg("m") = nb::none(),
         nb::arg("polarized") = nb::none()
  )
      .def_prop_rw(
          "freq",
          &WakeLrModeStruct::freq,
          &WakeLrModeStruct::set_freq,
          "Actual Frequency in Hz."
      )
      .def_prop_rw(
          "freq_in",
          &WakeLrModeStruct::freq_in,
          &WakeLrModeStruct::set_freq_in,
          "Input frequency in Hz."
      )
      .def_prop_rw(
          "R_over_Q",
          &WakeLrModeStruct::R_over_Q,
          &WakeLrModeStruct::set_R_over_Q,
          "Strength in V/C/m^(2*m_mode)."
      )
      .def_prop_rw(
          "Q",
          &WakeLrModeStruct::Q,
          &WakeLrModeStruct::set_Q,
          "Used for backwards compatability."
      )
      .def_prop_rw(
          "damp",
          &WakeLrModeStruct::damp,
          &WakeLrModeStruct::set_damp,
          "Damping factor = omega / 2 * Q = pi * freq / Q"
      )
      .def_prop_rw(
          "phi",
          &WakeLrModeStruct::phi,
          &WakeLrModeStruct::set_phi,
          "Phase in radians/2pi."
      )
      .def_prop_rw(
          "angle",
          &WakeLrModeStruct::angle,
          &WakeLrModeStruct::set_angle,
          "polarization angle (radians/2pi)."
      )
      .def_prop_rw(
          "b_sin",
          &WakeLrModeStruct::b_sin,
          &WakeLrModeStruct::set_b_sin,
          "non-skew sin-like component of the wake."
      )
      .def_prop_rw(
          "b_cos",
          &WakeLrModeStruct::b_cos,
          &WakeLrModeStruct::set_b_cos,
          "non-skew cos-like component of the wake."
      )
      .def_prop_rw(
          "a_sin",
          &WakeLrModeStruct::a_sin,
          &WakeLrModeStruct::set_a_sin,
          "skew sin-like component of the wake."
      )
      .def_prop_rw(
          "a_cos",
          &WakeLrModeStruct::a_cos,
          &WakeLrModeStruct::set_a_cos,
          "skew cos-like component of the wake."
      )
      .def_prop_rw(
          "m",
          &WakeLrModeStruct::m,
          &WakeLrModeStruct::set_m,
          "Mode order (1 = dipole, 2 = quad, etc.)"
      )
      .def_prop_rw(
          "polarized",
          &WakeLrModeStruct::polarized,
          &WakeLrModeStruct::set_polarized,
          "Polaraized mode?"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return WakeLrModeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = WakeLrModeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const WakeLrModeStruct &self, nb::dict &memo) { return WakeLrModeStruct(self); }
      )
      .def(
          "__eq__",
          [](const WakeLrModeStruct &self, const WakeLrModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeLrModeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_wake_lr_struct(nb::module_ &m, nb::class_<WakeLrStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>>(),
         nb::arg("file") = nb::none(),
         nb::arg("t_ref") = nb::none(),
         nb::arg("freq_spread") = nb::none(),
         nb::arg("amp_scale") = nb::none(),
         nb::arg("time_scale") = nb::none(),
         nb::arg("self_wake_on") = nb::none()
  )
      .def_prop_rw("file", &WakeLrStruct::file, &WakeLrStruct::set_file)
      .def_prop_ro("mode", &WakeLrStruct::mode, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "t_ref",
          &WakeLrStruct::t_ref,
          &WakeLrStruct::set_t_ref,
          "time reference value for computing the wake amplitude. This is used to prevent value "
          "overflow with long trains."
      )
      .def_prop_rw(
          "freq_spread",
          &WakeLrStruct::freq_spread,
          &WakeLrStruct::set_freq_spread,
          "Random frequency spread of long range modes."
      )
      .def_prop_rw(
          "amp_scale",
          &WakeLrStruct::amp_scale,
          &WakeLrStruct::set_amp_scale,
          "Wake amplitude scale factor."
      )
      .def_prop_rw(
          "time_scale",
          &WakeLrStruct::time_scale,
          &WakeLrStruct::set_time_scale,
          "time scale factor."
      )
      .def_prop_rw(
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
          [](const WakeLrStruct &self, nb::dict &memo) { return WakeLrStruct(self); }
      )
      .def(
          "__eq__",
          [](const WakeLrStruct &self, const WakeLrStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeLrStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D WakeLrStruct arrays are not used in structs/routines
  // 2D WakeLrStruct arrays are not used in structs/routines
  // 3D WakeLrStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_mode_struct
void init_wake_sr_mode_struct(nb::module_ &m, nb::class_<WakeSrModeStruct> &cls) {
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
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("amp") = nb::none(),
         nb::arg("damp") = nb::none(),
         nb::arg("k") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("b_sin") = nb::none(),
         nb::arg("b_cos") = nb::none(),
         nb::arg("a_sin") = nb::none(),
         nb::arg("a_cos") = nb::none(),
         nb::arg("polarization") = nb::none(),
         nb::arg("position_dependence") = nb::none()
  )
      .def_prop_rw("amp", &WakeSrModeStruct::amp, &WakeSrModeStruct::set_amp, "Amplitude")
      .def_prop_rw("damp", &WakeSrModeStruct::damp, &WakeSrModeStruct::set_damp, "Dampling factor.")
      .def_prop_rw("k", &WakeSrModeStruct::k, &WakeSrModeStruct::set_k, "k factor")
      .def_prop_rw(
          "phi",
          &WakeSrModeStruct::phi,
          &WakeSrModeStruct::set_phi,
          "Phase in radians/2pi"
      )
      .def_prop_rw(
          "b_sin",
          &WakeSrModeStruct::b_sin,
          &WakeSrModeStruct::set_b_sin,
          "non-skew (x) sin-like component of the wake"
      )
      .def_prop_rw(
          "b_cos",
          &WakeSrModeStruct::b_cos,
          &WakeSrModeStruct::set_b_cos,
          "non-skew (x) cos-like component of the wake"
      )
      .def_prop_rw(
          "a_sin",
          &WakeSrModeStruct::a_sin,
          &WakeSrModeStruct::set_a_sin,
          "skew (y) sin-like component of the wake"
      )
      .def_prop_rw(
          "a_cos",
          &WakeSrModeStruct::a_cos,
          &WakeSrModeStruct::set_a_cos,
          "skew (y) cos-like component of the wake"
      )
      .def_prop_rw(
          "polarization",
          &WakeSrModeStruct::polarization,
          &WakeSrModeStruct::set_polarization,
          "Transverse: none$, x_axis$, y_axis$. Not used for longitudinal."
      )
      .def_prop_rw(
          "position_dependence",
          &WakeSrModeStruct::position_dependence,
          &WakeSrModeStruct::set_position_dependence,
          "Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return WakeSrModeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = WakeSrModeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const WakeSrModeStruct &self, nb::dict &memo) { return WakeSrModeStruct(self); }
      )
      .def(
          "__eq__",
          [](const WakeSrModeStruct &self, const WakeSrModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeSrModeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_wake_sr_struct(nb::module_ &m, nb::class_<WakeSrStruct> &cls) {
  cls.def(
         "__init__",
         [](WakeSrStruct *self,
            std::optional<std::string> file,
            const WakeSrZLongStruct *z_long,
            std::optional<double> z_ref_long,
            std::optional<double> z_ref_trans,
            std::optional<double> z_max,
            std::optional<double> amp_scale,
            std::optional<double> z_scale,
            std::optional<bool> scale_with_length) {
           new (self) WakeSrStruct(
               file,
               ptr_to_opt_ref(z_long),
               z_ref_long,
               z_ref_trans,
               z_max,
               amp_scale,
               z_scale,
               scale_with_length
           );
         },
         nb::arg("file") = nb::none(),
         nb::arg("z_long") = nb::none(),
         nb::arg("z_ref_long") = nb::none(),
         nb::arg("z_ref_trans") = nb::none(),
         nb::arg("z_max") = nb::none(),
         nb::arg("amp_scale") = nb::none(),
         nb::arg("z_scale") = nb::none(),
         nb::arg("scale_with_length") = nb::none()
  )
      .def_prop_rw("file", &WakeSrStruct::file, &WakeSrStruct::set_file)
      .def_prop_rw(
          "z_long",
          &WakeSrStruct::z_long,
          &WakeSrStruct::set_z_long,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("long_wake", &WakeSrStruct::long_wake, nb::keep_alive<0, 1>())
      .def_prop_ro("trans_wake", &WakeSrStruct::trans_wake, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "z_ref_long",
          &WakeSrStruct::z_ref_long,
          &WakeSrStruct::set_z_ref_long,
          "z reference value for computing the wake amplitude."
      )
      .def_prop_rw(
          "z_ref_trans",
          &WakeSrStruct::z_ref_trans,
          &WakeSrStruct::set_z_ref_trans,
          "This is used to prevent value overflow with long bunches."
      )
      .def_prop_rw(
          "z_max",
          &WakeSrStruct::z_max,
          &WakeSrStruct::set_z_max,
          "Max allowable z value. 0-> ignore"
      )
      .def_prop_rw(
          "amp_scale",
          &WakeSrStruct::amp_scale,
          &WakeSrStruct::set_amp_scale,
          "Wake amplitude scale factor."
      )
      .def_prop_rw(
          "z_scale",
          &WakeSrStruct::z_scale,
          &WakeSrStruct::set_z_scale,
          "z-distance scale factor."
      )
      .def_prop_rw(
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
          [](const WakeSrStruct &self, nb::dict &memo) { return WakeSrStruct(self); }
      )
      .def(
          "__eq__",
          [](const WakeSrStruct &self, const WakeSrStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeSrStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D WakeSrStruct arrays are not used in structs/routines
  // 2D WakeSrStruct arrays are not used in structs/routines
  // 3D WakeSrStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_sr_z_long_struct
void init_wake_sr_z_long_struct(nb::module_ &m, nb::class_<WakeSrZLongStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("w") = nb::none(),
         nb::arg("fw") = nb::none(),
         nb::arg("fbunch") = nb::none(),
         nb::arg("w_out") = nb::none(),
         nb::arg("dz") = nb::none(),
         nb::arg("z0") = nb::none(),
         nb::arg("smoothing_sigma") = nb::none(),
         nb::arg("position_dependence") = nb::none(),
         nb::arg("time_based") = nb::none()
  )
      .def_prop_rw(
          "w",
          &WakeSrZLongStruct::w,
          &WakeSrZLongStruct::set_w,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Input single particle Wake. Indexed from 1."
      )
      .def_prop_rw(
          "fw",
          &WakeSrZLongStruct::fw,
          &WakeSrZLongStruct::set_fw,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Fourier transform of w."
      )
      .def_prop_rw(
          "fbunch",
          &WakeSrZLongStruct::fbunch,
          &WakeSrZLongStruct::set_fbunch,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Scratch space."
      )
      .def_prop_rw(
          "w_out",
          &WakeSrZLongStruct::w_out,
          &WakeSrZLongStruct::set_w_out,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Scratch space."
      )
      .def_prop_rw(
          "dz",
          &WakeSrZLongStruct::dz,
          &WakeSrZLongStruct::set_dz,
          "Distance between points. If zero there is no wake."
      )
      .def_prop_rw(
          "z0",
          &WakeSrZLongStruct::z0,
          &WakeSrZLongStruct::set_z0,
          "Wake extent is [-z0, z0]."
      )
      .def_prop_rw(
          "smoothing_sigma",
          &WakeSrZLongStruct::smoothing_sigma,
          &WakeSrZLongStruct::set_smoothing_sigma,
          "0 => No smoothing."
      )
      .def_prop_rw(
          "position_dependence",
          &WakeSrZLongStruct::position_dependence,
          &WakeSrZLongStruct::set_position_dependence,
          "Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$"
      )
      .def_prop_rw(
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
          [](const WakeSrZLongStruct &self, nb::dict &memo) { return WakeSrZLongStruct(self); }
      )
      .def(
          "__eq__",
          [](const WakeSrZLongStruct &self, const WakeSrZLongStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeSrZLongStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D WakeSrZLongStruct arrays are not used in structs/routines
  // 2D WakeSrZLongStruct arrays are not used in structs/routines
  // 3D WakeSrZLongStruct arrays are not used in structs/routines
}

// =============================================================================
// wake_struct
void init_wake_struct(nb::module_ &m, nb::class_<WakeStruct> &cls) {
  cls.def(
         "__init__",
         [](WakeStruct *self, const WakeSrStruct *sr, const WakeLrStruct *lr) {
           new (self) WakeStruct(ptr_to_opt_ref(sr), ptr_to_opt_ref(lr));
         },
         nb::arg("sr") = nb::none(),
         nb::arg("lr") = nb::none()
  )
      .def_prop_rw(
          "sr",
          &WakeStruct::sr,
          &WakeStruct::set_sr,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Short-range wake"
      )
      .def_prop_rw(
          "lr",
          &WakeStruct::lr,
          &WakeStruct::set_lr,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Long-range wake"
      )

      .def("__repr__", [](const WakeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const WakeStruct &self) {
            return WakeStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const WakeStruct &self, nb::dict &memo) { return WakeStruct(self); })
      .def(
          "__eq__",
          [](const WakeStruct &self, const WakeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const WakeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D WakeStruct arrays are not used in structs/routines
  // 2D WakeStruct arrays are not used in structs/routines
  // 3D WakeStruct arrays are not used in structs/routines
}

// =============================================================================
// wall3d_section_struct
void init_wall3d_section_struct(nb::module_ &m, nb::class_<Wall3dSectionStruct> &cls) {
  cls.def(
         "__init__",
         [](Wall3dSectionStruct *self,
            std::optional<std::string> name,
            std::optional<std::string> material,
            const PhotonReflectSurfaceStruct *surface,
            std::optional<int> type,
            std::optional<int> n_vertex_input,
            std::optional<int> ix_ele,
            std::optional<int> ix_branch,
            std::optional<int> vertices_state,
            std::optional<bool> patch_in_region,
            std::optional<double> thickness,
            std::optional<double> s,
            std::optional<std::vector<double>> r0,
            std::optional<double> dx0_ds,
            std::optional<double> dy0_ds,
            std::optional<std::vector<double>> x0_coef,
            std::optional<std::vector<double>> y0_coef,
            std::optional<double> dr_ds,
            std::optional<std::vector<double>> p1_coef,
            std::optional<std::vector<double>> p2_coef) {
           new (self) Wall3dSectionStruct(
               name,
               material,
               ptr_to_opt_ref(surface),
               type,
               n_vertex_input,
               ix_ele,
               ix_branch,
               vertices_state,
               patch_in_region,
               thickness,
               s,
               r0,
               dx0_ds,
               dy0_ds,
               x0_coef,
               y0_coef,
               dr_ds,
               p1_coef,
               p2_coef
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("material") = nb::none(),
         nb::arg("surface") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("n_vertex_input") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("vertices_state") = nb::none(),
         nb::arg("patch_in_region") = nb::none(),
         nb::arg("thickness") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("dx0_ds") = nb::none(),
         nb::arg("dy0_ds") = nb::none(),
         nb::arg("x0_coef") = nb::none(),
         nb::arg("y0_coef") = nb::none(),
         nb::arg("dr_ds") = nb::none(),
         nb::arg("p1_coef") = nb::none(),
         nb::arg("p2_coef") = nb::none()
  )
      .def_prop_rw(
          "name",
          &Wall3dSectionStruct::name,
          &Wall3dSectionStruct::set_name,
          "Identifying name"
      )
      .def_prop_rw(
          "material",
          &Wall3dSectionStruct::material,
          &Wall3dSectionStruct::set_material,
          "Material."
      )
      .def_prop_ro(
          "v",
          &Wall3dSectionStruct::v,
          nb::keep_alive<0, 1>(),
          "Array of vertices. Always stored relative."
      )
      .def_prop_rw(
          "surface",
          &Wall3dSectionStruct::surface,
          &Wall3dSectionStruct::set_surface,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Surface reflectivity tables."
      )
      .def_prop_rw(
          "type",
          &Wall3dSectionStruct::type,
          &Wall3dSectionStruct::set_type,
          "normal$, clear$, opaque$, wall_start$, wall_end$"
      )
      .def_prop_rw(
          "n_vertex_input",
          &Wall3dSectionStruct::n_vertex_input,
          &Wall3dSectionStruct::set_n_vertex_input,
          "Number of vertices specified by the user."
      )
      .def_prop_rw(
          "ix_ele",
          &Wall3dSectionStruct::ix_ele,
          &Wall3dSectionStruct::set_ix_ele,
          "index of lattice element containing section"
      )
      .def_prop_rw(
          "ix_branch",
          &Wall3dSectionStruct::ix_branch,
          &Wall3dSectionStruct::set_ix_branch,
          "Index of branch lattice element is in."
      )
      .def_prop_rw(
          "vertices_state",
          &Wall3dSectionStruct::vertices_state,
          &Wall3dSectionStruct::set_vertices_state,
          "absolute$, or shifted_to_relative$. If set to absolute$ on input, will be changed to "
          "shifted_to_relative$ by section initalizer."
      )
      .def_prop_rw(
          "patch_in_region",
          &Wall3dSectionStruct::patch_in_region,
          &Wall3dSectionStruct::set_patch_in_region,
          "Patch element exists between this section and previous one?"
      )
      .def_prop_rw(
          "thickness",
          &Wall3dSectionStruct::thickness,
          &Wall3dSectionStruct::set_thickness,
          "Material thickness."
      )
      .def_prop_rw(
          "s",
          &Wall3dSectionStruct::s,
          &Wall3dSectionStruct::set_s,
          "Longitudinal position"
      )
      .def_prop_rw(
          "r0",
          &Wall3dSectionStruct::r0,
          &Wall3dSectionStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Center of section Section-to-section spline interpolation of the center of the section"
      )
      .def_prop_rw(
          "dx0_ds",
          &Wall3dSectionStruct::dx0_ds,
          &Wall3dSectionStruct::set_dx0_ds,
          "Center of wall derivative"
      )
      .def_prop_rw(
          "dy0_ds",
          &Wall3dSectionStruct::dy0_ds,
          &Wall3dSectionStruct::set_dy0_ds,
          "Center of wall derivative"
      )
      .def_prop_rw(
          "x0_coef",
          &Wall3dSectionStruct::x0_coef,
          &Wall3dSectionStruct::set_x0_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spline coefs for x-center"
      )
      .def_prop_rw(
          "y0_coef",
          &Wall3dSectionStruct::y0_coef,
          &Wall3dSectionStruct::set_y0_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spline coefs for y-center Section-to_section spline interpolation of the wall."
      )
      .def_prop_rw(
          "dr_ds",
          &Wall3dSectionStruct::dr_ds,
          &Wall3dSectionStruct::set_dr_ds,
          "derivative of wall radius"
      )
      .def_prop_rw(
          "p1_coef",
          &Wall3dSectionStruct::p1_coef,
          &Wall3dSectionStruct::set_p1_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spline coefs for p0 function"
      )
      .def_prop_rw(
          "p2_coef",
          &Wall3dSectionStruct::p2_coef,
          &Wall3dSectionStruct::set_p2_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spline coefs for p1 function"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dSectionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dSectionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const Wall3dSectionStruct &self, nb::dict &memo) { return Wall3dSectionStruct(self); }
      )
      .def(
          "__eq__",
          [](const Wall3dSectionStruct &self, const Wall3dSectionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Wall3dSectionStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_wall3d_struct(nb::module_ &m, nb::class_<Wall3dStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>,
             std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("ix_wall3d") = nb::none(),
         nb::arg("n_link") = nb::none(),
         nb::arg("thickness") = nb::none(),
         nb::arg("clear_material") = nb::none(),
         nb::arg("opaque_material") = nb::none(),
         nb::arg("superimpose") = nb::none(),
         nb::arg("ele_anchor_pt") = nb::none()
  )
      .def_prop_rw("name", &Wall3dStruct::name, &Wall3dStruct::set_name)
      .def_prop_rw("type", &Wall3dStruct::type, &Wall3dStruct::set_type, "or mask_plate$")
      .def_prop_rw(
          "ix_wall3d",
          &Wall3dStruct::ix_wall3d,
          &Wall3dStruct::set_ix_wall3d,
          "Index in branch%wall3d(:) array."
      )
      .def_prop_rw(
          "n_link",
          &Wall3dStruct::n_link,
          &Wall3dStruct::set_n_link,
          "For memory management of ele%wall3d"
      )
      .def_prop_rw(
          "thickness",
          &Wall3dStruct::thickness,
          &Wall3dStruct::set_thickness,
          "For diffraction_plate elements"
      )
      .def_prop_rw(
          "clear_material",
          &Wall3dStruct::clear_material,
          &Wall3dStruct::set_clear_material
      )
      .def_prop_rw(
          "opaque_material",
          &Wall3dStruct::opaque_material,
          &Wall3dStruct::set_opaque_material
      )
      .def_prop_rw(
          "superimpose",
          &Wall3dStruct::superimpose,
          &Wall3dStruct::set_superimpose,
          "Can overlap another wall"
      )
      .def_prop_rw(
          "ele_anchor_pt",
          &Wall3dStruct::ele_anchor_pt,
          &Wall3dStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_prop_ro("section", &Wall3dStruct::section, nb::keep_alive<0, 1>(), "Indexed from 1.")
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const Wall3dStruct &self, nb::dict &memo) { return Wall3dStruct(self); }
      )
      .def(
          "__eq__",
          [](const Wall3dStruct &self, const Wall3dStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Wall3dStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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
void init_wall3d_vertex_struct(nb::module_ &m, nb::class_<Wall3dVertexStruct> &cls) {
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
             std::optional<int>>(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("radius_x") = nb::none(),
         nb::arg("radius_y") = nb::none(),
         nb::arg("tilt") = nb::none(),
         nb::arg("angle") = nb::none(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("type") = nb::none()
  )
      .def_prop_rw(
          "x",
          &Wall3dVertexStruct::x,
          &Wall3dVertexStruct::set_x,
          "Coordinates of the vertex."
      )
      .def_prop_rw(
          "y",
          &Wall3dVertexStruct::y,
          &Wall3dVertexStruct::set_y,
          "Coordinates of the vertex."
      )
      .def_prop_rw(
          "radius_x",
          &Wall3dVertexStruct::radius_x,
          &Wall3dVertexStruct::set_radius_x,
          "Radius of arc or ellipse x-axis half width. 0 => Straight line."
      )
      .def_prop_rw(
          "radius_y",
          &Wall3dVertexStruct::radius_y,
          &Wall3dVertexStruct::set_radius_y,
          "Ellipse y-axis half height."
      )
      .def_prop_rw(
          "tilt",
          &Wall3dVertexStruct::tilt,
          &Wall3dVertexStruct::set_tilt,
          "Tilt of ellipse"
      )
      .def_prop_rw(
          "angle",
          &Wall3dVertexStruct::angle,
          &Wall3dVertexStruct::set_angle,
          "Angle of (x, y) point."
      )
      .def_prop_rw("x0", &Wall3dVertexStruct::x0, &Wall3dVertexStruct::set_x0, "Center of ellipse")
      .def_prop_rw("y0", &Wall3dVertexStruct::y0, &Wall3dVertexStruct::set_y0, "Center of ellipse")
      .def_prop_rw(
          "type",
          &Wall3dVertexStruct::type,
          &Wall3dVertexStruct::set_type,
          "No longer used."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return Wall3dVertexStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Wall3dVertexStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const Wall3dVertexStruct &self, nb::dict &memo) { return Wall3dVertexStruct(self); }
      )
      .def(
          "__eq__",
          [](const Wall3dVertexStruct &self, const Wall3dVertexStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Wall3dVertexStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
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