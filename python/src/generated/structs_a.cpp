#include "pybmad/generated/structs_a.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// ac_kicker_freq_struct
void init_ac_kicker_freq_struct(
    py::module& m,
    py::class_<AcKickerFreqProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerFreqProxy.f (0D_NOT_real -
      .def_property("f", &AcKickerFreqProxy::f, &AcKickerFreqProxy::set_f)
      // AcKickerFreqProxy.amp (0D_NOT_real -
      .def_property("amp", &AcKickerFreqProxy::amp, &AcKickerFreqProxy::set_amp)
      // AcKickerFreqProxy.phi (0D_NOT_real -
      .def_property("phi", &AcKickerFreqProxy::phi, &AcKickerFreqProxy::set_phi)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AcKickerFreqProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AcKickerFreqProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerFreqProxy& self) {
            return AcKickerFreqProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerFreqProxy& self, py::dict& memo) {
            return AcKickerFreqProxy(self);
          })

      ;

  bind_FTypeArrayND<AcKickerFreqProxyArray1D>(m, "AcKickerFreqStructArray1D");
  bind_FTypeAlloc1D<AcKickerFreqProxyAlloc1D>(m, "AcKickerFreqStructAlloc1D");
  // 2D AcKickerFreqProxy arrays are not used in structs/routines
  // 3D AcKickerFreqProxy arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_struct
void init_ac_kicker_struct(py::module& m, py::class_<AcKickerProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerProxy.amp_vs_time (1D_ALLOC_type -
      .def_property_readonly("amp_vs_time", &AcKickerProxy::amp_vs_time)
      // AcKickerProxy.frequency (1D_ALLOC_type -
      .def_property_readonly("frequency", &AcKickerProxy::frequency)

      .def(
          "__repr__", [](const AcKickerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerProxy& self) {
            return AcKickerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerProxy& self, py::dict& memo) {
            return AcKickerProxy(self);
          })

      ;

  // 1D AcKickerProxy arrays are not used in structs/routines
  // 2D AcKickerProxy arrays are not used in structs/routines
  // 3D AcKickerProxy arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_time_struct
void init_ac_kicker_time_struct(
    py::module& m,
    py::class_<AcKickerTimeProxy>& cls) {
  cls.def(py::init<>())
      // AcKickerTimeProxy.amp (0D_NOT_real -
      .def_property("amp", &AcKickerTimeProxy::amp, &AcKickerTimeProxy::set_amp)
      // AcKickerTimeProxy.time (0D_NOT_real -
      .def_property(
          "time", &AcKickerTimeProxy::time, &AcKickerTimeProxy::set_time)
      // AcKickerTimeProxy.spline (0D_NOT_type -
      .def_property(
          "spline", &AcKickerTimeProxy::spline, &AcKickerTimeProxy::set_spline)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AcKickerTimeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AcKickerTimeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerTimeProxy& self) {
            return AcKickerTimeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AcKickerTimeProxy& self, py::dict& memo) {
            return AcKickerTimeProxy(self);
          })

      ;

  bind_FTypeArrayND<AcKickerTimeProxyArray1D>(m, "AcKickerTimeStructArray1D");
  bind_FTypeAlloc1D<AcKickerTimeProxyAlloc1D>(m, "AcKickerTimeStructAlloc1D");
  // 2D AcKickerTimeProxy arrays are not used in structs/routines
  // 3D AcKickerTimeProxy arrays are not used in structs/routines
}

// =============================================================================
// anormal_mode_struct
void init_anormal_mode_struct(
    py::module& m,
    py::class_<AnormalModeProxy>& cls) {
  cls.def(py::init<>())
      // AnormalModeProxy.emittance (0D_NOT_real - Beam emittance (unnormalized). Includes vertical photon opening angle.
      .def_property(
          "emittance",
          &AnormalModeProxy::emittance,
          &AnormalModeProxy::set_emittance)
      // AnormalModeProxy.emittance_no_vert (0D_NOT_real - Unnormalized beam emittance without the vertical photon opening angle taken into account.
      .def_property(
          "emittance_no_vert",
          &AnormalModeProxy::emittance_no_vert,
          &AnormalModeProxy::set_emittance_no_vert)
      // AnormalModeProxy.synch_int (1D_NOT_real - Synchrotron integrals
      .def_property_readonly("synch_int", &AnormalModeProxy::synch_int)
      // AnormalModeProxy.j_damp (0D_NOT_real - damping partition number
      .def_property(
          "j_damp", &AnormalModeProxy::j_damp, &AnormalModeProxy::set_j_damp)
      // AnormalModeProxy.alpha_damp (0D_NOT_real - damping per turn
      .def_property(
          "alpha_damp",
          &AnormalModeProxy::alpha_damp,
          &AnormalModeProxy::set_alpha_damp)
      // AnormalModeProxy.chrom (0D_NOT_real - Chromaticity
      .def_property(
          "chrom", &AnormalModeProxy::chrom, &AnormalModeProxy::set_chrom)
      // AnormalModeProxy.tune (0D_NOT_real - 'Fractional' tune in radians
      .def_property(
          "tune", &AnormalModeProxy::tune, &AnormalModeProxy::set_tune)

      .def(
          "__repr__",
          [](const AnormalModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AnormalModeProxy& self) {
            return AnormalModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AnormalModeProxy& self, py::dict& memo) {
            return AnormalModeProxy(self);
          })

      ;

  // 1D AnormalModeProxy arrays are not used in structs/routines
  // 2D AnormalModeProxy arrays are not used in structs/routines
  // 3D AnormalModeProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_param_struct
void init_aperture_param_struct(
    py::module& m,
    py::class_<ApertureParamProxy>& cls) {
  cls.def(py::init<>())
      // ApertureParamProxy.min_angle (0D_NOT_real -
      .def_property(
          "min_angle",
          &ApertureParamProxy::min_angle,
          &ApertureParamProxy::set_min_angle)
      // ApertureParamProxy.max_angle (0D_NOT_real -
      .def_property(
          "max_angle",
          &ApertureParamProxy::max_angle,
          &ApertureParamProxy::set_max_angle)
      // ApertureParamProxy.n_angle (0D_NOT_integer -
      .def_property(
          "n_angle",
          &ApertureParamProxy::n_angle,
          &ApertureParamProxy::set_n_angle)
      // ApertureParamProxy.n_turn (0D_NOT_integer - Number of turns a particle must survive.
      .def_property(
          "n_turn",
          &ApertureParamProxy::n_turn,
          &ApertureParamProxy::set_n_turn)
      // ApertureParamProxy.x_init (0D_NOT_real - Initial x coordinate to start with for theta_xy = 0.
      .def_property(
          "x_init",
          &ApertureParamProxy::x_init,
          &ApertureParamProxy::set_x_init)
      // ApertureParamProxy.y_init (0D_NOT_real - Initial y coordinate to start with for theta_xy = pi/2.
      .def_property(
          "y_init",
          &ApertureParamProxy::y_init,
          &ApertureParamProxy::set_y_init)
      // ApertureParamProxy.rel_accuracy (0D_NOT_real - Relative resolution of bracketed aperture.
      .def_property(
          "rel_accuracy",
          &ApertureParamProxy::rel_accuracy,
          &ApertureParamProxy::set_rel_accuracy)
      // ApertureParamProxy.abs_accuracy (0D_NOT_real - Absolute resolution of bracketed aperture (meters).
      .def_property(
          "abs_accuracy",
          &ApertureParamProxy::abs_accuracy,
          &ApertureParamProxy::set_abs_accuracy)
      // ApertureParamProxy.start_ele (0D_NOT_character - Element to start tracking at.
      .def_property(
          "start_ele",
          &ApertureParamProxy::start_ele,
          &ApertureParamProxy::set_start_ele)

      .def(
          "__repr__",
          [](const ApertureParamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureParamProxy& self) {
            return ApertureParamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ApertureParamProxy& self, py::dict& memo) {
            return ApertureParamProxy(self);
          })

      ;

  // 1D ApertureParamProxy arrays are not used in structs/routines
  // 2D ApertureParamProxy arrays are not used in structs/routines
  // 3D ApertureParamProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_point_struct
void init_aperture_point_struct(
    py::module& m,
    py::class_<AperturePointProxy>& cls) {
  cls.def(py::init<>())
      // AperturePointProxy.x (0D_NOT_real - (x,y) aperture point with respect to the reference orbit.
      .def_property("x", &AperturePointProxy::x, &AperturePointProxy::set_x)
      // AperturePointProxy.y (0D_NOT_real - (x,y) aperture point with respect to the reference orbit.
      .def_property("y", &AperturePointProxy::y, &AperturePointProxy::set_y)
      // AperturePointProxy.plane (0D_NOT_integer - plane determining loss
      .def_property(
          "plane", &AperturePointProxy::plane, &AperturePointProxy::set_plane)
      // AperturePointProxy.ix_ele (0D_NOT_integer - ele index particle lost at
      .def_property(
          "ix_ele",
          &AperturePointProxy::ix_ele,
          &AperturePointProxy::set_ix_ele)
      // AperturePointProxy.i_turn (0D_NOT_integer - turn particle lost at
      .def_property(
          "i_turn",
          &AperturePointProxy::i_turn,
          &AperturePointProxy::set_i_turn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return AperturePointProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const AperturePointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AperturePointProxy& self) {
            return AperturePointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AperturePointProxy& self, py::dict& memo) {
            return AperturePointProxy(self);
          })

      ;

  bind_FTypeArrayND<AperturePointProxyArray1D>(m, "AperturePointStructArray1D");
  bind_FTypeAlloc1D<AperturePointProxyAlloc1D>(m, "AperturePointStructAlloc1D");
  // 2D AperturePointProxy arrays are not used in structs/routines
  // 3D AperturePointProxy arrays are not used in structs/routines
}

// =============================================================================
// aperture_scan_struct
void init_aperture_scan_struct(
    py::module& m,
    py::class_<ApertureScanProxy>& cls) {
  cls.def(py::init<>())
      // ApertureScanProxy.point (1D_ALLOC_type - Set of aperture points at different angles.
      .def_property_readonly("point", &ApertureScanProxy::point)
      // ApertureScanProxy.ref_orb (0D_NOT_type - Ref orbit around which the scan is made.
      .def_property(
          "ref_orb",
          &ApertureScanProxy::ref_orb,
          &ApertureScanProxy::set_ref_orb)
      // ApertureScanProxy.pz_start (0D_NOT_real - Starting pz.
      .def_property(
          "pz_start",
          &ApertureScanProxy::pz_start,
          &ApertureScanProxy::set_pz_start)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ApertureScanProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ApertureScanProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureScanProxy& self) {
            return ApertureScanProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ApertureScanProxy& self, py::dict& memo) {
            return ApertureScanProxy(self);
          })

      ;

  bind_FTypeArrayND<ApertureScanProxyArray1D>(m, "ApertureScanStructArray1D");
  bind_FTypeAlloc1D<ApertureScanProxyAlloc1D>(m, "ApertureScanStructAlloc1D");
  // 2D ApertureScanProxy arrays are not used in structs/routines
  // 3D ApertureScanProxy arrays are not used in structs/routines
}

// =============================================================================
// all_encompassing_struct
void init_all_encompassing_struct(
    py::module& m,
    py::class_<AllEncompassingProxy>& cls) {
  cls.def(py::init<>())
      // AllEncompassingProxy.real_rp_0d (0D_NOT_real -
      .def_property(
          "real_rp_0d",
          &AllEncompassingProxy::real_rp_0d,
          &AllEncompassingProxy::set_real_rp_0d)
      // AllEncompassingProxy.real_rp_1d (1D_NOT_real -
      .def_property_readonly("real_rp_1d", &AllEncompassingProxy::real_rp_1d)
      // AllEncompassingProxy.real_rp_2d (2D_NOT_real -
      .def_property_readonly("real_rp_2d", &AllEncompassingProxy::real_rp_2d)
      // AllEncompassingProxy.real_rp_3d (3D_NOT_real -
      .def_property_readonly("real_rp_3d", &AllEncompassingProxy::real_rp_3d)
      // AllEncompassingProxy.real_rp_0d_ptr (0D_PTR_real -
      .def_property(
          "real_rp_0d_ptr",
          &AllEncompassingProxy::real_rp_0d_ptr,
          &AllEncompassingProxy::set_real_rp_0d_ptr)
      // AllEncompassingProxy.real_rp_1d_ptr (1D_PTR_real -
      .def_property_readonly(
          "real_rp_1d_ptr", &AllEncompassingProxy::real_rp_1d_ptr)
      // AllEncompassingProxy.real_rp_2d_ptr (2D_PTR_real -
      .def_property_readonly(
          "real_rp_2d_ptr", &AllEncompassingProxy::real_rp_2d_ptr)
      // AllEncompassingProxy.real_rp_3d_ptr (3D_PTR_real -
      .def_property_readonly(
          "real_rp_3d_ptr", &AllEncompassingProxy::real_rp_3d_ptr)
      // AllEncompassingProxy.real_rp_1d_alloc (1D_ALLOC_real -
      .def_property_readonly(
          "real_rp_1d_alloc", &AllEncompassingProxy::real_rp_1d_alloc)
      // AllEncompassingProxy.real_rp_2d_alloc (2D_ALLOC_real -
      .def_property_readonly(
          "real_rp_2d_alloc", &AllEncompassingProxy::real_rp_2d_alloc)
      // AllEncompassingProxy.real_rp_3d_alloc (3D_ALLOC_real - Real(dp)
      .def_property_readonly(
          "real_rp_3d_alloc", &AllEncompassingProxy::real_rp_3d_alloc)
      // AllEncompassingProxy.real_dp_0d (0D_NOT_real -
      .def_property(
          "real_dp_0d",
          &AllEncompassingProxy::real_dp_0d,
          &AllEncompassingProxy::set_real_dp_0d)
      // AllEncompassingProxy.real_dp_1d (1D_NOT_real -
      .def_property_readonly("real_dp_1d", &AllEncompassingProxy::real_dp_1d)
      // AllEncompassingProxy.real_dp_2d (2D_NOT_real -
      .def_property_readonly("real_dp_2d", &AllEncompassingProxy::real_dp_2d)
      // AllEncompassingProxy.real_dp_3d (3D_NOT_real -
      .def_property_readonly("real_dp_3d", &AllEncompassingProxy::real_dp_3d)
      // AllEncompassingProxy.real_dp_0d_ptr (0D_PTR_real -
      .def_property(
          "real_dp_0d_ptr",
          &AllEncompassingProxy::real_dp_0d_ptr,
          &AllEncompassingProxy::set_real_dp_0d_ptr)
      // AllEncompassingProxy.real_dp_1d_ptr (1D_PTR_real -
      .def_property_readonly(
          "real_dp_1d_ptr", &AllEncompassingProxy::real_dp_1d_ptr)
      // AllEncompassingProxy.real_dp_2d_ptr (2D_PTR_real -
      .def_property_readonly(
          "real_dp_2d_ptr", &AllEncompassingProxy::real_dp_2d_ptr)
      // AllEncompassingProxy.real_dp_3d_ptr (3D_PTR_real -
      .def_property_readonly(
          "real_dp_3d_ptr", &AllEncompassingProxy::real_dp_3d_ptr)
      // AllEncompassingProxy.real_dp_1d_alloc (1D_ALLOC_real -
      .def_property_readonly(
          "real_dp_1d_alloc", &AllEncompassingProxy::real_dp_1d_alloc)
      // AllEncompassingProxy.real_dp_2d_alloc (2D_ALLOC_real -
      .def_property_readonly(
          "real_dp_2d_alloc", &AllEncompassingProxy::real_dp_2d_alloc)
      // AllEncompassingProxy.real_dp_3d_alloc (3D_ALLOC_real - complex(dp)
      .def_property_readonly(
          "real_dp_3d_alloc", &AllEncompassingProxy::real_dp_3d_alloc)
      // AllEncompassingProxy.complex_dp_0d (0D_NOT_complex -
      .def_property(
          "complex_dp_0d",
          &AllEncompassingProxy::complex_dp_0d,
          &AllEncompassingProxy::set_complex_dp_0d)
      // AllEncompassingProxy.complex_dp_1d (1D_NOT_complex -
      .def_property_readonly(
          "complex_dp_1d", &AllEncompassingProxy::complex_dp_1d)
      // AllEncompassingProxy.complex_dp_2d (2D_NOT_complex -
      .def_property_readonly(
          "complex_dp_2d", &AllEncompassingProxy::complex_dp_2d)
      // AllEncompassingProxy.complex_dp_3d (3D_NOT_complex - TODO complex(dp), pointer :: complex_dp_0d_ptr
      .def_property_readonly(
          "complex_dp_3d", &AllEncompassingProxy::complex_dp_3d)
      // AllEncompassingProxy.complex_dp_1d_ptr (1D_PTR_complex -
      .def_property_readonly(
          "complex_dp_1d_ptr", &AllEncompassingProxy::complex_dp_1d_ptr)
      // AllEncompassingProxy.complex_dp_2d_ptr (2D_PTR_complex -
      .def_property_readonly(
          "complex_dp_2d_ptr", &AllEncompassingProxy::complex_dp_2d_ptr)
      // AllEncompassingProxy.complex_dp_3d_ptr (3D_PTR_complex -
      .def_property_readonly(
          "complex_dp_3d_ptr", &AllEncompassingProxy::complex_dp_3d_ptr)
      // AllEncompassingProxy.complex_dp_1d_alloc (1D_ALLOC_complex -
      .def_property_readonly(
          "complex_dp_1d_alloc", &AllEncompassingProxy::complex_dp_1d_alloc)
      // AllEncompassingProxy.complex_dp_2d_alloc (2D_ALLOC_complex -
      .def_property_readonly(
          "complex_dp_2d_alloc", &AllEncompassingProxy::complex_dp_2d_alloc)
      // AllEncompassingProxy.complex_dp_3d_alloc (3D_ALLOC_complex - Integer
      .def_property_readonly(
          "complex_dp_3d_alloc", &AllEncompassingProxy::complex_dp_3d_alloc)
      // AllEncompassingProxy.int_0d (0D_NOT_integer -
      .def_property(
          "int_0d",
          &AllEncompassingProxy::int_0d,
          &AllEncompassingProxy::set_int_0d)
      // AllEncompassingProxy.int_1d (1D_NOT_integer -
      .def_property_readonly("int_1d", &AllEncompassingProxy::int_1d)
      // AllEncompassingProxy.int_2d (2D_NOT_integer -
      .def_property_readonly("int_2d", &AllEncompassingProxy::int_2d)
      // AllEncompassingProxy.int_3d (3D_NOT_integer -
      .def_property_readonly("int_3d", &AllEncompassingProxy::int_3d)
      // AllEncompassingProxy.int_0d_ptr (0D_PTR_integer -
      .def_property(
          "int_0d_ptr",
          &AllEncompassingProxy::int_0d_ptr,
          &AllEncompassingProxy::set_int_0d_ptr)
      // AllEncompassingProxy.int_1d_ptr (1D_PTR_integer -
      .def_property_readonly("int_1d_ptr", &AllEncompassingProxy::int_1d_ptr)
      // AllEncompassingProxy.int_2d_ptr (2D_PTR_integer -
      .def_property_readonly("int_2d_ptr", &AllEncompassingProxy::int_2d_ptr)
      // AllEncompassingProxy.int_3d_ptr (3D_PTR_integer -
      .def_property_readonly("int_3d_ptr", &AllEncompassingProxy::int_3d_ptr)
      // AllEncompassingProxy.int_1d_alloc (1D_ALLOC_integer -
      .def_property_readonly(
          "int_1d_alloc", &AllEncompassingProxy::int_1d_alloc)
      // AllEncompassingProxy.int_2d_alloc (2D_ALLOC_integer -
      .def_property_readonly(
          "int_2d_alloc", &AllEncompassingProxy::int_2d_alloc)
      // AllEncompassingProxy.int_3d_alloc (3D_ALLOC_integer - Integer8
      .def_property_readonly(
          "int_3d_alloc", &AllEncompassingProxy::int_3d_alloc)
      // AllEncompassingProxy.int8_0d (0D_NOT_integer8 -
      .def_property(
          "int8_0d",
          &AllEncompassingProxy::int8_0d,
          &AllEncompassingProxy::set_int8_0d)
      // 1D_NOT_integer8 int8_1d proxy support missing
      // 2D_NOT_integer8 int8_2d proxy support missing
      // 3D_NOT_integer8 int8_3d proxy support missing
      // AllEncompassingProxy.int8_0d_ptr (0D_PTR_integer8 -
      .def_property(
          "int8_0d_ptr",
          &AllEncompassingProxy::int8_0d_ptr,
          &AllEncompassingProxy::set_int8_0d_ptr)
      // 1D_PTR_integer8 int8_1d_ptr proxy support missing
      // 2D_PTR_integer8 int8_2d_ptr proxy support missing
      // 3D_PTR_integer8 int8_3d_ptr proxy support missing
      // 1D_ALLOC_integer8 int8_1d_alloc proxy support missing
      // 2D_ALLOC_integer8 int8_2d_alloc proxy support missing
      // 3D_ALLOC_integer8 int8_3d_alloc proxy support missing
      // AllEncompassingProxy.logical_0d (0D_NOT_logical -
      .def_property(
          "logical_0d",
          &AllEncompassingProxy::logical_0d,
          &AllEncompassingProxy::set_logical_0d)
      // 1D_NOT_logical logical_1d proxy support missing
      // 2D_NOT_logical logical_2d proxy support missing
      // 3D_NOT_logical logical_3d proxy support missing
      // AllEncompassingProxy.logical_0d_ptr (0D_PTR_logical - logical, pointer :: logical_1d_ptr(:) logical, pointer :: logical_2d_ptr(:,:) logical, pointer :: logical_3d_ptr(:,:,:) logical, allocatable :: logical_1d_alloc(:) logical, allocatable :: logical_2d_alloc(:,:) logical, allocatable :: logical_3d_alloc(:,:,:) type
      .def_property(
          "logical_0d_ptr",
          &AllEncompassingProxy::logical_0d_ptr,
          &AllEncompassingProxy::set_logical_0d_ptr)
      // AllEncompassingProxy.type_0d (0D_NOT_type -
      .def_property(
          "type_0d",
          &AllEncompassingProxy::type_0d,
          &AllEncompassingProxy::set_type_0d)
      // AllEncompassingProxy.type_1d (1D_NOT_type -
      .def_property_readonly("type_1d", &AllEncompassingProxy::type_1d)
      // AllEncompassingProxy.type_2d (2D_NOT_type -
      .def_property_readonly("type_2d", &AllEncompassingProxy::type_2d)
      // AllEncompassingProxy.type_3d (3D_NOT_type -
      .def_property_readonly("type_3d", &AllEncompassingProxy::type_3d)
      // AllEncompassingProxy.type_0d_ptr (0D_PTR_type -
      .def_property(
          "type_0d_ptr",
          &AllEncompassingProxy::type_0d_ptr,
          &AllEncompassingProxy::set_type_0d_ptr)
      // AllEncompassingProxy.type_1d_ptr (1D_PTR_type -
      .def_property_readonly("type_1d_ptr", &AllEncompassingProxy::type_1d_ptr)
      // AllEncompassingProxy.type_2d_ptr (2D_PTR_type -
      .def_property_readonly("type_2d_ptr", &AllEncompassingProxy::type_2d_ptr)
      // AllEncompassingProxy.type_3d_ptr (3D_PTR_type -
      .def_property_readonly("type_3d_ptr", &AllEncompassingProxy::type_3d_ptr)
      // AllEncompassingProxy.type_1d_alloc (1D_ALLOC_type -
      .def_property_readonly(
          "type_1d_alloc", &AllEncompassingProxy::type_1d_alloc)
      // AllEncompassingProxy.type_2d_alloc (2D_ALLOC_type -
      .def_property_readonly(
          "type_2d_alloc", &AllEncompassingProxy::type_2d_alloc)
      // AllEncompassingProxy.type_3d_alloc (3D_ALLOC_type -
      .def_property_readonly(
          "type_3d_alloc", &AllEncompassingProxy::type_3d_alloc)

      .def(
          "__repr__",
          [](const AllEncompassingProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AllEncompassingProxy& self) {
            return AllEncompassingProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const AllEncompassingProxy& self, py::dict& memo) {
            return AllEncompassingProxy(self);
          })

      ;

  // 1D AllEncompassingProxy arrays are not used in structs/routines
  // 2D AllEncompassingProxy arrays are not used in structs/routines
  // 3D AllEncompassingProxy arrays are not used in structs/routines
}