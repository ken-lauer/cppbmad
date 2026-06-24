#include "pybmad/generated/structs_a.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// ac_kicker_freq_struct
void init_ac_kicker_freq_struct(py::module &m, py::class_<AcKickerFreqStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         py::arg("f") = py::none(),
         py::arg("amp") = py::none(),
         py::arg("phi") = py::none()
  )
      .def_property("f", &AcKickerFreqStruct::f, &AcKickerFreqStruct::set_f)
      .def_property("amp", &AcKickerFreqStruct::amp, &AcKickerFreqStruct::set_amp)
      .def_property("phi", &AcKickerFreqStruct::phi, &AcKickerFreqStruct::set_phi)
      .def_static(
          "new_array1d",
          [](int sz) { return AcKickerFreqStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AcKickerFreqStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const AcKickerFreqStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerFreqStruct &self) {
            return AcKickerFreqStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AcKickerFreqStruct &self, py::dict &memo) { return AcKickerFreqStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerFreqStruct &self, const AcKickerFreqStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerFreqStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<AcKickerFreqStructArray1D, AcKickerFreqStructAlloc1D>(
      m,
      "AcKickerFreqStructArray1D",
      "AcKickerFreqStructAlloc1D"
  );
  // 2D AcKickerFreqStruct arrays are not used in structs/routines
  // 3D AcKickerFreqStruct arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_struct
void init_ac_kicker_struct(py::module &m, py::class_<AcKickerStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly(
          "amp_vs_time",
          py::cpp_function(&AcKickerStruct::amp_vs_time, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "frequency",
          py::cpp_function(&AcKickerStruct::frequency, py::keep_alive<0, 1>())
      )

      .def("__repr__", [](const AcKickerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerStruct &self) {
            return AcKickerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AcKickerStruct &self, py::dict &memo) { return AcKickerStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerStruct &self, const AcKickerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D AcKickerStruct arrays are not used in structs/routines
  // 2D AcKickerStruct arrays are not used in structs/routines
  // 3D AcKickerStruct arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_time_struct
void init_ac_kicker_time_struct(py::module &m, py::class_<AcKickerTimeStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, optional_ref<const SplineStruct>>(),
         py::arg("amp") = py::none(),
         py::arg("time") = py::none(),
         py::arg("spline") = py::none()
  )
      .def_property("amp", &AcKickerTimeStruct::amp, &AcKickerTimeStruct::set_amp)
      .def_property("time", &AcKickerTimeStruct::time, &AcKickerTimeStruct::set_time)
      .def_property(
          "spline",
          py::cpp_function(&AcKickerTimeStruct::spline, py::keep_alive<0, 1>()),
          &AcKickerTimeStruct::set_spline
      )
      .def_static(
          "new_array1d",
          [](int sz) { return AcKickerTimeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AcKickerTimeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const AcKickerTimeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerTimeStruct &self) {
            return AcKickerTimeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AcKickerTimeStruct &self, py::dict &memo) { return AcKickerTimeStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerTimeStruct &self, const AcKickerTimeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerTimeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<AcKickerTimeStructArray1D, AcKickerTimeStructAlloc1D>(
      m,
      "AcKickerTimeStructArray1D",
      "AcKickerTimeStructAlloc1D"
  );
  // 2D AcKickerTimeStruct arrays are not used in structs/routines
  // 3D AcKickerTimeStruct arrays are not used in structs/routines
}

// =============================================================================
// anormal_mode_struct
void init_anormal_mode_struct(py::module &m, py::class_<AnormalModeStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("emittance") = py::none(),
         py::arg("emittance_no_vert") = py::none(),
         py::arg("synch_int") = py::none(),
         py::arg("j_damp") = py::none(),
         py::arg("alpha_damp") = py::none(),
         py::arg("chrom") = py::none(),
         py::arg("tune") = py::none()
  )
      .def_property(
          "emittance",
          &AnormalModeStruct::emittance,
          &AnormalModeStruct::set_emittance,
          "Beam emittance (unnormalized). Includes vertical photon opening angle."
      )
      .def_property(
          "emittance_no_vert",
          &AnormalModeStruct::emittance_no_vert,
          &AnormalModeStruct::set_emittance_no_vert,
          "Unnormalized beam emittance without the vertical photon opening angle taken into "
          "account."
      )
      .def_property(
          "synch_int",
          py::cpp_function(&AnormalModeStruct::synch_int, py::keep_alive<0, 1>()),
          &AnormalModeStruct::set_synch_int,
          "Synchrotron integrals"
      )
      .def_property(
          "j_damp",
          &AnormalModeStruct::j_damp,
          &AnormalModeStruct::set_j_damp,
          "damping partition number"
      )
      .def_property(
          "alpha_damp",
          &AnormalModeStruct::alpha_damp,
          &AnormalModeStruct::set_alpha_damp,
          "damping per turn"
      )
      .def_property(
          "chrom",
          &AnormalModeStruct::chrom,
          &AnormalModeStruct::set_chrom,
          "Chromaticity"
      )
      .def_property(
          "tune",
          &AnormalModeStruct::tune,
          &AnormalModeStruct::set_tune,
          "'Fractional' tune in radians"
      )

      .def("__repr__", [](const AnormalModeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AnormalModeStruct &self) {
            return AnormalModeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AnormalModeStruct &self, py::dict &memo) { return AnormalModeStruct(self); }
      )
      .def(
          "__eq__",
          [](const AnormalModeStruct &self, const AnormalModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AnormalModeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D AnormalModeStruct arrays are not used in structs/routines
  // 2D AnormalModeStruct arrays are not used in structs/routines
  // 3D AnormalModeStruct arrays are not used in structs/routines
}

// =============================================================================
// aperture_param_struct
void init_aperture_param_struct(py::module &m, py::class_<ApertureParamStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>>(),
         py::arg("min_angle") = py::none(),
         py::arg("max_angle") = py::none(),
         py::arg("n_angle") = py::none(),
         py::arg("n_turn") = py::none(),
         py::arg("x_init") = py::none(),
         py::arg("y_init") = py::none(),
         py::arg("rel_accuracy") = py::none(),
         py::arg("abs_accuracy") = py::none(),
         py::arg("start_ele") = py::none()
  )
      .def_property(
          "min_angle",
          &ApertureParamStruct::min_angle,
          &ApertureParamStruct::set_min_angle
      )
      .def_property(
          "max_angle",
          &ApertureParamStruct::max_angle,
          &ApertureParamStruct::set_max_angle
      )
      .def_property("n_angle", &ApertureParamStruct::n_angle, &ApertureParamStruct::set_n_angle)
      .def_property(
          "n_turn",
          &ApertureParamStruct::n_turn,
          &ApertureParamStruct::set_n_turn,
          "Number of turns a particle must survive."
      )
      .def_property(
          "x_init",
          &ApertureParamStruct::x_init,
          &ApertureParamStruct::set_x_init,
          "Initial x coordinate to start with for theta_xy = 0."
      )
      .def_property(
          "y_init",
          &ApertureParamStruct::y_init,
          &ApertureParamStruct::set_y_init,
          "Initial y coordinate to start with for theta_xy = pi/2."
      )
      .def_property(
          "rel_accuracy",
          &ApertureParamStruct::rel_accuracy,
          &ApertureParamStruct::set_rel_accuracy,
          "Relative resolution of bracketed aperture."
      )
      .def_property(
          "abs_accuracy",
          &ApertureParamStruct::abs_accuracy,
          &ApertureParamStruct::set_abs_accuracy,
          "Absolute resolution of bracketed aperture (meters)."
      )
      .def_property(
          "start_ele",
          &ApertureParamStruct::start_ele,
          &ApertureParamStruct::set_start_ele,
          "Element to start tracking at."
      )

      .def("__repr__", [](const ApertureParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureParamStruct &self) {
            return ApertureParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ApertureParamStruct &self, py::dict &memo) { return ApertureParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const ApertureParamStruct &self, const ApertureParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ApertureParamStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D ApertureParamStruct arrays are not used in structs/routines
  // 2D ApertureParamStruct arrays are not used in structs/routines
  // 3D ApertureParamStruct arrays are not used in structs/routines
}

// =============================================================================
// aperture_point_struct
void init_aperture_point_struct(py::module &m, py::class_<AperturePointStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("plane") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("i_turn") = py::none()
  )
      .def_property(
          "x",
          &AperturePointStruct::x,
          &AperturePointStruct::set_x,
          "(x,y) aperture point with respect to the reference orbit."
      )
      .def_property(
          "y",
          &AperturePointStruct::y,
          &AperturePointStruct::set_y,
          "(x,y) aperture point with respect to the reference orbit."
      )
      .def_property(
          "plane",
          &AperturePointStruct::plane,
          &AperturePointStruct::set_plane,
          "plane determining loss"
      )
      .def_property(
          "ix_ele",
          &AperturePointStruct::ix_ele,
          &AperturePointStruct::set_ix_ele,
          "ele index particle lost at"
      )
      .def_property(
          "i_turn",
          &AperturePointStruct::i_turn,
          &AperturePointStruct::set_i_turn,
          "turn particle lost at"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return AperturePointStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AperturePointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const AperturePointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AperturePointStruct &self) {
            return AperturePointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AperturePointStruct &self, py::dict &memo) { return AperturePointStruct(self); }
      )
      .def(
          "__eq__",
          [](const AperturePointStruct &self, const AperturePointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AperturePointStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<AperturePointStructArray1D, AperturePointStructAlloc1D>(
      m,
      "AperturePointStructArray1D",
      "AperturePointStructAlloc1D"
  );
  // 2D AperturePointStruct arrays are not used in structs/routines
  // 3D AperturePointStruct arrays are not used in structs/routines
}

// =============================================================================
// aperture_scan_struct
void init_aperture_scan_struct(py::module &m, py::class_<ApertureScanStruct> &cls) {
  cls.def(
         py::init<optional_ref<const CoordStruct>, std::optional<double>>(),
         py::arg("ref_orb") = py::none(),
         py::arg("pz_start") = py::none()
  )
      .def_property_readonly(
          "point",
          py::cpp_function(&ApertureScanStruct::point, py::keep_alive<0, 1>()),
          "Set of aperture points at different angles."
      )
      .def_property(
          "ref_orb",
          py::cpp_function(&ApertureScanStruct::ref_orb, py::keep_alive<0, 1>()),
          &ApertureScanStruct::set_ref_orb,
          "Ref orbit around which the scan is made."
      )
      .def_property(
          "pz_start",
          &ApertureScanStruct::pz_start,
          &ApertureScanStruct::set_pz_start,
          "Starting pz."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ApertureScanStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ApertureScanStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
      )

      .def("__repr__", [](const ApertureScanStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ApertureScanStruct &self) {
            return ApertureScanStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ApertureScanStruct &self, py::dict &memo) { return ApertureScanStruct(self); }
      )
      .def(
          "__eq__",
          [](const ApertureScanStruct &self, const ApertureScanStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ApertureScanStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<ApertureScanStructArray1D, ApertureScanStructAlloc1D>(
      m,
      "ApertureScanStructArray1D",
      "ApertureScanStructAlloc1D"
  );
  // 2D ApertureScanStruct arrays are not used in structs/routines
  // 3D ApertureScanStruct arrays are not used in structs/routines
}

// =============================================================================
// all_encompassing_struct
void init_all_encompassing_struct(py::module &m, py::class_<AllEncompassingStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<std::complex<double>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::vector<std::complex<double>>>>,
             std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>,
             std::optional<std::complex<double>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::vector<std::complex<double>>>>,
             std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>,
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::vector<std::vector<std::complex<double>>>>,
             std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>,
             std::optional<int>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<std::vector<int>>>,
             std::optional<std::vector<std::vector<std::vector<int>>>>,
             std::optional<int>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<std::vector<int>>>,
             std::optional<std::vector<std::vector<std::vector<int>>>>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<std::vector<int>>>,
             std::optional<std::vector<std::vector<std::vector<int>>>>,
             std::optional<int64_t>,
             std::optional<std::vector<int64_t>>,
             std::optional<std::vector<std::vector<int64_t>>>,
             std::optional<std::vector<std::vector<std::vector<int64_t>>>>,
             std::optional<int64_t>,
             std::optional<std::vector<int64_t>>,
             std::optional<std::vector<std::vector<int64_t>>>,
             std::optional<std::vector<std::vector<std::vector<int64_t>>>>,
             std::optional<std::vector<int64_t>>,
             std::optional<std::vector<std::vector<int64_t>>>,
             std::optional<std::vector<std::vector<std::vector<int64_t>>>>,
             std::optional<bool>,
             std::optional<std::vector<bool>>,
             std::optional<std::vector<std::vector<bool>>>,
             std::optional<std::vector<std::vector<std::vector<bool>>>>,
             std::optional<bool>,
             std::optional<std::vector<bool>>,
             optional_ref<const TestSubStruct>,
             optional_ref<const TestSubStruct>>(),
         py::arg("real_rp_0d") = py::none(),
         py::arg("real_rp_1d") = py::none(),
         py::arg("real_rp_2d") = py::none(),
         py::arg("real_rp_3d") = py::none(),
         py::arg("real_rp_0d_ptr") = py::none(),
         py::arg("real_rp_1d_ptr") = py::none(),
         py::arg("real_rp_2d_ptr") = py::none(),
         py::arg("real_rp_3d_ptr") = py::none(),
         py::arg("real_rp_1d_alloc") = py::none(),
         py::arg("real_rp_2d_alloc") = py::none(),
         py::arg("real_rp_3d_alloc") = py::none(),
         py::arg("real_dp_0d") = py::none(),
         py::arg("real_dp_1d") = py::none(),
         py::arg("real_dp_2d") = py::none(),
         py::arg("real_dp_3d") = py::none(),
         py::arg("real_dp_0d_ptr") = py::none(),
         py::arg("real_dp_1d_ptr") = py::none(),
         py::arg("real_dp_2d_ptr") = py::none(),
         py::arg("real_dp_3d_ptr") = py::none(),
         py::arg("real_dp_1d_alloc") = py::none(),
         py::arg("real_dp_2d_alloc") = py::none(),
         py::arg("real_dp_3d_alloc") = py::none(),
         py::arg("complex_dp_0d") = py::none(),
         py::arg("complex_dp_1d") = py::none(),
         py::arg("complex_dp_2d") = py::none(),
         py::arg("complex_dp_3d") = py::none(),
         py::arg("complex_dp_0d_ptr") = py::none(),
         py::arg("complex_dp_1d_ptr") = py::none(),
         py::arg("complex_dp_2d_ptr") = py::none(),
         py::arg("complex_dp_3d_ptr") = py::none(),
         py::arg("complex_dp_1d_alloc") = py::none(),
         py::arg("complex_dp_2d_alloc") = py::none(),
         py::arg("complex_dp_3d_alloc") = py::none(),
         py::arg("int_0d") = py::none(),
         py::arg("int_1d") = py::none(),
         py::arg("int_2d") = py::none(),
         py::arg("int_3d") = py::none(),
         py::arg("int_0d_ptr") = py::none(),
         py::arg("int_1d_ptr") = py::none(),
         py::arg("int_2d_ptr") = py::none(),
         py::arg("int_3d_ptr") = py::none(),
         py::arg("int_1d_alloc") = py::none(),
         py::arg("int_2d_alloc") = py::none(),
         py::arg("int_3d_alloc") = py::none(),
         py::arg("int8_0d") = py::none(),
         py::arg("int8_1d") = py::none(),
         py::arg("int8_2d") = py::none(),
         py::arg("int8_3d") = py::none(),
         py::arg("int8_0d_ptr") = py::none(),
         py::arg("int8_1d_ptr") = py::none(),
         py::arg("int8_2d_ptr") = py::none(),
         py::arg("int8_3d_ptr") = py::none(),
         py::arg("int8_1d_alloc") = py::none(),
         py::arg("int8_2d_alloc") = py::none(),
         py::arg("int8_3d_alloc") = py::none(),
         py::arg("logical_0d") = py::none(),
         py::arg("logical_1d") = py::none(),
         py::arg("logical_2d") = py::none(),
         py::arg("logical_3d") = py::none(),
         py::arg("logical_0d_ptr") = py::none(),
         py::arg("logical_1d_ptr") = py::none(),
         py::arg("type_0d") = py::none(),
         py::arg("type_0d_ptr") = py::none()
  )
      .def_property(
          "real_rp_0d",
          &AllEncompassingStruct::real_rp_0d,
          &AllEncompassingStruct::set_real_rp_0d
      )
      .def_property(
          "real_rp_1d",
          py::cpp_function(&AllEncompassingStruct::real_rp_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_1d
      )
      .def_property(
          "real_rp_2d",
          py::cpp_function(&AllEncompassingStruct::real_rp_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_2d
      )
      .def_property(
          "real_rp_3d",
          py::cpp_function(&AllEncompassingStruct::real_rp_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_3d
      )
      .def_property(
          "real_rp_0d_ptr",
          &AllEncompassingStruct::real_rp_0d_ptr,
          &AllEncompassingStruct::set_real_rp_0d_ptr
      )
      .def_property(
          "real_rp_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_rp_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_1d_ptr
      )
      .def_property(
          "real_rp_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_rp_2d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_2d_ptr
      )
      .def_property(
          "real_rp_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_rp_3d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_3d_ptr
      )
      .def_property(
          "real_rp_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_rp_1d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_1d_alloc
      )
      .def_property(
          "real_rp_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_rp_2d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_2d_alloc
      )
      .def_property(
          "real_rp_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_rp_3d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_rp_3d_alloc,
          "Real(dp)"
      )
      .def_property(
          "real_dp_0d",
          &AllEncompassingStruct::real_dp_0d,
          &AllEncompassingStruct::set_real_dp_0d
      )
      .def_property(
          "real_dp_1d",
          py::cpp_function(&AllEncompassingStruct::real_dp_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_1d
      )
      .def_property(
          "real_dp_2d",
          py::cpp_function(&AllEncompassingStruct::real_dp_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_2d
      )
      .def_property(
          "real_dp_3d",
          py::cpp_function(&AllEncompassingStruct::real_dp_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_3d
      )
      .def_property(
          "real_dp_0d_ptr",
          &AllEncompassingStruct::real_dp_0d_ptr,
          &AllEncompassingStruct::set_real_dp_0d_ptr
      )
      .def_property(
          "real_dp_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_dp_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_1d_ptr
      )
      .def_property(
          "real_dp_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_dp_2d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_2d_ptr
      )
      .def_property(
          "real_dp_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::real_dp_3d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_3d_ptr
      )
      .def_property(
          "real_dp_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_dp_1d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_1d_alloc
      )
      .def_property(
          "real_dp_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_dp_2d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_2d_alloc
      )
      .def_property(
          "real_dp_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::real_dp_3d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_real_dp_3d_alloc,
          "complex(dp)"
      )
      .def_property(
          "complex_dp_0d",
          &AllEncompassingStruct::complex_dp_0d,
          &AllEncompassingStruct::set_complex_dp_0d
      )
      .def_property(
          "complex_dp_1d",
          py::cpp_function(&AllEncompassingStruct::complex_dp_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_1d
      )
      .def_property(
          "complex_dp_2d",
          py::cpp_function(&AllEncompassingStruct::complex_dp_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_2d
      )
      .def_property(
          "complex_dp_3d",
          py::cpp_function(&AllEncompassingStruct::complex_dp_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_3d
      )
      .def_property(
          "complex_dp_0d_ptr",
          &AllEncompassingStruct::complex_dp_0d_ptr,
          &AllEncompassingStruct::set_complex_dp_0d_ptr
      )
      .def_property(
          "complex_dp_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::complex_dp_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_1d_ptr
      )
      .def_property(
          "complex_dp_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::complex_dp_2d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_2d_ptr
      )
      .def_property(
          "complex_dp_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::complex_dp_3d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_3d_ptr
      )
      .def_property(
          "complex_dp_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::complex_dp_1d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_1d_alloc
      )
      .def_property(
          "complex_dp_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::complex_dp_2d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_2d_alloc
      )
      .def_property(
          "complex_dp_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::complex_dp_3d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_complex_dp_3d_alloc,
          "Integer"
      )
      .def_property("int_0d", &AllEncompassingStruct::int_0d, &AllEncompassingStruct::set_int_0d)
      .def_property(
          "int_1d",
          py::cpp_function(&AllEncompassingStruct::int_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_1d
      )
      .def_property(
          "int_2d",
          py::cpp_function(&AllEncompassingStruct::int_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_2d
      )
      .def_property(
          "int_3d",
          py::cpp_function(&AllEncompassingStruct::int_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_3d
      )
      .def_property(
          "int_0d_ptr",
          &AllEncompassingStruct::int_0d_ptr,
          &AllEncompassingStruct::set_int_0d_ptr
      )
      .def_property(
          "int_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::int_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_1d_ptr
      )
      .def_property(
          "int_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::int_2d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_2d_ptr
      )
      .def_property(
          "int_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::int_3d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_3d_ptr
      )
      .def_property(
          "int_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::int_1d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_1d_alloc
      )
      .def_property(
          "int_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::int_2d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_2d_alloc
      )
      .def_property(
          "int_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::int_3d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int_3d_alloc,
          "Integer8"
      )
      .def_property("int8_0d", &AllEncompassingStruct::int8_0d, &AllEncompassingStruct::set_int8_0d)
      .def_property(
          "int8_1d",
          py::cpp_function(&AllEncompassingStruct::int8_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_1d
      )
      .def_property(
          "int8_2d",
          py::cpp_function(&AllEncompassingStruct::int8_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_2d
      )
      .def_property(
          "int8_3d",
          py::cpp_function(&AllEncompassingStruct::int8_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_3d
      )
      .def_property(
          "int8_0d_ptr",
          &AllEncompassingStruct::int8_0d_ptr,
          &AllEncompassingStruct::set_int8_0d_ptr
      )
      .def_property(
          "int8_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::int8_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_1d_ptr
      )
      .def_property(
          "int8_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::int8_2d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_2d_ptr
      )
      .def_property(
          "int8_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::int8_3d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_3d_ptr
      )
      .def_property(
          "int8_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::int8_1d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_1d_alloc
      )
      .def_property(
          "int8_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::int8_2d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_2d_alloc
      )
      .def_property(
          "int8_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::int8_3d_alloc, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_int8_3d_alloc,
          "logical"
      )
      .def_property(
          "logical_0d",
          &AllEncompassingStruct::logical_0d,
          &AllEncompassingStruct::set_logical_0d
      )
      .def_property(
          "logical_1d",
          py::cpp_function(&AllEncompassingStruct::logical_1d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_logical_1d
      )
      .def_property(
          "logical_2d",
          py::cpp_function(&AllEncompassingStruct::logical_2d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_logical_2d
      )
      .def_property(
          "logical_3d",
          py::cpp_function(&AllEncompassingStruct::logical_3d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_logical_3d
      )
      .def_property(
          "logical_0d_ptr",
          &AllEncompassingStruct::logical_0d_ptr,
          &AllEncompassingStruct::set_logical_0d_ptr
      )
      .def_property(
          "logical_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::logical_1d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_logical_1d_ptr,
          "logical, pointer :: logical_2d_ptr(:,:) logical, pointer :: logical_3d_ptr(:,:,:) "
          "logical, allocatable :: logical_1d_alloc(:) logical, allocatable :: "
          "logical_2d_alloc(:,:) logical, allocatable :: logical_3d_alloc(:,:,:) type"
      )
      .def_property(
          "type_0d",
          py::cpp_function(&AllEncompassingStruct::type_0d, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_type_0d
      )
      .def_property_readonly(
          "type_1d",
          py::cpp_function(&AllEncompassingStruct::type_1d, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_2d",
          py::cpp_function(&AllEncompassingStruct::type_2d, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_3d",
          py::cpp_function(&AllEncompassingStruct::type_3d, py::keep_alive<0, 1>())
      )
      .def_property(
          "type_0d_ptr",
          py::cpp_function(&AllEncompassingStruct::type_0d_ptr, py::keep_alive<0, 1>()),
          &AllEncompassingStruct::set_type_0d_ptr
      )
      .def_property_readonly(
          "type_1d_ptr",
          py::cpp_function(&AllEncompassingStruct::type_1d_ptr, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_2d_ptr",
          py::cpp_function(&AllEncompassingStruct::type_2d_ptr, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_3d_ptr",
          py::cpp_function(&AllEncompassingStruct::type_3d_ptr, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_1d_alloc",
          py::cpp_function(&AllEncompassingStruct::type_1d_alloc, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_2d_alloc",
          py::cpp_function(&AllEncompassingStruct::type_2d_alloc, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "type_3d_alloc",
          py::cpp_function(&AllEncompassingStruct::type_3d_alloc, py::keep_alive<0, 1>())
      )

      .def("__repr__", [](const AllEncompassingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AllEncompassingStruct &self) {
            return AllEncompassingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AllEncompassingStruct &self, py::dict &memo) {
            return AllEncompassingStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const AllEncompassingStruct &self, const AllEncompassingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const AllEncompassingStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D AllEncompassingStruct arrays are not used in structs/routines
  // 2D AllEncompassingStruct arrays are not used in structs/routines
  // 3D AllEncompassingStruct arrays are not used in structs/routines
}