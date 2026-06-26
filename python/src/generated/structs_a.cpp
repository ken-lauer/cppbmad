#include "pybmad/generated/structs_a.hpp"

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
// astra_lattice_param_struct
void init_astra_lattice_param_struct(nb::module_ &m, nb::class_<AstraLatticeParamStruct> &cls) {
  cls.def(nb::init<std::optional<int>>(), nb::arg("fieldmap_dimension") = nb::none())
      .def_prop_rw(
          "fieldmap_dimension",
          &AstraLatticeParamStruct::fieldmap_dimension,
          &AstraLatticeParamStruct::set_fieldmap_dimension,
          "Dimensions for field map. 1 or 3"
      )

      .def("__repr__", [](const AstraLatticeParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AstraLatticeParamStruct &self) {
            return AstraLatticeParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AstraLatticeParamStruct &self, nb::dict &memo) {
            return AstraLatticeParamStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const AstraLatticeParamStruct &self, const AstraLatticeParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AstraLatticeParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D AstraLatticeParamStruct arrays are not used in structs/routines
  // 2D AstraLatticeParamStruct arrays are not used in structs/routines
  // 3D AstraLatticeParamStruct arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_freq_struct
void init_ac_kicker_freq_struct(nb::module_ &m, nb::class_<AcKickerFreqStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("f") = nb::none(),
         nb::arg("amp") = nb::none(),
         nb::arg("phi") = nb::none()
  )
      .def_prop_rw("f", &AcKickerFreqStruct::f, &AcKickerFreqStruct::set_f)
      .def_prop_rw("amp", &AcKickerFreqStruct::amp, &AcKickerFreqStruct::set_amp)
      .def_prop_rw("phi", &AcKickerFreqStruct::phi, &AcKickerFreqStruct::set_phi)
      .def_static(
          "new_array1d",
          [](int sz) { return AcKickerFreqStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AcKickerFreqStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const AcKickerFreqStruct &self, nb::dict &memo) { return AcKickerFreqStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerFreqStruct &self, const AcKickerFreqStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerFreqStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_ac_kicker_struct(nb::module_ &m, nb::class_<AcKickerStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("amp_vs_time", &AcKickerStruct::amp_vs_time, nb::keep_alive<0, 1>())
      .def_prop_ro("frequency", &AcKickerStruct::frequency, nb::keep_alive<0, 1>())

      .def("__repr__", [](const AcKickerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AcKickerStruct &self) {
            return AcKickerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AcKickerStruct &self, nb::dict &memo) { return AcKickerStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerStruct &self, const AcKickerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D AcKickerStruct arrays are not used in structs/routines
  // 2D AcKickerStruct arrays are not used in structs/routines
  // 3D AcKickerStruct arrays are not used in structs/routines
}

// =============================================================================
// ac_kicker_time_struct
void init_ac_kicker_time_struct(nb::module_ &m, nb::class_<AcKickerTimeStruct> &cls) {
  cls.def(
         "__init__",
         [](AcKickerTimeStruct *self,
            std::optional<double> amp,
            std::optional<double> time,
            const SplineStruct *spline) {
           new (self) AcKickerTimeStruct(amp, time, ptr_to_opt_ref(spline));
         },
         nb::arg("amp") = nb::none(),
         nb::arg("time") = nb::none(),
         nb::arg("spline") = nb::none()
  )
      .def_prop_rw("amp", &AcKickerTimeStruct::amp, &AcKickerTimeStruct::set_amp)
      .def_prop_rw("time", &AcKickerTimeStruct::time, &AcKickerTimeStruct::set_time)
      .def_prop_rw(
          "spline",
          &AcKickerTimeStruct::spline,
          &AcKickerTimeStruct::set_spline,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return AcKickerTimeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AcKickerTimeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const AcKickerTimeStruct &self, nb::dict &memo) { return AcKickerTimeStruct(self); }
      )
      .def(
          "__eq__",
          [](const AcKickerTimeStruct &self, const AcKickerTimeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AcKickerTimeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_anormal_mode_struct(nb::module_ &m, nb::class_<AnormalModeStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("emittance") = nb::none(),
         nb::arg("emittance_no_vert") = nb::none(),
         nb::arg("synch_int") = nb::none(),
         nb::arg("j_damp") = nb::none(),
         nb::arg("alpha_damp") = nb::none(),
         nb::arg("chrom") = nb::none(),
         nb::arg("tune") = nb::none()
  )
      .def_prop_rw(
          "emittance",
          &AnormalModeStruct::emittance,
          &AnormalModeStruct::set_emittance,
          "Beam emittance (unnormalized). Includes vertical photon opening angle."
      )
      .def_prop_rw(
          "emittance_no_vert",
          &AnormalModeStruct::emittance_no_vert,
          &AnormalModeStruct::set_emittance_no_vert,
          "Unnormalized beam emittance without the vertical photon opening angle taken into "
          "account."
      )
      .def_prop_rw(
          "synch_int",
          &AnormalModeStruct::synch_int,
          &AnormalModeStruct::set_synch_int,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Synchrotron integrals"
      )
      .def_prop_rw(
          "j_damp",
          &AnormalModeStruct::j_damp,
          &AnormalModeStruct::set_j_damp,
          "damping partition number"
      )
      .def_prop_rw(
          "alpha_damp",
          &AnormalModeStruct::alpha_damp,
          &AnormalModeStruct::set_alpha_damp,
          "damping per turn"
      )
      .def_prop_rw(
          "chrom",
          &AnormalModeStruct::chrom,
          &AnormalModeStruct::set_chrom,
          "Chromaticity"
      )
      .def_prop_rw(
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
          [](const AnormalModeStruct &self, nb::dict &memo) { return AnormalModeStruct(self); }
      )
      .def(
          "__eq__",
          [](const AnormalModeStruct &self, const AnormalModeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AnormalModeStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D AnormalModeStruct arrays are not used in structs/routines
  // 2D AnormalModeStruct arrays are not used in structs/routines
  // 3D AnormalModeStruct arrays are not used in structs/routines
}

// =============================================================================
// aperture_param_struct
void init_aperture_param_struct(nb::module_ &m, nb::class_<ApertureParamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>>(),
         nb::arg("min_angle") = nb::none(),
         nb::arg("max_angle") = nb::none(),
         nb::arg("n_angle") = nb::none(),
         nb::arg("n_turn") = nb::none(),
         nb::arg("x_init") = nb::none(),
         nb::arg("y_init") = nb::none(),
         nb::arg("rel_accuracy") = nb::none(),
         nb::arg("abs_accuracy") = nb::none(),
         nb::arg("start_ele") = nb::none()
  )
      .def_prop_rw(
          "min_angle",
          &ApertureParamStruct::min_angle,
          &ApertureParamStruct::set_min_angle
      )
      .def_prop_rw(
          "max_angle",
          &ApertureParamStruct::max_angle,
          &ApertureParamStruct::set_max_angle
      )
      .def_prop_rw("n_angle", &ApertureParamStruct::n_angle, &ApertureParamStruct::set_n_angle)
      .def_prop_rw(
          "n_turn",
          &ApertureParamStruct::n_turn,
          &ApertureParamStruct::set_n_turn,
          "Number of turns a particle must survive."
      )
      .def_prop_rw(
          "x_init",
          &ApertureParamStruct::x_init,
          &ApertureParamStruct::set_x_init,
          "Initial x coordinate to start with for theta_xy = 0."
      )
      .def_prop_rw(
          "y_init",
          &ApertureParamStruct::y_init,
          &ApertureParamStruct::set_y_init,
          "Initial y coordinate to start with for theta_xy = pi/2."
      )
      .def_prop_rw(
          "rel_accuracy",
          &ApertureParamStruct::rel_accuracy,
          &ApertureParamStruct::set_rel_accuracy,
          "Relative resolution of bracketed aperture."
      )
      .def_prop_rw(
          "abs_accuracy",
          &ApertureParamStruct::abs_accuracy,
          &ApertureParamStruct::set_abs_accuracy,
          "Absolute resolution of bracketed aperture (meters)."
      )
      .def_prop_rw(
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
          [](const ApertureParamStruct &self, nb::dict &memo) { return ApertureParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const ApertureParamStruct &self, const ApertureParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ApertureParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ApertureParamStruct arrays are not used in structs/routines
  // 2D ApertureParamStruct arrays are not used in structs/routines
  // 3D ApertureParamStruct arrays are not used in structs/routines
}

// =============================================================================
// aperture_point_struct
void init_aperture_point_struct(nb::module_ &m, nb::class_<AperturePointStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("plane") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("i_turn") = nb::none()
  )
      .def_prop_rw(
          "x",
          &AperturePointStruct::x,
          &AperturePointStruct::set_x,
          "(x,y) aperture point with respect to the reference orbit."
      )
      .def_prop_rw(
          "y",
          &AperturePointStruct::y,
          &AperturePointStruct::set_y,
          "(x,y) aperture point with respect to the reference orbit."
      )
      .def_prop_rw(
          "plane",
          &AperturePointStruct::plane,
          &AperturePointStruct::set_plane,
          "plane determining loss"
      )
      .def_prop_rw(
          "ix_ele",
          &AperturePointStruct::ix_ele,
          &AperturePointStruct::set_ix_ele,
          "ele index particle lost at"
      )
      .def_prop_rw(
          "i_turn",
          &AperturePointStruct::i_turn,
          &AperturePointStruct::set_i_turn,
          "turn particle lost at"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return AperturePointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AperturePointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const AperturePointStruct &self, nb::dict &memo) { return AperturePointStruct(self); }
      )
      .def(
          "__eq__",
          [](const AperturePointStruct &self, const AperturePointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AperturePointStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_aperture_scan_struct(nb::module_ &m, nb::class_<ApertureScanStruct> &cls) {
  cls.def(
         "__init__",
         [](ApertureScanStruct *self, const CoordStruct *ref_orb, std::optional<double> pz_start) {
           new (self) ApertureScanStruct(ptr_to_opt_ref(ref_orb), pz_start);
         },
         nb::arg("ref_orb") = nb::none(),
         nb::arg("pz_start") = nb::none()
  )
      .def_prop_ro(
          "point",
          &ApertureScanStruct::point,
          nb::keep_alive<0, 1>(),
          "Set of aperture points at different angles."
      )
      .def_prop_rw(
          "ref_orb",
          &ApertureScanStruct::ref_orb,
          &ApertureScanStruct::set_ref_orb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Ref orbit around which the scan is made."
      )
      .def_prop_rw(
          "pz_start",
          &ApertureScanStruct::pz_start,
          &ApertureScanStruct::set_pz_start,
          "Starting pz."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ApertureScanStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ApertureScanStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const ApertureScanStruct &self, nb::dict &memo) { return ApertureScanStruct(self); }
      )
      .def(
          "__eq__",
          [](const ApertureScanStruct &self, const ApertureScanStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ApertureScanStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
// all_pointer_struct
void init_all_pointer_struct(nb::module_ &m, nb::class_<AllPointerStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<long double>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<int>>>(),
         nb::arg("r") = nb::none(),
         nb::arg("q") = nb::none(),
         nb::arg("i") = nb::none(),
         nb::arg("l") = nb::none(),
         nb::arg("r1") = nb::none(),
         nb::arg("i1") = nb::none()
  )
      .def_prop_rw("r", &AllPointerStruct::r, &AllPointerStruct::set_r)
      .def_prop_rw("q", &AllPointerStruct::q, &AllPointerStruct::set_q)
      .def_prop_rw("i", &AllPointerStruct::i, &AllPointerStruct::set_i)
      .def_prop_rw("l", &AllPointerStruct::l, &AllPointerStruct::set_l)
      .def_prop_rw(
          "r1",
          &AllPointerStruct::r1,
          &AllPointerStruct::set_r1,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "i1",
          &AllPointerStruct::i1,
          &AllPointerStruct::set_i1,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return AllPointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = AllPointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const AllPointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AllPointerStruct &self) {
            return AllPointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AllPointerStruct &self, nb::dict &memo) { return AllPointerStruct(self); }
      )
      .def(
          "__eq__",
          [](const AllPointerStruct &self, const AllPointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AllPointerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<AllPointerStructArray1D, AllPointerStructAlloc1D>(
      m,
      "AllPointerStructArray1D",
      "AllPointerStructAlloc1D"
  );
  // 2D AllPointerStruct arrays are not used in structs/routines
  // 3D AllPointerStruct arrays are not used in structs/routines
}

// =============================================================================
// all_encompassing_struct
void init_all_encompassing_struct(nb::module_ &m, nb::class_<AllEncompassingStruct> &cls) {
  cls.def(
         "__init__",
         [](AllEncompassingStruct *self,
            std::optional<double> real_rp_0d,
            std::optional<std::vector<double>> real_rp_1d,
            std::optional<std::vector<std::vector<double>>> real_rp_2d,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_rp_3d,
            std::optional<double> real_rp_0d_ptr,
            std::optional<std::vector<double>> real_rp_1d_ptr,
            std::optional<std::vector<std::vector<double>>> real_rp_2d_ptr,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_rp_3d_ptr,
            std::optional<std::vector<double>> real_rp_1d_alloc,
            std::optional<std::vector<std::vector<double>>> real_rp_2d_alloc,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_rp_3d_alloc,
            std::optional<double> real_dp_0d,
            std::optional<std::vector<double>> real_dp_1d,
            std::optional<std::vector<std::vector<double>>> real_dp_2d,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_dp_3d,
            std::optional<double> real_dp_0d_ptr,
            std::optional<std::vector<double>> real_dp_1d_ptr,
            std::optional<std::vector<std::vector<double>>> real_dp_2d_ptr,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_dp_3d_ptr,
            std::optional<std::vector<double>> real_dp_1d_alloc,
            std::optional<std::vector<std::vector<double>>> real_dp_2d_alloc,
            std::optional<std::vector<std::vector<std::vector<double>>>> real_dp_3d_alloc,
            std::optional<std::complex<double>> complex_dp_0d,
            std::optional<std::vector<std::complex<double>>> complex_dp_1d,
            std::optional<std::vector<std::vector<std::complex<double>>>> complex_dp_2d,
            std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>
                complex_dp_3d,
            std::optional<std::complex<double>> complex_dp_0d_ptr,
            std::optional<std::vector<std::complex<double>>> complex_dp_1d_ptr,
            std::optional<std::vector<std::vector<std::complex<double>>>> complex_dp_2d_ptr,
            std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>
                complex_dp_3d_ptr,
            std::optional<std::vector<std::complex<double>>> complex_dp_1d_alloc,
            std::optional<std::vector<std::vector<std::complex<double>>>> complex_dp_2d_alloc,
            std::optional<std::vector<std::vector<std::vector<std::complex<double>>>>>
                complex_dp_3d_alloc,
            std::optional<int> int_0d,
            std::optional<std::vector<int>> int_1d,
            std::optional<std::vector<std::vector<int>>> int_2d,
            std::optional<std::vector<std::vector<std::vector<int>>>> int_3d,
            std::optional<int> int_0d_ptr,
            std::optional<std::vector<int>> int_1d_ptr,
            std::optional<std::vector<std::vector<int>>> int_2d_ptr,
            std::optional<std::vector<std::vector<std::vector<int>>>> int_3d_ptr,
            std::optional<std::vector<int>> int_1d_alloc,
            std::optional<std::vector<std::vector<int>>> int_2d_alloc,
            std::optional<std::vector<std::vector<std::vector<int>>>> int_3d_alloc,
            std::optional<int64_t> int8_0d,
            std::optional<std::vector<int64_t>> int8_1d,
            std::optional<std::vector<std::vector<int64_t>>> int8_2d,
            std::optional<std::vector<std::vector<std::vector<int64_t>>>> int8_3d,
            std::optional<int64_t> int8_0d_ptr,
            std::optional<std::vector<int64_t>> int8_1d_ptr,
            std::optional<std::vector<std::vector<int64_t>>> int8_2d_ptr,
            std::optional<std::vector<std::vector<std::vector<int64_t>>>> int8_3d_ptr,
            std::optional<std::vector<int64_t>> int8_1d_alloc,
            std::optional<std::vector<std::vector<int64_t>>> int8_2d_alloc,
            std::optional<std::vector<std::vector<std::vector<int64_t>>>> int8_3d_alloc,
            std::optional<bool> logical_0d,
            std::optional<std::vector<bool>> logical_1d,
            std::optional<std::vector<std::vector<bool>>> logical_2d,
            std::optional<std::vector<std::vector<std::vector<bool>>>> logical_3d,
            std::optional<bool> logical_0d_ptr,
            std::optional<std::vector<bool>> logical_1d_ptr,
            const TestSubStruct *type_0d,
            const TestSubStruct *type_0d_ptr) {
           new (self) AllEncompassingStruct(
               real_rp_0d,
               real_rp_1d,
               real_rp_2d,
               real_rp_3d,
               real_rp_0d_ptr,
               real_rp_1d_ptr,
               real_rp_2d_ptr,
               real_rp_3d_ptr,
               real_rp_1d_alloc,
               real_rp_2d_alloc,
               real_rp_3d_alloc,
               real_dp_0d,
               real_dp_1d,
               real_dp_2d,
               real_dp_3d,
               real_dp_0d_ptr,
               real_dp_1d_ptr,
               real_dp_2d_ptr,
               real_dp_3d_ptr,
               real_dp_1d_alloc,
               real_dp_2d_alloc,
               real_dp_3d_alloc,
               complex_dp_0d,
               complex_dp_1d,
               complex_dp_2d,
               complex_dp_3d,
               complex_dp_0d_ptr,
               complex_dp_1d_ptr,
               complex_dp_2d_ptr,
               complex_dp_3d_ptr,
               complex_dp_1d_alloc,
               complex_dp_2d_alloc,
               complex_dp_3d_alloc,
               int_0d,
               int_1d,
               int_2d,
               int_3d,
               int_0d_ptr,
               int_1d_ptr,
               int_2d_ptr,
               int_3d_ptr,
               int_1d_alloc,
               int_2d_alloc,
               int_3d_alloc,
               int8_0d,
               int8_1d,
               int8_2d,
               int8_3d,
               int8_0d_ptr,
               int8_1d_ptr,
               int8_2d_ptr,
               int8_3d_ptr,
               int8_1d_alloc,
               int8_2d_alloc,
               int8_3d_alloc,
               logical_0d,
               logical_1d,
               logical_2d,
               logical_3d,
               logical_0d_ptr,
               logical_1d_ptr,
               ptr_to_opt_ref(type_0d),
               ptr_to_opt_ref(type_0d_ptr)
           );
         },
         nb::arg("real_rp_0d") = nb::none(),
         nb::arg("real_rp_1d") = nb::none(),
         nb::arg("real_rp_2d") = nb::none(),
         nb::arg("real_rp_3d") = nb::none(),
         nb::arg("real_rp_0d_ptr") = nb::none(),
         nb::arg("real_rp_1d_ptr") = nb::none(),
         nb::arg("real_rp_2d_ptr") = nb::none(),
         nb::arg("real_rp_3d_ptr") = nb::none(),
         nb::arg("real_rp_1d_alloc") = nb::none(),
         nb::arg("real_rp_2d_alloc") = nb::none(),
         nb::arg("real_rp_3d_alloc") = nb::none(),
         nb::arg("real_dp_0d") = nb::none(),
         nb::arg("real_dp_1d") = nb::none(),
         nb::arg("real_dp_2d") = nb::none(),
         nb::arg("real_dp_3d") = nb::none(),
         nb::arg("real_dp_0d_ptr") = nb::none(),
         nb::arg("real_dp_1d_ptr") = nb::none(),
         nb::arg("real_dp_2d_ptr") = nb::none(),
         nb::arg("real_dp_3d_ptr") = nb::none(),
         nb::arg("real_dp_1d_alloc") = nb::none(),
         nb::arg("real_dp_2d_alloc") = nb::none(),
         nb::arg("real_dp_3d_alloc") = nb::none(),
         nb::arg("complex_dp_0d") = nb::none(),
         nb::arg("complex_dp_1d") = nb::none(),
         nb::arg("complex_dp_2d") = nb::none(),
         nb::arg("complex_dp_3d") = nb::none(),
         nb::arg("complex_dp_0d_ptr") = nb::none(),
         nb::arg("complex_dp_1d_ptr") = nb::none(),
         nb::arg("complex_dp_2d_ptr") = nb::none(),
         nb::arg("complex_dp_3d_ptr") = nb::none(),
         nb::arg("complex_dp_1d_alloc") = nb::none(),
         nb::arg("complex_dp_2d_alloc") = nb::none(),
         nb::arg("complex_dp_3d_alloc") = nb::none(),
         nb::arg("int_0d") = nb::none(),
         nb::arg("int_1d") = nb::none(),
         nb::arg("int_2d") = nb::none(),
         nb::arg("int_3d") = nb::none(),
         nb::arg("int_0d_ptr") = nb::none(),
         nb::arg("int_1d_ptr") = nb::none(),
         nb::arg("int_2d_ptr") = nb::none(),
         nb::arg("int_3d_ptr") = nb::none(),
         nb::arg("int_1d_alloc") = nb::none(),
         nb::arg("int_2d_alloc") = nb::none(),
         nb::arg("int_3d_alloc") = nb::none(),
         nb::arg("int8_0d") = nb::none(),
         nb::arg("int8_1d") = nb::none(),
         nb::arg("int8_2d") = nb::none(),
         nb::arg("int8_3d") = nb::none(),
         nb::arg("int8_0d_ptr") = nb::none(),
         nb::arg("int8_1d_ptr") = nb::none(),
         nb::arg("int8_2d_ptr") = nb::none(),
         nb::arg("int8_3d_ptr") = nb::none(),
         nb::arg("int8_1d_alloc") = nb::none(),
         nb::arg("int8_2d_alloc") = nb::none(),
         nb::arg("int8_3d_alloc") = nb::none(),
         nb::arg("logical_0d") = nb::none(),
         nb::arg("logical_1d") = nb::none(),
         nb::arg("logical_2d") = nb::none(),
         nb::arg("logical_3d") = nb::none(),
         nb::arg("logical_0d_ptr") = nb::none(),
         nb::arg("logical_1d_ptr") = nb::none(),
         nb::arg("type_0d") = nb::none(),
         nb::arg("type_0d_ptr") = nb::none()
  )
      .def_prop_rw(
          "real_rp_0d",
          &AllEncompassingStruct::real_rp_0d,
          &AllEncompassingStruct::set_real_rp_0d
      )
      .def_prop_rw(
          "real_rp_1d",
          &AllEncompassingStruct::real_rp_1d,
          &AllEncompassingStruct::set_real_rp_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_2d",
          &AllEncompassingStruct::real_rp_2d,
          &AllEncompassingStruct::set_real_rp_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_3d",
          &AllEncompassingStruct::real_rp_3d,
          &AllEncompassingStruct::set_real_rp_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_0d_ptr",
          &AllEncompassingStruct::real_rp_0d_ptr,
          &AllEncompassingStruct::set_real_rp_0d_ptr
      )
      .def_prop_rw(
          "real_rp_1d_ptr",
          &AllEncompassingStruct::real_rp_1d_ptr,
          &AllEncompassingStruct::set_real_rp_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_2d_ptr",
          &AllEncompassingStruct::real_rp_2d_ptr,
          &AllEncompassingStruct::set_real_rp_2d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_3d_ptr",
          &AllEncompassingStruct::real_rp_3d_ptr,
          &AllEncompassingStruct::set_real_rp_3d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_1d_alloc",
          &AllEncompassingStruct::real_rp_1d_alloc,
          &AllEncompassingStruct::set_real_rp_1d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_2d_alloc",
          &AllEncompassingStruct::real_rp_2d_alloc,
          &AllEncompassingStruct::set_real_rp_2d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_rp_3d_alloc",
          &AllEncompassingStruct::real_rp_3d_alloc,
          &AllEncompassingStruct::set_real_rp_3d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Real(dp)"
      )
      .def_prop_rw(
          "real_dp_0d",
          &AllEncompassingStruct::real_dp_0d,
          &AllEncompassingStruct::set_real_dp_0d
      )
      .def_prop_rw(
          "real_dp_1d",
          &AllEncompassingStruct::real_dp_1d,
          &AllEncompassingStruct::set_real_dp_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_2d",
          &AllEncompassingStruct::real_dp_2d,
          &AllEncompassingStruct::set_real_dp_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_3d",
          &AllEncompassingStruct::real_dp_3d,
          &AllEncompassingStruct::set_real_dp_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_0d_ptr",
          &AllEncompassingStruct::real_dp_0d_ptr,
          &AllEncompassingStruct::set_real_dp_0d_ptr
      )
      .def_prop_rw(
          "real_dp_1d_ptr",
          &AllEncompassingStruct::real_dp_1d_ptr,
          &AllEncompassingStruct::set_real_dp_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_2d_ptr",
          &AllEncompassingStruct::real_dp_2d_ptr,
          &AllEncompassingStruct::set_real_dp_2d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_3d_ptr",
          &AllEncompassingStruct::real_dp_3d_ptr,
          &AllEncompassingStruct::set_real_dp_3d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_1d_alloc",
          &AllEncompassingStruct::real_dp_1d_alloc,
          &AllEncompassingStruct::set_real_dp_1d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_2d_alloc",
          &AllEncompassingStruct::real_dp_2d_alloc,
          &AllEncompassingStruct::set_real_dp_2d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "real_dp_3d_alloc",
          &AllEncompassingStruct::real_dp_3d_alloc,
          &AllEncompassingStruct::set_real_dp_3d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "complex(dp)"
      )
      .def_prop_rw(
          "complex_dp_0d",
          &AllEncompassingStruct::complex_dp_0d,
          &AllEncompassingStruct::set_complex_dp_0d
      )
      .def_prop_rw(
          "complex_dp_1d",
          &AllEncompassingStruct::complex_dp_1d,
          &AllEncompassingStruct::set_complex_dp_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_2d",
          &AllEncompassingStruct::complex_dp_2d,
          &AllEncompassingStruct::set_complex_dp_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_3d",
          &AllEncompassingStruct::complex_dp_3d,
          &AllEncompassingStruct::set_complex_dp_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_0d_ptr",
          &AllEncompassingStruct::complex_dp_0d_ptr,
          &AllEncompassingStruct::set_complex_dp_0d_ptr
      )
      .def_prop_rw(
          "complex_dp_1d_ptr",
          &AllEncompassingStruct::complex_dp_1d_ptr,
          &AllEncompassingStruct::set_complex_dp_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_2d_ptr",
          &AllEncompassingStruct::complex_dp_2d_ptr,
          &AllEncompassingStruct::set_complex_dp_2d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_3d_ptr",
          &AllEncompassingStruct::complex_dp_3d_ptr,
          &AllEncompassingStruct::set_complex_dp_3d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_1d_alloc",
          &AllEncompassingStruct::complex_dp_1d_alloc,
          &AllEncompassingStruct::set_complex_dp_1d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_2d_alloc",
          &AllEncompassingStruct::complex_dp_2d_alloc,
          &AllEncompassingStruct::set_complex_dp_2d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "complex_dp_3d_alloc",
          &AllEncompassingStruct::complex_dp_3d_alloc,
          &AllEncompassingStruct::set_complex_dp_3d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Integer"
      )
      .def_prop_rw("int_0d", &AllEncompassingStruct::int_0d, &AllEncompassingStruct::set_int_0d)
      .def_prop_rw(
          "int_1d",
          &AllEncompassingStruct::int_1d,
          &AllEncompassingStruct::set_int_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_2d",
          &AllEncompassingStruct::int_2d,
          &AllEncompassingStruct::set_int_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_3d",
          &AllEncompassingStruct::int_3d,
          &AllEncompassingStruct::set_int_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_0d_ptr",
          &AllEncompassingStruct::int_0d_ptr,
          &AllEncompassingStruct::set_int_0d_ptr
      )
      .def_prop_rw(
          "int_1d_ptr",
          &AllEncompassingStruct::int_1d_ptr,
          &AllEncompassingStruct::set_int_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_2d_ptr",
          &AllEncompassingStruct::int_2d_ptr,
          &AllEncompassingStruct::set_int_2d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_3d_ptr",
          &AllEncompassingStruct::int_3d_ptr,
          &AllEncompassingStruct::set_int_3d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_1d_alloc",
          &AllEncompassingStruct::int_1d_alloc,
          &AllEncompassingStruct::set_int_1d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_2d_alloc",
          &AllEncompassingStruct::int_2d_alloc,
          &AllEncompassingStruct::set_int_2d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int_3d_alloc",
          &AllEncompassingStruct::int_3d_alloc,
          &AllEncompassingStruct::set_int_3d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Integer8"
      )
      .def_prop_rw("int8_0d", &AllEncompassingStruct::int8_0d, &AllEncompassingStruct::set_int8_0d)
      .def_prop_rw(
          "int8_1d",
          &AllEncompassingStruct::int8_1d,
          &AllEncompassingStruct::set_int8_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_2d",
          &AllEncompassingStruct::int8_2d,
          &AllEncompassingStruct::set_int8_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_3d",
          &AllEncompassingStruct::int8_3d,
          &AllEncompassingStruct::set_int8_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_0d_ptr",
          &AllEncompassingStruct::int8_0d_ptr,
          &AllEncompassingStruct::set_int8_0d_ptr
      )
      .def_prop_rw(
          "int8_1d_ptr",
          &AllEncompassingStruct::int8_1d_ptr,
          &AllEncompassingStruct::set_int8_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_2d_ptr",
          &AllEncompassingStruct::int8_2d_ptr,
          &AllEncompassingStruct::set_int8_2d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_3d_ptr",
          &AllEncompassingStruct::int8_3d_ptr,
          &AllEncompassingStruct::set_int8_3d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_1d_alloc",
          &AllEncompassingStruct::int8_1d_alloc,
          &AllEncompassingStruct::set_int8_1d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_2d_alloc",
          &AllEncompassingStruct::int8_2d_alloc,
          &AllEncompassingStruct::set_int8_2d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "int8_3d_alloc",
          &AllEncompassingStruct::int8_3d_alloc,
          &AllEncompassingStruct::set_int8_3d_alloc,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "logical"
      )
      .def_prop_rw(
          "logical_0d",
          &AllEncompassingStruct::logical_0d,
          &AllEncompassingStruct::set_logical_0d
      )
      .def_prop_rw(
          "logical_1d",
          &AllEncompassingStruct::logical_1d,
          &AllEncompassingStruct::set_logical_1d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "logical_2d",
          &AllEncompassingStruct::logical_2d,
          &AllEncompassingStruct::set_logical_2d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "logical_3d",
          &AllEncompassingStruct::logical_3d,
          &AllEncompassingStruct::set_logical_3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "logical_0d_ptr",
          &AllEncompassingStruct::logical_0d_ptr,
          &AllEncompassingStruct::set_logical_0d_ptr
      )
      .def_prop_rw(
          "logical_1d_ptr",
          &AllEncompassingStruct::logical_1d_ptr,
          &AllEncompassingStruct::set_logical_1d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "logical, pointer :: logical_2d_ptr(:,:) logical, pointer :: logical_3d_ptr(:,:,:) "
          "logical, allocatable :: logical_1d_alloc(:) logical, allocatable :: "
          "logical_2d_alloc(:,:) logical, allocatable :: logical_3d_alloc(:,:,:) type"
      )
      .def_prop_rw(
          "type_0d",
          &AllEncompassingStruct::type_0d,
          &AllEncompassingStruct::set_type_0d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("type_1d", &AllEncompassingStruct::type_1d, nb::keep_alive<0, 1>())
      .def_prop_ro("type_2d", &AllEncompassingStruct::type_2d, nb::keep_alive<0, 1>())
      .def_prop_ro("type_3d", &AllEncompassingStruct::type_3d, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "type_0d_ptr",
          &AllEncompassingStruct::type_0d_ptr,
          &AllEncompassingStruct::set_type_0d_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("type_1d_ptr", &AllEncompassingStruct::type_1d_ptr, nb::keep_alive<0, 1>())
      .def_prop_ro("type_2d_ptr", &AllEncompassingStruct::type_2d_ptr, nb::keep_alive<0, 1>())
      .def_prop_ro("type_3d_ptr", &AllEncompassingStruct::type_3d_ptr, nb::keep_alive<0, 1>())
      .def_prop_ro("type_1d_alloc", &AllEncompassingStruct::type_1d_alloc, nb::keep_alive<0, 1>())
      .def_prop_ro("type_2d_alloc", &AllEncompassingStruct::type_2d_alloc, nb::keep_alive<0, 1>())
      .def_prop_ro("type_3d_alloc", &AllEncompassingStruct::type_3d_alloc, nb::keep_alive<0, 1>())

      .def("__repr__", [](const AllEncompassingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const AllEncompassingStruct &self) {
            return AllEncompassingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const AllEncompassingStruct &self, nb::dict &memo) {
            return AllEncompassingStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const AllEncompassingStruct &self, const AllEncompassingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const AllEncompassingStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D AllEncompassingStruct arrays are not used in structs/routines
  // 2D AllEncompassingStruct arrays are not used in structs/routines
  // 3D AllEncompassingStruct arrays are not used in structs/routines
}