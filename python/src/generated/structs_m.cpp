#include "pybmad/generated/structs_m.hpp"

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
// material_struct
void init_material_struct(nb::module_ &m, nb::class_<MaterialStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("species") = nb::none(),
         nb::arg("number") = nb::none(),
         nb::arg("density") = nb::none(),
         nb::arg("density_used") = nb::none(),
         nb::arg("area_density") = nb::none(),
         nb::arg("area_density_used") = nb::none(),
         nb::arg("radiation_length") = nb::none(),
         nb::arg("radiation_length_used") = nb::none()
  )
      .def_prop_rw("species", &MaterialStruct::species, &MaterialStruct::set_species)
      .def_prop_rw(
          "number",
          &MaterialStruct::number,
          &MaterialStruct::set_number,
          "Relative number"
      )
      .def_prop_rw("density", &MaterialStruct::density, &MaterialStruct::set_density)
      .def_prop_rw("density_used", &MaterialStruct::density_used, &MaterialStruct::set_density_used)
      .def_prop_rw("area_density", &MaterialStruct::area_density, &MaterialStruct::set_area_density)
      .def_prop_rw(
          "area_density_used",
          &MaterialStruct::area_density_used,
          &MaterialStruct::set_area_density_used
      )
      .def_prop_rw(
          "radiation_length",
          &MaterialStruct::radiation_length,
          &MaterialStruct::set_radiation_length
      )
      .def_prop_rw(
          "radiation_length_used",
          &MaterialStruct::radiation_length_used,
          &MaterialStruct::set_radiation_length_used
      )
      .def_static(
          "new_array1d",
          [](int sz) { return MaterialStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MaterialStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MaterialStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MaterialStruct &self) {
            return MaterialStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MaterialStruct &self, nb::dict &memo) { return MaterialStruct(self); }
      )
      .def(
          "__eq__",
          [](const MaterialStruct &self, const MaterialStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MaterialStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MaterialStructArray1D, MaterialStructAlloc1D>(
      m,
      "MaterialStructArray1D",
      "MaterialStructAlloc1D"
  );
  // 2D MaterialStruct arrays are not used in structs/routines
  // 3D MaterialStruct arrays are not used in structs/routines
}

// =============================================================================
// mode3_struct
void init_mode3_struct(nb::module_ &m, nb::class_<Mode3Struct> &cls) {
  cls.def(
         "__init__",
         [](Mode3Struct *self,
            std::optional<std::vector<std::vector<double>>> v,
            const TwissStruct *a,
            const TwissStruct *b,
            const TwissStruct *c,
            const TwissStruct *x,
            const TwissStruct *y) {
           new (self) Mode3Struct(
               v,
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(c),
               ptr_to_opt_ref(x),
               ptr_to_opt_ref(y)
           );
         },
         nb::arg("v") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("c") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none()
  )
      .def_prop_rw(
          "v",
          &Mode3Struct::v,
          &Mode3Struct::set_v,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "a",
          &Mode3Struct::a,
          &Mode3Struct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b",
          &Mode3Struct::b,
          &Mode3Struct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "c",
          &Mode3Struct::c,
          &Mode3Struct::set_c,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "x",
          &Mode3Struct::x,
          &Mode3Struct::set_x,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "y",
          &Mode3Struct::y,
          &Mode3Struct::set_y,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const Mode3Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mode3Struct &self) {
            return Mode3Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Mode3Struct &self, nb::dict &memo) { return Mode3Struct(self); }
      )
      .def(
          "__eq__",
          [](const Mode3Struct &self, const Mode3Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Mode3Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D Mode3Struct arrays are not used in structs/routines
  // 2D Mode3Struct arrays are not used in structs/routines
  // 3D Mode3Struct arrays are not used in structs/routines
}

// =============================================================================
// mode_info_struct
void init_mode_info_struct(nb::module_ &m, nb::class_<ModeInfoStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("stable") = nb::none(),
         nb::arg("tune") = nb::none(),
         nb::arg("emit") = nb::none(),
         nb::arg("chrom") = nb::none(),
         nb::arg("sigma") = nb::none(),
         nb::arg("sigmap") = nb::none()
  )
      .def_prop_rw(
          "stable",
          &ModeInfoStruct::stable,
          &ModeInfoStruct::set_stable,
          "Is the mode stable?"
      )
      .def_prop_rw(
          "tune",
          &ModeInfoStruct::tune,
          &ModeInfoStruct::set_tune,
          "'fractional' tune in radians"
      )
      .def_prop_rw(
          "emit",
          &ModeInfoStruct::emit,
          &ModeInfoStruct::set_emit,
          "Emittance (unnormalized)."
      )
      .def_prop_rw("chrom", &ModeInfoStruct::chrom, &ModeInfoStruct::set_chrom, "Chromaticity.")
      .def_prop_rw("sigma", &ModeInfoStruct::sigma, &ModeInfoStruct::set_sigma, "Beam size.")
      .def_prop_rw(
          "sigmap",
          &ModeInfoStruct::sigmap,
          &ModeInfoStruct::set_sigmap,
          "Beam divergence."
      )

      .def("__repr__", [](const ModeInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ModeInfoStruct &self) {
            return ModeInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ModeInfoStruct &self, nb::dict &memo) { return ModeInfoStruct(self); }
      )
      .def(
          "__eq__",
          [](const ModeInfoStruct &self, const ModeInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ModeInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ModeInfoStruct arrays are not used in structs/routines
  // 2D ModeInfoStruct arrays are not used in structs/routines
  // 3D ModeInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_all_info_struct
void init_multipass_all_info_struct(nb::module_ &m, nb::class_<MultipassAllInfoStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("lord", &MultipassAllInfoStruct::lord, nb::keep_alive<0, 1>(), "Array of lords")
      .def_prop_ro("branch", &MultipassAllInfoStruct::branch, nb::keep_alive<0, 1>())

      .def("__repr__", [](const MultipassAllInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassAllInfoStruct &self) {
            return MultipassAllInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassAllInfoStruct &self, nb::dict &memo) {
            return MultipassAllInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassAllInfoStruct &self, const MultipassAllInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassAllInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MultipassAllInfoStruct arrays are not used in structs/routines
  // 2D MultipassAllInfoStruct arrays are not used in structs/routines
  // 3D MultipassAllInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_branch_info_struct
void init_multipass_branch_info_struct(nb::module_ &m, nb::class_<MultipassBranchInfoStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("ele", &MultipassBranchInfoStruct::ele, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return MultipassBranchInfoStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MultipassBranchInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MultipassBranchInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassBranchInfoStruct &self) {
            return MultipassBranchInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassBranchInfoStruct &self, nb::dict &memo) {
            return MultipassBranchInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassBranchInfoStruct &self, const MultipassBranchInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassBranchInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MultipassBranchInfoStructArray1D, MultipassBranchInfoStructAlloc1D>(
      m,
      "MultipassBranchInfoStructArray1D",
      "MultipassBranchInfoStructAlloc1D"
  );
  // 2D MultipassBranchInfoStruct arrays are not used in structs/routines
  // 3D MultipassBranchInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_ele_info_struct
void init_multipass_ele_info_struct(nb::module_ &m, nb::class_<MultipassEleInfoStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<int>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<int>>>(),
         nb::arg("multipass") = nb::none(),
         nb::arg("ix_pass") = nb::none(),
         nb::arg("ix_lord") = nb::none(),
         nb::arg("ix_super") = nb::none()
  )
      .def_prop_rw(
          "multipass",
          &MultipassEleInfoStruct::multipass,
          &MultipassEleInfoStruct::set_multipass,
          "True if involved in multipass. False otherwise"
      )
      .def_prop_rw(
          "ix_pass",
          &MultipassEleInfoStruct::ix_pass,
          &MultipassEleInfoStruct::set_ix_pass,
          "Pass number"
      )
      .def_prop_rw(
          "ix_lord",
          &MultipassEleInfoStruct::ix_lord,
          &MultipassEleInfoStruct::set_ix_lord,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Pointers to lord(:) array"
      )
      .def_prop_rw(
          "ix_super",
          &MultipassEleInfoStruct::ix_super,
          &MultipassEleInfoStruct::set_ix_super,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Indexes to slave(ix_pass, super_slave%ix_ele) matrix"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return MultipassEleInfoStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MultipassEleInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MultipassEleInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassEleInfoStruct &self) {
            return MultipassEleInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassEleInfoStruct &self, nb::dict &memo) {
            return MultipassEleInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassEleInfoStruct &self, const MultipassEleInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassEleInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MultipassEleInfoStructArray1D, MultipassEleInfoStructAlloc1D>(
      m,
      "MultipassEleInfoStructArray1D",
      "MultipassEleInfoStructAlloc1D"
  );
  // 2D MultipassEleInfoStruct arrays are not used in structs/routines
  // 3D MultipassEleInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_lord_info_struct
void init_multipass_lord_info_struct(nb::module_ &m, nb::class_<MultipassLordInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](MultipassLordInfoStruct *self,
            const EleStruct *lord,
            std::optional<int> n_pass,
            std::optional<int> n_super_slave) {
           new (self) MultipassLordInfoStruct(ptr_to_opt_ref(lord), n_pass, n_super_slave);
         },
         nb::arg("lord") = nb::none(),
         nb::arg("n_pass") = nb::none(),
         nb::arg("n_super_slave") = nb::none()
  )
      .def_prop_rw(
          "lord",
          &MultipassLordInfoStruct::lord,
          &MultipassLordInfoStruct::set_lord,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Lord element"
      )
      .def_prop_rw(
          "n_pass",
          &MultipassLordInfoStruct::n_pass,
          &MultipassLordInfoStruct::set_n_pass,
          "Number of passes (= number of slaves)"
      )
      .def_prop_rw(
          "n_super_slave",
          &MultipassLordInfoStruct::n_super_slave,
          &MultipassLordInfoStruct::set_n_super_slave,
          "Number of super_slaves per super_lord."
      )
      .def_prop_ro(
          "super_lord",
          &MultipassLordInfoStruct::super_lord,
          nb::keep_alive<0, 1>(),
          "Super_lord list if they exist."
      )
      .def_prop_ro(
          "slave",
          &MultipassLordInfoStruct::slave,
          nb::keep_alive<0, 1>(),
          "Slaves list in tracking part."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return MultipassLordInfoStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MultipassLordInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MultipassLordInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassLordInfoStruct &self) {
            return MultipassLordInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassLordInfoStruct &self, nb::dict &memo) {
            return MultipassLordInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassLordInfoStruct &self, const MultipassLordInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassLordInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MultipassLordInfoStructArray1D, MultipassLordInfoStructAlloc1D>(
      m,
      "MultipassLordInfoStructArray1D",
      "MultipassLordInfoStructAlloc1D"
  );
  // 2D MultipassLordInfoStruct arrays are not used in structs/routines
  // 3D MultipassLordInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// multipole_cache_struct
void init_multipole_cache_struct(nb::module_ &m, nb::class_<MultipoleCacheStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>>(),
         nb::arg("a_pole_mag") = nb::none(),
         nb::arg("b_pole_mag") = nb::none(),
         nb::arg("a_kick_mag") = nb::none(),
         nb::arg("b_kick_mag") = nb::none(),
         nb::arg("ix_pole_mag_max") = nb::none(),
         nb::arg("ix_kick_mag_max") = nb::none(),
         nb::arg("mag_valid") = nb::none(),
         nb::arg("a_pole_elec") = nb::none(),
         nb::arg("b_pole_elec") = nb::none(),
         nb::arg("a_kick_elec") = nb::none(),
         nb::arg("b_kick_elec") = nb::none(),
         nb::arg("ix_pole_elec_max") = nb::none(),
         nb::arg("ix_kick_elec_max") = nb::none(),
         nb::arg("elec_valid") = nb::none()
  )
      .def_prop_rw(
          "a_pole_mag",
          &MultipoleCacheStruct::a_pole_mag,
          &MultipoleCacheStruct::set_a_pole_mag,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b_pole_mag",
          &MultipoleCacheStruct::b_pole_mag,
          &MultipoleCacheStruct::set_b_pole_mag,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "a_kick_mag",
          &MultipoleCacheStruct::a_kick_mag,
          &MultipoleCacheStruct::set_a_kick_mag,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b_kick_mag",
          &MultipoleCacheStruct::b_kick_mag,
          &MultipoleCacheStruct::set_b_kick_mag,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ix_pole_mag_max",
          &MultipoleCacheStruct::ix_pole_mag_max,
          &MultipoleCacheStruct::set_ix_pole_mag_max
      )
      .def_prop_rw(
          "ix_kick_mag_max",
          &MultipoleCacheStruct::ix_kick_mag_max,
          &MultipoleCacheStruct::set_ix_kick_mag_max
      )
      .def_prop_rw(
          "mag_valid",
          &MultipoleCacheStruct::mag_valid,
          &MultipoleCacheStruct::set_mag_valid,
          "From elseparator hkick and vkick."
      )
      .def_prop_rw(
          "a_pole_elec",
          &MultipoleCacheStruct::a_pole_elec,
          &MultipoleCacheStruct::set_a_pole_elec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b_pole_elec",
          &MultipoleCacheStruct::b_pole_elec,
          &MultipoleCacheStruct::set_b_pole_elec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "a_kick_elec",
          &MultipoleCacheStruct::a_kick_elec,
          &MultipoleCacheStruct::set_a_kick_elec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b_kick_elec",
          &MultipoleCacheStruct::b_kick_elec,
          &MultipoleCacheStruct::set_b_kick_elec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ix_pole_elec_max",
          &MultipoleCacheStruct::ix_pole_elec_max,
          &MultipoleCacheStruct::set_ix_pole_elec_max
      )
      .def_prop_rw(
          "ix_kick_elec_max",
          &MultipoleCacheStruct::ix_kick_elec_max,
          &MultipoleCacheStruct::set_ix_kick_elec_max
      )
      .def_prop_rw(
          "elec_valid",
          &MultipoleCacheStruct::elec_valid,
          &MultipoleCacheStruct::set_elec_valid
      )

      .def("__repr__", [](const MultipoleCacheStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipoleCacheStruct &self) {
            return MultipoleCacheStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipoleCacheStruct &self, nb::dict &memo) {
            return MultipoleCacheStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipoleCacheStruct &self, const MultipoleCacheStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipoleCacheStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MultipoleCacheStruct arrays are not used in structs/routines
  // 2D MultipoleCacheStruct arrays are not used in structs/routines
  // 3D MultipoleCacheStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_energy_struct
void init_mad_energy_struct(nb::module_ &m, nb::class_<MadEnergyStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("total") = nb::none(),
         nb::arg("beta") = nb::none(),
         nb::arg("gamma") = nb::none(),
         nb::arg("kinetic") = nb::none(),
         nb::arg("p0c") = nb::none(),
         nb::arg("particle") = nb::none()
  )
      .def_prop_rw("total", &MadEnergyStruct::total, &MadEnergyStruct::set_total)
      .def_prop_rw(
          "beta",
          &MadEnergyStruct::beta,
          &MadEnergyStruct::set_beta,
          "normalized velocity: v/c"
      )
      .def_prop_rw(
          "gamma",
          &MadEnergyStruct::gamma,
          &MadEnergyStruct::set_gamma,
          "relativistic factor: 1/sqrt(1-beta^2)"
      )
      .def_prop_rw(
          "kinetic",
          &MadEnergyStruct::kinetic,
          &MadEnergyStruct::set_kinetic,
          "kinetic energy"
      )
      .def_prop_rw("p0c", &MadEnergyStruct::p0c, &MadEnergyStruct::set_p0c, "particle momentum")
      .def_prop_rw(
          "particle",
          &MadEnergyStruct::particle,
          &MadEnergyStruct::set_particle,
          "particle species"
      )

      .def("__repr__", [](const MadEnergyStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadEnergyStruct &self) {
            return MadEnergyStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MadEnergyStruct &self, nb::dict &memo) { return MadEnergyStruct(self); }
      )
      .def(
          "__eq__",
          [](const MadEnergyStruct &self, const MadEnergyStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MadEnergyStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MadEnergyStruct arrays are not used in structs/routines
  // 2D MadEnergyStruct arrays are not used in structs/routines
  // 3D MadEnergyStruct arrays are not used in structs/routines
}

// =============================================================================
// mad_map_struct
void init_mad_map_struct(nb::module_ &m, nb::class_<MadMapStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>>(),
         nb::arg("k") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("t") = nb::none()
  )
      .def_prop_rw(
          "k",
          &MadMapStruct::k,
          &MadMapStruct::set_k,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th order map."
      )
      .def_prop_rw(
          "r",
          &MadMapStruct::r,
          &MadMapStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1st order map."
      )
      .def_prop_rw(
          "t",
          &MadMapStruct::t,
          &MadMapStruct::set_t,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "2nd order map."
      )

      .def("__repr__", [](const MadMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadMapStruct &self) {
            return MadMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MadMapStruct &self, nb::dict &memo) { return MadMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const MadMapStruct &self, const MadMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MadMapStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MadMapStruct arrays are not used in structs/routines
  // 2D MadMapStruct arrays are not used in structs/routines
  // 3D MadMapStruct arrays are not used in structs/routines
}

// =============================================================================
// mesh3d_struct
void init_mesh3d_struct(nb::module_ &m, nb::class_<Mesh3dStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<int>>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>>(),
         nb::arg("nlo") = nb::none(),
         nb::arg("nhi") = nb::none(),
         nb::arg("npad") = nb::none(),
         nb::arg("min") = nb::none(),
         nb::arg("max") = nb::none(),
         nb::arg("delta") = nb::none(),
         nb::arg("gamma") = nb::none(),
         nb::arg("charge") = nb::none(),
         nb::arg("rho") = nb::none(),
         nb::arg("phi") = nb::none()
  )
      .def_prop_rw(
          "nlo",
          &Mesh3dStruct::nlo,
          &Mesh3dStruct::set_nlo,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E)"
      )
      .def_prop_rw(
          "nhi",
          &Mesh3dStruct::nhi,
          &Mesh3dStruct::set_nhi,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E)"
      )
      .def_prop_rw(
          "npad",
          &Mesh3dStruct::npad,
          &Mesh3dStruct::set_npad,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Array padding for cyclic convolution"
      )
      .def_prop_rw(
          "min",
          &Mesh3dStruct::min,
          &Mesh3dStruct::set_min,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Minimim in each dimension"
      )
      .def_prop_rw(
          "max",
          &Mesh3dStruct::max,
          &Mesh3dStruct::set_max,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Maximum in each dimension"
      )
      .def_prop_rw(
          "delta",
          &Mesh3dStruct::delta,
          &Mesh3dStruct::set_delta,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Grid spacing"
      )
      .def_prop_rw("gamma", &Mesh3dStruct::gamma, &Mesh3dStruct::set_gamma, "Relativistic gamma")
      .def_prop_rw(
          "charge",
          &Mesh3dStruct::charge,
          &Mesh3dStruct::set_charge,
          "Total charge on mesh"
      )
      .def_prop_rw(
          "rho",
          &Mesh3dStruct::rho,
          &Mesh3dStruct::set_rho,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Charge density grid"
      )
      .def_prop_rw(
          "phi",
          &Mesh3dStruct::phi,
          &Mesh3dStruct::set_phi,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "electric potential grid"
      )
      // 4D_ALLOC_real efield proxy support missing
      // 4D_ALLOC_real bfield proxy support missing

      .def("__repr__", [](const Mesh3dStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mesh3dStruct &self) {
            return Mesh3dStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Mesh3dStruct &self, nb::dict &memo) { return Mesh3dStruct(self); }
      )
      .def(
          "__eq__",
          [](const Mesh3dStruct &self, const Mesh3dStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Mesh3dStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D Mesh3dStruct arrays are not used in structs/routines
  // 2D Mesh3dStruct arrays are not used in structs/routines
  // 3D Mesh3dStruct arrays are not used in structs/routines
}

// =============================================================================
// molecular_component_struct
void init_molecular_component_struct(nb::module_ &m, nb::class_<MolecularComponentStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("atom") = nb::none(),
         nb::arg("number") = nb::none()
  )
      .def_prop_rw("atom", &MolecularComponentStruct::atom, &MolecularComponentStruct::set_atom)
      .def_prop_rw(
          "number",
          &MolecularComponentStruct::number,
          &MolecularComponentStruct::set_number
      )
      .def_static(
          "new_array1d",
          [](int sz) { return MolecularComponentStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MolecularComponentStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MolecularComponentStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MolecularComponentStruct &self) {
            return MolecularComponentStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MolecularComponentStruct &self, nb::dict &memo) {
            return MolecularComponentStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MolecularComponentStruct &self, const MolecularComponentStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MolecularComponentStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MolecularComponentStructArray1D, MolecularComponentStructAlloc1D>(
      m,
      "MolecularComponentStructArray1D",
      "MolecularComponentStructAlloc1D"
  );
  // 2D MolecularComponentStruct arrays are not used in structs/routines
  // 3D MolecularComponentStruct arrays are not used in structs/routines
}

// =============================================================================
// momentum_aperture_struct
void init_momentum_aperture_struct(nb::module_ &m, nb::class_<MomentumApertureStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("s") = nb::none(),
         nb::arg("pos") = nb::none(),
         nb::arg("neg") = nb::none()
  )
      .def_prop_rw("s", &MomentumApertureStruct::s, &MomentumApertureStruct::set_s)
      .def_prop_rw("pos", &MomentumApertureStruct::pos, &MomentumApertureStruct::set_pos)
      .def_prop_rw("neg", &MomentumApertureStruct::neg, &MomentumApertureStruct::set_neg)
      .def_static(
          "new_array1d",
          [](int sz) { return MomentumApertureStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MomentumApertureStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MomentumApertureStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MomentumApertureStruct &self) {
            return MomentumApertureStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MomentumApertureStruct &self, nb::dict &memo) {
            return MomentumApertureStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MomentumApertureStruct &self, const MomentumApertureStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MomentumApertureStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MomentumApertureStructArray1D, MomentumApertureStructAlloc1D>(
      m,
      "MomentumApertureStructArray1D",
      "MomentumApertureStructAlloc1D"
  );
  // 2D MomentumApertureStruct arrays are not used in structs/routines
  // 3D MomentumApertureStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_region_branch_struct
void init_multipass_region_branch_struct(
    nb::module_ &m,
    nb::class_<MultipassRegionBranchStruct> &cls
) {
  cls.def(nb::init<>())
      .def_prop_ro("ele", &MultipassRegionBranchStruct::ele, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return MultipassRegionBranchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MultipassRegionBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MultipassRegionBranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassRegionBranchStruct &self) {
            return MultipassRegionBranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassRegionBranchStruct &self, nb::dict &memo) {
            return MultipassRegionBranchStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassRegionBranchStruct &self, const MultipassRegionBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassRegionBranchStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MultipassRegionBranchStructArray1D, MultipassRegionBranchStructAlloc1D>(
      m,
      "MultipassRegionBranchStructArray1D",
      "MultipassRegionBranchStructAlloc1D"
  );
  // 2D MultipassRegionBranchStruct arrays are not used in structs/routines
  // 3D MultipassRegionBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_region_ele_struct
void init_multipass_region_ele_struct(nb::module_ &m, nb::class_<MultipassRegionEleStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<bool>, std::optional<bool>>(),
         nb::arg("ix_region") = nb::none(),
         nb::arg("region_start_pt") = nb::none(),
         nb::arg("region_stop_pt") = nb::none()
  )
      .def_prop_rw(
          "ix_region",
          &MultipassRegionEleStruct::ix_region,
          &MultipassRegionEleStruct::set_ix_region
      )
      .def_prop_rw(
          "region_start_pt",
          &MultipassRegionEleStruct::region_start_pt,
          &MultipassRegionEleStruct::set_region_start_pt
      )
      .def_prop_rw(
          "region_stop_pt",
          &MultipassRegionEleStruct::region_stop_pt,
          &MultipassRegionEleStruct::set_region_stop_pt
      )
      .def_static(
          "new_array1d",
          [](int sz) { return MultipassRegionEleStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = MultipassRegionEleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const MultipassRegionEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassRegionEleStruct &self) {
            return MultipassRegionEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassRegionEleStruct &self, nb::dict &memo) {
            return MultipassRegionEleStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassRegionEleStruct &self, const MultipassRegionEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassRegionEleStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<MultipassRegionEleStructArray1D, MultipassRegionEleStructAlloc1D>(
      m,
      "MultipassRegionEleStructArray1D",
      "MultipassRegionEleStructAlloc1D"
  );
  // 2D MultipassRegionEleStruct arrays are not used in structs/routines
  // 3D MultipassRegionEleStruct arrays are not used in structs/routines
}

// =============================================================================
// multipass_region_lat_struct
void init_multipass_region_lat_struct(nb::module_ &m, nb::class_<MultipassRegionLatStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("branch", &MultipassRegionLatStruct::branch, nb::keep_alive<0, 1>())

      .def("__repr__", [](const MultipassRegionLatStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MultipassRegionLatStruct &self) {
            return MultipassRegionLatStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const MultipassRegionLatStruct &self, nb::dict &memo) {
            return MultipassRegionLatStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const MultipassRegionLatStruct &self, const MultipassRegionLatStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const MultipassRegionLatStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D MultipassRegionLatStruct arrays are not used in structs/routines
  // 2D MultipassRegionLatStruct arrays are not used in structs/routines
  // 3D MultipassRegionLatStruct arrays are not used in structs/routines
}