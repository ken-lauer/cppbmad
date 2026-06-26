#include "pybmad/generated/structs_r.hpp"

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
// rad_int1_struct
void init_rad_int1_struct(nb::module_ &m, nb::class_<RadInt1Struct> &cls) {
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
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("i0") = nb::none(),
         nb::arg("i1") = nb::none(),
         nb::arg("i2") = nb::none(),
         nb::arg("i3") = nb::none(),
         nb::arg("i4a") = nb::none(),
         nb::arg("i4b") = nb::none(),
         nb::arg("i4z") = nb::none(),
         nb::arg("i5a") = nb::none(),
         nb::arg("i5b") = nb::none(),
         nb::arg("i6b") = nb::none(),
         nb::arg("lin_i2_E4") = nb::none(),
         nb::arg("lin_i3_E7") = nb::none(),
         nb::arg("lin_i5a_E6") = nb::none(),
         nb::arg("lin_i5b_E6") = nb::none(),
         nb::arg("lin_norm_emit_a") = nb::none(),
         nb::arg("lin_norm_emit_b") = nb::none(),
         nb::arg("lin_sig_E") = nb::none(),
         nb::arg("n_steps") = nb::none()
  )
      .def_prop_rw("i0", &RadInt1Struct::i0, &RadInt1Struct::set_i0)
      .def_prop_rw("i1", &RadInt1Struct::i1, &RadInt1Struct::set_i1)
      .def_prop_rw("i2", &RadInt1Struct::i2, &RadInt1Struct::set_i2)
      .def_prop_rw("i3", &RadInt1Struct::i3, &RadInt1Struct::set_i3)
      .def_prop_rw("i4a", &RadInt1Struct::i4a, &RadInt1Struct::set_i4a)
      .def_prop_rw("i4b", &RadInt1Struct::i4b, &RadInt1Struct::set_i4b)
      .def_prop_rw("i4z", &RadInt1Struct::i4z, &RadInt1Struct::set_i4z)
      .def_prop_rw("i5a", &RadInt1Struct::i5a, &RadInt1Struct::set_i5a)
      .def_prop_rw("i5b", &RadInt1Struct::i5b, &RadInt1Struct::set_i5b)
      .def_prop_rw("i6b", &RadInt1Struct::i6b, &RadInt1Struct::set_i6b)
      .def_prop_rw("lin_i2_E4", &RadInt1Struct::lin_i2_E4, &RadInt1Struct::set_lin_i2_E4)
      .def_prop_rw("lin_i3_E7", &RadInt1Struct::lin_i3_E7, &RadInt1Struct::set_lin_i3_E7)
      .def_prop_rw("lin_i5a_E6", &RadInt1Struct::lin_i5a_E6, &RadInt1Struct::set_lin_i5a_E6)
      .def_prop_rw("lin_i5b_E6", &RadInt1Struct::lin_i5b_E6, &RadInt1Struct::set_lin_i5b_E6)
      .def_prop_rw(
          "lin_norm_emit_a",
          &RadInt1Struct::lin_norm_emit_a,
          &RadInt1Struct::set_lin_norm_emit_a,
          "Running sum"
      )
      .def_prop_rw(
          "lin_norm_emit_b",
          &RadInt1Struct::lin_norm_emit_b,
          &RadInt1Struct::set_lin_norm_emit_b,
          "Running sum"
      )
      .def_prop_rw(
          "lin_sig_E",
          &RadInt1Struct::lin_sig_E,
          &RadInt1Struct::set_lin_sig_E,
          "Running sum"
      )
      .def_prop_rw(
          "n_steps",
          &RadInt1Struct::n_steps,
          &RadInt1Struct::set_n_steps,
          "number of qromb steps needed"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RadInt1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RadInt1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const RadInt1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadInt1Struct &self) {
            return RadInt1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadInt1Struct &self, nb::dict &memo) { return RadInt1Struct(self); }
      )
      .def(
          "__eq__",
          [](const RadInt1Struct &self, const RadInt1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadInt1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<RadInt1StructArray1D, RadInt1StructAlloc1D>(
      m,
      "RadInt1StructArray1D",
      "RadInt1StructAlloc1D"
  );
  // 2D RadInt1Struct arrays are not used in structs/routines
  // 3D RadInt1Struct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_all_ele_struct
void init_rad_int_all_ele_struct(nb::module_ &m, nb::class_<RadIntAllEleStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro(
          "branch",
          &RadIntAllEleStruct::branch,
          nb::keep_alive<0, 1>(),
          "Array is indexed from 0"
      )

      .def("__repr__", [](const RadIntAllEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntAllEleStruct &self) {
            return RadIntAllEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntAllEleStruct &self, nb::dict &memo) { return RadIntAllEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const RadIntAllEleStruct &self, const RadIntAllEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadIntAllEleStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RadIntAllEleStruct arrays are not used in structs/routines
  // 2D RadIntAllEleStruct arrays are not used in structs/routines
  // 3D RadIntAllEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_branch_struct
void init_rad_int_branch_struct(nb::module_ &m, nb::class_<RadIntBranchStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro(
          "ele",
          &RadIntBranchStruct::ele,
          nb::keep_alive<0, 1>(),
          "Array is indexed from 0"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RadIntBranchStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RadIntBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const RadIntBranchStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntBranchStruct &self) {
            return RadIntBranchStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntBranchStruct &self, nb::dict &memo) { return RadIntBranchStruct(self); }
      )
      .def(
          "__eq__",
          [](const RadIntBranchStruct &self, const RadIntBranchStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadIntBranchStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<RadIntBranchStructArray1D, RadIntBranchStructAlloc1D>(
      m,
      "RadIntBranchStructArray1D",
      "RadIntBranchStructAlloc1D"
  );
  // 2D RadIntBranchStruct arrays are not used in structs/routines
  // 3D RadIntBranchStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_map_ele_struct
void init_rad_map_ele_struct(nb::module_ &m, nb::class_<RadMapEleStruct> &cls) {
  cls.def(
         "__init__",
         [](RadMapEleStruct *self,
            const RadMapStruct *rm0,
            const RadMapStruct *rm1,
            std::optional<bool> stale) {
           new (self) RadMapEleStruct(ptr_to_opt_ref(rm0), ptr_to_opt_ref(rm1), stale);
         },
         nb::arg("rm0") = nb::none(),
         nb::arg("rm1") = nb::none(),
         nb::arg("stale") = nb::none()
  )
      .def_prop_rw(
          "rm0",
          &RadMapEleStruct::rm0,
          &RadMapEleStruct::set_rm0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Upstream half and downstream half matrices for an element."
      )
      .def_prop_rw(
          "rm1",
          &RadMapEleStruct::rm1,
          &RadMapEleStruct::set_rm1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Upstream half and downstream half matrices for an element."
      )
      .def_prop_rw("stale", &RadMapEleStruct::stale, &RadMapEleStruct::set_stale)

      .def("__repr__", [](const RadMapEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapEleStruct &self) {
            return RadMapEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadMapEleStruct &self, nb::dict &memo) { return RadMapEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const RadMapEleStruct &self, const RadMapEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadMapEleStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RadMapEleStruct arrays are not used in structs/routines
  // 2D RadMapEleStruct arrays are not used in structs/routines
  // 3D RadMapEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_map_struct
void init_rad_map_struct(nb::module_ &m, nb::class_<RadMapStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>>(),
         nb::arg("ref_orb") = nb::none(),
         nb::arg("damp_dmat") = nb::none(),
         nb::arg("xfer_damp_vec") = nb::none(),
         nb::arg("xfer_damp_mat") = nb::none(),
         nb::arg("stoc_mat") = nb::none()
  )
      .def_prop_rw(
          "ref_orb",
          &RadMapStruct::ref_orb,
          &RadMapStruct::set_ref_orb,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Reference point around which damp_mat is calculated."
      )
      .def_prop_rw(
          "damp_dmat",
          &RadMapStruct::damp_dmat,
          &RadMapStruct::set_damp_dmat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "damp_correction = xfer_mat_with_damping - xfer_mat_without_damping."
      )
      .def_prop_rw(
          "xfer_damp_vec",
          &RadMapStruct::xfer_damp_vec,
          &RadMapStruct::set_xfer_damp_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Transfer map with damping 0th order vector."
      )
      .def_prop_rw(
          "xfer_damp_mat",
          &RadMapStruct::xfer_damp_mat,
          &RadMapStruct::set_xfer_damp_mat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1st order matrix: xfer_no_damp_mat + xfer_damp_correction."
      )
      .def_prop_rw(
          "stoc_mat",
          &RadMapStruct::stoc_mat,
          &RadMapStruct::set_stoc_mat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Stochastic variance or 'kick' (Cholesky decomposed) matrix."
      )

      .def("__repr__", [](const RadMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapStruct &self) {
            return RadMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadMapStruct &self, nb::dict &memo) { return RadMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const RadMapStruct &self, const RadMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadMapStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RadMapStruct arrays are not used in structs/routines
  // 2D RadMapStruct arrays are not used in structs/routines
  // 3D RadMapStruct arrays are not used in structs/routines
}

// =============================================================================
// ramper_lord_struct
void init_ramper_lord_struct(nb::module_ &m, nb::class_<RamperLordStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<int>, std::optional<double>>(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_con") = nb::none(),
         nb::arg("attrib_ptr") = nb::none()
  )
      .def_prop_rw("ix_ele", &RamperLordStruct::ix_ele, &RamperLordStruct::set_ix_ele, "Lord index")
      .def_prop_rw(
          "ix_con",
          &RamperLordStruct::ix_con,
          &RamperLordStruct::set_ix_con,
          "Index in lord%control%ramp(:) array"
      )
      .def_prop_rw(
          "attrib_ptr",
          &RamperLordStruct::attrib_ptr,
          &RamperLordStruct::set_attrib_ptr,
          "Pointer to attribute in this element."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RamperLordStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RamperLordStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const RamperLordStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RamperLordStruct &self) {
            return RamperLordStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RamperLordStruct &self, nb::dict &memo) { return RamperLordStruct(self); }
      )
      .def(
          "__eq__",
          [](const RamperLordStruct &self, const RamperLordStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RamperLordStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<RamperLordStructArray1D, RamperLordStructAlloc1D>(
      m,
      "RamperLordStructArray1D",
      "RamperLordStructAlloc1D"
  );
  // 2D RamperLordStruct arrays are not used in structs/routines
  // 3D RamperLordStruct arrays are not used in structs/routines
}

// =============================================================================
// resonance_h_struct
void init_resonance_h_struct(nb::module_ &m, nb::class_<ResonanceHStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<std::complex<double>>>(),
         nb::arg("id") = nb::none(),
         nb::arg("c_val") = nb::none()
  )
      .def_prop_rw(
          "id",
          &ResonanceHStruct::id,
          &ResonanceHStruct::set_id,
          "6 digit ID. EG: '003100'"
      )
      .def_prop_rw(
          "c_val",
          &ResonanceHStruct::c_val,
          &ResonanceHStruct::set_c_val,
          "Resonance value"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ResonanceHStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ResonanceHStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ResonanceHStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ResonanceHStruct &self) {
            return ResonanceHStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ResonanceHStruct &self, nb::dict &memo) { return ResonanceHStruct(self); }
      )
      .def(
          "__eq__",
          [](const ResonanceHStruct &self, const ResonanceHStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ResonanceHStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ResonanceHStructArray1D, ResonanceHStructAlloc1D>(
      m,
      "ResonanceHStructArray1D",
      "ResonanceHStructAlloc1D"
  );
  // 2D ResonanceHStruct arrays are not used in structs/routines
  // 3D ResonanceHStruct arrays are not used in structs/routines
}

// =============================================================================
// rf_ele_struct
void init_rf_ele_struct(nb::module_ &m, nb::class_<RfEleStruct> &cls) {
  cls.def(nb::init<std::optional<double>>(), nb::arg("ds_step") = nb::none())
      .def_prop_ro(
          "steps",
          &RfEleStruct::steps,
          nb::keep_alive<0, 1>(),
          "Energy stair step array indexed from zero."
      )
      .def_prop_rw(
          "ds_step",
          &RfEleStruct::ds_step,
          &RfEleStruct::set_ds_step,
          "length of a stair step."
      )

      .def("__repr__", [](const RfEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfEleStruct &self) {
            return RfEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RfEleStruct &self, nb::dict &memo) { return RfEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const RfEleStruct &self, const RfEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RfEleStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RfEleStruct arrays are not used in structs/routines
  // 2D RfEleStruct arrays are not used in structs/routines
  // 3D RfEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rf_stair_step_struct
void init_rf_stair_step_struct(nb::module_ &m, nb::class_<RfStairStepStruct> &cls) {
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
         nb::arg("E_tot0") = nb::none(),
         nb::arg("E_tot1") = nb::none(),
         nb::arg("p0c") = nb::none(),
         nb::arg("p1c") = nb::none(),
         nb::arg("scale") = nb::none(),
         nb::arg("time") = nb::none(),
         nb::arg("s0") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("ix_step") = nb::none()
  )
      .def_prop_rw(
          "E_tot0",
          &RfStairStepStruct::E_tot0,
          &RfStairStepStruct::set_E_tot0,
          "Reference energy in the drift region (before the kick point)."
      )
      .def_prop_rw(
          "E_tot1",
          &RfStairStepStruct::E_tot1,
          &RfStairStepStruct::set_E_tot1,
          "Reference energy after the kick point."
      )
      .def_prop_rw(
          "p0c",
          &RfStairStepStruct::p0c,
          &RfStairStepStruct::set_p0c,
          "Reference momentum in the drift region (before the kick point)."
      )
      .def_prop_rw(
          "p1c",
          &RfStairStepStruct::p1c,
          &RfStairStepStruct::set_p1c,
          "Reference momentum after the kick point."
      )
      .def_prop_rw(
          "scale",
          &RfStairStepStruct::scale,
          &RfStairStepStruct::set_scale,
          "Scale for multipole kick at the kick point. Sum over all steps will be 1."
      )
      .def_prop_rw(
          "time",
          &RfStairStepStruct::time,
          &RfStairStepStruct::set_time,
          "Reference particle time at the kick point with respect to beginning of element."
      )
      .def_prop_rw(
          "s0",
          &RfStairStepStruct::s0,
          &RfStairStepStruct::set_s0,
          "S-position at beginning of drift region relative to the beginning of the element."
      )
      .def_prop_rw(
          "s",
          &RfStairStepStruct::s,
          &RfStairStepStruct::set_s,
          "S-position at the kick point relative to the beginning of the element."
      )
      .def_prop_rw(
          "ix_step",
          &RfStairStepStruct::ix_step,
          &RfStairStepStruct::set_ix_step,
          "Step index in ele%rf%steps(:) array"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RfStairStepStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RfStairStepStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const RfStairStepStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfStairStepStruct &self) {
            return RfStairStepStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RfStairStepStruct &self, nb::dict &memo) { return RfStairStepStruct(self); }
      )
      .def(
          "__eq__",
          [](const RfStairStepStruct &self, const RfStairStepStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RfStairStepStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<RfStairStepStructArray1D, RfStairStepStructAlloc1D>(
      m,
      "RfStairStepStructArray1D",
      "RfStairStepStructAlloc1D"
  );
  // 2D RfStairStepStruct arrays are not used in structs/routines
  // 3D RfStairStepStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_cache1_struct
void init_rad_int_cache1_struct(nb::module_ &m, nb::class_<RadIntCache1Struct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<int>>(),
         nb::arg("n_pt") = nb::none(),
         nb::arg("cache_type") = nb::none()
  )
      .def_prop_ro("pt", &RadIntCache1Struct::pt, nb::keep_alive<0, 1>(), "pt(0:n_pt)")
      .def_prop_rw(
          "n_pt",
          &RadIntCache1Struct::n_pt,
          &RadIntCache1Struct::set_n_pt,
          "Upper bound of pt(0:n_pt)"
      )
      .def_prop_rw(
          "cache_type",
          &RadIntCache1Struct::cache_type,
          &RadIntCache1Struct::set_cache_type
      )

      .def("__repr__", [](const RadIntCache1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntCache1Struct &self) {
            return RadIntCache1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntCache1Struct &self, nb::dict &memo) { return RadIntCache1Struct(self); }
      )
      .def(
          "__eq__",
          [](const RadIntCache1Struct &self, const RadIntCache1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadIntCache1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RadIntCache1Struct arrays are not used in structs/routines
  // 2D RadIntCache1Struct arrays are not used in structs/routines
  // 3D RadIntCache1Struct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_info_struct
void init_rad_int_info_struct(nb::module_ &m, nb::class_<RadIntInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](RadIntInfoStruct *self,
            const BranchStruct *branch,
            const EleStruct *ele,
            const TwissStruct *a,
            const TwissStruct *b,
            const RadIntCache1Struct *cache_ele,
            std::optional<std::vector<double>> eta_a,
            std::optional<std::vector<double>> eta_b,
            std::optional<double> g,
            std::optional<double> g2,
            std::optional<double> g_x,
            std::optional<double> g_y,
            std::optional<double> dg2_x,
            std::optional<double> dg2_y) {
           new (self) RadIntInfoStruct(
               ptr_to_opt_ref(branch),
               ptr_to_opt_ref(ele),
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(cache_ele),
               eta_a,
               eta_b,
               g,
               g2,
               g_x,
               g_y,
               dg2_x,
               dg2_y
           );
         },
         nb::arg("branch") = nb::none(),
         nb::arg("ele") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("cache_ele") = nb::none(),
         nb::arg("eta_a") = nb::none(),
         nb::arg("eta_b") = nb::none(),
         nb::arg("g") = nb::none(),
         nb::arg("g2") = nb::none(),
         nb::arg("g_x") = nb::none(),
         nb::arg("g_y") = nb::none(),
         nb::arg("dg2_x") = nb::none(),
         nb::arg("dg2_y") = nb::none()
  )
      .def_prop_rw(
          "branch",
          &RadIntInfoStruct::branch,
          &RadIntInfoStruct::set_branch,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ele",
          &RadIntInfoStruct::ele,
          &RadIntInfoStruct::set_ele,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("orbit", &RadIntInfoStruct::orbit, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "a",
          &RadIntInfoStruct::a,
          &RadIntInfoStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "b",
          &RadIntInfoStruct::b,
          &RadIntInfoStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "cache_ele",
          &RadIntInfoStruct::cache_ele,
          &RadIntInfoStruct::set_cache_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "pointer to cache in use"
      )
      .def_prop_rw(
          "eta_a",
          &RadIntInfoStruct::eta_a,
          &RadIntInfoStruct::set_eta_a,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "eta_b",
          &RadIntInfoStruct::eta_b,
          &RadIntInfoStruct::set_eta_b,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "g",
          &RadIntInfoStruct::g,
          &RadIntInfoStruct::set_g,
          "bending strength (1/bending_radius)"
      )
      .def_prop_rw(
          "g2",
          &RadIntInfoStruct::g2,
          &RadIntInfoStruct::set_g2,
          "bending strength (1/bending_radius)"
      )
      .def_prop_rw(
          "g_x",
          &RadIntInfoStruct::g_x,
          &RadIntInfoStruct::set_g_x,
          "components in x-y plane"
      )
      .def_prop_rw(
          "g_y",
          &RadIntInfoStruct::g_y,
          &RadIntInfoStruct::set_g_y,
          "components in x-y plane"
      )
      .def_prop_rw("dg2_x", &RadIntInfoStruct::dg2_x, &RadIntInfoStruct::set_dg2_x)
      .def_prop_rw("dg2_y", &RadIntInfoStruct::dg2_y, &RadIntInfoStruct::set_dg2_y)

      .def("__repr__", [](const RadIntInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntInfoStruct &self) {
            return RadIntInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntInfoStruct &self, nb::dict &memo) { return RadIntInfoStruct(self); }
      )
      .def(
          "__eq__",
          [](const RadIntInfoStruct &self, const RadIntInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadIntInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RadIntInfoStruct arrays are not used in structs/routines
  // 2D RadIntInfoStruct arrays are not used in structs/routines
  // 3D RadIntInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_track_point_struct
void init_rad_int_track_point_struct(nb::module_ &m, nb::class_<RadIntTrackPointStruct> &cls) {
  cls.def(
         "__init__",
         [](RadIntTrackPointStruct *self,
            std::optional<double> s_body,
            std::optional<std::vector<std::vector<double>>> mat6,
            std::optional<std::vector<double>> vec0,
            const CoordStruct *ref_orb_in,
            const CoordStruct *ref_orb_out,
            std::optional<double> g_x0,
            std::optional<double> g_y0,
            std::optional<double> dgx_dx,
            std::optional<double> dgx_dy,
            std::optional<double> dgy_dx,
            std::optional<double> dgy_dy) {
           new (self) RadIntTrackPointStruct(
               s_body,
               mat6,
               vec0,
               ptr_to_opt_ref(ref_orb_in),
               ptr_to_opt_ref(ref_orb_out),
               g_x0,
               g_y0,
               dgx_dx,
               dgx_dy,
               dgy_dx,
               dgy_dy
           );
         },
         nb::arg("s_body") = nb::none(),
         nb::arg("mat6") = nb::none(),
         nb::arg("vec0") = nb::none(),
         nb::arg("ref_orb_in") = nb::none(),
         nb::arg("ref_orb_out") = nb::none(),
         nb::arg("g_x0") = nb::none(),
         nb::arg("g_y0") = nb::none(),
         nb::arg("dgx_dx") = nb::none(),
         nb::arg("dgx_dy") = nb::none(),
         nb::arg("dgy_dx") = nb::none(),
         nb::arg("dgy_dy") = nb::none()
  )
      .def_prop_rw("s_body", &RadIntTrackPointStruct::s_body, &RadIntTrackPointStruct::set_s_body)
      .def_prop_rw(
          "mat6",
          &RadIntTrackPointStruct::mat6,
          &RadIntTrackPointStruct::set_mat6,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "vec0",
          &RadIntTrackPointStruct::vec0,
          &RadIntTrackPointStruct::set_vec0,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ref_orb_in",
          &RadIntTrackPointStruct::ref_orb_in,
          &RadIntTrackPointStruct::set_ref_orb_in,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ref_orb_out",
          &RadIntTrackPointStruct::ref_orb_out,
          &RadIntTrackPointStruct::set_ref_orb_out,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "g_x0",
          &RadIntTrackPointStruct::g_x0,
          &RadIntTrackPointStruct::set_g_x0,
          "Additional g factors for bends."
      )
      .def_prop_rw(
          "g_y0",
          &RadIntTrackPointStruct::g_y0,
          &RadIntTrackPointStruct::set_g_y0,
          "Additional g factors for bends."
      )
      .def_prop_rw(
          "dgx_dx",
          &RadIntTrackPointStruct::dgx_dx,
          &RadIntTrackPointStruct::set_dgx_dx,
          "bending strength gradient"
      )
      .def_prop_rw(
          "dgx_dy",
          &RadIntTrackPointStruct::dgx_dy,
          &RadIntTrackPointStruct::set_dgx_dy,
          "bending strength gradient"
      )
      .def_prop_rw(
          "dgy_dx",
          &RadIntTrackPointStruct::dgy_dx,
          &RadIntTrackPointStruct::set_dgy_dx,
          "bending strength gradient"
      )
      .def_prop_rw(
          "dgy_dy",
          &RadIntTrackPointStruct::dgy_dy,
          &RadIntTrackPointStruct::set_dgy_dy,
          "bending strength gradient"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RadIntTrackPointStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RadIntTrackPointStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const RadIntTrackPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntTrackPointStruct &self) {
            return RadIntTrackPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntTrackPointStruct &self, nb::dict &memo) {
            return RadIntTrackPointStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const RadIntTrackPointStruct &self, const RadIntTrackPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RadIntTrackPointStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<RadIntTrackPointStructArray1D, RadIntTrackPointStructAlloc1D>(
      m,
      "RadIntTrackPointStructArray1D",
      "RadIntTrackPointStructAlloc1D"
  );
  // 2D RadIntTrackPointStruct arrays are not used in structs/routines
  // 3D RadIntTrackPointStruct arrays are not used in structs/routines
}

// =============================================================================
// random_state_struct
void init_random_state_struct(nb::module_ &m, nb::class_<RandomStateStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int64_t>,
             std::optional<int64_t>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<int>,
             std::optional<double>,
             std::optional<int64_t>,
             std::optional<std::vector<int64_t>>,
             std::optional<std::vector<double>>>(),
         nb::arg("ix") = nb::none(),
         nb::arg("iy") = nb::none(),
         nb::arg("number_stored") = nb::none(),
         nb::arg("h_saved") = nb::none(),
         nb::arg("engine") = nb::none(),
         nb::arg("seed") = nb::none(),
         nb::arg("am") = nb::none(),
         nb::arg("gauss_converter") = nb::none(),
         nb::arg("gauss_sigma_cut") = nb::none(),
         nb::arg("in_sobseq") = nb::none(),
         nb::arg("ix_sobseq") = nb::none(),
         nb::arg("x_sobseq") = nb::none()
  )
      .def_prop_rw("ix", &RandomStateStruct::ix, &RandomStateStruct::set_ix)
      .def_prop_rw("iy", &RandomStateStruct::iy, &RandomStateStruct::set_iy)
      .def_prop_rw(
          "number_stored",
          &RandomStateStruct::number_stored,
          &RandomStateStruct::set_number_stored
      )
      .def_prop_rw("h_saved", &RandomStateStruct::h_saved, &RandomStateStruct::set_h_saved)
      .def_prop_rw("engine", &RandomStateStruct::engine, &RandomStateStruct::set_engine, "Params")
      .def_prop_rw("seed", &RandomStateStruct::seed, &RandomStateStruct::set_seed)
      .def_prop_rw("am", &RandomStateStruct::am, &RandomStateStruct::set_am)
      .def_prop_rw(
          "gauss_converter",
          &RandomStateStruct::gauss_converter,
          &RandomStateStruct::set_gauss_converter
      )
      .def_prop_rw(
          "gauss_sigma_cut",
          &RandomStateStruct::gauss_sigma_cut,
          &RandomStateStruct::set_gauss_sigma_cut,
          "Only used if positive."
      )
      .def_prop_rw("in_sobseq", &RandomStateStruct::in_sobseq, &RandomStateStruct::set_in_sobseq)
      .def_prop_rw(
          "ix_sobseq",
          &RandomStateStruct::ix_sobseq,
          &RandomStateStruct::set_ix_sobseq,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "x_sobseq",
          &RandomStateStruct::x_sobseq,
          &RandomStateStruct::set_x_sobseq,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const RandomStateStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RandomStateStruct &self) {
            return RandomStateStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RandomStateStruct &self, nb::dict &memo) { return RandomStateStruct(self); }
      )
      .def(
          "__eq__",
          [](const RandomStateStruct &self, const RandomStateStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const RandomStateStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D RandomStateStruct arrays are not used in structs/routines
  // 2D RandomStateStruct arrays are not used in structs/routines
  // 3D RandomStateStruct arrays are not used in structs/routines
}