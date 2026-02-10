#include "pybmad/generated/structs_r.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// rad_int1_struct
void init_rad_int1_struct(py::module &m, py::class_<RadInt1Struct> &cls) {
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
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("i0") = py::none(),
         py::arg("i1") = py::none(),
         py::arg("i2") = py::none(),
         py::arg("i3") = py::none(),
         py::arg("i4a") = py::none(),
         py::arg("i4b") = py::none(),
         py::arg("i4z") = py::none(),
         py::arg("i5a") = py::none(),
         py::arg("i5b") = py::none(),
         py::arg("i6b") = py::none(),
         py::arg("lin_i2_E4") = py::none(),
         py::arg("lin_i3_E7") = py::none(),
         py::arg("lin_i5a_E6") = py::none(),
         py::arg("lin_i5b_E6") = py::none(),
         py::arg("lin_norm_emit_a") = py::none(),
         py::arg("lin_norm_emit_b") = py::none(),
         py::arg("lin_sig_E") = py::none(),
         py::arg("n_steps") = py::none()
  )
      .def_property("i0", &RadInt1Struct::i0, &RadInt1Struct::set_i0)
      .def_property("i1", &RadInt1Struct::i1, &RadInt1Struct::set_i1)
      .def_property("i2", &RadInt1Struct::i2, &RadInt1Struct::set_i2)
      .def_property("i3", &RadInt1Struct::i3, &RadInt1Struct::set_i3)
      .def_property("i4a", &RadInt1Struct::i4a, &RadInt1Struct::set_i4a)
      .def_property("i4b", &RadInt1Struct::i4b, &RadInt1Struct::set_i4b)
      .def_property("i4z", &RadInt1Struct::i4z, &RadInt1Struct::set_i4z)
      .def_property("i5a", &RadInt1Struct::i5a, &RadInt1Struct::set_i5a)
      .def_property("i5b", &RadInt1Struct::i5b, &RadInt1Struct::set_i5b)
      .def_property("i6b", &RadInt1Struct::i6b, &RadInt1Struct::set_i6b)
      .def_property("lin_i2_E4", &RadInt1Struct::lin_i2_E4, &RadInt1Struct::set_lin_i2_E4)
      .def_property("lin_i3_E7", &RadInt1Struct::lin_i3_E7, &RadInt1Struct::set_lin_i3_E7)
      .def_property("lin_i5a_E6", &RadInt1Struct::lin_i5a_E6, &RadInt1Struct::set_lin_i5a_E6)
      .def_property("lin_i5b_E6", &RadInt1Struct::lin_i5b_E6, &RadInt1Struct::set_lin_i5b_E6)
      .def_property(
          "lin_norm_emit_a",
          &RadInt1Struct::lin_norm_emit_a,
          &RadInt1Struct::set_lin_norm_emit_a,
          "Running sum"
      )
      .def_property(
          "lin_norm_emit_b",
          &RadInt1Struct::lin_norm_emit_b,
          &RadInt1Struct::set_lin_norm_emit_b,
          "Running sum"
      )
      .def_property(
          "lin_sig_E",
          &RadInt1Struct::lin_sig_E,
          &RadInt1Struct::set_lin_sig_E,
          "Running sum"
      )
      .def_property(
          "n_steps",
          &RadInt1Struct::n_steps,
          &RadInt1Struct::set_n_steps,
          "number of qromb steps needed"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RadInt1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RadInt1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const RadInt1Struct &self, py::dict &memo) { return RadInt1Struct(self); }
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
void init_rad_int_all_ele_struct(py::module &m, py::class_<RadIntAllEleStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly("branch", &RadIntAllEleStruct::branch, "Array is indexed from 0")

      .def("__repr__", [](const RadIntAllEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntAllEleStruct &self) {
            return RadIntAllEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadIntAllEleStruct &self, py::dict &memo) { return RadIntAllEleStruct(self); }
      )

      ;

  // 1D RadIntAllEleStruct arrays are not used in structs/routines
  // 2D RadIntAllEleStruct arrays are not used in structs/routines
  // 3D RadIntAllEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_int_branch_struct
void init_rad_int_branch_struct(py::module &m, py::class_<RadIntBranchStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly("ele", &RadIntBranchStruct::ele, "Array is indexed from 0")
      .def_static(
          "new_array1d",
          [](int sz) { return RadIntBranchStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RadIntBranchStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const RadIntBranchStruct &self, py::dict &memo) { return RadIntBranchStruct(self); }
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
void init_rad_map_ele_struct(py::module &m, py::class_<RadMapEleStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const RadMapStruct>,
             optional_ref<const RadMapStruct>,
             std::optional<bool>>(),
         py::arg("rm0") = py::none(),
         py::arg("rm1") = py::none(),
         py::arg("stale") = py::none()
  )
      .def_property(
          "rm0",
          &RadMapEleStruct::rm0,
          &RadMapEleStruct::set_rm0,
          "Upstream half and downstream half matrices for an element."
      )
      .def_property(
          "rm1",
          &RadMapEleStruct::rm1,
          &RadMapEleStruct::set_rm1,
          "Upstream half and downstream half matrices for an element."
      )
      .def_property("stale", &RadMapEleStruct::stale, &RadMapEleStruct::set_stale)

      .def("__repr__", [](const RadMapEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapEleStruct &self) {
            return RadMapEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RadMapEleStruct &self, py::dict &memo) { return RadMapEleStruct(self); }
      )

      ;

  // 1D RadMapEleStruct arrays are not used in structs/routines
  // 2D RadMapEleStruct arrays are not used in structs/routines
  // 3D RadMapEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rad_map_struct
void init_rad_map_struct(py::module &m, py::class_<RadMapStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("ref_orb") = py::none(),
         py::arg("damp_dmat") = py::none(),
         py::arg("xfer_damp_vec") = py::none(),
         py::arg("xfer_damp_mat") = py::none(),
         py::arg("stoc_mat") = py::none()
  )
      .def_property(
          "ref_orb",
          &RadMapStruct::ref_orb,
          &RadMapStruct::set_ref_orb,
          "Reference point around which damp_mat is calculated."
      )
      .def_property(
          "damp_dmat",
          &RadMapStruct::damp_dmat,
          &RadMapStruct::set_damp_dmat,
          "damp_correction = xfer_mat_with_damping - xfer_mat_without_damping."
      )
      .def_property(
          "xfer_damp_vec",
          &RadMapStruct::xfer_damp_vec,
          &RadMapStruct::set_xfer_damp_vec,
          "Transfer map with damping 0th order vector."
      )
      .def_property(
          "xfer_damp_mat",
          &RadMapStruct::xfer_damp_mat,
          &RadMapStruct::set_xfer_damp_mat,
          "1st order matrix: xfer_no_damp_mat + xfer_damp_correction."
      )
      .def_property(
          "stoc_mat",
          &RadMapStruct::stoc_mat,
          &RadMapStruct::set_stoc_mat,
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
          [](const RadMapStruct &self, py::dict &memo) { return RadMapStruct(self); }
      )

      ;

  // 1D RadMapStruct arrays are not used in structs/routines
  // 2D RadMapStruct arrays are not used in structs/routines
  // 3D RadMapStruct arrays are not used in structs/routines
}

// =============================================================================
// ramper_lord_struct
void init_ramper_lord_struct(py::module &m, py::class_<RamperLordStruct> &cls) {
  cls.def(
         py::init<std::optional<int>, std::optional<int>, std::optional<double>>(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_con") = py::none(),
         py::arg("attrib_ptr") = py::none()
  )
      .def_property(
          "ix_ele",
          &RamperLordStruct::ix_ele,
          &RamperLordStruct::set_ix_ele,
          "Lord index"
      )
      .def_property(
          "ix_con",
          &RamperLordStruct::ix_con,
          &RamperLordStruct::set_ix_con,
          "Index in lord%control%ramp(:) array"
      )
      .def_property(
          "attrib_ptr",
          &RamperLordStruct::attrib_ptr,
          &RamperLordStruct::set_attrib_ptr,
          "Pointer to attribute in this element."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RamperLordStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RamperLordStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const RamperLordStruct &self, py::dict &memo) { return RamperLordStruct(self); }
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
void init_resonance_h_struct(py::module &m, py::class_<ResonanceHStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<std::complex<double>>>(),
         py::arg("id") = py::none(),
         py::arg("c_val") = py::none()
  )
      .def_property(
          "id",
          &ResonanceHStruct::id,
          &ResonanceHStruct::set_id,
          "6 digit ID. EG: '003100'"
      )
      .def_property(
          "c_val",
          &ResonanceHStruct::c_val,
          &ResonanceHStruct::set_c_val,
          "Resonance value"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ResonanceHStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ResonanceHStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ResonanceHStruct &self, py::dict &memo) { return ResonanceHStruct(self); }
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
void init_rf_ele_struct(py::module &m, py::class_<RfEleStruct> &cls) {
  cls.def(py::init<std::optional<double>>(), py::arg("ds_step") = py::none())
      .def_property_readonly(
          "steps",
          &RfEleStruct::steps,
          "Energy stair step array indexed from zero."
      )
      .def_property(
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
          [](const RfEleStruct &self, py::dict &memo) { return RfEleStruct(self); }
      )

      ;

  // 1D RfEleStruct arrays are not used in structs/routines
  // 2D RfEleStruct arrays are not used in structs/routines
  // 3D RfEleStruct arrays are not used in structs/routines
}

// =============================================================================
// rf_stair_step_struct
void init_rf_stair_step_struct(py::module &m, py::class_<RfStairStepStruct> &cls) {
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
         py::arg("E_tot0") = py::none(),
         py::arg("E_tot1") = py::none(),
         py::arg("p0c") = py::none(),
         py::arg("p1c") = py::none(),
         py::arg("scale") = py::none(),
         py::arg("time") = py::none(),
         py::arg("s0") = py::none(),
         py::arg("s") = py::none(),
         py::arg("ix_step") = py::none()
  )
      .def_property(
          "E_tot0",
          &RfStairStepStruct::E_tot0,
          &RfStairStepStruct::set_E_tot0,
          "Reference energy in the drift region (before the kick point)."
      )
      .def_property(
          "E_tot1",
          &RfStairStepStruct::E_tot1,
          &RfStairStepStruct::set_E_tot1,
          "Reference energy after the kick point."
      )
      .def_property(
          "p0c",
          &RfStairStepStruct::p0c,
          &RfStairStepStruct::set_p0c,
          "Reference momentum in the drift region (before the kick point)."
      )
      .def_property(
          "p1c",
          &RfStairStepStruct::p1c,
          &RfStairStepStruct::set_p1c,
          "Reference momentum after the kick point."
      )
      .def_property(
          "scale",
          &RfStairStepStruct::scale,
          &RfStairStepStruct::set_scale,
          "Scale for multipole kick at the kick point. Sum over all steps will be 1."
      )
      .def_property(
          "time",
          &RfStairStepStruct::time,
          &RfStairStepStruct::set_time,
          "Reference particle time at the kick point with respect to beginning of element."
      )
      .def_property(
          "s0",
          &RfStairStepStruct::s0,
          &RfStairStepStruct::set_s0,
          "S-position at beginning of drift region relative to the beginning of the element."
      )
      .def_property(
          "s",
          &RfStairStepStruct::s,
          &RfStairStepStruct::set_s,
          "S-position at the kick point relative to the beginning of the element."
      )
      .def_property(
          "ix_step",
          &RfStairStepStruct::ix_step,
          &RfStairStepStruct::set_ix_step,
          "Step index in ele%rf%steps(:) array"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return RfStairStepStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = RfStairStepStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const RfStairStepStruct &self, py::dict &memo) { return RfStairStepStruct(self); }
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
// random_state_struct
void init_random_state_struct(py::module &m, py::class_<RandomStateStruct> &cls) {
  cls.def(
         py::init<
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
         py::arg("ix") = py::none(),
         py::arg("iy") = py::none(),
         py::arg("number_stored") = py::none(),
         py::arg("h_saved") = py::none(),
         py::arg("engine") = py::none(),
         py::arg("seed") = py::none(),
         py::arg("am") = py::none(),
         py::arg("gauss_converter") = py::none(),
         py::arg("gauss_sigma_cut") = py::none(),
         py::arg("in_sobseq") = py::none(),
         py::arg("ix_sobseq") = py::none(),
         py::arg("x_sobseq") = py::none()
  )
      .def_property("ix", &RandomStateStruct::ix, &RandomStateStruct::set_ix)
      .def_property("iy", &RandomStateStruct::iy, &RandomStateStruct::set_iy)
      .def_property(
          "number_stored",
          &RandomStateStruct::number_stored,
          &RandomStateStruct::set_number_stored
      )
      .def_property("h_saved", &RandomStateStruct::h_saved, &RandomStateStruct::set_h_saved)
      .def_property("engine", &RandomStateStruct::engine, &RandomStateStruct::set_engine, "Params")
      .def_property("seed", &RandomStateStruct::seed, &RandomStateStruct::set_seed)
      .def_property("am", &RandomStateStruct::am, &RandomStateStruct::set_am)
      .def_property(
          "gauss_converter",
          &RandomStateStruct::gauss_converter,
          &RandomStateStruct::set_gauss_converter
      )
      .def_property(
          "gauss_sigma_cut",
          &RandomStateStruct::gauss_sigma_cut,
          &RandomStateStruct::set_gauss_sigma_cut,
          "Only used if positive."
      )
      .def_property("in_sobseq", &RandomStateStruct::in_sobseq, &RandomStateStruct::set_in_sobseq)
      .def_property("ix_sobseq", &RandomStateStruct::ix_sobseq, &RandomStateStruct::set_ix_sobseq)
      .def_property("x_sobseq", &RandomStateStruct::x_sobseq, &RandomStateStruct::set_x_sobseq)

      .def("__repr__", [](const RandomStateStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RandomStateStruct &self) {
            return RandomStateStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const RandomStateStruct &self, py::dict &memo) { return RandomStateStruct(self); }
      )

      ;

  // 1D RandomStateStruct arrays are not used in structs/routines
  // 2D RandomStateStruct arrays are not used in structs/routines
  // 3D RandomStateStruct arrays are not used in structs/routines
}