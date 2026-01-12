#include "pybmad/generated/structs_r.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// rad_int1_struct
void init_rad_int1_struct(py::module& m, py::class_<RadInt1Proxy>& cls) {
  cls.def(py::init<>())
      // RadInt1Proxy.i0 (0D_NOT_real -
      .def_property("i0", &RadInt1Proxy::i0, &RadInt1Proxy::set_i0)
      // RadInt1Proxy.i1 (0D_NOT_real -
      .def_property("i1", &RadInt1Proxy::i1, &RadInt1Proxy::set_i1)
      // RadInt1Proxy.i2 (0D_NOT_real -
      .def_property("i2", &RadInt1Proxy::i2, &RadInt1Proxy::set_i2)
      // RadInt1Proxy.i3 (0D_NOT_real -
      .def_property("i3", &RadInt1Proxy::i3, &RadInt1Proxy::set_i3)
      // RadInt1Proxy.i4a (0D_NOT_real -
      .def_property("i4a", &RadInt1Proxy::i4a, &RadInt1Proxy::set_i4a)
      // RadInt1Proxy.i4b (0D_NOT_real -
      .def_property("i4b", &RadInt1Proxy::i4b, &RadInt1Proxy::set_i4b)
      // RadInt1Proxy.i4z (0D_NOT_real -
      .def_property("i4z", &RadInt1Proxy::i4z, &RadInt1Proxy::set_i4z)
      // RadInt1Proxy.i5a (0D_NOT_real -
      .def_property("i5a", &RadInt1Proxy::i5a, &RadInt1Proxy::set_i5a)
      // RadInt1Proxy.i5b (0D_NOT_real -
      .def_property("i5b", &RadInt1Proxy::i5b, &RadInt1Proxy::set_i5b)
      // RadInt1Proxy.i6b (0D_NOT_real -
      .def_property("i6b", &RadInt1Proxy::i6b, &RadInt1Proxy::set_i6b)
      // RadInt1Proxy.lin_i2_E4 (0D_NOT_real -
      .def_property(
          "lin_i2_E4", &RadInt1Proxy::lin_i2_E4, &RadInt1Proxy::set_lin_i2_E4)
      // RadInt1Proxy.lin_i3_E7 (0D_NOT_real -
      .def_property(
          "lin_i3_E7", &RadInt1Proxy::lin_i3_E7, &RadInt1Proxy::set_lin_i3_E7)
      // RadInt1Proxy.lin_i5a_E6 (0D_NOT_real -
      .def_property(
          "lin_i5a_E6",
          &RadInt1Proxy::lin_i5a_E6,
          &RadInt1Proxy::set_lin_i5a_E6)
      // RadInt1Proxy.lin_i5b_E6 (0D_NOT_real -
      .def_property(
          "lin_i5b_E6",
          &RadInt1Proxy::lin_i5b_E6,
          &RadInt1Proxy::set_lin_i5b_E6)
      // RadInt1Proxy.lin_norm_emit_a (0D_NOT_real - Running sum
      .def_property(
          "lin_norm_emit_a",
          &RadInt1Proxy::lin_norm_emit_a,
          &RadInt1Proxy::set_lin_norm_emit_a)
      // RadInt1Proxy.lin_norm_emit_b (0D_NOT_real - Running sum
      .def_property(
          "lin_norm_emit_b",
          &RadInt1Proxy::lin_norm_emit_b,
          &RadInt1Proxy::set_lin_norm_emit_b)
      // RadInt1Proxy.lin_sig_E (0D_NOT_real - Running sum
      .def_property(
          "lin_sig_E", &RadInt1Proxy::lin_sig_E, &RadInt1Proxy::set_lin_sig_E)
      // RadInt1Proxy.n_steps (0D_NOT_real - number of qromb steps needed
      .def_property(
          "n_steps", &RadInt1Proxy::n_steps, &RadInt1Proxy::set_n_steps)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return RadInt1ProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const RadInt1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadInt1Proxy& self) {
            return RadInt1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadInt1Proxy& self, py::dict& memo) {
            return RadInt1Proxy(self);
          })

      ;

  bind_FTypeArrayND<RadInt1ProxyArray1D>(m, "RadInt1StructArray1D");
  bind_FTypeAlloc1D<RadInt1ProxyAlloc1D>(m, "RadInt1StructAlloc1D");
  // 2D RadInt1Proxy arrays are not used in structs/routines
  // 3D RadInt1Proxy arrays are not used in structs/routines
}

// =============================================================================
// rad_int_all_ele_struct
void init_rad_int_all_ele_struct(
    py::module& m,
    py::class_<RadIntAllEleProxy>& cls) {
  cls.def(py::init<>())
      // RadIntAllEleProxy.branch (1D_ALLOC_type - Array is indexed from 0
      .def_property_readonly("branch", &RadIntAllEleProxy::branch)

      .def(
          "__repr__",
          [](const RadIntAllEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntAllEleProxy& self) {
            return RadIntAllEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadIntAllEleProxy& self, py::dict& memo) {
            return RadIntAllEleProxy(self);
          })

      ;

  // 1D RadIntAllEleProxy arrays are not used in structs/routines
  // 2D RadIntAllEleProxy arrays are not used in structs/routines
  // 3D RadIntAllEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_int_branch_struct
void init_rad_int_branch_struct(
    py::module& m,
    py::class_<RadIntBranchProxy>& cls) {
  cls.def(py::init<>())
      // RadIntBranchProxy.ele (1D_ALLOC_type - Array is indexed from 0
      .def_property_readonly("ele", &RadIntBranchProxy::ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return RadIntBranchProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RadIntBranchProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadIntBranchProxy& self) {
            return RadIntBranchProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadIntBranchProxy& self, py::dict& memo) {
            return RadIntBranchProxy(self);
          })

      ;

  bind_FTypeArrayND<RadIntBranchProxyArray1D>(m, "RadIntBranchStructArray1D");
  bind_FTypeAlloc1D<RadIntBranchProxyAlloc1D>(m, "RadIntBranchStructAlloc1D");
  // 2D RadIntBranchProxy arrays are not used in structs/routines
  // 3D RadIntBranchProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_map_ele_struct
void init_rad_map_ele_struct(py::module& m, py::class_<RadMapEleProxy>& cls) {
  cls.def(py::init<>())
      // RadMapEleProxy.rm0 (0D_NOT_type - Upstream half and downstream half matrices for an element.
      .def_property("rm0", &RadMapEleProxy::rm0, &RadMapEleProxy::set_rm0)
      // RadMapEleProxy.rm1 (0D_NOT_type - Upstream half and downstream half matrices for an element.
      .def_property("rm1", &RadMapEleProxy::rm1, &RadMapEleProxy::set_rm1)
      // RadMapEleProxy.stale (0D_NOT_logical -
      .def_property("stale", &RadMapEleProxy::stale, &RadMapEleProxy::set_stale)

      .def(
          "__repr__",
          [](const RadMapEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapEleProxy& self) {
            return RadMapEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadMapEleProxy& self, py::dict& memo) {
            return RadMapEleProxy(self);
          })

      ;

  // 1D RadMapEleProxy arrays are not used in structs/routines
  // 2D RadMapEleProxy arrays are not used in structs/routines
  // 3D RadMapEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rad_map_struct
void init_rad_map_struct(py::module& m, py::class_<RadMapProxy>& cls) {
  cls.def(py::init<>())
      // RadMapProxy.ref_orb (1D_NOT_real - Reference point around which damp_mat is calculated.
      .def_property_readonly("ref_orb", &RadMapProxy::ref_orb)
      // RadMapProxy.damp_dmat (2D_NOT_real - damp_correction = xfer_mat_with_damping - xfer_mat_without_damping.
      .def_property_readonly("damp_dmat", &RadMapProxy::damp_dmat)
      // RadMapProxy.xfer_damp_vec (1D_NOT_real - Transfer map with damping 0th order vector.
      .def_property_readonly("xfer_damp_vec", &RadMapProxy::xfer_damp_vec)
      // RadMapProxy.xfer_damp_mat (2D_NOT_real - 1st order matrix: xfer_no_damp_mat + xfer_damp_correction.
      .def_property_readonly("xfer_damp_mat", &RadMapProxy::xfer_damp_mat)
      // RadMapProxy.stoc_mat (2D_NOT_real - Stochastic variance or 'kick' (Cholesky decomposed) matrix.
      .def_property_readonly("stoc_mat", &RadMapProxy::stoc_mat)

      .def("__repr__", [](const RadMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RadMapProxy& self) {
            return RadMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RadMapProxy& self, py::dict& memo) {
            return RadMapProxy(self);
          })

      ;

  // 1D RadMapProxy arrays are not used in structs/routines
  // 2D RadMapProxy arrays are not used in structs/routines
  // 3D RadMapProxy arrays are not used in structs/routines
}

// =============================================================================
// ramper_lord_struct
void init_ramper_lord_struct(py::module& m, py::class_<RamperLordProxy>& cls) {
  cls.def(py::init<>())
      // RamperLordProxy.ix_ele (0D_NOT_integer - Lord index
      .def_property(
          "ix_ele", &RamperLordProxy::ix_ele, &RamperLordProxy::set_ix_ele)
      // RamperLordProxy.ix_con (0D_NOT_integer - Index in lord%control%ramp(:) array
      .def_property(
          "ix_con", &RamperLordProxy::ix_con, &RamperLordProxy::set_ix_con)
      // RamperLordProxy.attrib_ptr (0D_PTR_real - Pointer to attribute in this element.
      .def_property(
          "attrib_ptr",
          &RamperLordProxy::attrib_ptr,
          &RamperLordProxy::set_attrib_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return RamperLordProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RamperLordProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RamperLordProxy& self) {
            return RamperLordProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RamperLordProxy& self, py::dict& memo) {
            return RamperLordProxy(self);
          })

      ;

  bind_FTypeArrayND<RamperLordProxyArray1D>(m, "RamperLordStructArray1D");
  bind_FTypeAlloc1D<RamperLordProxyAlloc1D>(m, "RamperLordStructAlloc1D");
  // 2D RamperLordProxy arrays are not used in structs/routines
  // 3D RamperLordProxy arrays are not used in structs/routines
}

// =============================================================================
// resonance_h_struct
void init_resonance_h_struct(py::module& m, py::class_<ResonanceHProxy>& cls) {
  cls.def(py::init<>())
      // ResonanceHProxy.id (0D_NOT_character - 6 digit ID. EG: '003100'
      .def_property("id", &ResonanceHProxy::id, &ResonanceHProxy::set_id)
      // ResonanceHProxy.c_val (0D_NOT_complex - Resonance value
      .def_property(
          "c_val", &ResonanceHProxy::c_val, &ResonanceHProxy::set_c_val)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ResonanceHProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ResonanceHProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ResonanceHProxy& self) {
            return ResonanceHProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ResonanceHProxy& self, py::dict& memo) {
            return ResonanceHProxy(self);
          })

      ;

  bind_FTypeArrayND<ResonanceHProxyArray1D>(m, "ResonanceHStructArray1D");
  bind_FTypeAlloc1D<ResonanceHProxyAlloc1D>(m, "ResonanceHStructAlloc1D");
  // 2D ResonanceHProxy arrays are not used in structs/routines
  // 3D ResonanceHProxy arrays are not used in structs/routines
}

// =============================================================================
// rf_ele_struct
void init_rf_ele_struct(py::module& m, py::class_<RfEleProxy>& cls) {
  cls.def(py::init<>())
      // RfEleProxy.steps (1D_ALLOC_type - Energy stair step array indexed from zero.
      .def_property_readonly("steps", &RfEleProxy::steps)
      // RfEleProxy.ds_step (0D_NOT_real - length of a stair step.
      .def_property("ds_step", &RfEleProxy::ds_step, &RfEleProxy::set_ds_step)

      .def("__repr__", [](const RfEleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfEleProxy& self) {
            return RfEleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RfEleProxy& self, py::dict& memo) {
            return RfEleProxy(self);
          })

      ;

  // 1D RfEleProxy arrays are not used in structs/routines
  // 2D RfEleProxy arrays are not used in structs/routines
  // 3D RfEleProxy arrays are not used in structs/routines
}

// =============================================================================
// rf_stair_step_struct
void init_rf_stair_step_struct(
    py::module& m,
    py::class_<RfStairStepProxy>& cls) {
  cls.def(py::init<>())
      // RfStairStepProxy.E_tot0 (0D_NOT_real - Reference energy in the drift region (before the kick point).
      .def_property(
          "E_tot0", &RfStairStepProxy::E_tot0, &RfStairStepProxy::set_E_tot0)
      // RfStairStepProxy.E_tot1 (0D_NOT_real - Reference energy after the kick point.
      .def_property(
          "E_tot1", &RfStairStepProxy::E_tot1, &RfStairStepProxy::set_E_tot1)
      // RfStairStepProxy.p0c (0D_NOT_real - Reference momentum in the drift region (before the kick point).
      .def_property("p0c", &RfStairStepProxy::p0c, &RfStairStepProxy::set_p0c)
      // RfStairStepProxy.p1c (0D_NOT_real - Reference momentum after the kick point.
      .def_property("p1c", &RfStairStepProxy::p1c, &RfStairStepProxy::set_p1c)
      // RfStairStepProxy.scale (0D_NOT_real - Scale for multipole kick at the kick point. Sum over all steps will be 1.
      .def_property(
          "scale", &RfStairStepProxy::scale, &RfStairStepProxy::set_scale)
      // RfStairStepProxy.time (0D_NOT_real - Reference particle time at the kick point with respect to beginning of element.
      .def_property(
          "time", &RfStairStepProxy::time, &RfStairStepProxy::set_time)
      // RfStairStepProxy.s0 (0D_NOT_real - S-position at beginning of drift region relative to the beginning of the element.
      .def_property("s0", &RfStairStepProxy::s0, &RfStairStepProxy::set_s0)
      // RfStairStepProxy.s (0D_NOT_real - S-position at the kick point relative to the beginning of the element.
      .def_property("s", &RfStairStepProxy::s, &RfStairStepProxy::set_s)
      // RfStairStepProxy.ix_step (0D_NOT_integer - Step index in ele%rf%steps(:) array
      .def_property(
          "ix_step", &RfStairStepProxy::ix_step, &RfStairStepProxy::set_ix_step)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return RfStairStepProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const RfStairStepProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RfStairStepProxy& self) {
            return RfStairStepProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RfStairStepProxy& self, py::dict& memo) {
            return RfStairStepProxy(self);
          })

      ;

  bind_FTypeArrayND<RfStairStepProxyArray1D>(m, "RfStairStepStructArray1D");
  bind_FTypeAlloc1D<RfStairStepProxyAlloc1D>(m, "RfStairStepStructAlloc1D");
  // 2D RfStairStepProxy arrays are not used in structs/routines
  // 3D RfStairStepProxy arrays are not used in structs/routines
}

// =============================================================================
// random_state_struct
void init_random_state_struct(
    py::module& m,
    py::class_<RandomStateProxy>& cls) {
  cls.def(py::init<>())
      // RandomStateProxy.ix (0D_NOT_integer8 -
      .def_property("ix", &RandomStateProxy::ix, &RandomStateProxy::set_ix)
      // RandomStateProxy.iy (0D_NOT_integer8 -
      .def_property("iy", &RandomStateProxy::iy, &RandomStateProxy::set_iy)
      // RandomStateProxy.number_stored (0D_NOT_logical -
      .def_property(
          "number_stored",
          &RandomStateProxy::number_stored,
          &RandomStateProxy::set_number_stored)
      // RandomStateProxy.h_saved (0D_NOT_real -
      .def_property(
          "h_saved", &RandomStateProxy::h_saved, &RandomStateProxy::set_h_saved)
      // RandomStateProxy.engine (0D_NOT_integer - Params
      .def_property(
          "engine", &RandomStateProxy::engine, &RandomStateProxy::set_engine)
      // RandomStateProxy.seed (0D_NOT_integer -
      .def_property(
          "seed", &RandomStateProxy::seed, &RandomStateProxy::set_seed)
      // RandomStateProxy.am (0D_NOT_real -
      .def_property("am", &RandomStateProxy::am, &RandomStateProxy::set_am)
      // RandomStateProxy.gauss_converter (0D_NOT_integer -
      .def_property(
          "gauss_converter",
          &RandomStateProxy::gauss_converter,
          &RandomStateProxy::set_gauss_converter)
      // RandomStateProxy.gauss_sigma_cut (0D_NOT_real - Only used if positive.
      .def_property(
          "gauss_sigma_cut",
          &RandomStateProxy::gauss_sigma_cut,
          &RandomStateProxy::set_gauss_sigma_cut)
      // RandomStateProxy.in_sobseq (0D_NOT_integer8 -
      .def_property(
          "in_sobseq",
          &RandomStateProxy::in_sobseq,
          &RandomStateProxy::set_in_sobseq)
      // 1D_NOT_integer8 ix_sobseq proxy support missing
      // RandomStateProxy.x_sobseq (1D_NOT_real -
      .def_property_readonly("x_sobseq", &RandomStateProxy::x_sobseq)

      .def(
          "__repr__",
          [](const RandomStateProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const RandomStateProxy& self) {
            return RandomStateProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const RandomStateProxy& self, py::dict& memo) {
            return RandomStateProxy(self);
          })

      ;

  // 1D RandomStateProxy arrays are not used in structs/routines
  // 2D RandomStateProxy arrays are not used in structs/routines
  // 3D RandomStateProxy arrays are not used in structs/routines
}