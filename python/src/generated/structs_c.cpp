#include "pybmad/generated/structs_c.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// cartesian_map_struct
void init_cartesian_map_struct(py::module &m, py::class_<CartesianMapStruct> &cls) {
  cls.def(py::init<>())
      // CartesianMapStruct.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CartesianMapStruct::field_scale,
          &CartesianMapStruct::set_field_scale
      )
      // CartesianMapStruct.r0 (1D_NOT_real - Field origin offset.
      .def_property("r0", &CartesianMapStruct::r0, &CartesianMapStruct::set_r0)
      // CartesianMapStruct.master_parameter (0D_NOT_integer - Master parameter in ele%value(:)
      // array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CartesianMapStruct::master_parameter,
          &CartesianMapStruct::set_master_parameter
      )
      // CartesianMapStruct.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or
      // anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CartesianMapStruct::ele_anchor_pt,
          &CartesianMapStruct::set_ele_anchor_pt
      )
      // CartesianMapStruct.field_type (0D_NOT_integer - or electric$
      .def_property(
          "field_type",
          &CartesianMapStruct::field_type,
          &CartesianMapStruct::set_field_type
      )
      // CartesianMapStruct.ptr (0D_PTR_type -
      .def_property("ptr", &CartesianMapStruct::ptr, &CartesianMapStruct::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CartesianMapStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CartesianMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapStruct &self) {
            return CartesianMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapStruct &self, py::dict &memo) { return CartesianMapStruct(self); }
      )

      ;

  bind_FTypeArrayND<CartesianMapStructArray1D>(m, "CartesianMapStructArray1D");
  bind_FTypeAlloc1D<CartesianMapStructAlloc1D>(m, "CartesianMapStructAlloc1D");
  // 2D CartesianMapStruct arrays are not used in structs/routines
  // 3D CartesianMapStruct arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term1_struct
void init_cartesian_map_term1_struct(py::module &m, py::class_<CartesianMapTerm1Struct> &cls) {
  cls.def(py::init<>())
      // CartesianMapTerm1Struct.coef (0D_NOT_real -
      .def_property("coef", &CartesianMapTerm1Struct::coef, &CartesianMapTerm1Struct::set_coef)
      // CartesianMapTerm1Struct.kx (0D_NOT_real -
      .def_property("kx", &CartesianMapTerm1Struct::kx, &CartesianMapTerm1Struct::set_kx)
      // CartesianMapTerm1Struct.ky (0D_NOT_real -
      .def_property("ky", &CartesianMapTerm1Struct::ky, &CartesianMapTerm1Struct::set_ky)
      // CartesianMapTerm1Struct.kz (0D_NOT_real -
      .def_property("kz", &CartesianMapTerm1Struct::kz, &CartesianMapTerm1Struct::set_kz)
      // CartesianMapTerm1Struct.x0 (0D_NOT_real -
      .def_property("x0", &CartesianMapTerm1Struct::x0, &CartesianMapTerm1Struct::set_x0)
      // CartesianMapTerm1Struct.y0 (0D_NOT_real -
      .def_property("y0", &CartesianMapTerm1Struct::y0, &CartesianMapTerm1Struct::set_y0)
      // CartesianMapTerm1Struct.phi_z (0D_NOT_real -
      .def_property("phi_z", &CartesianMapTerm1Struct::phi_z, &CartesianMapTerm1Struct::set_phi_z)
      // CartesianMapTerm1Struct.family (0D_NOT_integer - family_x$, etc.
      .def_property(
          "family",
          &CartesianMapTerm1Struct::family,
          &CartesianMapTerm1Struct::set_family
      )
      // CartesianMapTerm1Struct.form (0D_NOT_integer - hyper_y$, etc.
      .def_property("form", &CartesianMapTerm1Struct::form, &CartesianMapTerm1Struct::set_form)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CartesianMapTerm1StructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CartesianMapTerm1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTerm1Struct &self) {
            return CartesianMapTerm1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapTerm1Struct &self, py::dict &memo) {
            return CartesianMapTerm1Struct(self);
          }
      )

      ;

  bind_FTypeArrayND<CartesianMapTerm1StructArray1D>(m, "CartesianMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CartesianMapTerm1StructAlloc1D>(m, "CartesianMapTerm1StructAlloc1D");
  // 2D CartesianMapTerm1Struct arrays are not used in structs/routines
  // 3D CartesianMapTerm1Struct arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term_struct
void init_cartesian_map_term_struct(py::module &m, py::class_<CartesianMapTermStruct> &cls) {
  cls.def(py::init<>())
      // CartesianMapTermStruct.file (0D_NOT_character - Input file name. Used also as ID for
      // instances.
      .def_property("file", &CartesianMapTermStruct::file, &CartesianMapTermStruct::set_file)
      // CartesianMapTermStruct.n_link (0D_NOT_integer - For memory management of %term
      .def_property("n_link", &CartesianMapTermStruct::n_link, &CartesianMapTermStruct::set_n_link)
      // CartesianMapTermStruct.term (1D_ALLOC_type -
      .def_property_readonly("term", &CartesianMapTermStruct::term)

      .def("__repr__", [](const CartesianMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTermStruct &self) {
            return CartesianMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapTermStruct &self, py::dict &memo) {
            return CartesianMapTermStruct(self);
          }
      )

      ;

  // 1D CartesianMapTermStruct arrays are not used in structs/routines
  // 2D CartesianMapTermStruct arrays are not used in structs/routines
  // 3D CartesianMapTermStruct arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_struct
void init_complex_taylor_struct(py::module &m, py::class_<ComplexTaylorStruct> &cls) {
  cls.def(py::init<>())
      // ComplexTaylorStruct.ref (0D_NOT_complex -
      .def_property("ref", &ComplexTaylorStruct::ref, &ComplexTaylorStruct::set_ref)
      // ComplexTaylorStruct.term (1D_PTR_type -
      .def_property_readonly("term", &ComplexTaylorStruct::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ComplexTaylorStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ComplexTaylorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorStruct &self) {
            return ComplexTaylorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ComplexTaylorStruct &self, py::dict &memo) { return ComplexTaylorStruct(self); }
      )

      ;

  bind_FTypeArrayND<ComplexTaylorStructArray1D>(m, "ComplexTaylorStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorStructAlloc1D>(m, "ComplexTaylorStructAlloc1D");
  // 2D ComplexTaylorStruct arrays are not used in structs/routines
  // 3D ComplexTaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_term_struct
void init_complex_taylor_term_struct(py::module &m, py::class_<ComplexTaylorTermStruct> &cls) {
  cls.def(py::init<>())
      // ComplexTaylorTermStruct.coef (0D_NOT_complex -
      .def_property("coef", &ComplexTaylorTermStruct::coef, &ComplexTaylorTermStruct::set_coef)
      // ComplexTaylorTermStruct.expn (1D_NOT_integer -
      .def_property("expn", &ComplexTaylorTermStruct::expn, &ComplexTaylorTermStruct::set_expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ComplexTaylorTermStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ComplexTaylorTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorTermStruct &self) {
            return ComplexTaylorTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ComplexTaylorTermStruct &self, py::dict &memo) {
            return ComplexTaylorTermStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<ComplexTaylorTermStructArray1D>(m, "ComplexTaylorTermStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorTermStructAlloc1D>(m, "ComplexTaylorTermStructAlloc1D");
  // 2D ComplexTaylorTermStruct arrays are not used in structs/routines
  // 3D ComplexTaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// control_ramp1_struct
void init_control_ramp1_struct(py::module &m, py::class_<ControlRamp1Struct> &cls) {
  cls.def(py::init<>())
      // ControlRamp1Struct.y_knot (1D_ALLOC_real -
      .def_property("y_knot", &ControlRamp1Struct::y_knot, &ControlRamp1Struct::set_y_knot)
      // ControlRamp1Struct.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlRamp1Struct::stack)
      // ControlRamp1Struct.attribute (0D_NOT_character - Name of attribute controlled. Set to
      // 'FIELD_OVERLAPS' for field overlaps.
      .def_property("attribute", &ControlRamp1Struct::attribute, &ControlRamp1Struct::set_attribute)
      // ControlRamp1Struct.slave_name (0D_NOT_character - Name of slave.
      .def_property(
          "slave_name",
          &ControlRamp1Struct::slave_name,
          &ControlRamp1Struct::set_slave_name
      )
      // ControlRamp1Struct.is_controller (0D_NOT_logical - Is the slave a controller? If so
      // bookkeeping is different.
      .def_property(
          "is_controller",
          &ControlRamp1Struct::is_controller,
          &ControlRamp1Struct::set_is_controller
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ControlRamp1StructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ControlRamp1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlRamp1Struct &self) {
            return ControlRamp1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlRamp1Struct &self, py::dict &memo) { return ControlRamp1Struct(self); }
      )

      ;

  bind_FTypeArrayND<ControlRamp1StructArray1D>(m, "ControlRamp1StructArray1D");
  bind_FTypeAlloc1D<ControlRamp1StructAlloc1D>(m, "ControlRamp1StructAlloc1D");
  // 2D ControlRamp1Struct arrays are not used in structs/routines
  // 3D ControlRamp1Struct arrays are not used in structs/routines
}

// =============================================================================
// control_struct
void init_control_struct(py::module &m, py::class_<ControlStruct> &cls) {
  cls.def(py::init<>())
      // ControlStruct.value (0D_NOT_real - Used by group, and overlay elements.
      .def_property("value", &ControlStruct::value, &ControlStruct::set_value)
      // ControlStruct.y_knot (1D_ALLOC_real -
      .def_property("y_knot", &ControlStruct::y_knot, &ControlStruct::set_y_knot)
      // ControlStruct.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlStruct::stack)
      // ControlStruct.slave (0D_NOT_type -
      .def_property("slave", &ControlStruct::slave, &ControlStruct::set_slave)
      // ControlStruct.lord (0D_NOT_type -
      .def_property("lord", &ControlStruct::lord, &ControlStruct::set_lord)
      // ControlStruct.slave_name (0D_NOT_character - Name of slave.
      .def_property("slave_name", &ControlStruct::slave_name, &ControlStruct::set_slave_name)
      // ControlStruct.attribute (0D_NOT_character - Name of attribute controlled. Set to
      // 'FIELD_OVERLAPS' for field overlaps. Set to 'INPUT' or 'OUTPUT' for feedback slaves.
      .def_property("attribute", &ControlStruct::attribute, &ControlStruct::set_attribute)
      // ControlStruct.ix_attrib (0D_NOT_integer - Index of attribute controlled. See note above!
      .def_property("ix_attrib", &ControlStruct::ix_attrib, &ControlStruct::set_ix_attrib)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ControlStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ControlStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlStruct &self) {
            return ControlStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlStruct &self, py::dict &memo) { return ControlStruct(self); }
      )

      ;

  bind_FTypeArrayND<ControlStructArray1D>(m, "ControlStructArray1D");
  bind_FTypeAlloc1D<ControlStructAlloc1D>(m, "ControlStructAlloc1D");
  // 2D ControlStruct arrays are not used in structs/routines
  // 3D ControlStruct arrays are not used in structs/routines
}

// =============================================================================
// control_var1_struct
void init_control_var1_struct(py::module &m, py::class_<ControlVar1Struct> &cls) {
  cls.def(py::init<>())
      // ControlVar1Struct.name (0D_NOT_character -
      .def_property("name", &ControlVar1Struct::name, &ControlVar1Struct::set_name)
      // ControlVar1Struct.value (0D_NOT_real -
      .def_property("value", &ControlVar1Struct::value, &ControlVar1Struct::set_value)
      // ControlVar1Struct.old_value (0D_NOT_real -
      .def_property("old_value", &ControlVar1Struct::old_value, &ControlVar1Struct::set_old_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ControlVar1StructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ControlVar1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlVar1Struct &self) {
            return ControlVar1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlVar1Struct &self, py::dict &memo) { return ControlVar1Struct(self); }
      )

      ;

  bind_FTypeArrayND<ControlVar1StructArray1D>(m, "ControlVar1StructArray1D");
  bind_FTypeAlloc1D<ControlVar1StructAlloc1D>(m, "ControlVar1StructAlloc1D");
  // 2D ControlVar1Struct arrays are not used in structs/routines
  // 3D ControlVar1Struct arrays are not used in structs/routines
}

// =============================================================================
// controller_struct
void init_controller_struct(py::module &m, py::class_<ControllerStruct> &cls) {
  cls.def(py::init<>())
      // ControllerStruct.var (1D_ALLOC_type -
      .def_property_readonly("var", &ControllerStruct::var)
      // ControllerStruct.ramp (1D_ALLOC_type - For ramper lord elements
      .def_property_readonly("ramp", &ControllerStruct::ramp)
      // ControllerStruct.ramper_lord (1D_ALLOC_type - Ramper lord info for this slave
      .def_property_readonly("ramper_lord", &ControllerStruct::ramper_lord)
      // ControllerStruct.x_knot (1D_ALLOC_real -
      .def_property("x_knot", &ControllerStruct::x_knot, &ControllerStruct::set_x_knot)

      .def("__repr__", [](const ControllerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControllerStruct &self) {
            return ControllerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControllerStruct &self, py::dict &memo) { return ControllerStruct(self); }
      )

      ;

  // 1D ControllerStruct arrays are not used in structs/routines
  // 2D ControllerStruct arrays are not used in structs/routines
  // 3D ControllerStruct arrays are not used in structs/routines
}

// =============================================================================
// coord_array_struct
void init_coord_array_struct(py::module &m, py::class_<CoordArrayStruct> &cls) {
  cls.def(py::init<>())
      // CoordArrayStruct.orbit (1D_ALLOC_type -
      .def_property_readonly("orbit", &CoordArrayStruct::orbit)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordArrayStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CoordArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordArrayStruct &self) {
            return CoordArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CoordArrayStruct &self, py::dict &memo) { return CoordArrayStruct(self); }
      )

      ;

  bind_FTypeArrayND<CoordArrayStructArray1D>(m, "CoordArrayStructArray1D");
  bind_FTypeAlloc1D<CoordArrayStructAlloc1D>(m, "CoordArrayStructAlloc1D");
  // 2D CoordArrayStruct arrays are not used in structs/routines
  // 3D CoordArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// coord_struct
void init_coord_struct(py::module &m, py::class_<CoordStruct> &cls) {
  cls.def(py::init<>())
      // CoordStruct.vec (1D_NOT_real - (x, px, y, py, z, pz). Generally phase space for charged
      // particles. See Bmad manual.
      .def_property("vec", &CoordStruct::vec, &CoordStruct::set_vec)
      // CoordStruct.s (0D_NOT_real - Longitudinal position
      .def_property("s", &CoordStruct::s, &CoordStruct::set_s)
      // CoordStruct.t (0D_NOT_real16 - Absolute time (not relative to reference). Note: Quad
      // precision!
      .def_property("t", &CoordStruct::t, &CoordStruct::set_t)
      // CoordStruct.spin (1D_NOT_real - Spin.
      .def_property("spin", &CoordStruct::spin, &CoordStruct::set_spin)
      // CoordStruct.field (1D_NOT_real - Photon E-field intensity (x,y).
      .def_property("field", &CoordStruct::field, &CoordStruct::set_field)
      // CoordStruct.phase (1D_NOT_real - Photon E-field phase (x,y). For charged particles,
      // phase(1) is RF phase.
      .def_property("phase", &CoordStruct::phase, &CoordStruct::set_phase)
      // CoordStruct.charge (0D_NOT_real - Macroparticle weight (which is different from particle
      // species charge). For some space charge calcs the weight is in Coulombs.
      .def_property("charge", &CoordStruct::charge, &CoordStruct::set_charge)
      // CoordStruct.dt_ref (0D_NOT_real - Used in: * time tracking for computing z. * by coherent
      // photons = path_length/c_light.
      .def_property("dt_ref", &CoordStruct::dt_ref, &CoordStruct::set_dt_ref)
      // CoordStruct.r (0D_NOT_real - For general use. Not used by Bmad.
      .def_property("r", &CoordStruct::r, &CoordStruct::set_r)
      // CoordStruct.p0c (0D_NOT_real - For non-photons: Reference momentum. For photons: Photon
      // momentum (not reference).
      .def_property("p0c", &CoordStruct::p0c, &CoordStruct::set_p0c)
      // CoordStruct.E_potential (0D_NOT_real - Potential energy.
      .def_property("E_potential", &CoordStruct::E_potential, &CoordStruct::set_E_potential)
      // CoordStruct.beta (0D_NOT_real - Velocity / c_light.
      .def_property("beta", &CoordStruct::beta, &CoordStruct::set_beta)
      // CoordStruct.ix_ele (0D_NOT_integer - Index of the lattice element the particle is in. May
      // be -1 if element is not associated with a lattice.
      .def_property("ix_ele", &CoordStruct::ix_ele, &CoordStruct::set_ix_ele)
      // CoordStruct.ix_branch (0D_NOT_integer - Index of the lattice branch the particle is in.
      .def_property("ix_branch", &CoordStruct::ix_branch, &CoordStruct::set_ix_branch)
      // CoordStruct.ix_turn (0D_NOT_integer - Turn index for multiturn tracking.
      .def_property("ix_turn", &CoordStruct::ix_turn, &CoordStruct::set_ix_turn)
      // CoordStruct.ix_user (0D_NOT_integer - For general use, not used by Bmad.
      .def_property("ix_user", &CoordStruct::ix_user, &CoordStruct::set_ix_user)
      // CoordStruct.state (0D_NOT_integer - alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc.
      .def_property("state", &CoordStruct::state, &CoordStruct::set_state)
      // CoordStruct.direction (0D_NOT_integer - +1 or -1. Sign of longitudinal direction of motion
      // (ds/dt). This is independent of the element orientation.
      .def_property("direction", &CoordStruct::direction, &CoordStruct::set_direction)
      // CoordStruct.time_dir (0D_NOT_integer - +1 or -1. Time direction. -1 => Traveling backwards
      // in time.
      .def_property("time_dir", &CoordStruct::time_dir, &CoordStruct::set_time_dir)
      // CoordStruct.species (0D_NOT_integer - positron$, proton$, etc.
      .def_property("species", &CoordStruct::species, &CoordStruct::set_species)
      // CoordStruct.location (0D_NOT_integer - upstream_end$, inside$, or downstream_end$
      .def_property("location", &CoordStruct::location, &CoordStruct::set_location)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CoordStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordStruct &self) {
            return CoordStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CoordStruct &self, py::dict &memo) { return CoordStruct(self); }
      )

      ;

  bind_FTypeArrayND<CoordStructArray1D>(m, "CoordStructArray1D");
  bind_FTypeAlloc1D<CoordStructAlloc1D>(m, "CoordStructAlloc1D");
  // 2D CoordStruct arrays are not used in structs/routines
  // 3D CoordStruct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_struct
void init_cylindrical_map_struct(py::module &m, py::class_<CylindricalMapStruct> &cls) {
  cls.def(py::init<>())
      // CylindricalMapStruct.m (0D_NOT_integer - Azimuthal Mode: varies as cos(m*phi -
      // theta0_azimuth)
      .def_property("m", &CylindricalMapStruct::m, &CylindricalMapStruct::set_m)
      // CylindricalMapStruct.harmonic (0D_NOT_integer - Harmonic of fundamental
      .def_property(
          "harmonic",
          &CylindricalMapStruct::harmonic,
          &CylindricalMapStruct::set_harmonic
      )
      // CylindricalMapStruct.phi0_fieldmap (0D_NOT_real - Mode oscillates as: twopi * (f * t +
      // phi0_fieldmap)
      .def_property(
          "phi0_fieldmap",
          &CylindricalMapStruct::phi0_fieldmap,
          &CylindricalMapStruct::set_phi0_fieldmap
      )
      // CylindricalMapStruct.theta0_azimuth (0D_NOT_real - Azimuthal ((x, y) plane) orientation of
      // mode.
      .def_property(
          "theta0_azimuth",
          &CylindricalMapStruct::theta0_azimuth,
          &CylindricalMapStruct::set_theta0_azimuth
      )
      // CylindricalMapStruct.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CylindricalMapStruct::field_scale,
          &CylindricalMapStruct::set_field_scale
      )
      // CylindricalMapStruct.master_parameter (0D_NOT_integer - Master parameter in ele%value(:)
      // array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CylindricalMapStruct::master_parameter,
          &CylindricalMapStruct::set_master_parameter
      )
      // CylindricalMapStruct.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or
      // anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CylindricalMapStruct::ele_anchor_pt,
          &CylindricalMapStruct::set_ele_anchor_pt
      )
      // CylindricalMapStruct.dz (0D_NOT_real - Distance between sampled field points.
      .def_property("dz", &CylindricalMapStruct::dz, &CylindricalMapStruct::set_dz)
      // CylindricalMapStruct.r0 (1D_NOT_real - Field origin offset.
      .def_property("r0", &CylindricalMapStruct::r0, &CylindricalMapStruct::set_r0)
      // CylindricalMapStruct.ptr (0D_PTR_type -
      .def_property("ptr", &CylindricalMapStruct::ptr, &CylindricalMapStruct::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CylindricalMapStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CylindricalMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapStruct &self) {
            return CylindricalMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapStruct &self, py::dict &memo) {
            return CylindricalMapStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<CylindricalMapStructArray1D>(m, "CylindricalMapStructArray1D");
  bind_FTypeAlloc1D<CylindricalMapStructAlloc1D>(m, "CylindricalMapStructAlloc1D");
  // 2D CylindricalMapStruct arrays are not used in structs/routines
  // 3D CylindricalMapStruct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term1_struct
void init_cylindrical_map_term1_struct(py::module &m, py::class_<CylindricalMapTerm1Struct> &cls) {
  cls.def(py::init<>())
      // CylindricalMapTerm1Struct.e_coef (0D_NOT_complex -
      .def_property(
          "e_coef",
          &CylindricalMapTerm1Struct::e_coef,
          &CylindricalMapTerm1Struct::set_e_coef
      )
      // CylindricalMapTerm1Struct.b_coef (0D_NOT_complex -
      .def_property(
          "b_coef",
          &CylindricalMapTerm1Struct::b_coef,
          &CylindricalMapTerm1Struct::set_b_coef
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CylindricalMapTerm1StructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const CylindricalMapTerm1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTerm1Struct &self) {
            return CylindricalMapTerm1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapTerm1Struct &self, py::dict &memo) {
            return CylindricalMapTerm1Struct(self);
          }
      )

      ;

  bind_FTypeArrayND<CylindricalMapTerm1StructArray1D>(m, "CylindricalMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CylindricalMapTerm1StructAlloc1D>(m, "CylindricalMapTerm1StructAlloc1D");
  // 2D CylindricalMapTerm1Struct arrays are not used in structs/routines
  // 3D CylindricalMapTerm1Struct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term_struct
void init_cylindrical_map_term_struct(py::module &m, py::class_<CylindricalMapTermStruct> &cls) {
  cls.def(py::init<>())
      // CylindricalMapTermStruct.file (0D_NOT_character - Input file name. Used also as ID for
      // instances.
      .def_property("file", &CylindricalMapTermStruct::file, &CylindricalMapTermStruct::set_file)
      // CylindricalMapTermStruct.n_link (0D_NOT_integer - For memory management of this structure
      .def_property(
          "n_link",
          &CylindricalMapTermStruct::n_link,
          &CylindricalMapTermStruct::set_n_link
      )
      // CylindricalMapTermStruct.term (1D_ALLOC_type -
      .def_property_readonly("term", &CylindricalMapTermStruct::term)

      .def("__repr__", [](const CylindricalMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTermStruct &self) {
            return CylindricalMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapTermStruct &self, py::dict &memo) {
            return CylindricalMapTermStruct(self);
          }
      )

      ;

  // 1D CylindricalMapTermStruct arrays are not used in structs/routines
  // 2D CylindricalMapTermStruct arrays are not used in structs/routines
  // 3D CylindricalMapTermStruct arrays are not used in structs/routines
}