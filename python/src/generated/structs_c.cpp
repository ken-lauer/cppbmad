#include "pybmad/generated/structs_c.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// cartesian_map_struct
void init_cartesian_map_struct(
    py::module& m,
    py::class_<CartesianMapProxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CartesianMapProxy::field_scale,
          &CartesianMapProxy::set_field_scale)
      // CartesianMapProxy.r0 (1D_NOT_real - Field origin offset.
      .def_property_readonly("r0", &CartesianMapProxy::r0)
      // CartesianMapProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CartesianMapProxy::master_parameter,
          &CartesianMapProxy::set_master_parameter)
      // CartesianMapProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CartesianMapProxy::ele_anchor_pt,
          &CartesianMapProxy::set_ele_anchor_pt)
      // CartesianMapProxy.field_type (0D_NOT_integer - or electric$
      .def_property(
          "field_type",
          &CartesianMapProxy::field_type,
          &CartesianMapProxy::set_field_type)
      // CartesianMapProxy.ptr (0D_PTR_type -
      .def_property("ptr", &CartesianMapProxy::ptr, &CartesianMapProxy::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CartesianMapProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CartesianMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapProxy& self) {
            return CartesianMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapProxy& self, py::dict& memo) {
            return CartesianMapProxy(self);
          })

      ;

  bind_FTypeArrayND<CartesianMapProxyArray1D>(m, "CartesianMapStructArray1D");
  bind_FTypeAlloc1D<CartesianMapProxyAlloc1D>(m, "CartesianMapStructAlloc1D");
  // 2D CartesianMapProxy arrays are not used in structs/routines
  // 3D CartesianMapProxy arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term1_struct
void init_cartesian_map_term1_struct(
    py::module& m,
    py::class_<CartesianMapTerm1Proxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapTerm1Proxy.coef (0D_NOT_real -
      .def_property(
          "coef",
          &CartesianMapTerm1Proxy::coef,
          &CartesianMapTerm1Proxy::set_coef)
      // CartesianMapTerm1Proxy.kx (0D_NOT_real -
      .def_property(
          "kx", &CartesianMapTerm1Proxy::kx, &CartesianMapTerm1Proxy::set_kx)
      // CartesianMapTerm1Proxy.ky (0D_NOT_real -
      .def_property(
          "ky", &CartesianMapTerm1Proxy::ky, &CartesianMapTerm1Proxy::set_ky)
      // CartesianMapTerm1Proxy.kz (0D_NOT_real -
      .def_property(
          "kz", &CartesianMapTerm1Proxy::kz, &CartesianMapTerm1Proxy::set_kz)
      // CartesianMapTerm1Proxy.x0 (0D_NOT_real -
      .def_property(
          "x0", &CartesianMapTerm1Proxy::x0, &CartesianMapTerm1Proxy::set_x0)
      // CartesianMapTerm1Proxy.y0 (0D_NOT_real -
      .def_property(
          "y0", &CartesianMapTerm1Proxy::y0, &CartesianMapTerm1Proxy::set_y0)
      // CartesianMapTerm1Proxy.phi_z (0D_NOT_real -
      .def_property(
          "phi_z",
          &CartesianMapTerm1Proxy::phi_z,
          &CartesianMapTerm1Proxy::set_phi_z)
      // CartesianMapTerm1Proxy.family (0D_NOT_integer - family_x$, etc.
      .def_property(
          "family",
          &CartesianMapTerm1Proxy::family,
          &CartesianMapTerm1Proxy::set_family)
      // CartesianMapTerm1Proxy.form (0D_NOT_integer - hyper_y$, etc.
      .def_property(
          "form",
          &CartesianMapTerm1Proxy::form,
          &CartesianMapTerm1Proxy::set_form)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CartesianMapTerm1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CartesianMapTerm1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTerm1Proxy& self) {
            return CartesianMapTerm1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapTerm1Proxy& self, py::dict& memo) {
            return CartesianMapTerm1Proxy(self);
          })

      ;

  bind_FTypeArrayND<CartesianMapTerm1ProxyArray1D>(
      m, "CartesianMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CartesianMapTerm1ProxyAlloc1D>(
      m, "CartesianMapTerm1StructAlloc1D");
  // 2D CartesianMapTerm1Proxy arrays are not used in structs/routines
  // 3D CartesianMapTerm1Proxy arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term_struct
void init_cartesian_map_term_struct(
    py::module& m,
    py::class_<CartesianMapTermProxy>& cls) {
  cls.def(py::init<>())
      // CartesianMapTermProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property(
          "file",
          &CartesianMapTermProxy::file,
          &CartesianMapTermProxy::set_file)
      // CartesianMapTermProxy.n_link (0D_NOT_integer - For memory management of %term
      .def_property(
          "n_link",
          &CartesianMapTermProxy::n_link,
          &CartesianMapTermProxy::set_n_link)
      // CartesianMapTermProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &CartesianMapTermProxy::term)

      .def(
          "__repr__",
          [](const CartesianMapTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTermProxy& self) {
            return CartesianMapTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CartesianMapTermProxy& self, py::dict& memo) {
            return CartesianMapTermProxy(self);
          })

      ;

  // 1D CartesianMapTermProxy arrays are not used in structs/routines
  // 2D CartesianMapTermProxy arrays are not used in structs/routines
  // 3D CartesianMapTermProxy arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_struct
void init_complex_taylor_struct(
    py::module& m,
    py::class_<ComplexTaylorProxy>& cls) {
  cls.def(py::init<>())
      // ComplexTaylorProxy.ref (0D_NOT_complex -
      .def_property(
          "ref", &ComplexTaylorProxy::ref, &ComplexTaylorProxy::set_ref)
      // ComplexTaylorProxy.term (1D_PTR_type -
      .def_property_readonly("term", &ComplexTaylorProxy::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ComplexTaylorProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ComplexTaylorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorProxy& self) {
            return ComplexTaylorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ComplexTaylorProxy& self, py::dict& memo) {
            return ComplexTaylorProxy(self);
          })

      ;

  bind_FTypeArrayND<ComplexTaylorProxyArray1D>(m, "ComplexTaylorStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorProxyAlloc1D>(m, "ComplexTaylorStructAlloc1D");
  // 2D ComplexTaylorProxy arrays are not used in structs/routines
  // 3D ComplexTaylorProxy arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_term_struct
void init_complex_taylor_term_struct(
    py::module& m,
    py::class_<ComplexTaylorTermProxy>& cls) {
  cls.def(py::init<>())
      // ComplexTaylorTermProxy.coef (0D_NOT_complex -
      .def_property(
          "coef",
          &ComplexTaylorTermProxy::coef,
          &ComplexTaylorTermProxy::set_coef)
      // ComplexTaylorTermProxy.expn (1D_NOT_integer -
      .def_property_readonly("expn", &ComplexTaylorTermProxy::expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ComplexTaylorTermProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ComplexTaylorTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorTermProxy& self) {
            return ComplexTaylorTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ComplexTaylorTermProxy& self, py::dict& memo) {
            return ComplexTaylorTermProxy(self);
          })

      ;

  bind_FTypeArrayND<ComplexTaylorTermProxyArray1D>(
      m, "ComplexTaylorTermStructArray1D");
  bind_FTypeAlloc1D<ComplexTaylorTermProxyAlloc1D>(
      m, "ComplexTaylorTermStructAlloc1D");
  // 2D ComplexTaylorTermProxy arrays are not used in structs/routines
  // 3D ComplexTaylorTermProxy arrays are not used in structs/routines
}

// =============================================================================
// control_ramp1_struct
void init_control_ramp1_struct(
    py::module& m,
    py::class_<ControlRamp1Proxy>& cls) {
  cls.def(py::init<>())
      // ControlRamp1Proxy.y_knot (1D_ALLOC_real -
      .def_property_readonly("y_knot", &ControlRamp1Proxy::y_knot)
      // ControlRamp1Proxy.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlRamp1Proxy::stack)
      // ControlRamp1Proxy.attribute (0D_NOT_character - Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps.
      .def_property(
          "attribute",
          &ControlRamp1Proxy::attribute,
          &ControlRamp1Proxy::set_attribute)
      // ControlRamp1Proxy.slave_name (0D_NOT_character - Name of slave.
      .def_property(
          "slave_name",
          &ControlRamp1Proxy::slave_name,
          &ControlRamp1Proxy::set_slave_name)
      // ControlRamp1Proxy.is_controller (0D_NOT_logical - Is the slave a controller? If so bookkeeping is different.
      .def_property(
          "is_controller",
          &ControlRamp1Proxy::is_controller,
          &ControlRamp1Proxy::set_is_controller)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ControlRamp1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ControlRamp1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlRamp1Proxy& self) {
            return ControlRamp1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlRamp1Proxy& self, py::dict& memo) {
            return ControlRamp1Proxy(self);
          })

      ;

  bind_FTypeArrayND<ControlRamp1ProxyArray1D>(m, "ControlRamp1StructArray1D");
  bind_FTypeAlloc1D<ControlRamp1ProxyAlloc1D>(m, "ControlRamp1StructAlloc1D");
  // 2D ControlRamp1Proxy arrays are not used in structs/routines
  // 3D ControlRamp1Proxy arrays are not used in structs/routines
}

// =============================================================================
// control_struct
void init_control_struct(py::module& m, py::class_<ControlProxy>& cls) {
  cls.def(py::init<>())
      // ControlProxy.value (0D_NOT_real - Used by group, and overlay elements.
      .def_property("value", &ControlProxy::value, &ControlProxy::set_value)
      // ControlProxy.y_knot (1D_ALLOC_real -
      .def_property_readonly("y_knot", &ControlProxy::y_knot)
      // ControlProxy.stack (1D_ALLOC_type - Evaluation stack
      .def_property_readonly("stack", &ControlProxy::stack)
      // ControlProxy.slave (0D_NOT_type -
      .def_property("slave", &ControlProxy::slave, &ControlProxy::set_slave)
      // ControlProxy.lord (0D_NOT_type -
      .def_property("lord", &ControlProxy::lord, &ControlProxy::set_lord)
      // ControlProxy.slave_name (0D_NOT_character - Name of slave.
      .def_property(
          "slave_name",
          &ControlProxy::slave_name,
          &ControlProxy::set_slave_name)
      // ControlProxy.attribute (0D_NOT_character - Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps. Set to 'INPUT' or 'OUTPUT' for feedback slaves.
      .def_property(
          "attribute", &ControlProxy::attribute, &ControlProxy::set_attribute)
      // ControlProxy.ix_attrib (0D_NOT_integer - Index of attribute controlled. See note above!
      .def_property(
          "ix_attrib", &ControlProxy::ix_attrib, &ControlProxy::set_ix_attrib)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ControlProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const ControlProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlProxy& self) {
            return ControlProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlProxy& self, py::dict& memo) {
            return ControlProxy(self);
          })

      ;

  bind_FTypeArrayND<ControlProxyArray1D>(m, "ControlStructArray1D");
  bind_FTypeAlloc1D<ControlProxyAlloc1D>(m, "ControlStructAlloc1D");
  // 2D ControlProxy arrays are not used in structs/routines
  // 3D ControlProxy arrays are not used in structs/routines
}

// =============================================================================
// control_var1_struct
void init_control_var1_struct(
    py::module& m,
    py::class_<ControlVar1Proxy>& cls) {
  cls.def(py::init<>())
      // ControlVar1Proxy.name (0D_NOT_character -
      .def_property(
          "name", &ControlVar1Proxy::name, &ControlVar1Proxy::set_name)
      // ControlVar1Proxy.value (0D_NOT_real -
      .def_property(
          "value", &ControlVar1Proxy::value, &ControlVar1Proxy::set_value)
      // ControlVar1Proxy.old_value (0D_NOT_real -
      .def_property(
          "old_value",
          &ControlVar1Proxy::old_value,
          &ControlVar1Proxy::set_old_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ControlVar1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ControlVar1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlVar1Proxy& self) {
            return ControlVar1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControlVar1Proxy& self, py::dict& memo) {
            return ControlVar1Proxy(self);
          })

      ;

  bind_FTypeArrayND<ControlVar1ProxyArray1D>(m, "ControlVar1StructArray1D");
  bind_FTypeAlloc1D<ControlVar1ProxyAlloc1D>(m, "ControlVar1StructAlloc1D");
  // 2D ControlVar1Proxy arrays are not used in structs/routines
  // 3D ControlVar1Proxy arrays are not used in structs/routines
}

// =============================================================================
// controller_struct
void init_controller_struct(py::module& m, py::class_<ControllerProxy>& cls) {
  cls.def(py::init<>())
      // ControllerProxy.var (1D_ALLOC_type -
      .def_property_readonly("var", &ControllerProxy::var)
      // ControllerProxy.ramp (1D_ALLOC_type - For ramper lord elements
      .def_property_readonly("ramp", &ControllerProxy::ramp)
      // ControllerProxy.ramper_lord (1D_ALLOC_type - Ramper lord info for this slave
      .def_property_readonly("ramper_lord", &ControllerProxy::ramper_lord)
      // ControllerProxy.x_knot (1D_ALLOC_real -
      .def_property_readonly("x_knot", &ControllerProxy::x_knot)

      .def(
          "__repr__",
          [](const ControllerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControllerProxy& self) {
            return ControllerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ControllerProxy& self, py::dict& memo) {
            return ControllerProxy(self);
          })

      ;

  // 1D ControllerProxy arrays are not used in structs/routines
  // 2D ControllerProxy arrays are not used in structs/routines
  // 3D ControllerProxy arrays are not used in structs/routines
}

// =============================================================================
// coord_array_struct
void init_coord_array_struct(py::module& m, py::class_<CoordArrayProxy>& cls) {
  cls.def(py::init<>())
      // CoordArrayProxy.orbit (1D_ALLOC_type -
      .def_property_readonly("orbit", &CoordArrayProxy::orbit)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordArrayProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CoordArrayProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordArrayProxy& self) {
            return CoordArrayProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CoordArrayProxy& self, py::dict& memo) {
            return CoordArrayProxy(self);
          })

      ;

  bind_FTypeArrayND<CoordArrayProxyArray1D>(m, "CoordArrayStructArray1D");
  bind_FTypeAlloc1D<CoordArrayProxyAlloc1D>(m, "CoordArrayStructAlloc1D");
  // 2D CoordArrayProxy arrays are not used in structs/routines
  // 3D CoordArrayProxy arrays are not used in structs/routines
}

// =============================================================================
// coord_struct
void init_coord_struct(py::module& m, py::class_<CoordProxy>& cls) {
  cls.def(py::init<>())
      // CoordProxy.vec (1D_NOT_real - (x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual.
      .def_property_readonly("vec", &CoordProxy::vec)
      // CoordProxy.s (0D_NOT_real - Longitudinal position
      .def_property("s", &CoordProxy::s, &CoordProxy::set_s)
      // CoordProxy.t (0D_NOT_real16 - Absolute time (not relative to reference). Note: Quad precision!
      .def_property("t", &CoordProxy::t, &CoordProxy::set_t)
      // CoordProxy.spin (1D_NOT_real - Spin.
      .def_property_readonly("spin", &CoordProxy::spin)
      // CoordProxy.field (1D_NOT_real - Photon E-field intensity (x,y).
      .def_property_readonly("field", &CoordProxy::field)
      // CoordProxy.phase (1D_NOT_real - Photon E-field phase (x,y). For charged particles, phase(1) is RF phase.
      .def_property_readonly("phase", &CoordProxy::phase)
      // CoordProxy.charge (0D_NOT_real - Macroparticle weight (which is different from particle species charge). For some space charge calcs the weight is in Coulombs.
      .def_property("charge", &CoordProxy::charge, &CoordProxy::set_charge)
      // CoordProxy.dt_ref (0D_NOT_real - Used in: * time tracking for computing z. * by coherent photons = path_length/c_light.
      .def_property("dt_ref", &CoordProxy::dt_ref, &CoordProxy::set_dt_ref)
      // CoordProxy.r (0D_NOT_real - For general use. Not used by Bmad.
      .def_property("r", &CoordProxy::r, &CoordProxy::set_r)
      // CoordProxy.p0c (0D_NOT_real - For non-photons: Reference momentum. For photons: Photon momentum (not reference).
      .def_property("p0c", &CoordProxy::p0c, &CoordProxy::set_p0c)
      // CoordProxy.E_potential (0D_NOT_real - Potential energy.
      .def_property(
          "E_potential", &CoordProxy::E_potential, &CoordProxy::set_E_potential)
      // CoordProxy.beta (0D_NOT_real - Velocity / c_light.
      .def_property("beta", &CoordProxy::beta, &CoordProxy::set_beta)
      // CoordProxy.ix_ele (0D_NOT_integer - Index of the lattice element the particle is in. May be -1 if element is not associated with a lattice.
      .def_property("ix_ele", &CoordProxy::ix_ele, &CoordProxy::set_ix_ele)
      // CoordProxy.ix_branch (0D_NOT_integer - Index of the lattice branch the particle is in.
      .def_property(
          "ix_branch", &CoordProxy::ix_branch, &CoordProxy::set_ix_branch)
      // CoordProxy.ix_turn (0D_NOT_integer - Turn index for multiturn tracking.
      .def_property("ix_turn", &CoordProxy::ix_turn, &CoordProxy::set_ix_turn)
      // CoordProxy.ix_user (0D_NOT_integer - For general use, not used by Bmad.
      .def_property("ix_user", &CoordProxy::ix_user, &CoordProxy::set_ix_user)
      // CoordProxy.state (0D_NOT_integer - alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc.
      .def_property("state", &CoordProxy::state, &CoordProxy::set_state)
      // CoordProxy.direction (0D_NOT_integer - +1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the element orientation.
      .def_property(
          "direction", &CoordProxy::direction, &CoordProxy::set_direction)
      // CoordProxy.time_dir (0D_NOT_integer - +1 or -1. Time direction. -1 => Traveling backwards in time.
      .def_property(
          "time_dir", &CoordProxy::time_dir, &CoordProxy::set_time_dir)
      // CoordProxy.species (0D_NOT_integer - positron$, proton$, etc.
      .def_property("species", &CoordProxy::species, &CoordProxy::set_species)
      // CoordProxy.location (0D_NOT_integer - upstream_end$, inside$, or downstream_end$
      .def_property(
          "location", &CoordProxy::location, &CoordProxy::set_location)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return CoordProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const CoordProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordProxy& self) {
            return CoordProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CoordProxy& self, py::dict& memo) {
            return CoordProxy(self);
          })

      ;

  bind_FTypeArrayND<CoordProxyArray1D>(m, "CoordStructArray1D");
  bind_FTypeAlloc1D<CoordProxyAlloc1D>(m, "CoordStructAlloc1D");
  // 2D CoordProxy arrays are not used in structs/routines
  // 3D CoordProxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_struct
void init_cylindrical_map_struct(
    py::module& m,
    py::class_<CylindricalMapProxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapProxy.m (0D_NOT_integer - Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)
      .def_property("m", &CylindricalMapProxy::m, &CylindricalMapProxy::set_m)
      // CylindricalMapProxy.harmonic (0D_NOT_integer - Harmonic of fundamental
      .def_property(
          "harmonic",
          &CylindricalMapProxy::harmonic,
          &CylindricalMapProxy::set_harmonic)
      // CylindricalMapProxy.phi0_fieldmap (0D_NOT_real - Mode oscillates as: twopi * (f * t + phi0_fieldmap)
      .def_property(
          "phi0_fieldmap",
          &CylindricalMapProxy::phi0_fieldmap,
          &CylindricalMapProxy::set_phi0_fieldmap)
      // CylindricalMapProxy.theta0_azimuth (0D_NOT_real - Azimuthal ((x, y) plane) orientation of mode.
      .def_property(
          "theta0_azimuth",
          &CylindricalMapProxy::theta0_azimuth,
          &CylindricalMapProxy::set_theta0_azimuth)
      // CylindricalMapProxy.field_scale (0D_NOT_real - Factor to scale the fields by
      .def_property(
          "field_scale",
          &CylindricalMapProxy::field_scale,
          &CylindricalMapProxy::set_field_scale)
      // CylindricalMapProxy.master_parameter (0D_NOT_integer - Master parameter in ele%value(:) array to use for scaling the field.
      .def_property(
          "master_parameter",
          &CylindricalMapProxy::master_parameter,
          &CylindricalMapProxy::set_master_parameter)
      // CylindricalMapProxy.ele_anchor_pt (0D_NOT_integer - anchor_beginning$, anchor_center$, or anchor_end$
      .def_property(
          "ele_anchor_pt",
          &CylindricalMapProxy::ele_anchor_pt,
          &CylindricalMapProxy::set_ele_anchor_pt)
      // CylindricalMapProxy.dz (0D_NOT_real - Distance between sampled field points.
      .def_property(
          "dz", &CylindricalMapProxy::dz, &CylindricalMapProxy::set_dz)
      // CylindricalMapProxy.r0 (1D_NOT_real - Field origin offset.
      .def_property_readonly("r0", &CylindricalMapProxy::r0)
      // CylindricalMapProxy.ptr (0D_PTR_type -
      .def_property(
          "ptr", &CylindricalMapProxy::ptr, &CylindricalMapProxy::set_ptr)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CylindricalMapProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CylindricalMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapProxy& self) {
            return CylindricalMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapProxy& self, py::dict& memo) {
            return CylindricalMapProxy(self);
          })

      ;

  bind_FTypeArrayND<CylindricalMapProxyArray1D>(
      m, "CylindricalMapStructArray1D");
  bind_FTypeAlloc1D<CylindricalMapProxyAlloc1D>(
      m, "CylindricalMapStructAlloc1D");
  // 2D CylindricalMapProxy arrays are not used in structs/routines
  // 3D CylindricalMapProxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term1_struct
void init_cylindrical_map_term1_struct(
    py::module& m,
    py::class_<CylindricalMapTerm1Proxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapTerm1Proxy.e_coef (0D_NOT_complex -
      .def_property(
          "e_coef",
          &CylindricalMapTerm1Proxy::e_coef,
          &CylindricalMapTerm1Proxy::set_e_coef)
      // CylindricalMapTerm1Proxy.b_coef (0D_NOT_complex -
      .def_property(
          "b_coef",
          &CylindricalMapTerm1Proxy::b_coef,
          &CylindricalMapTerm1Proxy::set_b_coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return CylindricalMapTerm1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const CylindricalMapTerm1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTerm1Proxy& self) {
            return CylindricalMapTerm1Proxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapTerm1Proxy& self, py::dict& memo) {
            return CylindricalMapTerm1Proxy(self);
          })

      ;

  bind_FTypeArrayND<CylindricalMapTerm1ProxyArray1D>(
      m, "CylindricalMapTerm1StructArray1D");
  bind_FTypeAlloc1D<CylindricalMapTerm1ProxyAlloc1D>(
      m, "CylindricalMapTerm1StructAlloc1D");
  // 2D CylindricalMapTerm1Proxy arrays are not used in structs/routines
  // 3D CylindricalMapTerm1Proxy arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term_struct
void init_cylindrical_map_term_struct(
    py::module& m,
    py::class_<CylindricalMapTermProxy>& cls) {
  cls.def(py::init<>())
      // CylindricalMapTermProxy.file (0D_NOT_character - Input file name. Used also as ID for instances.
      .def_property(
          "file",
          &CylindricalMapTermProxy::file,
          &CylindricalMapTermProxy::set_file)
      // CylindricalMapTermProxy.n_link (0D_NOT_integer - For memory management of this structure
      .def_property(
          "n_link",
          &CylindricalMapTermProxy::n_link,
          &CylindricalMapTermProxy::set_n_link)
      // CylindricalMapTermProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &CylindricalMapTermProxy::term)

      .def(
          "__repr__",
          [](const CylindricalMapTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTermProxy& self) {
            return CylindricalMapTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const CylindricalMapTermProxy& self, py::dict& memo) {
            return CylindricalMapTermProxy(self);
          })

      ;

  // 1D CylindricalMapTermProxy arrays are not used in structs/routines
  // 2D CylindricalMapTermProxy arrays are not used in structs/routines
  // 3D CylindricalMapTermProxy arrays are not used in structs/routines
}