#include "pybmad/generated/structs_e.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// ele_pointer_struct
void init_ele_pointer_struct(py::module& m, py::class_<ElePointerProxy>& cls) {
  cls.def(py::init<>())
      // ElePointerProxy.ele (0D_PTR_type -
      .def_property("ele", &ElePointerProxy::ele, &ElePointerProxy::set_ele)
      // ElePointerProxy.loc (0D_NOT_type -
      .def_property("loc", &ElePointerProxy::loc, &ElePointerProxy::set_loc)
      // ElePointerProxy.id (0D_NOT_integer - For general use. Not used by Bmad. In particular, used by Tao to designate universe ele is in.
      .def_property("id", &ElePointerProxy::id, &ElePointerProxy::set_id)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ElePointerProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ElePointerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ElePointerProxy& self) {
            return ElePointerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ElePointerProxy& self, py::dict& memo) {
            return ElePointerProxy(self);
          })

      ;

  bind_FTypeArrayND<ElePointerProxyArray1D>(m, "ElePointerStructArray1D");
  bind_FTypeAlloc1D<ElePointerProxyAlloc1D>(m, "ElePointerStructAlloc1D");
  // 2D ElePointerProxy arrays are not used in structs/routines
  // 3D ElePointerProxy arrays are not used in structs/routines
}

// =============================================================================
// ele_struct
void init_ele_struct(py::module& m, py::class_<EleProxy>& cls) {
  cls.def(py::init<>())
      // EleProxy.name (0D_NOT_character - name of element.
      .def_property("name", &EleProxy::name, &EleProxy::set_name)
      // EleProxy.type (0D_NOT_character - type name.
      .def_property("type", &EleProxy::type, &EleProxy::set_type)
      // EleProxy.alias (0D_NOT_character - Another name.
      .def_property("alias", &EleProxy::alias, &EleProxy::set_alias)
      // EleProxy.component_name (0D_NOT_character - Used by overlays, multipass patch, etc.
      .def_property(
          "component_name",
          &EleProxy::component_name,
          &EleProxy::set_component_name)
      // EleProxy.descrip (0D_PTR_character - Description string.
      .def_property("descrip", &EleProxy::descrip, &EleProxy::set_descrip)
      // EleProxy.a (0D_NOT_type - Twiss parameters at end of element
      .def_property("a", &EleProxy::a, &EleProxy::set_a)
      // EleProxy.b (0D_NOT_type - Twiss parameters at end of element
      .def_property("b", &EleProxy::b, &EleProxy::set_b)
      // EleProxy.z (0D_NOT_type - Twiss parameters at end of element
      .def_property("z", &EleProxy::z, &EleProxy::set_z)
      // EleProxy.x (0D_NOT_type - Projected dispersions.
      .def_property("x", &EleProxy::x, &EleProxy::set_x)
      // EleProxy.y (0D_NOT_type - Projected dispersions.
      .def_property("y", &EleProxy::y, &EleProxy::set_y)
      // EleProxy.ac_kick (0D_PTR_type - ac_kicker element parameters.
      .def_property("ac_kick", &EleProxy::ac_kick, &EleProxy::set_ac_kick)
      // EleProxy.bookkeeping_state (0D_NOT_type - Attribute bookkeeping
      .def_property(
          "bookkeeping_state",
          &EleProxy::bookkeeping_state,
          &EleProxy::set_bookkeeping_state)
      // EleProxy.branch (0D_PTR_type - Pointer to branch containing element.
      .def_property("branch", &EleProxy::branch, &EleProxy::set_branch)
      // EleProxy.control (0D_PTR_type - group & overlay variables.
      .def_property("control", &EleProxy::control, &EleProxy::set_control)
      // EleProxy.rf (0D_PTR_type - RF parameters.
      .def_property("rf", &EleProxy::rf, &EleProxy::set_rf)
      // EleProxy.lord (0D_PTR_type - Pointer to a slice lord.
      .def_property("lord", &EleProxy::lord, &EleProxy::set_lord)
      // EleProxy.ptc_fibre (0D_PTR_type - PTC track corresponding to this ele.
      .def_property("ptc_fibre", &EleProxy::ptc_fibre, &EleProxy::set_ptc_fibre)
      // EleProxy.floor (0D_NOT_type -
      .def_property("floor", &EleProxy::floor, &EleProxy::set_floor)
      // EleProxy.high_energy_space_charge (0D_PTR_type -
      .def_property(
          "high_energy_space_charge",
          &EleProxy::high_energy_space_charge,
          &EleProxy::set_high_energy_space_charge)
      // EleProxy.mode3 (0D_PTR_type - 6D normal mode structure.
      .def_property("mode3", &EleProxy::mode3, &EleProxy::set_mode3)
      // EleProxy.photon (0D_PTR_type -
      .def_property("photon", &EleProxy::photon, &EleProxy::set_photon)
      // EleProxy.rad_map (0D_PTR_type - Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map.
      .def_property("rad_map", &EleProxy::rad_map, &EleProxy::set_rad_map)
      // EleProxy.taylor (1D_NOT_type - Phase space Taylor map.
      .def_property_readonly("taylor", &EleProxy::taylor)
      // EleProxy.spin_taylor_ref_orb_in (1D_NOT_real -
      .def_property_readonly(
          "spin_taylor_ref_orb_in", &EleProxy::spin_taylor_ref_orb_in)
      // EleProxy.spin_taylor (1D_NOT_type - Quaternion Spin Taylor map.
      .def_property_readonly("spin_taylor", &EleProxy::spin_taylor)
      // EleProxy.wake (0D_PTR_type - Wakes
      .def_property("wake", &EleProxy::wake, &EleProxy::set_wake)
      // EleProxy.wall3d (1D_PTR_type - Chamber or capillary wall E/M field structs.
      .def_property_readonly("wall3d", &EleProxy::wall3d)
      // EleProxy.cartesian_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cartesian_map", &EleProxy::cartesian_map)
      // EleProxy.cylindrical_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cylindrical_map", &EleProxy::cylindrical_map)
      // EleProxy.gen_grad_map (1D_PTR_type - Used to define E/M fields.
      .def_property_readonly("gen_grad_map", &EleProxy::gen_grad_map)
      // EleProxy.grid_field (1D_PTR_type - Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference particle which is generally the zero orbit (non-zero, for example, in the second slice of a sliced wiggler).
      .def_property_readonly("grid_field", &EleProxy::grid_field)
      // EleProxy.map_ref_orb_in (0D_NOT_type - Entrance end transfer map ref orbit
      .def_property(
          "map_ref_orb_in",
          &EleProxy::map_ref_orb_in,
          &EleProxy::set_map_ref_orb_in)
      // EleProxy.map_ref_orb_out (0D_NOT_type - Exit end transfer map ref orbit
      .def_property(
          "map_ref_orb_out",
          &EleProxy::map_ref_orb_out,
          &EleProxy::set_map_ref_orb_out)
      // EleProxy.time_ref_orb_in (0D_NOT_type - Reference orbit at entrance end for ref_time calc.
      .def_property(
          "time_ref_orb_in",
          &EleProxy::time_ref_orb_in,
          &EleProxy::set_time_ref_orb_in)
      // EleProxy.time_ref_orb_out (0D_NOT_type - Reference orbit at exit end for ref_time calc.
      .def_property(
          "time_ref_orb_out",
          &EleProxy::time_ref_orb_out,
          &EleProxy::set_time_ref_orb_out)
      // EleProxy.value (1D_NOT_real - attribute values.
      .def_property_readonly("value", &EleProxy::value)
      // EleProxy.old_value (1D_NOT_real - Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit matrices is %map_ref_orb_in/out
      .def_property_readonly("old_value", &EleProxy::old_value)
      // EleProxy.spin_q (2D_NOT_real - 0th and 1st order Spin transport quaternion.
      .def_property_readonly("spin_q", &EleProxy::spin_q)
      // EleProxy.vec0 (1D_NOT_real - 0th order transport vector.
      .def_property_readonly("vec0", &EleProxy::vec0)
      // EleProxy.mat6 (2D_NOT_real - 1st order transport matrix.
      .def_property_readonly("mat6", &EleProxy::mat6)
      // EleProxy.c_mat (2D_NOT_real - 2x2 C coupling matrix
      .def_property_readonly("c_mat", &EleProxy::c_mat)
      // EleProxy.dc_mat_dpz (2D_NOT_real - d(c_mat)/dpz variation.
      .def_property_readonly("dc_mat_dpz", &EleProxy::dc_mat_dpz)
      // EleProxy.gamma_c (0D_NOT_real - gamma associated with C matrix
      .def_property("gamma_c", &EleProxy::gamma_c, &EleProxy::set_gamma_c)
      // EleProxy.s_start (0D_NOT_real - longitudinal ref position at entrance_end
      .def_property("s_start", &EleProxy::s_start, &EleProxy::set_s_start)
      // EleProxy.s (0D_NOT_real - longitudinal ref position at the exit end.
      .def_property("s", &EleProxy::s, &EleProxy::set_s)
      // EleProxy.ref_time (0D_NOT_real - Time ref particle passes exit end.
      .def_property("ref_time", &EleProxy::ref_time, &EleProxy::set_ref_time)
      // EleProxy.a_pole (1D_PTR_real - knl for multipole elements.
      .def_property_readonly("a_pole", &EleProxy::a_pole)
      // EleProxy.b_pole (1D_PTR_real - tilt for multipole elements.
      .def_property_readonly("b_pole", &EleProxy::b_pole)
      // EleProxy.a_pole_elec (1D_PTR_real - Electrostatic multipoles. ksnl for multipole elements.
      .def_property_readonly("a_pole_elec", &EleProxy::a_pole_elec)
      // EleProxy.b_pole_elec (1D_PTR_real - Electrostatic multipoles.
      .def_property_readonly("b_pole_elec", &EleProxy::b_pole_elec)
      // EleProxy.custom (1D_PTR_real - Custom attributes.
      .def_property_readonly("custom", &EleProxy::custom)
      // EleProxy.r (3D_PTR_real - For general use. Not used by Bmad.
      .def_property_readonly("r", &EleProxy::r)
      // EleProxy.key (0D_NOT_integer - Element class (quadrupole, etc.).
      .def_property("key", &EleProxy::key, &EleProxy::set_key)
      // EleProxy.sub_key (0D_NOT_integer - Records bend input type.
      .def_property("sub_key", &EleProxy::sub_key, &EleProxy::set_sub_key)
      // EleProxy.ix_ele (0D_NOT_integer - Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements.
      .def_property("ix_ele", &EleProxy::ix_ele, &EleProxy::set_ix_ele)
      // EleProxy.ix_branch (0D_NOT_integer - Index in lat%branch(:) array. Note: lat%ele => lat%branch(0).
      .def_property("ix_branch", &EleProxy::ix_branch, &EleProxy::set_ix_branch)
      // EleProxy.lord_status (0D_NOT_integer - Type of lord element this is. overlay_lord$, etc.
      .def_property(
          "lord_status", &EleProxy::lord_status, &EleProxy::set_lord_status)
      // EleProxy.n_slave (0D_NOT_integer - Number of slaves (except field overlap slaves) of this element.
      .def_property("n_slave", &EleProxy::n_slave, &EleProxy::set_n_slave)
      // EleProxy.n_slave_field (0D_NOT_integer - Number of field slaves of this element.
      .def_property(
          "n_slave_field",
          &EleProxy::n_slave_field,
          &EleProxy::set_n_slave_field)
      // EleProxy.ix1_slave (0D_NOT_integer - Pointer index to this element's slaves.
      .def_property("ix1_slave", &EleProxy::ix1_slave, &EleProxy::set_ix1_slave)
      // EleProxy.slave_status (0D_NOT_integer - Type of slave element this is. multipass_slave$, slice_slave$, etc.
      .def_property(
          "slave_status", &EleProxy::slave_status, &EleProxy::set_slave_status)
      // EleProxy.n_lord (0D_NOT_integer - Number of lords (except field overlap and ramper lords).
      .def_property("n_lord", &EleProxy::n_lord, &EleProxy::set_n_lord)
      // EleProxy.n_lord_field (0D_NOT_integer - Number of field lords of this element.
      .def_property(
          "n_lord_field", &EleProxy::n_lord_field, &EleProxy::set_n_lord_field)
      // EleProxy.n_lord_ramper (0D_NOT_integer - Number of ramper lords.
      .def_property(
          "n_lord_ramper",
          &EleProxy::n_lord_ramper,
          &EleProxy::set_n_lord_ramper)
      // EleProxy.ic1_lord (0D_NOT_integer - Pointer index to this element's lords.
      .def_property("ic1_lord", &EleProxy::ic1_lord, &EleProxy::set_ic1_lord)
      // EleProxy.ix_pointer (0D_NOT_integer - For general use. Not used by Bmad.
      .def_property(
          "ix_pointer", &EleProxy::ix_pointer, &EleProxy::set_ix_pointer)
      // EleProxy.ixx (0D_NOT_integer - Index for Bmad internal use.
      .def_property("ixx", &EleProxy::ixx, &EleProxy::set_ixx)
      // EleProxy.iyy (0D_NOT_integer - Index for Bmad internal use.
      .def_property("iyy", &EleProxy::iyy, &EleProxy::set_iyy)
      // EleProxy.izz (0D_NOT_integer - Index for Bmad internal use.
      .def_property("izz", &EleProxy::izz, &EleProxy::set_izz)
      // EleProxy.mat6_calc_method (0D_NOT_integer - taylor$, symp_lie_ptc$, etc.
      .def_property(
          "mat6_calc_method",
          &EleProxy::mat6_calc_method,
          &EleProxy::set_mat6_calc_method)
      // EleProxy.tracking_method (0D_NOT_integer - taylor$, linear$, etc.
      .def_property(
          "tracking_method",
          &EleProxy::tracking_method,
          &EleProxy::set_tracking_method)
      // EleProxy.spin_tracking_method (0D_NOT_integer - symp_lie_ptc$, etc.
      .def_property(
          "spin_tracking_method",
          &EleProxy::spin_tracking_method,
          &EleProxy::set_spin_tracking_method)
      // EleProxy.csr_method (0D_NOT_integer - or one_dim$ ('1_dim'), steady_state_3d$
      .def_property(
          "csr_method", &EleProxy::csr_method, &EleProxy::set_csr_method)
      // EleProxy.space_charge_method (0D_NOT_integer - slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$
      .def_property(
          "space_charge_method",
          &EleProxy::space_charge_method,
          &EleProxy::set_space_charge_method)
      // EleProxy.ptc_integration_type (0D_NOT_integer - drift_kick$, matrix_kick$, or ripken_kick$
      .def_property(
          "ptc_integration_type",
          &EleProxy::ptc_integration_type,
          &EleProxy::set_ptc_integration_type)
      // EleProxy.field_calc (0D_NOT_integer - no_field$, fieldmap$, refer_to_lords$, or custom$
      .def_property(
          "field_calc", &EleProxy::field_calc, &EleProxy::set_field_calc)
      // EleProxy.aperture_at (0D_NOT_integer - Aperture location: entrance_end$, ...
      .def_property(
          "aperture_at", &EleProxy::aperture_at, &EleProxy::set_aperture_at)
      // EleProxy.aperture_type (0D_NOT_integer - rectangular$, elliptical$, auto_aperture$, ...
      .def_property(
          "aperture_type",
          &EleProxy::aperture_type,
          &EleProxy::set_aperture_type)
      // EleProxy.ref_species (0D_NOT_integer - Reference species
      .def_property(
          "ref_species", &EleProxy::ref_species, &EleProxy::set_ref_species)
      // EleProxy.orientation (0D_NOT_integer - -1 -> Element is longitudinally reversed. +1 -> Normal.
      .def_property(
          "orientation", &EleProxy::orientation, &EleProxy::set_orientation)
      // EleProxy.symplectify (0D_NOT_logical - Symplectify mat6 matrices.
      .def_property(
          "symplectify", &EleProxy::symplectify, &EleProxy::set_symplectify)
      // EleProxy.mode_flip (0D_NOT_logical - Have the normal modes traded places?
      .def_property("mode_flip", &EleProxy::mode_flip, &EleProxy::set_mode_flip)
      // EleProxy.multipoles_on (0D_NOT_logical - For turning multipoles on/off
      .def_property(
          "multipoles_on",
          &EleProxy::multipoles_on,
          &EleProxy::set_multipoles_on)
      // EleProxy.scale_multipoles (0D_NOT_logical - Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the element?
      .def_property(
          "scale_multipoles",
          &EleProxy::scale_multipoles,
          &EleProxy::set_scale_multipoles)
      // EleProxy.taylor_map_includes_offsets (0D_NOT_logical - Taylor map calculated with element misalignments?
      .def_property(
          "taylor_map_includes_offsets",
          &EleProxy::taylor_map_includes_offsets,
          &EleProxy::set_taylor_map_includes_offsets)
      // EleProxy.field_master (0D_NOT_logical - Calculate strength from the field value?
      .def_property(
          "field_master", &EleProxy::field_master, &EleProxy::set_field_master)
      // EleProxy.is_on (0D_NOT_logical - For turning element on/off.
      .def_property("is_on", &EleProxy::is_on, &EleProxy::set_is_on)
      // EleProxy.logic (0D_NOT_logical - For general use. Not used by Bmad (except during lattice parsing).
      .def_property("logic", &EleProxy::logic, &EleProxy::set_logic)
      // EleProxy.bmad_logic (0D_NOT_logical - For Bmad internal use only.
      .def_property(
          "bmad_logic", &EleProxy::bmad_logic, &EleProxy::set_bmad_logic)
      // EleProxy.select (0D_NOT_logical - For Bmad internal use only.
      .def_property("select", &EleProxy::select, &EleProxy::set_select)
      // EleProxy.offset_moves_aperture (0D_NOT_logical - element offsets affects aperture? ! final :: ele_finalizer
      .def_property(
          "offset_moves_aperture",
          &EleProxy::offset_moves_aperture,
          &EleProxy::set_offset_moves_aperture)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EleProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const EleProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EleProxy& self) {
            return EleProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EleProxy& self, py::dict& memo) { return EleProxy(self); })

      ;

  bind_FTypeArrayND<EleProxyArray1D>(m, "EleStructArray1D");
  bind_FTypeAlloc1D<EleProxyAlloc1D>(m, "EleStructAlloc1D");
  // 2D EleProxy arrays are not used in structs/routines
  // 3D EleProxy arrays are not used in structs/routines
}

// =============================================================================
// ellipse_beam_init_struct
void init_ellipse_beam_init_struct(
    py::module& m,
    py::class_<EllipseBeamInitProxy>& cls) {
  cls.def(py::init<>())
      // EllipseBeamInitProxy.part_per_ellipse (0D_NOT_integer - number of particles per ellipse
      .def_property(
          "part_per_ellipse",
          &EllipseBeamInitProxy::part_per_ellipse,
          &EllipseBeamInitProxy::set_part_per_ellipse)
      // EllipseBeamInitProxy.n_ellipse (0D_NOT_integer - number of ellipses (>= 1)
      .def_property(
          "n_ellipse",
          &EllipseBeamInitProxy::n_ellipse,
          &EllipseBeamInitProxy::set_n_ellipse)
      // EllipseBeamInitProxy.sigma_cutoff (0D_NOT_real - sigma cutoff of the representation
      .def_property(
          "sigma_cutoff",
          &EllipseBeamInitProxy::sigma_cutoff,
          &EllipseBeamInitProxy::set_sigma_cutoff)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return EllipseBeamInitProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const EllipseBeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EllipseBeamInitProxy& self) {
            return EllipseBeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EllipseBeamInitProxy& self, py::dict& memo) {
            return EllipseBeamInitProxy(self);
          })

      ;

  bind_FTypeArrayND<EllipseBeamInitProxyArray1D>(
      m, "EllipseBeamInitStructArray1D");
  bind_FTypeAlloc1D<EllipseBeamInitProxyAlloc1D>(
      m, "EllipseBeamInitStructAlloc1D");
  // 2D EllipseBeamInitProxy arrays are not used in structs/routines
  // 3D EllipseBeamInitProxy arrays are not used in structs/routines
}

// =============================================================================
// em_field_struct
void init_em_field_struct(py::module& m, py::class_<EmFieldProxy>& cls) {
  cls.def(py::init<>())
      // EmFieldProxy.E (1D_NOT_real - electric field.
      .def_property_readonly("E", &EmFieldProxy::E)
      // EmFieldProxy.B (1D_NOT_real - magnetic field.
      .def_property_readonly("B", &EmFieldProxy::B)
      // EmFieldProxy.dE (2D_NOT_real - electric field gradient.
      .def_property_readonly("dE", &EmFieldProxy::dE)
      // EmFieldProxy.dB (2D_NOT_real - magnetic field gradient.
      .def_property_readonly("dB", &EmFieldProxy::dB)
      // EmFieldProxy.phi (0D_NOT_real - Electric scalar potential.
      .def_property("phi", &EmFieldProxy::phi, &EmFieldProxy::set_phi)
      // EmFieldProxy.phi_B (0D_NOT_real - Magnetic scalar potential.
      .def_property("phi_B", &EmFieldProxy::phi_B, &EmFieldProxy::set_phi_B)
      // EmFieldProxy.A (1D_NOT_real - Magnetic vector potential.
      .def_property_readonly("A", &EmFieldProxy::A)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmFieldProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const EmFieldProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmFieldProxy& self) {
            return EmFieldProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmFieldProxy& self, py::dict& memo) {
            return EmFieldProxy(self);
          })

      ;

  bind_FTypeArrayND<EmFieldProxyArray1D>(m, "EmFieldStructArray1D");
  bind_FTypeAlloc1D<EmFieldProxyAlloc1D>(m, "EmFieldStructAlloc1D");
  // 2D EmFieldProxy arrays are not used in structs/routines
  // 3D EmFieldProxy arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_struct
void init_em_taylor_struct(py::module& m, py::class_<EmTaylorProxy>& cls) {
  cls.def(py::init<>())
      // EmTaylorProxy.ref (0D_NOT_real -
      .def_property("ref", &EmTaylorProxy::ref, &EmTaylorProxy::set_ref)
      // EmTaylorProxy.term (1D_ALLOC_type -
      .def_property_readonly("term", &EmTaylorProxy::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmTaylorProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__", [](const EmTaylorProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorProxy& self) {
            return EmTaylorProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmTaylorProxy& self, py::dict& memo) {
            return EmTaylorProxy(self);
          })

      ;

  bind_FTypeArrayND<EmTaylorProxyArray1D>(m, "EmTaylorStructArray1D");
  bind_FTypeAlloc1D<EmTaylorProxyAlloc1D>(m, "EmTaylorStructAlloc1D");
  // 2D EmTaylorProxy arrays are not used in structs/routines
  // 3D EmTaylorProxy arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_term_struct
void init_em_taylor_term_struct(
    py::module& m,
    py::class_<EmTaylorTermProxy>& cls) {
  cls.def(py::init<>())
      // EmTaylorTermProxy.coef (0D_NOT_real -
      .def_property(
          "coef", &EmTaylorTermProxy::coef, &EmTaylorTermProxy::set_coef)
      // EmTaylorTermProxy.expn (1D_NOT_integer -
      .def_property_readonly("expn", &EmTaylorTermProxy::expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return EmTaylorTermProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const EmTaylorTermProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorTermProxy& self) {
            return EmTaylorTermProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const EmTaylorTermProxy& self, py::dict& memo) {
            return EmTaylorTermProxy(self);
          })

      ;

  bind_FTypeArrayND<EmTaylorTermProxyArray1D>(m, "EmTaylorTermStructArray1D");
  bind_FTypeAlloc1D<EmTaylorTermProxyAlloc1D>(m, "EmTaylorTermStructAlloc1D");
  // 2D EmTaylorTermProxy arrays are not used in structs/routines
  // 3D EmTaylorTermProxy arrays are not used in structs/routines
}

// =============================================================================
// expression_atom_struct
void init_expression_atom_struct(
    py::module& m,
    py::class_<ExpressionAtomProxy>& cls) {
  cls.def(py::init<>())
      // ExpressionAtomProxy.name (0D_NOT_character -
      .def_property(
          "name", &ExpressionAtomProxy::name, &ExpressionAtomProxy::set_name)
      // ExpressionAtomProxy.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name
      .def_property(
          "type", &ExpressionAtomProxy::type, &ExpressionAtomProxy::set_type)
      // ExpressionAtomProxy.value (0D_NOT_real -
      .def_property(
          "value", &ExpressionAtomProxy::value, &ExpressionAtomProxy::set_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ExpressionAtomProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ExpressionAtomProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionAtomProxy& self) {
            return ExpressionAtomProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ExpressionAtomProxy& self, py::dict& memo) {
            return ExpressionAtomProxy(self);
          })

      ;

  bind_FTypeArrayND<ExpressionAtomProxyArray1D>(
      m, "ExpressionAtomStructArray1D");
  bind_FTypeAlloc1D<ExpressionAtomProxyAlloc1D>(
      m, "ExpressionAtomStructAlloc1D");
  // 2D ExpressionAtomProxy arrays are not used in structs/routines
  // 3D ExpressionAtomProxy arrays are not used in structs/routines
}

// =============================================================================
// expression_tree_struct
void init_expression_tree_struct(
    py::module& m,
    py::class_<ExpressionTreeProxy>& cls) {
  cls.def(py::init<>())
      // ExpressionTreeProxy.name (0D_NOT_character -
      .def_property(
          "name", &ExpressionTreeProxy::name, &ExpressionTreeProxy::set_name)
      // ExpressionTreeProxy.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc.
      .def_property(
          "type", &ExpressionTreeProxy::type, &ExpressionTreeProxy::set_type)
      // ExpressionTreeProxy.value (0D_NOT_real -
      .def_property(
          "value", &ExpressionTreeProxy::value, &ExpressionTreeProxy::set_value)
      // ExpressionTreeProxy.node (1D_PTR_type - Child nodes. Note: Pointer used here since Ifort does not support allocatable.
      .def_property_readonly("node", &ExpressionTreeProxy::node)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return ExpressionTreeProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const ExpressionTreeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionTreeProxy& self) {
            return ExpressionTreeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ExpressionTreeProxy& self, py::dict& memo) {
            return ExpressionTreeProxy(self);
          })

      ;

  bind_FTypeArrayND<ExpressionTreeProxyArray1D>(
      m, "ExpressionTreeStructArray1D");
  bind_FTypeAlloc1D<ExpressionTreeProxyAlloc1D>(
      m, "ExpressionTreeStructAlloc1D");
  // 2D ExpressionTreeProxy arrays are not used in structs/routines
  // 3D ExpressionTreeProxy arrays are not used in structs/routines
}