#include "pybmad/generated/structs_e.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// ele_pointer_struct
void init_ele_pointer_struct(py::module &m, py::class_<ElePointerStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const EleStruct>,
             optional_ref<const LatEleLocStruct>,
             std::optional<int>>(),
         py::arg("ele") = py::none(),
         py::arg("loc") = py::none(),
         py::arg("id") = py::none()
  )
      // ElePointerStruct.ele (0D_PTR_type -
      .def_property("ele", &ElePointerStruct::ele, &ElePointerStruct::set_ele)
      // ElePointerStruct.loc (0D_NOT_type -
      .def_property("loc", &ElePointerStruct::loc, &ElePointerStruct::set_loc)
      // ElePointerStruct.id (0D_NOT_integer - For general use. Not used by Bmad. In particular,
      // used by Tao to designate universe ele is in.
      .def_property("id", &ElePointerStruct::id, &ElePointerStruct::set_id)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ElePointerStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ElePointerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ElePointerStruct &self) {
            return ElePointerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ElePointerStruct &self, py::dict &memo) { return ElePointerStruct(self); }
      )

      ;

  bind_FTypeArrayND<ElePointerStructArray1D>(m, "ElePointerStructArray1D");
  bind_FTypeAlloc1D<ElePointerStructAlloc1D>(m, "ElePointerStructAlloc1D");
  // 2D ElePointerStruct arrays are not used in structs/routines
  // 3D ElePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// ele_struct
void init_ele_struct(py::module &m, py::class_<EleStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             optional_ref<const std::string>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const TwissStruct>,
             optional_ref<const XyDispStruct>,
             optional_ref<const XyDispStruct>,
             optional_ref<const AcKickerStruct>,
             optional_ref<const BookkeepingStateStruct>,
             optional_ref<const BranchStruct>,
             optional_ref<const ControllerStruct>,
             optional_ref<const RfEleStruct>,
             optional_ref<const EleStruct>,
             optional_ref<const Fibre>,
             optional_ref<const FloorPositionStruct>,
             optional_ref<const HighEnergySpaceChargeStruct>,
             optional_ref<const Mode3Struct>,
             optional_ref<const PhotonElementStruct>,
             optional_ref<const RadMapEleStruct>,
             optional_ref<const std::vector<double>>,
             optional_ref<const WakeStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<double>>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<double>>>,
             optional_ref<const std::vector<std::vector<double>>>,
             optional_ref<const std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<std::vector<double>>>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("alias") = py::none(),
         py::arg("component_name") = py::none(),
         py::arg("descrip") = py::none(),
         py::arg("a") = py::none(),
         py::arg("b") = py::none(),
         py::arg("z") = py::none(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("ac_kick") = py::none(),
         py::arg("bookkeeping_state") = py::none(),
         py::arg("branch") = py::none(),
         py::arg("control") = py::none(),
         py::arg("rf") = py::none(),
         py::arg("lord") = py::none(),
         py::arg("ptc_fibre") = py::none(),
         py::arg("floor") = py::none(),
         py::arg("high_energy_space_charge") = py::none(),
         py::arg("mode3") = py::none(),
         py::arg("photon") = py::none(),
         py::arg("rad_map") = py::none(),
         py::arg("spin_taylor_ref_orb_in") = py::none(),
         py::arg("wake") = py::none(),
         py::arg("map_ref_orb_in") = py::none(),
         py::arg("map_ref_orb_out") = py::none(),
         py::arg("time_ref_orb_in") = py::none(),
         py::arg("time_ref_orb_out") = py::none(),
         py::arg("value") = py::none(),
         py::arg("old_value") = py::none(),
         py::arg("spin_q") = py::none(),
         py::arg("vec0") = py::none(),
         py::arg("mat6") = py::none(),
         py::arg("c_mat") = py::none(),
         py::arg("dc_mat_dpz") = py::none(),
         py::arg("gamma_c") = py::none(),
         py::arg("s_start") = py::none(),
         py::arg("s") = py::none(),
         py::arg("ref_time") = py::none(),
         py::arg("a_pole") = py::none(),
         py::arg("b_pole") = py::none(),
         py::arg("a_pole_elec") = py::none(),
         py::arg("b_pole_elec") = py::none(),
         py::arg("custom") = py::none(),
         py::arg("r") = py::none(),
         py::arg("key") = py::none(),
         py::arg("sub_key") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("lord_status") = py::none(),
         py::arg("n_slave") = py::none(),
         py::arg("n_slave_field") = py::none(),
         py::arg("ix1_slave") = py::none(),
         py::arg("slave_status") = py::none(),
         py::arg("n_lord") = py::none(),
         py::arg("n_lord_field") = py::none(),
         py::arg("n_lord_ramper") = py::none(),
         py::arg("ic1_lord") = py::none(),
         py::arg("ix_pointer") = py::none(),
         py::arg("ixx") = py::none(),
         py::arg("iyy") = py::none(),
         py::arg("izz") = py::none(),
         py::arg("mat6_calc_method") = py::none(),
         py::arg("tracking_method") = py::none(),
         py::arg("spin_tracking_method") = py::none(),
         py::arg("csr_method") = py::none(),
         py::arg("space_charge_method") = py::none(),
         py::arg("ptc_integration_type") = py::none(),
         py::arg("field_calc") = py::none(),
         py::arg("aperture_at") = py::none(),
         py::arg("aperture_type") = py::none(),
         py::arg("ref_species") = py::none(),
         py::arg("orientation") = py::none(),
         py::arg("symplectify") = py::none(),
         py::arg("mode_flip") = py::none(),
         py::arg("multipoles_on") = py::none(),
         py::arg("scale_multipoles") = py::none(),
         py::arg("taylor_map_includes_offsets") = py::none(),
         py::arg("field_master") = py::none(),
         py::arg("is_on") = py::none(),
         py::arg("logic") = py::none(),
         py::arg("bmad_logic") = py::none(),
         py::arg("select") = py::none(),
         py::arg("offset_moves_aperture") = py::none()
  )
      // EleStruct.name (0D_NOT_character - name of element.
      .def_property("name", &EleStruct::name, &EleStruct::set_name)
      // EleStruct.type (0D_NOT_character - type name.
      .def_property("type", &EleStruct::type, &EleStruct::set_type)
      // EleStruct.alias (0D_NOT_character - Another name.
      .def_property("alias", &EleStruct::alias, &EleStruct::set_alias)
      // EleStruct.component_name (0D_NOT_character - Used by overlays, multipass patch, etc.
      .def_property("component_name", &EleStruct::component_name, &EleStruct::set_component_name)
      // EleStruct.descrip (0D_PTR_character - Description string.
      .def_property("descrip", &EleStruct::descrip, &EleStruct::set_descrip)
      // EleStruct.a (0D_NOT_type - Twiss parameters at end of element
      .def_property("a", &EleStruct::a, &EleStruct::set_a)
      // EleStruct.b (0D_NOT_type - Twiss parameters at end of element
      .def_property("b", &EleStruct::b, &EleStruct::set_b)
      // EleStruct.z (0D_NOT_type - Twiss parameters at end of element
      .def_property("z", &EleStruct::z, &EleStruct::set_z)
      // EleStruct.x (0D_NOT_type - Projected dispersions.
      .def_property("x", &EleStruct::x, &EleStruct::set_x)
      // EleStruct.y (0D_NOT_type - Projected dispersions.
      .def_property("y", &EleStruct::y, &EleStruct::set_y)
      // EleStruct.ac_kick (0D_PTR_type - ac_kicker element parameters.
      .def_property("ac_kick", &EleStruct::ac_kick, &EleStruct::set_ac_kick)
      // EleStruct.bookkeeping_state (0D_NOT_type - Attribute bookkeeping
      .def_property(
          "bookkeeping_state",
          &EleStruct::bookkeeping_state,
          &EleStruct::set_bookkeeping_state
      )
      // EleStruct.branch (0D_PTR_type - Pointer to branch containing element.
      .def_property("branch", &EleStruct::branch, &EleStruct::set_branch)
      // EleStruct.control (0D_PTR_type - group & overlay variables.
      .def_property("control", &EleStruct::control, &EleStruct::set_control)
      // EleStruct.rf (0D_PTR_type - RF parameters.
      .def_property("rf", &EleStruct::rf, &EleStruct::set_rf)
      // EleStruct.lord (0D_PTR_type - Pointer to a slice lord.
      .def_property("lord", &EleStruct::lord, &EleStruct::set_lord)
      // EleStruct.ptc_fibre (0D_PTR_type - PTC track corresponding to this ele.
      .def_property("ptc_fibre", &EleStruct::ptc_fibre, &EleStruct::set_ptc_fibre)
      // EleStruct.floor (0D_NOT_type -
      .def_property("floor", &EleStruct::floor, &EleStruct::set_floor)
      // EleStruct.high_energy_space_charge (0D_PTR_type -
      .def_property(
          "high_energy_space_charge",
          &EleStruct::high_energy_space_charge,
          &EleStruct::set_high_energy_space_charge
      )
      // EleStruct.mode3 (0D_PTR_type - 6D normal mode structure.
      .def_property("mode3", &EleStruct::mode3, &EleStruct::set_mode3)
      // EleStruct.photon (0D_PTR_type -
      .def_property("photon", &EleStruct::photon, &EleStruct::set_photon)
      // EleStruct.rad_map (0D_PTR_type - Radiation kick parameters Note: The reference orbits for
      // spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor
      // maps can be with respect to the zero orbit independent of the orbital map.
      .def_property("rad_map", &EleStruct::rad_map, &EleStruct::set_rad_map)
      // EleStruct.taylor (1D_NOT_type - Phase space Taylor map.
      .def_property_readonly("taylor", &EleStruct::taylor)
      // EleStruct.spin_taylor_ref_orb_in (1D_NOT_real -
      .def_property(
          "spin_taylor_ref_orb_in",
          &EleStruct::spin_taylor_ref_orb_in,
          &EleStruct::set_spin_taylor_ref_orb_in
      )
      // EleStruct.spin_taylor (1D_NOT_type - Quaternion Spin Taylor map.
      .def_property_readonly("spin_taylor", &EleStruct::spin_taylor)
      // EleStruct.wake (0D_PTR_type - Wakes
      .def_property("wake", &EleStruct::wake, &EleStruct::set_wake)
      // EleStruct.wall3d (1D_PTR_type - Chamber or capillary wall E/M field structs.
      .def_property_readonly("wall3d", &EleStruct::wall3d)
      // EleStruct.cartesian_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cartesian_map", &EleStruct::cartesian_map)
      // EleStruct.cylindrical_map (1D_PTR_type - Used to define E/M fields
      .def_property_readonly("cylindrical_map", &EleStruct::cylindrical_map)
      // EleStruct.gen_grad_map (1D_PTR_type - Used to define E/M fields.
      .def_property_readonly("gen_grad_map", &EleStruct::gen_grad_map)
      // EleStruct.grid_field (1D_PTR_type - Used to define E/M fields. The difference between
      // map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the 1st order
      // spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference
      // particle which is generally the zero orbit (non-zero, for example, in the second slice of a
      // sliced wiggler).
      .def_property_readonly("grid_field", &EleStruct::grid_field)
      // EleStruct.map_ref_orb_in (0D_NOT_type - Entrance end transfer map ref orbit
      .def_property("map_ref_orb_in", &EleStruct::map_ref_orb_in, &EleStruct::set_map_ref_orb_in)
      // EleStruct.map_ref_orb_out (0D_NOT_type - Exit end transfer map ref orbit
      .def_property("map_ref_orb_out", &EleStruct::map_ref_orb_out, &EleStruct::set_map_ref_orb_out)
      // EleStruct.time_ref_orb_in (0D_NOT_type - Reference orbit at entrance end for ref_time calc.
      .def_property("time_ref_orb_in", &EleStruct::time_ref_orb_in, &EleStruct::set_time_ref_orb_in)
      // EleStruct.time_ref_orb_out (0D_NOT_type - Reference orbit at exit end for ref_time calc.
      .def_property(
          "time_ref_orb_out",
          &EleStruct::time_ref_orb_out,
          &EleStruct::set_time_ref_orb_out
      )
      // EleStruct.value (1D_NOT_real - attribute values.
      .def_property("value", &EleStruct::value, &EleStruct::set_value)
      // EleStruct.old_value (1D_NOT_real - Used to see if %value(:) array has changed. Note: The
      // reference orbit for spin/orbit matrices is %map_ref_orb_in/out
      .def_property("old_value", &EleStruct::old_value, &EleStruct::set_old_value)
      // EleStruct.spin_q (2D_NOT_real - 0th and 1st order Spin transport quaternion.
      .def_property("spin_q", &EleStruct::spin_q, &EleStruct::set_spin_q)
      // EleStruct.vec0 (1D_NOT_real - 0th order transport vector.
      .def_property("vec0", &EleStruct::vec0, &EleStruct::set_vec0)
      // EleStruct.mat6 (2D_NOT_real - 1st order transport matrix.
      .def_property("mat6", &EleStruct::mat6, &EleStruct::set_mat6)
      // EleStruct.c_mat (2D_NOT_real - 2x2 C coupling matrix
      .def_property("c_mat", &EleStruct::c_mat, &EleStruct::set_c_mat)
      // EleStruct.dc_mat_dpz (2D_NOT_real - d(c_mat)/dpz variation.
      .def_property("dc_mat_dpz", &EleStruct::dc_mat_dpz, &EleStruct::set_dc_mat_dpz)
      // EleStruct.gamma_c (0D_NOT_real - gamma associated with C matrix
      .def_property("gamma_c", &EleStruct::gamma_c, &EleStruct::set_gamma_c)
      // EleStruct.s_start (0D_NOT_real - longitudinal ref position at entrance_end
      .def_property("s_start", &EleStruct::s_start, &EleStruct::set_s_start)
      // EleStruct.s (0D_NOT_real - longitudinal ref position at the exit end.
      .def_property("s", &EleStruct::s, &EleStruct::set_s)
      // EleStruct.ref_time (0D_NOT_real - Time ref particle passes exit end.
      .def_property("ref_time", &EleStruct::ref_time, &EleStruct::set_ref_time)
      // EleStruct.a_pole (1D_PTR_real - knl for multipole elements.
      .def_property("a_pole", &EleStruct::a_pole, &EleStruct::set_a_pole)
      // EleStruct.b_pole (1D_PTR_real - tilt for multipole elements.
      .def_property("b_pole", &EleStruct::b_pole, &EleStruct::set_b_pole)
      // EleStruct.a_pole_elec (1D_PTR_real - Electrostatic multipoles. ksnl for multipole elements.
      .def_property("a_pole_elec", &EleStruct::a_pole_elec, &EleStruct::set_a_pole_elec)
      // EleStruct.b_pole_elec (1D_PTR_real - Electrostatic multipoles.
      .def_property("b_pole_elec", &EleStruct::b_pole_elec, &EleStruct::set_b_pole_elec)
      // EleStruct.custom (1D_PTR_real - Custom attributes.
      .def_property("custom", &EleStruct::custom, &EleStruct::set_custom)
      // EleStruct.r (3D_PTR_real - For general use. Not used by Bmad.
      .def_property("r", &EleStruct::r, &EleStruct::set_r)
      // EleStruct.key (0D_NOT_integer - Element class (quadrupole, etc.).
      .def_property("key", &EleStruct::key, &EleStruct::set_key)
      // EleStruct.sub_key (0D_NOT_integer - Records bend input type.
      .def_property("sub_key", &EleStruct::sub_key, &EleStruct::set_sub_key)
      // EleStruct.ix_ele (0D_NOT_integer - Index in branch ele(0:) array. Set to ix_slice_slave$ =
      // -2 for slice_slave$ elements.
      .def_property("ix_ele", &EleStruct::ix_ele, &EleStruct::set_ix_ele)
      // EleStruct.ix_branch (0D_NOT_integer - Index in lat%branch(:) array. Note: lat%ele =>
      // lat%branch(0).
      .def_property("ix_branch", &EleStruct::ix_branch, &EleStruct::set_ix_branch)
      // EleStruct.lord_status (0D_NOT_integer - Type of lord element this is. overlay_lord$, etc.
      .def_property("lord_status", &EleStruct::lord_status, &EleStruct::set_lord_status)
      // EleStruct.n_slave (0D_NOT_integer - Number of slaves (except field overlap slaves) of this
      // element.
      .def_property("n_slave", &EleStruct::n_slave, &EleStruct::set_n_slave)
      // EleStruct.n_slave_field (0D_NOT_integer - Number of field slaves of this element.
      .def_property("n_slave_field", &EleStruct::n_slave_field, &EleStruct::set_n_slave_field)
      // EleStruct.ix1_slave (0D_NOT_integer - Pointer index to this element's slaves.
      .def_property("ix1_slave", &EleStruct::ix1_slave, &EleStruct::set_ix1_slave)
      // EleStruct.slave_status (0D_NOT_integer - Type of slave element this is. multipass_slave$,
      // slice_slave$, etc.
      .def_property("slave_status", &EleStruct::slave_status, &EleStruct::set_slave_status)
      // EleStruct.n_lord (0D_NOT_integer - Number of lords (except field overlap and ramper lords).
      .def_property("n_lord", &EleStruct::n_lord, &EleStruct::set_n_lord)
      // EleStruct.n_lord_field (0D_NOT_integer - Number of field lords of this element.
      .def_property("n_lord_field", &EleStruct::n_lord_field, &EleStruct::set_n_lord_field)
      // EleStruct.n_lord_ramper (0D_NOT_integer - Number of ramper lords.
      .def_property("n_lord_ramper", &EleStruct::n_lord_ramper, &EleStruct::set_n_lord_ramper)
      // EleStruct.ic1_lord (0D_NOT_integer - Pointer index to this element's lords.
      .def_property("ic1_lord", &EleStruct::ic1_lord, &EleStruct::set_ic1_lord)
      // EleStruct.ix_pointer (0D_NOT_integer - For general use. Not used by Bmad.
      .def_property("ix_pointer", &EleStruct::ix_pointer, &EleStruct::set_ix_pointer)
      // EleStruct.ixx (0D_NOT_integer - Index for Bmad internal use.
      .def_property("ixx", &EleStruct::ixx, &EleStruct::set_ixx)
      // EleStruct.iyy (0D_NOT_integer - Index for Bmad internal use.
      .def_property("iyy", &EleStruct::iyy, &EleStruct::set_iyy)
      // EleStruct.izz (0D_NOT_integer - Index for Bmad internal use.
      .def_property("izz", &EleStruct::izz, &EleStruct::set_izz)
      // EleStruct.mat6_calc_method (0D_NOT_integer - taylor$, symp_lie_ptc$, etc.
      .def_property(
          "mat6_calc_method",
          &EleStruct::mat6_calc_method,
          &EleStruct::set_mat6_calc_method
      )
      // EleStruct.tracking_method (0D_NOT_integer - taylor$, linear$, etc.
      .def_property("tracking_method", &EleStruct::tracking_method, &EleStruct::set_tracking_method)
      // EleStruct.spin_tracking_method (0D_NOT_integer - symp_lie_ptc$, etc.
      .def_property(
          "spin_tracking_method",
          &EleStruct::spin_tracking_method,
          &EleStruct::set_spin_tracking_method
      )
      // EleStruct.csr_method (0D_NOT_integer - or one_dim$ ('1_dim'), steady_state_3d$
      .def_property("csr_method", &EleStruct::csr_method, &EleStruct::set_csr_method)
      // EleStruct.space_charge_method (0D_NOT_integer - slice$, slice_longitudinal$,
      // slice_transverse$, fft_3D$, cathode_fft_3d$
      .def_property(
          "space_charge_method",
          &EleStruct::space_charge_method,
          &EleStruct::set_space_charge_method
      )
      // EleStruct.ptc_integration_type (0D_NOT_integer - drift_kick$, matrix_kick$, or ripken_kick$
      .def_property(
          "ptc_integration_type",
          &EleStruct::ptc_integration_type,
          &EleStruct::set_ptc_integration_type
      )
      // EleStruct.field_calc (0D_NOT_integer - no_field$, fieldmap$, refer_to_lords$, or custom$
      .def_property("field_calc", &EleStruct::field_calc, &EleStruct::set_field_calc)
      // EleStruct.aperture_at (0D_NOT_integer - Aperture location: entrance_end$, ...
      .def_property("aperture_at", &EleStruct::aperture_at, &EleStruct::set_aperture_at)
      // EleStruct.aperture_type (0D_NOT_integer - rectangular$, elliptical$, auto_aperture$, ...
      .def_property("aperture_type", &EleStruct::aperture_type, &EleStruct::set_aperture_type)
      // EleStruct.ref_species (0D_NOT_integer - Reference species
      .def_property("ref_species", &EleStruct::ref_species, &EleStruct::set_ref_species)
      // EleStruct.orientation (0D_NOT_integer - -1 -> Element is longitudinally reversed. +1 ->
      // Normal.
      .def_property("orientation", &EleStruct::orientation, &EleStruct::set_orientation)
      // EleStruct.symplectify (0D_NOT_logical - Symplectify mat6 matrices.
      .def_property("symplectify", &EleStruct::symplectify, &EleStruct::set_symplectify)
      // EleStruct.mode_flip (0D_NOT_logical - Have the normal modes traded places?
      .def_property("mode_flip", &EleStruct::mode_flip, &EleStruct::set_mode_flip)
      // EleStruct.multipoles_on (0D_NOT_logical - For turning multipoles on/off
      .def_property("multipoles_on", &EleStruct::multipoles_on, &EleStruct::set_multipoles_on)
      // EleStruct.scale_multipoles (0D_NOT_logical - Are ab_multipoles within other elements (EG:
      // quads, etc.) scaled by the strength of the element?
      .def_property(
          "scale_multipoles",
          &EleStruct::scale_multipoles,
          &EleStruct::set_scale_multipoles
      )
      // EleStruct.taylor_map_includes_offsets (0D_NOT_logical - Taylor map calculated with element
      // misalignments?
      .def_property(
          "taylor_map_includes_offsets",
          &EleStruct::taylor_map_includes_offsets,
          &EleStruct::set_taylor_map_includes_offsets
      )
      // EleStruct.field_master (0D_NOT_logical - Calculate strength from the field value?
      .def_property("field_master", &EleStruct::field_master, &EleStruct::set_field_master)
      // EleStruct.is_on (0D_NOT_logical - For turning element on/off.
      .def_property("is_on", &EleStruct::is_on, &EleStruct::set_is_on)
      // EleStruct.logic (0D_NOT_logical - For general use. Not used by Bmad (except during lattice
      // parsing).
      .def_property("logic", &EleStruct::logic, &EleStruct::set_logic)
      // EleStruct.bmad_logic (0D_NOT_logical - For Bmad internal use only.
      .def_property("bmad_logic", &EleStruct::bmad_logic, &EleStruct::set_bmad_logic)
      // EleStruct.select (0D_NOT_logical - For Bmad internal use only.
      .def_property("select", &EleStruct::select, &EleStruct::set_select)
      // EleStruct.offset_moves_aperture (0D_NOT_logical - element offsets affects aperture? ! final
      // :: ele_finalizer
      .def_property(
          "offset_moves_aperture",
          &EleStruct::offset_moves_aperture,
          &EleStruct::set_offset_moves_aperture
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EleStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const EleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EleStruct &self) {
            return EleStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const EleStruct &self, py::dict &memo) { return EleStruct(self); })

      ;

  bind_FTypeArrayND<EleStructArray1D>(m, "EleStructArray1D");
  bind_FTypeAlloc1D<EleStructAlloc1D>(m, "EleStructAlloc1D");
  // 2D EleStruct arrays are not used in structs/routines
  // 3D EleStruct arrays are not used in structs/routines
}

// =============================================================================
// ellipse_beam_init_struct
void init_ellipse_beam_init_struct(py::module &m, py::class_<EllipseBeamInitStruct> &cls) {
  cls.def(
         py::init<std::optional<int>, std::optional<int>, std::optional<double>>(),
         py::arg("part_per_ellipse") = py::none(),
         py::arg("n_ellipse") = py::none(),
         py::arg("sigma_cutoff") = py::none()
  )
      // EllipseBeamInitStruct.part_per_ellipse (0D_NOT_integer - number of particles per ellipse
      .def_property(
          "part_per_ellipse",
          &EllipseBeamInitStruct::part_per_ellipse,
          &EllipseBeamInitStruct::set_part_per_ellipse
      )
      // EllipseBeamInitStruct.n_ellipse (0D_NOT_integer - number of ellipses (>= 1)
      .def_property(
          "n_ellipse",
          &EllipseBeamInitStruct::n_ellipse,
          &EllipseBeamInitStruct::set_n_ellipse
      )
      // EllipseBeamInitStruct.sigma_cutoff (0D_NOT_real - sigma cutoff of the representation
      .def_property(
          "sigma_cutoff",
          &EllipseBeamInitStruct::sigma_cutoff,
          &EllipseBeamInitStruct::set_sigma_cutoff
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EllipseBeamInitStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const EllipseBeamInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EllipseBeamInitStruct &self) {
            return EllipseBeamInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const EllipseBeamInitStruct &self, py::dict &memo) {
            return EllipseBeamInitStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<EllipseBeamInitStructArray1D>(m, "EllipseBeamInitStructArray1D");
  bind_FTypeAlloc1D<EllipseBeamInitStructAlloc1D>(m, "EllipseBeamInitStructAlloc1D");
  // 2D EllipseBeamInitStruct arrays are not used in structs/routines
  // 3D EllipseBeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// em_field_struct
void init_em_field_struct(py::module &m, py::class_<EmFieldStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<double>>>,
             optional_ref<const std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const std::vector<double>>>(),
         py::arg("E") = py::none(),
         py::arg("B") = py::none(),
         py::arg("dE") = py::none(),
         py::arg("dB") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("phi_B") = py::none(),
         py::arg("A") = py::none()
  )
      // EmFieldStruct.E (1D_NOT_real - electric field.
      .def_property("E", &EmFieldStruct::E, &EmFieldStruct::set_E)
      // EmFieldStruct.B (1D_NOT_real - magnetic field.
      .def_property("B", &EmFieldStruct::B, &EmFieldStruct::set_B)
      // EmFieldStruct.dE (2D_NOT_real - electric field gradient.
      .def_property("dE", &EmFieldStruct::dE, &EmFieldStruct::set_dE)
      // EmFieldStruct.dB (2D_NOT_real - magnetic field gradient.
      .def_property("dB", &EmFieldStruct::dB, &EmFieldStruct::set_dB)
      // EmFieldStruct.phi (0D_NOT_real - Electric scalar potential.
      .def_property("phi", &EmFieldStruct::phi, &EmFieldStruct::set_phi)
      // EmFieldStruct.phi_B (0D_NOT_real - Magnetic scalar potential.
      .def_property("phi_B", &EmFieldStruct::phi_B, &EmFieldStruct::set_phi_B)
      // EmFieldStruct.A (1D_NOT_real - Magnetic vector potential.
      .def_property("A", &EmFieldStruct::A, &EmFieldStruct::set_A)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmFieldStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const EmFieldStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmFieldStruct &self) {
            return EmFieldStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const EmFieldStruct &self, py::dict &memo) { return EmFieldStruct(self); }
      )

      ;

  bind_FTypeArrayND<EmFieldStructArray1D>(m, "EmFieldStructArray1D");
  bind_FTypeAlloc1D<EmFieldStructAlloc1D>(m, "EmFieldStructAlloc1D");
  // 2D EmFieldStruct arrays are not used in structs/routines
  // 3D EmFieldStruct arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_struct
void init_em_taylor_struct(py::module &m, py::class_<EmTaylorStruct> &cls) {
  cls.def(py::init<std::optional<double>>(), py::arg("ref") = py::none())
      // EmTaylorStruct.ref (0D_NOT_real -
      .def_property("ref", &EmTaylorStruct::ref, &EmTaylorStruct::set_ref)
      // EmTaylorStruct.term (1D_ALLOC_type -
      .def_property_readonly("term", &EmTaylorStruct::term)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmTaylorStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const EmTaylorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorStruct &self) {
            return EmTaylorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const EmTaylorStruct &self, py::dict &memo) { return EmTaylorStruct(self); }
      )

      ;

  bind_FTypeArrayND<EmTaylorStructArray1D>(m, "EmTaylorStructArray1D");
  bind_FTypeAlloc1D<EmTaylorStructAlloc1D>(m, "EmTaylorStructAlloc1D");
  // 2D EmTaylorStruct arrays are not used in structs/routines
  // 3D EmTaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_term_struct
void init_em_taylor_term_struct(py::module &m, py::class_<EmTaylorTermStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, optional_ref<const std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("expn") = py::none()
  )
      // EmTaylorTermStruct.coef (0D_NOT_real -
      .def_property("coef", &EmTaylorTermStruct::coef, &EmTaylorTermStruct::set_coef)
      // EmTaylorTermStruct.expn (1D_NOT_integer -
      .def_property("expn", &EmTaylorTermStruct::expn, &EmTaylorTermStruct::set_expn)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return EmTaylorTermStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const EmTaylorTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EmTaylorTermStruct &self) {
            return EmTaylorTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const EmTaylorTermStruct &self, py::dict &memo) { return EmTaylorTermStruct(self); }
      )

      ;

  bind_FTypeArrayND<EmTaylorTermStructArray1D>(m, "EmTaylorTermStructArray1D");
  bind_FTypeAlloc1D<EmTaylorTermStructAlloc1D>(m, "EmTaylorTermStructAlloc1D");
  // 2D EmTaylorTermStruct arrays are not used in structs/routines
  // 3D EmTaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// expression_atom_struct
void init_expression_atom_struct(py::module &m, py::class_<ExpressionAtomStruct> &cls) {
  cls.def(
         py::init<optional_ref<const std::string>, std::optional<int>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("value") = py::none()
  )
      // ExpressionAtomStruct.name (0D_NOT_character -
      .def_property("name", &ExpressionAtomStruct::name, &ExpressionAtomStruct::set_name)
      // ExpressionAtomStruct.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc. To convert to
      // string use: expression_op_name
      .def_property("type", &ExpressionAtomStruct::type, &ExpressionAtomStruct::set_type)
      // ExpressionAtomStruct.value (0D_NOT_real -
      .def_property("value", &ExpressionAtomStruct::value, &ExpressionAtomStruct::set_value)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ExpressionAtomStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ExpressionAtomStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionAtomStruct &self) {
            return ExpressionAtomStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ExpressionAtomStruct &self, py::dict &memo) {
            return ExpressionAtomStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<ExpressionAtomStructArray1D>(m, "ExpressionAtomStructArray1D");
  bind_FTypeAlloc1D<ExpressionAtomStructAlloc1D>(m, "ExpressionAtomStructAlloc1D");
  // 2D ExpressionAtomStruct arrays are not used in structs/routines
  // 3D ExpressionAtomStruct arrays are not used in structs/routines
}

// =============================================================================
// expression_tree_struct
void init_expression_tree_struct(py::module &m, py::class_<ExpressionTreeStruct> &cls) {
  cls.def(
         py::init<optional_ref<const std::string>, std::optional<int>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("value") = py::none()
  )
      // ExpressionTreeStruct.name (0D_NOT_character -
      .def_property("name", &ExpressionTreeStruct::name, &ExpressionTreeStruct::set_name)
      // ExpressionTreeStruct.type (0D_NOT_integer - plus$, minum$, sin$, cos$, etc.
      .def_property("type", &ExpressionTreeStruct::type, &ExpressionTreeStruct::set_type)
      // ExpressionTreeStruct.value (0D_NOT_real -
      .def_property("value", &ExpressionTreeStruct::value, &ExpressionTreeStruct::set_value)
      // ExpressionTreeStruct.node (1D_PTR_type - Child nodes. Note: Pointer used here since Ifort
      // does not support allocatable.
      .def_property_readonly("node", &ExpressionTreeStruct::node)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return ExpressionTreeStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const ExpressionTreeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExpressionTreeStruct &self) {
            return ExpressionTreeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ExpressionTreeStruct &self, py::dict &memo) {
            return ExpressionTreeStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<ExpressionTreeStructArray1D>(m, "ExpressionTreeStructArray1D");
  bind_FTypeAlloc1D<ExpressionTreeStructAlloc1D>(m, "ExpressionTreeStructAlloc1D");
  // 2D ExpressionTreeStruct arrays are not used in structs/routines
  // 3D ExpressionTreeStruct arrays are not used in structs/routines
}