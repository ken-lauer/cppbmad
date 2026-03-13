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
      .def_property(
          "ele",
          py::cpp_function(&ElePointerStruct::ele, py::keep_alive<0, 1>()),
          &ElePointerStruct::set_ele
      )
      .def_property(
          "loc",
          py::cpp_function(&ElePointerStruct::loc, py::keep_alive<0, 1>()),
          &ElePointerStruct::set_loc
      )
      .def_property(
          "id",
          &ElePointerStruct::id,
          &ElePointerStruct::set_id,
          "For general use. Not used by Bmad. In particular, used by Tao to designate universe ele "
          "is in."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ElePointerStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ElePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<ElePointerStructArray1D, ElePointerStructAlloc1D>(
      m,
      "ElePointerStructArray1D",
      "ElePointerStructAlloc1D"
  );
  // 2D ElePointerStruct arrays are not used in structs/routines
  // 3D ElePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// ele_struct
void init_ele_struct(py::module &m, py::class_<EleStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
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
             std::optional<std::vector<double>>,
             optional_ref<const WakeStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             optional_ref<const CoordStruct>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<std::vector<double>>>>,
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
      .def_property("name", &EleStruct::name, &EleStruct::set_name, "name of element.")
      .def_property("type", &EleStruct::type, &EleStruct::set_type, "type name.")
      .def_property("alias", &EleStruct::alias, &EleStruct::set_alias, "Another name.")
      .def_property(
          "component_name",
          &EleStruct::component_name,
          &EleStruct::set_component_name,
          "Used by overlays, multipass patch, etc."
      )
      .def_property("descrip", &EleStruct::descrip, &EleStruct::set_descrip, "Description string.")
      .def_property(
          "a",
          py::cpp_function(&EleStruct::a, py::keep_alive<0, 1>()),
          &EleStruct::set_a,
          "Twiss parameters at end of element"
      )
      .def_property(
          "b",
          py::cpp_function(&EleStruct::b, py::keep_alive<0, 1>()),
          &EleStruct::set_b,
          "Twiss parameters at end of element"
      )
      .def_property(
          "z",
          py::cpp_function(&EleStruct::z, py::keep_alive<0, 1>()),
          &EleStruct::set_z,
          "Twiss parameters at end of element"
      )
      .def_property(
          "x",
          py::cpp_function(&EleStruct::x, py::keep_alive<0, 1>()),
          &EleStruct::set_x,
          "Projected dispersions."
      )
      .def_property(
          "y",
          py::cpp_function(&EleStruct::y, py::keep_alive<0, 1>()),
          &EleStruct::set_y,
          "Projected dispersions."
      )
      .def_property(
          "ac_kick",
          py::cpp_function(&EleStruct::ac_kick, py::keep_alive<0, 1>()),
          &EleStruct::set_ac_kick,
          "ac_kicker element parameters."
      )
      .def_property(
          "bookkeeping_state",
          py::cpp_function(&EleStruct::bookkeeping_state, py::keep_alive<0, 1>()),
          &EleStruct::set_bookkeeping_state,
          "Attribute bookkeeping"
      )
      .def_property(
          "branch",
          py::cpp_function(&EleStruct::branch, py::keep_alive<0, 1>()),
          &EleStruct::set_branch,
          "Pointer to branch containing element."
      )
      .def_property(
          "control",
          py::cpp_function(&EleStruct::control, py::keep_alive<0, 1>()),
          &EleStruct::set_control,
          "group & overlay variables."
      )
      .def_property(
          "rf",
          py::cpp_function(&EleStruct::rf, py::keep_alive<0, 1>()),
          &EleStruct::set_rf,
          "RF parameters."
      )
      .def_property(
          "lord",
          py::cpp_function(&EleStruct::lord, py::keep_alive<0, 1>()),
          &EleStruct::set_lord,
          "Pointer to a slice lord."
      )
      .def_property(
          "ptc_fibre",
          py::cpp_function(&EleStruct::ptc_fibre, py::keep_alive<0, 1>()),
          &EleStruct::set_ptc_fibre,
          "PTC track corresponding to this ele. %floor is floor coord of lab coordinates at the "
          "downstream end. Notice that if ele%direction = -1, the lab coords have the z-axis "
          "antiparallel to the +s-direction."
      )
      .def_property(
          "floor",
          py::cpp_function(&EleStruct::floor, py::keep_alive<0, 1>()),
          &EleStruct::set_floor
      )
      .def_property(
          "high_energy_space_charge",
          py::cpp_function(&EleStruct::high_energy_space_charge, py::keep_alive<0, 1>()),
          &EleStruct::set_high_energy_space_charge
      )
      .def_property(
          "mode3",
          py::cpp_function(&EleStruct::mode3, py::keep_alive<0, 1>()),
          &EleStruct::set_mode3,
          "6D normal mode structure."
      )
      .def_property(
          "photon",
          py::cpp_function(&EleStruct::photon, py::keep_alive<0, 1>()),
          &EleStruct::set_photon
      )
      .def_property(
          "rad_map",
          py::cpp_function(&EleStruct::rad_map, py::keep_alive<0, 1>()),
          &EleStruct::set_rad_map,
          "Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are "
          "not necessarily the same. For example, Sprint spin Taylor maps can be with respect to "
          "the zero orbit independent of the orbital map."
      )
      .def_property_readonly(
          "taylor",
          py::cpp_function(&EleStruct::taylor, py::keep_alive<0, 1>()),
          "Phase space Taylor map."
      )
      .def_property(
          "spin_taylor_ref_orb_in",
          py::cpp_function(&EleStruct::spin_taylor_ref_orb_in, py::keep_alive<0, 1>()),
          &EleStruct::set_spin_taylor_ref_orb_in
      )
      .def_property_readonly(
          "spin_taylor",
          py::cpp_function(&EleStruct::spin_taylor, py::keep_alive<0, 1>()),
          "Quaternion Spin Taylor map."
      )
      .def_property(
          "wake",
          py::cpp_function(&EleStruct::wake, py::keep_alive<0, 1>()),
          &EleStruct::set_wake,
          "Wakes"
      )
      .def_property_readonly(
          "wall3d",
          py::cpp_function(&EleStruct::wall3d, py::keep_alive<0, 1>()),
          "Chamber or capillary wall E/M field structs."
      )
      .def_property_readonly(
          "cartesian_map",
          py::cpp_function(&EleStruct::cartesian_map, py::keep_alive<0, 1>()),
          "Used to define E/M fields"
      )
      .def_property_readonly(
          "cylindrical_map",
          py::cpp_function(&EleStruct::cylindrical_map, py::keep_alive<0, 1>()),
          "Used to define E/M fields"
      )
      .def_property_readonly(
          "gen_grad_map",
          py::cpp_function(&EleStruct::gen_grad_map, py::keep_alive<0, 1>()),
          "Used to define E/M fields."
      )
      .def_property_readonly(
          "grid_field",
          py::cpp_function(&EleStruct::grid_field, py::keep_alive<0, 1>()),
          "Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that "
          "map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, "
          "is non-zero while time_ref_orb follows the reference particle which is generally the "
          "zero orbit (non-zero, for example, in the second slice of a sliced wiggler)."
      )
      .def_property(
          "map_ref_orb_in",
          py::cpp_function(&EleStruct::map_ref_orb_in, py::keep_alive<0, 1>()),
          &EleStruct::set_map_ref_orb_in,
          "Entrance end transfer map ref orbit"
      )
      .def_property(
          "map_ref_orb_out",
          py::cpp_function(&EleStruct::map_ref_orb_out, py::keep_alive<0, 1>()),
          &EleStruct::set_map_ref_orb_out,
          "Exit end transfer map ref orbit"
      )
      .def_property(
          "time_ref_orb_in",
          py::cpp_function(&EleStruct::time_ref_orb_in, py::keep_alive<0, 1>()),
          &EleStruct::set_time_ref_orb_in,
          "Reference orbit at entrance end for ref_time calc."
      )
      .def_property(
          "time_ref_orb_out",
          py::cpp_function(&EleStruct::time_ref_orb_out, py::keep_alive<0, 1>()),
          &EleStruct::set_time_ref_orb_out,
          "Reference orbit at exit end for ref_time calc."
      )
      .def_property(
          "value",
          py::cpp_function(&EleStruct::value, py::keep_alive<0, 1>()),
          &EleStruct::set_value,
          "attribute values."
      )
      .def_property(
          "old_value",
          py::cpp_function(&EleStruct::old_value, py::keep_alive<0, 1>()),
          &EleStruct::set_old_value,
          "Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit "
          "matrices is %map_ref_orb_in/out"
      )
      .def_property(
          "spin_q",
          py::cpp_function(&EleStruct::spin_q, py::keep_alive<0, 1>()),
          &EleStruct::set_spin_q,
          "0th and 1st order Spin transport quaternion."
      )
      .def_property(
          "vec0",
          py::cpp_function(&EleStruct::vec0, py::keep_alive<0, 1>()),
          &EleStruct::set_vec0,
          "0th order transport vector."
      )
      .def_property(
          "mat6",
          py::cpp_function(&EleStruct::mat6, py::keep_alive<0, 1>()),
          &EleStruct::set_mat6,
          "1st order transport matrix."
      )
      .def_property(
          "c_mat",
          py::cpp_function(&EleStruct::c_mat, py::keep_alive<0, 1>()),
          &EleStruct::set_c_mat,
          "2x2 C coupling matrix"
      )
      .def_property(
          "dc_mat_dpz",
          py::cpp_function(&EleStruct::dc_mat_dpz, py::keep_alive<0, 1>()),
          &EleStruct::set_dc_mat_dpz,
          "d(c_mat)/dpz variation."
      )
      .def_property(
          "gamma_c",
          &EleStruct::gamma_c,
          &EleStruct::set_gamma_c,
          "gamma associated with C matrix"
      )
      .def_property(
          "s_start",
          &EleStruct::s_start,
          &EleStruct::set_s_start,
          "longitudinal ref position at entrance_end"
      )
      .def_property(
          "s",
          &EleStruct::s,
          &EleStruct::set_s,
          "longitudinal ref position at the exit end."
      )
      .def_property(
          "ref_time",
          &EleStruct::ref_time,
          &EleStruct::set_ref_time,
          "Time ref particle passes exit end."
      )
      .def_property(
          "a_pole",
          py::cpp_function(&EleStruct::a_pole, py::keep_alive<0, 1>()),
          &EleStruct::set_a_pole,
          "knl for multipole elements."
      )
      .def_property(
          "b_pole",
          py::cpp_function(&EleStruct::b_pole, py::keep_alive<0, 1>()),
          &EleStruct::set_b_pole,
          "tilt for multipole elements."
      )
      .def_property(
          "a_pole_elec",
          py::cpp_function(&EleStruct::a_pole_elec, py::keep_alive<0, 1>()),
          &EleStruct::set_a_pole_elec,
          "Electrostatic multipoles. ksnl for multipole elements."
      )
      .def_property(
          "b_pole_elec",
          py::cpp_function(&EleStruct::b_pole_elec, py::keep_alive<0, 1>()),
          &EleStruct::set_b_pole_elec,
          "Electrostatic multipoles."
      )
      .def_property(
          "custom",
          py::cpp_function(&EleStruct::custom, py::keep_alive<0, 1>()),
          &EleStruct::set_custom,
          "Custom attributes."
      )
      .def_property(
          "r",
          py::cpp_function(&EleStruct::r, py::keep_alive<0, 1>()),
          &EleStruct::set_r,
          "For general use. Not used by Bmad."
      )
      .def_property(
          "key",
          &EleStruct::key,
          &EleStruct::set_key,
          "Element class (quadrupole, etc.)."
      )
      .def_property(
          "sub_key",
          &EleStruct::sub_key,
          &EleStruct::set_sub_key,
          "Records bend input type."
      )
      .def_property(
          "ix_ele",
          &EleStruct::ix_ele,
          &EleStruct::set_ix_ele,
          "Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements."
      )
      .def_property(
          "ix_branch",
          &EleStruct::ix_branch,
          &EleStruct::set_ix_branch,
          "Index in lat%branch(:) array. Note: lat%ele => lat%branch(0)."
      )
      .def_property(
          "lord_status",
          &EleStruct::lord_status,
          &EleStruct::set_lord_status,
          "Type of lord element this is. overlay_lord$, etc."
      )
      .def_property(
          "n_slave",
          &EleStruct::n_slave,
          &EleStruct::set_n_slave,
          "Number of slaves (except field overlap slaves) of this element."
      )
      .def_property(
          "n_slave_field",
          &EleStruct::n_slave_field,
          &EleStruct::set_n_slave_field,
          "Number of field slaves of this element."
      )
      .def_property(
          "ix1_slave",
          &EleStruct::ix1_slave,
          &EleStruct::set_ix1_slave,
          "Pointer index to this element's slaves."
      )
      .def_property(
          "slave_status",
          &EleStruct::slave_status,
          &EleStruct::set_slave_status,
          "Type of slave element this is. multipass_slave$, slice_slave$, etc."
      )
      .def_property(
          "n_lord",
          &EleStruct::n_lord,
          &EleStruct::set_n_lord,
          "Number of lords (except field overlap and ramper lords)."
      )
      .def_property(
          "n_lord_field",
          &EleStruct::n_lord_field,
          &EleStruct::set_n_lord_field,
          "Number of field lords of this element."
      )
      .def_property(
          "n_lord_ramper",
          &EleStruct::n_lord_ramper,
          &EleStruct::set_n_lord_ramper,
          "Number of ramper lords."
      )
      .def_property(
          "ic1_lord",
          &EleStruct::ic1_lord,
          &EleStruct::set_ic1_lord,
          "Pointer index to this element's lords."
      )
      .def_property(
          "ix_pointer",
          &EleStruct::ix_pointer,
          &EleStruct::set_ix_pointer,
          "For general use. Not used by Bmad."
      )
      .def_property("ixx", &EleStruct::ixx, &EleStruct::set_ixx, "Index for Bmad internal use.")
      .def_property("iyy", &EleStruct::iyy, &EleStruct::set_iyy, "Index for Bmad internal use.")
      .def_property("izz", &EleStruct::izz, &EleStruct::set_izz, "Index for Bmad internal use.")
      .def_property(
          "mat6_calc_method",
          &EleStruct::mat6_calc_method,
          &EleStruct::set_mat6_calc_method,
          "taylor$, symp_lie_ptc$, etc."
      )
      .def_property(
          "tracking_method",
          &EleStruct::tracking_method,
          &EleStruct::set_tracking_method,
          "taylor$, linear$, etc."
      )
      .def_property(
          "spin_tracking_method",
          &EleStruct::spin_tracking_method,
          &EleStruct::set_spin_tracking_method,
          "symp_lie_ptc$, etc."
      )
      .def_property(
          "csr_method",
          &EleStruct::csr_method,
          &EleStruct::set_csr_method,
          "or one_dim$ ('1_dim'), steady_state_3d$"
      )
      .def_property(
          "space_charge_method",
          &EleStruct::space_charge_method,
          &EleStruct::set_space_charge_method,
          "slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$"
      )
      .def_property(
          "ptc_integration_type",
          &EleStruct::ptc_integration_type,
          &EleStruct::set_ptc_integration_type,
          "drift_kick$, matrix_kick$, or ripken_kick$"
      )
      .def_property(
          "field_calc",
          &EleStruct::field_calc,
          &EleStruct::set_field_calc,
          "no_field$, fieldmap$, refer_to_lords$, or custom$"
      )
      .def_property(
          "aperture_at",
          &EleStruct::aperture_at,
          &EleStruct::set_aperture_at,
          "Aperture location: entrance_end$, ..."
      )
      .def_property(
          "aperture_type",
          &EleStruct::aperture_type,
          &EleStruct::set_aperture_type,
          "rectangular$, elliptical$, auto_aperture$, ..."
      )
      .def_property(
          "ref_species",
          &EleStruct::ref_species,
          &EleStruct::set_ref_species,
          "Reference species"
      )
      .def_property(
          "orientation",
          &EleStruct::orientation,
          &EleStruct::set_orientation,
          "-1 -> Element is longitudinally reversed. +1 -> Normal."
      )
      .def_property(
          "symplectify",
          &EleStruct::symplectify,
          &EleStruct::set_symplectify,
          "Symplectify mat6 matrices."
      )
      .def_property(
          "mode_flip",
          &EleStruct::mode_flip,
          &EleStruct::set_mode_flip,
          "Have the normal modes traded places?"
      )
      .def_property(
          "multipoles_on",
          &EleStruct::multipoles_on,
          &EleStruct::set_multipoles_on,
          "For turning multipoles on/off"
      )
      .def_property(
          "scale_multipoles",
          &EleStruct::scale_multipoles,
          &EleStruct::set_scale_multipoles,
          "Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the "
          "element?"
      )
      .def_property(
          "taylor_map_includes_offsets",
          &EleStruct::taylor_map_includes_offsets,
          &EleStruct::set_taylor_map_includes_offsets,
          "Taylor map calculated with element misalignments?"
      )
      .def_property(
          "field_master",
          &EleStruct::field_master,
          &EleStruct::set_field_master,
          "Calculate strength from the field value?"
      )
      .def_property(
          "is_on",
          &EleStruct::is_on,
          &EleStruct::set_is_on,
          "For turning element on/off."
      )
      .def_property(
          "logic",
          &EleStruct::logic,
          &EleStruct::set_logic,
          "For general use. Not used by Bmad (except during lattice parsing)."
      )
      .def_property(
          "bmad_logic",
          &EleStruct::bmad_logic,
          &EleStruct::set_bmad_logic,
          "For Bmad internal use only."
      )
      .def_property(
          "select",
          &EleStruct::select,
          &EleStruct::set_select,
          "For Bmad internal use only."
      )
      .def_property(
          "offset_moves_aperture",
          &EleStruct::offset_moves_aperture,
          &EleStruct::set_offset_moves_aperture,
          "element offsets affects aperture? ! final :: ele_finalizer"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EleStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<EleStructArray1D, EleStructAlloc1D>(
      m,
      "EleStructArray1D",
      "EleStructAlloc1D"
  );
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
      .def_property(
          "part_per_ellipse",
          &EllipseBeamInitStruct::part_per_ellipse,
          &EllipseBeamInitStruct::set_part_per_ellipse,
          "number of particles per ellipse"
      )
      .def_property(
          "n_ellipse",
          &EllipseBeamInitStruct::n_ellipse,
          &EllipseBeamInitStruct::set_n_ellipse,
          "number of ellipses (>= 1)"
      )
      .def_property(
          "sigma_cutoff",
          &EllipseBeamInitStruct::sigma_cutoff,
          &EllipseBeamInitStruct::set_sigma_cutoff,
          "sigma cutoff of the representation"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EllipseBeamInitStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EllipseBeamInitStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<EllipseBeamInitStructArray1D, EllipseBeamInitStructAlloc1D>(
      m,
      "EllipseBeamInitStructArray1D",
      "EllipseBeamInitStructAlloc1D"
  );
  // 2D EllipseBeamInitStruct arrays are not used in structs/routines
  // 3D EllipseBeamInitStruct arrays are not used in structs/routines
}

// =============================================================================
// em_field_struct
void init_em_field_struct(py::module &m, py::class_<EmFieldStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         py::arg("E") = py::none(),
         py::arg("B") = py::none(),
         py::arg("dE") = py::none(),
         py::arg("dB") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("phi_B") = py::none(),
         py::arg("A") = py::none()
  )
      .def_property(
          "E",
          py::cpp_function(&EmFieldStruct::E, py::keep_alive<0, 1>()),
          &EmFieldStruct::set_E,
          "electric field."
      )
      .def_property(
          "B",
          py::cpp_function(&EmFieldStruct::B, py::keep_alive<0, 1>()),
          &EmFieldStruct::set_B,
          "magnetic field."
      )
      .def_property(
          "dE",
          py::cpp_function(&EmFieldStruct::dE, py::keep_alive<0, 1>()),
          &EmFieldStruct::set_dE,
          "electric field gradient."
      )
      .def_property(
          "dB",
          py::cpp_function(&EmFieldStruct::dB, py::keep_alive<0, 1>()),
          &EmFieldStruct::set_dB,
          "magnetic field gradient."
      )
      .def_property(
          "phi",
          &EmFieldStruct::phi,
          &EmFieldStruct::set_phi,
          "Electric scalar potential."
      )
      .def_property(
          "phi_B",
          &EmFieldStruct::phi_B,
          &EmFieldStruct::set_phi_B,
          "Magnetic scalar potential."
      )
      .def_property(
          "A",
          py::cpp_function(&EmFieldStruct::A, py::keep_alive<0, 1>()),
          &EmFieldStruct::set_A,
          "Magnetic vector potential."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EmFieldStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmFieldStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<EmFieldStructArray1D, EmFieldStructAlloc1D>(
      m,
      "EmFieldStructArray1D",
      "EmFieldStructAlloc1D"
  );
  // 2D EmFieldStruct arrays are not used in structs/routines
  // 3D EmFieldStruct arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_struct
void init_em_taylor_struct(py::module &m, py::class_<EmTaylorStruct> &cls) {
  cls.def(py::init<std::optional<double>>(), py::arg("ref") = py::none())
      .def_property("ref", &EmTaylorStruct::ref, &EmTaylorStruct::set_ref)
      .def_property_readonly(
          "term",
          py::cpp_function(&EmTaylorStruct::term, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EmTaylorStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<EmTaylorStructArray1D, EmTaylorStructAlloc1D>(
      m,
      "EmTaylorStructArray1D",
      "EmTaylorStructAlloc1D"
  );
  // 2D EmTaylorStruct arrays are not used in structs/routines
  // 3D EmTaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// em_taylor_term_struct
void init_em_taylor_term_struct(py::module &m, py::class_<EmTaylorTermStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("expn") = py::none()
  )
      .def_property("coef", &EmTaylorTermStruct::coef, &EmTaylorTermStruct::set_coef)
      .def_property(
          "expn",
          py::cpp_function(&EmTaylorTermStruct::expn, py::keep_alive<0, 1>()),
          &EmTaylorTermStruct::set_expn
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EmTaylorTermStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<EmTaylorTermStructArray1D, EmTaylorTermStructAlloc1D>(
      m,
      "EmTaylorTermStructArray1D",
      "EmTaylorTermStructAlloc1D"
  );
  // 2D EmTaylorTermStruct arrays are not used in structs/routines
  // 3D EmTaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// expression_atom_struct
void init_expression_atom_struct(py::module &m, py::class_<ExpressionAtomStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("value") = py::none()
  )
      .def_property("name", &ExpressionAtomStruct::name, &ExpressionAtomStruct::set_name)
      .def_property(
          "type",
          &ExpressionAtomStruct::type,
          &ExpressionAtomStruct::set_type,
          "plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name"
      )
      .def_property("value", &ExpressionAtomStruct::value, &ExpressionAtomStruct::set_value)
      .def_static(
          "new_array1d",
          [](int sz) { return ExpressionAtomStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ExpressionAtomStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<ExpressionAtomStructArray1D, ExpressionAtomStructAlloc1D>(
      m,
      "ExpressionAtomStructArray1D",
      "ExpressionAtomStructAlloc1D"
  );
  // 2D ExpressionAtomStruct arrays are not used in structs/routines
  // 3D ExpressionAtomStruct arrays are not used in structs/routines
}

// =============================================================================
// expression_tree_struct
void init_expression_tree_struct(py::module &m, py::class_<ExpressionTreeStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("type") = py::none(),
         py::arg("value") = py::none()
  )
      .def_property("name", &ExpressionTreeStruct::name, &ExpressionTreeStruct::set_name)
      .def_property(
          "type",
          &ExpressionTreeStruct::type,
          &ExpressionTreeStruct::set_type,
          "plus$, minum$, sin$, cos$, etc."
      )
      .def_property("value", &ExpressionTreeStruct::value, &ExpressionTreeStruct::set_value)
      .def_property_readonly(
          "node",
          py::cpp_function(&ExpressionTreeStruct::node, py::keep_alive<0, 1>()),
          "Child nodes. Note: Pointer used here since Ifort does not support allocatable."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ExpressionTreeStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ExpressionTreeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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

  bind_1d_type_array_pair<ExpressionTreeStructArray1D, ExpressionTreeStructAlloc1D>(
      m,
      "ExpressionTreeStructArray1D",
      "ExpressionTreeStructAlloc1D"
  );
  // 2D ExpressionTreeStruct arrays are not used in structs/routines
  // 3D ExpressionTreeStruct arrays are not used in structs/routines
}