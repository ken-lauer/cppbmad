#include "pybmad/generated/structs_e.hpp"

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
// ele_attribute_struct
void init_ele_attribute_struct(nb::module_ &m, nb::class_<EleAttributeStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<double>>(),
         nb::arg("name") = nb::none(),
         nb::arg("state") = nb::none(),
         nb::arg("kind") = nb::none(),
         nb::arg("units") = nb::none(),
         nb::arg("ix_attrib") = nb::none(),
         nb::arg("value") = nb::none()
  )
      .def_prop_rw("name", &EleAttributeStruct::name, &EleAttributeStruct::set_name)
      .def_prop_rw(
          "state",
          &EleAttributeStruct::state,
          &EleAttributeStruct::set_state,
          "See above."
      )
      .def_prop_rw(
          "kind",
          &EleAttributeStruct::kind,
          &EleAttributeStruct::set_kind,
          "Is_switch$, is_real$, etc. See attribute_type routine."
      )
      .def_prop_rw(
          "units",
          &EleAttributeStruct::units,
          &EleAttributeStruct::set_units,
          "EG: 'T*m'."
      )
      .def_prop_rw(
          "ix_attrib",
          &EleAttributeStruct::ix_attrib,
          &EleAttributeStruct::set_ix_attrib,
          "Attribute index. Frequently will be where in the ele%value(:) array the attribute is."
      )
      .def_prop_rw(
          "value",
          &EleAttributeStruct::value,
          &EleAttributeStruct::set_value,
          "Used by type_ele."
      )

      .def("__repr__", [](const EleAttributeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EleAttributeStruct &self) {
            return EleAttributeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const EleAttributeStruct &self, nb::dict &memo) { return EleAttributeStruct(self); }
      )
      .def(
          "__eq__",
          [](const EleAttributeStruct &self, const EleAttributeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EleAttributeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D EleAttributeStruct arrays are not used in structs/routines
  // 2D EleAttributeStruct arrays are not used in structs/routines
  // 3D EleAttributeStruct arrays are not used in structs/routines
}

// =============================================================================
// ele_pointer_struct
void init_ele_pointer_struct(nb::module_ &m, nb::class_<ElePointerStruct> &cls) {
  cls.def(
         "__init__",
         [](ElePointerStruct *self,
            const EleStruct *ele,
            const LatEleLocStruct *loc,
            std::optional<int> id) {
           new (self) ElePointerStruct(ptr_to_opt_ref(ele), ptr_to_opt_ref(loc), id);
         },
         nb::arg("ele") = nb::none(),
         nb::arg("loc") = nb::none(),
         nb::arg("id") = nb::none()
  )
      .def_prop_rw(
          "ele",
          &ElePointerStruct::ele,
          &ElePointerStruct::set_ele,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "loc",
          &ElePointerStruct::loc,
          &ElePointerStruct::set_loc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "id",
          &ElePointerStruct::id,
          &ElePointerStruct::set_id,
          "For general use. Not used by Bmad. In particular, used by Tao to designate universe ele "
          "is in."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ElePointerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ElePointerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const ElePointerStruct &self, nb::dict &memo) { return ElePointerStruct(self); }
      )
      .def(
          "__eq__",
          [](const ElePointerStruct &self, const ElePointerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ElePointerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<ElePointerStructArray1D, ElePointerStructAlloc1D>(
      m,
      "ElePointerStructArray1D",
      "ElePointerStructAlloc1D"
  );
  bind_FTypeArrayND<ElePointerStructArray2D>(m, "ElePointerStructArray2D");
  // 3D ElePointerStruct arrays are not used in structs/routines
}

// =============================================================================
// ele_struct
void init_ele_struct(nb::module_ &m, nb::class_<EleStruct> &cls) {
  cls.def(
         "__init__",
         [](EleStruct *self,
            std::optional<std::string> name,
            std::optional<std::string> type,
            std::optional<std::string> alias,
            std::optional<std::string> component_name,
            std::optional<std::string> descrip,
            const TwissStruct *a,
            const TwissStruct *b,
            const TwissStruct *z,
            const XyDispStruct *x,
            const XyDispStruct *y,
            const AcKickerStruct *ac_kick,
            const BookkeepingStateStruct *bookkeeping_state,
            const BranchStruct *branch,
            const ControllerStruct *control,
            const ConverterStruct *converter,
            const RfEleStruct *rf,
            const FoilStruct *foil,
            const EleStruct *lord,
            const FloorPositionStruct *floor,
            const HighEnergySpaceChargeStruct *high_energy_space_charge,
            const Mode3Struct *mode3,
            const PhotonElementStruct *photon,
            const RadMapEleStruct *rad_map,
            std::optional<std::vector<double>> spin_taylor_ref_orb_in,
            const WakeStruct *wake,
            const CoordStruct *map_ref_orb_in,
            const CoordStruct *map_ref_orb_out,
            const CoordStruct *time_ref_orb_in,
            const CoordStruct *time_ref_orb_out,
            std::optional<std::vector<double>> value,
            std::optional<std::vector<double>> old_value,
            std::optional<std::vector<std::vector<double>>> spin_q,
            std::optional<std::vector<double>> vec0,
            std::optional<std::vector<std::vector<double>>> mat6,
            std::optional<std::vector<std::vector<double>>> c_mat,
            std::optional<std::vector<std::vector<double>>> dc_mat_dpz,
            std::optional<double> gamma_c,
            std::optional<double> s_start,
            std::optional<double> s,
            std::optional<double> ref_time,
            std::optional<std::vector<double>> a_pole,
            std::optional<std::vector<double>> b_pole,
            std::optional<std::vector<double>> a_pole_elec,
            std::optional<std::vector<double>> b_pole_elec,
            std::optional<std::vector<double>> custom,
            std::optional<std::vector<std::vector<std::vector<double>>>> r,
            std::optional<int> key,
            std::optional<int> sub_key,
            std::optional<int> ix_ele,
            std::optional<int> ix_branch,
            std::optional<int> lord_status,
            std::optional<int> n_slave,
            std::optional<int> n_slave_field,
            std::optional<int> ix1_slave,
            std::optional<int> slave_status,
            std::optional<int> n_lord,
            std::optional<int> n_lord_field,
            std::optional<int> n_lord_ramper,
            std::optional<int> ic1_lord,
            std::optional<int> ix_pointer,
            std::optional<int> ixx,
            std::optional<int> iyy,
            std::optional<int> izz,
            std::optional<int> mat6_calc_method,
            std::optional<int> tracking_method,
            std::optional<int> spin_tracking_method,
            std::optional<int> csr_method,
            std::optional<int> space_charge_method,
            std::optional<int> ptc_integration_type,
            std::optional<int> field_calc,
            std::optional<int> aperture_at,
            std::optional<int> aperture_type,
            std::optional<int> ref_species,
            std::optional<int> orientation,
            std::optional<bool> symplectify,
            std::optional<bool> mode_flip,
            std::optional<bool> multipoles_on,
            std::optional<bool> scale_multipoles,
            std::optional<bool> taylor_map_includes_offsets,
            std::optional<bool> field_master,
            std::optional<bool> is_on,
            std::optional<bool> logic,
            std::optional<bool> bmad_logic,
            std::optional<bool> select,
            std::optional<bool> offset_moves_aperture) {
           new (self) EleStruct(
               name,
               type,
               alias,
               component_name,
               descrip,
               ptr_to_opt_ref(a),
               ptr_to_opt_ref(b),
               ptr_to_opt_ref(z),
               ptr_to_opt_ref(x),
               ptr_to_opt_ref(y),
               ptr_to_opt_ref(ac_kick),
               ptr_to_opt_ref(bookkeeping_state),
               ptr_to_opt_ref(branch),
               ptr_to_opt_ref(control),
               ptr_to_opt_ref(converter),
               ptr_to_opt_ref(rf),
               ptr_to_opt_ref(foil),
               ptr_to_opt_ref(lord),
               ptr_to_opt_ref(floor),
               ptr_to_opt_ref(high_energy_space_charge),
               ptr_to_opt_ref(mode3),
               ptr_to_opt_ref(photon),
               ptr_to_opt_ref(rad_map),
               spin_taylor_ref_orb_in,
               ptr_to_opt_ref(wake),
               ptr_to_opt_ref(map_ref_orb_in),
               ptr_to_opt_ref(map_ref_orb_out),
               ptr_to_opt_ref(time_ref_orb_in),
               ptr_to_opt_ref(time_ref_orb_out),
               value,
               old_value,
               spin_q,
               vec0,
               mat6,
               c_mat,
               dc_mat_dpz,
               gamma_c,
               s_start,
               s,
               ref_time,
               a_pole,
               b_pole,
               a_pole_elec,
               b_pole_elec,
               custom,
               r,
               key,
               sub_key,
               ix_ele,
               ix_branch,
               lord_status,
               n_slave,
               n_slave_field,
               ix1_slave,
               slave_status,
               n_lord,
               n_lord_field,
               n_lord_ramper,
               ic1_lord,
               ix_pointer,
               ixx,
               iyy,
               izz,
               mat6_calc_method,
               tracking_method,
               spin_tracking_method,
               csr_method,
               space_charge_method,
               ptc_integration_type,
               field_calc,
               aperture_at,
               aperture_type,
               ref_species,
               orientation,
               symplectify,
               mode_flip,
               multipoles_on,
               scale_multipoles,
               taylor_map_includes_offsets,
               field_master,
               is_on,
               logic,
               bmad_logic,
               select,
               offset_moves_aperture
           );
         },
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("alias") = nb::none(),
         nb::arg("component_name") = nb::none(),
         nb::arg("descrip") = nb::none(),
         nb::arg("a") = nb::none(),
         nb::arg("b") = nb::none(),
         nb::arg("z") = nb::none(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("ac_kick") = nb::none(),
         nb::arg("bookkeeping_state") = nb::none(),
         nb::arg("branch") = nb::none(),
         nb::arg("control") = nb::none(),
         nb::arg("converter") = nb::none(),
         nb::arg("rf") = nb::none(),
         nb::arg("foil") = nb::none(),
         nb::arg("lord") = nb::none(),
         nb::arg("floor") = nb::none(),
         nb::arg("high_energy_space_charge") = nb::none(),
         nb::arg("mode3") = nb::none(),
         nb::arg("photon") = nb::none(),
         nb::arg("rad_map") = nb::none(),
         nb::arg("spin_taylor_ref_orb_in") = nb::none(),
         nb::arg("wake") = nb::none(),
         nb::arg("map_ref_orb_in") = nb::none(),
         nb::arg("map_ref_orb_out") = nb::none(),
         nb::arg("time_ref_orb_in") = nb::none(),
         nb::arg("time_ref_orb_out") = nb::none(),
         nb::arg("value") = nb::none(),
         nb::arg("old_value") = nb::none(),
         nb::arg("spin_q") = nb::none(),
         nb::arg("vec0") = nb::none(),
         nb::arg("mat6") = nb::none(),
         nb::arg("c_mat") = nb::none(),
         nb::arg("dc_mat_dpz") = nb::none(),
         nb::arg("gamma_c") = nb::none(),
         nb::arg("s_start") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("ref_time") = nb::none(),
         nb::arg("a_pole") = nb::none(),
         nb::arg("b_pole") = nb::none(),
         nb::arg("a_pole_elec") = nb::none(),
         nb::arg("b_pole_elec") = nb::none(),
         nb::arg("custom") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("key") = nb::none(),
         nb::arg("sub_key") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("lord_status") = nb::none(),
         nb::arg("n_slave") = nb::none(),
         nb::arg("n_slave_field") = nb::none(),
         nb::arg("ix1_slave") = nb::none(),
         nb::arg("slave_status") = nb::none(),
         nb::arg("n_lord") = nb::none(),
         nb::arg("n_lord_field") = nb::none(),
         nb::arg("n_lord_ramper") = nb::none(),
         nb::arg("ic1_lord") = nb::none(),
         nb::arg("ix_pointer") = nb::none(),
         nb::arg("ixx") = nb::none(),
         nb::arg("iyy") = nb::none(),
         nb::arg("izz") = nb::none(),
         nb::arg("mat6_calc_method") = nb::none(),
         nb::arg("tracking_method") = nb::none(),
         nb::arg("spin_tracking_method") = nb::none(),
         nb::arg("csr_method") = nb::none(),
         nb::arg("space_charge_method") = nb::none(),
         nb::arg("ptc_integration_type") = nb::none(),
         nb::arg("field_calc") = nb::none(),
         nb::arg("aperture_at") = nb::none(),
         nb::arg("aperture_type") = nb::none(),
         nb::arg("ref_species") = nb::none(),
         nb::arg("orientation") = nb::none(),
         nb::arg("symplectify") = nb::none(),
         nb::arg("mode_flip") = nb::none(),
         nb::arg("multipoles_on") = nb::none(),
         nb::arg("scale_multipoles") = nb::none(),
         nb::arg("taylor_map_includes_offsets") = nb::none(),
         nb::arg("field_master") = nb::none(),
         nb::arg("is_on") = nb::none(),
         nb::arg("logic") = nb::none(),
         nb::arg("bmad_logic") = nb::none(),
         nb::arg("select") = nb::none(),
         nb::arg("offset_moves_aperture") = nb::none()
  )
      .def_prop_rw("name", &EleStruct::name, &EleStruct::set_name, "name of element.")
      .def_prop_rw("type", &EleStruct::type, &EleStruct::set_type, "type name.")
      .def_prop_rw("alias", &EleStruct::alias, &EleStruct::set_alias, "Another name.")
      .def_prop_rw(
          "component_name",
          &EleStruct::component_name,
          &EleStruct::set_component_name,
          "Used by overlays, multipass patch, etc."
      )
      .def_prop_rw("descrip", &EleStruct::descrip, &EleStruct::set_descrip, "Description string.")
      .def_prop_rw(
          "a",
          &EleStruct::a,
          &EleStruct::set_a,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Twiss parameters at end of element"
      )
      .def_prop_rw(
          "b",
          &EleStruct::b,
          &EleStruct::set_b,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Twiss parameters at end of element"
      )
      .def_prop_rw(
          "z",
          &EleStruct::z,
          &EleStruct::set_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Twiss parameters at end of element"
      )
      .def_prop_rw(
          "x",
          &EleStruct::x,
          &EleStruct::set_x,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Projected dispersions."
      )
      .def_prop_rw(
          "y",
          &EleStruct::y,
          &EleStruct::set_y,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Projected dispersions."
      )
      .def_prop_rw(
          "ac_kick",
          &EleStruct::ac_kick,
          &EleStruct::set_ac_kick,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "ac_kicker element parameters."
      )
      .def_prop_rw(
          "bookkeeping_state",
          &EleStruct::bookkeeping_state,
          &EleStruct::set_bookkeeping_state,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Attribute bookkeeping"
      )
      .def_prop_rw(
          "branch",
          &EleStruct::branch,
          &EleStruct::set_branch,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Pointer to branch containing element."
      )
      .def_prop_rw(
          "control",
          &EleStruct::control,
          &EleStruct::set_control,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "group & overlay variables."
      )
      .def_prop_rw(
          "converter",
          &EleStruct::converter,
          &EleStruct::set_converter,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "EG: Positron converter in linac."
      )
      .def_prop_rw(
          "rf",
          &EleStruct::rf,
          &EleStruct::set_rf,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "RF parameters."
      )
      .def_prop_rw(
          "foil",
          &EleStruct::foil,
          &EleStruct::set_foil,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "lord",
          &EleStruct::lord,
          &EleStruct::set_lord,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Pointer to a slice lord."
      )
      .def_prop_rw(
          "floor",
          &EleStruct::floor,
          &EleStruct::set_floor,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "high_energy_space_charge",
          &EleStruct::high_energy_space_charge,
          &EleStruct::set_high_energy_space_charge,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "mode3",
          &EleStruct::mode3,
          &EleStruct::set_mode3,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "6D normal mode structure."
      )
      .def_prop_rw(
          "photon",
          &EleStruct::photon,
          &EleStruct::set_photon,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      // 0D_ALLOC_type multipole_cache proxy support missing
      .def_prop_rw(
          "rad_map",
          &EleStruct::rad_map,
          &EleStruct::set_rad_map,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are "
          "not necessarily the same. For example, Sprint spin Taylor maps can be with respect to "
          "the zero orbit independent of the orbital map."
      )
      .def_prop_ro("taylor", &EleStruct::taylor, nb::keep_alive<0, 1>(), "Phase space Taylor map.")
      .def_prop_rw(
          "spin_taylor_ref_orb_in",
          &EleStruct::spin_taylor_ref_orb_in,
          &EleStruct::set_spin_taylor_ref_orb_in,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro(
          "spin_taylor",
          &EleStruct::spin_taylor,
          nb::keep_alive<0, 1>(),
          "Quaternion Spin Taylor map."
      )
      .def_prop_rw(
          "wake",
          &EleStruct::wake,
          &EleStruct::set_wake,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Wakes"
      )
      .def_prop_ro(
          "wall3d",
          &EleStruct::wall3d,
          nb::keep_alive<0, 1>(),
          "Chamber or capillary wall E/M field structs."
      )
      .def_prop_ro(
          "cartesian_map",
          &EleStruct::cartesian_map,
          nb::keep_alive<0, 1>(),
          "Used to define E/M fields"
      )
      .def_prop_ro(
          "cylindrical_map",
          &EleStruct::cylindrical_map,
          nb::keep_alive<0, 1>(),
          "Used to define E/M fields"
      )
      .def_prop_ro(
          "gen_grad_map",
          &EleStruct::gen_grad_map,
          nb::keep_alive<0, 1>(),
          "Used to define E/M fields."
      )
      .def_prop_ro(
          "grid_field",
          &EleStruct::grid_field,
          nb::keep_alive<0, 1>(),
          "Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that "
          "map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, "
          "is non-zero while time_ref_orb follows the reference particle which is generally the "
          "zero orbit (non-zero, for example, in the second slice of a sliced wiggler)."
      )
      .def_prop_rw(
          "map_ref_orb_in",
          &EleStruct::map_ref_orb_in,
          &EleStruct::set_map_ref_orb_in,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Entrance end transfer map ref orbit"
      )
      .def_prop_rw(
          "map_ref_orb_out",
          &EleStruct::map_ref_orb_out,
          &EleStruct::set_map_ref_orb_out,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Exit end transfer map ref orbit"
      )
      .def_prop_rw(
          "time_ref_orb_in",
          &EleStruct::time_ref_orb_in,
          &EleStruct::set_time_ref_orb_in,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Reference orbit at entrance end for ref_time calc."
      )
      .def_prop_rw(
          "time_ref_orb_out",
          &EleStruct::time_ref_orb_out,
          &EleStruct::set_time_ref_orb_out,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Reference orbit at exit end for ref_time calc."
      )
      .def_prop_rw(
          "value",
          &EleStruct::value,
          &EleStruct::set_value,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "attribute values."
      )
      .def_prop_rw(
          "old_value",
          &EleStruct::old_value,
          &EleStruct::set_old_value,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit "
          "matrices is %map_ref_orb_in/out"
      )
      .def_prop_rw(
          "spin_q",
          &EleStruct::spin_q,
          &EleStruct::set_spin_q,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th and 1st order Spin transport quaternion."
      )
      .def_prop_rw(
          "vec0",
          &EleStruct::vec0,
          &EleStruct::set_vec0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th order transport vector."
      )
      .def_prop_rw(
          "mat6",
          &EleStruct::mat6,
          &EleStruct::set_mat6,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1st order transport matrix."
      )
      .def_prop_rw(
          "c_mat",
          &EleStruct::c_mat,
          &EleStruct::set_c_mat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "2x2 C coupling matrix"
      )
      .def_prop_rw(
          "dc_mat_dpz",
          &EleStruct::dc_mat_dpz,
          &EleStruct::set_dc_mat_dpz,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "d(c_mat)/dpz variation."
      )
      .def_prop_rw(
          "gamma_c",
          &EleStruct::gamma_c,
          &EleStruct::set_gamma_c,
          "gamma associated with C matrix"
      )
      .def_prop_rw(
          "s_start",
          &EleStruct::s_start,
          &EleStruct::set_s_start,
          "longitudinal ref position at entrance_end"
      )
      .def_prop_rw(
          "s",
          &EleStruct::s,
          &EleStruct::set_s,
          "longitudinal ref position at the exit end."
      )
      .def_prop_rw(
          "ref_time",
          &EleStruct::ref_time,
          &EleStruct::set_ref_time,
          "Time ref particle passes exit end."
      )
      .def_prop_rw(
          "a_pole",
          &EleStruct::a_pole,
          &EleStruct::set_a_pole,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "knl for multipole elements."
      )
      .def_prop_rw(
          "b_pole",
          &EleStruct::b_pole,
          &EleStruct::set_b_pole,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "tilt for multipole elements."
      )
      .def_prop_rw(
          "a_pole_elec",
          &EleStruct::a_pole_elec,
          &EleStruct::set_a_pole_elec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Electrostatic multipoles. ksnl for multipole elements."
      )
      .def_prop_rw(
          "b_pole_elec",
          &EleStruct::b_pole_elec,
          &EleStruct::set_b_pole_elec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Electrostatic multipoles."
      )
      .def_prop_rw(
          "custom",
          &EleStruct::custom,
          &EleStruct::set_custom,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Custom attributes."
      )
      .def_prop_rw(
          "r",
          &EleStruct::r,
          &EleStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "For general use. Not used by Bmad."
      )
      .def_prop_rw("key", &EleStruct::key, &EleStruct::set_key, "Element class (quadrupole, etc.).")
      .def_prop_rw(
          "sub_key",
          &EleStruct::sub_key,
          &EleStruct::set_sub_key,
          "Records bend input type."
      )
      .def_prop_rw(
          "ix_ele",
          &EleStruct::ix_ele,
          &EleStruct::set_ix_ele,
          "Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements."
      )
      .def_prop_rw(
          "ix_branch",
          &EleStruct::ix_branch,
          &EleStruct::set_ix_branch,
          "Index in lat%branch(:) array. Note: lat%ele => lat%branch(0)."
      )
      .def_prop_rw(
          "lord_status",
          &EleStruct::lord_status,
          &EleStruct::set_lord_status,
          "Type of lord element this is. overlay_lord$, etc."
      )
      .def_prop_rw(
          "n_slave",
          &EleStruct::n_slave,
          &EleStruct::set_n_slave,
          "Number of slaves (except field overlap slaves) of this element."
      )
      .def_prop_rw(
          "n_slave_field",
          &EleStruct::n_slave_field,
          &EleStruct::set_n_slave_field,
          "Number of field slaves of this element."
      )
      .def_prop_rw(
          "ix1_slave",
          &EleStruct::ix1_slave,
          &EleStruct::set_ix1_slave,
          "Pointer index to this element's slaves."
      )
      .def_prop_rw(
          "slave_status",
          &EleStruct::slave_status,
          &EleStruct::set_slave_status,
          "Type of slave element this is. multipass_slave$, slice_slave$, etc."
      )
      .def_prop_rw(
          "n_lord",
          &EleStruct::n_lord,
          &EleStruct::set_n_lord,
          "Number of lords (except field overlap and ramper lords)."
      )
      .def_prop_rw(
          "n_lord_field",
          &EleStruct::n_lord_field,
          &EleStruct::set_n_lord_field,
          "Number of field lords of this element."
      )
      .def_prop_rw(
          "n_lord_ramper",
          &EleStruct::n_lord_ramper,
          &EleStruct::set_n_lord_ramper,
          "Number of ramper lords."
      )
      .def_prop_rw(
          "ic1_lord",
          &EleStruct::ic1_lord,
          &EleStruct::set_ic1_lord,
          "Pointer index to this element's lords."
      )
      .def_prop_rw(
          "ix_pointer",
          &EleStruct::ix_pointer,
          &EleStruct::set_ix_pointer,
          "For general use. Not used by Bmad."
      )
      .def_prop_rw("ixx", &EleStruct::ixx, &EleStruct::set_ixx, "Index for Bmad internal use.")
      .def_prop_rw("iyy", &EleStruct::iyy, &EleStruct::set_iyy, "Index for Bmad internal use.")
      .def_prop_rw("izz", &EleStruct::izz, &EleStruct::set_izz, "Index for Bmad internal use.")
      .def_prop_rw(
          "mat6_calc_method",
          &EleStruct::mat6_calc_method,
          &EleStruct::set_mat6_calc_method,
          "taylor$, symp_lie_ptc$, etc."
      )
      .def_prop_rw(
          "tracking_method",
          &EleStruct::tracking_method,
          &EleStruct::set_tracking_method,
          "taylor$, linear$, etc."
      )
      .def_prop_rw(
          "spin_tracking_method",
          &EleStruct::spin_tracking_method,
          &EleStruct::set_spin_tracking_method,
          "symp_lie_ptc$, etc."
      )
      .def_prop_rw(
          "csr_method",
          &EleStruct::csr_method,
          &EleStruct::set_csr_method,
          "or one_dim$ ('1_dim'), steady_state_3d$"
      )
      .def_prop_rw(
          "space_charge_method",
          &EleStruct::space_charge_method,
          &EleStruct::set_space_charge_method,
          "slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$"
      )
      .def_prop_rw(
          "ptc_integration_type",
          &EleStruct::ptc_integration_type,
          &EleStruct::set_ptc_integration_type,
          "drift_kick$, matrix_kick$, or ripken_kick$"
      )
      .def_prop_rw(
          "field_calc",
          &EleStruct::field_calc,
          &EleStruct::set_field_calc,
          "no_field$, fieldmap$, refer_to_lords$, or custom$"
      )
      .def_prop_rw(
          "aperture_at",
          &EleStruct::aperture_at,
          &EleStruct::set_aperture_at,
          "Aperture location: entrance_end$, ..."
      )
      .def_prop_rw(
          "aperture_type",
          &EleStruct::aperture_type,
          &EleStruct::set_aperture_type,
          "rectangular$, elliptical$, auto_aperture$, ..."
      )
      .def_prop_rw(
          "ref_species",
          &EleStruct::ref_species,
          &EleStruct::set_ref_species,
          "Reference species"
      )
      .def_prop_rw(
          "orientation",
          &EleStruct::orientation,
          &EleStruct::set_orientation,
          "-1 -> Element is longitudinally reversed. +1 -> Normal."
      )
      .def_prop_rw(
          "symplectify",
          &EleStruct::symplectify,
          &EleStruct::set_symplectify,
          "Symplectify mat6 matrices."
      )
      .def_prop_rw(
          "mode_flip",
          &EleStruct::mode_flip,
          &EleStruct::set_mode_flip,
          "Have the normal modes traded places?"
      )
      .def_prop_rw(
          "multipoles_on",
          &EleStruct::multipoles_on,
          &EleStruct::set_multipoles_on,
          "For turning multipoles on/off"
      )
      .def_prop_rw(
          "scale_multipoles",
          &EleStruct::scale_multipoles,
          &EleStruct::set_scale_multipoles,
          "Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the "
          "element?"
      )
      .def_prop_rw(
          "taylor_map_includes_offsets",
          &EleStruct::taylor_map_includes_offsets,
          &EleStruct::set_taylor_map_includes_offsets,
          "Taylor map calculated with element misalignments?"
      )
      .def_prop_rw(
          "field_master",
          &EleStruct::field_master,
          &EleStruct::set_field_master,
          "Calculate strength from the field value?"
      )
      .def_prop_rw("is_on", &EleStruct::is_on, &EleStruct::set_is_on, "For turning element on/off.")
      .def_prop_rw(
          "logic",
          &EleStruct::logic,
          &EleStruct::set_logic,
          "For general use. Not used by Bmad (except during lattice parsing)."
      )
      .def_prop_rw(
          "bmad_logic",
          &EleStruct::bmad_logic,
          &EleStruct::set_bmad_logic,
          "For Bmad internal use only."
      )
      .def_prop_rw(
          "select",
          &EleStruct::select,
          &EleStruct::set_select,
          "For Bmad internal use only."
      )
      .def_prop_rw(
          "offset_moves_aperture",
          &EleStruct::offset_moves_aperture,
          &EleStruct::set_offset_moves_aperture,
          "element offsets affects aperture? ! final :: ele_finalizer"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EleStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const EleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const EleStruct &self) {
            return EleStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const EleStruct &self, nb::dict &memo) { return EleStruct(self); })
      .def(
          "__eq__",
          [](const EleStruct &self, const EleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EleStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

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
void init_ellipse_beam_init_struct(nb::module_ &m, nb::class_<EllipseBeamInitStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<int>, std::optional<double>>(),
         nb::arg("part_per_ellipse") = nb::none(),
         nb::arg("n_ellipse") = nb::none(),
         nb::arg("sigma_cutoff") = nb::none()
  )
      .def_prop_rw(
          "part_per_ellipse",
          &EllipseBeamInitStruct::part_per_ellipse,
          &EllipseBeamInitStruct::set_part_per_ellipse,
          "number of particles per ellipse"
      )
      .def_prop_rw(
          "n_ellipse",
          &EllipseBeamInitStruct::n_ellipse,
          &EllipseBeamInitStruct::set_n_ellipse,
          "number of ellipses (>= 1)"
      )
      .def_prop_rw(
          "sigma_cutoff",
          &EllipseBeamInitStruct::sigma_cutoff,
          &EllipseBeamInitStruct::set_sigma_cutoff,
          "sigma cutoff of the representation"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EllipseBeamInitStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EllipseBeamInitStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const EllipseBeamInitStruct &self, nb::dict &memo) {
            return EllipseBeamInitStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const EllipseBeamInitStruct &self, const EllipseBeamInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EllipseBeamInitStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_em_field_struct(nb::module_ &m, nb::class_<EmFieldStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         nb::arg("E") = nb::none(),
         nb::arg("B") = nb::none(),
         nb::arg("dE") = nb::none(),
         nb::arg("dB") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("phi_B") = nb::none(),
         nb::arg("A") = nb::none()
  )
      .def_prop_rw(
          "E",
          &EmFieldStruct::E,
          &EmFieldStruct::set_E,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "electric field."
      )
      .def_prop_rw(
          "B",
          &EmFieldStruct::B,
          &EmFieldStruct::set_B,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "magnetic field."
      )
      .def_prop_rw(
          "dE",
          &EmFieldStruct::dE,
          &EmFieldStruct::set_dE,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "electric field gradient."
      )
      .def_prop_rw(
          "dB",
          &EmFieldStruct::dB,
          &EmFieldStruct::set_dB,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "magnetic field gradient."
      )
      .def_prop_rw(
          "phi",
          &EmFieldStruct::phi,
          &EmFieldStruct::set_phi,
          "Electric scalar potential."
      )
      .def_prop_rw(
          "phi_B",
          &EmFieldStruct::phi_B,
          &EmFieldStruct::set_phi_B,
          "Magnetic scalar potential."
      )
      .def_prop_rw(
          "A",
          &EmFieldStruct::A,
          &EmFieldStruct::set_A,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Magnetic vector potential."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EmFieldStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmFieldStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const EmFieldStruct &self, nb::dict &memo) { return EmFieldStruct(self); }
      )
      .def(
          "__eq__",
          [](const EmFieldStruct &self, const EmFieldStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EmFieldStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
void init_em_taylor_struct(nb::module_ &m, nb::class_<EmTaylorStruct> &cls) {
  cls.def(nb::init<std::optional<double>>(), nb::arg("ref") = nb::none())
      .def_prop_rw("ref", &EmTaylorStruct::ref, &EmTaylorStruct::set_ref)
      .def_prop_ro("term", &EmTaylorStruct::term, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return EmTaylorStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const EmTaylorStruct &self, nb::dict &memo) { return EmTaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const EmTaylorStruct &self, const EmTaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EmTaylorStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
void init_em_taylor_term_struct(nb::module_ &m, nb::class_<EmTaylorTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("expn") = nb::none()
  )
      .def_prop_rw("coef", &EmTaylorTermStruct::coef, &EmTaylorTermStruct::set_coef)
      .def_prop_rw(
          "expn",
          &EmTaylorTermStruct::expn,
          &EmTaylorTermStruct::set_expn,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return EmTaylorTermStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = EmTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const EmTaylorTermStruct &self, nb::dict &memo) { return EmTaylorTermStruct(self); }
      )
      .def(
          "__eq__",
          [](const EmTaylorTermStruct &self, const EmTaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const EmTaylorTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
void init_expression_atom_struct(nb::module_ &m, nb::class_<ExpressionAtomStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>, std::optional<double>>(),
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("value") = nb::none()
  )
      .def_prop_rw("name", &ExpressionAtomStruct::name, &ExpressionAtomStruct::set_name)
      .def_prop_rw(
          "type",
          &ExpressionAtomStruct::type,
          &ExpressionAtomStruct::set_type,
          "plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name"
      )
      .def_prop_rw("value", &ExpressionAtomStruct::value, &ExpressionAtomStruct::set_value)
      .def_static(
          "new_array1d",
          [](int sz) { return ExpressionAtomStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ExpressionAtomStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const ExpressionAtomStruct &self, nb::dict &memo) {
            return ExpressionAtomStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ExpressionAtomStruct &self, const ExpressionAtomStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ExpressionAtomStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_expression_tree_struct(nb::module_ &m, nb::class_<ExpressionTreeStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>, std::optional<double>>(),
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("value") = nb::none()
  )
      .def_prop_rw("name", &ExpressionTreeStruct::name, &ExpressionTreeStruct::set_name)
      .def_prop_rw(
          "type",
          &ExpressionTreeStruct::type,
          &ExpressionTreeStruct::set_type,
          "plus$, minum$, sin$, cos$, etc."
      )
      .def_prop_rw("value", &ExpressionTreeStruct::value, &ExpressionTreeStruct::set_value)
      .def_prop_ro(
          "node",
          &ExpressionTreeStruct::node,
          nb::keep_alive<0, 1>(),
          "Child nodes. Note: Pointer used here since Ifort does not support allocatable."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ExpressionTreeStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ExpressionTreeStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const ExpressionTreeStruct &self, nb::dict &memo) {
            return ExpressionTreeStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ExpressionTreeStruct &self, const ExpressionTreeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ExpressionTreeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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

// =============================================================================
// extra_parsing_info_struct
void init_extra_parsing_info_struct(nb::module_ &m, nb::class_<ExtraParsingInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](ExtraParsingInfoStruct *self,
            const RandomStateStruct *ran_state,
            std::optional<int> ran_seed,
            std::optional<bool> undeterministic_ran_function_called,
            std::optional<bool> d_orb_set,
            std::optional<bool> max_aperture_limit_set,
            std::optional<bool> default_ds_step_set,
            std::optional<bool> significant_length_set,
            std::optional<bool> rel_tol_tracking_set,
            std::optional<bool> abs_tol_tracking_set,
            std::optional<bool> rel_tol_adaptive_tracking_set,
            std::optional<bool> abs_tol_adaptive_tracking_set,
            std::optional<bool> init_ds_adaptive_tracking_set,
            std::optional<bool> min_ds_adaptive_tracking_set,
            std::optional<bool> fatal_ds_adaptive_tracking_set,
            std::optional<bool> synch_rad_scale_set,
            std::optional<bool> autoscale_amp_abs_tol_set,
            std::optional<bool> autoscale_amp_rel_tol_set,
            std::optional<bool> autoscale_phase_tol_set,
            std::optional<bool> rf_phase_below_transition_ref_set,
            std::optional<bool> electric_dipole_moment_set,
            std::optional<bool> taylor_order_set,
            std::optional<bool> runge_kutta_order_set,
            std::optional<bool> default_integ_order_set,
            std::optional<bool> sr_wakes_on_set,
            std::optional<bool> lr_wakes_on_set,
            std::optional<bool> high_energy_space_charge_on_set,
            std::optional<bool> csr_and_space_charge_on_set,
            std::optional<bool> spin_tracking_on_set,
            std::optional<bool> spin_sokolov_ternov_flipping_on_set,
            std::optional<bool> radiation_damping_on_set,
            std::optional<bool> radiation_zero_average_set,
            std::optional<bool> radiation_fluctuations_on_set,
            std::optional<bool> conserve_taylor_maps_set,
            std::optional<bool> absolute_time_tracking_set,
            std::optional<bool> absolute_time_ref_shift_set,
            std::optional<bool> convert_to_kinetic_momentum_set,
            std::optional<bool> aperture_limit_on_set,
            std::optional<bool> normalize_twiss_set,
            std::optional<bool> sad_eps_scale_set,
            std::optional<bool> sad_amp_max_set,
            std::optional<bool> sad_n_div_max_set,
            std::optional<bool> max_num_runge_kutta_step_set,
            std::optional<bool> spin_n0_direction_user_set_set,
            std::optional<bool> debug_set,
            std::optional<bool> ds_track_step_set,
            std::optional<bool> dt_track_step_set,
            std::optional<bool> cathode_strength_cutoff_set,
            std::optional<bool> sc_rel_tol_tracking_set,
            std::optional<bool> sc_abs_tol_tracking_set,
            std::optional<bool> beam_chamber_height_set,
            std::optional<bool> lsc_sigma_cutoff_set,
            std::optional<bool> particle_sigma_cutoff_set,
            std::optional<bool> space_charge_mesh_size_set,
            std::optional<bool> csr3d_mesh_size_set,
            std::optional<bool> n_bin_set,
            std::optional<bool> particle_bin_span_set,
            std::optional<bool> n_shield_images_set,
            std::optional<bool> sc_min_in_bin_set,
            std::optional<bool> lsc_kick_transverse_dependence_set,
            std::optional<bool> sc_debug_set,
            std::optional<bool> diagnostic_output_file_set,
            std::optional<bool> old_integrator_set,
            std::optional<bool> use_orientation_patches_set,
            std::optional<bool> print_info_messages_set,
            std::optional<bool> max_fringe_order_set,
            std::optional<bool> exact_model_set,
            std::optional<bool> exact_misalign_set,
            std::optional<bool> vertical_kick_set,
            std::optional<bool> cut_factor_set,
            std::optional<bool> translate_patch_drift_time_set) {
           new (self) ExtraParsingInfoStruct(
               ptr_to_opt_ref(ran_state),
               ran_seed,
               undeterministic_ran_function_called,
               d_orb_set,
               max_aperture_limit_set,
               default_ds_step_set,
               significant_length_set,
               rel_tol_tracking_set,
               abs_tol_tracking_set,
               rel_tol_adaptive_tracking_set,
               abs_tol_adaptive_tracking_set,
               init_ds_adaptive_tracking_set,
               min_ds_adaptive_tracking_set,
               fatal_ds_adaptive_tracking_set,
               synch_rad_scale_set,
               autoscale_amp_abs_tol_set,
               autoscale_amp_rel_tol_set,
               autoscale_phase_tol_set,
               rf_phase_below_transition_ref_set,
               electric_dipole_moment_set,
               taylor_order_set,
               runge_kutta_order_set,
               default_integ_order_set,
               sr_wakes_on_set,
               lr_wakes_on_set,
               high_energy_space_charge_on_set,
               csr_and_space_charge_on_set,
               spin_tracking_on_set,
               spin_sokolov_ternov_flipping_on_set,
               radiation_damping_on_set,
               radiation_zero_average_set,
               radiation_fluctuations_on_set,
               conserve_taylor_maps_set,
               absolute_time_tracking_set,
               absolute_time_ref_shift_set,
               convert_to_kinetic_momentum_set,
               aperture_limit_on_set,
               normalize_twiss_set,
               sad_eps_scale_set,
               sad_amp_max_set,
               sad_n_div_max_set,
               max_num_runge_kutta_step_set,
               spin_n0_direction_user_set_set,
               debug_set,
               ds_track_step_set,
               dt_track_step_set,
               cathode_strength_cutoff_set,
               sc_rel_tol_tracking_set,
               sc_abs_tol_tracking_set,
               beam_chamber_height_set,
               lsc_sigma_cutoff_set,
               particle_sigma_cutoff_set,
               space_charge_mesh_size_set,
               csr3d_mesh_size_set,
               n_bin_set,
               particle_bin_span_set,
               n_shield_images_set,
               sc_min_in_bin_set,
               lsc_kick_transverse_dependence_set,
               sc_debug_set,
               diagnostic_output_file_set,
               old_integrator_set,
               use_orientation_patches_set,
               print_info_messages_set,
               max_fringe_order_set,
               exact_model_set,
               exact_misalign_set,
               vertical_kick_set,
               cut_factor_set,
               translate_patch_drift_time_set
           );
         },
         nb::arg("ran_state") = nb::none(),
         nb::arg("ran_seed") = nb::none(),
         nb::arg("undeterministic_ran_function_called") = nb::none(),
         nb::arg("d_orb_set") = nb::none(),
         nb::arg("max_aperture_limit_set") = nb::none(),
         nb::arg("default_ds_step_set") = nb::none(),
         nb::arg("significant_length_set") = nb::none(),
         nb::arg("rel_tol_tracking_set") = nb::none(),
         nb::arg("abs_tol_tracking_set") = nb::none(),
         nb::arg("rel_tol_adaptive_tracking_set") = nb::none(),
         nb::arg("abs_tol_adaptive_tracking_set") = nb::none(),
         nb::arg("init_ds_adaptive_tracking_set") = nb::none(),
         nb::arg("min_ds_adaptive_tracking_set") = nb::none(),
         nb::arg("fatal_ds_adaptive_tracking_set") = nb::none(),
         nb::arg("synch_rad_scale_set") = nb::none(),
         nb::arg("autoscale_amp_abs_tol_set") = nb::none(),
         nb::arg("autoscale_amp_rel_tol_set") = nb::none(),
         nb::arg("autoscale_phase_tol_set") = nb::none(),
         nb::arg("rf_phase_below_transition_ref_set") = nb::none(),
         nb::arg("electric_dipole_moment_set") = nb::none(),
         nb::arg("taylor_order_set") = nb::none(),
         nb::arg("runge_kutta_order_set") = nb::none(),
         nb::arg("default_integ_order_set") = nb::none(),
         nb::arg("sr_wakes_on_set") = nb::none(),
         nb::arg("lr_wakes_on_set") = nb::none(),
         nb::arg("high_energy_space_charge_on_set") = nb::none(),
         nb::arg("csr_and_space_charge_on_set") = nb::none(),
         nb::arg("spin_tracking_on_set") = nb::none(),
         nb::arg("spin_sokolov_ternov_flipping_on_set") = nb::none(),
         nb::arg("radiation_damping_on_set") = nb::none(),
         nb::arg("radiation_zero_average_set") = nb::none(),
         nb::arg("radiation_fluctuations_on_set") = nb::none(),
         nb::arg("conserve_taylor_maps_set") = nb::none(),
         nb::arg("absolute_time_tracking_set") = nb::none(),
         nb::arg("absolute_time_ref_shift_set") = nb::none(),
         nb::arg("convert_to_kinetic_momentum_set") = nb::none(),
         nb::arg("aperture_limit_on_set") = nb::none(),
         nb::arg("normalize_twiss_set") = nb::none(),
         nb::arg("sad_eps_scale_set") = nb::none(),
         nb::arg("sad_amp_max_set") = nb::none(),
         nb::arg("sad_n_div_max_set") = nb::none(),
         nb::arg("max_num_runge_kutta_step_set") = nb::none(),
         nb::arg("spin_n0_direction_user_set_set") = nb::none(),
         nb::arg("debug_set") = nb::none(),
         nb::arg("ds_track_step_set") = nb::none(),
         nb::arg("dt_track_step_set") = nb::none(),
         nb::arg("cathode_strength_cutoff_set") = nb::none(),
         nb::arg("sc_rel_tol_tracking_set") = nb::none(),
         nb::arg("sc_abs_tol_tracking_set") = nb::none(),
         nb::arg("beam_chamber_height_set") = nb::none(),
         nb::arg("lsc_sigma_cutoff_set") = nb::none(),
         nb::arg("particle_sigma_cutoff_set") = nb::none(),
         nb::arg("space_charge_mesh_size_set") = nb::none(),
         nb::arg("csr3d_mesh_size_set") = nb::none(),
         nb::arg("n_bin_set") = nb::none(),
         nb::arg("particle_bin_span_set") = nb::none(),
         nb::arg("n_shield_images_set") = nb::none(),
         nb::arg("sc_min_in_bin_set") = nb::none(),
         nb::arg("lsc_kick_transverse_dependence_set") = nb::none(),
         nb::arg("sc_debug_set") = nb::none(),
         nb::arg("diagnostic_output_file_set") = nb::none(),
         nb::arg("old_integrator_set") = nb::none(),
         nb::arg("use_orientation_patches_set") = nb::none(),
         nb::arg("print_info_messages_set") = nb::none(),
         nb::arg("max_fringe_order_set") = nb::none(),
         nb::arg("exact_model_set") = nb::none(),
         nb::arg("exact_misalign_set") = nb::none(),
         nb::arg("vertical_kick_set") = nb::none(),
         nb::arg("cut_factor_set") = nb::none(),
         nb::arg("translate_patch_drift_time_set") = nb::none()
  )
      .def_prop_rw(
          "ran_state",
          &ExtraParsingInfoStruct::ran_state,
          &ExtraParsingInfoStruct::set_ran_state,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "ran_seed",
          &ExtraParsingInfoStruct::ran_seed,
          &ExtraParsingInfoStruct::set_ran_seed
      )
      .def_prop_rw(
          "undeterministic_ran_function_called",
          &ExtraParsingInfoStruct::undeterministic_ran_function_called,
          &ExtraParsingInfoStruct::set_undeterministic_ran_function_called,
          "Used with bmad_com"
      )
      .def_prop_rw(
          "d_orb_set",
          &ExtraParsingInfoStruct::d_orb_set,
          &ExtraParsingInfoStruct::set_d_orb_set
      )
      .def_prop_rw(
          "max_aperture_limit_set",
          &ExtraParsingInfoStruct::max_aperture_limit_set,
          &ExtraParsingInfoStruct::set_max_aperture_limit_set
      )
      .def_prop_rw(
          "default_ds_step_set",
          &ExtraParsingInfoStruct::default_ds_step_set,
          &ExtraParsingInfoStruct::set_default_ds_step_set
      )
      .def_prop_rw(
          "significant_length_set",
          &ExtraParsingInfoStruct::significant_length_set,
          &ExtraParsingInfoStruct::set_significant_length_set
      )
      .def_prop_rw(
          "rel_tol_tracking_set",
          &ExtraParsingInfoStruct::rel_tol_tracking_set,
          &ExtraParsingInfoStruct::set_rel_tol_tracking_set
      )
      .def_prop_rw(
          "abs_tol_tracking_set",
          &ExtraParsingInfoStruct::abs_tol_tracking_set,
          &ExtraParsingInfoStruct::set_abs_tol_tracking_set
      )
      .def_prop_rw(
          "rel_tol_adaptive_tracking_set",
          &ExtraParsingInfoStruct::rel_tol_adaptive_tracking_set,
          &ExtraParsingInfoStruct::set_rel_tol_adaptive_tracking_set
      )
      .def_prop_rw(
          "abs_tol_adaptive_tracking_set",
          &ExtraParsingInfoStruct::abs_tol_adaptive_tracking_set,
          &ExtraParsingInfoStruct::set_abs_tol_adaptive_tracking_set
      )
      .def_prop_rw(
          "init_ds_adaptive_tracking_set",
          &ExtraParsingInfoStruct::init_ds_adaptive_tracking_set,
          &ExtraParsingInfoStruct::set_init_ds_adaptive_tracking_set
      )
      .def_prop_rw(
          "min_ds_adaptive_tracking_set",
          &ExtraParsingInfoStruct::min_ds_adaptive_tracking_set,
          &ExtraParsingInfoStruct::set_min_ds_adaptive_tracking_set
      )
      .def_prop_rw(
          "fatal_ds_adaptive_tracking_set",
          &ExtraParsingInfoStruct::fatal_ds_adaptive_tracking_set,
          &ExtraParsingInfoStruct::set_fatal_ds_adaptive_tracking_set
      )
      .def_prop_rw(
          "synch_rad_scale_set",
          &ExtraParsingInfoStruct::synch_rad_scale_set,
          &ExtraParsingInfoStruct::set_synch_rad_scale_set
      )
      .def_prop_rw(
          "autoscale_amp_abs_tol_set",
          &ExtraParsingInfoStruct::autoscale_amp_abs_tol_set,
          &ExtraParsingInfoStruct::set_autoscale_amp_abs_tol_set
      )
      .def_prop_rw(
          "autoscale_amp_rel_tol_set",
          &ExtraParsingInfoStruct::autoscale_amp_rel_tol_set,
          &ExtraParsingInfoStruct::set_autoscale_amp_rel_tol_set
      )
      .def_prop_rw(
          "autoscale_phase_tol_set",
          &ExtraParsingInfoStruct::autoscale_phase_tol_set,
          &ExtraParsingInfoStruct::set_autoscale_phase_tol_set
      )
      .def_prop_rw(
          "rf_phase_below_transition_ref_set",
          &ExtraParsingInfoStruct::rf_phase_below_transition_ref_set,
          &ExtraParsingInfoStruct::set_rf_phase_below_transition_ref_set
      )
      .def_prop_rw(
          "electric_dipole_moment_set",
          &ExtraParsingInfoStruct::electric_dipole_moment_set,
          &ExtraParsingInfoStruct::set_electric_dipole_moment_set
      )
      .def_prop_rw(
          "taylor_order_set",
          &ExtraParsingInfoStruct::taylor_order_set,
          &ExtraParsingInfoStruct::set_taylor_order_set
      )
      .def_prop_rw(
          "runge_kutta_order_set",
          &ExtraParsingInfoStruct::runge_kutta_order_set,
          &ExtraParsingInfoStruct::set_runge_kutta_order_set
      )
      .def_prop_rw(
          "default_integ_order_set",
          &ExtraParsingInfoStruct::default_integ_order_set,
          &ExtraParsingInfoStruct::set_default_integ_order_set
      )
      .def_prop_rw(
          "sr_wakes_on_set",
          &ExtraParsingInfoStruct::sr_wakes_on_set,
          &ExtraParsingInfoStruct::set_sr_wakes_on_set
      )
      .def_prop_rw(
          "lr_wakes_on_set",
          &ExtraParsingInfoStruct::lr_wakes_on_set,
          &ExtraParsingInfoStruct::set_lr_wakes_on_set
      )
      .def_prop_rw(
          "high_energy_space_charge_on_set",
          &ExtraParsingInfoStruct::high_energy_space_charge_on_set,
          &ExtraParsingInfoStruct::set_high_energy_space_charge_on_set
      )
      .def_prop_rw(
          "csr_and_space_charge_on_set",
          &ExtraParsingInfoStruct::csr_and_space_charge_on_set,
          &ExtraParsingInfoStruct::set_csr_and_space_charge_on_set
      )
      .def_prop_rw(
          "spin_tracking_on_set",
          &ExtraParsingInfoStruct::spin_tracking_on_set,
          &ExtraParsingInfoStruct::set_spin_tracking_on_set
      )
      .def_prop_rw(
          "spin_sokolov_ternov_flipping_on_set",
          &ExtraParsingInfoStruct::spin_sokolov_ternov_flipping_on_set,
          &ExtraParsingInfoStruct::set_spin_sokolov_ternov_flipping_on_set
      )
      .def_prop_rw(
          "radiation_damping_on_set",
          &ExtraParsingInfoStruct::radiation_damping_on_set,
          &ExtraParsingInfoStruct::set_radiation_damping_on_set
      )
      .def_prop_rw(
          "radiation_zero_average_set",
          &ExtraParsingInfoStruct::radiation_zero_average_set,
          &ExtraParsingInfoStruct::set_radiation_zero_average_set
      )
      .def_prop_rw(
          "radiation_fluctuations_on_set",
          &ExtraParsingInfoStruct::radiation_fluctuations_on_set,
          &ExtraParsingInfoStruct::set_radiation_fluctuations_on_set
      )
      .def_prop_rw(
          "conserve_taylor_maps_set",
          &ExtraParsingInfoStruct::conserve_taylor_maps_set,
          &ExtraParsingInfoStruct::set_conserve_taylor_maps_set
      )
      .def_prop_rw(
          "absolute_time_tracking_set",
          &ExtraParsingInfoStruct::absolute_time_tracking_set,
          &ExtraParsingInfoStruct::set_absolute_time_tracking_set
      )
      .def_prop_rw(
          "absolute_time_ref_shift_set",
          &ExtraParsingInfoStruct::absolute_time_ref_shift_set,
          &ExtraParsingInfoStruct::set_absolute_time_ref_shift_set
      )
      .def_prop_rw(
          "convert_to_kinetic_momentum_set",
          &ExtraParsingInfoStruct::convert_to_kinetic_momentum_set,
          &ExtraParsingInfoStruct::set_convert_to_kinetic_momentum_set
      )
      .def_prop_rw(
          "aperture_limit_on_set",
          &ExtraParsingInfoStruct::aperture_limit_on_set,
          &ExtraParsingInfoStruct::set_aperture_limit_on_set
      )
      .def_prop_rw(
          "normalize_twiss_set",
          &ExtraParsingInfoStruct::normalize_twiss_set,
          &ExtraParsingInfoStruct::set_normalize_twiss_set
      )
      .def_prop_rw(
          "sad_eps_scale_set",
          &ExtraParsingInfoStruct::sad_eps_scale_set,
          &ExtraParsingInfoStruct::set_sad_eps_scale_set
      )
      .def_prop_rw(
          "sad_amp_max_set",
          &ExtraParsingInfoStruct::sad_amp_max_set,
          &ExtraParsingInfoStruct::set_sad_amp_max_set
      )
      .def_prop_rw(
          "sad_n_div_max_set",
          &ExtraParsingInfoStruct::sad_n_div_max_set,
          &ExtraParsingInfoStruct::set_sad_n_div_max_set
      )
      .def_prop_rw(
          "max_num_runge_kutta_step_set",
          &ExtraParsingInfoStruct::max_num_runge_kutta_step_set,
          &ExtraParsingInfoStruct::set_max_num_runge_kutta_step_set
      )
      .def_prop_rw(
          "spin_n0_direction_user_set_set",
          &ExtraParsingInfoStruct::spin_n0_direction_user_set_set,
          &ExtraParsingInfoStruct::set_spin_n0_direction_user_set_set
      )
      .def_prop_rw(
          "debug_set",
          &ExtraParsingInfoStruct::debug_set,
          &ExtraParsingInfoStruct::set_debug_set,
          "Used with space_charge_com"
      )
      .def_prop_rw(
          "ds_track_step_set",
          &ExtraParsingInfoStruct::ds_track_step_set,
          &ExtraParsingInfoStruct::set_ds_track_step_set
      )
      .def_prop_rw(
          "dt_track_step_set",
          &ExtraParsingInfoStruct::dt_track_step_set,
          &ExtraParsingInfoStruct::set_dt_track_step_set
      )
      .def_prop_rw(
          "cathode_strength_cutoff_set",
          &ExtraParsingInfoStruct::cathode_strength_cutoff_set,
          &ExtraParsingInfoStruct::set_cathode_strength_cutoff_set
      )
      .def_prop_rw(
          "sc_rel_tol_tracking_set",
          &ExtraParsingInfoStruct::sc_rel_tol_tracking_set,
          &ExtraParsingInfoStruct::set_sc_rel_tol_tracking_set,
          "For: space_charge_com%rel_tol_tracking"
      )
      .def_prop_rw(
          "sc_abs_tol_tracking_set",
          &ExtraParsingInfoStruct::sc_abs_tol_tracking_set,
          &ExtraParsingInfoStruct::set_sc_abs_tol_tracking_set,
          "For: space_charge_com%abs_tol_tracking"
      )
      .def_prop_rw(
          "beam_chamber_height_set",
          &ExtraParsingInfoStruct::beam_chamber_height_set,
          &ExtraParsingInfoStruct::set_beam_chamber_height_set
      )
      .def_prop_rw(
          "lsc_sigma_cutoff_set",
          &ExtraParsingInfoStruct::lsc_sigma_cutoff_set,
          &ExtraParsingInfoStruct::set_lsc_sigma_cutoff_set
      )
      .def_prop_rw(
          "particle_sigma_cutoff_set",
          &ExtraParsingInfoStruct::particle_sigma_cutoff_set,
          &ExtraParsingInfoStruct::set_particle_sigma_cutoff_set
      )
      .def_prop_rw(
          "space_charge_mesh_size_set",
          &ExtraParsingInfoStruct::space_charge_mesh_size_set,
          &ExtraParsingInfoStruct::set_space_charge_mesh_size_set
      )
      .def_prop_rw(
          "csr3d_mesh_size_set",
          &ExtraParsingInfoStruct::csr3d_mesh_size_set,
          &ExtraParsingInfoStruct::set_csr3d_mesh_size_set
      )
      .def_prop_rw(
          "n_bin_set",
          &ExtraParsingInfoStruct::n_bin_set,
          &ExtraParsingInfoStruct::set_n_bin_set
      )
      .def_prop_rw(
          "particle_bin_span_set",
          &ExtraParsingInfoStruct::particle_bin_span_set,
          &ExtraParsingInfoStruct::set_particle_bin_span_set
      )
      .def_prop_rw(
          "n_shield_images_set",
          &ExtraParsingInfoStruct::n_shield_images_set,
          &ExtraParsingInfoStruct::set_n_shield_images_set
      )
      .def_prop_rw(
          "sc_min_in_bin_set",
          &ExtraParsingInfoStruct::sc_min_in_bin_set,
          &ExtraParsingInfoStruct::set_sc_min_in_bin_set
      )
      .def_prop_rw(
          "lsc_kick_transverse_dependence_set",
          &ExtraParsingInfoStruct::lsc_kick_transverse_dependence_set,
          &ExtraParsingInfoStruct::set_lsc_kick_transverse_dependence_set
      )
      .def_prop_rw(
          "sc_debug_set",
          &ExtraParsingInfoStruct::sc_debug_set,
          &ExtraParsingInfoStruct::set_sc_debug_set
      )
      .def_prop_rw(
          "diagnostic_output_file_set",
          &ExtraParsingInfoStruct::diagnostic_output_file_set,
          &ExtraParsingInfoStruct::set_diagnostic_output_file_set,
          "Used with ptc_com"
      )
      .def_prop_rw(
          "old_integrator_set",
          &ExtraParsingInfoStruct::old_integrator_set,
          &ExtraParsingInfoStruct::set_old_integrator_set
      )
      .def_prop_rw(
          "use_orientation_patches_set",
          &ExtraParsingInfoStruct::use_orientation_patches_set,
          &ExtraParsingInfoStruct::set_use_orientation_patches_set
      )
      .def_prop_rw(
          "print_info_messages_set",
          &ExtraParsingInfoStruct::print_info_messages_set,
          &ExtraParsingInfoStruct::set_print_info_messages_set
      )
      .def_prop_rw(
          "max_fringe_order_set",
          &ExtraParsingInfoStruct::max_fringe_order_set,
          &ExtraParsingInfoStruct::set_max_fringe_order_set
      )
      .def_prop_rw(
          "exact_model_set",
          &ExtraParsingInfoStruct::exact_model_set,
          &ExtraParsingInfoStruct::set_exact_model_set
      )
      .def_prop_rw(
          "exact_misalign_set",
          &ExtraParsingInfoStruct::exact_misalign_set,
          &ExtraParsingInfoStruct::set_exact_misalign_set
      )
      .def_prop_rw(
          "vertical_kick_set",
          &ExtraParsingInfoStruct::vertical_kick_set,
          &ExtraParsingInfoStruct::set_vertical_kick_set
      )
      .def_prop_rw(
          "cut_factor_set",
          &ExtraParsingInfoStruct::cut_factor_set,
          &ExtraParsingInfoStruct::set_cut_factor_set
      )
      .def_prop_rw(
          "translate_patch_drift_time_set",
          &ExtraParsingInfoStruct::translate_patch_drift_time_set,
          &ExtraParsingInfoStruct::set_translate_patch_drift_time_set
      )

      .def("__repr__", [](const ExtraParsingInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ExtraParsingInfoStruct &self) {
            return ExtraParsingInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ExtraParsingInfoStruct &self, nb::dict &memo) {
            return ExtraParsingInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ExtraParsingInfoStruct &self, const ExtraParsingInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ExtraParsingInfoStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D ExtraParsingInfoStruct arrays are not used in structs/routines
  // 2D ExtraParsingInfoStruct arrays are not used in structs/routines
  // 3D ExtraParsingInfoStruct arrays are not used in structs/routines
}