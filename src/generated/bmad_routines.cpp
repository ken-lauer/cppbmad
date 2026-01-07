#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/bmad_routines.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"
#include "json.hpp"

using namespace Bmad;

using json = nlohmann::json;
Bmad::AbMultipoleKick Bmad::ab_multipole_kick(
    double a,
    double b,
    int n,
    int ref_species,
    int ele_orientation,
    CoordProxy& coord,
    std::optional<int> pole_type,
    std::optional<double> scale) {
  double _kx{};
  double _ky{};
  FixedArray2D<Real, 2, 2> dk;
  double _dk_vec[2 * 2];
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  double scale_lvalue;
  auto* _scale{&scale_lvalue};
  if (scale.has_value()) {
    scale_lvalue = scale.value();
  } else {
    _scale = nullptr;
  }
  fortran_ab_multipole_kick(
      /* double& */ a,
      /* double& */ b,
      /* int& */ n,
      /* int& */ ref_species,
      /* int& */ ele_orientation,
      /* void* */ coord.get_fortran_ptr(),
      /* double& */ _kx,
      /* double& */ _ky,
      /* double* */ _dk_vec,
      /* int* */ _pole_type,
      /* double* */ _scale);
  vec_to_matrix(_dk_vec, dk);
  return AbMultipoleKick{_kx, _ky, dk};
}
void Bmad::ab_multipole_kicks(
    RealAlloc1D& an,
    RealAlloc1D& bn,
    int ix_pole_max,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<int> pole_type,
    std::optional<double> scale,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  double scale_lvalue;
  auto* _scale{&scale_lvalue};
  if (scale.has_value()) {
    scale_lvalue = scale.value();
  } else {
    _scale = nullptr;
  }
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_ab_multipole_kicks(
      /* void* */ an.get_fortran_ptr(),
      /* void* */ bn.get_fortran_ptr(),
      /* int& */ ix_pole_max,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* double* */ _scale,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::absolute_photon_position(CoordProxy& e_orb, CoordProxy& photon_orb) {
  fortran_absolute_photon_position(
      /* void* */ e_orb.get_fortran_ptr(),
      /* void* */ photon_orb.get_fortran_ptr());
}
void Bmad::absolute_time_tracking(EleProxy& ele, bool& is_abs_time) {
  fortran_absolute_time_tracking(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ is_abs_time);
}
void Bmad::ac_kicker_amp(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<double> true_time,
    double& ac_amp) {
  double true_time_lvalue;
  auto* _true_time{&true_time_lvalue};
  if (true_time.has_value()) {
    true_time_lvalue = true_time.value();
  } else {
    _true_time = nullptr;
  }
  fortran_ac_kicker_amp(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _true_time,
      /* double& */ ac_amp);
}
Bmad::ActionToXyz Bmad::action_to_xyz(
    LatProxy& ring,
    int ix,
    FixedArray1D<Real, 6> J) {
  auto* _J = J.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 6> _X;
  bool _err_flag{};
  fortran_action_to_xyz(
      /* void* */ ring.get_fortran_ptr(),
      /* int& */ ix,
      /* double* */ _J,
      /* double* */ _X.data(),
      /* bool& */ _err_flag);
  return ActionToXyz{_X, _err_flag};
}
void Bmad::add_lattice_control_structs(
    EleProxy& ele,
    std::optional<int> n_add_slave,
    std::optional<int> n_add_lord,
    std::optional<int> n_add_slave_field,
    std::optional<int> n_add_lord_field,
    std::optional<bool> add_at_end) {
  int n_add_slave_lvalue;
  auto* _n_add_slave{&n_add_slave_lvalue};
  if (n_add_slave.has_value()) {
    n_add_slave_lvalue = n_add_slave.value();
  } else {
    _n_add_slave = nullptr;
  }
  int n_add_lord_lvalue;
  auto* _n_add_lord{&n_add_lord_lvalue};
  if (n_add_lord.has_value()) {
    n_add_lord_lvalue = n_add_lord.value();
  } else {
    _n_add_lord = nullptr;
  }
  int n_add_slave_field_lvalue;
  auto* _n_add_slave_field{&n_add_slave_field_lvalue};
  if (n_add_slave_field.has_value()) {
    n_add_slave_field_lvalue = n_add_slave_field.value();
  } else {
    _n_add_slave_field = nullptr;
  }
  int n_add_lord_field_lvalue;
  auto* _n_add_lord_field{&n_add_lord_field_lvalue};
  if (n_add_lord_field.has_value()) {
    n_add_lord_field_lvalue = n_add_lord_field.value();
  } else {
    _n_add_lord_field = nullptr;
  }
  bool add_at_end_lvalue;
  auto* _add_at_end{&add_at_end_lvalue};
  if (add_at_end.has_value()) {
    add_at_end_lvalue = add_at_end.value();
  } else {
    _add_at_end = nullptr;
  }
  fortran_add_lattice_control_structs(
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _n_add_slave,
      /* int* */ _n_add_lord,
      /* int* */ _n_add_slave_field,
      /* int* */ _n_add_lord_field,
      /* bool* */ _add_at_end);
}
Bmad::AddSuperimpose Bmad::add_superimpose(
    LatProxy& lat,
    EleProxy& super_ele_in,
    int ix_branch,
    std::optional<bool> save_null_drift,
    std::optional<bool> create_jumbo_slave,
    std::optional<int> ix_insert,
    std::optional<bool> mangle_slave_names,
    std::optional<bool> wrap) {
  bool _err_flag{};
  EleProxy _super_ele_out;
  bool save_null_drift_lvalue;
  auto* _save_null_drift{&save_null_drift_lvalue};
  if (save_null_drift.has_value()) {
    save_null_drift_lvalue = save_null_drift.value();
  } else {
    _save_null_drift = nullptr;
  }
  bool create_jumbo_slave_lvalue;
  auto* _create_jumbo_slave{&create_jumbo_slave_lvalue};
  if (create_jumbo_slave.has_value()) {
    create_jumbo_slave_lvalue = create_jumbo_slave.value();
  } else {
    _create_jumbo_slave = nullptr;
  }
  int ix_insert_lvalue;
  auto* _ix_insert{&ix_insert_lvalue};
  if (ix_insert.has_value()) {
    ix_insert_lvalue = ix_insert.value();
  } else {
    _ix_insert = nullptr;
  }
  bool mangle_slave_names_lvalue;
  auto* _mangle_slave_names{&mangle_slave_names_lvalue};
  if (mangle_slave_names.has_value()) {
    mangle_slave_names_lvalue = mangle_slave_names.value();
  } else {
    _mangle_slave_names = nullptr;
  }
  bool wrap_lvalue;
  auto* _wrap{&wrap_lvalue};
  if (wrap.has_value()) {
    wrap_lvalue = wrap.value();
  } else {
    _wrap = nullptr;
  }
  fortran_add_superimpose(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ super_ele_in.get_fortran_ptr(),
      /* int& */ ix_branch,
      /* bool& */ _err_flag,
      /* void* */ _super_ele_out.get_fortran_ptr(),
      /* bool* */ _save_null_drift,
      /* bool* */ _create_jumbo_slave,
      /* int* */ _ix_insert,
      /* bool* */ _mangle_slave_names,
      /* bool* */ _wrap);
  return AddSuperimpose{_err_flag, std::move(_super_ele_out)};
}
void Bmad::add_this_multipass(
    LatProxy& lat,
    LatEleLocProxyAlloc1D& m_slaves,
    optional_ref<EleProxy> lord_in) {
  // intent=inout allocatable type array
  auto* _lord_in = lord_in.has_value() ? lord_in->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  fortran_add_this_multipass(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ m_slaves.get_fortran_ptr(),
      /* void* */ _lord_in);
}
void Bmad::add_this_taylor_term(
    EleProxy& ele,
    int& i_out,
    double& coef,
    FixedArray1D<Int, 6> expn) {
  auto* _expn = expn.data(); // CppWrapperGeneralArgument
  fortran_add_this_taylor_term(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ i_out,
      /* double& */ coef,
      /* int* */ _expn);
}
void Bmad::adjust_super_slave_names(
    LatProxy& lat,
    int& ix1_lord,
    int& ix2_lord,
    optional_ref<bool> first_time) {
  auto* _first_time =
      first_time.has_value() ? &first_time->get() : nullptr; // inout, optional
  fortran_adjust_super_slave_names(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ ix1_lord,
      /* int& */ ix2_lord,
      /* bool* */ _first_time);
}
void Bmad::allocate_branch_array(LatProxy& lat, int upper_bound) {
  fortran_allocate_branch_array(
      /* void* */ lat.get_fortran_ptr(), /* int& */ upper_bound);
}
void Bmad::allocate_lat_ele_array(
    LatProxy& lat,
    std::optional<int> upper_bound,
    std::optional<int> ix_branch,
    std::optional<bool> do_ramper_slave_setup) {
  int upper_bound_lvalue;
  auto* _upper_bound{&upper_bound_lvalue};
  if (upper_bound.has_value()) {
    upper_bound_lvalue = upper_bound.value();
  } else {
    _upper_bound = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool do_ramper_slave_setup_lvalue;
  auto* _do_ramper_slave_setup{&do_ramper_slave_setup_lvalue};
  if (do_ramper_slave_setup.has_value()) {
    do_ramper_slave_setup_lvalue = do_ramper_slave_setup.value();
  } else {
    _do_ramper_slave_setup = nullptr;
  }
  fortran_allocate_lat_ele_array(
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _upper_bound,
      /* int* */ _ix_branch,
      /* bool* */ _do_ramper_slave_setup);
}
void Bmad::angle_between_polars(
    SpinPolarProxy& polar1,
    SpinPolarProxy& polar2,
    double& angle) {
  fortran_angle_between_polars(
      /* void* */ polar1.get_fortran_ptr(),
      /* void* */ polar2.get_fortran_ptr(),
      /* double& */ angle);
}
void Bmad::angle_to_canonical_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type) {
  const char* _coord_type =
      coord_type.has_value() ? coord_type->c_str() : nullptr;
  fortran_angle_to_canonical_coords(
      /* void* */ orbit.get_fortran_ptr(), /* const char* */ _coord_type);
}
void Bmad::aperture_bookkeeper(EleProxy& ele) {
  fortran_aperture_bookkeeper(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::apply_all_rampers(LatProxy& lat) {
  bool _err_flag{};
  fortran_apply_all_rampers(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::apply_energy_kick(
    double dE,
    CoordProxy& orbit,
    FixedArray1D<Real, 2> ddE_dr,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  auto* _ddE_dr = ddE_dr.data(); // CppWrapperGeneralArgument
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_apply_energy_kick(
      /* double& */ dE,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _ddE_dr,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::apply_patch_to_ptc_fibre(EleProxy& ele) {
  fortran_apply_patch_to_ptc_fibre(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::apply_rampers_to_slave(EleProxy& slave) {
  bool _err_flag{};
  fortran_apply_rampers_to_slave(
      /* void* */ slave.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::array_re_str(
    RealAlloc1D& arr,
    optional_ref<std::string> parens_in,
    std::string& str_out) {
  // intent=inout allocatable general array
  const char* _parens_in =
      parens_in.has_value() ? parens_in->get().c_str() : nullptr;
  auto _str_out = str_out.c_str(); // ptr, inout, required
  fortran_array_re_str(
      /* void* */ arr.get_fortran_ptr(),
      /* const char* */ _parens_in,
      /* const char* */ _str_out);
}
void Bmad::astra_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double& field_value) {
  fortran_astra_max_field_reference(
      /* void* */ pt0.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ field_value);
}
void Bmad::at_this_ele_end(int now_at, int where_at, bool& is_at_this_end) {
  fortran_at_this_ele_end(
      /* int& */ now_at, /* int& */ where_at, /* bool& */ is_at_this_end);
}
void Bmad::attribute_bookkeeper(
    EleProxy& ele,
    std::optional<bool> force_bookkeeping) {
  bool force_bookkeeping_lvalue;
  auto* _force_bookkeeping{&force_bookkeeping_lvalue};
  if (force_bookkeeping.has_value()) {
    force_bookkeeping_lvalue = force_bookkeeping.value();
  } else {
    _force_bookkeeping = nullptr;
  }
  fortran_attribute_bookkeeper(
      /* void* */ ele.get_fortran_ptr(), /* bool* */ _force_bookkeeping);
}
Bmad::AttributeFree1 Bmad::attribute_free(
    int ix_ele,
    std::string attrib_name,
    LatProxy& lat,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto* _except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto* _dependent_attribs_free{&dependent_attribs_free_lvalue};
  if (dependent_attribs_free.has_value()) {
    dependent_attribs_free_lvalue = dependent_attribs_free.value();
  } else {
    _dependent_attribs_free = nullptr;
  }
  int _why_not_free{};
  bool _free{};
  fortran_attribute_free1(
      /* int& */ ix_ele,
      /* const char* */ _attrib_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _err_print_flag,
      /* bool* */ _except_overlay,
      /* bool* */ _dependent_attribs_free,
      /* int& */ _why_not_free,
      /* bool& */ _free);
  return AttributeFree1{_why_not_free, _free};
}
Bmad::AttributeFree2 Bmad::attribute_free(
    EleProxy& ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto* _except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto* _dependent_attribs_free{&dependent_attribs_free_lvalue};
  if (dependent_attribs_free.has_value()) {
    dependent_attribs_free_lvalue = dependent_attribs_free.value();
  } else {
    _dependent_attribs_free = nullptr;
  }
  int _why_not_free{};
  bool _free{};
  fortran_attribute_free2(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_name,
      /* bool* */ _err_print_flag,
      /* bool* */ _except_overlay,
      /* bool* */ _dependent_attribs_free,
      /* int& */ _why_not_free,
      /* bool& */ _free);
  return AttributeFree2{_why_not_free, _free};
}
Bmad::AttributeFree3 Bmad::attribute_free(
    int ix_ele,
    int ix_branch,
    std::string attrib_name,
    LatProxy& lat,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto* _except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto* _dependent_attribs_free{&dependent_attribs_free_lvalue};
  if (dependent_attribs_free.has_value()) {
    dependent_attribs_free_lvalue = dependent_attribs_free.value();
  } else {
    _dependent_attribs_free = nullptr;
  }
  int _why_not_free{};
  bool _free{};
  fortran_attribute_free3(
      /* int& */ ix_ele,
      /* int& */ ix_branch,
      /* const char* */ _attrib_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _err_print_flag,
      /* bool* */ _except_overlay,
      /* bool* */ _dependent_attribs_free,
      /* int& */ _why_not_free,
      /* bool& */ _free);
  return AttributeFree3{_why_not_free, _free};
}
Bmad::AttributeIndex1 Bmad::attribute_index(
    EleProxy& ele,
    std::string name,
    std::optional<bool> can_abbreviate,
    std::optional<bool> print_error) {
  auto _name = name.c_str();
  char _full_name[4096];
  bool can_abbreviate_lvalue;
  auto* _can_abbreviate{&can_abbreviate_lvalue};
  if (can_abbreviate.has_value()) {
    can_abbreviate_lvalue = can_abbreviate.value();
  } else {
    _can_abbreviate = nullptr;
  }
  bool print_error_lvalue;
  auto* _print_error{&print_error_lvalue};
  if (print_error.has_value()) {
    print_error_lvalue = print_error.value();
  } else {
    _print_error = nullptr;
  }
  int _attrib_index{};
  fortran_attribute_index1(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _name,
      /* const char* */ _full_name,
      /* bool* */ _can_abbreviate,
      /* bool* */ _print_error,
      /* int& */ _attrib_index);
  return AttributeIndex1{_full_name, _attrib_index};
}
Bmad::AttributeIndex2 Bmad::attribute_index(
    int key,
    std::string name,
    std::optional<bool> can_abbreviate,
    std::optional<bool> print_error) {
  auto _name = name.c_str();
  char _full_name[4096];
  bool can_abbreviate_lvalue;
  auto* _can_abbreviate{&can_abbreviate_lvalue};
  if (can_abbreviate.has_value()) {
    can_abbreviate_lvalue = can_abbreviate.value();
  } else {
    _can_abbreviate = nullptr;
  }
  bool print_error_lvalue;
  auto* _print_error{&print_error_lvalue};
  if (print_error.has_value()) {
    print_error_lvalue = print_error.value();
  } else {
    _print_error = nullptr;
  }
  int _attrib_index{};
  fortran_attribute_index2(
      /* int& */ key,
      /* const char* */ _name,
      /* const char* */ _full_name,
      /* bool* */ _can_abbreviate,
      /* bool* */ _print_error,
      /* int& */ _attrib_index);
  return AttributeIndex2{_full_name, _attrib_index};
}
std::string Bmad::attribute_name(
    int key,
    int ix_att,
    std::optional<bool> show_private) {
  bool show_private_lvalue;
  auto* _show_private{&show_private_lvalue};
  if (show_private.has_value()) {
    show_private_lvalue = show_private.value();
  } else {
    _show_private = nullptr;
  }
  char _attrib_name[4096];
  fortran_attribute_name1(
      /* int& */ key,
      /* int& */ ix_att,
      /* bool* */ _show_private,
      /* const char* */ _attrib_name);
  return _attrib_name;
}
std::string Bmad::attribute_name(
    EleProxy& ele,
    int ix_att,
    std::optional<bool> show_private) {
  bool show_private_lvalue;
  auto* _show_private{&show_private_lvalue};
  if (show_private.has_value()) {
    show_private_lvalue = show_private.value();
  } else {
    _show_private = nullptr;
  }
  char _attrib_name[4096];
  fortran_attribute_name2(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ ix_att,
      /* bool* */ _show_private,
      /* const char* */ _attrib_name);
  return _attrib_name;
}
int Bmad::attribute_type(std::string attrib_name, optional_ref<EleProxy> ele) {
  auto _attrib_name = attrib_name.c_str();
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  int _attrib_type{};
  fortran_attribute_type(
      /* const char* */ _attrib_name,
      /* void* */ _ele,
      /* int& */ _attrib_type);
  return _attrib_type;
}
std::string Bmad::attribute_units(
    std::string attrib_name,
    std::optional<std::string> unrecognized_units) {
  auto _attrib_name = attrib_name.c_str();
  const char* _unrecognized_units =
      unrecognized_units.has_value() ? unrecognized_units->c_str() : nullptr;
  char _attrib_units[4096];
  fortran_attribute_units(
      /* const char* */ _attrib_name,
      /* const char* */ _unrecognized_units,
      /* const char* */ _attrib_units);
  return _attrib_units;
}
bool Bmad::autoscale_phase_and_amp(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> scale_phase,
    std::optional<bool> scale_amp,
    std::optional<bool> call_bookkeeper) {
  bool _err_flag{};
  bool scale_phase_lvalue;
  auto* _scale_phase{&scale_phase_lvalue};
  if (scale_phase.has_value()) {
    scale_phase_lvalue = scale_phase.value();
  } else {
    _scale_phase = nullptr;
  }
  bool scale_amp_lvalue;
  auto* _scale_amp{&scale_amp_lvalue};
  if (scale_amp.has_value()) {
    scale_amp_lvalue = scale_amp.value();
  } else {
    _scale_amp = nullptr;
  }
  bool call_bookkeeper_lvalue;
  auto* _call_bookkeeper{&call_bookkeeper_lvalue};
  if (call_bookkeeper.has_value()) {
    call_bookkeeper_lvalue = call_bookkeeper.value();
  } else {
    _call_bookkeeper = nullptr;
  }
  fortran_autoscale_phase_and_amp(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _scale_phase,
      /* bool* */ _scale_amp,
      /* bool* */ _call_bookkeeper);
  return _err_flag;
}
void Bmad::average_twiss(
    double frac1,
    TwissProxy& twiss1,
    TwissProxy& twiss2,
    TwissProxy& ave_twiss) {
  fortran_average_twiss(
      /* double& */ frac1,
      /* void* */ twiss1.get_fortran_ptr(),
      /* void* */ twiss2.get_fortran_ptr(),
      /* void* */ ave_twiss.get_fortran_ptr());
}
Bmad::BbiKick Bmad::bbi_kick(double x, double y, FixedArray1D<Real, 2> sigma) {
  auto* _sigma = sigma.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 2> _nk;
  FixedArray2D<Real, 2, 2> dnk;
  double _dnk_vec[2 * 2];
  fortran_bbi_kick(
      /* double& */ x,
      /* double& */ y,
      /* double* */ _sigma,
      /* double* */ _nk.data(),
      /* double* */ _dnk_vec);
  vec_to_matrix(_dnk_vec, dnk);
  return BbiKick{_nk, dnk};
}
RealAlloc1D Bmad::bbi_slice_calc(EleProxy& ele, int n_slice) {
  // intent=out allocatable general array
  auto z_slice{RealAlloc1D()};
  fortran_bbi_slice_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ n_slice,
      /* void* */ z_slice.get_fortran_ptr());
  return std::move(z_slice);
}
FixedArray2D<Real, 6, 6> Bmad::beam_envelope_ibs(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species) {
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  FixedArray2D<Real, 6, 6> ibs_mat;
  double _ibs_mat_vec[6 * 6];
  fortran_beam_envelope_ibs(
      /* double* */ _sigma_mat_vec,
      /* double* */ _ibs_mat_vec,
      /* bool& */ tail_cut,
      /* double& */ tau,
      /* double& */ energy,
      /* double& */ n_part,
      /* int& */ species);
  vec_to_matrix(_ibs_mat_vec, ibs_mat);
  return ibs_mat;
}
void Bmad::beam_equal_beam(BeamProxy& beam1, BeamProxy& beam2) {
  fortran_beam_equal_beam(
      /* void* */ beam1.get_fortran_ptr(), /* void* */ beam2.get_fortran_ptr());
}
void Bmad::beam_init_setup(
    BeamInitProxy& beam_init_in,
    EleProxy& ele,
    int species,
    optional_ref<NormalModesProxy> modes,
    std::optional<bool> err_flag,
    BeamInitProxy& beam_init_set) {
  auto* _modes = modes.has_value() ? modes->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  bool err_flag_lvalue;
  auto* _err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_beam_init_setup(
      /* void* */ beam_init_in.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ species,
      /* void* */ _modes,
      /* bool* */ _err_flag,
      /* void* */ beam_init_set.get_fortran_ptr());
}
Bmad::BeamTilts Bmad::beam_tilts(FixedArray2D<Real, 6, 6> S) {
  double _S_vec[6 * 6];
  matrix_to_vec(S, _S_vec);
  double _angle_xy{};
  double _angle_xz{};
  double _angle_yz{};
  double _angle_xpz{};
  double _angle_ypz{};
  fortran_beam_tilts(
      /* double* */ _S_vec,
      /* double& */ _angle_xy,
      /* double& */ _angle_xz,
      /* double& */ _angle_yz,
      /* double& */ _angle_xpz,
      /* double& */ _angle_ypz);
  return BeamTilts{_angle_xy, _angle_xz, _angle_yz, _angle_xpz, _angle_ypz};
}
void Bmad::bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<bool> track_spin) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  bool track_spin_lvalue;
  auto* _track_spin{&track_spin_lvalue};
  if (track_spin.has_value()) {
    track_spin_lvalue = track_spin.value();
  } else {
    _track_spin = nullptr;
  }
  fortran_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* bool* */ _track_spin);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
EmFieldProxy Bmad::bend_exact_multipole_field(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential) {
  EmFieldProxy _field;
  bool calc_dfield_lvalue;
  auto* _calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto* _calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_bend_exact_multipole_field(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ local_ref_frame,
      /* void* */ _field.get_fortran_ptr(),
      /* bool* */ _calc_dfield,
      /* bool* */ _calc_potential);
  return std::move(_field);
}
void Bmad::bend_length_has_been_set(EleProxy& ele, bool& is_set) {
  fortran_bend_length_has_been_set(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ is_set);
}
double Bmad::bend_photon_e_rel_init(std::optional<double> r_in) {
  double r_in_lvalue;
  auto* _r_in{&r_in_lvalue};
  if (r_in.has_value()) {
    r_in_lvalue = r_in.value();
  } else {
    _r_in = nullptr;
  }
  double _E_rel{};
  fortran_bend_photon_e_rel_init(/* double* */ _r_in, /* double& */ _E_rel);
  return _E_rel;
}
double Bmad::bend_photon_energy_integ_prob(
    double E_photon,
    double g_bend,
    double gamma) {
  double _integ_prob{};
  fortran_bend_photon_energy_integ_prob(
      /* double& */ E_photon,
      /* double& */ g_bend,
      /* double& */ gamma,
      /* double& */ _integ_prob);
  return _integ_prob;
}
double Bmad::bend_photon_energy_normalized_probability(double E_rel) {
  double _prob{};
  fortran_bend_photon_energy_normalized_probability(
      /* double& */ E_rel, /* double& */ _prob);
  return _prob;
}
CoordProxy Bmad::bend_photon_init(
    double g_bend_x,
    double g_bend_y,
    double gamma,
    std::optional<double> E_min,
    std::optional<double> E_max,
    std::optional<double> E_integ_prob,
    std::optional<double> vert_angle_min,
    std::optional<double> vert_angle_max,
    std::optional<bool> vert_angle_symmetric,
    std::optional<double> emit_probability) {
  CoordProxy _orbit;
  double E_min_lvalue;
  auto* _E_min{&E_min_lvalue};
  if (E_min.has_value()) {
    E_min_lvalue = E_min.value();
  } else {
    _E_min = nullptr;
  }
  double E_max_lvalue;
  auto* _E_max{&E_max_lvalue};
  if (E_max.has_value()) {
    E_max_lvalue = E_max.value();
  } else {
    _E_max = nullptr;
  }
  double E_integ_prob_lvalue;
  auto* _E_integ_prob{&E_integ_prob_lvalue};
  if (E_integ_prob.has_value()) {
    E_integ_prob_lvalue = E_integ_prob.value();
  } else {
    _E_integ_prob = nullptr;
  }
  double vert_angle_min_lvalue;
  auto* _vert_angle_min{&vert_angle_min_lvalue};
  if (vert_angle_min.has_value()) {
    vert_angle_min_lvalue = vert_angle_min.value();
  } else {
    _vert_angle_min = nullptr;
  }
  double vert_angle_max_lvalue;
  auto* _vert_angle_max{&vert_angle_max_lvalue};
  if (vert_angle_max.has_value()) {
    vert_angle_max_lvalue = vert_angle_max.value();
  } else {
    _vert_angle_max = nullptr;
  }
  bool vert_angle_symmetric_lvalue;
  auto* _vert_angle_symmetric{&vert_angle_symmetric_lvalue};
  if (vert_angle_symmetric.has_value()) {
    vert_angle_symmetric_lvalue = vert_angle_symmetric.value();
  } else {
    _vert_angle_symmetric = nullptr;
  }
  double emit_probability_lvalue;
  auto* _emit_probability{&emit_probability_lvalue};
  if (emit_probability.has_value()) {
    emit_probability_lvalue = emit_probability.value();
  } else {
    _emit_probability = nullptr;
  }
  fortran_bend_photon_init(
      /* double& */ g_bend_x,
      /* double& */ g_bend_y,
      /* double& */ gamma,
      /* void* */ _orbit.get_fortran_ptr(),
      /* double* */ _E_min,
      /* double* */ _E_max,
      /* double* */ _E_integ_prob,
      /* double* */ _vert_angle_min,
      /* double* */ _vert_angle_max,
      /* bool* */ _vert_angle_symmetric,
      /* double* */ _emit_probability);
  return std::move(_orbit);
}
CoordProxy Bmad::bend_photon_polarization_init(
    double g_bend_x,
    double g_bend_y,
    double E_rel,
    double gamma_phi) {
  CoordProxy _orbit;
  fortran_bend_photon_polarization_init(
      /* double& */ g_bend_x,
      /* double& */ g_bend_y,
      /* double& */ E_rel,
      /* double& */ gamma_phi,
      /* void* */ _orbit.get_fortran_ptr());
  return std::move(_orbit);
}
double Bmad::bend_photon_vert_angle_init(
    double E_rel,
    double gamma,
    std::optional<double> r_in,
    std::optional<bool> invert) {
  double r_in_lvalue;
  auto* _r_in{&r_in_lvalue};
  if (r_in.has_value()) {
    r_in_lvalue = r_in.value();
  } else {
    _r_in = nullptr;
  }
  bool invert_lvalue;
  auto* _invert{&invert_lvalue};
  if (invert.has_value()) {
    invert_lvalue = invert.value();
  } else {
    _invert = nullptr;
  }
  double _phi{};
  fortran_bend_photon_vert_angle_init(
      /* double& */ E_rel,
      /* double& */ gamma,
      /* double* */ _r_in,
      /* bool* */ _invert,
      /* double& */ _phi);
  return _phi;
}
FixedArray2D<Real, 3, 3> Bmad::bend_shift(
    FloorPositionProxy& position1,
    double g,
    double delta_s,
    std::optional<double> ref_tilt,
    FloorPositionProxy& position2) {
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  double ref_tilt_lvalue;
  auto* _ref_tilt{&ref_tilt_lvalue};
  if (ref_tilt.has_value()) {
    ref_tilt_lvalue = ref_tilt.value();
  } else {
    _ref_tilt = nullptr;
  }
  fortran_bend_shift(
      /* void* */ position1.get_fortran_ptr(),
      /* double& */ g,
      /* double& */ delta_s,
      /* double* */ _w_mat_vec,
      /* double* */ _ref_tilt,
      /* void* */ position2.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
double Bmad::bend_vert_angle_integ_prob(
    double vert_angle,
    double E_rel,
    double gamma) {
  double _integ_prob{};
  fortran_bend_vert_angle_integ_prob(
      /* double& */ vert_angle,
      /* double& */ E_rel,
      /* double& */ gamma,
      /* double& */ _integ_prob);
  return _integ_prob;
}
double Bmad::bl_via_vlassov(
    double current,
    double alpha,
    double Energy,
    double sigma_p,
    double Vrf,
    double omega,
    double U0,
    double circ,
    double R,
    double L) {
  double _sigma_z{};
  fortran_bl_via_vlassov(
      /* double& */ current,
      /* double& */ alpha,
      /* double& */ Energy,
      /* double& */ sigma_p,
      /* double& */ Vrf,
      /* double& */ omega,
      /* double& */ U0,
      /* double& */ circ,
      /* double& */ R,
      /* double& */ L,
      /* double& */ _sigma_z);
  return _sigma_z;
}
Bmad::BmadParser Bmad::bmad_parser(
    std::string lat_file,
    std::optional<bool> make_mats6,
    std::optional<std::string> use_line) {
  auto _lat_file = lat_file.c_str();
  LatProxy _lat;
  bool make_mats6_lvalue;
  auto* _make_mats6{&make_mats6_lvalue};
  if (make_mats6.has_value()) {
    make_mats6_lvalue = make_mats6.value();
  } else {
    _make_mats6 = nullptr;
  }
  bool _digested_read_ok{};
  const char* _use_line = use_line.has_value() ? use_line->c_str() : nullptr;
  bool _err_flag{};
  LatProxy _parse_lat;
  fortran_bmad_parser(
      /* const char* */ _lat_file,
      /* void* */ _lat.get_fortran_ptr(),
      /* bool* */ _make_mats6,
      /* bool& */ _digested_read_ok,
      /* const char* */ _use_line,
      /* bool& */ _err_flag,
      /* void* */ _parse_lat.get_fortran_ptr());
  return BmadParser{
      std::move(_lat), _digested_read_ok, _err_flag, std::move(_parse_lat)};
}
void Bmad::bmad_parser2(
    std::string lat_file,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> orbit,
    std::optional<bool> make_mats6,
    optional_ref<bool> err_flag,
    optional_ref<LatProxy> parse_lat) {
  auto _lat_file = lat_file.c_str();
  // intent=in allocatable type array
  auto* _orbit = orbit.has_value() ? orbit->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  bool make_mats6_lvalue;
  auto* _make_mats6{&make_mats6_lvalue};
  if (make_mats6.has_value()) {
    make_mats6_lvalue = make_mats6.value();
  } else {
    _make_mats6 = nullptr;
  }
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  auto* _parse_lat = parse_lat.has_value() ? parse_lat->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  fortran_bmad_parser2(
      /* const char* */ _lat_file,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _orbit,
      /* bool* */ _make_mats6,
      /* bool* */ _err_flag,
      /* void* */ _parse_lat);
}
void Bmad::bmad_patch_parameters_to_ptc(
    FixedArray1D<Real, 3> ang,
    FixedArray2D<Real, 3, 3> exi) {
  auto* _ang = ang.data(); // CppWrapperGeneralArgument
  double _exi_vec[3 * 3];
  matrix_to_vec(exi, _exi_vec);
  fortran_bmad_patch_parameters_to_ptc(
      /* double* */ _ang, /* double* */ _exi_vec);
  vec_to_matrix(_exi_vec, exi);
}
void Bmad::bp_set_ran_status() {
  fortran_bp_set_ran_status();
}
void Bmad::branch_equal_branch(BranchProxy& branch1, BranchProxy& branch2) {
  fortran_branch_equal_branch(
      /* void* */ branch1.get_fortran_ptr(),
      /* void* */ branch2.get_fortran_ptr());
}
void Bmad::branch_name(BranchProxy& branch, std::string& name) {
  auto _name = name.c_str(); // ptr, inout, required
  fortran_branch_name(
      /* void* */ branch.get_fortran_ptr(), /* const char* */ _name);
}
void Bmad::branch_to_ptc_m_u(BranchProxy& branch) {
  fortran_branch_to_ptc_m_u(/* void* */ branch.get_fortran_ptr());
}
void Bmad::bunch_equal_bunch(BunchProxy& bunch1, BunchProxy& bunch2) {
  fortran_bunch_equal_bunch(
      /* void* */ bunch1.get_fortran_ptr(),
      /* void* */ bunch2.get_fortran_ptr());
}
FixedArray2D<Real, 2, 2> Bmad::c_to_cbar(EleProxy& ele) {
  FixedArray2D<Real, 2, 2> cbar_mat;
  double _cbar_mat_vec[2 * 2];
  fortran_c_to_cbar(
      /* void* */ ele.get_fortran_ptr(), /* double* */ _cbar_mat_vec);
  vec_to_matrix(_cbar_mat_vec, cbar_mat);
  return cbar_mat;
}
void Bmad::calc_bunch_params(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    bool error,
    std::optional<bool> print_err,
    std::optional<FixedArray2D<Real, 6, 6>> n_mat,
    std::optional<bool> is_time_coords,
    optional_ref<EleProxy> ele) {
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  double _n_mat_vec[6 * 6];
  const double* _n_mat = nullptr;
  if (n_mat.has_value()) {
    matrix_to_vec(n_mat.value(), _n_mat_vec);
    _n_mat = _n_mat_vec;
  }
  bool is_time_coords_lvalue;
  auto* _is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_calc_bunch_params(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* bool& */ error,
      /* bool* */ _print_err,
      /* double* */ _n_mat_vec,
      /* bool* */ _is_time_coords,
      /* void* */ _ele);
}
void Bmad::calc_bunch_params_slice(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    int plane,
    double slice_center,
    double slice_spread,
    bool err,
    std::optional<bool> print_err,
    std::optional<bool> is_time_coords,
    optional_ref<EleProxy> ele) {
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool is_time_coords_lvalue;
  auto* _is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_calc_bunch_params_slice(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* int& */ plane,
      /* double& */ slice_center,
      /* double& */ slice_spread,
      /* bool& */ err,
      /* bool* */ _print_err,
      /* bool* */ _is_time_coords,
      /* void* */ _ele);
}
void Bmad::calc_bunch_params_z_slice(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    FixedArray1D<Real, 2> slice_bounds,
    bool err,
    std::optional<bool> print_err,
    std::optional<bool> is_time_coords,
    optional_ref<EleProxy> ele) {
  auto* _slice_bounds = slice_bounds.data(); // CppWrapperGeneralArgument
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool is_time_coords_lvalue;
  auto* _is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_calc_bunch_params_z_slice(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* double* */ _slice_bounds,
      /* bool& */ err,
      /* bool* */ _print_err,
      /* bool* */ _is_time_coords,
      /* void* */ _ele);
}
BunchParamsProxy Bmad::calc_bunch_sigma_matrix_etc(
    CoordProxyAlloc1D& particle,
    RealAlloc1D& charge,
    optional_ref<bool> is_time_coords,
    optional_ref<EleProxy> ele) {
  // intent=in allocatable type array
  // intent=in allocatable general array
  BunchParamsProxy _bunch_params;
  auto* _is_time_coords = is_time_coords.has_value()
      ? &is_time_coords->get()
      : nullptr; // inout, optional
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_calc_bunch_sigma_matrix_etc(
      /* void* */ particle.get_fortran_ptr(),
      /* void* */ charge.get_fortran_ptr(),
      /* void* */ _bunch_params.get_fortran_ptr(),
      /* bool* */ _is_time_coords,
      /* void* */ _ele);
  return std::move(_bunch_params);
}
Bmad::CalcEmittancesAndTwissFromSigmaMatrix Bmad::
    calc_emittances_and_twiss_from_sigma_matrix(
        FixedArray2D<Real, 6, 6> sigma_mat,
        std::optional<bool> print_err) {
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  BunchParamsProxy _bunch_params;
  bool _error{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  FixedArray2D<Real, 6, 6> n_mat;
  double _n_mat_vec[6 * 6];
  fortran_calc_emittances_and_twiss_from_sigma_matrix(
      /* double* */ _sigma_mat_vec,
      /* void* */ _bunch_params.get_fortran_ptr(),
      /* bool& */ _error,
      /* bool* */ _print_err,
      /* double* */ _n_mat_vec);
  vec_to_matrix(_n_mat_vec, n_mat);
  return CalcEmittancesAndTwissFromSigmaMatrix{
      std::move(_bunch_params), _error, n_mat};
}
BunchParamsProxy Bmad::calc_spin_params(BunchProxy& bunch) {
  BunchParamsProxy _bunch_params;
  fortran_calc_spin_params(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ _bunch_params.get_fortran_ptr());
  return std::move(_bunch_params);
}
EleProxy Bmad::calc_super_slave_key(
    EleProxy& lord1,
    EleProxy& lord2,
    std::optional<bool> create_jumbo_slave) {
  EleProxy _slave;
  bool create_jumbo_slave_lvalue;
  auto* _create_jumbo_slave{&create_jumbo_slave_lvalue};
  if (create_jumbo_slave.has_value()) {
    create_jumbo_slave_lvalue = create_jumbo_slave.value();
  } else {
    _create_jumbo_slave = nullptr;
  }
  fortran_calc_super_slave_key(
      /* void* */ lord1.get_fortran_ptr(),
      /* void* */ lord2.get_fortran_ptr(),
      /* void* */ _slave.get_fortran_ptr(),
      /* bool* */ _create_jumbo_slave);
  return std::move(_slave);
}
Bmad::CalcWallRadius Bmad::calc_wall_radius(
    Wall3dVertexProxyAlloc1D& v,
    double cos_ang,
    double sin_ang) {
  // intent=in allocatable type array
  double _r_wall{};
  double _dr_dtheta{};
  int _ix_vertex{};
  fortran_calc_wall_radius(
      /* void* */ v.get_fortran_ptr(),
      /* double& */ cos_ang,
      /* double& */ sin_ang,
      /* double& */ _r_wall,
      /* double& */ _dr_dtheta,
      /* int& */ _ix_vertex);
  return CalcWallRadius{_r_wall, _dr_dtheta, _ix_vertex};
}
void Bmad::calc_z_tune(BranchProxy& branch) {
  fortran_calc_z_tune(/* void* */ branch.get_fortran_ptr());
}
void Bmad::canonical_to_angle_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type) {
  const char* _coord_type =
      coord_type.has_value() ? coord_type->c_str() : nullptr;
  fortran_canonical_to_angle_coords(
      /* void* */ orbit.get_fortran_ptr(), /* const char* */ _coord_type);
}
FixedArray2D<Real, 2, 2> Bmad::cbar_to_c(
    FixedArray2D<Real, 2, 2> cbar_mat,
    TwissProxy& a,
    TwissProxy& b) {
  double _cbar_mat_vec[2 * 2];
  matrix_to_vec(cbar_mat, _cbar_mat_vec);
  FixedArray2D<Real, 2, 2> c_mat;
  double _c_mat_vec[2 * 2];
  fortran_cbar_to_c(
      /* double* */ _cbar_mat_vec,
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* double* */ _c_mat_vec);
  vec_to_matrix(_c_mat_vec, c_mat);
  return c_mat;
}
void Bmad::check_aperture_limit(
    CoordProxy& orb,
    EleProxy& ele,
    int particle_at,
    LatParamProxy& param,
    optional_ref<CoordProxy> old_orb,
    std::optional<bool> check_momentum) {
  auto* _old_orb = old_orb.has_value() ? old_orb->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  bool check_momentum_lvalue;
  auto* _check_momentum{&check_momentum_lvalue};
  if (check_momentum.has_value()) {
    check_momentum_lvalue = check_momentum.value();
  } else {
    _check_momentum = nullptr;
  }
  fortran_check_aperture_limit(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _old_orb,
      /* bool* */ _check_momentum);
}
bool Bmad::check_controller_controls(
    int ele_key,
    ControlProxyAlloc1D& contrl,
    std::string name) {
  // intent=in allocatable type array
  auto _name = name.c_str();
  bool _err{};
  fortran_check_controller_controls(
      /* int& */ ele_key,
      /* void* */ contrl.get_fortran_ptr(),
      /* const char* */ _name,
      /* bool& */ _err);
  return _err;
}
void Bmad::check_for_superimpose_problem(
    BranchProxy& branch,
    EleProxy& super_ele,
    bool& err_flag,
    optional_ref<EleProxy> ref_ele,
    bool& wrap) {
  auto* _ref_ele = ref_ele.has_value() ? ref_ele->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  fortran_check_for_superimpose_problem(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ super_ele.get_fortran_ptr(),
      /* bool& */ err_flag,
      /* void* */ _ref_ele,
      /* bool& */ wrap);
}
Bmad::CheckIfSInBounds Bmad::check_if_s_in_bounds(
    BranchProxy& branch,
    double s,
    std::optional<bool> print_err) {
  bool _err_flag{};
  double _translated_s{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_check_if_s_in_bounds(
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ s,
      /* bool& */ _err_flag,
      /* double& */ _translated_s,
      /* bool* */ _print_err);
  return CheckIfSInBounds{_err_flag, _translated_s};
}
Bmad::ChooseQuadsForSetTune Bmad::choose_quads_for_set_tune(
    BranchProxy& branch,
    std::optional<std::string> mask) {
  // intent=out allocatable general array
  auto dk1{RealAlloc1D()};
  // intent=out allocatable type array
  auto eles{ElePointerProxyAlloc1D()};
  const char* _mask = mask.has_value() ? mask->c_str() : nullptr;
  bool _err_flag{};
  fortran_choose_quads_for_set_tune(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ dk1.get_fortran_ptr(),
      /* void* */ eles.get_fortran_ptr(),
      /* const char* */ _mask,
      /* bool& */ _err_flag);
  return ChooseQuadsForSetTune{std::move(dk1), std::move(eles), _err_flag};
}
Bmad::ChromCalc Bmad::chrom_calc(
    LatProxy& lat,
    double& delta_e,
    std::optional<double> pz,
    std::optional<int> ix_branch,
    optional_ref<CoordProxy> orb0) {
  double _chrom_a{};
  double _chrom_b{};
  bool _err_flag{};
  double pz_lvalue;
  auto* _pz{&pz_lvalue};
  if (pz.has_value()) {
    pz_lvalue = pz.value();
  } else {
    _pz = nullptr;
  }
  LatProxy _low_E_lat;
  LatProxy _high_E_lat;
  // intent=out allocatable type array
  auto low_E_orb{CoordProxyAlloc1D()};
  // intent=out allocatable type array
  auto high_E_orb{CoordProxyAlloc1D()};
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  auto* _orb0 = orb0.has_value() ? orb0->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  fortran_chrom_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ delta_e,
      /* double& */ _chrom_a,
      /* double& */ _chrom_b,
      /* bool& */ _err_flag,
      /* double* */ _pz,
      /* void* */ _low_E_lat.get_fortran_ptr(),
      /* void* */ _high_E_lat.get_fortran_ptr(),
      /* void* */ low_E_orb.get_fortran_ptr(),
      /* void* */ high_E_orb.get_fortran_ptr(),
      /* int* */ _ix_branch,
      /* void* */ _orb0);
  return ChromCalc{
      _chrom_a,
      _chrom_b,
      _err_flag,
      std::move(_low_E_lat),
      std::move(_high_E_lat),
      std::move(low_E_orb),
      std::move(high_E_orb)};
}
bool Bmad::chrom_tune(
    LatProxy& lat,
    double& delta_e,
    double target_x,
    double target_y,
    double err_tol) {
  bool _err_flag{};
  fortran_chrom_tune(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ delta_e,
      /* double& */ target_x,
      /* double& */ target_y,
      /* double& */ err_tol,
      /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::classical_radius(int species, double& radius) {
  fortran_classical_radius(/* int& */ species, /* double& */ radius);
}
LatProxy Bmad::clear_lat_1turn_mats() {
  LatProxy _lat;
  fortran_clear_lat_1turn_mats(/* void* */ _lat.get_fortran_ptr());
  return std::move(_lat);
}
void Bmad::clear_taylor_maps_from_elements(LatProxy& lat) {
  fortran_clear_taylor_maps_from_elements(/* void* */ lat.get_fortran_ptr());
}
bool Bmad::closed_orbit_calc(
    LatProxy& lat,
    CoordProxyAlloc1D& closed_orb,
    std::optional<int> i_dim,
    std::optional<int> direction,
    std::optional<int> ix_branch,
    std::optional<bool> print_err) {
  // intent=inout allocatable type array
  int i_dim_lvalue;
  auto* _i_dim{&i_dim_lvalue};
  if (i_dim.has_value()) {
    i_dim_lvalue = i_dim.value();
  } else {
    _i_dim = nullptr;
  }
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_closed_orbit_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ closed_orb.get_fortran_ptr(),
      /* int* */ _i_dim,
      /* int* */ _direction,
      /* int* */ _ix_branch,
      /* bool& */ _err_flag,
      /* bool* */ _print_err);
  return _err_flag;
}
Bmad::ClosedOrbitFromTracking Bmad::closed_orbit_from_tracking(
    LatProxy& lat,
    int i_dim,
    optional_ref<RealAlloc1D> eps_rel,
    optional_ref<RealAlloc1D> eps_abs,
    optional_ref<CoordProxy> init_guess) {
  // intent=out allocatable type array
  auto closed_orb{CoordProxyAlloc1D()};
  // intent=in allocatable general array
  auto* _eps_rel = eps_rel.has_value() ? eps_rel->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  // intent=in allocatable general array
  auto* _eps_abs = eps_abs.has_value() ? eps_abs->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  auto* _init_guess = init_guess.has_value()
      ? init_guess->get().get_fortran_ptr()
      : nullptr; // input, optional
  bool _err_flag{};
  fortran_closed_orbit_from_tracking(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ closed_orb.get_fortran_ptr(),
      /* int& */ i_dim,
      /* void* */ _eps_rel,
      /* void* */ _eps_abs,
      /* void* */ _init_guess,
      /* bool& */ _err_flag);
  return ClosedOrbitFromTracking{std::move(closed_orb), _err_flag};
}
void Bmad::cmplx_re_str(std::complex<double>& cmp, std::string& str_out) {
  auto _str_out = str_out.c_str(); // ptr, inout, required
  fortran_cmplx_re_str(
      /* std::complex<double>& */ cmp, /* const char* */ _str_out);
}
bool Bmad::combine_consecutive_elements(LatProxy& lat) {
  bool _error{};
  fortran_combine_consecutive_elements(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ _error);
  return _error;
}
void Bmad::complex_taylor_clean(ComplexTaylorProxy& complex_taylor) {
  fortran_complex_taylor_clean(/* void* */ complex_taylor.get_fortran_ptr());
}
void Bmad::complex_taylor_coef(
    ComplexTaylorProxy& complex_taylor,
    IntAlloc1D& exp,
    std::complex<double>& coef) {
  // intent=in allocatable general array
  fortran_complex_taylor_coef1(
      /* void* */ complex_taylor.get_fortran_ptr(),
      /* void* */ exp.get_fortran_ptr(),
      /* std::complex<double>& */ coef);
}
void Bmad::complex_taylor_coef(
    ComplexTaylorProxy& complex_taylor,
    std::optional<int> i1,
    std::optional<int> i2,
    std::optional<int> i3,
    std::optional<int> i4,
    std::optional<int> i5,
    std::optional<int> i6,
    std::optional<int> i7,
    std::optional<int> i8,
    std::optional<int> i9,
    std::complex<double>& coef) {
  int i1_lvalue;
  auto* _i1{&i1_lvalue};
  if (i1.has_value()) {
    i1_lvalue = i1.value();
  } else {
    _i1 = nullptr;
  }
  int i2_lvalue;
  auto* _i2{&i2_lvalue};
  if (i2.has_value()) {
    i2_lvalue = i2.value();
  } else {
    _i2 = nullptr;
  }
  int i3_lvalue;
  auto* _i3{&i3_lvalue};
  if (i3.has_value()) {
    i3_lvalue = i3.value();
  } else {
    _i3 = nullptr;
  }
  int i4_lvalue;
  auto* _i4{&i4_lvalue};
  if (i4.has_value()) {
    i4_lvalue = i4.value();
  } else {
    _i4 = nullptr;
  }
  int i5_lvalue;
  auto* _i5{&i5_lvalue};
  if (i5.has_value()) {
    i5_lvalue = i5.value();
  } else {
    _i5 = nullptr;
  }
  int i6_lvalue;
  auto* _i6{&i6_lvalue};
  if (i6.has_value()) {
    i6_lvalue = i6.value();
  } else {
    _i6 = nullptr;
  }
  int i7_lvalue;
  auto* _i7{&i7_lvalue};
  if (i7.has_value()) {
    i7_lvalue = i7.value();
  } else {
    _i7 = nullptr;
  }
  int i8_lvalue;
  auto* _i8{&i8_lvalue};
  if (i8.has_value()) {
    i8_lvalue = i8.value();
  } else {
    _i8 = nullptr;
  }
  int i9_lvalue;
  auto* _i9{&i9_lvalue};
  if (i9.has_value()) {
    i9_lvalue = i9.value();
  } else {
    _i9 = nullptr;
  }
  fortran_complex_taylor_coef2(
      /* void* */ complex_taylor.get_fortran_ptr(),
      /* int* */ _i1,
      /* int* */ _i2,
      /* int* */ _i3,
      /* int* */ _i4,
      /* int* */ _i5,
      /* int* */ _i6,
      /* int* */ _i7,
      /* int* */ _i8,
      /* int* */ _i9,
      /* std::complex<double>& */ coef);
}
void Bmad::complex_taylor_equal_complex_taylor(
    ComplexTaylorProxy& complex_taylor1,
    ComplexTaylorProxy& complex_taylor2) {
  fortran_complex_taylor_equal_complex_taylor(
      /* void* */ complex_taylor1.get_fortran_ptr(),
      /* void* */ complex_taylor2.get_fortran_ptr());
}
int Bmad::complex_taylor_exponent_index(FixedArray1D<Int, 6> expn) {
  auto* _expn = expn.data(); // CppWrapperGeneralArgument
  int _index{};
  fortran_complex_taylor_exponent_index(/* int* */ _expn, /* int& */ _index);
  return _index;
}
ComplexTaylorProxyAlloc1D Bmad::complex_taylor_make_unit() {
  // intent=out allocatable type array
  auto complex_taylor{ComplexTaylorProxyAlloc1D()};
  fortran_complex_taylor_make_unit(
      /* void* */ complex_taylor.get_fortran_ptr());
  return std::move(complex_taylor);
}
Bmad::ComplexTaylorToMat6 Bmad::complex_taylor_to_mat6(
    FixedArray1D<ComplexTaylorProxy, 6> a_complex_taylor,
    ComplexAlloc1D& r_in) {
  // intent=in allocatable general array
  FixedArray1D<Complex, 6> _vec0;
  FixedArray2D<Complex, 6, 6> mat6;
  std::complex<double> _mat6_vec[6 * 6];
  // intent=out allocatable general array
  auto r_out{ComplexAlloc1D()};
  fortran_complex_taylor_to_mat6(
      /* void* */ a_complex_taylor.data(),
      /* void* */ r_in.get_fortran_ptr(),
      /* std::complex<double>* */ _vec0.data(),
      /* std::complex<double>* */ _mat6_vec,
      /* void* */ r_out.get_fortran_ptr());
  vec_to_matrix(_mat6_vec, mat6);
  return ComplexTaylorToMat6{_vec0, mat6, std::move(r_out)};
}
void Bmad::complex_taylors_equal_complex_taylors(
    ComplexTaylorProxyAlloc1D& complex_taylor1,
    ComplexTaylorProxyAlloc1D& complex_taylor2) {
  // intent=inout allocatable type array
  // intent=in allocatable type array
  fortran_complex_taylors_equal_complex_taylors(
      /* void* */ complex_taylor1.get_fortran_ptr(),
      /* void* */ complex_taylor2.get_fortran_ptr());
}
void Bmad::compute_slave_coupler(EleProxy& slave) {
  fortran_compute_slave_coupler(/* void* */ slave.get_fortran_ptr());
}
void Bmad::concat_ele_taylor(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    bool err_flag,
    optional_ref<TaylorProxyAlloc1D> spin_taylor) {
  // intent=in allocatable type array
  // intent=in allocatable type array
  auto* _spin_taylor = spin_taylor.has_value()
      ? spin_taylor->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_concat_ele_taylor(
      /* void* */ orb_taylor.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ err_flag,
      /* void* */ _spin_taylor);
}
void Bmad::concat_taylor(
    TaylorProxyAlloc1D& taylor1,
    TaylorProxyAlloc1D& taylor2,
    TaylorProxyAlloc1D& taylor3) {
  // intent=in allocatable type array
  // intent=in allocatable type array
  // intent=in allocatable type array
  fortran_concat_taylor(
      /* void* */ taylor1.get_fortran_ptr(),
      /* void* */ taylor2.get_fortran_ptr(),
      /* void* */ taylor3.get_fortran_ptr());
}
FixedArray2D<Real, 6, 6> Bmad::concat_transfer_mat(
    FixedArray2D<Real, 6, 6> mat_1,
    FixedArray1D<Real, 6> vec_1,
    FixedArray2D<Real, 6, 6> mat_0,
    FixedArray1D<Real, 6> vec_0,
    FixedArray1D<Real, 6> vec_out) {
  double _mat_1_vec[6 * 6];
  matrix_to_vec(mat_1, _mat_1_vec);
  auto* _vec_1 = vec_1.data(); // CppWrapperGeneralArgument
  double _mat_0_vec[6 * 6];
  matrix_to_vec(mat_0, _mat_0_vec);
  auto* _vec_0 = vec_0.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 6, 6> mat_out;
  double _mat_out_vec[6 * 6];
  auto* _vec_out = vec_out.data(); // CppWrapperGeneralArgument
  fortran_concat_transfer_mat(
      /* double* */ _mat_1_vec,
      /* double* */ _vec_1,
      /* double* */ _mat_0_vec,
      /* double* */ _vec_0,
      /* double* */ _mat_out_vec,
      /* double* */ _vec_out);
  vec_to_matrix(_mat_out_vec, mat_out);
  return mat_out;
}
void Bmad::control_bookkeeper(
    LatProxy& lat,
    optional_ref<EleProxy> ele,
    std::optional<bool> err_flag) {
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool err_flag_lvalue;
  auto* _err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_control_bookkeeper(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _ele,
      /* bool* */ _err_flag);
}
void Bmad::convert_bend_exact_multipole(
    double g,
    int out_type,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> an,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> bn) {
  auto* _an = an.data(); // CppWrapperGeneralArgument
  auto* _bn = bn.data(); // CppWrapperGeneralArgument
  fortran_convert_bend_exact_multipole(
      /* double& */ g,
      /* int& */ out_type,
      /* double* */ _an,
      /* double* */ _bn);
}
Bmad::ConvertCoords Bmad::convert_coords(
    std::string in_type_str,
    CoordProxy& coord_in,
    EleProxy& ele) {
  auto _in_type_str = in_type_str.c_str();
  char _out_type_str[4096];
  CoordProxy _coord_out;
  bool _err_flag{};
  fortran_convert_coords(
      /* const char* */ _in_type_str,
      /* void* */ coord_in.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _out_type_str,
      /* void* */ _coord_out.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return ConvertCoords{_out_type_str, std::move(_coord_out), _err_flag};
}
EmFieldProxy Bmad::convert_field_ele_to_lab(
    EleProxy& ele,
    double s_here,
    bool forward_transform,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential) {
  EmFieldProxy _field;
  bool calc_dfield_lvalue;
  auto* _calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto* _calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_convert_field_ele_to_lab(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ s_here,
      /* bool& */ forward_transform,
      /* void* */ _field.get_fortran_ptr(),
      /* bool* */ _calc_dfield,
      /* bool* */ _calc_potential);
  return std::move(_field);
}
void Bmad::convert_local_cartesian_to_local_curvilinear(
    double& x,
    double& z,
    double& g,
    double& xout,
    double& sout) {
  fortran_convert_local_cartesian_to_local_curvilinear(
      /* double& */ x,
      /* double& */ z,
      /* double& */ g,
      /* double& */ xout,
      /* double& */ sout);
}
void Bmad::convert_local_curvilinear_to_local_cartesian(
    double& x,
    double& s,
    double& g,
    double& xout,
    double& zout) {
  fortran_convert_local_curvilinear_to_local_cartesian(
      /* double& */ x,
      /* double& */ s,
      /* double& */ g,
      /* double& */ xout,
      /* double& */ zout);
}
void Bmad::convert_particle_coordinates_s_to_t(
    CoordProxy& particle,
    double s_body,
    int orientation) {
  fortran_convert_particle_coordinates_s_to_t(
      /* void* */ particle.get_fortran_ptr(),
      /* double& */ s_body,
      /* int& */ orientation);
}
double Bmad::convert_particle_coordinates_t_to_s(
    CoordProxy& particle,
    EleProxy& ele,
    std::optional<bool> use_downstream_p0c) {
  double _s_body{};
  bool use_downstream_p0c_lvalue;
  auto* _use_downstream_p0c{&use_downstream_p0c_lvalue};
  if (use_downstream_p0c.has_value()) {
    use_downstream_p0c_lvalue = use_downstream_p0c.value();
  } else {
    _use_downstream_p0c = nullptr;
  }
  fortran_convert_particle_coordinates_t_to_s(
      /* void* */ particle.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _s_body,
      /* bool* */ _use_downstream_p0c);
  return _s_body;
}
Bmad::ConvertPcTo Bmad::convert_pc_to(double pc, int particle) {
  double _E_tot{};
  double _gamma{};
  double _kinetic{};
  double _beta{};
  double _brho{};
  double _beta1{};
  bool _err_flag{};
  fortran_convert_pc_to(
      /* double& */ pc,
      /* int& */ particle,
      /* double& */ _E_tot,
      /* double& */ _gamma,
      /* double& */ _kinetic,
      /* double& */ _beta,
      /* double& */ _brho,
      /* double& */ _beta1,
      /* bool& */ _err_flag);
  return ConvertPcTo{_E_tot, _gamma, _kinetic, _beta, _brho, _beta1, _err_flag};
}
Bmad::ConvertTotalEnergyTo Bmad::convert_total_energy_to(
    double E_tot,
    int particle,
    std::optional<bool> print_err) {
  double _gamma{};
  double _kinetic{};
  double _beta{};
  double _pc{};
  double _brho{};
  double _beta1{};
  bool _err_flag{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_convert_total_energy_to(
      /* double& */ E_tot,
      /* int& */ particle,
      /* double& */ _gamma,
      /* double& */ _kinetic,
      /* double& */ _beta,
      /* double& */ _pc,
      /* double& */ _brho,
      /* double& */ _beta1,
      /* bool& */ _err_flag,
      /* bool* */ _print_err);
  return ConvertTotalEnergyTo{
      _gamma, _kinetic, _beta, _pc, _brho, _beta1, _err_flag};
}
Bmad::ConverterDistributionParser Bmad::converter_distribution_parser(
    EleProxy& ele) {
  char _delim[4096];
  bool _delim_found{};
  bool _err_flag{};
  fortran_converter_distribution_parser(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* bool& */ _err_flag);
  return ConverterDistributionParser{_delim, _delim_found, _err_flag};
}
CoordProxy Bmad::coord_equal_coord(CoordProxy& coord2) {
  CoordProxy _coord1;
  fortran_coord_equal_coord(
      /* void* */ _coord1.get_fortran_ptr(),
      /* void* */ coord2.get_fortran_ptr());
  return std::move(_coord1);
}
std::string Bmad::coord_state_name(
    int coord_state,
    optional_ref<bool> one_word) {
  auto* _one_word =
      one_word.has_value() ? &one_word->get() : nullptr; // inout, optional
  char _state_str[4096];
  fortran_coord_state_name(
      /* int& */ coord_state,
      /* bool* */ _one_word,
      /* const char* */ _state_str);
  return _state_str;
}
void Bmad::coords_body_to_local(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& local_position) {
  double _w_mat_vec[3 * 3];
  const double* _w_mat = nullptr;
  if (w_mat.has_value()) {
    matrix_to_vec(w_mat.value(), _w_mat_vec);
    _w_mat = _w_mat_vec;
  }
  bool calculate_angles_lvalue;
  auto* _calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  fortran_coords_body_to_local(
      /* void* */ body_position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _w_mat_vec,
      /* bool* */ _calculate_angles,
      /* void* */ local_position.get_fortran_ptr());
}
void Bmad::coords_body_to_rel_exit(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& rel_exit) {
  double _w_mat_vec[3 * 3];
  const double* _w_mat = nullptr;
  if (w_mat.has_value()) {
    matrix_to_vec(w_mat.value(), _w_mat_vec);
    _w_mat = _w_mat_vec;
  }
  bool calculate_angles_lvalue;
  auto* _calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  fortran_coords_body_to_rel_exit(
      /* void* */ body_position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _w_mat_vec,
      /* bool* */ _calculate_angles,
      /* void* */ rel_exit.get_fortran_ptr());
}
bool Bmad::coords_curvilinear_to_floor(
    FixedArray1D<Real, 3> xys,
    BranchProxy& branch,
    FloorPositionProxy& global) {
  auto* _xys = xys.data(); // CppWrapperGeneralArgument
  bool _err_flag{};
  fortran_coords_curvilinear_to_floor(
      /* double* */ _xys,
      /* void* */ branch.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ global.get_fortran_ptr());
  return _err_flag;
}
Bmad::CoordsFloorToCurvilinear Bmad::coords_floor_to_curvilinear(
    FloorPositionProxy& floor_coords,
    EleProxy& ele0,
    FloorPositionProxy& local_coords) {
  EleProxy _ele1;
  int _status{};
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  fortran_coords_floor_to_curvilinear(
      /* void* */ floor_coords.get_fortran_ptr(),
      /* void* */ ele0.get_fortran_ptr(),
      /* void* */ _ele1.get_fortran_ptr(),
      /* int& */ _status,
      /* double* */ _w_mat_vec,
      /* void* */ local_coords.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsFloorToCurvilinear{std::move(_ele1), _status, w_mat};
}
Bmad::CoordsFloorToLocalCurvilinear Bmad::coords_floor_to_local_curvilinear(
    FloorPositionProxy& global_position,
    EleProxy& ele,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position) {
  int _status{};
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  int relative_to_lvalue;
  auto* _relative_to{&relative_to_lvalue};
  if (relative_to.has_value()) {
    relative_to_lvalue = relative_to.value();
  } else {
    _relative_to = nullptr;
  }
  fortran_coords_floor_to_local_curvilinear(
      /* void* */ global_position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ _status,
      /* double* */ _w_mat_vec,
      /* int* */ _relative_to,
      /* void* */ local_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsFloorToLocalCurvilinear{_status, w_mat};
}
void Bmad::coords_floor_to_relative(
    FloorPositionProxy& floor0,
    FloorPositionProxy& global_position,
    std::optional<bool> calculate_angles,
    std::optional<bool> is_delta_position,
    FloorPositionProxy& local_position) {
  bool calculate_angles_lvalue;
  auto* _calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  bool is_delta_position_lvalue;
  auto* _is_delta_position{&is_delta_position_lvalue};
  if (is_delta_position.has_value()) {
    is_delta_position_lvalue = is_delta_position.value();
  } else {
    _is_delta_position = nullptr;
  }
  fortran_coords_floor_to_relative(
      /* void* */ floor0.get_fortran_ptr(),
      /* void* */ global_position.get_fortran_ptr(),
      /* bool* */ _calculate_angles,
      /* bool* */ _is_delta_position,
      /* void* */ local_position.get_fortran_ptr());
}
void Bmad::coords_local_curvilinear_to_body(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& body_position) {
  double _w_mat_vec[3 * 3];
  const double* _w_mat = nullptr;
  if (w_mat.has_value()) {
    matrix_to_vec(w_mat.value(), _w_mat_vec);
    _w_mat = _w_mat_vec;
  }
  bool calculate_angles_lvalue;
  auto* _calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  fortran_coords_local_curvilinear_to_body(
      /* void* */ local_position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _w_mat_vec,
      /* bool* */ _calculate_angles,
      /* void* */ body_position.get_fortran_ptr());
}
FixedArray2D<Real, 3, 3> Bmad::coords_local_curvilinear_to_floor(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<bool> in_body_frame,
    std::optional<bool> calculate_angles,
    std::optional<int> relative_to,
    FloorPositionProxy& global_position) {
  bool in_body_frame_lvalue;
  auto* _in_body_frame{&in_body_frame_lvalue};
  if (in_body_frame.has_value()) {
    in_body_frame_lvalue = in_body_frame.value();
  } else {
    _in_body_frame = nullptr;
  }
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  bool calculate_angles_lvalue;
  auto* _calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  int relative_to_lvalue;
  auto* _relative_to{&relative_to_lvalue};
  if (relative_to.has_value()) {
    relative_to_lvalue = relative_to.value();
  } else {
    _relative_to = nullptr;
  }
  fortran_coords_local_curvilinear_to_floor(
      /* void* */ local_position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _in_body_frame,
      /* double* */ _w_mat_vec,
      /* bool* */ _calculate_angles,
      /* int* */ _relative_to,
      /* void* */ global_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
void Bmad::coords_relative_to_floor(
    FloorPositionProxy& floor0,
    FixedArray1D<Real, 3> dr,
    optional_ref<double> theta,
    optional_ref<double> phi,
    optional_ref<double> psi,
    FloorPositionProxy& floor1) {
  auto* _dr = dr.data(); // CppWrapperGeneralArgument
  auto* _theta = theta.has_value() ? &theta->get() : nullptr; // inout, optional
  auto* _phi = phi.has_value() ? &phi->get() : nullptr; // inout, optional
  auto* _psi = psi.has_value() ? &psi->get() : nullptr; // inout, optional
  fortran_coords_relative_to_floor(
      /* void* */ floor0.get_fortran_ptr(),
      /* double* */ _dr,
      /* double* */ _theta,
      /* double* */ _phi,
      /* double* */ _psi,
      /* void* */ floor1.get_fortran_ptr());
}
void Bmad::coulombfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& res) {
  fortran_coulombfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ res);
}
void Bmad::create_concatenated_wall3d(LatProxy& lat, bool& err) {
  fortran_create_concatenated_wall3d(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ err);
}
Bmad::CreateElementSlice Bmad::create_element_slice(
    EleProxy& ele_in,
    double l_slice,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end,
    optional_ref<EleProxy> old_slice,
    optional_ref<CoordProxy> orb_in) {
  EleProxy _sliced_ele;
  bool _err_flag{};
  auto* _old_slice = old_slice.has_value() ? old_slice->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  auto* _orb_in = orb_in.has_value() ? orb_in->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  fortran_create_element_slice(
      /* void* */ _sliced_ele.get_fortran_ptr(),
      /* void* */ ele_in.get_fortran_ptr(),
      /* double& */ l_slice,
      /* double& */ offset,
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ include_upstream_end,
      /* bool& */ include_downstream_end,
      /* bool& */ _err_flag,
      /* void* */ _old_slice,
      /* void* */ _orb_in);
  return CreateElementSlice{std::move(_sliced_ele), _err_flag};
}
bool Bmad::create_field_overlap(
    LatProxy& lat,
    std::string lord_name,
    std::string slave_name) {
  auto _lord_name = lord_name.c_str();
  auto _slave_name = slave_name.c_str();
  bool _err_flag{};
  fortran_create_field_overlap(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _lord_name,
      /* const char* */ _slave_name,
      /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::create_girder(
    LatProxy& lat,
    int ix_girder,
    ControlProxyAlloc1D& contrl,
    EleProxy& girder_info,
    bool& err_flag) {
  // intent=in allocatable type array
  fortran_create_girder(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ ix_girder,
      /* void* */ contrl.get_fortran_ptr(),
      /* void* */ girder_info.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::create_group(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err) {
  // intent=in allocatable type array
  fortran_create_group(
      /* void* */ lord.get_fortran_ptr(),
      /* void* */ contrl.get_fortran_ptr(),
      /* bool& */ err);
}
void Bmad::create_lat_ele_nametable(LatProxy& lat, NametableProxy& nametable) {
  fortran_create_lat_ele_nametable(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ nametable.get_fortran_ptr());
}
void Bmad::create_overlay(
    EleProxy& lord,
    ControlProxyAlloc1D& contrl,
    bool err) {
  // intent=in allocatable type array
  fortran_create_overlay(
      /* void* */ lord.get_fortran_ptr(),
      /* void* */ contrl.get_fortran_ptr(),
      /* bool& */ err);
}
Bmad::CreatePlanarWigglerModel Bmad::create_planar_wiggler_model(
    EleProxy& wiggler_in,
    std::optional<bool> print_err) {
  LatProxy _lat;
  bool _err_flag{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_create_planar_wiggler_model(
      /* void* */ wiggler_in.get_fortran_ptr(),
      /* void* */ _lat.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _print_err);
  return CreatePlanarWigglerModel{std::move(_lat), _err_flag};
}
void Bmad::create_ramper(
    EleProxy& lord,
    ControlProxyAlloc1D& contrl,
    bool err) {
  // intent=in allocatable type array
  fortran_create_ramper(
      /* void* */ lord.get_fortran_ptr(),
      /* void* */ contrl.get_fortran_ptr(),
      /* bool& */ err);
}
void Bmad::create_sol_quad_model(EleProxy& sol_quad, LatProxy& lat) {
  fortran_create_sol_quad_model(
      /* void* */ sol_quad.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr());
}
void Bmad::create_unique_ele_names(LatProxy& lat, int key, std::string suffix) {
  auto _suffix = suffix.c_str();
  fortran_create_unique_ele_names(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ key,
      /* const char* */ _suffix);
}
CartesianMapProxy Bmad::create_wiggler_cartesian_map(EleProxy& ele) {
  CartesianMapProxy _cart_map;
  fortran_create_wiggler_cartesian_map(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _cart_map.get_fortran_ptr());
  return std::move(_cart_map);
}
void Bmad::crystal_attribute_bookkeeper(EleProxy& ele) {
  fortran_crystal_attribute_bookkeeper(/* void* */ ele.get_fortran_ptr());
}
void Bmad::crystal_h_misalign(
    EleProxy& ele,
    CoordProxy& orbit,
    FixedArray1D<Real, 3> h_vec) {
  auto* _h_vec = h_vec.data(); // CppWrapperGeneralArgument
  fortran_crystal_h_misalign(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _h_vec);
}
bool Bmad::crystal_type_to_crystal_params(EleProxy& ele) {
  bool _err_flag{};
  fortran_crystal_type_to_crystal_params(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
int Bmad::custom_attribute_ubound_index(int ele_class) {
  int _ix_ubound{};
  fortran_custom_attribute_ubound_index(
      /* int& */ ele_class, /* int& */ _ix_ubound);
  return _ix_ubound;
}
void Bmad::damping_matrix_d(
    double& gamma,
    double& g_tot,
    double& B0,
    double& B1,
    double& delta,
    int& species,
    FixedArray2D<Real, 6, 6> mat) {
  double _mat_vec[6 * 6];
  matrix_to_vec(mat, _mat_vec);
  fortran_damping_matrix_d(
      /* double& */ gamma,
      /* double& */ g_tot,
      /* double& */ B0,
      /* double& */ B1,
      /* double& */ delta,
      /* int& */ species,
      /* double* */ _mat_vec);
  vec_to_matrix(_mat_vec, mat);
}
void Bmad::deallocate_ele_pointers(
    EleProxy& ele,
    std::optional<bool> nullify_only,
    std::optional<bool> nullify_branch,
    std::optional<bool> dealloc_poles) {
  bool nullify_only_lvalue;
  auto* _nullify_only{&nullify_only_lvalue};
  if (nullify_only.has_value()) {
    nullify_only_lvalue = nullify_only.value();
  } else {
    _nullify_only = nullptr;
  }
  bool nullify_branch_lvalue;
  auto* _nullify_branch{&nullify_branch_lvalue};
  if (nullify_branch.has_value()) {
    nullify_branch_lvalue = nullify_branch.value();
  } else {
    _nullify_branch = nullptr;
  }
  bool dealloc_poles_lvalue;
  auto* _dealloc_poles{&dealloc_poles_lvalue};
  if (dealloc_poles.has_value()) {
    dealloc_poles_lvalue = dealloc_poles.value();
  } else {
    _dealloc_poles = nullptr;
  }
  fortran_deallocate_ele_pointers(
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _nullify_only,
      /* bool* */ _nullify_branch,
      /* bool* */ _dealloc_poles);
}
void Bmad::deallocate_expression_tree(ExpressionTreeProxy& tree) {
  fortran_deallocate_expression_tree(/* void* */ tree.get_fortran_ptr());
}
void Bmad::deallocate_lat_pointers(LatProxy& lat) {
  fortran_deallocate_lat_pointers(/* void* */ lat.get_fortran_ptr());
}
void Bmad::default_tracking_species(LatParamProxy& param, int& species) {
  fortran_default_tracking_species(
      /* void* */ param.get_fortran_ptr(), /* int& */ species);
}
FixedArray1D<Int, 2> Bmad::detector_pixel_pt(CoordProxy& orbit, EleProxy& ele) {
  FixedArray1D<Int, 2> _ix_pix;
  fortran_detector_pixel_pt(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _ix_pix.data());
  return _ix_pix;
}
void Bmad::diffraction_plate_or_mask_hit_spot(
    EleProxy& ele,
    CoordProxy& orbit,
    int& ix_section) {
  fortran_diffraction_plate_or_mask_hit_spot(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ ix_section);
}
void Bmad::diffusion_matrix_b(
    double& gamma,
    double& g_tot,
    int& species,
    FixedArray2D<Real, 6, 6> mat) {
  double _mat_vec[6 * 6];
  matrix_to_vec(mat, _mat_vec);
  fortran_diffusion_matrix_b(
      /* double& */ gamma,
      /* double& */ g_tot,
      /* int& */ species,
      /* double* */ _mat_vec);
  vec_to_matrix(_mat_vec, mat);
}
bool Bmad::distance_to_aperture(
    CoordProxy& orbit,
    int particle_at,
    EleProxy& ele,
    double& dist) {
  bool _no_aperture_here{};
  fortran_distance_to_aperture(
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _no_aperture_here,
      /* double& */ dist);
  return _no_aperture_here;
}
bool Bmad::do_mode_flip(EleProxy& ele) {
  bool _err_flag{};
  fortran_do_mode_flip(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::dpc_given_de(double& pc_old, double& mass, double& dE, double& dpc) {
  fortran_dpc_given_de(
      /* double& */ pc_old,
      /* double& */ mass,
      /* double& */ dE,
      /* double& */ dpc);
}
void Bmad::drift_and_pipe_track_methods_adjustment(LatProxy& lat) {
  fortran_drift_and_pipe_track_methods_adjustment(
      /* void* */ lat.get_fortran_ptr());
}
void Bmad::drift_multipass_name_correction(LatProxy& lat) {
  fortran_drift_multipass_name_correction(/* void* */ lat.get_fortran_ptr());
}
void Bmad::drift_orbit_time(
    CoordProxy& orbit,
    double beta0,
    std::optional<double> delta_s,
    std::optional<double> delta_t) {
  double delta_s_lvalue;
  auto* _delta_s{&delta_s_lvalue};
  if (delta_s.has_value()) {
    delta_s_lvalue = delta_s.value();
  } else {
    _delta_s = nullptr;
  }
  double delta_t_lvalue;
  auto* _delta_t{&delta_t_lvalue};
  if (delta_t.has_value()) {
    delta_t_lvalue = delta_t.value();
  } else {
    _delta_t = nullptr;
  }
  fortran_drift_orbit_time(
      /* void* */ orbit.get_fortran_ptr(),
      /* double& */ beta0,
      /* double* */ _delta_s,
      /* double* */ _delta_t);
}
void Bmad::drift_particle_to_s(CoordProxy& p, double s, BranchProxy& branch) {
  fortran_drift_particle_to_s(
      /* void* */ p.get_fortran_ptr(),
      /* double& */ s,
      /* void* */ branch.get_fortran_ptr());
}
void Bmad::drift_particle_to_t(CoordProxy& p, double t, BranchProxy& branch) {
  fortran_drift_particle_to_t(
      /* void* */ p.get_fortran_ptr(),
      /* double& */ t,
      /* void* */ branch.get_fortran_ptr());
}
double Bmad::dspline_len(
    double s_chord0,
    double s_chord1,
    SplineProxy& spline,
    std::optional<double> dtheta_ref) {
  double dtheta_ref_lvalue;
  auto* _dtheta_ref{&dtheta_ref_lvalue};
  if (dtheta_ref.has_value()) {
    dtheta_ref_lvalue = dtheta_ref.value();
  } else {
    _dtheta_ref = nullptr;
  }
  double _dlen{};
  fortran_dspline_len(
      /* double& */ s_chord0,
      /* double& */ s_chord1,
      /* void* */ spline.get_fortran_ptr(),
      /* double* */ _dtheta_ref,
      /* double& */ _dlen);
  return _dlen;
}
AperturePointProxy Bmad::dynamic_aperture_point(
    BranchProxy& branch,
    EleProxy& ele0,
    CoordProxy& orb0,
    double theta_xy,
    ApertureParamProxy& ap_param,
    std::optional<bool> check_xy_init) {
  AperturePointProxy _ap_point;
  bool check_xy_init_lvalue;
  auto* _check_xy_init{&check_xy_init_lvalue};
  if (check_xy_init.has_value()) {
    check_xy_init_lvalue = check_xy_init.value();
  } else {
    _check_xy_init = nullptr;
  }
  fortran_dynamic_aperture_point(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ ele0.get_fortran_ptr(),
      /* void* */ orb0.get_fortran_ptr(),
      /* double& */ theta_xy,
      /* void* */ ap_param.get_fortran_ptr(),
      /* void* */ _ap_point.get_fortran_ptr(),
      /* bool* */ _check_xy_init);
  return std::move(_ap_point);
}
ApertureScanProxyAlloc1D Bmad::dynamic_aperture_scan(
    ApertureParamProxy& aperture_param,
    RealAlloc1D& pz_start,
    LatProxy& lat,
    std::optional<bool> print_timing) {
  // intent=out allocatable type array
  auto aperture_scan{ApertureScanProxyAlloc1D()};
  // intent=in allocatable general array
  bool print_timing_lvalue;
  auto* _print_timing{&print_timing_lvalue};
  if (print_timing.has_value()) {
    print_timing_lvalue = print_timing.value();
  } else {
    _print_timing = nullptr;
  }
  fortran_dynamic_aperture_scan(
      /* void* */ aperture_scan.get_fortran_ptr(),
      /* void* */ aperture_param.get_fortran_ptr(),
      /* void* */ pz_start.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _print_timing);
  return std::move(aperture_scan);
}
void Bmad::e_accel_field(
    EleProxy& ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking,
    double& field) {
  bool bmad_standard_tracking_lvalue;
  auto* _bmad_standard_tracking{&bmad_standard_tracking_lvalue};
  if (bmad_standard_tracking.has_value()) {
    bmad_standard_tracking_lvalue = bmad_standard_tracking.value();
  } else {
    _bmad_standard_tracking = nullptr;
  }
  fortran_e_accel_field(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ voltage_or_gradient,
      /* bool* */ _bmad_standard_tracking,
      /* double& */ field);
}
double Bmad::e_crit_photon(double gamma, double g_bend) {
  double _E_crit{};
  fortran_e_crit_photon(
      /* double& */ gamma, /* double& */ g_bend, /* double& */ _E_crit);
  return _E_crit;
}
Bmad::EigenDecomp6mat Bmad::eigen_decomp_6mat(FixedArray2D<Real, 6, 6> mat) {
  double _mat_vec[6 * 6];
  matrix_to_vec(mat, _mat_vec);
  FixedArray1D<Complex, 6> _eval;
  FixedArray2D<Complex, 6, 6> evec;
  std::complex<double> _evec_vec[6 * 6];
  bool _err_flag{};
  FixedArray1D<Real, 3> _tunes;
  fortran_eigen_decomp_6mat(
      /* double* */ _mat_vec,
      /* std::complex<double>* */ _eval.data(),
      /* std::complex<double>* */ _evec_vec,
      /* bool& */ _err_flag,
      /* double* */ _tunes.data());
  vec_to_matrix(_evec_vec, evec);
  return EigenDecomp6mat{_eval, evec, _err_flag, _tunes};
}
void Bmad::ele_compute_ref_energy_and_time(
    EleProxy& ele0,
    EleProxy& ele,
    LatParamProxy& param,
    bool err_flag) {
  fortran_ele_compute_ref_energy_and_time(
      /* void* */ ele0.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::ele_equal_ele(EleProxy& ele_out, EleProxy& ele_in) {
  fortran_ele_equal_ele(
      /* void* */ ele_out.get_fortran_ptr(),
      /* void* */ ele_in.get_fortran_ptr());
}
EleProxy Bmad::ele_equals_ele(EleProxy& ele_in, bool update_nametable) {
  EleProxy _ele_out;
  fortran_ele_equals_ele(
      /* void* */ _ele_out.get_fortran_ptr(),
      /* void* */ ele_in.get_fortran_ptr(),
      /* bool& */ update_nametable);
  return std::move(_ele_out);
}
void Bmad::ele_finalizer(EleProxy& ele) {
  fortran_ele_finalizer(/* void* */ ele.get_fortran_ptr());
}
void Bmad::ele_full_name(
    EleProxy& ele,
    std::optional<std::string> template_,
    std::string& str) {
  const char* _template_ = template_.has_value() ? template_->c_str() : nullptr;
  auto _str = str.c_str(); // ptr, inout, required
  fortran_ele_full_name(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _template_,
      /* const char* */ _str);
}
FloorPositionProxy Bmad::ele_geometry(
    FloorPositionProxy& floor_start,
    EleProxy& ele,
    std::optional<double> len_scale,
    std::optional<bool> ignore_patch_err) {
  FloorPositionProxy _floor_end;
  double len_scale_lvalue;
  auto* _len_scale{&len_scale_lvalue};
  if (len_scale.has_value()) {
    len_scale_lvalue = len_scale.value();
  } else {
    _len_scale = nullptr;
  }
  bool ignore_patch_err_lvalue;
  auto* _ignore_patch_err{&ignore_patch_err_lvalue};
  if (ignore_patch_err.has_value()) {
    ignore_patch_err_lvalue = ignore_patch_err.value();
  } else {
    _ignore_patch_err = nullptr;
  }
  fortran_ele_geometry(
      /* void* */ floor_start.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _floor_end.get_fortran_ptr(),
      /* double* */ _len_scale,
      /* bool* */ _ignore_patch_err);
  return std::move(_floor_end);
}
void Bmad::ele_geometry_with_misalignments(
    EleProxy& ele,
    std::optional<double> len_scale,
    FloorPositionProxy& floor) {
  double len_scale_lvalue;
  auto* _len_scale{&len_scale_lvalue};
  if (len_scale.has_value()) {
    len_scale_lvalue = len_scale.value();
  } else {
    _len_scale = nullptr;
  }
  fortran_ele_geometry_with_misalignments(
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _len_scale,
      /* void* */ floor.get_fortran_ptr());
}
void Bmad::ele_has_constant_ds_dt_ref(EleProxy& ele, bool& is_const) {
  fortran_ele_has_constant_ds_dt_ref(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ is_const);
}
EleProxy Bmad::ele_has_nonzero_kick(bool& has_kick) {
  EleProxy _ele;
  fortran_ele_has_nonzero_kick(
      /* void* */ _ele.get_fortran_ptr(), /* bool& */ has_kick);
  return std::move(_ele);
}
void Bmad::ele_has_nonzero_offset(EleProxy& ele, bool& has_offset) {
  fortran_ele_has_nonzero_offset(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ has_offset);
}
bool Bmad::ele_is_monitor(EleProxy& ele, std::optional<bool> print_warning) {
  bool print_warning_lvalue;
  auto* _print_warning{&print_warning_lvalue};
  if (print_warning.has_value()) {
    print_warning_lvalue = print_warning.value();
  } else {
    _print_warning = nullptr;
  }
  bool _is_monitor{};
  fortran_ele_is_monitor(
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _print_warning,
      /* bool& */ _is_monitor);
  return _is_monitor;
}
void Bmad::ele_loc(EleProxy& ele, LatEleLocProxy& loc) {
  fortran_ele_loc(
      /* void* */ ele.get_fortran_ptr(), /* void* */ loc.get_fortran_ptr());
}
void Bmad::ele_loc_name(
    EleProxy& ele,
    std::optional<bool> show_branch0,
    std::optional<std::string> parens,
    std::string& str) {
  bool show_branch0_lvalue;
  auto* _show_branch0{&show_branch0_lvalue};
  if (show_branch0.has_value()) {
    show_branch0_lvalue = show_branch0.value();
  } else {
    _show_branch0 = nullptr;
  }
  const char* _parens = parens.has_value() ? parens->c_str() : nullptr;
  auto _str = str.c_str(); // ptr, inout, required
  fortran_ele_loc_name(
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _show_branch0,
      /* const char* */ _parens,
      /* const char* */ _str);
}
Bmad::EleMisalignmentLSCalc Bmad::ele_misalignment_l_s_calc(EleProxy& ele) {
  FixedArray1D<Real, 3> _L_mis;
  FixedArray2D<Real, 3, 3> S_mis;
  double _S_mis_vec[3 * 3];
  fortran_ele_misalignment_l_s_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _L_mis.data(),
      /* double* */ _S_mis_vec);
  vec_to_matrix(_S_mis_vec, S_mis);
  return EleMisalignmentLSCalc{_L_mis, S_mis};
}
void Bmad::ele_nametable_index(EleProxy& ele, int& ix_nt) {
  fortran_ele_nametable_index(
      /* void* */ ele.get_fortran_ptr(), /* int& */ ix_nt);
}
LatEleOrderProxy Bmad::ele_order_calc(LatProxy& lat) {
  LatEleOrderProxy _order;
  fortran_ele_order_calc(
      /* void* */ lat.get_fortran_ptr(), /* void* */ _order.get_fortran_ptr());
  return std::move(_order);
}
void Bmad::ele_reference_energy_correction(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_ele_reference_energy_correction(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ particle_at,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::ele_rf_step_index(
    double E_ref,
    double s_rel,
    EleProxy& ele,
    int& ix_step) {
  fortran_ele_rf_step_index(
      /* double& */ E_ref,
      /* double& */ s_rel,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ ix_step);
}
Bmad::EleToPtcMagneticBnAn Bmad::ele_to_ptc_magnetic_bn_an(EleProxy& ele) {
  // intent=out allocatable general array
  auto bn{RealAlloc1D()};
  // intent=out allocatable general array
  auto an{RealAlloc1D()};
  int _n_max{};
  fortran_ele_to_ptc_magnetic_bn_an(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ bn.get_fortran_ptr(),
      /* void* */ an.get_fortran_ptr(),
      /* int& */ _n_max);
  return EleToPtcMagneticBnAn{std::move(bn), std::move(an), _n_max};
}
void Bmad::ele_to_spin_taylor(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orb0) {
  fortran_ele_to_spin_taylor(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orb0.get_fortran_ptr());
}
Bmad::EleToTaylor Bmad::ele_to_taylor(
    EleProxy& ele,
    optional_ref<CoordProxy> orb0,
    std::optional<bool> taylor_map_includes_offsets,
    std::optional<bool> include_damping) {
  auto* _orb0 = orb0.has_value() ? orb0->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  bool taylor_map_includes_offsets_lvalue;
  auto* _taylor_map_includes_offsets{&taylor_map_includes_offsets_lvalue};
  if (taylor_map_includes_offsets.has_value()) {
    taylor_map_includes_offsets_lvalue = taylor_map_includes_offsets.value();
  } else {
    _taylor_map_includes_offsets = nullptr;
  }
  bool include_damping_lvalue;
  auto* _include_damping{&include_damping_lvalue};
  if (include_damping.has_value()) {
    include_damping_lvalue = include_damping.value();
  } else {
    _include_damping = nullptr;
  }
  // Output-only type array
  auto orbital_taylor = TaylorProxyArray1D::allocate(6, 1);

  // Output-only type array
  auto spin_taylor = TaylorProxyArray1D::allocate(4, 1);

  fortran_ele_to_taylor(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _orb0,
      /* bool* */ _taylor_map_includes_offsets,
      /* bool* */ _include_damping,
      /* void* */ orbital_taylor.get_fortran_ptr(),
      /* void* */ spin_taylor.get_fortran_ptr());
  return EleToTaylor{
      std::move(std::move(orbital_taylor)), std::move(std::move(spin_taylor))};
}
void Bmad::ele_unique_name(
    EleProxy& ele,
    LatEleOrderProxy& order,
    std::string& unique_name) {
  auto _unique_name = unique_name.c_str(); // ptr, inout, required
  fortran_ele_unique_name(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ order.get_fortran_ptr(),
      /* const char* */ _unique_name);
}
void Bmad::ele_value_has_changed(
    EleProxy& ele,
    IntAlloc1D& list,
    RealAlloc1D& abs_tol,
    bool set_old,
    bool& has_changed) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  fortran_ele_value_has_changed(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ list.get_fortran_ptr(),
      /* void* */ abs_tol.get_fortran_ptr(),
      /* bool& */ set_old,
      /* bool& */ has_changed);
}
void Bmad::ele_vec_equal_ele_vec(EleProxyAlloc1D& ele1, EleProxyAlloc1D& ele2) {
  // intent=inout allocatable type array
  // intent=in allocatable type array
  fortran_ele_vec_equal_ele_vec(
      /* void* */ ele1.get_fortran_ptr(), /* void* */ ele2.get_fortran_ptr());
}
Bmad::ElecMultipoleField Bmad::elec_multipole_field(
    double a,
    double b,
    int n,
    CoordProxy& coord) {
  double _Ex{};
  double _Ey{};
  FixedArray2D<Real, 2, 2> dE;
  double _dE_vec[2 * 2];
  bool _compute_dE{};
  fortran_elec_multipole_field(
      /* double& */ a,
      /* double& */ b,
      /* int& */ n,
      /* void* */ coord.get_fortran_ptr(),
      /* double& */ _Ex,
      /* double& */ _Ey,
      /* double* */ _dE_vec,
      /* bool& */ _compute_dE);
  vec_to_matrix(_dE_vec, dE);
  return ElecMultipoleField{_Ex, _Ey, dE, _compute_dE};
}
Bmad::ElementAtSBranch Bmad::element_at_s(
    BranchProxy& branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err) {
  bool _err_flag{};
  double _s_eff{};
  CoordProxy _position;
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _ix_ele{};
  fortran_element_at_s_branch(
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ s,
      /* bool& */ choose_max,
      /* bool& */ _err_flag,
      /* double& */ _s_eff,
      /* void* */ _position.get_fortran_ptr(),
      /* bool* */ _print_err,
      /* int& */ _ix_ele);
  return ElementAtSBranch{_err_flag, _s_eff, std::move(_position), _ix_ele};
}
Bmad::ElementAtSLat Bmad::element_at_s(
    LatProxy& lat,
    double s,
    bool choose_max,
    std::optional<int> ix_branch,
    std::optional<bool> print_err) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  double _s_eff{};
  CoordProxy _position;
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _ix_ele{};
  fortran_element_at_s_lat(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ s,
      /* bool& */ choose_max,
      /* int* */ _ix_branch,
      /* bool& */ _err_flag,
      /* double& */ _s_eff,
      /* void* */ _position.get_fortran_ptr(),
      /* bool* */ _print_err,
      /* int& */ _ix_ele);
  return ElementAtSLat{_err_flag, _s_eff, std::move(_position), _ix_ele};
}
void Bmad::element_slice_iterator(
    EleProxy& ele,
    LatParamProxy& param,
    int i_slice,
    int n_slice_tot,
    EleProxy& sliced_ele,
    std::optional<double> s_start,
    std::optional<double> s_end) {
  double s_start_lvalue;
  auto* _s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto* _s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  fortran_element_slice_iterator(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ i_slice,
      /* int& */ n_slice_tot,
      /* void* */ sliced_ele.get_fortran_ptr(),
      /* double* */ _s_start,
      /* double* */ _s_end);
}
void Bmad::ellipinc_test() {
  fortran_ellipinc_test();
}
Bmad::EmFieldCalc Bmad::em_field_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_pos,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential,
    std::optional<bool> use_overlap,
    std::optional<bool> grid_allow_s_out_of_bounds,
    std::optional<double> rf_time,
    optional_ref<ElePointerProxyAlloc1D> used_eles,
    std::optional<bool> print_err,
    optional_ref<EleProxy> original_ele) {
  EmFieldProxy _field;
  bool calc_dfield_lvalue;
  auto* _calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool _err_flag{};
  bool calc_potential_lvalue;
  auto* _calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  bool use_overlap_lvalue;
  auto* _use_overlap{&use_overlap_lvalue};
  if (use_overlap.has_value()) {
    use_overlap_lvalue = use_overlap.value();
  } else {
    _use_overlap = nullptr;
  }
  bool grid_allow_s_out_of_bounds_lvalue;
  auto* _grid_allow_s_out_of_bounds{&grid_allow_s_out_of_bounds_lvalue};
  if (grid_allow_s_out_of_bounds.has_value()) {
    grid_allow_s_out_of_bounds_lvalue = grid_allow_s_out_of_bounds.value();
  } else {
    _grid_allow_s_out_of_bounds = nullptr;
  }
  double rf_time_lvalue;
  auto* _rf_time{&rf_time_lvalue};
  if (rf_time.has_value()) {
    rf_time_lvalue = rf_time.value();
  } else {
    _rf_time = nullptr;
  }
  // intent=in allocatable type array
  auto* _used_eles = used_eles.has_value() ? used_eles->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  auto* _original_ele = original_ele.has_value()
      ? original_ele->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_em_field_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ s_pos,
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ local_ref_frame,
      /* void* */ _field.get_fortran_ptr(),
      /* bool* */ _calc_dfield,
      /* bool& */ _err_flag,
      /* bool* */ _calc_potential,
      /* bool* */ _use_overlap,
      /* bool* */ _grid_allow_s_out_of_bounds,
      /* double* */ _rf_time,
      /* void* */ _used_eles,
      /* bool* */ _print_err,
      /* void* */ _original_ele);
  return EmFieldCalc{std::move(_field), _err_flag};
}
EmFieldProxy Bmad::em_field_derivatives(
    EleProxy& ele,
    LatParamProxy& param,
    double& s_pos,
    CoordProxy& orbit,
    bool& local_ref_frame,
    optional_ref<bool> grid_allow_s_out_of_bounds,
    optional_ref<double> rf_time) {
  EmFieldProxy _dfield;
  auto* _grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds.has_value()
      ? &grid_allow_s_out_of_bounds->get()
      : nullptr; // inout, optional
  auto* _rf_time =
      rf_time.has_value() ? &rf_time->get() : nullptr; // inout, optional
  fortran_em_field_derivatives(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ s_pos,
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ local_ref_frame,
      /* void* */ _dfield.get_fortran_ptr(),
      /* bool* */ _grid_allow_s_out_of_bounds,
      /* double* */ _rf_time);
  return std::move(_dfield);
}
FixedArray1D<Real, 10> Bmad::em_field_kick_vector_time(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orbit,
    bool err_flag,
    std::optional<bool> print_err,
    optional_ref<EmFieldProxy> extra_field) {
  FixedArray1D<Real, 10> _dvec_dt;
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  auto* _extra_field = extra_field.has_value()
      ? extra_field->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_em_field_kick_vector_time(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ rf_time,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _dvec_dt.data(),
      /* bool& */ err_flag,
      /* bool* */ _print_err,
      /* void* */ _extra_field);
  return _dvec_dt;
}
void Bmad::em_field_plus_em_field(
    EmFieldProxy& field1,
    EmFieldProxy& field2,
    EmFieldProxy& field_tot) {
  fortran_em_field_plus_em_field(
      /* void* */ field1.get_fortran_ptr(),
      /* void* */ field2.get_fortran_ptr(),
      /* void* */ field_tot.get_fortran_ptr());
}
void Bmad::em_taylor_equal_em_taylor(
    EmTaylorProxy& em_taylor1,
    EmTaylorProxy& em_taylor2) {
  fortran_em_taylor_equal_em_taylor(
      /* void* */ em_taylor1.get_fortran_ptr(),
      /* void* */ em_taylor2.get_fortran_ptr());
}
void Bmad::em_taylors_equal_em_taylors(
    EmTaylorProxyAlloc1D& em_taylor1,
    EmTaylorProxyAlloc1D& em_taylor2) {
  // intent=inout allocatable type array
  // intent=in allocatable type array
  fortran_em_taylors_equal_em_taylors(
      /* void* */ em_taylor1.get_fortran_ptr(),
      /* void* */ em_taylor2.get_fortran_ptr());
}
Bmad::Emit6d Bmad::emit_6d(
    EleProxy& ele_ref,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit) {
  NormalModesProxy _mode;
  FixedArray2D<Real, 6, 6> sigma_mat;
  double _sigma_mat_vec[6 * 6];
  // intent=in allocatable type array
  auto* _closed_orbit = closed_orbit.has_value()
      ? closed_orbit->get().get_fortran_ptr()
      : nullptr; // input, optional
  RadIntAllEleProxy _rad_int_by_ele;
  fortran_emit_6d(
      /* void* */ ele_ref.get_fortran_ptr(),
      /* bool& */ include_opening_angle,
      /* void* */ _mode.get_fortran_ptr(),
      /* double* */ _sigma_mat_vec,
      /* void* */ _closed_orbit,
      /* void* */ _rad_int_by_ele.get_fortran_ptr());
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  return Emit6d{std::move(_mode), sigma_mat, std::move(_rad_int_by_ele)};
}
void Bmad::entering_element(
    CoordProxy& orbit,
    int particle_at,
    bool& is_entering) {
  fortran_entering_element(
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ particle_at,
      /* bool& */ is_entering);
}
void Bmad::envelope_radints(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    FixedArray1D<Real, 3> alpha,
    FixedArray1D<Real, 3> emit) {
  std::complex<double> _Lambda_vec[6 * 6];
  matrix_to_vec(Lambda, _Lambda_vec);
  std::complex<double> _Theta_vec[6 * 6];
  matrix_to_vec(Theta, _Theta_vec);
  std::complex<double> _Iota_vec[6 * 6];
  matrix_to_vec(Iota, _Iota_vec);
  auto* _alpha = alpha.data(); // CppWrapperGeneralArgument
  auto* _emit = emit.data(); // CppWrapperGeneralArgument
  fortran_envelope_radints(
      /* std::complex<double>* */ _Lambda_vec,
      /* std::complex<double>* */ _Theta_vec,
      /* std::complex<double>* */ _Iota_vec,
      /* double* */ _alpha,
      /* double* */ _emit);
  vec_to_matrix(_Lambda_vec, Lambda);
  vec_to_matrix(_Theta_vec, Theta);
  vec_to_matrix(_Iota_vec, Iota);
}
Bmad::EnvelopeRadintsIbs Bmad::envelope_radints_ibs(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    EleProxyAlloc1D& eles,
    NormalModesProxy& mode,
    bool tail_cut,
    double npart,
    int species) {
  std::complex<double> _Lambda_vec[6 * 6];
  matrix_to_vec(Lambda, _Lambda_vec);
  std::complex<double> _Theta_vec[6 * 6];
  matrix_to_vec(Theta, _Theta_vec);
  std::complex<double> _Iota_vec[6 * 6];
  matrix_to_vec(Iota, _Iota_vec);
  // intent=in allocatable type array
  FixedArray1D<Real, 3> _alpha;
  FixedArray1D<Real, 3> _emit;
  fortran_envelope_radints_ibs(
      /* std::complex<double>* */ _Lambda_vec,
      /* std::complex<double>* */ _Theta_vec,
      /* std::complex<double>* */ _Iota_vec,
      /* void* */ eles.get_fortran_ptr(),
      /* double* */ _alpha.data(),
      /* double* */ _emit.data(),
      /* void* */ mode.get_fortran_ptr(),
      /* bool& */ tail_cut,
      /* double& */ npart,
      /* int& */ species);
  return EnvelopeRadintsIbs{_alpha, _emit};
}
void Bmad::eq_ac_kicker(AcKickerProxy& f1, AcKickerProxy& f2, bool& is_eq) {
  fortran_eq_ac_kicker(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_ac_kicker_freq(
    AcKickerFreqProxy& f1,
    AcKickerFreqProxy& f2,
    bool& is_eq) {
  fortran_eq_ac_kicker_freq(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_ac_kicker_time(
    AcKickerTimeProxy& f1,
    AcKickerTimeProxy& f2,
    bool& is_eq) {
  fortran_eq_ac_kicker_time(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_anormal_mode(
    AnormalModeProxy& f1,
    AnormalModeProxy& f2,
    bool& is_eq) {
  fortran_eq_anormal_mode(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_aperture_param(
    ApertureParamProxy& f1,
    ApertureParamProxy& f2,
    bool& is_eq) {
  fortran_eq_aperture_param(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_aperture_point(
    AperturePointProxy& f1,
    AperturePointProxy& f2,
    bool& is_eq) {
  fortran_eq_aperture_point(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_aperture_scan(
    ApertureScanProxy& f1,
    ApertureScanProxy& f2,
    bool& is_eq) {
  fortran_eq_aperture_scan(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_beam(BeamProxy& f1, BeamProxy& f2, bool& is_eq) {
  fortran_eq_beam(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_beam_init(BeamInitProxy& f1, BeamInitProxy& f2, bool& is_eq) {
  fortran_eq_beam_init(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_bmad_common(
    BmadCommonProxy& f1,
    BmadCommonProxy& f2,
    bool& is_eq) {
  fortran_eq_bmad_common(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_bookkeeping_state(
    BookkeepingStateProxy& f1,
    BookkeepingStateProxy& f2,
    bool& is_eq) {
  fortran_eq_bookkeeping_state(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_bpm_phase_coupling(
    BpmPhaseCouplingProxy& f1,
    BpmPhaseCouplingProxy& f2,
    bool& is_eq) {
  fortran_eq_bpm_phase_coupling(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_branch(BranchProxy& f1, BranchProxy& f2, bool& is_eq) {
  fortran_eq_branch(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_bunch(BunchProxy& f1, BunchProxy& f2, bool& is_eq) {
  fortran_eq_bunch(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_bunch_params(
    BunchParamsProxy& f1,
    BunchParamsProxy& f2,
    bool& is_eq) {
  fortran_eq_bunch_params(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map(
    CartesianMapProxy& f1,
    CartesianMapProxy& f2,
    bool& is_eq) {
  fortran_eq_cartesian_map(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map_term(
    CartesianMapTermProxy& f1,
    CartesianMapTermProxy& f2,
    bool& is_eq) {
  fortran_eq_cartesian_map_term(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map_term1(
    CartesianMapTerm1Proxy& f1,
    CartesianMapTerm1Proxy& f2,
    bool& is_eq) {
  fortran_eq_cartesian_map_term1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_complex_taylor(
    ComplexTaylorProxy& f1,
    ComplexTaylorProxy& f2,
    bool& is_eq) {
  fortran_eq_complex_taylor(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_complex_taylor_term(
    ComplexTaylorTermProxy& f1,
    ComplexTaylorTermProxy& f2,
    bool& is_eq) {
  fortran_eq_complex_taylor_term(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_control(ControlProxy& f1, ControlProxy& f2, bool& is_eq) {
  fortran_eq_control(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_control_ramp1(
    ControlRamp1Proxy& f1,
    ControlRamp1Proxy& f2,
    bool& is_eq) {
  fortran_eq_control_ramp1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_control_var1(
    ControlVar1Proxy& f1,
    ControlVar1Proxy& f2,
    bool& is_eq) {
  fortran_eq_control_var1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_controller(
    ControllerProxy& f1,
    ControllerProxy& f2,
    bool& is_eq) {
  fortran_eq_controller(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_coord(CoordProxy& f1, CoordProxy& f2, bool& is_eq) {
  fortran_eq_coord(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_coord_array(
    CoordArrayProxy& f1,
    CoordArrayProxy& f2,
    bool& is_eq) {
  fortran_eq_coord_array(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map(
    CylindricalMapProxy& f1,
    CylindricalMapProxy& f2,
    bool& is_eq) {
  fortran_eq_cylindrical_map(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map_term(
    CylindricalMapTermProxy& f1,
    CylindricalMapTermProxy& f2,
    bool& is_eq) {
  fortran_eq_cylindrical_map_term(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map_term1(
    CylindricalMapTerm1Proxy& f1,
    CylindricalMapTerm1Proxy& f2,
    bool& is_eq) {
  fortran_eq_cylindrical_map_term1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_ele(EleProxy& f1, EleProxy& f2, bool& is_eq) {
  fortran_eq_ele(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_ellipse_beam_init(
    EllipseBeamInitProxy& f1,
    EllipseBeamInitProxy& f2,
    bool& is_eq) {
  fortran_eq_ellipse_beam_init(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_em_field(EmFieldProxy& f1, EmFieldProxy& f2, bool& is_eq) {
  fortran_eq_em_field(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_em_taylor(EmTaylorProxy& f1, EmTaylorProxy& f2, bool& is_eq) {
  fortran_eq_em_taylor(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_em_taylor_term(
    EmTaylorTermProxy& f1,
    EmTaylorTermProxy& f2,
    bool& is_eq) {
  fortran_eq_em_taylor_term(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_expression_atom(
    ExpressionAtomProxy& f1,
    ExpressionAtomProxy& f2,
    bool& is_eq) {
  fortran_eq_expression_atom(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_floor_position(
    FloorPositionProxy& f1,
    FloorPositionProxy& f2,
    bool& is_eq) {
  fortran_eq_floor_position(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_gen_grad1(GenGrad1Proxy& f1, GenGrad1Proxy& f2, bool& is_eq) {
  fortran_eq_gen_grad1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_gen_grad_map(
    GenGradMapProxy& f1,
    GenGradMapProxy& f2,
    bool& is_eq) {
  fortran_eq_gen_grad_map(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_grid_beam_init(
    GridBeamInitProxy& f1,
    GridBeamInitProxy& f2,
    bool& is_eq) {
  fortran_eq_grid_beam_init(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_grid_field(GridFieldProxy& f1, GridFieldProxy& f2, bool& is_eq) {
  fortran_eq_grid_field(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_grid_field_pt(
    GridFieldPtProxy& f1,
    GridFieldPtProxy& f2,
    bool& is_eq) {
  fortran_eq_grid_field_pt(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_grid_field_pt1(
    GridFieldPt1Proxy& f1,
    GridFieldPt1Proxy& f2,
    bool& is_eq) {
  fortran_eq_grid_field_pt1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_high_energy_space_charge(
    HighEnergySpaceChargeProxy& f1,
    HighEnergySpaceChargeProxy& f2,
    bool& is_eq) {
  fortran_eq_high_energy_space_charge(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_interval1_coef(
    Interval1CoefProxy& f1,
    Interval1CoefProxy& f2,
    bool& is_eq) {
  fortran_eq_interval1_coef(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_kv_beam_init(
    KvBeamInitProxy& f1,
    KvBeamInitProxy& f2,
    bool& is_eq) {
  fortran_eq_kv_beam_init(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_lat(LatProxy& f1, LatProxy& f2, bool& is_eq) {
  fortran_eq_lat(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_lat_ele_loc(LatEleLocProxy& f1, LatEleLocProxy& f2, bool& is_eq) {
  fortran_eq_lat_ele_loc(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_lat_param(LatParamProxy& f1, LatParamProxy& f2, bool& is_eq) {
  fortran_eq_lat_param(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_linac_normal_mode(
    LinacNormalModeProxy& f1,
    LinacNormalModeProxy& f2,
    bool& is_eq) {
  fortran_eq_linac_normal_mode(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_mode3(Mode3Proxy& f1, Mode3Proxy& f2, bool& is_eq) {
  fortran_eq_mode3(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_mode_info(ModeInfoProxy& f1, ModeInfoProxy& f2, bool& is_eq) {
  fortran_eq_mode_info(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_normal_modes(
    NormalModesProxy& f1,
    NormalModesProxy& f2,
    bool& is_eq) {
  fortran_eq_normal_modes(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_photon_element(
    PhotonElementProxy& f1,
    PhotonElementProxy& f2,
    bool& is_eq) {
  fortran_eq_photon_element(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_photon_material(
    PhotonMaterialProxy& f1,
    PhotonMaterialProxy& f2,
    bool& is_eq) {
  fortran_eq_photon_material(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_photon_reflect_surface(
    PhotonReflectSurfaceProxy& f1,
    PhotonReflectSurfaceProxy& f2,
    bool& is_eq) {
  fortran_eq_photon_reflect_surface(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_photon_reflect_table(
    PhotonReflectTableProxy& f1,
    PhotonReflectTableProxy& f2,
    bool& is_eq) {
  fortran_eq_photon_reflect_table(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_photon_target(
    PhotonTargetProxy& f1,
    PhotonTargetProxy& f2,
    bool& is_eq) {
  fortran_eq_photon_target(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_pixel_detec(
    PixelDetecProxy& f1,
    PixelDetecProxy& f2,
    bool& is_eq) {
  fortran_eq_pixel_detec(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_pixel_pt(PixelPtProxy& f1, PixelPtProxy& f2, bool& is_eq) {
  fortran_eq_pixel_pt(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_pre_tracker(
    PreTrackerProxy& f1,
    PreTrackerProxy& f2,
    bool& is_eq) {
  fortran_eq_pre_tracker(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_rad_int1(RadInt1Proxy& f1, RadInt1Proxy& f2, bool& is_eq) {
  fortran_eq_rad_int1(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_rad_int_all_ele(
    RadIntAllEleProxy& f1,
    RadIntAllEleProxy& f2,
    bool& is_eq) {
  fortran_eq_rad_int_all_ele(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_rad_int_branch(
    RadIntBranchProxy& f1,
    RadIntBranchProxy& f2,
    bool& is_eq) {
  fortran_eq_rad_int_branch(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_rad_map(RadMapProxy& f1, RadMapProxy& f2, bool& is_eq) {
  fortran_eq_rad_map(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_rad_map_ele(RadMapEleProxy& f1, RadMapEleProxy& f2, bool& is_eq) {
  fortran_eq_rad_map_ele(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_ramper_lord(
    RamperLordProxy& f1,
    RamperLordProxy& f2,
    bool& is_eq) {
  fortran_eq_ramper_lord(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_space_charge_common(
    SpaceChargeCommonProxy& f1,
    SpaceChargeCommonProxy& f2,
    bool& is_eq) {
  fortran_eq_space_charge_common(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_spin_polar(SpinPolarProxy& f1, SpinPolarProxy& f2, bool& is_eq) {
  fortran_eq_spin_polar(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_spline(SplineProxy& f1, SplineProxy& f2, bool& is_eq) {
  fortran_eq_spline(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_strong_beam(
    StrongBeamProxy& f1,
    StrongBeamProxy& f2,
    bool& is_eq) {
  fortran_eq_strong_beam(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_curvature(
    SurfaceCurvatureProxy& f1,
    SurfaceCurvatureProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_curvature(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_displacement(
    SurfaceDisplacementProxy& f1,
    SurfaceDisplacementProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_displacement(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_displacement_pt(
    SurfaceDisplacementPtProxy& f1,
    SurfaceDisplacementPtProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_displacement_pt(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_h_misalign(
    SurfaceHMisalignProxy& f1,
    SurfaceHMisalignProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_h_misalign(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_h_misalign_pt(
    SurfaceHMisalignPtProxy& f1,
    SurfaceHMisalignPtProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_h_misalign_pt(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_segmented(
    SurfaceSegmentedProxy& f1,
    SurfaceSegmentedProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_segmented(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_surface_segmented_pt(
    SurfaceSegmentedPtProxy& f1,
    SurfaceSegmentedPtProxy& f2,
    bool& is_eq) {
  fortran_eq_surface_segmented_pt(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_target_point(
    TargetPointProxy& f1,
    TargetPointProxy& f2,
    bool& is_eq) {
  fortran_eq_target_point(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_taylor(TaylorProxy& f1, TaylorProxy& f2, bool& is_eq) {
  fortran_eq_taylor(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_taylor_term(
    TaylorTermProxy& f1,
    TaylorTermProxy& f2,
    bool& is_eq) {
  fortran_eq_taylor_term(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_track(TrackProxy& f1, TrackProxy& f2, bool& is_eq) {
  fortran_eq_track(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_track_point(
    TrackPointProxy& f1,
    TrackPointProxy& f2,
    bool& is_eq) {
  fortran_eq_track_point(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_twiss(TwissProxy& f1, TwissProxy& f2, bool& is_eq) {
  fortran_eq_twiss(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake(WakeProxy& f1, WakeProxy& f2, bool& is_eq) {
  fortran_eq_wake(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake_lr(WakeLrProxy& f1, WakeLrProxy& f2, bool& is_eq) {
  fortran_eq_wake_lr(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake_lr_mode(
    WakeLrModeProxy& f1,
    WakeLrModeProxy& f2,
    bool& is_eq) {
  fortran_eq_wake_lr_mode(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake_sr(WakeSrProxy& f1, WakeSrProxy& f2, bool& is_eq) {
  fortran_eq_wake_sr(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake_sr_mode(
    WakeSrModeProxy& f1,
    WakeSrModeProxy& f2,
    bool& is_eq) {
  fortran_eq_wake_sr_mode(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wake_sr_z_long(
    WakeSrZLongProxy& f1,
    WakeSrZLongProxy& f2,
    bool& is_eq) {
  fortran_eq_wake_sr_z_long(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wall3d(Wall3dProxy& f1, Wall3dProxy& f2, bool& is_eq) {
  fortran_eq_wall3d(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wall3d_section(
    Wall3dSectionProxy& f1,
    Wall3dSectionProxy& f2,
    bool& is_eq) {
  fortran_eq_wall3d_section(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_wall3d_vertex(
    Wall3dVertexProxy& f1,
    Wall3dVertexProxy& f2,
    bool& is_eq) {
  fortran_eq_wall3d_vertex(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::eq_xy_disp(XyDispProxy& f1, XyDispProxy& f2, bool& is_eq) {
  fortran_eq_xy_disp(
      /* void* */ f1.get_fortran_ptr(),
      /* void* */ f2.get_fortran_ptr(),
      /* bool& */ is_eq);
}
void Bmad::equal_sign_here(EleProxy& ele, std::string& delim, bool& is_here) {
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_equal_sign_here(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ is_here);
}
void Bmad::equivalent_taylor_attributes(
    EleProxy& ele_taylor,
    EleProxy& ele2,
    bool& equiv) {
  fortran_equivalent_taylor_attributes(
      /* void* */ ele_taylor.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool& */ equiv);
}
void Bmad::etdiv(
    double& A,
    double& B,
    double& C,
    double& D,
    double& E,
    double& F) {
  fortran_etdiv(
      /* double& */ A,
      /* double& */ B,
      /* double& */ C,
      /* double& */ D,
      /* double& */ E,
      /* double& */ F);
}
Bmad::EvaluateArrayIndex Bmad::evaluate_array_index(
    std::string delim_list1,
    std::string delim_list2) {
  bool _err_flag{};
  auto _delim_list1 = delim_list1.c_str();
  char _word2[4096];
  auto _delim_list2 = delim_list2.c_str();
  char _delim2[4096];
  int _this_index{};
  fortran_evaluate_array_index(
      /* bool& */ _err_flag,
      /* const char* */ _delim_list1,
      /* const char* */ _word2,
      /* const char* */ _delim_list2,
      /* const char* */ _delim2,
      /* int& */ _this_index);
  return EvaluateArrayIndex{_err_flag, _word2, _delim2, _this_index};
}
Bmad::EvaluateLogical Bmad::evaluate_logical(std::string word) {
  auto _word = word.c_str();
  int _iostat{};
  bool _this_logic{};
  fortran_evaluate_logical(
      /* const char* */ _word, /* int& */ _iostat, /* bool& */ _this_logic);
  return EvaluateLogical{_iostat, _this_logic};
}
void Bmad::exact_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_exact_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
double Bmad::exp_bessi0(double t, double B1, double B2) {
  double _func_retval__{};
  fortran_exp_bessi0(
      /* double& */ t,
      /* double& */ B1,
      /* double& */ B2,
      /* double& */ _func_retval__);
  return _func_retval__;
}
void Bmad::expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string& delim,
    bool& delim_found,
    bool& is_ok) {
  auto _delim_list = delim_list.c_str();
  auto _ele_name = ele_name.c_str();
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_expect_one_of(
      /* const char* */ _delim_list,
      /* bool& */ check_input_delim,
      /* const char* */ _ele_name,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ is_ok);
}
Bmad::ExpectThis Bmad::expect_this(
    std::string expecting,
    bool check_delim,
    bool call_check,
    std::string err_str,
    EleProxy& ele) {
  auto _expecting = expecting.c_str();
  auto _err_str = err_str.c_str();
  char _delim[4096];
  bool _delim_found{};
  bool _is_ok{};
  fortran_expect_this(
      /* const char* */ _expecting,
      /* bool& */ check_delim,
      /* bool& */ call_check,
      /* const char* */ _err_str,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* bool& */ _is_ok);
  return ExpectThis{_delim, _delim_found, _is_ok};
}
std::string Bmad::expression_stack_to_string(
    ExpressionAtomProxyAlloc1D& stack,
    std::optional<bool> polish) {
  // intent=in allocatable type array
  bool polish_lvalue;
  auto* _polish{&polish_lvalue};
  if (polish.has_value()) {
    polish_lvalue = polish.value();
  } else {
    _polish = nullptr;
  }
  char _str[4096];
  fortran_expression_stack_to_string(
      /* void* */ stack.get_fortran_ptr(),
      /* bool* */ _polish,
      /* const char* */ _str);
  return _str;
}
Bmad::ExpressionStackValue Bmad::expression_stack_value(
    ExpressionAtomProxyAlloc1D& stack,
    optional_ref<ControlVar1ProxyAlloc1D> var,
    std::optional<bool> use_old) {
  // intent=in allocatable type array
  bool _err_flag{};
  char _err_str[4096];
  // intent=in allocatable type array
  auto* _var = var.has_value() ? var->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool use_old_lvalue;
  auto* _use_old{&use_old_lvalue};
  if (use_old.has_value()) {
    use_old_lvalue = use_old.value();
  } else {
    _use_old = nullptr;
  }
  double _value{};
  fortran_expression_stack_value(
      /* void* */ stack.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* const char* */ _err_str,
      /* void* */ _var,
      /* bool* */ _use_old,
      /* double& */ _value);
  return ExpressionStackValue{_err_flag, _err_str, _value};
}
Bmad::ExpressionStringToStack Bmad::expression_string_to_stack(
    std::string string) {
  auto _string = string.c_str();
  // intent=out allocatable type array
  auto stack{ExpressionAtomProxyAlloc1D()};
  int _n_stack{};
  bool _err_flag{};
  char _err_str[4096];
  fortran_expression_string_to_stack(
      /* const char* */ _string,
      /* void* */ stack.get_fortran_ptr(),
      /* int& */ _n_stack,
      /* bool& */ _err_flag,
      /* const char* */ _err_str);
  return ExpressionStringToStack{
      std::move(stack), _n_stack, _err_flag, _err_str};
}
Bmad::ExpressionStringToTree Bmad::expression_string_to_tree(
    std::string string,
    ExpressionTreeProxy& root_tree) {
  auto _string = string.c_str();
  bool _err_flag{};
  char _err_str[4096];
  fortran_expression_string_to_tree(
      /* const char* */ _string,
      /* void* */ root_tree.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* const char* */ _err_str);
  return ExpressionStringToTree{_err_flag, _err_str};
}
std::string Bmad::expression_tree_to_string(
    ExpressionTreeProxy& tree,
    std::optional<bool> include_root,
    std::optional<int> n_node,
    optional_ref<ExpressionTreeProxy> parent) {
  bool include_root_lvalue;
  auto* _include_root{&include_root_lvalue};
  if (include_root.has_value()) {
    include_root_lvalue = include_root.value();
  } else {
    _include_root = nullptr;
  }
  int n_node_lvalue;
  auto* _n_node{&n_node_lvalue};
  if (n_node.has_value()) {
    n_node_lvalue = n_node.value();
  } else {
    _n_node = nullptr;
  }
  auto* _parent = parent.has_value() ? parent->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  char _str_out[4096];
  fortran_expression_tree_to_string(
      /* void* */ tree.get_fortran_ptr(),
      /* bool* */ _include_root,
      /* int* */ _n_node,
      /* void* */ _parent,
      /* const char* */ _str_out);
  return _str_out;
}
Bmad::ExpressionValue Bmad::expression_value(
    std::string expression,
    optional_ref<ControlVar1ProxyAlloc1D> var,
    std::optional<bool> use_old) {
  auto _expression = expression.c_str();
  bool _err_flag{};
  char _err_str[4096];
  // intent=in allocatable type array
  auto* _var = var.has_value() ? var->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool use_old_lvalue;
  auto* _use_old{&use_old_lvalue};
  if (use_old.has_value()) {
    use_old_lvalue = use_old.value();
  } else {
    _use_old = nullptr;
  }
  double _value{};
  fortran_expression_value(
      /* const char* */ _expression,
      /* bool& */ _err_flag,
      /* const char* */ _err_str,
      /* void* */ _var,
      /* bool* */ _use_old,
      /* double& */ _value);
  return ExpressionValue{_err_flag, _err_str, _value};
}
int Bmad::fft1(RealAlloc1D& a, RealAlloc1D& b, int n, int isn) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  int _ierr{};
  fortran_fft1(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* int& */ n,
      /* int& */ isn,
      /* int& */ _ierr);
  return _ierr;
}
bool Bmad::field_attribute_free(EleProxy& ele, std::string attrib_name) {
  auto _attrib_name = attrib_name.c_str();
  bool _free{};
  fortran_field_attribute_free(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_name,
      /* bool& */ _free);
  return _free;
}
void Bmad::finalize_reflectivity_table(
    PhotonReflectTableProxy& table,
    bool in_degrees) {
  fortran_finalize_reflectivity_table(
      /* void* */ table.get_fortran_ptr(), /* bool& */ in_degrees);
}
Bmad::FindElementEnds Bmad::find_element_ends(
    EleProxy& ele,
    std::optional<int> ix_multipass) {
  EleProxy _ele1;
  EleProxy _ele2;
  int ix_multipass_lvalue;
  auto* _ix_multipass{&ix_multipass_lvalue};
  if (ix_multipass.has_value()) {
    ix_multipass_lvalue = ix_multipass.value();
  } else {
    _ix_multipass = nullptr;
  }
  fortran_find_element_ends(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _ele1.get_fortran_ptr(),
      /* void* */ _ele2.get_fortran_ptr(),
      /* int* */ _ix_multipass);
  return FindElementEnds{std::move(_ele1), std::move(_ele2)};
}
double Bmad::find_fwhm(double bound, FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  double _fwhm{};
  fortran_find_fwhm(
      /* double& */ bound, /* double* */ _args, /* double& */ _fwhm);
  return _fwhm;
}
Bmad::FindMatchingFieldmap Bmad::find_matching_fieldmap(
    std::string file_name,
    EleProxy& ele,
    int fm_type,
    std::optional<bool> ignore_slaves) {
  auto _file_name = file_name.c_str();
  EleProxy _match_ele;
  int _ix_field{};
  bool ignore_slaves_lvalue;
  auto* _ignore_slaves{&ignore_slaves_lvalue};
  if (ignore_slaves.has_value()) {
    ignore_slaves_lvalue = ignore_slaves.value();
  } else {
    _ignore_slaves = nullptr;
  }
  fortran_find_matching_fieldmap(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ fm_type,
      /* void* */ _match_ele.get_fortran_ptr(),
      /* int& */ _ix_field,
      /* bool* */ _ignore_slaves);
  return FindMatchingFieldmap{std::move(_match_ele), _ix_field};
}
double Bmad::find_normalization(
    double bound,
    double p0,
    FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  double _pnrml{};
  fortran_find_normalization(
      /* double& */ bound,
      /* double& */ p0,
      /* double* */ _args,
      /* double& */ _pnrml);
  return _pnrml;
}
Bmad::FloorAnglesToWMat Bmad::floor_angles_to_w_mat(
    double theta,
    double phi,
    double psi) {
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  FixedArray2D<Real, 3, 3> w_mat_inv;
  double _w_mat_inv_vec[3 * 3];
  fortran_floor_angles_to_w_mat(
      /* double& */ theta,
      /* double& */ phi,
      /* double& */ psi,
      /* double* */ _w_mat_vec,
      /* double* */ _w_mat_inv_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
  vec_to_matrix(_w_mat_inv_vec, w_mat_inv);
  return FloorAnglesToWMat{w_mat, w_mat_inv};
}
Bmad::FloorWMatToAngles Bmad::floor_w_mat_to_angles(
    FixedArray2D<Real, 3, 3> w_mat,
    optional_ref<FloorPositionProxy> floor0) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  double _theta{};
  double _phi{};
  double _psi{};
  auto* _floor0 = floor0.has_value() ? floor0->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  fortran_floor_w_mat_to_angles(
      /* double* */ _w_mat_vec,
      /* double& */ _theta,
      /* double& */ _phi,
      /* double& */ _psi,
      /* void* */ _floor0);
  return FloorWMatToAngles{_theta, _phi, _psi};
}
ComplexTaylorProxy Bmad::form_complex_taylor(
    TaylorProxy& re_taylor,
    TaylorProxy& im_taylor) {
  ComplexTaylorProxy _complex_taylor;
  fortran_form_complex_taylor(
      /* void* */ re_taylor.get_fortran_ptr(),
      /* void* */ im_taylor.get_fortran_ptr(),
      /* void* */ _complex_taylor.get_fortran_ptr());
  return std::move(_complex_taylor);
}
Bmad::FormDigestedBmadFileName Bmad::form_digested_bmad_file_name(
    std::string lat_file,
    std::optional<std::string> use_line) {
  auto _lat_file = lat_file.c_str();
  char _digested_file[4096];
  char _full_lat_file[4096];
  const char* _use_line = use_line.has_value() ? use_line->c_str() : nullptr;
  fortran_form_digested_bmad_file_name(
      /* const char* */ _lat_file,
      /* const char* */ _digested_file,
      /* const char* */ _full_lat_file,
      /* const char* */ _use_line);
  return FormDigestedBmadFileName{_digested_file, _full_lat_file};
}
void Bmad::fringe_here(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    bool& is_here) {
  fortran_fringe_here(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ particle_at,
      /* bool& */ is_here);
}
FixedArray1D<Real, 3> Bmad::g_bend_from_em_field(
    FixedArray1D<Real, 3> b,
    FixedArray1D<Real, 3> e,
    CoordProxy& orbit) {
  auto* _b = b.data(); // CppWrapperGeneralArgument
  auto* _e = e.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _g_bend;
  fortran_g_bend_from_em_field(
      /* double* */ _b,
      /* double* */ _e,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _g_bend.data());
  return _g_bend;
}
Bmad::GBendingStrengthFromEmField Bmad::g_bending_strength_from_em_field(
    EleProxy& ele,
    LatParamProxy& param,
    double s_rel,
    CoordProxy& orbit,
    bool local_ref_frame) {
  FixedArray1D<Real, 3> _g;
  FixedArray2D<Real, 3, 3> dg;
  double _dg_vec[3 * 3];
  fortran_g_bending_strength_from_em_field(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ s_rel,
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ local_ref_frame,
      /* double* */ _g.data(),
      /* double* */ _dg_vec);
  vec_to_matrix(_dg_vec, dg);
  return GBendingStrengthFromEmField{_g, dg};
}
void Bmad::g_integrals_calc(LatProxy& lat) {
  fortran_g_integrals_calc(/* void* */ lat.get_fortran_ptr());
}
void Bmad::gamma_ref(EleProxy& ele, double& gamma) {
  fortran_gamma_ref(/* void* */ ele.get_fortran_ptr(), /* double& */ gamma);
}
EmTaylorProxyArray1D Bmad::gen_grad1_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    int iz) {
  // Output-only type array
  auto em_taylor = EmTaylorProxyArray1D::allocate(3, 1);

  fortran_gen_grad1_to_em_taylor(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ gen_grad.get_fortran_ptr(),
      /* int& */ iz,
      /* void* */ em_taylor.get_fortran_ptr());
  return std::move(std::move(em_taylor));
}
EmTaylorProxyArray1D Bmad::gen_grad_at_s_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    double s_pos) {
  // Output-only type array
  auto em_taylor = EmTaylorProxyArray1D::allocate(3, 1);

  fortran_gen_grad_at_s_to_em_taylor(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ gen_grad.get_fortran_ptr(),
      /* double& */ s_pos,
      /* void* */ em_taylor.get_fortran_ptr());
  return std::move(std::move(em_taylor));
}
void Bmad::gen_grad_field(
    RealAlloc1D& deriv,
    GenGrad1Proxy& gg,
    double& rho,
    double& theta,
    FixedArray1D<Real, 3> field) {
  // intent=inout allocatable general array
  auto* _field = field.data(); // CppWrapperGeneralArgument
  fortran_gen_grad_field(
      /* void* */ deriv.get_fortran_ptr(),
      /* void* */ gg.get_fortran_ptr(),
      /* double& */ rho,
      /* double& */ theta,
      /* double* */ _field);
}
double Bmad::get_bl_from_fwhm(double bound, FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  double _sigma{};
  fortran_get_bl_from_fwhm(
      /* double& */ bound, /* double* */ _args, /* double& */ _sigma);
  return _sigma;
}
void Bmad::get_called_file(
    std::string& delim,
    std::string& call_file,
    bool& err) {
  auto _delim = delim.c_str(); // ptr, inout, required
  auto _call_file = call_file.c_str(); // ptr, inout, required
  fortran_get_called_file(
      /* const char* */ _delim, /* const char* */ _call_file, /* bool& */ err);
}
Bmad::GetEmitFromSigmaMat Bmad::get_emit_from_sigma_mat(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<FixedArray2D<Real, 6, 6>> Nmat) {
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  FixedArray1D<Real, 3> _normal;
  double _Nmat_vec[6 * 6];
  const double* _Nmat = nullptr;
  if (Nmat.has_value()) {
    matrix_to_vec(Nmat.value(), _Nmat_vec);
    _Nmat = _Nmat_vec;
  }
  bool _err_flag{};
  fortran_get_emit_from_sigma_mat(
      /* double* */ _sigma_mat_vec,
      /* double* */ _normal.data(),
      /* double* */ _Nmat_vec,
      /* bool& */ _err_flag);
  return GetEmitFromSigmaMat{_normal, _err_flag};
}
void Bmad::get_next_word(
    std::string word,
    int ix_word,
    std::string delim_list,
    std::string delim,
    bool delim_found,
    std::optional<bool> upper_case_word,
    std::optional<bool> call_check,
    std::optional<bool> err_flag) {
  auto _word = word.c_str();
  auto _delim_list = delim_list.c_str();
  auto _delim = delim.c_str();
  bool upper_case_word_lvalue;
  auto* _upper_case_word{&upper_case_word_lvalue};
  if (upper_case_word.has_value()) {
    upper_case_word_lvalue = upper_case_word.value();
  } else {
    _upper_case_word = nullptr;
  }
  bool call_check_lvalue;
  auto* _call_check{&call_check_lvalue};
  if (call_check.has_value()) {
    call_check_lvalue = call_check.value();
  } else {
    _call_check = nullptr;
  }
  bool err_flag_lvalue;
  auto* _err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_get_next_word(
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim_list,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool* */ _upper_case_word,
      /* bool* */ _call_check,
      /* bool* */ _err_flag);
}
Bmad::GetSlaveList Bmad::get_slave_list(EleProxy& lord) {
  // intent=out allocatable type array
  auto slaves{ElePointerProxyAlloc1D()};
  int _n_slave{};
  fortran_get_slave_list(
      /* void* */ lord.get_fortran_ptr(),
      /* void* */ slaves.get_fortran_ptr(),
      /* int& */ _n_slave);
  return GetSlaveList{std::move(slaves), _n_slave};
}
void Bmad::gpt_field_grid_scaling(
    EleProxy& ele,
    int& dimensions,
    double& field_scale,
    double& ref_time) {
  fortran_gpt_field_grid_scaling(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ dimensions,
      /* double& */ field_scale,
      /* double& */ ref_time);
}
void Bmad::gpt_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double& field_value) {
  fortran_gpt_max_field_reference(
      /* void* */ pt0.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ field_value);
}
Bmad::GptToParticleBunch Bmad::gpt_to_particle_bunch(
    std::string gpt_file,
    EleProxy& ele) {
  auto _gpt_file = gpt_file.c_str();
  BunchProxy _bunch;
  bool _err_flag{};
  fortran_gpt_to_particle_bunch(
      /* const char* */ _gpt_file,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _bunch.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return GptToParticleBunch{std::move(_bunch), _err_flag};
}
void Bmad::gradient_shift_sr_wake(
    EleProxy& ele,
    LatParamProxy& param,
    double& grad_shift) {
  fortran_gradient_shift_sr_wake(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ grad_shift);
}
GridFieldPt1Proxy Bmad::grid_field_interpolate(
    EleProxy& ele,
    CoordProxy& orbit,
    GridFieldProxy& grid,
    bool err_flag,
    double x1,
    std::optional<double> x2,
    std::optional<double> x3,
    std::optional<bool> allow_s_out_of_bounds,
    std::optional<bool> print_err) {
  GridFieldPt1Proxy _g_field;
  double x2_lvalue;
  auto* _x2{&x2_lvalue};
  if (x2.has_value()) {
    x2_lvalue = x2.value();
  } else {
    _x2 = nullptr;
  }
  double x3_lvalue;
  auto* _x3{&x3_lvalue};
  if (x3.has_value()) {
    x3_lvalue = x3.value();
  } else {
    _x3 = nullptr;
  }
  bool allow_s_out_of_bounds_lvalue;
  auto* _allow_s_out_of_bounds{&allow_s_out_of_bounds_lvalue};
  if (allow_s_out_of_bounds.has_value()) {
    allow_s_out_of_bounds_lvalue = allow_s_out_of_bounds.value();
  } else {
    _allow_s_out_of_bounds = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_grid_field_interpolate(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ grid.get_fortran_ptr(),
      /* void* */ _g_field.get_fortran_ptr(),
      /* bool& */ err_flag,
      /* double& */ x1,
      /* double* */ _x2,
      /* double* */ _x3,
      /* bool* */ _allow_s_out_of_bounds,
      /* bool* */ _print_err);
  return std::move(_g_field);
}
void Bmad::hard_multipole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_hard_multipole_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::has_attribute(EleProxy& ele, std::string& attrib, bool& has_it) {
  auto _attrib = attrib.c_str(); // ptr, inout, required
  fortran_has_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib,
      /* bool& */ has_it);
}
bool Bmad::has_curvature(PhotonElementProxy& phot_ele) {
  bool _curved{};
  fortran_has_curvature(
      /* void* */ phot_ele.get_fortran_ptr(), /* bool& */ _curved);
  return _curved;
}
bool Bmad::has_orientation_attributes(EleProxy& ele) {
  bool _has_attribs{};
  fortran_has_orientation_attributes(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _has_attribs);
  return _has_attribs;
}
void Bmad::hdf5_write_beam(
    std::string& file_name,
    BunchProxyAlloc1D& bunches,
    bool& append,
    bool& error,
    optional_ref<LatProxy> lat,
    optional_ref<bool> alive_only) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  // intent=inout allocatable type array
  auto* _lat = lat.has_value() ? lat->get().get_fortran_ptr()
                               : nullptr; // input, optional
  auto* _alive_only =
      alive_only.has_value() ? &alive_only->get() : nullptr; // inout, optional
  fortran_hdf5_write_beam(
      /* const char* */ _file_name,
      /* void* */ bunches.get_fortran_ptr(),
      /* bool& */ append,
      /* bool& */ error,
      /* void* */ _lat,
      /* bool* */ _alive_only);
}
void Bmad::hdf5_write_grid_field(
    std::string& file_name,
    EleProxy& ele,
    GridFieldProxyAlloc1D& g_field,
    bool& err_flag) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  // intent=inout allocatable type array
  fortran_hdf5_write_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ g_field.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::hwang_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_hwang_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool& tail_cut,
    double& tau,
    double& energy,
    double& n_part,
    int& species,
    FixedArray2D<Real, 6, 6> ibs_mat) {
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  double _ibs_mat_vec[6 * 6];
  matrix_to_vec(ibs_mat, _ibs_mat_vec);
  fortran_ibs_matrix_c(
      /* double* */ _sigma_mat_vec,
      /* bool& */ tail_cut,
      /* double& */ tau,
      /* double& */ energy,
      /* double& */ n_part,
      /* int& */ species,
      /* double* */ _ibs_mat_vec);
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  vec_to_matrix(_ibs_mat_vec, ibs_mat);
}
void Bmad::igfcoulombfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res) {
  fortran_igfcoulombfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res);
}
void Bmad::igfexfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res) {
  fortran_igfexfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res);
}
void Bmad::igfeyfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res) {
  fortran_igfeyfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res);
}
void Bmad::igfezfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res) {
  fortran_igfezfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res);
}
void Bmad::init_attribute_name1(
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state,
    std::optional<bool> override) {
  auto _name = name.c_str();
  int attrib_state_lvalue;
  auto* _attrib_state{&attrib_state_lvalue};
  if (attrib_state.has_value()) {
    attrib_state_lvalue = attrib_state.value();
  } else {
    _attrib_state = nullptr;
  }
  bool override_lvalue;
  auto* _override{&override_lvalue};
  if (override.has_value()) {
    override_lvalue = override.value();
  } else {
    _override = nullptr;
  }
  fortran_init_attribute_name1(
      /* int& */ ix_key,
      /* int& */ ix_attrib,
      /* const char* */ _name,
      /* int* */ _attrib_state,
      /* bool* */ _override);
}
void Bmad::init_attribute_name_array() {
  fortran_init_attribute_name_array();
}
Bmad::InitBeamDistribution Bmad::init_beam_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    optional_ref<NormalModesProxy> modes,
    std::optional<bool> print_p0c_shift_warning,
    optional_ref<bool> conserve_momentum) {
  BeamProxy _beam;
  bool _err_flag{};
  auto* _modes = modes.has_value() ? modes->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  BeamInitProxy _beam_init_set;
  bool print_p0c_shift_warning_lvalue;
  auto* _print_p0c_shift_warning{&print_p0c_shift_warning_lvalue};
  if (print_p0c_shift_warning.has_value()) {
    print_p0c_shift_warning_lvalue = print_p0c_shift_warning.value();
  } else {
    _print_p0c_shift_warning = nullptr;
  }
  auto* _conserve_momentum = conserve_momentum.has_value()
      ? &conserve_momentum->get()
      : nullptr; // inout, optional
  fortran_init_beam_distribution(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* void* */ _beam.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _modes,
      /* void* */ _beam_init_set.get_fortran_ptr(),
      /* bool* */ _print_p0c_shift_warning,
      /* bool* */ _conserve_momentum);
  return InitBeamDistribution{
      std::move(_beam), _err_flag, std::move(_beam_init_set)};
}
void Bmad::init_bmad() {
  fortran_init_bmad();
}
void Bmad::init_bmad_parser_common(optional_ref<LatProxy> lat) {
  auto* _lat = lat.has_value() ? lat->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_init_bmad_parser_common(/* void* */ _lat);
}
Bmad::InitBunchDistribution Bmad::init_bunch_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    int ix_bunch,
    optional_ref<NormalModesProxy> modes,
    std::optional<bool> print_p0c_shift_warning,
    optional_ref<bool> conserve_momentum) {
  BunchProxy _bunch;
  bool _err_flag{};
  auto* _modes = modes.has_value() ? modes->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  BeamInitProxy _beam_init_used;
  bool print_p0c_shift_warning_lvalue;
  auto* _print_p0c_shift_warning{&print_p0c_shift_warning_lvalue};
  if (print_p0c_shift_warning.has_value()) {
    print_p0c_shift_warning_lvalue = print_p0c_shift_warning.value();
  } else {
    _print_p0c_shift_warning = nullptr;
  }
  auto* _conserve_momentum = conserve_momentum.has_value()
      ? &conserve_momentum->get()
      : nullptr; // inout, optional
  fortran_init_bunch_distribution(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* int& */ ix_bunch,
      /* void* */ _bunch.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _modes,
      /* void* */ _beam_init_used.get_fortran_ptr(),
      /* bool* */ _print_p0c_shift_warning,
      /* bool* */ _conserve_momentum);
  return InitBunchDistribution{
      std::move(_bunch), _err_flag, std::move(_beam_init_used)};
}
void Bmad::init_complex_taylor_series(
    ComplexTaylorProxy& complex_taylor,
    int n_term,
    std::optional<bool> save) {
  bool save_lvalue;
  auto* _save{&save_lvalue};
  if (save.has_value()) {
    save_lvalue = save.value();
  } else {
    _save = nullptr;
  }
  fortran_init_complex_taylor_series(
      /* void* */ complex_taylor.get_fortran_ptr(),
      /* int& */ n_term,
      /* bool* */ _save);
}
void Bmad::init_coord(
    CoordProxy& orb,
    FixedArray1D<Real, 6> vec,
    optional_ref<EleProxy> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin,
    std::optional<double> s_pos,
    std::optional<bool> random_on) {
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  int element_end_lvalue;
  auto* _element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto* _particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto* _E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto* _t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto* _shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  double* _spin = spin.has_value() ? spin.value().data() : nullptr;
  double s_pos_lvalue;
  auto* _s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  bool random_on_lvalue;
  auto* _random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_coord1(
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _vec,
      /* void* */ _ele,
      /* int* */ _element_end,
      /* int* */ _particle,
      /* int* */ _direction,
      /* double* */ _E_photon,
      /* double* */ _t_offset,
      /* bool* */ _shift_vec6,
      /* double* */ _spin,
      /* double* */ _s_pos,
      /* bool* */ _random_on);
}
CoordProxy Bmad::init_coord(
    CoordProxy& orb_in,
    optional_ref<EleProxy> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin,
    std::optional<double> s_pos,
    std::optional<bool> random_on) {
  CoordProxy _orb_out;
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  int element_end_lvalue;
  auto* _element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto* _particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto* _E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto* _t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto* _shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  double* _spin = spin.has_value() ? spin.value().data() : nullptr;
  double s_pos_lvalue;
  auto* _s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  bool random_on_lvalue;
  auto* _random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_coord2(
      /* void* */ _orb_out.get_fortran_ptr(),
      /* void* */ orb_in.get_fortran_ptr(),
      /* void* */ _ele,
      /* int* */ _element_end,
      /* int* */ _particle,
      /* int* */ _direction,
      /* double* */ _E_photon,
      /* double* */ _t_offset,
      /* bool* */ _shift_vec6,
      /* double* */ _spin,
      /* double* */ _s_pos,
      /* bool* */ _random_on);
  return std::move(_orb_out);
}
void Bmad::init_coord(
    CoordProxy& orb,
    optional_ref<EleProxy> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin) {
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  int element_end_lvalue;
  auto* _element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto* _particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto* _E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto* _t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto* _shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  double* _spin = spin.has_value() ? spin.value().data() : nullptr;
  fortran_init_coord3(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ _ele,
      /* int* */ _element_end,
      /* int* */ _particle,
      /* int* */ _direction,
      /* double* */ _E_photon,
      /* double* */ _t_offset,
      /* bool* */ _shift_vec6,
      /* double* */ _spin);
}
void Bmad::init_custom(LatProxy& lat) {
  fortran_init_custom(/* void* */ lat.get_fortran_ptr());
}
EleProxy Bmad::init_ele(
    std::optional<int> key,
    std::optional<int> sub_key,
    std::optional<int> ix_ele,
    optional_ref<BranchProxy> branch) {
  EleProxy _ele;
  int key_lvalue;
  auto* _key{&key_lvalue};
  if (key.has_value()) {
    key_lvalue = key.value();
  } else {
    _key = nullptr;
  }
  int sub_key_lvalue;
  auto* _sub_key{&sub_key_lvalue};
  if (sub_key.has_value()) {
    sub_key_lvalue = sub_key.value();
  } else {
    _sub_key = nullptr;
  }
  int ix_ele_lvalue;
  auto* _ix_ele{&ix_ele_lvalue};
  if (ix_ele.has_value()) {
    ix_ele_lvalue = ix_ele.value();
  } else {
    _ix_ele = nullptr;
  }
  auto* _branch = branch.has_value() ? branch->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  fortran_init_ele(
      /* void* */ _ele.get_fortran_ptr(),
      /* int* */ _key,
      /* int* */ _sub_key,
      /* int* */ _ix_ele,
      /* void* */ _branch);
  return std::move(_ele);
}
void Bmad::init_em_taylor_series(
    EmTaylorProxy& em_taylor,
    int n_term,
    std::optional<bool> save_old) {
  bool save_old_lvalue;
  auto* _save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  fortran_init_em_taylor_series(
      /* void* */ em_taylor.get_fortran_ptr(),
      /* int& */ n_term,
      /* bool* */ _save_old);
}
LatProxy Bmad::init_lat(
    std::optional<int> n,
    std::optional<bool> init_beginning_ele) {
  LatProxy _lat;
  int n_lvalue;
  auto* _n{&n_lvalue};
  if (n.has_value()) {
    n_lvalue = n.value();
  } else {
    _n = nullptr;
  }
  bool init_beginning_ele_lvalue;
  auto* _init_beginning_ele{&init_beginning_ele_lvalue};
  if (init_beginning_ele.has_value()) {
    init_beginning_ele_lvalue = init_beginning_ele.value();
  } else {
    _init_beginning_ele = nullptr;
  }
  fortran_init_lat(
      /* void* */ _lat.get_fortran_ptr(),
      /* int* */ _n,
      /* bool* */ _init_beginning_ele);
  return std::move(_lat);
}
void Bmad::init_multipole_cache(EleProxy& ele) {
  fortran_init_multipole_cache(/* void* */ ele.get_fortran_ptr());
}
CoordProxy Bmad::init_photon_from_a_photon_init_ele(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> random_on) {
  CoordProxy _orbit;
  bool random_on_lvalue;
  auto* _random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_photon_from_a_photon_init_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _orbit.get_fortran_ptr(),
      /* bool* */ _random_on);
  return std::move(_orbit);
}
Bmad::InitPhotonIntegProb Bmad::init_photon_integ_prob(
    double gamma,
    double g,
    double E_min,
    double E_max,
    std::optional<double> vert_angle_min,
    std::optional<double> vert_angle_max,
    std::optional<bool> vert_angle_symmetric,
    std::optional<double> energy_integ_prob) {
  double vert_angle_min_lvalue;
  auto* _vert_angle_min{&vert_angle_min_lvalue};
  if (vert_angle_min.has_value()) {
    vert_angle_min_lvalue = vert_angle_min.value();
  } else {
    _vert_angle_min = nullptr;
  }
  double vert_angle_max_lvalue;
  auto* _vert_angle_max{&vert_angle_max_lvalue};
  if (vert_angle_max.has_value()) {
    vert_angle_max_lvalue = vert_angle_max.value();
  } else {
    _vert_angle_max = nullptr;
  }
  bool vert_angle_symmetric_lvalue;
  auto* _vert_angle_symmetric{&vert_angle_symmetric_lvalue};
  if (vert_angle_symmetric.has_value()) {
    vert_angle_symmetric_lvalue = vert_angle_symmetric.value();
  } else {
    _vert_angle_symmetric = nullptr;
  }
  double energy_integ_prob_lvalue;
  auto* _energy_integ_prob{&energy_integ_prob_lvalue};
  if (energy_integ_prob.has_value()) {
    energy_integ_prob_lvalue = energy_integ_prob.value();
  } else {
    _energy_integ_prob = nullptr;
  }
  double _E_photon{};
  double _integ_prob{};
  fortran_init_photon_integ_prob(
      /* double& */ gamma,
      /* double& */ g,
      /* double& */ E_min,
      /* double& */ E_max,
      /* double* */ _vert_angle_min,
      /* double* */ _vert_angle_max,
      /* bool* */ _vert_angle_symmetric,
      /* double* */ _energy_integ_prob,
      /* double& */ _E_photon,
      /* double& */ _integ_prob);
  return InitPhotonIntegProb{_E_photon, _integ_prob};
}
BunchProxy Bmad::init_spin_distribution(
    BeamInitProxy& beam_init,
    EleProxy& ele) {
  BunchProxy _bunch;
  fortran_init_spin_distribution(
      /* void* */ beam_init.get_fortran_ptr(),
      /* void* */ _bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr());
  return std::move(_bunch);
}
void Bmad::init_surface_segment(PhotonElementProxy& phot, int& ix, int& iy) {
  fortran_init_surface_segment(
      /* void* */ phot.get_fortran_ptr(), /* int& */ ix, /* int& */ iy);
}
void Bmad::init_taylor_series(
    TaylorProxy& bmad_taylor,
    int n_term,
    std::optional<bool> save_old) {
  bool save_old_lvalue;
  auto* _save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  fortran_init_taylor_series(
      /* void* */ bmad_taylor.get_fortran_ptr(),
      /* int& */ n_term,
      /* bool* */ _save_old);
}
WakeProxy Bmad::init_wake(
    int n_sr_long,
    int n_sr_trans,
    int n_sr_z,
    int n_lr_mode,
    std::optional<bool> always_allocate) {
  WakeProxy _wake;
  bool always_allocate_lvalue;
  auto* _always_allocate{&always_allocate_lvalue};
  if (always_allocate.has_value()) {
    always_allocate_lvalue = always_allocate.value();
  } else {
    _always_allocate = nullptr;
  }
  fortran_init_wake(
      /* void* */ _wake.get_fortran_ptr(),
      /* int& */ n_sr_long,
      /* int& */ n_sr_trans,
      /* int& */ n_sr_z,
      /* int& */ n_lr_mode,
      /* bool* */ _always_allocate);
  return std::move(_wake);
}
void Bmad::insert_element(
    LatProxy& lat,
    EleProxy& insert_ele,
    int ix_ele,
    std::optional<int> ix_branch,
    optional_ref<CoordProxyAlloc1D> orbit) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  // intent=inout allocatable type array
  auto* _orbit = orbit.has_value() ? orbit->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  fortran_insert_element(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ insert_ele.get_fortran_ptr(),
      /* int& */ ix_ele,
      /* int* */ _ix_branch,
      /* void* */ _orbit);
}
void Bmad::integrand_base(double t, RealAlloc1D& args, double& func_retval__) {
  // intent=inout allocatable general array
  fortran_integrand_base(
      /* double& */ t,
      /* void* */ args.get_fortran_ptr(),
      /* double& */ func_retval__);
}
double Bmad::integrate_psi(
    double bound,
    double p0,
    FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  double _result{};
  fortran_integrate_psi(
      /* double& */ bound,
      /* double& */ p0,
      /* double* */ _args,
      /* double& */ _result);
  return _result;
}
void Bmad::integrated_mats(
    EleProxyAlloc1D& eles,
    CoordProxyAlloc1D& coos,
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    NormalModesProxy& mode) {
  // intent=inout allocatable type array
  // intent=inout allocatable type array
  std::complex<double> _Lambda_vec[6 * 6];
  matrix_to_vec(Lambda, _Lambda_vec);
  std::complex<double> _Theta_vec[6 * 6];
  matrix_to_vec(Theta, _Theta_vec);
  std::complex<double> _Iota_vec[6 * 6];
  matrix_to_vec(Iota, _Iota_vec);
  fortran_integrated_mats(
      /* void* */ eles.get_fortran_ptr(),
      /* void* */ coos.get_fortran_ptr(),
      /* std::complex<double>* */ _Lambda_vec,
      /* std::complex<double>* */ _Theta_vec,
      /* std::complex<double>* */ _Iota_vec,
      /* void* */ mode.get_fortran_ptr());
  vec_to_matrix(_Lambda_vec, Lambda);
  vec_to_matrix(_Theta_vec, Theta);
  vec_to_matrix(_Iota_vec, Iota);
}
void Bmad::integration_timer(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start,
    CoordProxy& orb_max,
    double& tol) {
  fortran_integration_timer_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ start.get_fortran_ptr(),
      /* void* */ orb_max.get_fortran_ptr(),
      /* double& */ tol);
}
FixedArray1D<Real, 3> Bmad::ion_kick(
    CoordProxy& orbit,
    FixedArray1D<Real, 2> r_beam,
    double n_beam_part,
    TwissProxy& a_twiss,
    TwissProxy& b_twiss,
    double sig_ee) {
  auto* _r_beam = r_beam.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _kick;
  fortran_ion_kick(
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _r_beam,
      /* double& */ n_beam_part,
      /* void* */ a_twiss.get_fortran_ptr(),
      /* void* */ b_twiss.get_fortran_ptr(),
      /* double& */ sig_ee,
      /* double* */ _kick.data());
  return _kick;
}
bool Bmad::is_attribute(int ix_attrib, int which) {
  bool _is_attrib{};
  fortran_is_attribute(
      /* int& */ ix_attrib, /* int& */ which, /* bool& */ _is_attrib);
  return _is_attrib;
}
void Bmad::key_name_to_key_index(
    std::string key_str,
    std::optional<bool> abbrev_allowed,
    int& key_index) {
  auto _key_str = key_str.c_str();
  bool abbrev_allowed_lvalue;
  auto* _abbrev_allowed{&abbrev_allowed_lvalue};
  if (abbrev_allowed.has_value()) {
    abbrev_allowed_lvalue = abbrev_allowed.value();
  } else {
    _abbrev_allowed = nullptr;
  }
  fortran_key_name_to_key_index(
      /* const char* */ _key_str,
      /* bool* */ _abbrev_allowed,
      /* int& */ key_index);
}
Bmad::KickVectorCalc Bmad::kick_vector_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_body,
    CoordProxy& orbit,
    optional_ref<bool> print_err) {
  FixedArray1D<Real, 11> _dr_ds;
  bool _err{};
  auto* _print_err =
      print_err.has_value() ? &print_err->get() : nullptr; // inout, optional
  fortran_kick_vector_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ s_body,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _dr_ds.data(),
      /* bool& */ _err,
      /* bool* */ _print_err);
  return KickVectorCalc{_dr_ds, _err};
}
void Bmad::kill_complex_taylor(ComplexTaylorProxyAlloc1D& complex_taylor) {
  // intent=inout allocatable type array
  fortran_kill_complex_taylor(/* void* */ complex_taylor.get_fortran_ptr());
}
void Bmad::kill_ptc_layouts(LatProxy& lat) {
  fortran_kill_ptc_layouts(/* void* */ lat.get_fortran_ptr());
}
void Bmad::kill_taylor(TaylorProxyAlloc1D& bmad_taylor) {
  // intent=inout allocatable type array
  fortran_kill_taylor(/* void* */ bmad_taylor.get_fortran_ptr());
}
std::string Bmad::kind_name(int this_kind) {
  char _kind_str[4096];
  fortran_kind_name(/* int* */ &this_kind, /* const char* */ _kind_str);
  return _kind_str;
}
bool Bmad::knot_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x_pt,
    int interpolation,
    double& y_pt) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  bool _err_flag{};
  fortran_knot_interpolate(
      /* void* */ x_knot.get_fortran_ptr(),
      /* void* */ y_knot.get_fortran_ptr(),
      /* double& */ x_pt,
      /* int& */ interpolation,
      /* bool& */ _err_flag,
      /* double& */ y_pt);
  return _err_flag;
}
void Bmad::knots_to_string(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    std::string& str) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  auto _str = str.c_str(); // ptr, inout, required
  fortran_knots_to_string(
      /* void* */ x_knot.get_fortran_ptr(),
      /* void* */ y_knot.get_fortran_ptr(),
      /* const char* */ _str);
}
void Bmad::lafun(double& x, double& y, double& z, double& res) {
  fortran_lafun(
      /* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
bool Bmad::lat_compute_ref_energy_and_time(LatProxy& lat) {
  bool _err_flag{};
  fortran_lat_compute_ref_energy_and_time(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::lat_ele_locator(
    std::string loc_str,
    LatProxy& lat,
    ElePointerProxyAlloc1D& eles,
    int& n_loc,
    std::optional<bool> above_ubound_is_err,
    std::optional<int> ix_dflt_branch,
    std::optional<bool> order_by_index,
    std::optional<bool> append_eles) {
  auto _loc_str = loc_str.c_str();
  // intent=inout allocatable type array
  bool _err{};
  bool above_ubound_is_err_lvalue;
  auto* _above_ubound_is_err{&above_ubound_is_err_lvalue};
  if (above_ubound_is_err.has_value()) {
    above_ubound_is_err_lvalue = above_ubound_is_err.value();
  } else {
    _above_ubound_is_err = nullptr;
  }
  int ix_dflt_branch_lvalue;
  auto* _ix_dflt_branch{&ix_dflt_branch_lvalue};
  if (ix_dflt_branch.has_value()) {
    ix_dflt_branch_lvalue = ix_dflt_branch.value();
  } else {
    _ix_dflt_branch = nullptr;
  }
  bool order_by_index_lvalue;
  auto* _order_by_index{&order_by_index_lvalue};
  if (order_by_index.has_value()) {
    order_by_index_lvalue = order_by_index.value();
  } else {
    _order_by_index = nullptr;
  }
  bool append_eles_lvalue;
  auto* _append_eles{&append_eles_lvalue};
  if (append_eles.has_value()) {
    append_eles_lvalue = append_eles.value();
  } else {
    _append_eles = nullptr;
  }
  fortran_lat_ele_locator(
      /* const char* */ _loc_str,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ eles.get_fortran_ptr(),
      /* int& */ n_loc,
      /* bool& */ _err,
      /* bool* */ _above_ubound_is_err,
      /* int* */ _ix_dflt_branch,
      /* bool* */ _order_by_index,
      /* bool* */ _append_eles);
  return _err;
}
void Bmad::lat_equal_lat(LatProxy& lat_out, LatProxy& lat_in) {
  fortran_lat_equal_lat(
      /* void* */ lat_out.get_fortran_ptr(),
      /* void* */ lat_in.get_fortran_ptr());
}
void Bmad::lat_geometry(LatProxy& lat) {
  fortran_lat_geometry(/* void* */ lat.get_fortran_ptr());
}
bool Bmad::lat_make_mat6(
    LatProxy& lat,
    std::optional<int> ix_ele,
    optional_ref<CoordProxyAlloc1D> ref_orb,
    std::optional<int> ix_branch) {
  int ix_ele_lvalue;
  auto* _ix_ele{&ix_ele_lvalue};
  if (ix_ele.has_value()) {
    ix_ele_lvalue = ix_ele.value();
  } else {
    _ix_ele = nullptr;
  }
  // intent=in allocatable type array
  auto* _ref_orb = ref_orb.has_value() ? ref_orb->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  fortran_lat_make_mat6(
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _ix_ele,
      /* void* */ _ref_orb,
      /* int* */ _ix_branch,
      /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::lat_sanity_check(LatProxy& lat) {
  bool _err_flag{};
  fortran_lat_sanity_check(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::lat_to_ptc_layout(LatProxy& lat) {
  fortran_lat_to_ptc_layout(/* void* */ lat.get_fortran_ptr());
}
void Bmad::lat_vec_equal_lat_vec(LatProxyAlloc1D& lat1, LatProxyAlloc1D& lat2) {
  // intent=inout allocatable type array
  // intent=in allocatable type array
  fortran_lat_vec_equal_lat_vec(
      /* void* */ lat1.get_fortran_ptr(), /* void* */ lat2.get_fortran_ptr());
}
bool Bmad::lattice_bookkeeper(LatProxy& lat) {
  bool _err_flag{};
  fortran_lattice_bookkeeper(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::lcavity_rf_step_setup(EleProxy& ele) {
  fortran_lcavity_rf_step_setup(/* void* */ ele.get_fortran_ptr());
}
void Bmad::linear_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_linear_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
Bmad::LinearCoef Bmad::linear_coef(ExpressionAtomProxyAlloc1D& stack) {
  // intent=in allocatable type array
  bool _err_flag{};
  double _coef{};
  fortran_linear_coef(
      /* void* */ stack.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* double& */ _coef);
  return LinearCoef{_err_flag, _coef};
}
TaylorProxyArray1D Bmad::linear_to_spin_taylor(FixedArray2D<Real, 4, 7> q_map) {
  double _q_map_vec[4 * 7];
  matrix_to_vec(q_map, _q_map_vec);
  // Output-only type array
  auto spin_taylor = TaylorProxyArray1D::allocate(4, 1);

  fortran_linear_to_spin_taylor(
      /* double* */ _q_map_vec, /* void* */ spin_taylor.get_fortran_ptr());
  return std::move(std::move(spin_taylor));
}
Bmad::LoadParseLine Bmad::load_parse_line(std::string action, int ix_start) {
  auto _action = action.c_str();
  bool _end_of_file{};
  bool _err_flag{};
  fortran_load_parse_line(
      /* const char* */ _action,
      /* int& */ ix_start,
      /* bool& */ _end_of_file,
      /* bool& */ _err_flag);
  return LoadParseLine{_end_of_file, _err_flag};
}
void Bmad::lord_edge_aligned(
    EleProxy& slave,
    int slave_edge,
    EleProxy& lord,
    bool& is_aligned) {
  fortran_lord_edge_aligned(
      /* void* */ slave.get_fortran_ptr(),
      /* int& */ slave_edge,
      /* void* */ lord.get_fortran_ptr(),
      /* bool& */ is_aligned);
}
void Bmad::low_energy_z_correction(
    CoordProxy& orbit,
    EleProxy& ele,
    double ds,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    double& dz) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_low_energy_z_correction(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ ds,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* double& */ dz);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
MadMapProxy Bmad::mad_add_offsets_and_multipoles(EleProxy& ele) {
  MadMapProxy _map;
  fortran_mad_add_offsets_and_multipoles(
      /* void* */ ele.get_fortran_ptr(), /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_concat_map2(MadMapProxy& map1, MadMapProxy& map2) {
  MadMapProxy _map3;
  fortran_mad_concat_map2(
      /* void* */ map1.get_fortran_ptr(),
      /* void* */ map2.get_fortran_ptr(),
      /* void* */ _map3.get_fortran_ptr());
  return std::move(_map3);
}
MadMapProxy Bmad::mad_drift(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_drift(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_elsep(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_elsep(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
TaylorProxyAlloc1D Bmad::mad_map_to_taylor(
    MadMapProxy& map,
    MadEnergyProxy& energy) {
  // intent=out allocatable type array
  auto taylor{TaylorProxyAlloc1D()};
  fortran_mad_map_to_taylor(
      /* void* */ map.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ taylor.get_fortran_ptr());
  return std::move(taylor);
}
MadMapProxy Bmad::mad_quadrupole(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_quadrupole(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_rfcavity(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_rfcavity(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_sbend(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_sbend(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_sbend_body(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_sbend_body(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_sbend_fringe(
    EleProxy& ele,
    MadEnergyProxy& energy,
    bool into) {
  MadMapProxy _map;
  fortran_mad_sbend_fringe(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* bool& */ into,
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_sextupole(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_sextupole(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapProxy Bmad::mad_solenoid(EleProxy& ele, MadEnergyProxy& energy) {
  MadMapProxy _map;
  fortran_mad_solenoid(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
Bmad::MadTmfoc Bmad::mad_tmfoc(double el, double sk1) {
  double _c{};
  double _s{};
  double _d{};
  double _f{};
  fortran_mad_tmfoc(
      /* double& */ el,
      /* double& */ sk1,
      /* double& */ _c,
      /* double& */ _s,
      /* double& */ _d,
      /* double& */ _f);
  return MadTmfoc{_c, _s, _d, _f};
}
void Bmad::mad_tmsymm(FixedArray3D<Real, 6, 6, 6> te) {
  double _te_vec[6 * 6 * 6];
  tensor_to_vec(te, _te_vec);
  fortran_mad_tmsymm(/* double* */ _te_vec);
  vec_to_tensor(_te_vec, te);
}
void Bmad::mad_tmtilt(MadMapProxy& map, double tilt) {
  fortran_mad_tmtilt(/* void* */ map.get_fortran_ptr(), /* double& */ tilt);
}
CoordProxy Bmad::mad_track1(CoordProxy& c0, MadMapProxy& map) {
  CoordProxy _c1;
  fortran_mad_track1(
      /* void* */ c0.get_fortran_ptr(),
      /* void* */ map.get_fortran_ptr(),
      /* void* */ _c1.get_fortran_ptr());
  return std::move(_c1);
}
void Bmad::make_g2_mats(
    TwissProxy& twiss,
    FixedArray2D<Real, 2, 2> g2_mat,
    FixedArray2D<Real, 2, 2> g2_inv_mat) {
  double _g2_mat_vec[2 * 2];
  matrix_to_vec(g2_mat, _g2_mat_vec);
  double _g2_inv_mat_vec[2 * 2];
  matrix_to_vec(g2_inv_mat, _g2_inv_mat_vec);
  fortran_make_g2_mats(
      /* void* */ twiss.get_fortran_ptr(),
      /* double* */ _g2_mat_vec,
      /* double* */ _g2_inv_mat_vec);
  vec_to_matrix(_g2_mat_vec, g2_mat);
  vec_to_matrix(_g2_inv_mat_vec, g2_inv_mat);
}
Bmad::MakeGMats Bmad::make_g_mats(EleProxy& ele) {
  FixedArray2D<Real, 4, 4> g_mat;
  double _g_mat_vec[4 * 4];
  FixedArray2D<Real, 4, 4> g_inv_mat;
  double _g_inv_mat_vec[4 * 4];
  fortran_make_g_mats(
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _g_mat_vec,
      /* double* */ _g_inv_mat_vec);
  vec_to_matrix(_g_mat_vec, g_mat);
  vec_to_matrix(_g_inv_mat_vec, g_inv_mat);
  return MakeGMats{g_mat, g_inv_mat};
}
Bmad::MakeHvbp Bmad::make_hvbp(FixedArray2D<Real, 6, 6> N) {
  double _N_vec[6 * 6];
  matrix_to_vec(N, _N_vec);
  FixedArray2D<Real, 6, 6> B;
  double _B_vec[6 * 6];
  FixedArray2D<Real, 6, 6> V;
  double _V_vec[6 * 6];
  FixedArray2D<Real, 6, 6> H;
  double _H_vec[6 * 6];
  FixedArray2D<Real, 6, 6> Vbar;
  double _Vbar_vec[6 * 6];
  FixedArray2D<Real, 6, 6> Hbar;
  double _Hbar_vec[6 * 6];
  fortran_make_hvbp(
      /* double* */ _N_vec,
      /* double* */ _B_vec,
      /* double* */ _V_vec,
      /* double* */ _H_vec,
      /* double* */ _Vbar_vec,
      /* double* */ _Hbar_vec);
  vec_to_matrix(_B_vec, B);
  vec_to_matrix(_V_vec, V);
  vec_to_matrix(_H_vec, H);
  vec_to_matrix(_Vbar_vec, Vbar);
  vec_to_matrix(_Hbar_vec, Hbar);
  return MakeHvbp{B, V, H, Vbar, Hbar};
}
LatProxy Bmad::make_hybrid_lat(
    LatProxy& lat_in,
    std::optional<bool> use_taylor,
    optional_ref<CoordArrayProxyAlloc1D> orb0_arr) {
  LatProxy _lat_out;
  bool use_taylor_lvalue;
  auto* _use_taylor{&use_taylor_lvalue};
  if (use_taylor.has_value()) {
    use_taylor_lvalue = use_taylor.value();
  } else {
    _use_taylor = nullptr;
  }
  // intent=in allocatable type array
  auto* _orb0_arr = orb0_arr.has_value() ? orb0_arr->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  fortran_make_hybrid_lat(
      /* void* */ lat_in.get_fortran_ptr(),
      /* void* */ _lat_out.get_fortran_ptr(),
      /* bool* */ _use_taylor,
      /* void* */ _orb0_arr);
  return std::move(_lat_out);
}
Bmad::MakeMadMap Bmad::make_mad_map(EleProxy& ele, LatParamProxy& param) {
  MadEnergyProxy _energy;
  MadMapProxy _map;
  fortran_make_mad_map(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return MakeMadMap{std::move(_energy), std::move(_map)};
}
Bmad::MakeMat6 Bmad::make_mat6(
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> start_orb) {
  auto* _start_orb = start_orb.has_value() ? start_orb->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  CoordProxy _end_orb;
  bool _err_flag{};
  fortran_make_mat6(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _start_orb,
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return MakeMat6{std::move(_end_orb), _err_flag};
}
Bmad::MakeMat6Bmad Bmad::make_mat6_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb) {
  CoordProxy _end_orb;
  bool _err{};
  fortran_make_mat6_bmad(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool& */ _err);
  return MakeMat6Bmad{std::move(_end_orb), _err};
}
Bmad::MakeMat6BmadPhoton Bmad::make_mat6_bmad_photon(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb) {
  CoordProxy _end_orb;
  bool _err{};
  fortran_make_mat6_bmad_photon(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool& */ _err);
  return MakeMat6BmadPhoton{std::move(_end_orb), _err};
}
void Bmad::make_mat6_high_energy_space_charge(
    EleProxy& ele,
    LatParamProxy& param) {
  fortran_make_mat6_high_energy_space_charge(
      /* void* */ ele.get_fortran_ptr(), /* void* */ param.get_fortran_ptr());
}
CoordProxy Bmad::make_mat6_mad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& c0) {
  CoordProxy _c1;
  fortran_make_mat6_mad(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ c0.get_fortran_ptr(),
      /* void* */ _c1.get_fortran_ptr());
  return std::move(_c1);
}
CoordProxy Bmad::make_mat6_symp_lie_ptc(EleProxy& ele, CoordProxy& start_orb) {
  CoordProxy _end_orb;
  fortran_make_mat6_symp_lie_ptc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr());
  return std::move(_end_orb);
}
CoordProxy Bmad::make_mat6_taylor(
    EleProxy& ele,
    CoordProxy& start_orb,
    optional_ref<bool> err_flag) {
  CoordProxy _end_orb;
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  fortran_make_mat6_taylor(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool* */ _err_flag);
  return std::move(_end_orb);
}
Bmad::MakeMat6Tracking Bmad::make_mat6_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb,
    std::optional<bool> spin_only) {
  CoordProxy _end_orb;
  bool _err_flag{};
  bool spin_only_lvalue;
  auto* _spin_only{&spin_only_lvalue};
  if (spin_only.has_value()) {
    spin_only_lvalue = spin_only.value();
  } else {
    _spin_only = nullptr;
  }
  fortran_make_mat6_tracking(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _spin_only);
  return MakeMat6Tracking{std::move(_end_orb), _err_flag};
}
Bmad::MakeN Bmad::make_n(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<FixedArray1D<Real, 3>> abz_tunes) {
  double _t6_vec[6 * 6];
  matrix_to_vec(t6, _t6_vec);
  FixedArray2D<Real, 6, 6> N;
  double _N_vec[6 * 6];
  bool _err_flag{};
  double* _abz_tunes =
      abz_tunes.has_value() ? abz_tunes.value().data() : nullptr;
  FixedArray1D<Real, 3> _tunes_out;
  FixedArray2D<Real, 6, 6> U;
  double _U_vec[6 * 6];
  fortran_make_n(
      /* double* */ _t6_vec,
      /* double* */ _N_vec,
      /* bool& */ _err_flag,
      /* double* */ _abz_tunes,
      /* double* */ _tunes_out.data(),
      /* double* */ _U_vec);
  vec_to_matrix(_N_vec, N);
  vec_to_matrix(_U_vec, U);
  return MakeN{N, _err_flag, _tunes_out, U};
}
Bmad::MakePbrh Bmad::make_pbrh(
    FixedArray2D<Real, 6, 6> M,
    FixedArray1D<Real, 3> abz_tunes) {
  double _M_vec[6 * 6];
  matrix_to_vec(M, _M_vec);
  FixedArray2D<Complex, 6, 6> P;
  std::complex<double> _P_vec[6 * 6];
  FixedArray2D<Complex, 6, 6> Bp;
  std::complex<double> _Bp_vec[6 * 6];
  FixedArray2D<Complex, 6, 6> R;
  std::complex<double> _R_vec[6 * 6];
  FixedArray2D<Complex, 6, 6> H;
  std::complex<double> _H_vec[6 * 6];
  auto* _abz_tunes = abz_tunes.data(); // CppWrapperGeneralArgument
  fortran_make_pbrh(
      /* double* */ _M_vec,
      /* std::complex<double>* */ _P_vec,
      /* std::complex<double>* */ _Bp_vec,
      /* std::complex<double>* */ _R_vec,
      /* std::complex<double>* */ _H_vec,
      /* double* */ _abz_tunes);
  vec_to_matrix(_P_vec, P);
  vec_to_matrix(_Bp_vec, Bp);
  vec_to_matrix(_R_vec, R);
  vec_to_matrix(_H_vec, H);
  return MakePbrh{P, Bp, R, H};
}
Bmad::MakeSmatFromAbc Bmad::make_smat_from_abc(
    FixedArray2D<Real, 6, 6> t6,
    NormalModesProxy& mode) {
  double _t6_vec[6 * 6];
  matrix_to_vec(t6, _t6_vec);
  FixedArray2D<Real, 6, 6> sigma_mat;
  double _sigma_mat_vec[6 * 6];
  bool _err_flag{};
  FixedArray2D<Real, 6, 6> Nout;
  double _Nout_vec[6 * 6];
  fortran_make_smat_from_abc(
      /* double* */ _t6_vec,
      /* void* */ mode.get_fortran_ptr(),
      /* double* */ _sigma_mat_vec,
      /* bool& */ _err_flag,
      /* double* */ _Nout_vec);
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  vec_to_matrix(_Nout_vec, Nout);
  return MakeSmatFromAbc{sigma_mat, _err_flag, Nout};
}
void Bmad::make_unit_mad_map(MadMapProxy& map) {
  fortran_make_unit_mad_map(/* void* */ map.get_fortran_ptr());
}
void Bmad::make_v(
    FixedArray2D<Real, 6, 6> M,
    FixedArray2D<Complex, 6, 6> V,
    FixedArray1D<Real, 3> abz_tunes) {
  double _M_vec[6 * 6];
  matrix_to_vec(M, _M_vec);
  std::complex<double> _V_vec[6 * 6];
  matrix_to_vec(V, _V_vec);
  auto* _abz_tunes = abz_tunes.data(); // CppWrapperGeneralArgument
  fortran_make_v(
      /* double* */ _M_vec,
      /* std::complex<double>* */ _V_vec,
      /* double* */ _abz_tunes);
  vec_to_matrix(_M_vec, M);
  vec_to_matrix(_V_vec, V);
}
Bmad::MakeVMats Bmad::make_v_mats(EleProxy& ele) {
  FixedArray2D<Real, 4, 4> v_mat;
  double _v_mat_vec[4 * 4];
  FixedArray2D<Real, 4, 4> v_inv_mat;
  double _v_inv_mat_vec[4 * 4];
  fortran_make_v_mats(
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _v_mat_vec,
      /* double* */ _v_inv_mat_vec);
  vec_to_matrix(_v_mat_vec, v_mat);
  vec_to_matrix(_v_inv_mat_vec, v_inv_mat);
  return MakeVMats{v_mat, v_inv_mat};
}
void Bmad::makeup_control_slave(
    LatProxy& lat,
    EleProxy& slave,
    bool& err_flag) {
  fortran_makeup_control_slave(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ slave.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::makeup_group_lord(LatProxy& lat, EleProxy& lord, bool& err_flag) {
  fortran_makeup_group_lord(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ lord.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::makeup_multipass_slave(
    LatProxy& lat,
    EleProxy& slave,
    bool& err_flag) {
  fortran_makeup_multipass_slave(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ slave.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::makeup_super_slave(LatProxy& lat, EleProxy& slave, bool& err_flag) {
  fortran_makeup_super_slave(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ slave.get_fortran_ptr(),
      /* bool& */ err_flag);
}
bool Bmad::makeup_super_slave1(
    EleProxy& slave,
    EleProxy& lord,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end) {
  bool _err_flag{};
  fortran_makeup_super_slave1(
      /* void* */ slave.get_fortran_ptr(),
      /* void* */ lord.get_fortran_ptr(),
      /* double& */ offset,
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ include_upstream_end,
      /* bool& */ include_downstream_end,
      /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::map1_inverse(
    SpinOrbitMap1Proxy& map1,
    SpinOrbitMap1Proxy& inv_map1) {
  fortran_map1_inverse(
      /* void* */ map1.get_fortran_ptr(),
      /* void* */ inv_map1.get_fortran_ptr());
}
SpinOrbitMap1Proxy Bmad::map1_make_unit() {
  SpinOrbitMap1Proxy _map1;
  fortran_map1_make_unit(/* void* */ _map1.get_fortran_ptr());
  return std::move(_map1);
}
void Bmad::map1_times_map1(
    SpinOrbitMap1Proxy& map2,
    SpinOrbitMap1Proxy& map1,
    SpinOrbitMap1Proxy& map_out) {
  fortran_map1_times_map1(
      /* void* */ map2.get_fortran_ptr(),
      /* void* */ map1.get_fortran_ptr(),
      /* void* */ map_out.get_fortran_ptr());
}
TaylorProxyArray1D Bmad::map_to_angle_coords(
    FixedArray1D<TaylorProxy, 6> t_canon) {
  // Output-only type array
  auto t_angle = TaylorProxyArray1D::allocate(6, 1);

  fortran_map_to_angle_coords(
      /* void* */ t_canon.data(), /* void* */ t_angle.get_fortran_ptr());
  return std::move(std::move(t_angle));
}
void Bmad::mark_patch_regions(BranchProxy& branch) {
  fortran_mark_patch_regions(/* void* */ branch.get_fortran_ptr());
}
void Bmad::master_parameter_value(
    int master_parameter,
    EleProxy& ele,
    double& value) {
  fortran_master_parameter_value(
      /* int& */ master_parameter,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ value);
}
FixedArray2D<Real, 4, 4> Bmad::mat4_multipole(
    double knl,
    double tilt,
    int& n,
    CoordProxy& orbit) {
  FixedArray2D<Real, 4, 4> kick_mat;
  double _kick_mat_vec[4 * 4];
  fortran_mat4_multipole(
      /* double& */ knl,
      /* double& */ tilt,
      /* int& */ n,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _kick_mat_vec);
  vec_to_matrix(_kick_mat_vec, kick_mat);
  return kick_mat;
}
void Bmad::mat6_add_offsets(EleProxy& ele, LatParamProxy& param) {
  fortran_mat6_add_offsets(
      /* void* */ ele.get_fortran_ptr(), /* void* */ param.get_fortran_ptr());
}
void Bmad::mat6_add_pitch(
    double x_pitch_tot,
    double y_pitch_tot,
    int orientation,
    FixedArray2D<Real, 6, 6> mat6) {
  double _mat6_vec[6 * 6];
  matrix_to_vec(mat6, _mat6_vec);
  fortran_mat6_add_pitch(
      /* double& */ x_pitch_tot,
      /* double& */ y_pitch_tot,
      /* int& */ orientation,
      /* double* */ _mat6_vec);
  vec_to_matrix(_mat6_vec, mat6);
}
ComplexTaylorProxyArray1D Bmad::mat6_to_complex_taylor(
    FixedArray1D<Complex, 6> vec0,
    FixedArray2D<Complex, 6, 6> mat6) {
  auto* _vec0 = vec0.data(); // CppWrapperGeneralArgument
  std::complex<double> _mat6_vec[6 * 6];
  matrix_to_vec(mat6, _mat6_vec);
  // Output-only type array
  auto complex_taylor = ComplexTaylorProxyArray1D::allocate(6, 1);

  fortran_mat6_to_complex_taylor(
      /* std::complex<double>* */ _vec0,
      /* std::complex<double>* */ _mat6_vec,
      /* void* */ complex_taylor.get_fortran_ptr());
  return std::move(std::move(complex_taylor));
}
Bmad::MatSympDecouple Bmad::mat_symp_decouple(
    FixedArray2D<Real, 4, 4> t0,
    FixedArray2D<Real, 4, 4> U,
    FixedArray2D<Real, 4, 4> V,
    FixedArray2D<Real, 4, 4> Ubar,
    FixedArray2D<Real, 4, 4> Vbar,
    FixedArray2D<Real, 4, 4> G,
    bool type_out) {
  double _t0_vec[4 * 4];
  matrix_to_vec(t0, _t0_vec);
  int _stat{};
  double _U_vec[4 * 4];
  matrix_to_vec(U, _U_vec);
  double _V_vec[4 * 4];
  matrix_to_vec(V, _V_vec);
  double _Ubar_vec[4 * 4];
  matrix_to_vec(Ubar, _Ubar_vec);
  double _Vbar_vec[4 * 4];
  matrix_to_vec(Vbar, _Vbar_vec);
  double _G_vec[4 * 4];
  matrix_to_vec(G, _G_vec);
  TwissProxy _twiss1;
  TwissProxy _twiss2;
  double _gamma{};
  fortran_mat_symp_decouple(
      /* double* */ _t0_vec,
      /* int& */ _stat,
      /* double* */ _U_vec,
      /* double* */ _V_vec,
      /* double* */ _Ubar_vec,
      /* double* */ _Vbar_vec,
      /* double* */ _G_vec,
      /* void* */ _twiss1.get_fortran_ptr(),
      /* void* */ _twiss2.get_fortran_ptr(),
      /* double& */ _gamma,
      /* bool& */ type_out);
  vec_to_matrix(_U_vec, U);
  vec_to_matrix(_V_vec, V);
  vec_to_matrix(_Ubar_vec, Ubar);
  vec_to_matrix(_Vbar_vec, Vbar);
  vec_to_matrix(_G_vec, G);
  return MatSympDecouple{_stat, std::move(_twiss1), std::move(_twiss2), _gamma};
}
Bmad::MatchEleToMat6 Bmad::match_ele_to_mat6(
    EleProxy& ele,
    CoordProxy& start_orb,
    std::optional<bool> include_delta_time,
    std::optional<bool> set_trombone) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  FixedArray1D<Real, 6> _vec0;
  bool _err_flag{};
  bool include_delta_time_lvalue;
  auto* _include_delta_time{&include_delta_time_lvalue};
  if (include_delta_time.has_value()) {
    include_delta_time_lvalue = include_delta_time.value();
  } else {
    _include_delta_time = nullptr;
  }
  bool set_trombone_lvalue;
  auto* _set_trombone{&set_trombone_lvalue};
  if (set_trombone.has_value()) {
    set_trombone_lvalue = set_trombone.value();
  } else {
    _set_trombone = nullptr;
  }
  fortran_match_ele_to_mat6(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ start_orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* double* */ _vec0.data(),
      /* bool& */ _err_flag,
      /* bool* */ _include_delta_time,
      /* bool* */ _set_trombone);
  vec_to_matrix(_mat6_vec, mat6);
  return MatchEleToMat6{mat6, _vec0, _err_flag};
}
void Bmad::mexp(double x, int m, double& this_exp) {
  fortran_mexp(/* double& */ x, /* int& */ m, /* double& */ this_exp);
}
int Bmad::mfft1(
    RealAlloc1D& a,
    RealAlloc1D& b,
    IntAlloc1D& n,
    int ndim,
    int isn) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  // intent=in allocatable general array
  int _ierr{};
  fortran_mfft1(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* void* */ n.get_fortran_ptr(),
      /* int& */ ndim,
      /* int& */ isn,
      /* int& */ _ierr);
  return _ierr;
}
void Bmad::momentum_compaction(BranchProxy& branch, double& mom_comp) {
  fortran_momentum_compaction(
      /* void* */ branch.get_fortran_ptr(), /* double& */ mom_comp);
}
Bmad::MultiTurnTrackingAnalysis Bmad::multi_turn_tracking_analysis(
    CoordProxyAlloc1D& track,
    int i_dim) {
  // intent=in allocatable type array
  CoordProxy _track0;
  EleProxy _ele;
  bool _stable{};
  double _growth_rate{};
  double _chi{};
  bool _err_flag{};
  fortran_multi_turn_tracking_analysis(
      /* void* */ track.get_fortran_ptr(),
      /* int& */ i_dim,
      /* void* */ _track0.get_fortran_ptr(),
      /* void* */ _ele.get_fortran_ptr(),
      /* bool& */ _stable,
      /* double& */ _growth_rate,
      /* double& */ _chi,
      /* bool& */ _err_flag);
  return MultiTurnTrackingAnalysis{
      std::move(_track0),
      std::move(_ele),
      _stable,
      _growth_rate,
      _chi,
      _err_flag};
}
bool Bmad::multilayer_type_to_multilayer_params(EleProxy& ele) {
  bool _err_flag{};
  fortran_multilayer_type_to_multilayer_params(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::multipass_chain(
    EleProxy& ele,
    int ix_pass,
    int n_links,
    optional_ref<ElePointerProxyAlloc1D> chain_ele,
    std::optional<bool> use_super_lord) {
  // intent=in allocatable type array
  auto* _chain_ele = chain_ele.has_value() ? chain_ele->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  bool use_super_lord_lvalue;
  auto* _use_super_lord{&use_super_lord_lvalue};
  if (use_super_lord.has_value()) {
    use_super_lord_lvalue = use_super_lord.value();
  } else {
    _use_super_lord = nullptr;
  }
  fortran_multipass_chain(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ ix_pass,
      /* int& */ n_links,
      /* void* */ _chain_ele,
      /* bool* */ _use_super_lord);
}
Bmad::Multipole1AbToKt Bmad::multipole1_ab_to_kt(double an, double bn, int n) {
  double _knl{};
  double _tn{};
  fortran_multipole1_ab_to_kt(
      /* double& */ an,
      /* double& */ bn,
      /* int& */ n,
      /* double& */ _knl,
      /* double& */ _tn);
  return Multipole1AbToKt{_knl, _tn};
}
Bmad::Multipole1KtToAb Bmad::multipole1_kt_to_ab(
    double knl,
    double knsl,
    double tn,
    int n) {
  double _an{};
  double _bn{};
  fortran_multipole1_kt_to_ab(
      /* double& */ knl,
      /* double& */ knsl,
      /* double& */ tn,
      /* int& */ n,
      /* double& */ _an,
      /* double& */ _bn);
  return Multipole1KtToAb{_an, _bn};
}
Bmad::MultipoleAbToKt Bmad::multipole_ab_to_kt(
    RealAlloc1D& an,
    RealAlloc1D& bn) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=out allocatable general array
  auto knl{RealAlloc1D()};
  // intent=out allocatable general array
  auto tn{RealAlloc1D()};
  fortran_multipole_ab_to_kt(
      /* void* */ an.get_fortran_ptr(),
      /* void* */ bn.get_fortran_ptr(),
      /* void* */ knl.get_fortran_ptr(),
      /* void* */ tn.get_fortran_ptr());
  return MultipoleAbToKt{std::move(knl), std::move(tn)};
}
Bmad::MultipoleEleToAb Bmad::multipole_ele_to_ab(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type,
    std::optional<int> include_kicks,
    std::optional<bool> original) {
  int _ix_pole_max{};
  FixedArray1D<Real, Bmad::N_POLE_MAXX> _a;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> _b;
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  int include_kicks_lvalue;
  auto* _include_kicks{&include_kicks_lvalue};
  if (include_kicks.has_value()) {
    include_kicks_lvalue = include_kicks.value();
  } else {
    _include_kicks = nullptr;
  }
  double _b1{};
  bool original_lvalue;
  auto* _original{&original_lvalue};
  if (original.has_value()) {
    original_lvalue = original.value();
  } else {
    _original = nullptr;
  }
  fortran_multipole_ele_to_ab(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ use_ele_tilt,
      /* int& */ _ix_pole_max,
      /* double* */ _a.data(),
      /* double* */ _b.data(),
      /* int* */ _pole_type,
      /* int* */ _include_kicks,
      /* double& */ _b1,
      /* bool* */ _original);
  return MultipoleEleToAb{_ix_pole_max, _a, _b, _b1};
}
Bmad::MultipoleEleToKt Bmad::multipole_ele_to_kt(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type,
    std::optional<int> include_kicks) {
  int _ix_pole_max{};
  // intent=out allocatable general array
  auto knl{RealAlloc1D()};
  // intent=out allocatable general array
  auto tilt{RealAlloc1D()};
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  int include_kicks_lvalue;
  auto* _include_kicks{&include_kicks_lvalue};
  if (include_kicks.has_value()) {
    include_kicks_lvalue = include_kicks.value();
  } else {
    _include_kicks = nullptr;
  }
  fortran_multipole_ele_to_kt(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ use_ele_tilt,
      /* int& */ _ix_pole_max,
      /* void* */ knl.get_fortran_ptr(),
      /* void* */ tilt.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* int* */ _include_kicks);
  return MultipoleEleToKt{_ix_pole_max, std::move(knl), std::move(tilt)};
}
EleProxy Bmad::multipole_init(int who, std::optional<bool> zero) {
  EleProxy _ele;
  bool zero_lvalue;
  auto* _zero{&zero_lvalue};
  if (zero.has_value()) {
    zero_lvalue = zero.value();
  } else {
    _zero = nullptr;
  }
  fortran_multipole_init(
      /* void* */ _ele.get_fortran_ptr(), /* int& */ who, /* bool* */ _zero);
  return std::move(_ele);
}
void Bmad::multipole_kick(
    double knl,
    double tilt,
    int n,
    int ref_species,
    int ele_orientation,
    CoordProxy& coord,
    std::optional<int> pole_type,
    std::optional<bool> ref_orb_offset) {
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  bool ref_orb_offset_lvalue;
  auto* _ref_orb_offset{&ref_orb_offset_lvalue};
  if (ref_orb_offset.has_value()) {
    ref_orb_offset_lvalue = ref_orb_offset.value();
  } else {
    _ref_orb_offset = nullptr;
  }
  fortran_multipole_kick(
      /* double& */ knl,
      /* double& */ tilt,
      /* int& */ n,
      /* int& */ ref_species,
      /* int& */ ele_orientation,
      /* void* */ coord.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* bool* */ _ref_orb_offset);
}
FixedArray2D<Real, 6, 6> Bmad::multipole_kick_mat(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    int ref_species,
    EleProxy& ele,
    CoordProxy& orbit,
    double factor) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  fortran_multipole_kick_mat(
      /* void* */ knl.get_fortran_ptr(),
      /* void* */ tilt.get_fortran_ptr(),
      /* int& */ ref_species,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* double& */ factor,
      /* double* */ _mat6_vec);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::multipole_kicks(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<int> pole_type,
    std::optional<bool> ref_orb_offset) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  int pole_type_lvalue;
  auto* _pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  bool ref_orb_offset_lvalue;
  auto* _ref_orb_offset{&ref_orb_offset_lvalue};
  if (ref_orb_offset.has_value()) {
    ref_orb_offset_lvalue = ref_orb_offset.value();
  } else {
    _ref_orb_offset = nullptr;
  }
  fortran_multipole_kicks(
      /* void* */ knl.get_fortran_ptr(),
      /* void* */ tilt.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* bool* */ _ref_orb_offset);
}
Bmad::MultipoleKtToAb Bmad::multipole_kt_to_ab(
    RealAlloc1D& knl,
    RealAlloc1D& knsl,
    RealAlloc1D& tn) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=out allocatable general array
  auto an{RealAlloc1D()};
  // intent=out allocatable general array
  auto bn{RealAlloc1D()};
  fortran_multipole_kt_to_ab(
      /* void* */ knl.get_fortran_ptr(),
      /* void* */ knsl.get_fortran_ptr(),
      /* void* */ tn.get_fortran_ptr(),
      /* void* */ an.get_fortran_ptr(),
      /* void* */ bn.get_fortran_ptr());
  return MultipoleKtToAb{std::move(an), std::move(bn)};
}
void Bmad::multipole_spin_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_multipole_spin_tracking(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::mytan(double& y, double& x, double& arg) {
  fortran_mytan(/* double& */ y, /* double& */ x, /* double& */ arg);
}
int Bmad::n_attrib_string_max_len() {
  int _max_len{};
  fortran_n_attrib_string_max_len(/* int& */ _max_len);
  return _max_len;
}
void Bmad::new_control(
    LatProxy& lat,
    int ix_ele,
    std::optional<std::string> ele_name) {
  const char* _ele_name = ele_name.has_value() ? ele_name->c_str() : nullptr;
  fortran_new_control(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ ix_ele,
      /* const char* */ _ele_name);
}
int Bmad::nint_chk(double re_val) {
  int _int_val{};
  fortran_nint_chk(/* double& */ re_val, /* int& */ _int_val);
  return _int_val;
}
void Bmad::normal_form_complex_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool& rf_on,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> F,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> L,
    std::optional<FixedArray1D<TaylorProxy, 6>> A,
    std::optional<FixedArray1D<TaylorProxy, 6>> A_inverse,
    optional_ref<int> order) {
  auto* _F = F.has_value() ? F->data() : nullptr; // input, optional
  auto* _L = L.has_value() ? L->data() : nullptr; // input, optional
  auto* _A = A.has_value() ? A->data() : nullptr; // input, optional
  auto* _A_inverse =
      A_inverse.has_value() ? A_inverse->data() : nullptr; // input, optional
  auto* _order = order.has_value() ? &order->get() : nullptr; // inout, optional
  fortran_normal_form_complex_taylors(
      /* void* */ one_turn_taylor.data(),
      /* bool& */ rf_on,
      /* void* */ _F,
      /* void* */ _L,
      /* void* */ _A,
      /* void* */ _A_inverse,
      /* int* */ _order);
}
Bmad::NormalFormTaylors Bmad::normal_form_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool rf_on) {
  // Output-only type array
  auto dhdj = TaylorProxyArray1D::allocate(6, 1);

  // Output-only type array
  auto A = TaylorProxyArray1D::allocate(6, 1);

  // Output-only type array
  auto A_inverse = TaylorProxyArray1D::allocate(6, 1);

  fortran_normal_form_taylors(
      /* void* */ one_turn_taylor.data(),
      /* bool& */ rf_on,
      /* void* */ dhdj.get_fortran_ptr(),
      /* void* */ A.get_fortran_ptr(),
      /* void* */ A_inverse.get_fortran_ptr());
  return NormalFormTaylors{
      std::move(std::move(dhdj)),
      std::move(std::move(A)),
      std::move(std::move(A_inverse))};
}
Bmad::NormalMode3Calc Bmad::normal_mode3_calc(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<bool> above_transition,
    std::optional<FixedArray1D<Real, 3>> abz_tunes) {
  double _t6_vec[6 * 6];
  matrix_to_vec(t6, _t6_vec);
  FixedArray1D<Real, 3> _tune;
  FixedArray2D<Real, 6, 6> B;
  double _B_vec[6 * 6];
  FixedArray2D<Real, 6, 6> HV;
  double _HV_vec[6 * 6];
  bool above_transition_lvalue;
  auto* _above_transition{&above_transition_lvalue};
  if (above_transition.has_value()) {
    above_transition_lvalue = above_transition.value();
  } else {
    _above_transition = nullptr;
  }
  double* _abz_tunes =
      abz_tunes.has_value() ? abz_tunes.value().data() : nullptr;
  fortran_normal_mode3_calc(
      /* double* */ _t6_vec,
      /* double* */ _tune.data(),
      /* double* */ _B_vec,
      /* double* */ _HV_vec,
      /* bool* */ _above_transition,
      /* double* */ _abz_tunes);
  vec_to_matrix(_t6_vec, t6);
  vec_to_matrix(_B_vec, B);
  vec_to_matrix(_HV_vec, HV);
  return NormalMode3Calc{_tune, B, HV};
}
void Bmad::normal_mode_dispersion(EleProxy& ele, std::optional<bool> reverse) {
  bool reverse_lvalue;
  auto* _reverse{&reverse_lvalue};
  if (reverse.has_value()) {
    reverse_lvalue = reverse.value();
  } else {
    _reverse = nullptr;
  }
  fortran_normal_mode_dispersion(
      /* void* */ ele.get_fortran_ptr(), /* bool* */ _reverse);
}
bool Bmad::normalize_evecs(FixedArray2D<Complex, 6, 6> evec) {
  std::complex<double> _evec_vec[6 * 6];
  matrix_to_vec(evec, _evec_vec);
  bool _err_flag{};
  fortran_normalize_evecs(
      /* std::complex<double>* */ _evec_vec, /* bool& */ _err_flag);
  vec_to_matrix(_evec_vec, evec);
  return _err_flag;
}
void Bmad::num_field_eles(EleProxy& ele, int& n_field_ele) {
  fortran_num_field_eles(
      /* void* */ ele.get_fortran_ptr(), /* int& */ n_field_ele);
}
void Bmad::num_lords(EleProxy& slave, int lord_type, int& num) {
  fortran_num_lords(
      /* void* */ slave.get_fortran_ptr(),
      /* int& */ lord_type,
      /* int& */ num);
}
Bmad::OdeintBmad Bmad::odeint_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    double s1_body,
    double s2_body,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  bool _err_flag{};
  TrackProxy _track;
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_odeint_bmad(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ s1_body,
      /* double& */ s2_body,
      /* bool& */ _err_flag,
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return OdeintBmad{_err_flag, std::move(_track)};
}
Bmad::OdeintBmadTime Bmad::odeint_bmad_time(
    CoordProxy& orb,
    EleProxy& ele,
    LatParamProxy& param,
    int t_dir,
    double& rf_time,
    optional_ref<TrackProxy> track,
    std::optional<double> t_end,
    optional_ref<EmFieldProxy> extra_field) {
  bool _err_flag{};
  auto* _track = track.has_value() ? track->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  double t_end_lvalue;
  auto* _t_end{&t_end_lvalue};
  if (t_end.has_value()) {
    t_end_lvalue = t_end.value();
  } else {
    _t_end = nullptr;
  }
  double _dt_step{};
  auto* _extra_field = extra_field.has_value()
      ? extra_field->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_odeint_bmad_time(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ t_dir,
      /* double& */ rf_time,
      /* bool& */ _err_flag,
      /* void* */ _track,
      /* double* */ _t_end,
      /* double& */ _dt_step,
      /* void* */ _extra_field);
  return OdeintBmadTime{_err_flag, _dt_step};
}
Bmad::OffsetParticle Bmad::offset_particle(
    EleProxy& ele,
    bool set,
    CoordProxy& orbit,
    std::optional<bool> set_tilt,
    std::optional<bool> set_hvkicks,
    std::optional<int> drift_to_edge,
    std::optional<double> s_pos,
    std::optional<bool> set_spin,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    optional_ref<double> time) {
  bool set_tilt_lvalue;
  auto* _set_tilt{&set_tilt_lvalue};
  if (set_tilt.has_value()) {
    set_tilt_lvalue = set_tilt.value();
  } else {
    _set_tilt = nullptr;
  }
  bool set_hvkicks_lvalue;
  auto* _set_hvkicks{&set_hvkicks_lvalue};
  if (set_hvkicks.has_value()) {
    set_hvkicks_lvalue = set_hvkicks.value();
  } else {
    _set_hvkicks = nullptr;
  }
  int drift_to_edge_lvalue;
  auto* _drift_to_edge{&drift_to_edge_lvalue};
  if (drift_to_edge.has_value()) {
    drift_to_edge_lvalue = drift_to_edge.value();
  } else {
    _drift_to_edge = nullptr;
  }
  double s_pos_lvalue;
  auto* _s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  double _s_out{};
  bool set_spin_lvalue;
  auto* _set_spin{&set_spin_lvalue};
  if (set_spin.has_value()) {
    set_spin_lvalue = set_spin.value();
  } else {
    _set_spin = nullptr;
  }
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  FixedArray1D<Real, 4> _spin_qrot;
  auto* _time = time.has_value() ? &time->get() : nullptr; // inout, optional
  fortran_offset_particle(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ set,
      /* void* */ orbit.get_fortran_ptr(),
      /* bool* */ _set_tilt,
      /* bool* */ _set_hvkicks,
      /* int* */ _drift_to_edge,
      /* double* */ _s_pos,
      /* double& */ _s_out,
      /* bool* */ _set_spin,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* double* */ _spin_qrot.data(),
      /* double* */ _time);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return OffsetParticle{_s_out, _spin_qrot};
}
void Bmad::offset_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    std::optional<bool> offset_position_only,
    std::optional<FixedArray2D<Real, 3, 3>> rot_mat) {
  bool offset_position_only_lvalue;
  auto* _offset_position_only{&offset_position_only_lvalue};
  if (offset_position_only.has_value()) {
    offset_position_only_lvalue = offset_position_only.value();
  } else {
    _offset_position_only = nullptr;
  }
  double _rot_mat_vec[3 * 3];
  const double* _rot_mat = nullptr;
  if (rot_mat.has_value()) {
    matrix_to_vec(rot_mat.value(), _rot_mat_vec);
    _rot_mat = _rot_mat_vec;
  }
  fortran_offset_photon(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ set,
      /* bool* */ _offset_position_only,
      /* double* */ _rot_mat_vec);
}
FixedArray2D<Real, 4, 4> Bmad::one_turn_mat_at_ele(
    EleProxy& ele,
    double phi_a,
    double phi_b) {
  FixedArray2D<Real, 4, 4> mat4;
  double _mat4_vec[4 * 4];
  fortran_one_turn_mat_at_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ phi_a,
      /* double& */ phi_b,
      /* double* */ _mat4_vec);
  vec_to_matrix(_mat4_vec, mat4);
  return mat4;
}
Bmad::OpenBinaryFile Bmad::open_binary_file(
    std::string file_name,
    std::string action,
    std::string r_name) {
  auto _file_name = file_name.c_str();
  auto _action = action.c_str();
  int _iu{};
  auto _r_name = r_name.c_str();
  int _iver{};
  bool _is_ok{};
  fortran_open_binary_file(
      /* const char* */ _file_name,
      /* const char* */ _action,
      /* int& */ _iu,
      /* const char* */ _r_name,
      /* int& */ _iver,
      /* bool& */ _is_ok);
  return OpenBinaryFile{_iu, _iver, _is_ok};
}
Bmad::OrbitAmplitudeCalc Bmad::orbit_amplitude_calc(
    EleProxy& ele,
    CoordProxy& orb) {
  double _amp_a{};
  double _amp_b{};
  double _amp_na{};
  double _amp_nb{};
  fortran_orbit_amplitude_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orb.get_fortran_ptr(),
      /* double& */ _amp_a,
      /* double& */ _amp_b,
      /* double& */ _amp_na,
      /* double& */ _amp_nb);
  return OrbitAmplitudeCalc{_amp_a, _amp_b, _amp_na, _amp_nb};
}
void Bmad::orbit_reference_energy_correction(
    CoordProxy& orbit,
    double p0c_new,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_orbit_reference_energy_correction(
      /* void* */ orbit.get_fortran_ptr(),
      /* double& */ p0c_new,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::orbit_to_floor_phase_space(
    CoordProxy& orbit,
    EleProxy& ele,
    FixedArray1D<Real, 6> floor_phase_space) {
  auto* _floor_phase_space =
      floor_phase_space.data(); // CppWrapperGeneralArgument
  fortran_orbit_to_floor_phase_space(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _floor_phase_space);
}
void Bmad::orbit_to_local_curvilinear(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<int> z_direction,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position) {
  int z_direction_lvalue;
  auto* _z_direction{&z_direction_lvalue};
  if (z_direction.has_value()) {
    z_direction_lvalue = z_direction.value();
  } else {
    _z_direction = nullptr;
  }
  int relative_to_lvalue;
  auto* _relative_to{&relative_to_lvalue};
  if (relative_to.has_value()) {
    relative_to_lvalue = relative_to.value();
  } else {
    _relative_to = nullptr;
  }
  fortran_orbit_to_local_curvilinear(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _z_direction,
      /* int* */ _relative_to,
      /* void* */ local_position.get_fortran_ptr());
}
LatParamProxy Bmad::orbit_too_large(
    CoordProxy& orbit,
    std::optional<bool> check_momentum,
    bool& is_too_large) {
  LatParamProxy _param;
  bool check_momentum_lvalue;
  auto* _check_momentum{&check_momentum_lvalue};
  if (check_momentum.has_value()) {
    check_momentum_lvalue = check_momentum.value();
  } else {
    _check_momentum = nullptr;
  }
  fortran_orbit_too_large(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ _param.get_fortran_ptr(),
      /* bool* */ _check_momentum,
      /* bool& */ is_too_large);
  return std::move(_param);
}
Bmad::OrderEvecsByNSimilarity Bmad::order_evecs_by_n_similarity(
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray2D<Real, 6, 6> Nmat) {
  FixedArray2D<Complex, 6, 6> evec;
  std::complex<double> _evec_vec[6 * 6];
  auto* _eval = eval.data(); // CppWrapperGeneralArgument
  auto* _mat_tunes = mat_tunes.data(); // CppWrapperGeneralArgument
  double _Nmat_vec[6 * 6];
  matrix_to_vec(Nmat, _Nmat_vec);
  bool _err_flag{};
  fortran_order_evecs_by_n_similarity(
      /* std::complex<double>* */ _evec_vec,
      /* std::complex<double>* */ _eval,
      /* double* */ _mat_tunes,
      /* double* */ _Nmat_vec,
      /* bool& */ _err_flag);
  vec_to_matrix(_evec_vec, evec);
  return OrderEvecsByNSimilarity{evec, _err_flag};
}
void Bmad::order_evecs_by_plane_dominance(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    std::optional<FixedArray1D<Real, 3>> mat_tunes) {
  std::complex<double> _evec_vec[6 * 6];
  matrix_to_vec(evec, _evec_vec);
  auto* _eval = eval.data(); // CppWrapperGeneralArgument
  double* _mat_tunes =
      mat_tunes.has_value() ? mat_tunes.value().data() : nullptr;
  fortran_order_evecs_by_plane_dominance(
      /* std::complex<double>* */ _evec_vec,
      /* std::complex<double>* */ _eval,
      /* double* */ _mat_tunes);
  vec_to_matrix(_evec_vec, evec);
}
bool Bmad::order_evecs_by_tune(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray1D<Real, 3> abz_tunes) {
  std::complex<double> _evec_vec[6 * 6];
  matrix_to_vec(evec, _evec_vec);
  auto* _eval = eval.data(); // CppWrapperGeneralArgument
  auto* _mat_tunes = mat_tunes.data(); // CppWrapperGeneralArgument
  auto* _abz_tunes = abz_tunes.data(); // CppWrapperGeneralArgument
  bool _err_flag{};
  fortran_order_evecs_by_tune(
      /* std::complex<double>* */ _evec_vec,
      /* std::complex<double>* */ _eval,
      /* double* */ _mat_tunes,
      /* double* */ _abz_tunes,
      /* bool& */ _err_flag);
  vec_to_matrix(_evec_vec, evec);
  return _err_flag;
}
void Bmad::order_particles_in_z(BunchProxy& bunch) {
  fortran_order_particles_in_z(/* void* */ bunch.get_fortran_ptr());
}
void Bmad::order_super_lord_slaves(LatProxy& lat, int ix_lord) {
  fortran_order_super_lord_slaves(
      /* void* */ lat.get_fortran_ptr(), /* int& */ ix_lord);
}
void Bmad::osc_alloc_freespace_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad) {
  auto* _nlo = nlo.data(); // CppWrapperGeneralArgument
  auto* _nhi = nhi.data(); // CppWrapperGeneralArgument
  auto* _npad = npad.data(); // CppWrapperGeneralArgument
  fortran_osc_alloc_freespace_array(
      /* int* */ _nlo, /* int* */ _nhi, /* int* */ _npad);
}
void Bmad::osc_alloc_image_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad) {
  auto* _nlo = nlo.data(); // CppWrapperGeneralArgument
  auto* _nhi = nhi.data(); // CppWrapperGeneralArgument
  auto* _npad = npad.data(); // CppWrapperGeneralArgument
  fortran_osc_alloc_image_array(
      /* int* */ _nlo, /* int* */ _nhi, /* int* */ _npad);
}
void Bmad::osc_alloc_rectpipe_arrays(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad) {
  auto* _nlo = nlo.data(); // CppWrapperGeneralArgument
  auto* _nhi = nhi.data(); // CppWrapperGeneralArgument
  auto* _npad = npad.data(); // CppWrapperGeneralArgument
  fortran_osc_alloc_rectpipe_arrays(
      /* int* */ _nlo, /* int* */ _nhi, /* int* */ _npad);
}
void Bmad::osc_getgrnpipe(
    double& gam,
    double& a,
    double& b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad) {
  auto* _delta = delta.data(); // CppWrapperGeneralArgument
  auto* _umin = umin.data(); // CppWrapperGeneralArgument
  auto* _npad = npad.data(); // CppWrapperGeneralArgument
  fortran_osc_getgrnpipe(
      /* double& */ gam,
      /* double& */ a,
      /* double& */ b,
      /* double* */ _delta,
      /* double* */ _umin,
      /* int* */ _npad);
}
void Bmad::osc_read_rectpipe_grn() {
  fortran_osc_read_rectpipe_grn();
}
void Bmad::osc_write_rectpipe_grn(
    double& apipe,
    double& bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double& gamma) {
  auto* _delta = delta.data(); // CppWrapperGeneralArgument
  auto* _umin = umin.data(); // CppWrapperGeneralArgument
  auto* _umax = umax.data(); // CppWrapperGeneralArgument
  auto* _nlo = nlo.data(); // CppWrapperGeneralArgument
  auto* _nhi = nhi.data(); // CppWrapperGeneralArgument
  fortran_osc_write_rectpipe_grn(
      /* double& */ apipe,
      /* double& */ bpipe,
      /* double* */ _delta,
      /* double* */ _umin,
      /* double* */ _umax,
      /* int* */ _nlo,
      /* int* */ _nhi,
      /* double& */ gamma);
}
void Bmad::parse_cartesian_map(
    CartesianMapProxy& ct_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parse_cartesian_map(
      /* void* */ ct_map.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
void Bmad::parse_cylindrical_map(
    CylindricalMapProxy& cl_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _cl_map = &cl_map; // input, required, pointer
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parse_cylindrical_map(
      /* void* */ &cl_map,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
void Bmad::parse_gen_grad_map(
    GenGradMapProxy& gg_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _gg_map = &gg_map; // input, required, pointer
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parse_gen_grad_map(
      /* void* */ &gg_map,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
void Bmad::parse_grid_field(
    GridFieldProxy& g_field,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _g_field = &g_field; // input, required, pointer
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parse_grid_field(
      /* void* */ &g_field,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
void Bmad::parse_integer_list(
    std::string& err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    bool& exact_size,
    std::string& delim,
    bool& delim_found,
    optional_ref<std::string> open_delim,
    optional_ref<std::string> separator,
    optional_ref<std::string> close_delim,
    optional_ref<int> default_value,
    bool& is_ok) {
  auto _err_str = err_str.c_str(); // ptr, inout, required
  // intent=inout allocatable general array
  auto _delim = delim.c_str(); // ptr, inout, required
  const char* _open_delim =
      open_delim.has_value() ? open_delim->get().c_str() : nullptr;
  const char* _separator =
      separator.has_value() ? separator->get().c_str() : nullptr;
  const char* _close_delim =
      close_delim.has_value() ? close_delim->get().c_str() : nullptr;
  auto* _default_value = default_value.has_value() ? &default_value->get()
                                                   : nullptr; // inout, optional
  fortran_parse_integer_list(
      /* const char* */ _err_str,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ int_array.get_fortran_ptr(),
      /* bool& */ exact_size,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* const char* */ _open_delim,
      /* const char* */ _separator,
      /* const char* */ _close_delim,
      /* int* */ _default_value,
      /* bool& */ is_ok);
}
Bmad::ParseIntegerList2 Bmad::parse_integer_list2(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    optional_ref<int> num_expected,
    optional_ref<std::string> open_delim,
    optional_ref<std::string> separator,
    optional_ref<std::string> close_delim,
    optional_ref<int> default_value) {
  auto _err_str = err_str.c_str();
  // intent=inout allocatable general array
  int _num_found{};
  char _delim[4096];
  bool _delim_found{};
  auto* _num_expected = num_expected.has_value() ? &num_expected->get()
                                                 : nullptr; // inout, optional
  const char* _open_delim =
      open_delim.has_value() ? open_delim->get().c_str() : nullptr;
  const char* _separator =
      separator.has_value() ? separator->get().c_str() : nullptr;
  const char* _close_delim =
      close_delim.has_value() ? close_delim->get().c_str() : nullptr;
  auto* _default_value = default_value.has_value() ? &default_value->get()
                                                   : nullptr; // inout, optional
  bool _is_ok{};
  fortran_parse_integer_list2(
      /* const char* */ _err_str,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ int_array.get_fortran_ptr(),
      /* int& */ _num_found,
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* int* */ _num_expected,
      /* const char* */ _open_delim,
      /* const char* */ _separator,
      /* const char* */ _close_delim,
      /* int* */ _default_value,
      /* bool& */ _is_ok);
  return ParseIntegerList2{_num_found, _delim, _delim_found, _is_ok};
}
Bmad::ParseRealList Bmad::parse_real_list(
    LatProxy& lat,
    std::string err_str,
    bool exact_size,
    std::optional<std::string> open_delim,
    std::optional<std::string> separator,
    std::optional<std::string> close_delim,
    std::optional<double> default_value) {
  auto _err_str = err_str.c_str();
  // intent=out allocatable general array
  auto real_array{RealAlloc1D()};
  char _delim[4096];
  bool _delim_found{};
  const char* _open_delim =
      open_delim.has_value() ? open_delim->c_str() : nullptr;
  const char* _separator = separator.has_value() ? separator->c_str() : nullptr;
  const char* _close_delim =
      close_delim.has_value() ? close_delim->c_str() : nullptr;
  double default_value_lvalue;
  auto* _default_value{&default_value_lvalue};
  if (default_value.has_value()) {
    default_value_lvalue = default_value.value();
  } else {
    _default_value = nullptr;
  }
  int _num_found{};
  bool _is_ok{};
  fortran_parse_real_list(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _err_str,
      /* void* */ real_array.get_fortran_ptr(),
      /* bool& */ exact_size,
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* const char* */ _open_delim,
      /* const char* */ _separator,
      /* const char* */ _close_delim,
      /* double* */ _default_value,
      /* int& */ _num_found,
      /* bool& */ _is_ok);
  return ParseRealList{
      std::move(real_array), _delim, _delim_found, _num_found, _is_ok};
}
Bmad::ParseRealList2 Bmad::parse_real_list2(
    LatProxy& lat,
    std::string err_str,
    RealAlloc1D& real_array,
    optional_ref<int> num_expected,
    optional_ref<std::string> open_brace,
    optional_ref<std::string> separator,
    optional_ref<std::string> close_brace,
    optional_ref<double> default_value,
    optional_ref<bool> single_value) {
  auto _err_str = err_str.c_str();
  // intent=inout allocatable general array
  int _num_found{};
  char _delim[4096];
  bool _delim_found{};
  auto* _num_expected = num_expected.has_value() ? &num_expected->get()
                                                 : nullptr; // inout, optional
  const char* _open_brace =
      open_brace.has_value() ? open_brace->get().c_str() : nullptr;
  const char* _separator =
      separator.has_value() ? separator->get().c_str() : nullptr;
  const char* _close_brace =
      close_brace.has_value() ? close_brace->get().c_str() : nullptr;
  auto* _default_value = default_value.has_value() ? &default_value->get()
                                                   : nullptr; // inout, optional
  auto* _single_value = single_value.has_value() ? &single_value->get()
                                                 : nullptr; // inout, optional
  bool _is_ok{};
  fortran_parse_real_list2(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _err_str,
      /* void* */ real_array.get_fortran_ptr(),
      /* int& */ _num_found,
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* int* */ _num_expected,
      /* const char* */ _open_brace,
      /* const char* */ _separator,
      /* const char* */ _close_brace,
      /* double* */ _default_value,
      /* bool* */ _single_value,
      /* bool& */ _is_ok);
  return ParseRealList2{_num_found, _delim, _delim_found, _is_ok};
}
void Bmad::parser_add_constant(
    std::string& word,
    LatProxy& lat,
    bool& redef_is_error) {
  auto _word = word.c_str(); // ptr, inout, required
  fortran_parser_add_constant(
      /* const char* */ _word,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ redef_is_error);
}
void Bmad::parser_call_check(
    std::string& word,
    int& ix_word,
    std::string& delim,
    bool& delim_found,
    bool& call_found,
    optional_ref<bool> err_flag) {
  auto _word = word.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  fortran_parser_call_check(
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ call_found,
      /* bool* */ _err_flag);
}
Bmad::ParserFastComplexRead Bmad::parser_fast_complex_read(
    EleProxy& ele,
    std::string err_str) {
  // intent=out allocatable general array
  auto cmplx_vec{ComplexAlloc1D()};
  char _delim[4096];
  auto _err_str = err_str.c_str();
  bool _is_ok{};
  fortran_parser_fast_complex_read(
      /* void* */ cmplx_vec.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* const char* */ _err_str,
      /* bool& */ _is_ok);
  return ParserFastComplexRead{std::move(cmplx_vec), _delim, _is_ok};
}
void Bmad::parser_fast_integer_read(
    IntAlloc1D& int_vec,
    EleProxy& ele,
    std::string& delim_wanted,
    std::string& err_str,
    bool& is_ok) {
  // intent=inout allocatable general array
  auto _delim_wanted = delim_wanted.c_str(); // ptr, inout, required
  auto _err_str = err_str.c_str(); // ptr, inout, required
  fortran_parser_fast_integer_read(
      /* void* */ int_vec.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim_wanted,
      /* const char* */ _err_str,
      /* bool& */ is_ok);
}
Bmad::ParserFastRealRead Bmad::parser_fast_real_read(
    EleProxy& ele,
    std::string end_delims,
    std::string err_str,
    std::optional<bool> exact_size) {
  // intent=out allocatable general array
  auto real_vec{RealAlloc1D()};
  auto _end_delims = end_delims.c_str();
  char _delim[4096];
  auto _err_str = err_str.c_str();
  bool exact_size_lvalue;
  auto* _exact_size{&exact_size_lvalue};
  if (exact_size.has_value()) {
    exact_size_lvalue = exact_size.value();
  } else {
    _exact_size = nullptr;
  }
  int _n_real{};
  bool _is_ok{};
  fortran_parser_fast_real_read(
      /* void* */ real_vec.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _end_delims,
      /* const char* */ _delim,
      /* const char* */ _err_str,
      /* bool* */ _exact_size,
      /* int& */ _n_real,
      /* bool& */ _is_ok);
  return ParserFastRealRead{std::move(real_vec), _delim, _n_real, _is_ok};
}
void Bmad::parser_file_stack(
    std::string& how,
    optional_ref<std::string> file_name_in,
    optional_ref<bool> finished,
    optional_ref<bool> err,
    optional_ref<bool> open_file,
    optional_ref<bool> abort_on_open_error) {
  auto _how = how.c_str(); // ptr, inout, required
  const char* _file_name_in =
      file_name_in.has_value() ? file_name_in->get().c_str() : nullptr;
  auto* _finished =
      finished.has_value() ? &finished->get() : nullptr; // inout, optional
  auto* _err = err.has_value() ? &err->get() : nullptr; // inout, optional
  auto* _open_file =
      open_file.has_value() ? &open_file->get() : nullptr; // inout, optional
  auto* _abort_on_open_error = abort_on_open_error.has_value()
      ? &abort_on_open_error->get()
      : nullptr; // inout, optional
  fortran_parser_file_stack(
      /* const char* */ _how,
      /* const char* */ _file_name_in,
      /* bool* */ _finished,
      /* bool* */ _err,
      /* bool* */ _open_file,
      /* bool* */ _abort_on_open_error);
}
void Bmad::parser_get_integer(
    int& int_val,
    std::string& word,
    int& ix_word,
    std::string& delim,
    bool& delim_found,
    bool& err,
    optional_ref<std::string> str1,
    optional_ref<std::string> str2) {
  auto _word = word.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  const char* _str1 = str1.has_value() ? str1->get().c_str() : nullptr;
  const char* _str2 = str2.has_value() ? str2->get().c_str() : nullptr;
  fortran_parser_get_integer(
      /* int& */ int_val,
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err,
      /* const char* */ _str1,
      /* const char* */ _str2);
}
void Bmad::parser_get_logical(
    std::string& attrib_name,
    bool& this_logic,
    std::string& ele_name,
    std::string& delim,
    bool& delim_found,
    bool& err) {
  auto _attrib_name = attrib_name.c_str(); // ptr, inout, required
  auto _ele_name = ele_name.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parser_get_logical(
      /* const char* */ _attrib_name,
      /* bool& */ this_logic,
      /* const char* */ _ele_name,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err);
}
void Bmad::parser_identify_fork_to_element(LatProxy& lat) {
  fortran_parser_identify_fork_to_element(/* void* */ lat.get_fortran_ptr());
}
void Bmad::parser_init_custom_elements(LatProxy& lat) {
  fortran_parser_init_custom_elements(/* void* */ lat.get_fortran_ptr());
}
void Bmad::parser_print_line(LatProxy& lat, bool& end_of_file) {
  fortran_parser_print_line(
      /* void* */ lat.get_fortran_ptr(), /* bool& */ end_of_file);
}
void Bmad::parser_read_lr_wake(
    EleProxy& ele,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parser_read_lr_wake(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
void Bmad::parser_read_old_format_lr_wake(
    EleProxy& ele,
    std::string lr_file_name) {
  auto _lr_file_name = lr_file_name.c_str();
  fortran_parser_read_old_format_lr_wake(
      /* void* */ ele.get_fortran_ptr(), /* const char* */ _lr_file_name);
}
void Bmad::parser_read_old_format_sr_wake(
    EleProxy& ele,
    std::string sr_file_name) {
  auto _sr_file_name = sr_file_name.c_str();
  fortran_parser_read_old_format_sr_wake(
      /* void* */ ele.get_fortran_ptr(), /* const char* */ _sr_file_name);
}
void Bmad::parser_read_sr_wake(
    EleProxy& ele,
    std::string& delim,
    bool& delim_found,
    bool& err_flag) {
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_parser_read_sr_wake(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err_flag);
}
ControlProxy Bmad::parser_transfer_control_struct(
    ControlProxy& con_in,
    EleProxy& lord,
    int ix_var) {
  ControlProxy _con_out;
  fortran_parser_transfer_control_struct(
      /* void* */ con_in.get_fortran_ptr(),
      /* void* */ _con_out.get_fortran_ptr(),
      /* void* */ lord.get_fortran_ptr(),
      /* int& */ ix_var);
  return std::move(_con_out);
}
void Bmad::particle_in_global_frame(
    CoordProxy& orb,
    BranchProxy& branch,
    std::optional<bool> in_time_coordinates,
    std::optional<bool> in_body_frame,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat_out,
    CoordProxy& particle) {
  bool in_time_coordinates_lvalue;
  auto* _in_time_coordinates{&in_time_coordinates_lvalue};
  if (in_time_coordinates.has_value()) {
    in_time_coordinates_lvalue = in_time_coordinates.value();
  } else {
    _in_time_coordinates = nullptr;
  }
  bool in_body_frame_lvalue;
  auto* _in_body_frame{&in_body_frame_lvalue};
  if (in_body_frame.has_value()) {
    in_body_frame_lvalue = in_body_frame.value();
  } else {
    _in_body_frame = nullptr;
  }
  double _w_mat_out_vec[3 * 3];
  const double* _w_mat_out = nullptr;
  if (w_mat_out.has_value()) {
    matrix_to_vec(w_mat_out.value(), _w_mat_out_vec);
    _w_mat_out = _w_mat_out_vec;
  }
  fortran_particle_in_global_frame(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* bool* */ _in_time_coordinates,
      /* bool* */ _in_body_frame,
      /* double* */ _w_mat_out_vec,
      /* void* */ particle.get_fortran_ptr());
  if (w_mat_out.has_value())
    vec_to_matrix(_w_mat_out_vec, w_mat_out.value());
}
void Bmad::particle_is_moving_backwards(
    CoordProxy& orbit,
    bool& is_moving_backwards) {
  fortran_particle_is_moving_backwards(
      /* void* */ orbit.get_fortran_ptr(), /* bool& */ is_moving_backwards);
}
void Bmad::particle_is_moving_forward(
    CoordProxy& orbit,
    std::optional<int> dir,
    bool& is_moving_forward) {
  int dir_lvalue;
  auto* _dir{&dir_lvalue};
  if (dir.has_value()) {
    dir_lvalue = dir.value();
  } else {
    _dir = nullptr;
  }
  fortran_particle_is_moving_forward(
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _dir,
      /* bool& */ is_moving_forward);
}
void Bmad::particle_rf_time(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<bool> reference_active_edge,
    std::optional<double> s_rel,
    std::optional<bool> time_coords,
    std::optional<double> rf_freq,
    std::optional<bool> abs_time,
    long double& time) {
  bool reference_active_edge_lvalue;
  auto* _reference_active_edge{&reference_active_edge_lvalue};
  if (reference_active_edge.has_value()) {
    reference_active_edge_lvalue = reference_active_edge.value();
  } else {
    _reference_active_edge = nullptr;
  }
  double s_rel_lvalue;
  auto* _s_rel{&s_rel_lvalue};
  if (s_rel.has_value()) {
    s_rel_lvalue = s_rel.value();
  } else {
    _s_rel = nullptr;
  }
  bool time_coords_lvalue;
  auto* _time_coords{&time_coords_lvalue};
  if (time_coords.has_value()) {
    time_coords_lvalue = time_coords.value();
  } else {
    _time_coords = nullptr;
  }
  double rf_freq_lvalue;
  auto* _rf_freq{&rf_freq_lvalue};
  if (rf_freq.has_value()) {
    rf_freq_lvalue = rf_freq.value();
  } else {
    _rf_freq = nullptr;
  }
  bool abs_time_lvalue;
  auto* _abs_time{&abs_time_lvalue};
  if (abs_time.has_value()) {
    abs_time_lvalue = abs_time.value();
  } else {
    _abs_time = nullptr;
  }
  fortran_particle_rf_time(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _reference_active_edge,
      /* double* */ _s_rel,
      /* bool* */ _time_coords,
      /* double* */ _rf_freq,
      /* bool* */ _abs_time,
      /* long double& */ time);
}
void Bmad::patch_flips_propagation_direction(
    double x_pitch,
    double y_pitch,
    bool& is_flip) {
  fortran_patch_flips_propagation_direction(
      /* double& */ x_pitch, /* double& */ y_pitch, /* bool& */ is_flip);
}
void Bmad::patch_length(
    EleProxy& patch,
    std::optional<int> ref_coords,
    double& length) {
  int ref_coords_lvalue;
  auto* _ref_coords{&ref_coords_lvalue};
  if (ref_coords.has_value()) {
    ref_coords_lvalue = ref_coords.value();
  } else {
    _ref_coords = nullptr;
  }
  fortran_patch_length(
      /* void* */ patch.get_fortran_ptr(),
      /* int* */ _ref_coords,
      /* double& */ length);
}
Bmad::PhotonAbsorptionAndPhaseShift Bmad::photon_absorption_and_phase_shift(
    std::string material,
    double Energy) {
  auto _material = material.c_str();
  double _absorption{};
  double _phase_shift{};
  bool _err_flag{};
  fortran_photon_absorption_and_phase_shift(
      /* const char* */ _material,
      /* double& */ Energy,
      /* double& */ _absorption,
      /* double& */ _phase_shift,
      /* bool& */ _err_flag);
  return PhotonAbsorptionAndPhaseShift{_absorption, _phase_shift, _err_flag};
}
void Bmad::photon_add_to_detector_statistics(
    CoordProxy& orbit0,
    CoordProxy& orbit,
    EleProxy& ele,
    optional_ref<int> ix_pt,
    optional_ref<int> iy_pt,
    optional_ref<PixelPtProxy> pixel_pt) {
  auto* _ix_pt = ix_pt.has_value() ? &ix_pt->get() : nullptr; // inout, optional
  auto* _iy_pt = iy_pt.has_value() ? &iy_pt->get() : nullptr; // inout, optional
  auto* _pixel_pt = pixel_pt.has_value() ? pixel_pt->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  fortran_photon_add_to_detector_statistics(
      /* void* */ orbit0.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _ix_pt,
      /* int* */ _iy_pt,
      /* void* */ _pixel_pt);
}
Bmad::PhotonReflection Bmad::photon_reflection(
    double graze_angle_in,
    double energy,
    PhotonReflectSurfaceProxy& surface) {
  double _graze_angle_out{};
  double _phi_out{};
  fortran_photon_reflection(
      /* double& */ graze_angle_in,
      /* double& */ energy,
      /* void* */ surface.get_fortran_ptr(),
      /* double& */ _graze_angle_out,
      /* double& */ _phi_out);
  return PhotonReflection{_graze_angle_out, _phi_out};
}
PhotonReflectSurfaceProxy Bmad::photon_reflection_std_surface_init() {
  PhotonReflectSurfaceProxy _surface;
  fortran_photon_reflection_std_surface_init(
      /* void* */ _surface.get_fortran_ptr());
  return std::move(_surface);
}
Bmad::PhotonReflectivity Bmad::photon_reflectivity(
    double angle,
    double energy,
    PhotonReflectSurfaceProxy& surface) {
  double _p_reflect{};
  double _rel_p_specular{};
  fortran_photon_reflectivity(
      /* double& */ angle,
      /* double& */ energy,
      /* void* */ surface.get_fortran_ptr(),
      /* double& */ _p_reflect,
      /* double& */ _rel_p_specular);
  return PhotonReflectivity{_p_reflect, _rel_p_specular};
}
TargetPointProxy Bmad::photon_target_corner_calc(
    EleProxy& aperture_ele,
    double& x_lim,
    double& y_lim,
    double& z_lim,
    EleProxy& source_ele) {
  TargetPointProxy _corner;
  fortran_photon_target_corner_calc(
      /* void* */ aperture_ele.get_fortran_ptr(),
      /* double& */ x_lim,
      /* double& */ y_lim,
      /* double& */ z_lim,
      /* void* */ source_ele.get_fortran_ptr(),
      /* void* */ _corner.get_fortran_ptr());
  return std::move(_corner);
}
void Bmad::photon_target_setup(EleProxy& ele) {
  fortran_photon_target_setup(/* void* */ ele.get_fortran_ptr());
}
int Bmad::photon_type(EleProxy& ele) {
  int _e_type{};
  fortran_photon_type(/* void* */ ele.get_fortran_ptr(), /* int& */ _e_type);
  return _e_type;
}
void Bmad::physical_ele_end(
    int track_end,
    CoordProxy& orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end,
    int& physical_end) {
  bool return_stream_end_lvalue;
  auto* _return_stream_end{&return_stream_end_lvalue};
  if (return_stream_end.has_value()) {
    return_stream_end_lvalue = return_stream_end.value();
  } else {
    _return_stream_end = nullptr;
  }
  fortran_physical_ele_end(
      /* int& */ track_end,
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ ele_orientation,
      /* bool* */ _return_stream_end,
      /* int& */ physical_end);
}
void Bmad::point_photon_emission(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    int direction,
    double max_target_area,
    std::optional<FixedArray2D<Real, 3, 3>> w_to_surface) {
  double _w_to_surface_vec[3 * 3];
  const double* _w_to_surface = nullptr;
  if (w_to_surface.has_value()) {
    matrix_to_vec(w_to_surface.value(), _w_to_surface_vec);
    _w_to_surface = _w_to_surface_vec;
  }
  fortran_point_photon_emission(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ direction,
      /* double& */ max_target_area,
      /* double* */ _w_to_surface_vec);
}
BranchProxy Bmad::pointer_to_branch(EleProxy& ele) {
  BranchProxy _branch_ptr;
  fortran_pointer_to_branch_given_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _branch_ptr.get_fortran_ptr());
  return std::move(_branch_ptr);
}
BranchProxy Bmad::pointer_to_branch(
    std::string branch_name,
    LatProxy& lat,
    std::optional<bool> parameter_is_branch0,
    std::optional<int> blank_branch) {
  auto _branch_name = branch_name.c_str();
  bool parameter_is_branch0_lvalue;
  auto* _parameter_is_branch0{&parameter_is_branch0_lvalue};
  if (parameter_is_branch0.has_value()) {
    parameter_is_branch0_lvalue = parameter_is_branch0.value();
  } else {
    _parameter_is_branch0 = nullptr;
  }
  int blank_branch_lvalue;
  auto* _blank_branch{&blank_branch_lvalue};
  if (blank_branch.has_value()) {
    blank_branch_lvalue = blank_branch.value();
  } else {
    _blank_branch = nullptr;
  }
  BranchProxy _branch_ptr;
  fortran_pointer_to_branch_given_name(
      /* const char* */ _branch_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _parameter_is_branch0,
      /* int* */ _blank_branch,
      /* void* */ _branch_ptr.get_fortran_ptr());
  return std::move(_branch_ptr);
}
EleProxy Bmad::pointer_to_ele(
    LatProxy& lat,
    int ix_ele,
    std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  EleProxy _ele_ptr;
  fortran_pointer_to_ele1(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ ix_ele,
      /* int* */ _ix_branch,
      /* void* */ _ele_ptr.get_fortran_ptr());
  return std::move(_ele_ptr);
}
EleProxy Bmad::pointer_to_ele(LatProxy& lat, LatEleLocProxy& ele_loc) {
  EleProxy _ele_ptr;
  fortran_pointer_to_ele2(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ ele_loc.get_fortran_ptr(),
      /* void* */ _ele_ptr.get_fortran_ptr());
  return std::move(_ele_ptr);
}
EleProxy Bmad::pointer_to_ele(LatProxy& lat, std::string ele_name) {
  auto _ele_name = ele_name.c_str();
  EleProxy _ele_ptr;
  fortran_pointer_to_ele3(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _ele_name,
      /* void* */ _ele_ptr.get_fortran_ptr());
  return std::move(_ele_ptr);
}
EleProxy Bmad::pointer_to_ele(LatProxy& lat, EleProxy& foreign_ele) {
  EleProxy _ele_ptr;
  fortran_pointer_to_ele4(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ foreign_ele.get_fortran_ptr(),
      /* void* */ _ele_ptr.get_fortran_ptr());
  return std::move(_ele_ptr);
}
Bmad::PointerToElementAtS Bmad::pointer_to_element_at_s(
    BranchProxy& branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err) {
  bool _err_flag{};
  double _s_eff{};
  CoordProxy _position;
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  EleProxy _ele;
  fortran_pointer_to_element_at_s(
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ s,
      /* bool& */ choose_max,
      /* bool& */ _err_flag,
      /* double& */ _s_eff,
      /* void* */ _position.get_fortran_ptr(),
      /* bool* */ _print_err,
      /* void* */ _ele.get_fortran_ptr());
  return PointerToElementAtS{
      _err_flag, _s_eff, std::move(_position), std::move(_ele)};
}
double Bmad::pointer_to_field_ele(
    EleProxy& ele,
    int ix_field_ele,
    EleProxy& field_ele) {
  double _dz_offset{};
  auto _field_ele = &field_ele; // input, required, pointer
  fortran_pointer_to_field_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ ix_field_ele,
      /* double& */ _dz_offset,
      /* void* */ &field_ele);
  return _dz_offset;
}
int Bmad::pointer_to_girder(EleProxy& ele, EleProxy& girder) {
  int _ix_slave_back{};
  auto _girder = &girder; // input, required, pointer
  fortran_pointer_to_girder(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ _ix_slave_back,
      /* void* */ &girder);
  return _ix_slave_back;
}
Bmad::PointerToLord Bmad::pointer_to_lord(
    EleProxy& slave,
    int ix_lord,
    std::optional<int> lord_type,
    EleProxy& lord_ptr) {
  ControlProxy _control;
  int _ix_slave_back{};
  int lord_type_lvalue;
  auto* _lord_type{&lord_type_lvalue};
  if (lord_type.has_value()) {
    lord_type_lvalue = lord_type.value();
  } else {
    _lord_type = nullptr;
  }
  int _ix_control{};
  int _ix_ic{};
  auto _lord_ptr = &lord_ptr; // input, required, pointer
  fortran_pointer_to_lord(
      /* void* */ slave.get_fortran_ptr(),
      /* int& */ ix_lord,
      /* void* */ _control.get_fortran_ptr(),
      /* int& */ _ix_slave_back,
      /* int* */ _lord_type,
      /* int& */ _ix_control,
      /* int& */ _ix_ic,
      /* void* */ &lord_ptr);
  return PointerToLord{
      std::move(_control), _ix_slave_back, _ix_control, _ix_ic};
}
Bmad::PointerToMultipassLord Bmad::pointer_to_multipass_lord(
    EleProxy& ele,
    EleProxy& multi_lord) {
  int _ix_pass{};
  EleProxy _super_lord;
  auto _multi_lord = &multi_lord; // input, required, pointer
  fortran_pointer_to_multipass_lord(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ _ix_pass,
      /* void* */ _super_lord.get_fortran_ptr(),
      /* void* */ &multi_lord);
  return PointerToMultipassLord{_ix_pass, std::move(_super_lord)};
}
void Bmad::pointer_to_next_ele(
    EleProxy& this_ele,
    std::optional<int> offset,
    std::optional<bool> skip_beginning,
    std::optional<bool> follow_fork,
    EleProxy& next_ele) {
  int offset_lvalue;
  auto* _offset{&offset_lvalue};
  if (offset.has_value()) {
    offset_lvalue = offset.value();
  } else {
    _offset = nullptr;
  }
  bool skip_beginning_lvalue;
  auto* _skip_beginning{&skip_beginning_lvalue};
  if (skip_beginning.has_value()) {
    skip_beginning_lvalue = skip_beginning.value();
  } else {
    _skip_beginning = nullptr;
  }
  bool follow_fork_lvalue;
  auto* _follow_fork{&follow_fork_lvalue};
  if (follow_fork.has_value()) {
    follow_fork_lvalue = follow_fork.value();
  } else {
    _follow_fork = nullptr;
  }
  auto _next_ele = &next_ele; // input, required, pointer
  fortran_pointer_to_next_ele(
      /* void* */ this_ele.get_fortran_ptr(),
      /* int* */ _offset,
      /* bool* */ _skip_beginning,
      /* bool* */ _follow_fork,
      /* void* */ &next_ele);
}
Bmad::PointerToSlave Bmad::pointer_to_slave(
    EleProxy& lord,
    int ix_slave,
    std::optional<int> lord_type) {
  ControlProxy _control;
  int lord_type_lvalue;
  auto* _lord_type{&lord_type_lvalue};
  if (lord_type.has_value()) {
    lord_type_lvalue = lord_type.value();
  } else {
    _lord_type = nullptr;
  }
  int _ix_lord_back{};
  int _ix_control{};
  int _ix_ic{};
  EleProxy _slave_ptr;
  fortran_pointer_to_slave(
      /* void* */ lord.get_fortran_ptr(),
      /* int& */ ix_slave,
      /* void* */ _control.get_fortran_ptr(),
      /* int* */ _lord_type,
      /* int& */ _ix_lord_back,
      /* int& */ _ix_control,
      /* int& */ _ix_ic,
      /* void* */ _slave_ptr.get_fortran_ptr());
  return PointerToSlave{
      std::move(_control),
      _ix_lord_back,
      _ix_control,
      _ix_ic,
      std::move(_slave_ptr)};
}
Bmad::PointerToSuperLord Bmad::pointer_to_super_lord(
    EleProxy& slave,
    std::optional<int> lord_type,
    EleProxy& lord_ptr) {
  ControlProxy _control;
  int _ix_slave_back{};
  int _ix_control{};
  int _ix_ic{};
  int lord_type_lvalue;
  auto* _lord_type{&lord_type_lvalue};
  if (lord_type.has_value()) {
    lord_type_lvalue = lord_type.value();
  } else {
    _lord_type = nullptr;
  }
  auto _lord_ptr = &lord_ptr; // input, required, pointer
  fortran_pointer_to_super_lord(
      /* void* */ slave.get_fortran_ptr(),
      /* void* */ _control.get_fortran_ptr(),
      /* int& */ _ix_slave_back,
      /* int& */ _ix_control,
      /* int& */ _ix_ic,
      /* int* */ _lord_type,
      /* void* */ &lord_ptr);
  return PointerToSuperLord{
      std::move(_control), _ix_slave_back, _ix_control, _ix_ic};
}
SurfaceDisplacementPtProxy Bmad::pointer_to_surface_displacement_pt(
    EleProxy& ele,
    bool nearest,
    double& x,
    double& y,
    optional_ref<int> ix,
    optional_ref<int> iy,
    std::optional<bool> extend_grid,
    optional_ref<double> xx,
    optional_ref<double> yy) {
  auto* _ix = ix.has_value() ? &ix->get() : nullptr; // inout, optional
  auto* _iy = iy.has_value() ? &iy->get() : nullptr; // inout, optional
  bool extend_grid_lvalue;
  auto* _extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  auto* _xx = xx.has_value() ? &xx->get() : nullptr; // inout, optional
  auto* _yy = yy.has_value() ? &yy->get() : nullptr; // inout, optional
  SurfaceDisplacementPtProxy _pt;
  fortran_pointer_to_surface_displacement_pt(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ nearest,
      /* double& */ x,
      /* double& */ y,
      /* int* */ _ix,
      /* int* */ _iy,
      /* bool* */ _extend_grid,
      /* double* */ _xx,
      /* double* */ _yy,
      /* void* */ _pt.get_fortran_ptr());
  return std::move(_pt);
}
SurfaceSegmentedPtProxy Bmad::pointer_to_surface_segmented_pt(
    EleProxy& ele,
    bool nearest,
    double& x,
    double& y,
    optional_ref<int> ix,
    optional_ref<int> iy,
    std::optional<bool> extend_grid,
    optional_ref<double> xx,
    optional_ref<double> yy) {
  auto* _ix = ix.has_value() ? &ix->get() : nullptr; // inout, optional
  auto* _iy = iy.has_value() ? &iy->get() : nullptr; // inout, optional
  bool extend_grid_lvalue;
  auto* _extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  auto* _xx = xx.has_value() ? &xx->get() : nullptr; // inout, optional
  auto* _yy = yy.has_value() ? &yy->get() : nullptr; // inout, optional
  SurfaceSegmentedPtProxy _pt;
  fortran_pointer_to_surface_segmented_pt(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ nearest,
      /* double& */ x,
      /* double& */ y,
      /* int* */ _ix,
      /* int* */ _iy,
      /* bool* */ _extend_grid,
      /* double* */ _xx,
      /* double* */ _yy,
      /* void* */ _pt.get_fortran_ptr());
  return std::move(_pt);
}
double Bmad::pointer_to_wake_ele(EleProxy& ele, EleProxy& wake_ele) {
  double _delta_s{};
  auto _wake_ele = &wake_ele; // input, required, pointer
  fortran_pointer_to_wake_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _delta_s,
      /* void* */ &wake_ele);
  return _delta_s;
}
Bmad::PointerToWall3d Bmad::pointer_to_wall3d(
    EleProxy& ele,
    std::optional<int> ix_wall) {
  int ix_wall_lvalue;
  auto* _ix_wall{&ix_wall_lvalue};
  if (ix_wall.has_value()) {
    ix_wall_lvalue = ix_wall.value();
  } else {
    _ix_wall = nullptr;
  }
  double _ds_offset{};
  bool _is_branch_wall{};
  Wall3dProxy _wall3d;
  fortran_pointer_to_wall3d(
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _ix_wall,
      /* double& */ _ds_offset,
      /* bool& */ _is_branch_wall,
      /* void* */ _wall3d.get_fortran_ptr());
  return PointerToWall3d{_ds_offset, _is_branch_wall, std::move(_wall3d)};
}
void Bmad::polar_to_spinor(
    SpinPolarProxy& polar,
    FixedArray1D<Complex, 2> spinor) {
  auto* _spinor = spinor.data(); // CppWrapperGeneralArgument
  fortran_polar_to_spinor(
      /* void* */ polar.get_fortran_ptr(), /* std::complex<double>* */ _spinor);
}
void Bmad::polar_to_vec(SpinPolarProxy& polar, FixedArray1D<Real, 3> vec) {
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  fortran_polar_to_vec(/* void* */ polar.get_fortran_ptr(), /* double* */ _vec);
}
Bmad::ProjectEmitToXyz Bmad::project_emit_to_xyz(
    LatProxy& ring,
    int ix,
    NormalModesProxy& mode) {
  double _sigma_x{};
  double _sigma_y{};
  double _sigma_z{};
  fortran_project_emit_to_xyz(
      /* void* */ ring.get_fortran_ptr(),
      /* int& */ ix,
      /* void* */ mode.get_fortran_ptr(),
      /* double& */ _sigma_x,
      /* double& */ _sigma_y,
      /* double& */ _sigma_z);
  return ProjectEmitToXyz{_sigma_x, _sigma_y, _sigma_z};
}
double Bmad::psi_prime_sca(double t, double p, FixedArray1D<Real, 8> args) {
  double _dpdt{};
  auto* _args = args.data(); // CppWrapperGeneralArgument
  fortran_psi_prime_sca(
      /* double& */ t,
      /* double& */ p,
      /* double& */ _dpdt,
      /* double* */ _args);
  return _dpdt;
}
void Bmad::ptc_bookkeeper(LatProxy& lat) {
  fortran_ptc_bookkeeper(/* void* */ lat.get_fortran_ptr());
}
CoordProxyAlloc1D Bmad::ptc_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<bool> radiation_damping_on) {
  // intent=out allocatable type array
  auto closed_orbit{CoordProxyAlloc1D()};
  bool radiation_damping_on_lvalue;
  auto* _radiation_damping_on{&radiation_damping_on_lvalue};
  if (radiation_damping_on.has_value()) {
    radiation_damping_on_lvalue = radiation_damping_on.value();
  } else {
    _radiation_damping_on = nullptr;
  }
  fortran_ptc_closed_orbit_calc(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ closed_orbit.get_fortran_ptr(),
      /* bool* */ _radiation_damping_on);
  return std::move(closed_orbit);
}
Bmad::PtcEmitCalc Bmad::ptc_emit_calc(
    EleProxy& ele,
    FixedArray2D<Real, 6, 6> sigma_mat) {
  NormalModesProxy _norm_mode;
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  CoordProxy _closed_orb;
  fortran_ptc_emit_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _norm_mode.get_fortran_ptr(),
      /* double* */ _sigma_mat_vec,
      /* void* */ _closed_orb.get_fortran_ptr());
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  return PtcEmitCalc{std::move(_norm_mode), std::move(_closed_orb)};
}
void Bmad::ptc_layouts_resplit(
    double dKL_max,
    double l_max,
    bool l_max_drift_only,
    double bend_dorb,
    double sex_dx,
    std::optional<bool> even,
    std::optional<FixedArray1D<Int, 2>> crossover,
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler) {
  bool even_lvalue;
  auto* _even{&even_lvalue};
  if (even.has_value()) {
    even_lvalue = even.value();
  } else {
    _even = nullptr;
  }
  int* _crossover = crossover.has_value() ? crossover.value().data() : nullptr;
  int* _crossover_wiggler = crossover_wiggler.has_value()
      ? crossover_wiggler.value().data()
      : nullptr;
  fortran_ptc_layouts_resplit(
      /* double& */ dKL_max,
      /* double& */ l_max,
      /* bool& */ l_max_drift_only,
      /* double& */ bend_dorb,
      /* double& */ sex_dx,
      /* bool* */ _even,
      /* int* */ _crossover,
      /* int* */ _crossover_wiggler);
}
void Bmad::ptc_one_turn_mat_and_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<double> pz) {
  double pz_lvalue;
  auto* _pz{&pz_lvalue};
  if (pz.has_value()) {
    pz_lvalue = pz.value();
  } else {
    _pz = nullptr;
  }
  fortran_ptc_one_turn_mat_and_closed_orbit_calc(
      /* void* */ branch.get_fortran_ptr(), /* double* */ _pz);
}
void Bmad::ptc_ran_seed_put(int iseed) {
  fortran_ptc_ran_seed_put(/* int& */ iseed);
}
void Bmad::ptc_set_rf_state_for_c_normal(bool nocavity) {
  fortran_ptc_set_rf_state_for_c_normal(/* bool& */ nocavity);
}
void Bmad::ptc_set_taylor_order_if_needed() {
  fortran_ptc_set_taylor_order_if_needed();
}
Bmad::PtcSpinCalc Bmad::ptc_spin_calc(
    EleProxy& ele,
    FixedArray2D<Real, 6, 6> sigma_mat) {
  NormalModesProxy _norm_mode;
  double _sigma_mat_vec[6 * 6];
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  CoordProxy _closed_orb;
  fortran_ptc_spin_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _norm_mode.get_fortran_ptr(),
      /* double* */ _sigma_mat_vec,
      /* void* */ _closed_orb.get_fortran_ptr());
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  return PtcSpinCalc{std::move(_norm_mode), std::move(_closed_orb)};
}
Bmad::PtcTrackAll Bmad::ptc_track_all(
    BranchProxy& branch,
    CoordProxyAlloc1D& orbit) {
  // intent=inout allocatable type array
  int _track_state{};
  bool _err_flag{};
  fortran_ptc_track_all(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ _track_state,
      /* bool& */ _err_flag);
  return PtcTrackAll{_track_state, _err_flag};
}
bool Bmad::ptc_transfer_map_with_spin(
    BranchProxy& branch,
    FixedArray1D<TaylorProxy, 6> t_map,
    FixedArray1D<TaylorProxy, 4> s_map,
    CoordProxy& orb0,
    std::optional<int> ix1,
    std::optional<int> ix2,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start) {
  bool _err_flag{};
  int ix1_lvalue;
  auto* _ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto* _ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  bool one_turn_lvalue;
  auto* _one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto* _unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  fortran_ptc_transfer_map_with_spin(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ t_map.data(),
      /* void* */ s_map.data(),
      /* void* */ orb0.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* int* */ _ix1,
      /* int* */ _ix2,
      /* bool* */ _one_turn,
      /* bool* */ _unit_start);
  return _err_flag;
}
FixedArray2D<Real, 6, 6> Bmad::pwd_mat(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> t6,
    double inductance,
    double sig_z) {
  double _t6_vec[6 * 6];
  matrix_to_vec(t6, _t6_vec);
  FixedArray2D<Real, 6, 6> t6_pwd;
  double _t6_pwd_vec[6 * 6];
  fortran_pwd_mat(
      /* void* */ lat.get_fortran_ptr(),
      /* double* */ _t6_vec,
      /* double& */ inductance,
      /* double& */ sig_z,
      /* double* */ _t6_pwd_vec);
  vec_to_matrix(_t6_pwd_vec, t6_pwd);
  return t6_pwd;
}
Bmad::Rad1DampAndStocMats Bmad::rad1_damp_and_stoc_mats(
    EleProxy& ele,
    bool include_opening_angle,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double g2_tol,
    double g3_tol,
    optional_ref<EleProxy> ele0) {
  RadMapProxy _rad_map;
  bool _err_flag{};
  auto* _ele0 = ele0.has_value() ? ele0->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  RadInt1Proxy _rad_int1;
  fortran_rad1_damp_and_stoc_mats(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ include_opening_angle,
      /* void* */ orb_in.get_fortran_ptr(),
      /* void* */ orb_out.get_fortran_ptr(),
      /* void* */ _rad_map.get_fortran_ptr(),
      /* double& */ g2_tol,
      /* double& */ g3_tol,
      /* bool& */ _err_flag,
      /* void* */ _ele0,
      /* void* */ _rad_int1.get_fortran_ptr());
  return Rad1DampAndStocMats{
      std::move(_rad_map), _err_flag, std::move(_rad_int1)};
}
Bmad::RadDampAndStocMats Bmad::rad_damp_and_stoc_mats(
    EleProxy& ele1,
    EleProxy& ele2,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit) {
  RadMapProxy _rmap;
  NormalModesProxy _mode;
  FixedArray2D<Real, 6, 6> xfer_nodamp_mat;
  double _xfer_nodamp_mat_vec[6 * 6];
  bool _err_flag{};
  // intent=in allocatable type array
  auto* _closed_orbit = closed_orbit.has_value()
      ? closed_orbit->get().get_fortran_ptr()
      : nullptr; // input, optional
  RadIntBranchProxy _rad_int_branch;
  fortran_rad_damp_and_stoc_mats(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool& */ include_opening_angle,
      /* void* */ _rmap.get_fortran_ptr(),
      /* void* */ _mode.get_fortran_ptr(),
      /* double* */ _xfer_nodamp_mat_vec,
      /* bool& */ _err_flag,
      /* void* */ _closed_orbit,
      /* void* */ _rad_int_branch.get_fortran_ptr());
  vec_to_matrix(_xfer_nodamp_mat_vec, xfer_nodamp_mat);
  return RadDampAndStocMats{
      std::move(_rmap),
      std::move(_mode),
      xfer_nodamp_mat,
      _err_flag,
      std::move(_rad_int_branch)};
}
FixedArray1D<Real, 2> Bmad::rad_g_integrals(
    EleProxy& ele,
    int where,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double& int_g2,
    double& int_g3,
    double g_tol,
    double g2_tol,
    double g3_tol) {
  FixedArray1D<Real, 2> _int_g;
  fortran_rad_g_integrals(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ where,
      /* void* */ orb_in.get_fortran_ptr(),
      /* void* */ orb_out.get_fortran_ptr(),
      /* double* */ _int_g.data(),
      /* double& */ int_g2,
      /* double& */ int_g3,
      /* double& */ g_tol,
      /* double& */ g2_tol,
      /* double& */ g3_tol);
  return _int_g;
}
Bmad::RadiationIntegrals Bmad::radiation_integrals(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    optional_ref<int> ix_cache,
    std::optional<int> ix_branch) {
  // intent=in allocatable type array
  NormalModesProxy _mode;
  auto* _ix_cache =
      ix_cache.has_value() ? &ix_cache->get() : nullptr; // inout, optional
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  RadIntAllEleProxy _rad_int_by_ele;
  fortran_radiation_integrals(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ _mode.get_fortran_ptr(),
      /* int* */ _ix_cache,
      /* int* */ _ix_branch,
      /* void* */ _rad_int_by_ele.get_fortran_ptr());
  return RadiationIntegrals{std::move(_mode), std::move(_rad_int_by_ele)};
}
bool Bmad::radiation_map_setup(
    EleProxy& ele,
    optional_ref<CoordProxy> ref_orbit_in) {
  bool _err_flag{};
  auto* _ref_orbit_in = ref_orbit_in.has_value()
      ? ref_orbit_in->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_radiation_map_setup(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _ref_orbit_in);
  return _err_flag;
}
void Bmad::ramper_slave_setup(LatProxy& lat, std::optional<bool> force_setup) {
  bool force_setup_lvalue;
  auto* _force_setup{&force_setup_lvalue};
  if (force_setup.has_value()) {
    force_setup_lvalue = force_setup.value();
  } else {
    _force_setup = nullptr;
  }
  fortran_ramper_slave_setup(
      /* void* */ lat.get_fortran_ptr(), /* bool* */ _force_setup);
}
bool Bmad::ramper_value(
    EleProxy& ramper,
    ControlRamp1Proxy& r1,
    double& value) {
  bool _err_flag{};
  fortran_ramper_value(
      /* void* */ ramper.get_fortran_ptr(),
      /* void* */ r1.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* double& */ value);
  return _err_flag;
}
bool Bmad::randomize_lr_wake_frequencies(EleProxy& ele) {
  bool _set_done{};
  fortran_randomize_lr_wake_frequencies(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _set_done);
  return _set_done;
}
void Bmad::rchomp(double& rel, int& plc, std::string& out) {
  auto _out = out.c_str(); // ptr, inout, required
  fortran_rchomp(/* double& */ rel, /* int& */ plc, /* const char* */ _out);
}
void Bmad::re_allocate_eles(
    ElePointerProxyAlloc1D& eles,
    int n,
    std::optional<bool> save_old,
    std::optional<bool> exact) {
  // intent=inout allocatable type array
  bool save_old_lvalue;
  auto* _save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_allocate_eles(
      /* void* */ eles.get_fortran_ptr(),
      /* int& */ n,
      /* bool* */ _save_old,
      /* bool* */ _exact);
}
void Bmad::re_allocate(
    Wall3dSectionProxyAlloc1D& section,
    int n,
    optional_ref<bool> exact) {
  // intent=inout allocatable type array
  auto* _exact = exact.has_value() ? &exact->get() : nullptr; // inout, optional
  fortran_re_allocate_wall3d_section_array(
      /* void* */ section.get_fortran_ptr(), /* int& */ n, /* bool* */ _exact);
}
void Bmad::re_allocate(
    Wall3dVertexProxyAlloc1D& v,
    int n,
    optional_ref<bool> exact) {
  // intent=inout allocatable type array
  auto* _exact = exact.has_value() ? &exact->get() : nullptr; // inout, optional
  fortran_re_allocate_wall3d_vertex_array(
      /* void* */ v.get_fortran_ptr(), /* int& */ n, /* bool* */ _exact);
}
void Bmad::re_associate_node_array(
    ExpressionTreeProxy& tree,
    int n,
    std::optional<bool> exact) {
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_associate_node_array(
      /* void* */ tree.get_fortran_ptr(), /* int& */ n, /* bool* */ _exact);
}
void Bmad::re_str(long double& rel, std::string& str_out) {
  auto _str_out = str_out.c_str(); // ptr, inout, required
  fortran_re_str_qp(/* long double& */ rel, /* const char* */ _str_out);
}
void Bmad::re_str(double& rel, std::string& str_out) {
  auto _str_out = str_out.c_str(); // ptr, inout, required
  fortran_re_str_rp(/* double& */ rel, /* const char* */ _str_out);
}
Bmad::ReadBeamAscii Bmad::read_beam_ascii(
    std::string file_name,
    BeamInitProxy& beam_init) {
  auto _file_name = file_name.c_str();
  BeamProxy _beam;
  bool _err_flag{};
  fortran_read_beam_ascii(
      /* const char* */ _file_name,
      /* void* */ _beam.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return ReadBeamAscii{std::move(_beam), _err_flag};
}
Bmad::ReadBeamFile Bmad::read_beam_file(
    std::string file_name,
    BeamInitProxy& beam_init,
    optional_ref<EleProxy> ele,
    std::optional<bool> print_mom_shift_warning,
    optional_ref<bool> conserve_momentum) {
  auto _file_name = file_name.c_str();
  BeamProxy _beam;
  bool _err_flag{};
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool print_mom_shift_warning_lvalue;
  auto* _print_mom_shift_warning{&print_mom_shift_warning_lvalue};
  if (print_mom_shift_warning.has_value()) {
    print_mom_shift_warning_lvalue = print_mom_shift_warning.value();
  } else {
    _print_mom_shift_warning = nullptr;
  }
  auto* _conserve_momentum = conserve_momentum.has_value()
      ? &conserve_momentum->get()
      : nullptr; // inout, optional
  fortran_read_beam_file(
      /* const char* */ _file_name,
      /* void* */ _beam.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _ele,
      /* bool* */ _print_mom_shift_warning,
      /* bool* */ _conserve_momentum);
  return ReadBeamFile{std::move(_beam), _err_flag};
}
void Bmad::read_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_read_binary_cartesian_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cart_map.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::read_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_read_binary_cylindrical_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cl_map.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::read_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_read_binary_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ g_field.get_fortran_ptr(),
      /* bool& */ err_flag);
}
PhotonReflectSurfaceProxy Bmad::read_surface_reflection_file(
    std::string file_name) {
  auto _file_name = file_name.c_str();
  PhotonReflectSurfaceProxy _surface;
  fortran_read_surface_reflection_file(
      /* const char* */ _file_name, /* void* */ _surface.get_fortran_ptr());
  return std::move(_surface);
}
void Bmad::reallocate_beam(
    BeamProxy& beam,
    int n_bunch,
    std::optional<int> n_particle,
    optional_ref<bool> extend) {
  int n_particle_lvalue;
  auto* _n_particle{&n_particle_lvalue};
  if (n_particle.has_value()) {
    n_particle_lvalue = n_particle.value();
  } else {
    _n_particle = nullptr;
  }
  auto* _extend =
      extend.has_value() ? &extend->get() : nullptr; // inout, optional
  fortran_reallocate_beam(
      /* void* */ beam.get_fortran_ptr(),
      /* int& */ n_bunch,
      /* int* */ _n_particle,
      /* bool* */ _extend);
}
void Bmad::reallocate_bp_com_const() {
  fortran_reallocate_bp_com_const();
}
BunchProxy Bmad::reallocate_bunch(int n_particle, std::optional<bool> save) {
  BunchProxy _bunch;
  bool save_lvalue;
  auto* _save{&save_lvalue};
  if (save.has_value()) {
    save_lvalue = save.value();
  } else {
    _save = nullptr;
  }
  fortran_reallocate_bunch(
      /* void* */ _bunch.get_fortran_ptr(),
      /* int& */ n_particle,
      /* bool* */ _save);
  return std::move(_bunch);
}
void Bmad::reallocate_control(LatProxy& lat, int n) {
  fortran_reallocate_control(/* void* */ lat.get_fortran_ptr(), /* int& */ n);
}
void Bmad::reallocate_coord(
    CoordArrayProxyAlloc1D& coord_array,
    LatProxy& lat) {
  // intent=inout allocatable type array
  fortran_reallocate_coord_array(
      /* void* */ coord_array.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr());
}
void Bmad::reallocate_coord(
    CoordProxyAlloc1D& coord,
    LatProxy& lat,
    std::optional<int> ix_branch) {
  // intent=inout allocatable type array
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  fortran_reallocate_coord_lat(
      /* void* */ coord.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _ix_branch);
}
void Bmad::reallocate_coord(CoordProxyAlloc1D& coord, int n_coord) {
  // intent=inout allocatable type array
  fortran_reallocate_coord_n(
      /* void* */ coord.get_fortran_ptr(), /* int& */ n_coord);
}
void Bmad::reallocate_expression_stack(
    ExpressionAtomProxyAlloc1D& stack,
    int n,
    std::optional<bool> exact) {
  // intent=inout allocatable type array
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_reallocate_expression_stack(
      /* void* */ stack.get_fortran_ptr(), /* int& */ n, /* bool* */ _exact);
}
void Bmad::rel_tracking_charge_to_mass(
    CoordProxy& orbit,
    int ref_species,
    double& rel_charge) {
  fortran_rel_tracking_charge_to_mass(
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ ref_species,
      /* double& */ rel_charge);
}
void Bmad::relative_mode_flip(
    EleProxy& ele1,
    EleProxy& ele2,
    bool& func_retval__) {
  fortran_relative_mode_flip(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool& */ func_retval__);
}
void Bmad::release_rad_int_cache(int& ix_cache) {
  fortran_release_rad_int_cache(/* int& */ ix_cache);
}
Bmad::RemoveConstantTaylor Bmad::remove_constant_taylor(
    TaylorProxyAlloc1D& taylor_in,
    bool remove_higher_order_terms) {
  // intent=in allocatable type array
  // intent=out allocatable type array
  auto taylor_out{TaylorProxyAlloc1D()};
  // intent=out allocatable general array
  auto c0{RealAlloc1D()};
  fortran_remove_constant_taylor(
      /* void* */ taylor_in.get_fortran_ptr(),
      /* void* */ taylor_out.get_fortran_ptr(),
      /* void* */ c0.get_fortran_ptr(),
      /* bool& */ remove_higher_order_terms);
  return RemoveConstantTaylor{std::move(taylor_out), std::move(c0)};
}
BunchProxy Bmad::remove_dead_from_bunch(BunchProxy& bunch_in) {
  BunchProxy _bunch_out;
  fortran_remove_dead_from_bunch(
      /* void* */ bunch_in.get_fortran_ptr(),
      /* void* */ _bunch_out.get_fortran_ptr());
  return std::move(_bunch_out);
}
void Bmad::remove_eles_from_lat(
    LatProxy& lat,
    std::optional<bool> check_sanity) {
  bool check_sanity_lvalue;
  auto* _check_sanity{&check_sanity_lvalue};
  if (check_sanity.has_value()) {
    check_sanity_lvalue = check_sanity.value();
  } else {
    _check_sanity = nullptr;
  }
  fortran_remove_eles_from_lat(
      /* void* */ lat.get_fortran_ptr(), /* bool* */ _check_sanity);
}
void Bmad::remove_lord_slave_link(EleProxy& lord, EleProxy& slave) {
  fortran_remove_lord_slave_link(
      /* void* */ lord.get_fortran_ptr(), /* void* */ slave.get_fortran_ptr());
}
LatProxy Bmad::reverse_lat(
    LatProxy& lat_in,
    std::optional<bool> track_antiparticle) {
  LatProxy _lat_rev;
  bool track_antiparticle_lvalue;
  auto* _track_antiparticle{&track_antiparticle_lvalue};
  if (track_antiparticle.has_value()) {
    track_antiparticle_lvalue = track_antiparticle.value();
  } else {
    _track_antiparticle = nullptr;
  }
  fortran_reverse_lat(
      /* void* */ lat_in.get_fortran_ptr(),
      /* void* */ _lat_rev.get_fortran_ptr(),
      /* bool* */ _track_antiparticle);
  return std::move(_lat_rev);
}
void Bmad::rf_coupler_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    double phase,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_rf_coupler_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* double& */ phase,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::rf_is_on(
    BranchProxy& branch,
    std::optional<int> ix_ele1,
    std::optional<int> ix_ele2,
    bool& is_on) {
  int ix_ele1_lvalue;
  auto* _ix_ele1{&ix_ele1_lvalue};
  if (ix_ele1.has_value()) {
    ix_ele1_lvalue = ix_ele1.value();
  } else {
    _ix_ele1 = nullptr;
  }
  int ix_ele2_lvalue;
  auto* _ix_ele2{&ix_ele2_lvalue};
  if (ix_ele2.has_value()) {
    ix_ele2_lvalue = ix_ele2.value();
  } else {
    _ix_ele2 = nullptr;
  }
  fortran_rf_is_on(
      /* void* */ branch.get_fortran_ptr(),
      /* int* */ _ix_ele1,
      /* int* */ _ix_ele2,
      /* bool& */ is_on);
}
void Bmad::rf_ref_time_offset(
    EleProxy& ele,
    std::optional<double> ds,
    double& time) {
  double ds_lvalue;
  auto* _ds{&ds_lvalue};
  if (ds.has_value()) {
    ds_lvalue = ds.value();
  } else {
    _ds = nullptr;
  }
  fortran_rf_ref_time_offset(
      /* void* */ ele.get_fortran_ptr(), /* double* */ _ds, /* double& */ time);
}
void Bmad::rfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& a,
    double& b,
    double& hz,
    int& i,
    int& j,
    double& res) {
  fortran_rfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ a,
      /* double& */ b,
      /* double& */ hz,
      /* int& */ i,
      /* int& */ j,
      /* double& */ res);
}
void Bmad::rk_adaptive_time_step(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orb,
    int& t_dir,
    double& rf_time,
    double& dt_try,
    double& dt_did,
    double& dt_next,
    bool& err_flag,
    optional_ref<EmFieldProxy> extra_field) {
  auto* _extra_field = extra_field.has_value()
      ? extra_field->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_rk_adaptive_time_step(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orb.get_fortran_ptr(),
      /* int& */ t_dir,
      /* double& */ rf_time,
      /* double& */ dt_try,
      /* double& */ dt_did,
      /* double& */ dt_next,
      /* bool& */ err_flag,
      /* void* */ _extra_field);
}
FixedArray1D<Real, 10> Bmad::rk_time_step1(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orb,
    double dt,
    CoordProxy& new_orb,
    std::optional<FixedArray1D<Real, 10>> dr_dt,
    bool& err_flag,
    optional_ref<bool> print_err,
    optional_ref<EmFieldProxy> extra_field) {
  FixedArray1D<Real, 10> _r_err;
  double* _dr_dt = dr_dt.has_value() ? dr_dt.value().data() : nullptr;
  auto* _print_err =
      print_err.has_value() ? &print_err->get() : nullptr; // inout, optional
  auto* _extra_field = extra_field.has_value()
      ? extra_field->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_rk_time_step1(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ rf_time,
      /* void* */ orb.get_fortran_ptr(),
      /* double& */ dt,
      /* void* */ new_orb.get_fortran_ptr(),
      /* double* */ _r_err.data(),
      /* double* */ _dr_dt,
      /* bool& */ err_flag,
      /* bool* */ _print_err,
      /* void* */ _extra_field);
  return _r_err;
}
void Bmad::rotate3(
    FixedArray1D<Real, 3> vec,
    double& angle,
    FixedArray1D<Real, 3> rvec) {
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  auto* _rvec = rvec.data(); // CppWrapperGeneralArgument
  fortran_rotate3(/* double* */ _vec, /* double& */ angle, /* double* */ _rvec);
}
void Bmad::rotate_em_field(
    EmFieldProxy& field,
    FixedArray2D<Real, 3, 3> w_mat,
    FixedArray2D<Real, 3, 3> w_inv,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  double _w_inv_vec[3 * 3];
  matrix_to_vec(w_inv, _w_inv_vec);
  bool calc_dfield_lvalue;
  auto* _calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto* _calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_rotate_em_field(
      /* void* */ field.get_fortran_ptr(),
      /* double* */ _w_mat_vec,
      /* double* */ _w_inv_vec,
      /* bool* */ _calc_dfield,
      /* bool* */ _calc_potential);
}
void Bmad::rotate_field_zx(EmFieldProxy& field, double& theta) {
  fortran_rotate_field_zx(
      /* void* */ field.get_fortran_ptr(), /* double& */ theta);
}
void Bmad::rotate_for_curved_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    FixedArray2D<Real, 3, 3> rot_mat) {
  double _rot_mat_vec[3 * 3];
  matrix_to_vec(rot_mat, _rot_mat_vec);
  fortran_rotate_for_curved_surface(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ set,
      /* double* */ _rot_mat_vec);
  vec_to_matrix(_rot_mat_vec, rot_mat);
}
FixedArray1D<Real, 4> Bmad::rotate_spin(
    FixedArray1D<Real, 3> rot_vec,
    FixedArray1D<Real, 3> spin) {
  auto* _rot_vec = rot_vec.data(); // CppWrapperGeneralArgument
  auto* _spin = spin.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _qrot;
  fortran_rotate_spin(
      /* double* */ _rot_vec, /* double* */ _spin, /* double* */ _qrot.data());
  return _qrot;
}
void Bmad::rotate_spin_a_step(
    CoordProxy& orbit,
    EmFieldProxy& field,
    EleProxy& ele,
    double ds) {
  fortran_rotate_spin_a_step(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ field.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ ds);
}
void Bmad::rotate_spin_given_field(
    CoordProxy& orbit,
    int sign_z_vel,
    std::optional<FixedArray1D<Real, 3>> BL,
    std::optional<FixedArray1D<Real, 3>> EL,
    std::optional<FixedArray1D<Real, 4>> qrot) {
  double* _BL = BL.has_value() ? BL.value().data() : nullptr;
  double* _EL = EL.has_value() ? EL.value().data() : nullptr;
  double* _qrot = qrot.has_value() ? qrot.value().data() : nullptr;
  fortran_rotate_spin_given_field(
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ sign_z_vel,
      /* double* */ _BL,
      /* double* */ _EL,
      /* double* */ _qrot);
}
void Bmad::s_body_calc(CoordProxy& orbit, EleProxy& ele, double& s_body) {
  fortran_s_body_calc(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ s_body);
}
void Bmad::s_calc(LatProxy& lat) {
  fortran_s_calc(/* void* */ lat.get_fortran_ptr());
}
void Bmad::sad_mult_hard_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sad_mult_hard_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::sad_soft_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sad_soft_bend_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orb.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::save_a_beam_step(
    EleProxy& ele,
    BeamProxy& beam,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks,
    std::optional<double> s_body,
    std::optional<bool> is_time_coords) {
  // intent=in allocatable type array
  auto* _bunch_tracks = bunch_tracks.has_value()
      ? bunch_tracks->get().get_fortran_ptr()
      : nullptr; // input, optional
  double s_body_lvalue;
  auto* _s_body{&s_body_lvalue};
  if (s_body.has_value()) {
    s_body_lvalue = s_body.value();
  } else {
    _s_body = nullptr;
  }
  bool is_time_coords_lvalue;
  auto* _is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  fortran_save_a_beam_step(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ beam.get_fortran_ptr(),
      /* void* */ _bunch_tracks,
      /* double* */ _s_body,
      /* bool* */ _is_time_coords);
}
void Bmad::save_a_bunch_step(
    EleProxy& ele,
    BunchProxy& bunch,
    optional_ref<BunchTrackProxy> bunch_track,
    std::optional<double> s_body,
    std::optional<bool> is_time_coords) {
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  double s_body_lvalue;
  auto* _s_body{&s_body_lvalue};
  if (s_body.has_value()) {
    s_body_lvalue = s_body.value();
  } else {
    _s_body = nullptr;
  }
  bool is_time_coords_lvalue;
  auto* _is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  fortran_save_a_bunch_step(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ _bunch_track,
      /* double* */ _s_body,
      /* bool* */ _is_time_coords);
}
void Bmad::save_a_step(
    TrackProxy& track,
    EleProxy& ele,
    LatParamProxy& param,
    bool local_ref_frame,
    CoordProxy& orb,
    double s_rel,
    std::optional<bool> save_field,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<double> rf_time,
    optional_ref<StrongBeamProxy> strong_beam) {
  bool save_field_lvalue;
  auto* _save_field{&save_field_lvalue};
  if (save_field.has_value()) {
    save_field_lvalue = save_field.value();
  } else {
    _save_field = nullptr;
  }
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  double rf_time_lvalue;
  auto* _rf_time{&rf_time_lvalue};
  if (rf_time.has_value()) {
    rf_time_lvalue = rf_time.value();
  } else {
    _rf_time = nullptr;
  }
  auto* _strong_beam = strong_beam.has_value()
      ? strong_beam->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_save_a_step(
      /* void* */ track.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ local_ref_frame,
      /* void* */ orb.get_fortran_ptr(),
      /* double& */ s_rel,
      /* bool* */ _save_field,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* double* */ _rf_time,
      /* void* */ _strong_beam);
}
void Bmad::sbend_body_with_k1_map(
    EleProxy& ele,
    double dg,
    double b1,
    LatParamProxy& param,
    int n_step,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sbend_body_with_k1_map(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ dg,
      /* double& */ b1,
      /* void* */ param.get_fortran_ptr(),
      /* int& */ n_step,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
double Bmad::sc_adaptive_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool& include_image,
    double t_now,
    double& dt_step,
    EmFieldProxyAlloc1D& sc_field) {
  double _dt_next{};
  // intent=in allocatable type array
  fortran_sc_adaptive_step(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ include_image,
      /* double& */ t_now,
      /* double& */ dt_step,
      /* double& */ _dt_next,
      /* void* */ sc_field.get_fortran_ptr());
  return _dt_next;
}
int Bmad::sc_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool& include_image,
    double t_end,
    EmFieldProxyAlloc1D& sc_field) {
  // intent=in allocatable type array
  int _n_emit{};
  fortran_sc_step(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ include_image,
      /* double& */ t_end,
      /* void* */ sc_field.get_fortran_ptr(),
      /* int& */ _n_emit);
  return _n_emit;
}
CoordProxy Bmad::set_active_fixer(
    EleProxy& fixer,
    std::optional<bool> turn_on) {
  bool turn_on_lvalue;
  auto* _turn_on{&turn_on_lvalue};
  if (turn_on.has_value()) {
    turn_on_lvalue = turn_on.value();
  } else {
    _turn_on = nullptr;
  }
  CoordProxy _orbit;
  fortran_set_active_fixer(
      /* void* */ fixer.get_fortran_ptr(),
      /* bool* */ _turn_on,
      /* void* */ _orbit.get_fortran_ptr());
  return std::move(_orbit);
}
bool Bmad::set_custom_attribute_name(
    std::string custom_name,
    std::optional<int> custom_index) {
  auto _custom_name = custom_name.c_str();
  bool _err_flag{};
  int custom_index_lvalue;
  auto* _custom_index{&custom_index_lvalue};
  if (custom_index.has_value()) {
    custom_index_lvalue = custom_index.value();
  } else {
    _custom_index = nullptr;
  }
  fortran_set_custom_attribute_name(
      /* const char* */ _custom_name,
      /* bool& */ _err_flag,
      /* int* */ _custom_index);
  return _err_flag;
}
Bmad::SetEleAttribute Bmad::set_ele_attribute(
    EleProxy& ele,
    std::string set_string,
    std::optional<bool> err_print_flag,
    std::optional<bool> set_lords) {
  auto _set_string = set_string.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool set_lords_lvalue;
  auto* _set_lords{&set_lords_lvalue};
  if (set_lords.has_value()) {
    set_lords_lvalue = set_lords.value();
  } else {
    _set_lords = nullptr;
  }
  int _err_id{};
  fortran_set_ele_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _set_string,
      /* bool& */ _err_flag,
      /* bool* */ _err_print_flag,
      /* bool* */ _set_lords,
      /* int& */ _err_id);
  return SetEleAttribute{_err_flag, _err_id};
}
void Bmad::set_ele_defaults(EleProxy& ele, std::optional<bool> do_allocate) {
  bool do_allocate_lvalue;
  auto* _do_allocate{&do_allocate_lvalue};
  if (do_allocate.has_value()) {
    do_allocate_lvalue = do_allocate.value();
  } else {
    _do_allocate = nullptr;
  }
  fortran_set_ele_defaults(
      /* void* */ ele.get_fortran_ptr(), /* bool* */ _do_allocate);
}
void Bmad::set_ele_name(EleProxy& ele, std::string name) {
  auto _name = name.c_str();
  fortran_set_ele_name(
      /* void* */ ele.get_fortran_ptr(), /* const char* */ _name);
}
bool Bmad::set_ele_real_attribute(
    EleProxy& ele,
    std::string attrib_name,
    double value,
    std::optional<bool> err_print_flag) {
  auto _attrib_name = attrib_name.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  fortran_set_ele_real_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_name,
      /* double& */ value,
      /* bool& */ _err_flag,
      /* bool* */ _err_print_flag);
  return _err_flag;
}
Bmad::SetEleStatusStale Bmad::set_ele_status_stale() {
  EleProxy _ele;
  int _status_group{};
  bool _set_slaves{};
  fortran_set_ele_status_stale(
      /* void* */ _ele.get_fortran_ptr(),
      /* int& */ _status_group,
      /* bool& */ _set_slaves);
  return SetEleStatusStale{std::move(_ele), _status_group, _set_slaves};
}
void Bmad::set_flags_for_changed_attribute(
    EleProxy& ele,
    int& attrib,
    std::optional<bool> set_dependent) {
  bool set_dependent_lvalue;
  auto* _set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_integer_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ attrib,
      /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(
    LatProxy& lat,
    std::optional<bool> set_dependent) {
  bool set_dependent_lvalue;
  auto* _set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_lat_attribute(
      /* void* */ lat.get_fortran_ptr(), /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(
    EleProxy& ele,
    bool& attrib,
    std::optional<bool> set_dependent) {
  bool set_dependent_lvalue;
  auto* _set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_logical_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ attrib,
      /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(
    EleProxy& ele,
    optional_ref<double> attrib,
    std::optional<bool> set_dependent) {
  auto* _attrib =
      attrib.has_value() ? &attrib->get() : nullptr; // inout, optional
  bool set_dependent_lvalue;
  auto* _set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_real_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _attrib,
      /* bool* */ _set_dependent);
}
void Bmad::set_fringe_on_off(double& fringe_at, int ele_end, int on_or_off) {
  fortran_set_fringe_on_off(
      /* double& */ fringe_at, /* int& */ ele_end, /* int& */ on_or_off);
}
void Bmad::set_lords_status_stale(
    EleProxy& ele,
    int stat_group,
    std::optional<bool> control_bookkeeping,
    std::optional<int> flag) {
  bool control_bookkeeping_lvalue;
  auto* _control_bookkeeping{&control_bookkeeping_lvalue};
  if (control_bookkeeping.has_value()) {
    control_bookkeeping_lvalue = control_bookkeeping.value();
  } else {
    _control_bookkeeping = nullptr;
  }
  int flag_lvalue;
  auto* _flag{&flag_lvalue};
  if (flag.has_value()) {
    flag_lvalue = flag.value();
  } else {
    _flag = nullptr;
  }
  fortran_set_lords_status_stale(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ stat_group,
      /* bool* */ _control_bookkeeping,
      /* int* */ _flag);
}
void Bmad::set_on_off(
    int key,
    LatProxy& lat,
    int switch_,
    optional_ref<CoordProxyAlloc1D> orb,
    std::optional<bool> use_ref_orb,
    std::optional<int> ix_branch,
    optional_ref<RealAlloc1D> saved_values,
    std::optional<std::string> attribute,
    std::optional<int> set_val) {
  // intent=in allocatable type array
  auto* _orb = orb.has_value() ? orb->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool use_ref_orb_lvalue;
  auto* _use_ref_orb{&use_ref_orb_lvalue};
  if (use_ref_orb.has_value()) {
    use_ref_orb_lvalue = use_ref_orb.value();
  } else {
    _use_ref_orb = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  // intent=inout allocatable general array
  auto* _saved_values = saved_values.has_value()
      ? saved_values->get().get_fortran_ptr()
      : nullptr; // input, optional
  const char* _attribute = attribute.has_value() ? attribute->c_str() : nullptr;
  int set_val_lvalue;
  auto* _set_val{&set_val_lvalue};
  if (set_val.has_value()) {
    set_val_lvalue = set_val.value();
  } else {
    _set_val = nullptr;
  }
  fortran_set_on_off(
      /* int& */ key,
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ switch_,
      /* void* */ _orb,
      /* bool* */ _use_ref_orb,
      /* int* */ _ix_branch,
      /* void* */ _saved_values,
      /* const char* */ _attribute,
      /* int* */ _set_val);
}
CoordProxyAlloc1D Bmad::set_orbit_to_zero(
    int n1,
    int n2,
    std::optional<int> ix_noset) {
  // intent=out allocatable type array
  auto orbit{CoordProxyAlloc1D()};
  int ix_noset_lvalue;
  auto* _ix_noset{&ix_noset_lvalue};
  if (ix_noset.has_value()) {
    ix_noset_lvalue = ix_noset.value();
  } else {
    _ix_noset = nullptr;
  }
  fortran_set_orbit_to_zero(
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ n1,
      /* int& */ n2,
      /* int* */ _ix_noset);
  return std::move(orbit);
}
void Bmad::set_ptc(
    std::optional<double> e_tot,
    std::optional<int> particle,
    std::optional<int> taylor_order,
    std::optional<int> integ_order,
    std::optional<int> n_step,
    std::optional<bool> no_cavity,
    std::optional<bool> force_init) {
  double e_tot_lvalue;
  auto* _e_tot{&e_tot_lvalue};
  if (e_tot.has_value()) {
    e_tot_lvalue = e_tot.value();
  } else {
    _e_tot = nullptr;
  }
  int particle_lvalue;
  auto* _particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int taylor_order_lvalue;
  auto* _taylor_order{&taylor_order_lvalue};
  if (taylor_order.has_value()) {
    taylor_order_lvalue = taylor_order.value();
  } else {
    _taylor_order = nullptr;
  }
  int integ_order_lvalue;
  auto* _integ_order{&integ_order_lvalue};
  if (integ_order.has_value()) {
    integ_order_lvalue = integ_order.value();
  } else {
    _integ_order = nullptr;
  }
  int n_step_lvalue;
  auto* _n_step{&n_step_lvalue};
  if (n_step.has_value()) {
    n_step_lvalue = n_step.value();
  } else {
    _n_step = nullptr;
  }
  bool no_cavity_lvalue;
  auto* _no_cavity{&no_cavity_lvalue};
  if (no_cavity.has_value()) {
    no_cavity_lvalue = no_cavity.value();
  } else {
    _no_cavity = nullptr;
  }
  bool force_init_lvalue;
  auto* _force_init{&force_init_lvalue};
  if (force_init.has_value()) {
    force_init_lvalue = force_init.value();
  } else {
    _force_init = nullptr;
  }
  fortran_set_ptc(
      /* double* */ _e_tot,
      /* int* */ _particle,
      /* int* */ _taylor_order,
      /* int* */ _integ_order,
      /* int* */ _n_step,
      /* bool* */ _no_cavity,
      /* bool* */ _force_init);
}
bool Bmad::set_ptc_base_state(std::string component, bool set_val) {
  auto _component = component.c_str();
  bool _old_val{};
  fortran_set_ptc_base_state(
      /* const char* */ _component, /* bool& */ set_val, /* bool& */ _old_val);
  return _old_val;
}
void Bmad::set_ptc_com_pointers() {
  fortran_set_ptc_com_pointers();
}
void Bmad::set_ptc_quiet(int channel, bool set, int& old_val) {
  fortran_set_ptc_quiet(
      /* int& */ channel, /* bool& */ set, /* int& */ old_val);
}
void Bmad::set_ptc_verbose(bool& on) {
  fortran_set_ptc_verbose(/* bool& */ on);
}
void Bmad::set_pwd_ele(
    LatProxy& lat,
    NormalModesProxy& mode0,
    double inductance) {
  fortran_set_pwd_ele(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ mode0.get_fortran_ptr(),
      /* double& */ inductance);
}
BookkeepingStateProxy Bmad::set_status_flags(int stat) {
  BookkeepingStateProxy _bookkeeping_state;
  fortran_set_status_flags(
      /* void* */ _bookkeeping_state.get_fortran_ptr(), /* int& */ stat);
  return std::move(_bookkeeping_state);
}
void Bmad::set_tune(
    double phi_a_set,
    double phi_b_set,
    RealAlloc1D& dk1,
    ElePointerProxyAlloc1D& eles,
    BranchProxy& branch,
    CoordProxyAlloc1D& orb,
    std::optional<bool> print_err,
    bool& ok) {
  // intent=in allocatable general array
  // intent=in allocatable type array
  // intent=inout allocatable type array
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_tune(
      /* double& */ phi_a_set,
      /* double& */ phi_b_set,
      /* void* */ dk1.get_fortran_ptr(),
      /* void* */ eles.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ orb.get_fortran_ptr(),
      /* bool* */ _print_err,
      /* bool& */ ok);
}
void Bmad::set_twiss(
    BranchProxy& branch,
    EleProxy& twiss_ele,
    int ix_ele,
    bool match_deta_ds,
    bool err_flag,
    std::optional<bool> print_err) {
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_twiss(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ twiss_ele.get_fortran_ptr(),
      /* int& */ ix_ele,
      /* bool& */ match_deta_ds,
      /* bool& */ err_flag,
      /* bool* */ _print_err);
}
bool Bmad::set_z_tune(
    BranchProxy& branch,
    double z_tune,
    std::optional<bool> print_err) {
  bool _ok{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_z_tune(
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ z_tune,
      /* bool& */ _ok,
      /* bool* */ _print_err);
  return _ok;
}
void Bmad::settable_dep_var_bookkeeping(EleProxy& ele) {
  fortran_settable_dep_var_bookkeeping(/* void* */ ele.get_fortran_ptr());
}
void Bmad::setup_high_energy_space_charge_calc(
    bool calc_on,
    BranchProxy& branch,
    double n_part,
    NormalModesProxy& mode,
    optional_ref<CoordProxyAlloc1D> closed_orb) {
  // intent=in allocatable type array
  auto* _closed_orb = closed_orb.has_value()
      ? closed_orb->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_setup_high_energy_space_charge_calc(
      /* bool& */ calc_on,
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ n_part,
      /* void* */ mode.get_fortran_ptr(),
      /* void* */ _closed_orb);
}
FixedArray2D<Real, 6, 6> Bmad::sigma_mat_ptc_to_bmad(
    FixedArray2D<Real, 6, 6> sigma_mat_ptc,
    double beta0) {
  double _sigma_mat_ptc_vec[6 * 6];
  matrix_to_vec(sigma_mat_ptc, _sigma_mat_ptc_vec);
  FixedArray2D<Real, 6, 6> sigma_mat_bmad;
  double _sigma_mat_bmad_vec[6 * 6];
  fortran_sigma_mat_ptc_to_bmad(
      /* double* */ _sigma_mat_ptc_vec,
      /* double& */ beta0,
      /* double* */ _sigma_mat_bmad_vec);
  vec_to_matrix(_sigma_mat_bmad_vec, sigma_mat_bmad);
  return sigma_mat_bmad;
}
void Bmad::significant_difference(
    double value1,
    double value2,
    std::optional<double> abs_tol,
    std::optional<double> rel_tol,
    bool& is_different) {
  double abs_tol_lvalue;
  auto* _abs_tol{&abs_tol_lvalue};
  if (abs_tol.has_value()) {
    abs_tol_lvalue = abs_tol.value();
  } else {
    _abs_tol = nullptr;
  }
  double rel_tol_lvalue;
  auto* _rel_tol{&rel_tol_lvalue};
  if (rel_tol.has_value()) {
    rel_tol_lvalue = rel_tol.value();
  } else {
    _rel_tol = nullptr;
  }
  fortran_significant_difference(
      /* double& */ value1,
      /* double& */ value2,
      /* double* */ _abs_tol,
      /* double* */ _rel_tol,
      /* bool& */ is_different);
}
void Bmad::skip_ele_blender(EleProxy& ele, bool& skip) {
  fortran_skip_ele_blender(/* void* */ ele.get_fortran_ptr(), /* bool& */ skip);
}
bool Bmad::slice_lattice(
    LatProxy& lat,
    std::string ele_list,
    std::optional<bool> do_bookkeeping) {
  auto _ele_list = ele_list.c_str();
  bool _error{};
  bool do_bookkeeping_lvalue;
  auto* _do_bookkeeping{&do_bookkeeping_lvalue};
  if (do_bookkeeping.has_value()) {
    do_bookkeeping_lvalue = do_bookkeeping.value();
  } else {
    _do_bookkeeping = nullptr;
  }
  fortran_slice_lattice(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _ele_list,
      /* bool& */ _error,
      /* bool* */ _do_bookkeeping);
  return _error;
}
void Bmad::soft_quadrupole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_soft_quadrupole_edge_kick(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* int& */ particle_at,
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::sol_quad_mat6_calc(
    double& ks_in,
    double& k1_in,
    double tilt,
    double length,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sol_quad_mat6_calc(
      /* double& */ ks_in,
      /* double& */ k1_in,
      /* double& */ tilt,
      /* double& */ length,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
double Bmad::solve_psi_adaptive(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  double _p1{};
  fortran_solve_psi_adaptive(
      /* double& */ t0,
      /* double& */ t1,
      /* double& */ p0,
      /* double* */ _args,
      /* double& */ _p1);
  return _p1;
}
Bmad::SolvePsiFixedSteps Bmad::solve_psi_fixed_steps(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args) {
  auto* _args = args.data(); // CppWrapperGeneralArgument
  // intent=out allocatable general array
  auto t{RealAlloc1D()};
  // intent=out allocatable general array
  auto p{RealAlloc1D()};
  fortran_solve_psi_fixed_steps(
      /* double& */ t0,
      /* double& */ t1,
      /* double& */ p0,
      /* double* */ _args,
      /* void* */ t.get_fortran_ptr(),
      /* void* */ p.get_fortran_ptr());
  return SolvePsiFixedSteps{std::move(t), std::move(p)};
}
ComplexTaylorProxy Bmad::sort_complex_taylor_terms(
    ComplexTaylorProxy& complex_taylor_in) {
  ComplexTaylorProxy _complex_taylor_sorted;
  fortran_sort_complex_taylor_terms(
      /* void* */ complex_taylor_in.get_fortran_ptr(),
      /* void* */ _complex_taylor_sorted.get_fortran_ptr());
  return std::move(_complex_taylor_sorted);
}
bool Bmad::spin_dn_dpz_from_mat8(
    FixedArray2D<Real, 8, 8> mat_1turn,
    std::optional<FixedArray2D<Real, 3, 3>> dn_dpz_partial,
    FixedArray1D<Real, 3> dn_dpz) {
  double _mat_1turn_vec[8 * 8];
  matrix_to_vec(mat_1turn, _mat_1turn_vec);
  double _dn_dpz_partial_vec[3 * 3];
  const double* _dn_dpz_partial = nullptr;
  if (dn_dpz_partial.has_value()) {
    matrix_to_vec(dn_dpz_partial.value(), _dn_dpz_partial_vec);
    _dn_dpz_partial = _dn_dpz_partial_vec;
  }
  bool _error{};
  auto* _dn_dpz = dn_dpz.data(); // CppWrapperGeneralArgument
  fortran_spin_dn_dpz_from_mat8(
      /* double* */ _mat_1turn_vec,
      /* double* */ _dn_dpz_partial_vec,
      /* bool& */ _error,
      /* double* */ _dn_dpz);
  return _error;
}
bool Bmad::spin_dn_dpz_from_qmap(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> q_map,
    FixedArray2D<Real, 3, 3> dn_dpz_partial,
    FixedArray2D<Real, 3, 3> dn_dpz_partial2,
    std::optional<FixedArray1D<Real, 3>> n0,
    FixedArray1D<Real, 3> dn_dpz) {
  double _orb_mat_vec[6 * 6];
  matrix_to_vec(orb_mat, _orb_mat_vec);
  double _q_map_vec[4 * 7];
  matrix_to_vec(q_map, _q_map_vec);
  double _dn_dpz_partial_vec[3 * 3];
  matrix_to_vec(dn_dpz_partial, _dn_dpz_partial_vec);
  double _dn_dpz_partial2_vec[3 * 3];
  matrix_to_vec(dn_dpz_partial2, _dn_dpz_partial2_vec);
  bool _error{};
  double* _n0 = n0.has_value() ? n0.value().data() : nullptr;
  auto* _dn_dpz = dn_dpz.data(); // CppWrapperGeneralArgument
  fortran_spin_dn_dpz_from_qmap(
      /* double* */ _orb_mat_vec,
      /* double* */ _q_map_vec,
      /* double* */ _dn_dpz_partial_vec,
      /* double* */ _dn_dpz_partial2_vec,
      /* bool& */ _error,
      /* double* */ _n0,
      /* double* */ _dn_dpz);
  return _error;
}
void Bmad::spin_map1_normalize(FixedArray2D<Real, 4, 7> spin1) {
  double _spin1_vec[4 * 7];
  matrix_to_vec(spin1, _spin1_vec);
  fortran_spin_map1_normalize(/* double* */ _spin1_vec);
  vec_to_matrix(_spin1_vec, spin1);
}
Bmad::SpinMat8ResonanceStrengths Bmad::spin_mat8_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 6, 6> mat8) {
  auto* _orb_evec = orb_evec.data(); // CppWrapperGeneralArgument
  double _mat8_vec[6 * 6];
  matrix_to_vec(mat8, _mat8_vec);
  double _xi_sum{};
  double _xi_diff{};
  fortran_spin_mat8_resonance_strengths(
      /* std::complex<double>* */ _orb_evec,
      /* double* */ _mat8_vec,
      /* double& */ _xi_sum,
      /* double& */ _xi_diff);
  return SpinMat8ResonanceStrengths{_xi_sum, _xi_diff};
}
Bmad::SpinMatToEigen Bmad::spin_mat_to_eigen(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> spin_map) {
  double _orb_mat_vec[6 * 6];
  matrix_to_vec(orb_mat, _orb_mat_vec);
  double _spin_map_vec[4 * 7];
  matrix_to_vec(spin_map, _spin_map_vec);
  FixedArray1D<Complex, 6> _orb_eval;
  FixedArray2D<Complex, 6, 6> orb_evec;
  std::complex<double> _orb_evec_vec[6 * 6];
  FixedArray1D<Real, 3> _n0;
  FixedArray2D<Complex, 6, 3> spin_evec;
  std::complex<double> _spin_evec_vec[6 * 3];
  bool _error{};
  fortran_spin_mat_to_eigen(
      /* double* */ _orb_mat_vec,
      /* double* */ _spin_map_vec,
      /* std::complex<double>* */ _orb_eval.data(),
      /* std::complex<double>* */ _orb_evec_vec,
      /* double* */ _n0.data(),
      /* std::complex<double>* */ _spin_evec_vec,
      /* bool& */ _error);
  vec_to_matrix(_orb_evec_vec, orb_evec);
  vec_to_matrix(_spin_evec_vec, spin_evec);
  return SpinMatToEigen{_orb_eval, orb_evec, _n0, spin_evec, _error};
}
void Bmad::spin_omega(
    EmFieldProxy& field,
    CoordProxy& orbit,
    int& sign_z_vel,
    optional_ref<bool> phase_space_coords,
    FixedArray1D<Real, 3> omega) {
  auto* _phase_space_coords = phase_space_coords.has_value()
      ? &phase_space_coords->get()
      : nullptr; // inout, optional
  auto* _omega = omega.data(); // CppWrapperGeneralArgument
  fortran_spin_omega(
      /* void* */ field.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ sign_z_vel,
      /* bool* */ _phase_space_coords,
      /* double* */ _omega);
}
Bmad::SpinQuatResonanceStrengths Bmad::spin_quat_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 4, 7> spin_q) {
  auto* _orb_evec = orb_evec.data(); // CppWrapperGeneralArgument
  double _spin_q_vec[4 * 7];
  matrix_to_vec(spin_q, _spin_q_vec);
  double _xi_sum{};
  double _xi_diff{};
  fortran_spin_quat_resonance_strengths(
      /* std::complex<double>* */ _orb_evec,
      /* double* */ _spin_q_vec,
      /* double& */ _xi_sum,
      /* double& */ _xi_diff);
  return SpinQuatResonanceStrengths{_xi_sum, _xi_diff};
}
void Bmad::spin_taylor_to_linear(
    FixedArray1D<TaylorProxy, 4> spin_taylor,
    bool normalize,
    FixedArray1D<Real, 6> dref_orb,
    bool is_on,
    FixedArray2D<Real, 4, 7> spin_map1) {
  auto* _dref_orb = dref_orb.data(); // CppWrapperGeneralArgument
  double _spin_map1_vec[4 * 7];
  matrix_to_vec(spin_map1, _spin_map1_vec);
  fortran_spin_taylor_to_linear(
      /* void* */ spin_taylor.data(),
      /* bool& */ normalize,
      /* double* */ _dref_orb,
      /* bool& */ is_on,
      /* double* */ _spin_map1_vec);
  vec_to_matrix(_spin_map1_vec, spin_map1);
}
void Bmad::spinor_to_polar(
    FixedArray1D<Complex, 2> spinor,
    SpinPolarProxy& polar) {
  auto* _spinor = spinor.data(); // CppWrapperGeneralArgument
  fortran_spinor_to_polar(
      /* std::complex<double>* */ _spinor, /* void* */ polar.get_fortran_ptr());
}
void Bmad::spinor_to_vec(
    FixedArray1D<Complex, 2> spinor,
    FixedArray1D<Real, 3> vec) {
  auto* _spinor = spinor.data(); // CppWrapperGeneralArgument
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  fortran_spinor_to_vec(
      /* std::complex<double>* */ _spinor, /* double* */ _vec);
}
void Bmad::spline_fit_orbit(
    CoordProxy& start_orb,
    CoordProxy& end_orb,
    FixedArray1D<Real, 4> spline_x,
    FixedArray1D<Real, 4> spline_y) {
  auto* _spline_x = spline_x.data(); // CppWrapperGeneralArgument
  auto* _spline_y = spline_y.data(); // CppWrapperGeneralArgument
  fortran_spline_fit_orbit(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ end_orb.get_fortran_ptr(),
      /* double* */ _spline_x,
      /* double* */ _spline_y);
}
Bmad::SplitLat Bmad::split_lat(
    LatProxy& lat,
    double s_split,
    int ix_branch,
    std::optional<bool> add_suffix,
    std::optional<bool> check_sanity,
    std::optional<bool> save_null_drift,
    std::optional<bool> choose_max,
    std::optional<int> ix_insert) {
  int _ix_split{};
  bool _split_done{};
  bool add_suffix_lvalue;
  auto* _add_suffix{&add_suffix_lvalue};
  if (add_suffix.has_value()) {
    add_suffix_lvalue = add_suffix.value();
  } else {
    _add_suffix = nullptr;
  }
  bool check_sanity_lvalue;
  auto* _check_sanity{&check_sanity_lvalue};
  if (check_sanity.has_value()) {
    check_sanity_lvalue = check_sanity.value();
  } else {
    _check_sanity = nullptr;
  }
  bool save_null_drift_lvalue;
  auto* _save_null_drift{&save_null_drift_lvalue};
  if (save_null_drift.has_value()) {
    save_null_drift_lvalue = save_null_drift.value();
  } else {
    _save_null_drift = nullptr;
  }
  bool _err_flag{};
  bool choose_max_lvalue;
  auto* _choose_max{&choose_max_lvalue};
  if (choose_max.has_value()) {
    choose_max_lvalue = choose_max.value();
  } else {
    _choose_max = nullptr;
  }
  int ix_insert_lvalue;
  auto* _ix_insert{&ix_insert_lvalue};
  if (ix_insert.has_value()) {
    ix_insert_lvalue = ix_insert.value();
  } else {
    _ix_insert = nullptr;
  }
  fortran_split_lat(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ s_split,
      /* int& */ ix_branch,
      /* int& */ _ix_split,
      /* bool& */ _split_done,
      /* bool* */ _add_suffix,
      /* bool* */ _check_sanity,
      /* bool* */ _save_null_drift,
      /* bool& */ _err_flag,
      /* bool* */ _choose_max,
      /* int* */ _ix_insert);
  return SplitLat{_ix_split, _split_done, _err_flag};
}
void Bmad::sprint_spin_taylor_map(
    EleProxy& ele,
    std::optional<FixedArray1D<Real, 6>> start_orbit) {
  double* _start_orbit =
      start_orbit.has_value() ? start_orbit.value().data() : nullptr;
  fortran_sprint_spin_taylor_map(
      /* void* */ ele.get_fortran_ptr(), /* double* */ _start_orbit);
}
void Bmad::sr_longitudinal_wake_particle(EleProxy& ele, CoordProxy& orbit) {
  fortran_sr_longitudinal_wake_particle(
      /* void* */ ele.get_fortran_ptr(), /* void* */ orbit.get_fortran_ptr());
}
void Bmad::sr_transverse_wake_particle(EleProxy& ele, CoordProxy& orbit) {
  fortran_sr_transverse_wake_particle(
      /* void* */ ele.get_fortran_ptr(), /* void* */ orbit.get_fortran_ptr());
}
void Bmad::sr_z_long_wake(EleProxy& ele, BunchProxy& bunch, double z_ave) {
  fortran_sr_z_long_wake(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ bunch.get_fortran_ptr(),
      /* double& */ z_ave);
}
SummationRdtProxy Bmad::srdt_calc(
    LatProxy& lat,
    int order,
    std::optional<int> n_slices_gen_opt,
    std::optional<int> n_slices_sxt_opt,
    optional_ref<SummationRdtProxyAlloc1D> per_ele_out) {
  SummationRdtProxy _srdt_sums;
  int n_slices_gen_opt_lvalue;
  auto* _n_slices_gen_opt{&n_slices_gen_opt_lvalue};
  if (n_slices_gen_opt.has_value()) {
    n_slices_gen_opt_lvalue = n_slices_gen_opt.value();
  } else {
    _n_slices_gen_opt = nullptr;
  }
  int n_slices_sxt_opt_lvalue;
  auto* _n_slices_sxt_opt{&n_slices_sxt_opt_lvalue};
  if (n_slices_sxt_opt.has_value()) {
    n_slices_sxt_opt_lvalue = n_slices_sxt_opt.value();
  } else {
    _n_slices_sxt_opt = nullptr;
  }
  // intent=inout allocatable type array
  auto* _per_ele_out = per_ele_out.has_value()
      ? per_ele_out->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_srdt_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _srdt_sums.get_fortran_ptr(),
      /* int& */ order,
      /* int* */ _n_slices_gen_opt,
      /* int* */ _n_slices_sxt_opt,
      /* void* */ _per_ele_out);
  return std::move(_srdt_sums);
}
RealAlloc1D Bmad::srdt_lsq_solution(
    LatProxy& lat,
    IntAlloc1D& var_indexes,
    std::optional<int> n_slices_gen_opt,
    std::optional<int> n_slices_sxt_opt,
    std::optional<double> chrom_set_x_opt,
    std::optional<double> chrom_set_y_opt,
    std::optional<FixedArray1D<Real, 10>> weight_in) {
  // intent=in allocatable general array
  // intent=out allocatable general array
  auto ls_soln{RealAlloc1D()};
  int n_slices_gen_opt_lvalue;
  auto* _n_slices_gen_opt{&n_slices_gen_opt_lvalue};
  if (n_slices_gen_opt.has_value()) {
    n_slices_gen_opt_lvalue = n_slices_gen_opt.value();
  } else {
    _n_slices_gen_opt = nullptr;
  }
  int n_slices_sxt_opt_lvalue;
  auto* _n_slices_sxt_opt{&n_slices_sxt_opt_lvalue};
  if (n_slices_sxt_opt.has_value()) {
    n_slices_sxt_opt_lvalue = n_slices_sxt_opt.value();
  } else {
    _n_slices_sxt_opt = nullptr;
  }
  double chrom_set_x_opt_lvalue;
  auto* _chrom_set_x_opt{&chrom_set_x_opt_lvalue};
  if (chrom_set_x_opt.has_value()) {
    chrom_set_x_opt_lvalue = chrom_set_x_opt.value();
  } else {
    _chrom_set_x_opt = nullptr;
  }
  double chrom_set_y_opt_lvalue;
  auto* _chrom_set_y_opt{&chrom_set_y_opt_lvalue};
  if (chrom_set_y_opt.has_value()) {
    chrom_set_y_opt_lvalue = chrom_set_y_opt.value();
  } else {
    _chrom_set_y_opt = nullptr;
  }
  double* _weight_in =
      weight_in.has_value() ? weight_in.value().data() : nullptr;
  fortran_srdt_lsq_solution(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ var_indexes.get_fortran_ptr(),
      /* void* */ ls_soln.get_fortran_ptr(),
      /* int* */ _n_slices_gen_opt,
      /* int* */ _n_slices_sxt_opt,
      /* double* */ _chrom_set_x_opt,
      /* double* */ _chrom_set_y_opt,
      /* double* */ _weight_in);
  return std::move(ls_soln);
}
bool Bmad::start_branch_at(
    LatProxy& lat,
    std::string ele_start,
    bool move_end_marker) {
  auto _ele_start = ele_start.c_str();
  bool _error{};
  fortran_start_branch_at(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _ele_start,
      /* bool& */ move_end_marker,
      /* bool& */ _error);
  return _error;
}
void Bmad::stream_ele_end(
    int physical_end,
    int ele_orientation,
    int& stream_end) {
  fortran_stream_ele_end(
      /* int& */ physical_end,
      /* int& */ ele_orientation,
      /* int& */ stream_end);
}
std::string Bmad::string_attrib(std::string attrib_name, EleProxy& ele) {
  auto _attrib_name = attrib_name.c_str();
  char _attrib_value[4096];
  fortran_string_attrib(
      /* const char* */ _attrib_name,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_value);
  return _attrib_value;
}
Bmad::StrongBeamSigmaCalc Bmad::strong_beam_sigma_calc(
    EleProxy& ele,
    double s_pos) {
  FixedArray1D<Real, 2> _sigma;
  double _bbi_const{};
  FixedArray1D<Real, 2> _dsigma_ds;
  fortran_strong_beam_sigma_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ s_pos,
      /* double* */ _sigma.data(),
      /* double& */ _bbi_const,
      /* double* */ _dsigma_ds.data());
  return StrongBeamSigmaCalc{_sigma, _bbi_const, _dsigma_ds};
}
void Bmad::strong_beam_strength(EleProxy& ele, double& strength) {
  fortran_strong_beam_strength(
      /* void* */ ele.get_fortran_ptr(), /* double& */ strength);
}
void Bmad::surface_grid_displacement(
    EleProxy& ele,
    double& x,
    double& y,
    bool err_flag,
    double z,
    std::optional<FixedArray1D<Real, 2>> dz_dxy,
    std::optional<bool> extend_grid) {
  double* _dz_dxy = dz_dxy.has_value() ? dz_dxy.value().data() : nullptr;
  bool extend_grid_lvalue;
  auto* _extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  fortran_surface_grid_displacement(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ x,
      /* double& */ y,
      /* bool& */ err_flag,
      /* double& */ z,
      /* double* */ _dz_dxy,
      /* bool* */ _extend_grid);
}
TrackProxy Bmad::symp_lie_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<bool> offset_ele) {
  TrackProxy _track;
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  bool offset_ele_lvalue;
  auto* _offset_ele{&offset_ele_lvalue};
  if (offset_ele.has_value()) {
    offset_ele_lvalue = offset_ele.value();
  } else {
    _offset_ele = nullptr;
  }
  fortran_symp_lie_bmad(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* bool* */ _offset_ele);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return std::move(_track);
}
Bmad::T6ToB123 Bmad::t6_to_b123(
    FixedArray2D<Real, 6, 6> t6,
    FixedArray1D<Real, 3> abz_tunes) {
  double _t6_vec[6 * 6];
  matrix_to_vec(t6, _t6_vec);
  auto* _abz_tunes = abz_tunes.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 6, 6> B1;
  double _B1_vec[6 * 6];
  FixedArray2D<Real, 6, 6> B2;
  double _B2_vec[6 * 6];
  FixedArray2D<Real, 6, 6> B3;
  double _B3_vec[6 * 6];
  bool _err_flag{};
  fortran_t6_to_b123(
      /* double* */ _t6_vec,
      /* double* */ _abz_tunes,
      /* double* */ _B1_vec,
      /* double* */ _B2_vec,
      /* double* */ _B3_vec,
      /* bool& */ _err_flag);
  vec_to_matrix(_B1_vec, B1);
  vec_to_matrix(_B2_vec, B2);
  vec_to_matrix(_B3_vec, B3);
  return T6ToB123{B1, B2, B3, _err_flag};
}
void Bmad::taper_mag_strengths(
    LatProxy& lat,
    optional_ref<LatProxy> ref_lat,
    std::optional<std::string> except,
    optional_ref<bool> err_flag) {
  auto* _ref_lat = ref_lat.has_value() ? ref_lat->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  const char* _except = except.has_value() ? except->c_str() : nullptr;
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  fortran_taper_mag_strengths(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _ref_lat,
      /* const char* */ _except,
      /* bool* */ _err_flag);
}
void Bmad::target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double& y_min,
    double& y_max,
    double& phi_min,
    double& phi_max,
    std::optional<bool> initial) {
  auto* _r_corner1 = r_corner1.data(); // CppWrapperGeneralArgument
  auto* _r_corner2 = r_corner2.data(); // CppWrapperGeneralArgument
  bool initial_lvalue;
  auto* _initial{&initial_lvalue};
  if (initial.has_value()) {
    initial_lvalue = initial.value();
  } else {
    _initial = nullptr;
  }
  fortran_target_min_max_calc(
      /* double* */ _r_corner1,
      /* double* */ _r_corner2,
      /* double& */ y_min,
      /* double& */ y_max,
      /* double& */ phi_min,
      /* double& */ phi_max,
      /* bool* */ _initial);
}
Bmad::TargetRotMats Bmad::target_rot_mats(FixedArray1D<Real, 3> r_center) {
  auto* _r_center = r_center.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 3, 3> w_to_target;
  double _w_to_target_vec[3 * 3];
  FixedArray2D<Real, 3, 3> w_to_ele;
  double _w_to_ele_vec[3 * 3];
  fortran_target_rot_mats(
      /* double* */ _r_center,
      /* double* */ _w_to_target_vec,
      /* double* */ _w_to_ele_vec);
  vec_to_matrix(_w_to_target_vec, w_to_target);
  vec_to_matrix(_w_to_ele_vec, w_to_ele);
  return TargetRotMats{w_to_target, w_to_ele};
}
void Bmad::taylor_equal_taylor(TaylorProxy& taylor1, TaylorProxy& taylor2) {
  fortran_taylor_equal_taylor(
      /* void* */ taylor1.get_fortran_ptr(),
      /* void* */ taylor2.get_fortran_ptr());
}
Bmad::TaylorInverse Bmad::taylor_inverse(TaylorProxyAlloc1D& taylor_in) {
  // intent=in allocatable type array
  // intent=out allocatable type array
  auto taylor_inv{TaylorProxyAlloc1D()};
  bool _err{};
  fortran_taylor_inverse(
      /* void* */ taylor_in.get_fortran_ptr(),
      /* void* */ taylor_inv.get_fortran_ptr(),
      /* bool& */ _err);
  return TaylorInverse{std::move(taylor_inv), _err};
}
bool Bmad::taylor_propagate1(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> ref_in,
    optional_ref<TaylorProxyAlloc1D> spin_taylor) {
  // intent=inout allocatable type array
  bool _err_flag{};
  auto* _ref_in = ref_in.has_value() ? ref_in->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  // intent=inout allocatable type array
  auto* _spin_taylor = spin_taylor.has_value()
      ? spin_taylor->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_taylor_propagate1(
      /* void* */ orb_taylor.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _ref_in,
      /* void* */ _spin_taylor);
  return _err_flag;
}
MadMapProxy Bmad::taylor_to_mad_map(
    TaylorProxyAlloc1D& taylor,
    MadEnergyProxy& energy) {
  // intent=in allocatable type array
  MadMapProxy _map;
  fortran_taylor_to_mad_map(
      /* void* */ taylor.get_fortran_ptr(),
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
void Bmad::taylors_equal_taylors(
    TaylorProxyAlloc1D& taylor1,
    TaylorProxyAlloc1D& taylor2) {
  // intent=inout allocatable type array
  // intent=in allocatable type array
  fortran_taylors_equal_taylors(
      /* void* */ taylor1.get_fortran_ptr(),
      /* void* */ taylor2.get_fortran_ptr());
}
void Bmad::tilt_coords(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  // intent=inout allocatable general array
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_tilt_coords(
      /* double& */ tilt_val,
      /* void* */ coord.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::tilt_coords_photon(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat) {
  // intent=inout allocatable general array
  double _w_mat_vec[3 * 3];
  const double* _w_mat = nullptr;
  if (w_mat.has_value()) {
    matrix_to_vec(w_mat.value(), _w_mat_vec);
    _w_mat = _w_mat_vec;
  }
  fortran_tilt_coords_photon(
      /* double& */ tilt_val,
      /* void* */ coord.get_fortran_ptr(),
      /* double* */ _w_mat_vec);
  if (w_mat.has_value())
    vec_to_matrix(_w_mat_vec, w_mat.value());
}
void Bmad::tilt_mat6(FixedArray2D<Real, 6, 6> mat6, double tilt) {
  double _mat6_vec[6 * 6];
  matrix_to_vec(mat6, _mat6_vec);
  fortran_tilt_mat6(/* double* */ _mat6_vec, /* double& */ tilt);
  vec_to_matrix(_mat6_vec, mat6);
}
Bmad::ToEtaReading Bmad::to_eta_reading(
    RealAlloc1D& eta_actual,
    EleProxy& ele,
    int axis,
    bool add_noise) {
  // intent=in allocatable general array
  double _reading{};
  bool _err{};
  fortran_to_eta_reading(
      /* void* */ eta_actual.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ axis,
      /* bool& */ add_noise,
      /* double& */ _reading,
      /* bool& */ _err);
  return ToEtaReading{_reading, _err};
}
void Bmad::to_fieldmap_coords(
    EleProxy& ele,
    CoordProxy& local_orb,
    double s_body,
    int ele_anchor_pt,
    FixedArray1D<Real, 3> r0,
    bool curved_ref_frame,
    double& x,
    double& y,
    double& z,
    double& cos_ang,
    double& sin_ang,
    bool err_flag) {
  auto* _r0 = r0.data(); // CppWrapperGeneralArgument
  fortran_to_fieldmap_coords(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ local_orb.get_fortran_ptr(),
      /* double& */ s_body,
      /* int& */ ele_anchor_pt,
      /* double* */ _r0,
      /* bool& */ curved_ref_frame,
      /* double& */ x,
      /* double& */ y,
      /* double& */ z,
      /* double& */ cos_ang,
      /* double& */ sin_ang,
      /* bool& */ err_flag);
}
Bmad::ToOrbitReading Bmad::to_orbit_reading(
    CoordProxy& orb,
    EleProxy& ele,
    int axis,
    bool add_noise) {
  double _reading{};
  bool _err{};
  fortran_to_orbit_reading(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ axis,
      /* bool& */ add_noise,
      /* double& */ _reading,
      /* bool& */ _err);
  return ToOrbitReading{_reading, _err};
}
Bmad::ToPhaseAndCouplingReading Bmad::to_phase_and_coupling_reading(
    EleProxy& ele,
    bool add_noise) {
  BpmPhaseCouplingProxy _reading;
  bool _err{};
  fortran_to_phase_and_coupling_reading(
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ add_noise,
      /* void* */ _reading.get_fortran_ptr(),
      /* bool& */ _err);
  return ToPhaseAndCouplingReading{std::move(_reading), _err};
}
CoordProxy Bmad::to_photon_angle_coords(CoordProxy& orb_in, EleProxy& ele) {
  CoordProxy _orb_out;
  fortran_to_photon_angle_coords(
      /* void* */ orb_in.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _orb_out.get_fortran_ptr());
  return std::move(_orb_out);
}
CoordProxy Bmad::to_surface_coords(CoordProxy& lab_orbit, EleProxy& ele) {
  CoordProxy _surface_orbit;
  fortran_to_surface_coords(
      /* void* */ lab_orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _surface_orbit.get_fortran_ptr());
  return std::move(_surface_orbit);
}
double Bmad::touschek_lifetime(NormalModesProxy& mode, LatProxy& lat) {
  double _Tl{};
  fortran_touschek_lifetime(
      /* void* */ mode.get_fortran_ptr(),
      /* double& */ _Tl,
      /* void* */ lat.get_fortran_ptr());
  return _Tl;
}
double Bmad::touschek_rate1(
    NormalModesProxy& mode,
    LatProxy& lat,
    std::optional<int> ix,
    std::optional<double> s) {
  double _rate{};
  int ix_lvalue;
  auto* _ix{&ix_lvalue};
  if (ix.has_value()) {
    ix_lvalue = ix.value();
  } else {
    _ix = nullptr;
  }
  double s_lvalue;
  auto* _s{&s_lvalue};
  if (s.has_value()) {
    s_lvalue = s.value();
  } else {
    _s = nullptr;
  }
  fortran_touschek_rate1(
      /* void* */ mode.get_fortran_ptr(),
      /* double& */ _rate,
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _ix,
      /* double* */ _s);
  return _rate;
}
void Bmad::touschek_rate1_zap(
    NormalModesProxy& mode,
    double& rate,
    LatProxy& lat,
    optional_ref<int> ix,
    optional_ref<double> s) {
  auto* _ix = ix.has_value() ? &ix->get() : nullptr; // inout, optional
  auto* _s = s.has_value() ? &s->get() : nullptr; // inout, optional
  fortran_touschek_rate1_zap(
      /* void* */ mode.get_fortran_ptr(),
      /* double& */ rate,
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _ix,
      /* double* */ _s);
}
Bmad::Track1 Bmad::track1(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<TrackProxy> track,
    std::optional<bool> ignore_radiation,
    std::optional<bool> make_map1,
    std::optional<bool> init_to_edge) {
  CoordProxy _end_orb;
  auto* _track = track.has_value() ? track->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  bool _err_flag{};
  bool ignore_radiation_lvalue;
  auto* _ignore_radiation{&ignore_radiation_lvalue};
  if (ignore_radiation.has_value()) {
    ignore_radiation_lvalue = ignore_radiation.value();
  } else {
    _ignore_radiation = nullptr;
  }
  bool make_map1_lvalue;
  auto* _make_map1{&make_map1_lvalue};
  if (make_map1.has_value()) {
    make_map1_lvalue = make_map1.value();
  } else {
    _make_map1 = nullptr;
  }
  bool init_to_edge_lvalue;
  auto* _init_to_edge{&init_to_edge_lvalue};
  if (init_to_edge.has_value()) {
    init_to_edge_lvalue = init_to_edge.value();
  } else {
    _init_to_edge = nullptr;
  }
  fortran_track1(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* void* */ _track,
      /* bool& */ _err_flag,
      /* bool* */ _ignore_radiation,
      /* bool* */ _make_map1,
      /* bool* */ _init_to_edge);
  return Track1{std::move(_end_orb), _err_flag};
}
bool Bmad::track1_beam(
    BeamProxy& beam,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid,
    std::optional<int> direction) {
  bool _err{};
  // intent=in allocatable type array
  auto* _centroid = centroid.has_value() ? centroid->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  fortran_track1_beam(
      /* void* */ beam.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err,
      /* void* */ _centroid,
      /* int* */ _direction);
  return _err;
}
Bmad::Track1Bmad Bmad::track1_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  bool _err_flag{};
  TrackProxy _track;
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_bmad(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return Track1Bmad{_err_flag, std::move(_track)};
}
bool Bmad::track1_bmad_photon(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param) {
  bool _err_flag{};
  fortran_track1_bmad_photon(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::track1_bunch(
    BunchProxy& bunch,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid,
    std::optional<int> direction,
    optional_ref<BunchTrackProxy> bunch_track) {
  bool _err{};
  // intent=in allocatable type array
  auto* _centroid = centroid.has_value() ? centroid->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track1_bunch(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err,
      /* void* */ _centroid,
      /* int* */ _direction,
      /* void* */ _bunch_track);
  return _err;
}
bool Bmad::track1_bunch_csr(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start,
    std::optional<double> s_end,
    optional_ref<BunchTrackProxy> bunch_track) {
  // intent=in allocatable type array
  bool _err{};
  double s_start_lvalue;
  auto* _s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto* _s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track1_bunch_csr(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ centroid.get_fortran_ptr(),
      /* bool& */ _err,
      /* double* */ _s_start,
      /* double* */ _s_end,
      /* void* */ _bunch_track);
  return _err;
}
bool Bmad::track1_bunch_csr3d(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start,
    std::optional<double> s_end,
    optional_ref<BunchTrackProxy> bunch_track) {
  // intent=in allocatable type array
  bool _err{};
  double s_start_lvalue;
  auto* _s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto* _s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track1_bunch_csr3d(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ centroid.get_fortran_ptr(),
      /* bool& */ _err,
      /* double* */ _s_start,
      /* double* */ _s_end,
      /* void* */ _bunch_track);
  return _err;
}
void Bmad::track1_bunch_hom(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<int> direction,
    optional_ref<BunchTrackProxy> bunch_track) {
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track1_bunch_hom(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _direction,
      /* void* */ _bunch_track);
}
bool Bmad::track1_bunch_space_charge(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<bool> track_to_same_s,
    optional_ref<BunchTrackProxy> bunch_track) {
  bool _err{};
  bool track_to_same_s_lvalue;
  auto* _track_to_same_s{&track_to_same_s_lvalue};
  if (track_to_same_s.has_value()) {
    track_to_same_s_lvalue = track_to_same_s.value();
  } else {
    _track_to_same_s = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track1_bunch_space_charge(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err,
      /* bool* */ _track_to_same_s,
      /* void* */ _bunch_track);
  return _err;
}
void Bmad::track1_crystal(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_crystal(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_diffraction_plate_or_mask(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_diffraction_plate_or_mask(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_high_energy_space_charge(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_high_energy_space_charge(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_lens(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit) {
  fortran_track1_lens(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_linear(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param) {
  fortran_track1_linear(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr());
}
void Bmad::track1_lr_wake(BunchProxy& bunch, EleProxy& ele) {
  fortran_track1_lr_wake(
      /* void* */ bunch.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
void Bmad::track1_mad(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param) {
  fortran_track1_mad(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr());
}
void Bmad::track1_mirror(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_mirror(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_mosaic_crystal(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_mosaic_crystal(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_multilayer_mirror(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_multilayer_mirror(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_radiation(CoordProxy& orbit, EleProxy& ele, int edge) {
  fortran_track1_radiation(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ edge);
}
void Bmad::track1_radiation_center(
    CoordProxy& orbit,
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> rad_damp,
    std::optional<bool> rad_fluct) {
  bool rad_damp_lvalue;
  auto* _rad_damp{&rad_damp_lvalue};
  if (rad_damp.has_value()) {
    rad_damp_lvalue = rad_damp.value();
  } else {
    _rad_damp = nullptr;
  }
  bool rad_fluct_lvalue;
  auto* _rad_fluct{&rad_fluct_lvalue};
  if (rad_fluct.has_value()) {
    rad_fluct_lvalue = rad_fluct.value();
  } else {
    _rad_fluct = nullptr;
  }
  fortran_track1_radiation_center(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool* */ _rad_damp,
      /* bool* */ _rad_fluct);
}
Bmad::Track1RungeKutta Bmad::track1_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  bool _err_flag{};
  TrackProxy _track;
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_runge_kutta(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return Track1RungeKutta{_err_flag, std::move(_track)};
}
void Bmad::track1_sample(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit) {
  fortran_track1_sample(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
Bmad::Track1Spin Bmad::track1_spin(
    CoordProxy& start_orb,
    LatParamProxy& param,
    optional_ref<bool> make_quaternion) {
  EleProxy _ele;
  CoordProxy _end_orb;
  auto* _make_quaternion = make_quaternion.has_value()
      ? &make_quaternion->get()
      : nullptr; // inout, optional
  fortran_track1_spin(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ _ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* bool* */ _make_quaternion);
  return Track1Spin{std::move(_ele), std::move(_end_orb)};
}
CoordProxy Bmad::track1_spin_integration(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param) {
  CoordProxy _end_orb;
  fortran_track1_spin_integration(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr());
  return std::move(_end_orb);
}
CoordProxy Bmad::track1_spin_taylor(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param) {
  CoordProxy _end_orb;
  fortran_track1_spin_taylor(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr());
  return std::move(_end_orb);
}
void Bmad::track1_sr_wake(BunchProxy& bunch, EleProxy& ele) {
  fortran_track1_sr_wake(
      /* void* */ bunch.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
TrackProxy Bmad::track1_symp_lie_ptc(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param) {
  TrackProxy _track;
  fortran_track1_symp_lie_ptc(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _track.get_fortran_ptr());
  return std::move(_track);
}
FixedArray2D<Real, 6, 6> Bmad::track1_taylor(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<FixedArray1D<TaylorProxy, 6>> taylor,
    std::optional<bool> make_matrix) {
  auto* _taylor =
      taylor.has_value() ? taylor->data() : nullptr; // input, optional
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_taylor(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _taylor,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::Track1TimeRungeKutta Bmad::track1_time_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<double> t_end,
    optional_ref<double> dt_step) {
  bool _err_flag{};
  TrackProxy _track;
  double t_end_lvalue;
  auto* _t_end{&t_end_lvalue};
  if (t_end.has_value()) {
    t_end_lvalue = t_end.value();
  } else {
    _t_end = nullptr;
  }
  auto* _dt_step =
      dt_step.has_value() ? &dt_step->get() : nullptr; // inout, optional
  fortran_track1_time_runge_kutta(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _t_end,
      /* double* */ _dt_step);
  return Track1TimeRungeKutta{_err_flag, std::move(_track)};
}
Bmad::TrackABeambeam Bmad::track_a_beambeam(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  TrackProxy _track;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_beambeam(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* void* */ _track.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return TrackABeambeam{std::move(_track), mat6};
}
void Bmad::track_a_bend(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_bend(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_bend_photon(CoordProxy& orb, EleProxy& ele, double length) {
  fortran_track_a_bend_photon(
      /* void* */ orb.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ length);
}
void Bmad::track_a_capillary(CoordProxy& orb, EleProxy& ele) {
  fortran_track_a_capillary(
      /* void* */ orb.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
FixedArray2D<Real, 6, 6> Bmad::track_a_converter(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_converter(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
FixedArray2D<Real, 6, 6> Bmad::track_a_crab_cavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_crab_cavity(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_drift(
    CoordProxy& orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<int> ele_orientation,
    std::optional<bool> include_ref_motion,
    optional_ref<double> time) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  int ele_orientation_lvalue;
  auto* _ele_orientation{&ele_orientation_lvalue};
  if (ele_orientation.has_value()) {
    ele_orientation_lvalue = ele_orientation.value();
  } else {
    _ele_orientation = nullptr;
  }
  bool include_ref_motion_lvalue;
  auto* _include_ref_motion{&include_ref_motion_lvalue};
  if (include_ref_motion.has_value()) {
    include_ref_motion_lvalue = include_ref_motion.value();
  } else {
    _include_ref_motion = nullptr;
  }
  auto* _time = time.has_value() ? &time->get() : nullptr; // inout, optional
  fortran_track_a_drift(
      /* void* */ orb.get_fortran_ptr(),
      /* double& */ length,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix,
      /* int* */ _ele_orientation,
      /* bool* */ _include_ref_motion,
      /* double* */ _time);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_drift_photon(
    CoordProxy& orb,
    double length,
    bool phase_relative_to_ref) {
  fortran_track_a_drift_photon(
      /* void* */ orb.get_fortran_ptr(),
      /* double& */ length,
      /* bool& */ phase_relative_to_ref);
}
FixedArray2D<Real, 6, 6> Bmad::track_a_foil(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_foil(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_gkicker(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_gkicker(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_lcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_lcavity(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_lcavity_old(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_lcavity_old(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
FixedArray2D<Real, 6, 6> Bmad::track_a_mask(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_mask(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
FixedArray2D<Real, 6, 6> Bmad::track_a_match(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag,
    std::optional<bool> make_matrix) {
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_match(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool* */ _err_flag,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::TrackAPatch Bmad::track_a_patch(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit,
    std::optional<bool> track_spin,
    std::optional<bool> make_matrix) {
  bool drift_to_exit_lvalue;
  auto* _drift_to_exit{&drift_to_exit_lvalue};
  if (drift_to_exit.has_value()) {
    drift_to_exit_lvalue = drift_to_exit.value();
  } else {
    _drift_to_exit = nullptr;
  }
  double _s_ent{};
  double _ds_ref{};
  bool track_spin_lvalue;
  auto* _track_spin{&track_spin_lvalue};
  if (track_spin.has_value()) {
    track_spin_lvalue = track_spin.value();
  } else {
    _track_spin = nullptr;
  }
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_patch(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool* */ _drift_to_exit,
      /* double& */ _s_ent,
      /* double& */ _ds_ref,
      /* bool* */ _track_spin,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return TrackAPatch{_s_ent, _ds_ref, mat6};
}
void Bmad::track_a_patch_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit,
    std::optional<bool> use_z_pos) {
  bool drift_to_exit_lvalue;
  auto* _drift_to_exit{&drift_to_exit_lvalue};
  if (drift_to_exit.has_value()) {
    drift_to_exit_lvalue = drift_to_exit.value();
  } else {
    _drift_to_exit = nullptr;
  }
  bool use_z_pos_lvalue;
  auto* _use_z_pos{&use_z_pos_lvalue};
  if (use_z_pos.has_value()) {
    use_z_pos_lvalue = use_z_pos.value();
  } else {
    _use_z_pos = nullptr;
  }
  fortran_track_a_patch_photon(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool* */ _drift_to_exit,
      /* bool* */ _use_z_pos);
}
FixedArray2D<Real, 6, 6> Bmad::track_a_pickup(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag,
    std::optional<bool> make_matrix) {
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_pickup(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool* */ _err_flag,
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
FixedArray2D<Real, 6, 6> Bmad::track_a_quadrupole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_quadrupole(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
FixedArray2D<Real, 6, 6> Bmad::track_a_rfcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_rfcavity(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_sad_mult(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_sad_mult(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
FixedArray2D<Real, 6, 6> Bmad::track_a_sol_quad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_sol_quad(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_thick_multipole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix) {
  double _mat6_vec[6 * 6];
  const double* _mat6 = nullptr;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6 = _mat6_vec;
  }
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_thick_multipole(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
FixedArray2D<Real, 6, 6> Bmad::track_a_wiggler(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix) {
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  bool make_matrix_lvalue;
  auto* _make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_wiggler(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _mat6_vec,
      /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::TrackAZeroLengthElement Bmad::track_a_zero_length_element(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param) {
  bool _err_flag{};
  TrackProxy _track;
  fortran_track_a_zero_length_element(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _track.get_fortran_ptr());
  return TrackAZeroLengthElement{_err_flag, std::move(_track)};
}
Bmad::TrackAll Bmad::track_all(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    std::optional<int> ix_branch,
    std::optional<bool> init_lost) {
  // intent=inout allocatable type array
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  bool _err_flag{};
  // intent=out allocatable type array
  auto orbit0{CoordProxyAlloc1D()};
  bool init_lost_lvalue;
  auto* _init_lost{&init_lost_lvalue};
  if (init_lost.has_value()) {
    init_lost_lvalue = init_lost.value();
  } else {
    _init_lost = nullptr;
  }
  fortran_track_all(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _ix_branch,
      /* int& */ _track_state,
      /* bool& */ _err_flag,
      /* void* */ orbit0.get_fortran_ptr(),
      /* bool* */ _init_lost);
  return TrackAll{_track_state, _err_flag, std::move(orbit0)};
}
bool Bmad::track_beam(
    LatProxy& lat,
    BeamProxy& beam,
    optional_ref<EleProxy> ele1,
    optional_ref<EleProxy> ele2,
    optional_ref<CoordProxyAlloc1D> centroid,
    std::optional<int> direction,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks) {
  auto* _ele1 = ele1.has_value() ? ele1->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  auto* _ele2 = ele2.has_value() ? ele2->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  bool _err{};
  // intent=in allocatable type array
  auto* _centroid = centroid.has_value() ? centroid->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  // intent=inout allocatable type array
  auto* _bunch_tracks = bunch_tracks.has_value()
      ? bunch_tracks->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track_beam(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ beam.get_fortran_ptr(),
      /* void* */ _ele1,
      /* void* */ _ele2,
      /* bool& */ _err,
      /* void* */ _centroid,
      /* int* */ _direction,
      /* void* */ _bunch_tracks);
  return _err;
}
bool Bmad::track_bunch(
    LatProxy& lat,
    BunchProxy& bunch,
    optional_ref<EleProxy> ele1,
    optional_ref<EleProxy> ele2,
    optional_ref<CoordProxyAlloc1D> centroid,
    std::optional<int> direction,
    optional_ref<BunchTrackProxy> bunch_track) {
  auto* _ele1 = ele1.has_value() ? ele1->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  auto* _ele2 = ele2.has_value() ? ele2->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  bool _err{};
  // intent=in allocatable type array
  auto* _centroid = centroid.has_value() ? centroid->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  int direction_lvalue;
  auto* _direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto* _bunch_track = bunch_track.has_value()
      ? bunch_track->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track_bunch(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ _ele1,
      /* void* */ _ele2,
      /* bool& */ _err,
      /* void* */ _centroid,
      /* int* */ _direction,
      /* void* */ _bunch_track);
  return _err;
}
void Bmad::track_bunch_time(
    BunchProxy& bunch,
    BranchProxy& branch,
    double t_end,
    double s_end,
    optional_ref<RealAlloc1D> dt_step,
    optional_ref<EmFieldProxyAlloc1D> extra_field) {
  // intent=inout allocatable general array
  auto* _dt_step = dt_step.has_value() ? dt_step->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  // intent=in allocatable type array
  auto* _extra_field = extra_field.has_value()
      ? extra_field->get().get_fortran_ptr()
      : nullptr; // input, optional
  fortran_track_bunch_time(
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ t_end,
      /* double& */ s_end,
      /* void* */ _dt_step,
      /* void* */ _extra_field);
}
void Bmad::track_bunch_to_s(BunchProxy& bunch, double s, BranchProxy& branch) {
  fortran_track_bunch_to_s(
      /* void* */ bunch.get_fortran_ptr(),
      /* double& */ s,
      /* void* */ branch.get_fortran_ptr());
}
void Bmad::track_bunch_to_t(
    BunchProxy& bunch,
    double t_target,
    BranchProxy& branch) {
  fortran_track_bunch_to_t(
      /* void* */ bunch.get_fortran_ptr(),
      /* double& */ t_target,
      /* void* */ branch.get_fortran_ptr());
}
ComplexAlloc1D Bmad::track_complex_taylor(
    ComplexAlloc1D& start_orb,
    ComplexTaylorProxyAlloc1D& complex_taylor) {
  // intent=in allocatable general array
  // intent=in allocatable type array
  // intent=out allocatable general array
  auto end_orb{ComplexAlloc1D()};
  fortran_track_complex_taylor(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ complex_taylor.get_fortran_ptr(),
      /* void* */ end_orb.get_fortran_ptr());
  return std::move(end_orb);
}
Bmad::TrackFromSToS Bmad::track_from_s_to_s(
    LatProxy& lat,
    double s_start,
    double s_end,
    CoordProxy& orbit_start,
    std::optional<int> ix_branch,
    std::optional<int> ix_ele_end) {
  CoordProxy _orbit_end;
  // intent=out allocatable type array
  auto all_orb{CoordProxyAlloc1D()};
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  int ix_ele_end_lvalue;
  auto* _ix_ele_end{&ix_ele_end_lvalue};
  if (ix_ele_end.has_value()) {
    ix_ele_end_lvalue = ix_ele_end.value();
  } else {
    _ix_ele_end = nullptr;
  }
  fortran_track_from_s_to_s(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ s_start,
      /* double& */ s_end,
      /* void* */ orbit_start.get_fortran_ptr(),
      /* void* */ _orbit_end.get_fortran_ptr(),
      /* void* */ all_orb.get_fortran_ptr(),
      /* int* */ _ix_branch,
      /* int& */ _track_state,
      /* int* */ _ix_ele_end);
  return TrackFromSToS{std::move(_orbit_end), std::move(all_orb), _track_state};
}
int Bmad::track_many(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    int ix_start,
    int ix_end,
    int direction,
    std::optional<int> ix_branch) {
  // intent=inout allocatable type array
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  fortran_track_many(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ ix_start,
      /* int& */ ix_end,
      /* int& */ direction,
      /* int* */ _ix_branch,
      /* int& */ _track_state);
  return _track_state;
}
FixedArray2D<Real, 3, 3> Bmad::track_to_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    LatParamProxy& param) {
  FixedArray2D<Real, 3, 3> w_surface;
  double _w_surface_vec[3 * 3];
  fortran_track_to_surface(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double* */ _w_surface_vec);
  vec_to_matrix(_w_surface_vec, w_surface);
  return w_surface;
}
Bmad::TrackUntilDead Bmad::track_until_dead(
    CoordProxy& start_orb,
    LatProxy& lat) {
  CoordProxy _end_orb;
  TrackProxy _track;
  fortran_track_until_dead(
      /* void* */ start_orb.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _end_orb.get_fortran_ptr(),
      /* void* */ _track.get_fortran_ptr());
  return TrackUntilDead{std::move(_end_orb), std::move(_track)};
}
Bmad::TrackingRadMapSetup Bmad::tracking_rad_map_setup(
    EleProxy& ele,
    double tollerance,
    int ref_edge) {
  RadMapProxy _rad_map;
  bool _err_flag{};
  fortran_tracking_rad_map_setup(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ tollerance,
      /* int& */ ref_edge,
      /* void* */ _rad_map.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return TrackingRadMapSetup{std::move(_rad_map), _err_flag};
}
AcKickerProxy Bmad::transfer_ac_kick(AcKickerProxy& ac_in) {
  auto _ac_in = &ac_in; // input, required, pointer
  AcKickerProxy _ac_out;
  fortran_transfer_ac_kick(
      /* void* */ &ac_in, /* void* */ _ac_out.get_fortran_ptr());
  return std::move(_ac_out);
}
BranchProxy Bmad::transfer_branch(BranchProxy& branch1) {
  BranchProxy _branch2;
  fortran_transfer_branch(
      /* void* */ branch1.get_fortran_ptr(),
      /* void* */ _branch2.get_fortran_ptr());
  return std::move(_branch2);
}
BranchProxy Bmad::transfer_branch_parameters(BranchProxy& branch_in) {
  BranchProxy _branch_out;
  fortran_transfer_branch_parameters(
      /* void* */ branch_in.get_fortran_ptr(),
      /* void* */ _branch_out.get_fortran_ptr());
  return std::move(_branch_out);
}
BranchProxyAlloc1D Bmad::transfer_branches(BranchProxyAlloc1D& branch1) {
  // intent=in allocatable type array
  // intent=out allocatable type array
  auto branch2{BranchProxyAlloc1D()};
  fortran_transfer_branches(
      /* void* */ branch1.get_fortran_ptr(),
      /* void* */ branch2.get_fortran_ptr());
  return std::move(branch2);
}
EleProxy Bmad::transfer_ele(
    EleProxy& ele1,
    std::optional<bool> nullify_pointers) {
  EleProxy _ele2;
  bool nullify_pointers_lvalue;
  auto* _nullify_pointers{&nullify_pointers_lvalue};
  if (nullify_pointers.has_value()) {
    nullify_pointers_lvalue = nullify_pointers.value();
  } else {
    _nullify_pointers = nullptr;
  }
  fortran_transfer_ele(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ _ele2.get_fortran_ptr(),
      /* bool* */ _nullify_pointers);
  return std::move(_ele2);
}
EleProxy Bmad::transfer_ele_taylor(
    EleProxy& ele_in,
    std::optional<int> taylor_order) {
  EleProxy _ele_out;
  int taylor_order_lvalue;
  auto* _taylor_order{&taylor_order_lvalue};
  if (taylor_order.has_value()) {
    taylor_order_lvalue = taylor_order.value();
  } else {
    _taylor_order = nullptr;
  }
  fortran_transfer_ele_taylor(
      /* void* */ ele_in.get_fortran_ptr(),
      /* void* */ _ele_out.get_fortran_ptr(),
      /* int* */ _taylor_order);
  return std::move(_ele_out);
}
EleProxyAlloc1D Bmad::transfer_eles(EleProxyAlloc1D& ele1) {
  // intent=in allocatable type array
  // intent=out allocatable type array
  auto ele2{EleProxyAlloc1D()};
  fortran_transfer_eles(
      /* void* */ ele1.get_fortran_ptr(), /* void* */ ele2.get_fortran_ptr());
  return std::move(ele2);
}
EleProxy Bmad::transfer_fieldmap(EleProxy& ele_in, int who) {
  EleProxy _ele_out;
  fortran_transfer_fieldmap(
      /* void* */ ele_in.get_fortran_ptr(),
      /* void* */ _ele_out.get_fortran_ptr(),
      /* int& */ who);
  return std::move(_ele_out);
}
bool Bmad::transfer_fixer_params(
    EleProxy& fixer,
    bool to_stored,
    optional_ref<CoordProxy> orbit,
    std::optional<std::string> who) {
  auto* _orbit = orbit.has_value() ? orbit->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  const char* _who = who.has_value() ? who->c_str() : nullptr;
  bool _is_ok{};
  fortran_transfer_fixer_params(
      /* void* */ fixer.get_fortran_ptr(),
      /* bool& */ to_stored,
      /* void* */ _orbit,
      /* const char* */ _who,
      /* bool& */ _is_ok);
  return _is_ok;
}
LatProxy Bmad::transfer_lat(LatProxy& lat1) {
  LatProxy _lat2;
  fortran_transfer_lat(
      /* void* */ lat1.get_fortran_ptr(), /* void* */ _lat2.get_fortran_ptr());
  return std::move(_lat2);
}
LatProxy Bmad::transfer_lat_parameters(LatProxy& lat_in) {
  LatProxy _lat_out;
  fortran_transfer_lat_parameters(
      /* void* */ lat_in.get_fortran_ptr(),
      /* void* */ _lat_out.get_fortran_ptr());
  return std::move(_lat_out);
}
bool Bmad::transfer_map_calc(
    LatProxy& lat,
    TaylorProxyAlloc1D& orb_map,
    std::optional<int> ix1,
    std::optional<int> ix2,
    optional_ref<CoordProxy> ref_orb,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start,
    std::optional<bool> concat_if_possible,
    optional_ref<TaylorProxyAlloc1D> spin_map) {
  // intent=inout allocatable type array
  bool _err_flag{};
  int ix1_lvalue;
  auto* _ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto* _ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  auto* _ref_orb = ref_orb.has_value() ? ref_orb->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto* _one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto* _unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  bool concat_if_possible_lvalue;
  auto* _concat_if_possible{&concat_if_possible_lvalue};
  if (concat_if_possible.has_value()) {
    concat_if_possible_lvalue = concat_if_possible.value();
  } else {
    _concat_if_possible = nullptr;
  }
  // intent=inout allocatable type array
  auto* _spin_map = spin_map.has_value() ? spin_map->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  fortran_transfer_map_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orb_map.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* int* */ _ix1,
      /* int* */ _ix2,
      /* void* */ _ref_orb,
      /* int* */ _ix_branch,
      /* bool* */ _one_turn,
      /* bool* */ _unit_start,
      /* bool* */ _concat_if_possible,
      /* void* */ _spin_map);
  return _err_flag;
}
Bmad::TransferMapFromSToS Bmad::transfer_map_from_s_to_s(
    LatProxy& lat,
    TaylorProxyAlloc1D& t_map,
    std::optional<double> s1,
    std::optional<double> s2,
    optional_ref<CoordProxy> ref_orb_in,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start,
    std::optional<bool> concat_if_possible,
    optional_ref<TaylorProxyAlloc1D> spin_map) {
  // intent=inout allocatable type array
  double s1_lvalue;
  auto* _s1{&s1_lvalue};
  if (s1.has_value()) {
    s1_lvalue = s1.value();
  } else {
    _s1 = nullptr;
  }
  double s2_lvalue;
  auto* _s2{&s2_lvalue};
  if (s2.has_value()) {
    s2_lvalue = s2.value();
  } else {
    _s2 = nullptr;
  }
  auto* _ref_orb_in = ref_orb_in.has_value()
      ? ref_orb_in->get().get_fortran_ptr()
      : nullptr; // input, optional
  CoordProxy _ref_orb_out;
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto* _one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto* _unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  bool _err_flag{};
  bool concat_if_possible_lvalue;
  auto* _concat_if_possible{&concat_if_possible_lvalue};
  if (concat_if_possible.has_value()) {
    concat_if_possible_lvalue = concat_if_possible.value();
  } else {
    _concat_if_possible = nullptr;
  }
  // intent=inout allocatable type array
  auto* _spin_map = spin_map.has_value() ? spin_map->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  fortran_transfer_map_from_s_to_s(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ t_map.get_fortran_ptr(),
      /* double* */ _s1,
      /* double* */ _s2,
      /* void* */ _ref_orb_in,
      /* void* */ _ref_orb_out.get_fortran_ptr(),
      /* int* */ _ix_branch,
      /* bool* */ _one_turn,
      /* bool* */ _unit_start,
      /* bool& */ _err_flag,
      /* bool* */ _concat_if_possible,
      /* void* */ _spin_map);
  return TransferMapFromSToS{std::move(_ref_orb_out), _err_flag};
}
FixedArray2D<Real, 2, 2> Bmad::transfer_mat2_from_twiss(
    TwissProxy& twiss1,
    TwissProxy& twiss2) {
  FixedArray2D<Real, 2, 2> mat;
  double _mat_vec[2 * 2];
  fortran_transfer_mat2_from_twiss(
      /* void* */ twiss1.get_fortran_ptr(),
      /* void* */ twiss2.get_fortran_ptr(),
      /* double* */ _mat_vec);
  vec_to_matrix(_mat_vec, mat);
  return mat;
}
FixedArray2D<Real, 6, 6> Bmad::transfer_mat_from_twiss(
    EleProxy& ele1,
    EleProxy& ele2,
    FixedArray1D<Real, 6> orb1,
    FixedArray1D<Real, 6> orb2) {
  auto* _orb1 = orb1.data(); // CppWrapperGeneralArgument
  auto* _orb2 = orb2.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 6, 6> m;
  double _m_vec[6 * 6];
  fortran_transfer_mat_from_twiss(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* double* */ _orb1,
      /* double* */ _orb2,
      /* double* */ _m_vec);
  vec_to_matrix(_m_vec, m);
  return m;
}
void Bmad::transfer_matrix_calc(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> xfer_mat,
    std::optional<FixedArray1D<Real, 6>> xfer_vec,
    std::optional<int> ix1,
    std::optional<int> ix2,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn) {
  double _xfer_mat_vec[6 * 6];
  matrix_to_vec(xfer_mat, _xfer_mat_vec);
  double* _xfer_vec = xfer_vec.has_value() ? xfer_vec.value().data() : nullptr;
  int ix1_lvalue;
  auto* _ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto* _ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto* _one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  fortran_transfer_matrix_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* double* */ _xfer_mat_vec,
      /* double* */ _xfer_vec,
      /* int* */ _ix1,
      /* int* */ _ix2,
      /* int* */ _ix_branch,
      /* bool* */ _one_turn);
  vec_to_matrix(_xfer_mat_vec, xfer_mat);
}
EleProxy Bmad::transfer_twiss(EleProxy& ele_in, std::optional<bool> reverse) {
  EleProxy _ele_out;
  bool reverse_lvalue;
  auto* _reverse{&reverse_lvalue};
  if (reverse.has_value()) {
    reverse_lvalue = reverse.value();
  } else {
    _reverse = nullptr;
  }
  fortran_transfer_twiss(
      /* void* */ ele_in.get_fortran_ptr(),
      /* void* */ _ele_out.get_fortran_ptr(),
      /* bool* */ _reverse);
  return std::move(_ele_out);
}
WakeProxy Bmad::transfer_wake(WakeProxy& wake_in) {
  auto _wake_in = &wake_in; // input, required, pointer
  WakeProxy _wake_out;
  fortran_transfer_wake(
      /* void* */ &wake_in, /* void* */ _wake_out.get_fortran_ptr());
  return std::move(_wake_out);
}
ComplexTaylorProxyAlloc1D Bmad::truncate_complex_taylor_to_order(
    ComplexTaylorProxyAlloc1D& complex_taylor_in,
    int order) {
  // intent=in allocatable type array
  // intent=out allocatable type array
  auto complex_taylor_out{ComplexTaylorProxyAlloc1D()};
  fortran_truncate_complex_taylor_to_order(
      /* void* */ complex_taylor_in.get_fortran_ptr(),
      /* int& */ order,
      /* void* */ complex_taylor_out.get_fortran_ptr());
  return std::move(complex_taylor_out);
}
Bmad::Twiss1Propagate Bmad::twiss1_propagate(
    TwissProxy& twiss1,
    FixedArray2D<Real, 2, 2> mat2,
    int ele_key,
    double length) {
  double _mat2_vec[2 * 2];
  matrix_to_vec(mat2, _mat2_vec);
  TwissProxy _twiss2;
  bool _err{};
  fortran_twiss1_propagate(
      /* void* */ twiss1.get_fortran_ptr(),
      /* double* */ _mat2_vec,
      /* int& */ ele_key,
      /* double& */ length,
      /* void* */ _twiss2.get_fortran_ptr(),
      /* bool& */ _err);
  return Twiss1Propagate{std::move(_twiss2), _err};
}
FixedArray1D<Real, 3> Bmad::twiss3_at_start(
    LatProxy& lat,
    bool& err_flag,
    std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  FixedArray1D<Real, 3> _tune3;
  fortran_twiss3_at_start(
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ err_flag,
      /* int* */ _ix_branch,
      /* double* */ _tune3.data());
  return _tune3;
}
void Bmad::twiss3_from_twiss2(EleProxy& ele) {
  fortran_twiss3_from_twiss2(/* void* */ ele.get_fortran_ptr());
}
void Bmad::twiss3_propagate1(EleProxy& ele1, EleProxy& ele2, bool& err_flag) {
  fortran_twiss3_propagate1(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::twiss3_propagate_all(LatProxy& lat, std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  fortran_twiss3_propagate_all(
      /* void* */ lat.get_fortran_ptr(), /* int* */ _ix_branch);
}
int Bmad::twiss_and_track(
    LatProxy& lat,
    CoordArrayProxyAlloc1D& orb_array,
    std::optional<bool> print_err,
    std::optional<bool> calc_chrom) {
  // intent=inout allocatable type array
  int _status{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool calc_chrom_lvalue;
  auto* _calc_chrom{&calc_chrom_lvalue};
  if (calc_chrom.has_value()) {
    calc_chrom_lvalue = calc_chrom.value();
  } else {
    _calc_chrom = nullptr;
  }
  fortran_twiss_and_track_all(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orb_array.get_fortran_ptr(),
      /* int& */ _status,
      /* bool* */ _print_err,
      /* bool* */ _calc_chrom);
  return _status;
}
bool Bmad::twiss_and_track_at_s(
    LatProxy& lat,
    double s,
    optional_ref<EleProxy> ele_at_s,
    optional_ref<CoordProxyAlloc1D> orb,
    optional_ref<CoordProxy> orb_at_s,
    std::optional<int> ix_branch,
    std::optional<bool> use_last,
    std::optional<bool> compute_floor_coords) {
  auto* _ele_at_s = ele_at_s.has_value() ? ele_at_s->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  // intent=in allocatable type array
  auto* _orb = orb.has_value() ? orb->get().get_fortran_ptr()
                               : nullptr; // input, optional
  auto* _orb_at_s = orb_at_s.has_value() ? orb_at_s->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err{};
  bool use_last_lvalue;
  auto* _use_last{&use_last_lvalue};
  if (use_last.has_value()) {
    use_last_lvalue = use_last.value();
  } else {
    _use_last = nullptr;
  }
  bool compute_floor_coords_lvalue;
  auto* _compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  fortran_twiss_and_track_at_s(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ s,
      /* void* */ _ele_at_s,
      /* void* */ _orb,
      /* void* */ _orb_at_s,
      /* int* */ _ix_branch,
      /* bool& */ _err,
      /* bool* */ _use_last,
      /* bool* */ _compute_floor_coords);
  return _err;
}
int Bmad::twiss_and_track(
    LatProxy& lat,
    CoordProxyAlloc1D& orb,
    std::optional<int> ix_branch,
    std::optional<bool> print_err,
    std::optional<bool> calc_chrom,
    optional_ref<CoordProxy> orb_start) {
  // intent=inout allocatable type array
  int _status{};
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool calc_chrom_lvalue;
  auto* _calc_chrom{&calc_chrom_lvalue};
  if (calc_chrom.has_value()) {
    calc_chrom_lvalue = calc_chrom.value();
  } else {
    _calc_chrom = nullptr;
  }
  auto* _orb_start = orb_start.has_value() ? orb_start->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  fortran_twiss_and_track_branch(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ orb.get_fortran_ptr(),
      /* int& */ _status,
      /* int* */ _ix_branch,
      /* bool* */ _print_err,
      /* bool* */ _calc_chrom,
      /* void* */ _orb_start);
  return _status;
}
Bmad::TwissAndTrackFromSToS Bmad::twiss_and_track_from_s_to_s(
    BranchProxy& branch,
    CoordProxy& orbit_start,
    double s_end,
    optional_ref<EleProxy> ele_start,
    std::optional<bool> compute_floor_coords,
    std::optional<bool> compute_twiss) {
  CoordProxy _orbit_end;
  auto* _ele_start = ele_start.has_value() ? ele_start->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  EleProxy _ele_end;
  bool _err{};
  bool compute_floor_coords_lvalue;
  auto* _compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  bool compute_twiss_lvalue;
  auto* _compute_twiss{&compute_twiss_lvalue};
  if (compute_twiss.has_value()) {
    compute_twiss_lvalue = compute_twiss.value();
  } else {
    _compute_twiss = nullptr;
  }
  fortran_twiss_and_track_from_s_to_s(
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ orbit_start.get_fortran_ptr(),
      /* double& */ s_end,
      /* void* */ _orbit_end.get_fortran_ptr(),
      /* void* */ _ele_start,
      /* void* */ _ele_end.get_fortran_ptr(),
      /* bool& */ _err,
      /* bool* */ _compute_floor_coords,
      /* bool* */ _compute_twiss);
  return TwissAndTrackFromSToS{
      std::move(_orbit_end), std::move(_ele_end), _err};
}
Bmad::TwissAndTrackIntraEle Bmad::twiss_and_track_intra_ele(
    EleProxy& ele,
    LatParamProxy& param,
    double l_start,
    double l_end,
    bool track_upstream_end,
    bool track_downstream_end,
    optional_ref<CoordProxy> orbit_start,
    optional_ref<EleProxy> ele_start,
    optional_ref<EleProxy> ele_end,
    std::optional<bool> compute_floor_coords,
    std::optional<bool> compute_twiss,
    std::optional<bool> reuse_ele_end) {
  auto* _orbit_start = orbit_start.has_value()
      ? orbit_start->get().get_fortran_ptr()
      : nullptr; // input, optional
  CoordProxy _orbit_end;
  auto* _ele_start = ele_start.has_value() ? ele_start->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  auto* _ele_end = ele_end.has_value() ? ele_end->get().get_fortran_ptr()
                                       : nullptr; // input, optional
  bool _err{};
  bool compute_floor_coords_lvalue;
  auto* _compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  bool compute_twiss_lvalue;
  auto* _compute_twiss{&compute_twiss_lvalue};
  if (compute_twiss.has_value()) {
    compute_twiss_lvalue = compute_twiss.value();
  } else {
    _compute_twiss = nullptr;
  }
  bool reuse_ele_end_lvalue;
  auto* _reuse_ele_end{&reuse_ele_end_lvalue};
  if (reuse_ele_end.has_value()) {
    reuse_ele_end_lvalue = reuse_ele_end.value();
  } else {
    _reuse_ele_end = nullptr;
  }
  fortran_twiss_and_track_intra_ele(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ l_start,
      /* double& */ l_end,
      /* bool& */ track_upstream_end,
      /* bool& */ track_downstream_end,
      /* void* */ _orbit_start,
      /* void* */ _orbit_end.get_fortran_ptr(),
      /* void* */ _ele_start,
      /* void* */ _ele_end,
      /* bool& */ _err,
      /* bool* */ _compute_floor_coords,
      /* bool* */ _compute_twiss,
      /* bool* */ _reuse_ele_end);
  return TwissAndTrackIntraEle{std::move(_orbit_end), _err};
}
Bmad::TwissAtElement Bmad::twiss_at_element(EleProxy& ele) {
  EleProxy _start;
  EleProxy _end;
  EleProxy _average;
  fortran_twiss_at_element(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _start.get_fortran_ptr(),
      /* void* */ _end.get_fortran_ptr(),
      /* void* */ _average.get_fortran_ptr());
  return TwissAtElement{
      std::move(_start), std::move(_end), std::move(_average)};
}
int Bmad::twiss_at_start(
    LatProxy& lat,
    std::optional<int> ix_branch,
    std::optional<bool> type_out) {
  int _status{};
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool type_out_lvalue;
  auto* _type_out{&type_out_lvalue};
  if (type_out.has_value()) {
    type_out_lvalue = type_out.value();
  } else {
    _type_out = nullptr;
  }
  fortran_twiss_at_start(
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ _status,
      /* int* */ _ix_branch,
      /* bool* */ _type_out);
  return _status;
}
Bmad::TwissFromTracking Bmad::twiss_from_tracking(
    LatProxy& lat,
    CoordProxy& ref_orb0,
    optional_ref<RealAlloc1D> d_orb) {
  double _symp_err{};
  bool _err_flag{};
  // intent=in allocatable general array
  auto* _d_orb = d_orb.has_value() ? d_orb->get().get_fortran_ptr()
                                   : nullptr; // input, optional
  fortran_twiss_from_tracking(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ ref_orb0.get_fortran_ptr(),
      /* double& */ _symp_err,
      /* bool& */ _err_flag,
      /* void* */ _d_orb);
  return TwissFromTracking{_symp_err, _err_flag};
}
bool Bmad::twiss_propagate1(
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> forward) {
  bool _err_flag{};
  bool forward_lvalue;
  auto* _forward{&forward_lvalue};
  if (forward.has_value()) {
    forward_lvalue = forward.value();
  } else {
    _forward = nullptr;
  }
  fortran_twiss_propagate1(
      /* void* */ ele1.get_fortran_ptr(),
      /* void* */ ele2.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _forward);
  return _err_flag;
}
bool Bmad::twiss_propagate_all(
    LatProxy& lat,
    std::optional<int> ix_branch,
    std::optional<int> ie_start,
    std::optional<int> ie_end) {
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  int ie_start_lvalue;
  auto* _ie_start{&ie_start_lvalue};
  if (ie_start.has_value()) {
    ie_start_lvalue = ie_start.value();
  } else {
    _ie_start = nullptr;
  }
  int ie_end_lvalue;
  auto* _ie_end{&ie_end_lvalue};
  if (ie_end.has_value()) {
    ie_end_lvalue = ie_end.value();
  } else {
    _ie_end = nullptr;
  }
  fortran_twiss_propagate_all(
      /* void* */ lat.get_fortran_ptr(),
      /* int* */ _ix_branch,
      /* bool& */ _err_flag,
      /* int* */ _ie_start,
      /* int* */ _ie_end);
  return _err_flag;
}
FixedArray2D<Real, 2, 2> Bmad::twiss_to_1_turn_mat(
    TwissProxy& twiss,
    double phi) {
  FixedArray2D<Real, 2, 2> mat2;
  double _mat2_vec[2 * 2];
  fortran_twiss_to_1_turn_mat(
      /* void* */ twiss.get_fortran_ptr(),
      /* double& */ phi,
      /* double* */ _mat2_vec);
  vec_to_matrix(_mat2_vec, mat2);
  return mat2;
}
void Bmad::type_coord(CoordProxy& coord) {
  fortran_type_coord(/* void* */ coord.get_fortran_ptr());
}
void Bmad::type_expression_tree(
    ExpressionTreeProxy& tree,
    std::optional<int> indent) {
  int indent_lvalue;
  auto* _indent{&indent_lvalue};
  if (indent.has_value()) {
    indent_lvalue = indent.value();
  } else {
    _indent = nullptr;
  }
  fortran_type_expression_tree(
      /* void* */ tree.get_fortran_ptr(), /* int* */ _indent);
}
void Bmad::update_ele_from_fibre(EleProxy& ele) {
  fortran_update_ele_from_fibre(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::update_fibre_from_ele(EleProxy& ele) {
  bool _survey_needed{};
  fortran_update_fibre_from_ele(
      /* void* */ ele.get_fortran_ptr(), /* bool& */ _survey_needed);
  return _survey_needed;
}
void Bmad::update_floor_angles(
    FloorPositionProxy& floor,
    optional_ref<FloorPositionProxy> floor0) {
  auto* _floor0 = floor0.has_value() ? floor0->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  fortran_update_floor_angles(
      /* void* */ floor.get_fortran_ptr(), /* void* */ _floor0);
}
void Bmad::valid_field_calc(EleProxy& ele, int field_calc, bool& is_valid) {
  fortran_valid_field_calc(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ field_calc,
      /* bool& */ is_valid);
}
void Bmad::valid_fringe_type(EleProxy& ele, int fringe_type, bool& is_valid) {
  fortran_valid_fringe_type(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ fringe_type,
      /* bool& */ is_valid);
}
void Bmad::valid_mat6_calc_method(
    EleProxy& ele,
    int species,
    int mat6_calc_method,
    bool& is_valid) {
  fortran_valid_mat6_calc_method(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ species,
      /* int& */ mat6_calc_method,
      /* bool& */ is_valid);
}
void Bmad::valid_spin_tracking_method(
    EleProxy& ele,
    int spin_tracking_method,
    bool& is_valid) {
  fortran_valid_spin_tracking_method(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ spin_tracking_method,
      /* bool& */ is_valid);
}
void Bmad::valid_tracking_method(
    EleProxy& ele,
    int species,
    int tracking_method,
    bool& is_valid) {
  fortran_valid_tracking_method(
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ species,
      /* int& */ tracking_method,
      /* bool& */ is_valid);
}
bool Bmad::value_of_attribute(
    EleProxy& ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<double> err_value,
    double& value) {
  auto _attrib_name = attrib_name.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto* _err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  double err_value_lvalue;
  auto* _err_value{&err_value_lvalue};
  if (err_value.has_value()) {
    err_value_lvalue = err_value.value();
  } else {
    _err_value = nullptr;
  }
  fortran_value_of_attribute(
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_name,
      /* bool& */ _err_flag,
      /* bool* */ _err_print_flag,
      /* double* */ _err_value,
      /* double& */ value);
  return _err_flag;
}
void Bmad::value_to_line(
    std::string& line,
    double& value,
    std::string& str,
    std::string& typ,
    optional_ref<bool> ignore_if_zero,
    optional_ref<bool> use_comma) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _str = str.c_str(); // ptr, inout, required
  auto _typ = typ.c_str(); // ptr, inout, required
  auto* _ignore_if_zero = ignore_if_zero.has_value()
      ? &ignore_if_zero->get()
      : nullptr; // inout, optional
  auto* _use_comma =
      use_comma.has_value() ? &use_comma->get() : nullptr; // inout, optional
  fortran_value_to_line(
      /* const char* */ _line,
      /* double& */ value,
      /* const char* */ _str,
      /* const char* */ _typ,
      /* bool* */ _ignore_if_zero,
      /* bool* */ _use_comma);
}
void Bmad::vec_to_polar(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    SpinPolarProxy& polar) {
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  double phase_lvalue;
  auto* _phase{&phase_lvalue};
  if (phase.has_value()) {
    phase_lvalue = phase.value();
  } else {
    _phase = nullptr;
  }
  fortran_vec_to_polar(
      /* double* */ _vec,
      /* double* */ _phase,
      /* void* */ polar.get_fortran_ptr());
}
void Bmad::vec_to_spinor(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    FixedArray1D<Complex, 2> spinor) {
  auto* _vec = vec.data(); // CppWrapperGeneralArgument
  double phase_lvalue;
  auto* _phase{&phase_lvalue};
  if (phase.has_value()) {
    phase_lvalue = phase.value();
  } else {
    _phase = nullptr;
  }
  auto* _spinor = spinor.data(); // CppWrapperGeneralArgument
  fortran_vec_to_spinor(
      /* double* */ _vec,
      /* double* */ _phase,
      /* std::complex<double>* */ _spinor);
}
bool Bmad::verify_valid_name(
    std::string name,
    int ix_name,
    std::optional<bool> pure_name,
    std::optional<bool> include_wild) {
  auto _name = name.c_str();
  bool pure_name_lvalue;
  auto* _pure_name{&pure_name_lvalue};
  if (pure_name.has_value()) {
    pure_name_lvalue = pure_name.value();
  } else {
    _pure_name = nullptr;
  }
  bool include_wild_lvalue;
  auto* _include_wild{&include_wild_lvalue};
  if (include_wild.has_value()) {
    include_wild_lvalue = include_wild.value();
  } else {
    _include_wild = nullptr;
  }
  bool _is_valid{};
  fortran_verify_valid_name(
      /* const char* */ _name,
      /* int& */ ix_name,
      /* bool* */ _pure_name,
      /* bool* */ _include_wild,
      /* bool& */ _is_valid);
  return _is_valid;
}
void Bmad::w_mat_for_bend_angle(
    double angle,
    double ref_tilt,
    std::optional<FixedArray1D<Real, 3>> r_vec,
    FixedArray2D<Real, 3, 3> w_mat) {
  double* _r_vec = r_vec.has_value() ? r_vec.value().data() : nullptr;
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  fortran_w_mat_for_bend_angle(
      /* double& */ angle,
      /* double& */ ref_tilt,
      /* double* */ _r_vec,
      /* double* */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
}
void Bmad::w_mat_for_tilt(
    double tilt,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat) {
  bool return_inverse_lvalue;
  auto* _return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  fortran_w_mat_for_tilt(
      /* double& */ tilt,
      /* bool* */ _return_inverse,
      /* double* */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
}
void Bmad::w_mat_for_x_pitch(
    double x_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat) {
  bool return_inverse_lvalue;
  auto* _return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  fortran_w_mat_for_x_pitch(
      /* double& */ x_pitch,
      /* bool* */ _return_inverse,
      /* double* */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
}
void Bmad::w_mat_for_y_pitch(
    double y_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat) {
  bool return_inverse_lvalue;
  auto* _return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  fortran_w_mat_for_y_pitch(
      /* double& */ y_pitch,
      /* bool* */ _return_inverse,
      /* double* */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
}
Bmad::Wall3dDRadius Bmad::wall3d_d_radius(
    RealAlloc1D& position,
    EleProxy& ele,
    std::optional<int> ix_wall) {
  // intent=in allocatable general array
  int ix_wall_lvalue;
  auto* _ix_wall{&ix_wall_lvalue};
  if (ix_wall.has_value()) {
    ix_wall_lvalue = ix_wall.value();
  } else {
    _ix_wall = nullptr;
  }
  FixedArray1D<Real, 3> _perp;
  int _ix_section{};
  bool _no_wall_here{};
  FixedArray1D<Real, 3> _origin;
  double _radius_wall{};
  bool _err_flag{};
  double _d_radius{};
  fortran_wall3d_d_radius(
      /* void* */ position.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _ix_wall,
      /* double* */ _perp.data(),
      /* int& */ _ix_section,
      /* bool& */ _no_wall_here,
      /* double* */ _origin.data(),
      /* double& */ _radius_wall,
      /* bool& */ _err_flag,
      /* double& */ _d_radius);
  return Wall3dDRadius{
      _perp,
      _ix_section,
      _no_wall_here,
      _origin,
      _radius_wall,
      _err_flag,
      _d_radius};
}
bool Bmad::wall3d_initializer(Wall3dProxy& wall3d) {
  bool _err{};
  fortran_wall3d_initializer(
      /* void* */ wall3d.get_fortran_ptr(), /* bool& */ _err);
  return _err;
}
bool Bmad::wall3d_section_initializer(Wall3dSectionProxy& section) {
  bool _err{};
  fortran_wall3d_section_initializer(
      /* void* */ section.get_fortran_ptr(), /* bool& */ _err);
  return _err;
}
FixedArray1D<Real, 6> Bmad::wall3d_to_position(
    CoordProxy& orbit,
    EleProxy& ele) {
  FixedArray1D<Real, 6> _position;
  fortran_wall3d_to_position(
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double* */ _position.data());
  return _position;
}
void Bmad::word_to_value(
    std::string& word,
    LatProxy& lat,
    double& value,
    bool& err_flag,
    optional_ref<EleProxy> ele) {
  auto _word = word.c_str(); // ptr, inout, required
  auto* _ele = ele.has_value() ? ele->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_word_to_value(
      /* const char* */ _word,
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ value,
      /* bool& */ err_flag,
      /* void* */ _ele);
}
void Bmad::write_ascii_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file,
    std::optional<bool> alive_only) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto* _new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  bool alive_only_lvalue;
  auto* _alive_only{&alive_only_lvalue};
  if (alive_only.has_value()) {
    alive_only_lvalue = alive_only.value();
  } else {
    _alive_only = nullptr;
  }
  fortran_write_ascii_beam_file(
      /* const char* */ _file_name,
      /* void* */ beam.get_fortran_ptr(),
      /* bool* */ _new_file,
      /* bool* */ _alive_only);
}
void Bmad::write_astra_bend(
    int& iu,
    double& strength,
    int& id,
    FixedArray1D<Real, 2> d1,
    FixedArray1D<Real, 2> d2,
    FixedArray1D<Real, 2> d3,
    FixedArray1D<Real, 2> d4) {
  auto* _d1 = d1.data(); // CppWrapperGeneralArgument
  auto* _d2 = d2.data(); // CppWrapperGeneralArgument
  auto* _d3 = d3.data(); // CppWrapperGeneralArgument
  auto* _d4 = d4.data(); // CppWrapperGeneralArgument
  fortran_write_astra_bend(
      /* int& */ iu,
      /* double& */ strength,
      /* int& */ id,
      /* double* */ _d1,
      /* double* */ _d2,
      /* double* */ _d3,
      /* double* */ _d4);
}
Bmad::WriteAstraFieldGridFile Bmad::write_astra_field_grid_file(
    int astra_file_unit,
    EleProxy& ele,
    std::optional<double> dz) {
  double _maxfield{};
  double dz_lvalue;
  auto* _dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  bool _err{};
  fortran_write_astra_field_grid_file(
      /* int& */ astra_file_unit,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* double* */ _dz,
      /* bool& */ _err);
  return WriteAstraFieldGridFile{_maxfield, _err};
}
Bmad::WriteAstraFieldGridFile3d Bmad::write_astra_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz) {
  auto _base_filename = base_filename.c_str();
  double _maxfield{};
  double dz_lvalue;
  auto* _dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  bool _err{};
  fortran_write_astra_field_grid_file_3d(
      /* const char* */ _base_filename,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* double* */ _dz,
      /* bool& */ _err);
  return WriteAstraFieldGridFile3d{_maxfield, _err};
}
void Bmad::write_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file,
    std::optional<int> file_format,
    optional_ref<LatProxy> lat,
    std::optional<bool> alive_only) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto* _new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  int file_format_lvalue;
  auto* _file_format{&file_format_lvalue};
  if (file_format.has_value()) {
    file_format_lvalue = file_format.value();
  } else {
    _file_format = nullptr;
  }
  auto* _lat = lat.has_value() ? lat->get().get_fortran_ptr()
                               : nullptr; // input, optional
  bool alive_only_lvalue;
  auto* _alive_only{&alive_only_lvalue};
  if (alive_only.has_value()) {
    alive_only_lvalue = alive_only.value();
  } else {
    _alive_only = nullptr;
  }
  fortran_write_beam_file(
      /* const char* */ _file_name,
      /* void* */ beam.get_fortran_ptr(),
      /* bool* */ _new_file,
      /* int* */ _file_format,
      /* void* */ _lat,
      /* bool* */ _alive_only);
}
void Bmad::write_beam_floor_positions(
    std::string file_name,
    BeamProxy& beam,
    EleProxy& ele,
    std::optional<bool> new_file) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto* _new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  fortran_write_beam_floor_positions(
      /* const char* */ _file_name,
      /* void* */ beam.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _new_file);
}
void Bmad::write_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_write_binary_cartesian_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cart_map.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::write_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_write_binary_cylindrical_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cl_map.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::write_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_write_binary_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ g_field.get_fortran_ptr(),
      /* bool& */ err_flag);
}
void Bmad::write_blender_ele(
    int& iu,
    EleProxy& ele,
    optional_ref<bool> old_format) {
  auto* _old_format =
      old_format.has_value() ? &old_format->get() : nullptr; // inout, optional
  fortran_write_blender_ele(
      /* int& */ iu,
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _old_format);
}
void Bmad::write_blender_lat_layout(std::string& file_name, LatProxy& lat) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_write_blender_lat_layout(
      /* const char* */ _file_name, /* void* */ lat.get_fortran_ptr());
}
bool Bmad::write_bmad_lattice_file(
    std::string bmad_file,
    LatProxy& lat,
    std::optional<int> output_form,
    optional_ref<CoordProxy> orbit0) {
  auto _bmad_file = bmad_file.c_str();
  bool _err{};
  int output_form_lvalue;
  auto* _output_form{&output_form_lvalue};
  if (output_form.has_value()) {
    output_form_lvalue = output_form.value();
  } else {
    _output_form = nullptr;
  }
  auto* _orbit0 = orbit0.has_value() ? orbit0->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  fortran_write_bmad_lattice_file(
      /* const char* */ _bmad_file,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err,
      /* int* */ _output_form,
      /* void* */ _orbit0);
  return _err;
}
Bmad::WriteGptFieldGridFile1d Bmad::write_gpt_field_grid_file_1d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dz) {
  double _maxfield{};
  double _ref_time{};
  double dz_lvalue;
  auto* _dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  bool _err{};
  fortran_write_gpt_field_grid_file_1d(
      /* int& */ gpt_file_unit,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* double& */ _ref_time,
      /* double* */ _dz,
      /* bool& */ _err);
  return WriteGptFieldGridFile1d{_maxfield, _ref_time, _err};
}
Bmad::WriteGptFieldGridFile2d Bmad::write_gpt_field_grid_file_2d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dr,
    std::optional<double> dz,
    std::optional<double> r_max) {
  double _maxfield{};
  double _ref_time{};
  double dr_lvalue;
  auto* _dr{&dr_lvalue};
  if (dr.has_value()) {
    dr_lvalue = dr.value();
  } else {
    _dr = nullptr;
  }
  double dz_lvalue;
  auto* _dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  double r_max_lvalue;
  auto* _r_max{&r_max_lvalue};
  if (r_max.has_value()) {
    r_max_lvalue = r_max.value();
  } else {
    _r_max = nullptr;
  }
  bool _err{};
  fortran_write_gpt_field_grid_file_2d(
      /* int& */ gpt_file_unit,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* double& */ _ref_time,
      /* double* */ _dr,
      /* double* */ _dz,
      /* double* */ _r_max,
      /* bool& */ _err);
  return WriteGptFieldGridFile2d{_maxfield, _ref_time, _err};
}
Bmad::WriteGptFieldGridFile3d Bmad::write_gpt_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz) {
  auto _base_filename = base_filename.c_str();
  double _maxfield{};
  double _ref_time{};
  double dz_lvalue;
  auto* _dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  bool _err{};
  fortran_write_gpt_field_grid_file_3d(
      /* const char* */ _base_filename,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* double& */ _ref_time,
      /* double* */ _dz,
      /* bool& */ _err);
  return WriteGptFieldGridFile3d{_maxfield, _ref_time, _err};
}
void Bmad::write_lat_line(
    std::string& line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split,
    std::optional<bool> scibmad) {
  auto _line = line.c_str(); // ptr, inout, required
  bool do_split_lvalue;
  auto* _do_split{&do_split_lvalue};
  if (do_split.has_value()) {
    do_split_lvalue = do_split.value();
  } else {
    _do_split = nullptr;
  }
  bool scibmad_lvalue;
  auto* _scibmad{&scibmad_lvalue};
  if (scibmad.has_value()) {
    scibmad_lvalue = scibmad.value();
  } else {
    _scibmad = nullptr;
  }
  fortran_write_lat_line(
      /* const char* */ _line,
      /* int& */ iu,
      /* bool& */ end_is_neigh,
      /* bool* */ _do_split,
      /* bool* */ _scibmad);
}
bool Bmad::write_lattice_in_elegant_format(
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch) {
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto* _ref_orbit = ref_orbit.has_value() ? ref_orbit->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto* _use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto* _include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto* _dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err{};
  fortran_write_lattice_in_elegant_format(
      /* const char* */ _out_file_name,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _ref_orbit,
      /* bool* */ _use_matrix_model,
      /* bool* */ _include_apertures,
      /* double* */ _dr12_drift_max,
      /* int* */ _ix_branch,
      /* bool& */ _err);
  return _err;
}
bool Bmad::write_lattice_in_foreign_format(
    std::string out_type,
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch) {
  auto _out_type = out_type.c_str();
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto* _ref_orbit = ref_orbit.has_value() ? ref_orbit->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto* _use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto* _include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto* _dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err{};
  fortran_write_lattice_in_foreign_format(
      /* const char* */ _out_type,
      /* const char* */ _out_file_name,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _ref_orbit,
      /* bool* */ _use_matrix_model,
      /* bool* */ _include_apertures,
      /* double* */ _dr12_drift_max,
      /* int* */ _ix_branch,
      /* bool& */ _err);
  return _err;
}
bool Bmad::write_lattice_in_mad_format(
    std::string out_type,
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch) {
  auto _out_type = out_type.c_str();
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto* _ref_orbit = ref_orbit.has_value() ? ref_orbit->get().get_fortran_ptr()
                                           : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto* _use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto* _include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto* _dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto* _ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err{};
  fortran_write_lattice_in_mad_format(
      /* const char* */ _out_type,
      /* const char* */ _out_file_name,
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ _ref_orbit,
      /* bool* */ _use_matrix_model,
      /* bool* */ _include_apertures,
      /* double* */ _dr12_drift_max,
      /* int* */ _ix_branch,
      /* bool& */ _err);
  return _err;
}
void Bmad::write_lattice_in_sad_format(
    std::string& out_file_name,
    LatProxy& lat,
    optional_ref<bool> include_apertures,
    optional_ref<int> ix_branch,
    optional_ref<bool> err) {
  auto _out_file_name = out_file_name.c_str(); // ptr, inout, required
  auto* _include_apertures = include_apertures.has_value()
      ? &include_apertures->get()
      : nullptr; // inout, optional
  auto* _ix_branch =
      ix_branch.has_value() ? &ix_branch->get() : nullptr; // inout, optional
  auto* _err = err.has_value() ? &err->get() : nullptr; // inout, optional
  fortran_write_lattice_in_sad_format(
      /* const char* */ _out_file_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _include_apertures,
      /* int* */ _ix_branch,
      /* bool* */ _err);
}
Bmad::WriteLatticeInScibmad Bmad::write_lattice_in_scibmad(LatProxy& lat) {
  char _scibmad_file[4096];
  bool _err_flag{};
  fortran_write_lattice_in_scibmad(
      /* const char* */ _scibmad_file,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err_flag);
  return WriteLatticeInScibmad{_scibmad_file, _err_flag};
}
void Bmad::write_line_element(
    std::string& line,
    int& iu,
    EleProxy& ele,
    LatProxy& lat) {
  auto _line = line.c_str(); // ptr, inout, required
  fortran_write_line_element(
      /* const char* */ _line,
      /* int& */ iu,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr());
}
Bmad::WriteOpalFieldGridFile Bmad::write_opal_field_grid_file(
    int opal_file_unit,
    EleProxy& ele,
    LatParamProxy& param) {
  double _maxfield{};
  bool _err{};
  fortran_write_opal_field_grid_file(
      /* int& */ opal_file_unit,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* bool& */ _err);
  return WriteOpalFieldGridFile{_maxfield, _err};
}
bool Bmad::write_opal_lattice_file(int opal_file_unit, LatProxy& lat) {
  bool _err{};
  fortran_write_opal_lattice_file(
      /* int& */ opal_file_unit,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err);
  return _err;
}
bool Bmad::write_time_particle_distribution(
    int time_file_unit,
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<std::string> style,
    optional_ref<BranchProxy> branch,
    std::optional<std::string> format) {
  const char* _style = style.has_value() ? style->c_str() : nullptr;
  auto* _branch = branch.has_value() ? branch->get().get_fortran_ptr()
                                     : nullptr; // input, optional
  const char* _format = format.has_value() ? format->c_str() : nullptr;
  bool _err{};
  fortran_write_time_particle_distribution(
      /* int& */ time_file_unit,
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _style,
      /* void* */ _branch,
      /* const char* */ _format,
      /* bool& */ _err);
  return _err;
}
void Bmad::xlafun(double& x, double& y, double& z, double& res) {
  fortran_xlafun(
      /* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
int Bmad::xraylib_nist_compound(std::string name) {
  auto _name = name.c_str();
  int _indx{};
  fortran_xraylib_nist_compound(/* const char* */ _name, /* int& */ _indx);
  return _indx;
}
void Bmad::ylafun(double& x, double& y, double& z, double& res) {
  fortran_ylafun(
      /* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
Bmad::ZAtSurface Bmad::z_at_surface(
    EleProxy& ele,
    double& x,
    double& y,
    std::optional<bool> extend_grid) {
  bool _err_flag{};
  bool extend_grid_lvalue;
  auto* _extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  FixedArray1D<Real, 2> _dz_dxy;
  double _z{};
  fortran_z_at_surface(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ x,
      /* double& */ y,
      /* bool& */ _err_flag,
      /* bool* */ _extend_grid,
      /* double* */ _dz_dxy.data(),
      /* double& */ _z);
  return ZAtSurface{_err_flag, _dz_dxy, _z};
}
EleProxy Bmad::zero_ele_kicks() {
  EleProxy _ele;
  fortran_zero_ele_kicks(/* void* */ _ele.get_fortran_ptr());
  return std::move(_ele);
}
EleProxy Bmad::zero_ele_offsets() {
  EleProxy _ele;
  fortran_zero_ele_offsets(/* void* */ _ele.get_fortran_ptr());
  return std::move(_ele);
}
void Bmad::zero_lr_wakes_in_lat(LatProxy& lat) {
  fortran_zero_lr_wakes_in_lat(/* void* */ lat.get_fortran_ptr());
}
void Bmad::zlafun(double& x, double& y, double& z, double& res) {
  fortran_zlafun(
      /* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}