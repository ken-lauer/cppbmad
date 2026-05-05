#include "bmad/generated/bmad_routines.hpp"

#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/json.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

using json = nlohmann::json;
Bmad::AbMultipoleKick Bmad::ab_multipole_kick(
    double a,
    double b,
    int n,
    int ref_species,
    int ele_orientation,
    CoordStruct &coord,
    std::optional<int> pole_type,
    std::optional<double> scale
) {
  double _kx{};
  double _ky{};
  // dk: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _dk_desc;
  _dk_desc.rank = 2;
  FixedArray2D<Real, 2, 2> dk;
  double _dk_vec[2 * 2];
  _dk_desc.data_ptr = _dk_vec;
  _dk_desc.dims[0] = 2;
  _dk_desc.dims[1] = 2;
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  double scale_lvalue;
  auto *_scale{&scale_lvalue};
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
      /* Bmad::array_descriptor_t& */ _dk_desc,
      /* int* */ _pole_type,
      /* double* */ _scale
  );
  vec_to_matrix(_dk_vec, dk);
  return AbMultipoleKick{_kx, _ky, dk};
}
void Bmad::ab_multipole_kicks(
    FArray1D<Real> &an,
    FArray1D<Real> &bn,
    int ix_pole_max,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<int> pole_type,
    std::optional<double> scale,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // an: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _an_desc;
  _an_desc.rank = 1;
  _an_desc.data_ptr = an.data();
  _an_desc.dims[0] = an.size();
  // bn: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _bn_desc;
  _bn_desc.rank = 1;
  _bn_desc.data_ptr = bn.data();
  _bn_desc.dims[0] = bn.size();
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  double scale_lvalue;
  auto *_scale{&scale_lvalue};
  if (scale.has_value()) {
    scale_lvalue = scale.value();
  } else {
    _scale = nullptr;
  }
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_ab_multipole_kicks(
      /* Bmad::array_descriptor_t& */ _an_desc,
      /* Bmad::array_descriptor_t& */ _bn_desc,
      /* int& */ ix_pole_max,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* double* */ _scale,
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* bool* */ _make_matrix
  );
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::absolute_photon_position(CoordStruct &e_orb, CoordStruct &photon_orb) {
  fortran_absolute_photon_position(/* void* */ e_orb.get_fortran_ptr(),
                                   /* void* */ photon_orb.get_fortran_ptr());
}
bool Bmad::absolute_time_tracking(EleStruct &ele) {
  bool _is_abs_time{};
  fortran_absolute_time_tracking(/* void* */ ele.get_fortran_ptr(), /* bool& */ _is_abs_time);
  return _is_abs_time;
}
double Bmad::ac_kicker_amp(EleStruct &ele, CoordStruct &orbit, std::optional<double> true_time) {
  double true_time_lvalue;
  auto *_true_time{&true_time_lvalue};
  if (true_time.has_value()) {
    true_time_lvalue = true_time.value();
  } else {
    _true_time = nullptr;
  }
  double _ac_amp{};
  fortran_ac_kicker_amp(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr(),
                        /* double* */ _true_time,
                        /* double& */ _ac_amp);
  return _ac_amp;
}
Bmad::ActionToXyz Bmad::action_to_xyz(LatStruct &ring, int ix, FixedArray1D<Real, 6> J) {
  // J: in NOT (CppWrapperGeneralArgumentArray) (['1:6'])
  Bmad::array_descriptor_t _J_desc;
  _J_desc.rank = 1;
  _J_desc.data_ptr = J.data();
  _J_desc.dims[0] = J.size();
  // X: out NOT (CppWrapperGeneralArgumentArray) (['1:6'])
  Bmad::array_descriptor_t _X_desc;
  _X_desc.rank = 1;
  FixedArray1D<Real, 6> _X;
  _X_desc.data_ptr = _X.data();
  _X_desc.dims[0] = _X.size();
  bool _err_flag{};
  fortran_action_to_xyz(/* void* */ ring.get_fortran_ptr(),
                        /* int& */ ix,
                        /* Bmad::array_descriptor_t& */ _J_desc,
                        /* Bmad::array_descriptor_t& */ _X_desc,
                        /* bool& */ _err_flag);
  return ActionToXyz{_X, _err_flag};
}
void Bmad::add_lattice_control_structs(
    EleStruct &ele,
    std::optional<int> n_add_slave,
    std::optional<int> n_add_lord,
    std::optional<int> n_add_slave_field,
    std::optional<int> n_add_lord_field,
    std::optional<bool> add_at_end
) {
  int n_add_slave_lvalue;
  auto *_n_add_slave{&n_add_slave_lvalue};
  if (n_add_slave.has_value()) {
    n_add_slave_lvalue = n_add_slave.value();
  } else {
    _n_add_slave = nullptr;
  }
  int n_add_lord_lvalue;
  auto *_n_add_lord{&n_add_lord_lvalue};
  if (n_add_lord.has_value()) {
    n_add_lord_lvalue = n_add_lord.value();
  } else {
    _n_add_lord = nullptr;
  }
  int n_add_slave_field_lvalue;
  auto *_n_add_slave_field{&n_add_slave_field_lvalue};
  if (n_add_slave_field.has_value()) {
    n_add_slave_field_lvalue = n_add_slave_field.value();
  } else {
    _n_add_slave_field = nullptr;
  }
  int n_add_lord_field_lvalue;
  auto *_n_add_lord_field{&n_add_lord_field_lvalue};
  if (n_add_lord_field.has_value()) {
    n_add_lord_field_lvalue = n_add_lord_field.value();
  } else {
    _n_add_lord_field = nullptr;
  }
  bool add_at_end_lvalue;
  auto *_add_at_end{&add_at_end_lvalue};
  if (add_at_end.has_value()) {
    add_at_end_lvalue = add_at_end.value();
  } else {
    _add_at_end = nullptr;
  }
  fortran_add_lattice_control_structs(/* void* */ ele.get_fortran_ptr(),
                                      /* int* */ _n_add_slave,
                                      /* int* */ _n_add_lord,
                                      /* int* */ _n_add_slave_field,
                                      /* int* */ _n_add_lord_field,
                                      /* bool* */ _add_at_end);
}
Bmad::AddSuperimpose Bmad::add_superimpose(
    LatStruct &lat,
    EleStruct &super_ele_in,
    int ix_branch,
    std::optional<bool> save_null_drift,
    std::optional<bool> create_jumbo_slave,
    std::optional<int> ix_insert,
    std::optional<bool> mangle_slave_names,
    std::optional<bool> wrap
) {
  bool _err_flag{};
  void *_super_ele_out;
  bool save_null_drift_lvalue;
  auto *_save_null_drift{&save_null_drift_lvalue};
  if (save_null_drift.has_value()) {
    save_null_drift_lvalue = save_null_drift.value();
  } else {
    _save_null_drift = nullptr;
  }
  bool create_jumbo_slave_lvalue;
  auto *_create_jumbo_slave{&create_jumbo_slave_lvalue};
  if (create_jumbo_slave.has_value()) {
    create_jumbo_slave_lvalue = create_jumbo_slave.value();
  } else {
    _create_jumbo_slave = nullptr;
  }
  int ix_insert_lvalue;
  auto *_ix_insert{&ix_insert_lvalue};
  if (ix_insert.has_value()) {
    ix_insert_lvalue = ix_insert.value();
  } else {
    _ix_insert = nullptr;
  }
  bool mangle_slave_names_lvalue;
  auto *_mangle_slave_names{&mangle_slave_names_lvalue};
  if (mangle_slave_names.has_value()) {
    mangle_slave_names_lvalue = mangle_slave_names.value();
  } else {
    _mangle_slave_names = nullptr;
  }
  bool wrap_lvalue;
  auto *_wrap{&wrap_lvalue};
  if (wrap.has_value()) {
    wrap_lvalue = wrap.value();
  } else {
    _wrap = nullptr;
  }
  fortran_add_superimpose(/* void* */ lat.get_fortran_ptr(),
                          /* void* */ super_ele_in.get_fortran_ptr(),
                          /* int& */ ix_branch,
                          /* bool& */ _err_flag,
                          /* void* */ &_super_ele_out,
                          /* bool* */ _save_null_drift,
                          /* bool* */ _create_jumbo_slave,
                          /* int* */ _ix_insert,
                          /* bool* */ _mangle_slave_names,
                          /* bool* */ _wrap);
  return AddSuperimpose{
      _err_flag,
      std::move((_super_ele_out ? std::make_optional<EleStruct>(_super_ele_out) : std::nullopt))
  };
}
void Bmad::add_this_multipass(
    LatStruct &lat,
    LatEleLocStructArray1D m_slaves,
    optional_ref<EleStruct> lord_in
) {
  // m_slaves: LatEleLocStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _m_slaves_desc;
  _m_slaves_desc.rank = 1;
  _m_slaves_desc.data_ptr = m_slaves.data();
  _m_slaves_desc.dims[0] = m_slaves.size();
  _m_slaves_desc.strides[0] = 1;
  auto *_lord_in =
      lord_in.has_value() ? lord_in->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_add_this_multipass(/* void* */ lat.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _m_slaves_desc,
                             /* void* */ _lord_in);
}
void Bmad::add_this_name_to_list(
    EleStruct &ele,
    CharacterAlloc1D &names,
    IntAlloc1D &an_indexx,
    int n_names,
    int ix_match,
    bool has_been_added,
    ElePointerStructAlloc1D named_eles
) {
  // intent=inout character array container
  // intent=inout allocatable general array
  // intent=inout allocatable type array
  fortran_add_this_name_to_list(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ names.get_fortran_ptr(),
                                /* void* */ an_indexx.get_fortran_ptr(),
                                /* int& */ n_names,
                                /* int& */ ix_match,
                                /* bool& */ has_been_added,
                                /* void* */ named_eles.get_fortran_ptr());
}
void Bmad::add_this_taylor_term(EleStruct &ele, int i_out, double coef, FixedArray1D<Int, 6> expn) {
  // expn: inout NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _expn_desc;
  _expn_desc.rank = 1;
  _expn_desc.data_ptr = expn.data();
  _expn_desc.dims[0] = expn.size();
  fortran_add_this_taylor_term(/* void* */ ele.get_fortran_ptr(),
                               /* int& */ i_out,
                               /* double& */ coef,
                               /* Bmad::array_descriptor_t& */ _expn_desc);
}
void Bmad::adjust_super_slave_names(
    LatStruct &lat,
    int ix1_lord,
    int ix2_lord,
    std::optional<bool> first_time
) {
  bool first_time_lvalue;
  auto *_first_time{&first_time_lvalue};
  if (first_time.has_value()) {
    first_time_lvalue = first_time.value();
  } else {
    _first_time = nullptr;
  }
  fortran_adjust_super_slave_names(/* void* */ lat.get_fortran_ptr(),
                                   /* int& */ ix1_lord,
                                   /* int& */ ix2_lord,
                                   /* bool* */ _first_time);
}
void Bmad::allocate_branch_array(LatStruct &lat, int upper_bound) {
  fortran_allocate_branch_array(/* void* */ lat.get_fortran_ptr(), /* int& */ upper_bound);
}
void Bmad::allocate_grid_field(GridFieldStructArray1D g_field, int n_gf) {
  // g_field: GridFieldStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _g_field_desc;
  _g_field_desc.rank = 1;
  _g_field_desc.data_ptr = g_field.data();
  _g_field_desc.dims[0] = g_field.size();
  _g_field_desc.strides[0] = 1;
  fortran_allocate_grid_field(/* Bmad::array_descriptor_t& */ _g_field_desc, /* int& */ n_gf);
}
void Bmad::allocate_lat_ele_array(
    LatStruct &lat,
    std::optional<int> upper_bound,
    std::optional<int> ix_branch,
    std::optional<bool> do_ramper_slave_setup
) {
  int upper_bound_lvalue;
  auto *_upper_bound{&upper_bound_lvalue};
  if (upper_bound.has_value()) {
    upper_bound_lvalue = upper_bound.value();
  } else {
    _upper_bound = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool do_ramper_slave_setup_lvalue;
  auto *_do_ramper_slave_setup{&do_ramper_slave_setup_lvalue};
  if (do_ramper_slave_setup.has_value()) {
    do_ramper_slave_setup_lvalue = do_ramper_slave_setup.value();
  } else {
    _do_ramper_slave_setup = nullptr;
  }
  fortran_allocate_lat_ele_array(/* void* */ lat.get_fortran_ptr(),
                                 /* int* */ _upper_bound,
                                 /* int* */ _ix_branch,
                                 /* bool* */ _do_ramper_slave_setup);
}
double Bmad::angle_between_polars(SpinPolarStruct &polar1, SpinPolarStruct &polar2) {
  double _angle{};
  fortran_angle_between_polars(/* void* */ polar1.get_fortran_ptr(),
                               /* void* */ polar2.get_fortran_ptr(),
                               /* double& */ _angle);
  return _angle;
}
void Bmad::angle_to_canonical_coords(CoordStruct &orbit, std::optional<std::string> coord_type) {
  const char *_coord_type = coord_type.has_value() ? coord_type->c_str() : nullptr;
  fortran_angle_to_canonical_coords(/* void* */ orbit.get_fortran_ptr(),
                                    /* const char* */ _coord_type);
}
void Bmad::aperture_bookkeeper(EleStruct &ele) {
  fortran_aperture_bookkeeper(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::apply_all_rampers(LatStruct &lat) {
  bool _err_flag{};
  fortran_apply_all_rampers(/* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::apply_energy_kick(
    double dE,
    CoordStruct &orbit,
    FixedArray1D<Real, 2> ddE_dr,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // ddE_dr: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _ddE_dr_desc;
  _ddE_dr_desc.rank = 1;
  _ddE_dr_desc.data_ptr = ddE_dr.data();
  _ddE_dr_desc.dims[0] = ddE_dr.size();
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_apply_energy_kick(
      /* double& */ dE,
      /* void* */ orbit.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _ddE_dr_desc,
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* bool* */ _make_matrix
  );
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::apply_patch_to_ptc_fibre(EleStruct &ele) {
  fortran_apply_patch_to_ptc_fibre(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::apply_rampers_to_slave(EleStruct &slave) {
  bool _err_flag{};
  fortran_apply_rampers_to_slave(/* void* */ slave.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::array_re_str(
    FArray1D<Real> &arr,
    std::string str_out,
    std::optional<std::string> parens_in
) {
  // arr: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_desc;
  _arr_desc.rank = 1;
  _arr_desc.data_ptr = arr.data();
  _arr_desc.dims[0] = arr.size();
  const char *_parens_in = parens_in.has_value() ? parens_in->c_str() : nullptr;
  auto _str_out = str_out.c_str();
  fortran_array_re_str(
      /* Bmad::array_descriptor_t& */ _arr_desc,
      /* const char* */ _parens_in,
      /* const char* */ _str_out
  );
}
void Bmad::astra_max_field_reference(GridFieldPt1Struct &pt0, EleStruct &ele, double field_value) {
  fortran_astra_max_field_reference(/* void* */ pt0.get_fortran_ptr(),
                                    /* void* */ ele.get_fortran_ptr(),
                                    /* double& */ field_value);
}
bool Bmad::at_this_ele_end(int now_at, int where_at) {
  bool _is_at_this_end{};
  fortran_at_this_ele_end(/* int& */ now_at, /* int& */ where_at, /* bool& */ _is_at_this_end);
  return _is_at_this_end;
}
void Bmad::attribute_bookkeeper(EleStruct &ele, std::optional<bool> force_bookkeeping) {
  bool force_bookkeeping_lvalue;
  auto *_force_bookkeeping{&force_bookkeeping_lvalue};
  if (force_bookkeeping.has_value()) {
    force_bookkeeping_lvalue = force_bookkeeping.value();
  } else {
    _force_bookkeeping = nullptr;
  }
  fortran_attribute_bookkeeper(/* void* */ ele.get_fortran_ptr(), /* bool* */ _force_bookkeeping);
}
Bmad::AttributeFree1 Bmad::attribute_free(
    int ix_ele,
    std::string attrib_name,
    LatStruct &lat,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free
) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto *_except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto *_dependent_attribs_free{&dependent_attribs_free_lvalue};
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
      /* bool& */ _free
  );
  return AttributeFree1{_why_not_free, _free};
}
Bmad::AttributeFree2 Bmad::attribute_free(
    EleStruct &ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free
) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto *_except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto *_dependent_attribs_free{&dependent_attribs_free_lvalue};
  if (dependent_attribs_free.has_value()) {
    dependent_attribs_free_lvalue = dependent_attribs_free.value();
  } else {
    _dependent_attribs_free = nullptr;
  }
  int _why_not_free{};
  bool _free{};
  fortran_attribute_free2(/* void* */ ele.get_fortran_ptr(),
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
    LatStruct &lat,
    std::optional<bool> err_print_flag,
    std::optional<bool> except_overlay,
    std::optional<bool> dependent_attribs_free
) {
  auto _attrib_name = attrib_name.c_str();
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool except_overlay_lvalue;
  auto *_except_overlay{&except_overlay_lvalue};
  if (except_overlay.has_value()) {
    except_overlay_lvalue = except_overlay.value();
  } else {
    _except_overlay = nullptr;
  }
  bool dependent_attribs_free_lvalue;
  auto *_dependent_attribs_free{&dependent_attribs_free_lvalue};
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
      /* bool& */ _free
  );
  return AttributeFree3{_why_not_free, _free};
}
Bmad::AttributeIndex1 Bmad::attribute_index(
    EleStruct &ele,
    std::string name,
    std::optional<bool> can_abbreviate,
    std::optional<bool> print_error
) {
  auto _name = name.c_str();
  char _full_name[4096];
  bool can_abbreviate_lvalue;
  auto *_can_abbreviate{&can_abbreviate_lvalue};
  if (can_abbreviate.has_value()) {
    can_abbreviate_lvalue = can_abbreviate.value();
  } else {
    _can_abbreviate = nullptr;
  }
  bool print_error_lvalue;
  auto *_print_error{&print_error_lvalue};
  if (print_error.has_value()) {
    print_error_lvalue = print_error.value();
  } else {
    _print_error = nullptr;
  }
  int _attrib_index{};
  fortran_attribute_index1(/* void* */ ele.get_fortran_ptr(),
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
    std::optional<bool> print_error
) {
  auto _name = name.c_str();
  char _full_name[4096];
  bool can_abbreviate_lvalue;
  auto *_can_abbreviate{&can_abbreviate_lvalue};
  if (can_abbreviate.has_value()) {
    can_abbreviate_lvalue = can_abbreviate.value();
  } else {
    _can_abbreviate = nullptr;
  }
  bool print_error_lvalue;
  auto *_print_error{&print_error_lvalue};
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
      /* int& */ _attrib_index
  );
  return AttributeIndex2{_full_name, _attrib_index};
}
std::string Bmad::attribute_name(int key, int ix_att, std::optional<bool> show_private) {
  bool show_private_lvalue;
  auto *_show_private{&show_private_lvalue};
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
      /* const char* */ _attrib_name
  );
  return _attrib_name;
}
std::string Bmad::attribute_name(EleStruct &ele, int ix_att, std::optional<bool> show_private) {
  bool show_private_lvalue;
  auto *_show_private{&show_private_lvalue};
  if (show_private.has_value()) {
    show_private_lvalue = show_private.value();
  } else {
    _show_private = nullptr;
  }
  char _attrib_name[4096];
  fortran_attribute_name2(/* void* */ ele.get_fortran_ptr(),
                          /* int& */ ix_att,
                          /* bool* */ _show_private,
                          /* const char* */ _attrib_name);
  return _attrib_name;
}
int Bmad::attribute_type(std::string attrib_name, optional_ref<EleStruct> ele) {
  auto _attrib_name = attrib_name.c_str();
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  int _attrib_type{};
  fortran_attribute_type(/* const char* */ _attrib_name, /* void* */ _ele, /* int& */ _attrib_type);
  return _attrib_type;
}
std::string
Bmad::attribute_units(std::string attrib_name, std::optional<std::string> unrecognized_units) {
  auto _attrib_name = attrib_name.c_str();
  const char *_unrecognized_units =
      unrecognized_units.has_value() ? unrecognized_units->c_str() : nullptr;
  char _attrib_units[4096];
  fortran_attribute_units(
      /* const char* */ _attrib_name,
      /* const char* */ _unrecognized_units,
      /* const char* */ _attrib_units
  );
  return _attrib_units;
}
bool Bmad::autoscale_phase_and_amp(
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> scale_phase,
    std::optional<bool> scale_amp,
    std::optional<bool> call_bookkeeper
) {
  bool _err_flag{};
  bool scale_phase_lvalue;
  auto *_scale_phase{&scale_phase_lvalue};
  if (scale_phase.has_value()) {
    scale_phase_lvalue = scale_phase.value();
  } else {
    _scale_phase = nullptr;
  }
  bool scale_amp_lvalue;
  auto *_scale_amp{&scale_amp_lvalue};
  if (scale_amp.has_value()) {
    scale_amp_lvalue = scale_amp.value();
  } else {
    _scale_amp = nullptr;
  }
  bool call_bookkeeper_lvalue;
  auto *_call_bookkeeper{&call_bookkeeper_lvalue};
  if (call_bookkeeper.has_value()) {
    call_bookkeeper_lvalue = call_bookkeeper.value();
  } else {
    _call_bookkeeper = nullptr;
  }
  fortran_autoscale_phase_and_amp(/* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* bool& */ _err_flag,
                                  /* bool* */ _scale_phase,
                                  /* bool* */ _scale_amp,
                                  /* bool* */ _call_bookkeeper);
  return _err_flag;
}
TwissStruct Bmad::average_twiss(double frac1, TwissStruct &twiss1, TwissStruct &twiss2) {
  TwissStruct _ave_twiss;
  fortran_average_twiss(
      /* double& */ frac1,
      /* void* */ twiss1.get_fortran_ptr(),
      /* void* */ twiss2.get_fortran_ptr(),
      /* void* */ _ave_twiss.get_fortran_ptr()
  );
  return std::move(_ave_twiss);
}
Bmad::BbiKick
Bmad::bbi_kick(double x, double y, FixedArray1D<Real, 2> sigma, std::optional<bool> linear_kick) {
  // sigma: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _sigma_desc;
  _sigma_desc.rank = 1;
  _sigma_desc.data_ptr = sigma.data();
  _sigma_desc.dims[0] = sigma.size();
  // nk: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _nk_desc;
  _nk_desc.rank = 1;
  FixedArray1D<Real, 2> _nk;
  _nk_desc.data_ptr = _nk.data();
  _nk_desc.dims[0] = _nk.size();
  // dnk: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _dnk_desc;
  _dnk_desc.rank = 2;
  FixedArray2D<Real, 2, 2> dnk;
  double _dnk_vec[2 * 2];
  _dnk_desc.data_ptr = _dnk_vec;
  _dnk_desc.dims[0] = 2;
  _dnk_desc.dims[1] = 2;
  bool linear_kick_lvalue;
  auto *_linear_kick{&linear_kick_lvalue};
  if (linear_kick.has_value()) {
    linear_kick_lvalue = linear_kick.value();
  } else {
    _linear_kick = nullptr;
  }
  fortran_bbi_kick(
      /* double& */ x,
      /* double& */ y,
      /* Bmad::array_descriptor_t& */ _sigma_desc,
      /* Bmad::array_descriptor_t& */ _nk_desc,
      /* Bmad::array_descriptor_t& */ _dnk_desc,
      /* bool* */ _linear_kick
  );
  vec_to_matrix(_dnk_vec, dnk);
  return BbiKick{_nk, dnk};
}
void Bmad::bbi_slice_calc(EleStruct &ele, int n_slice, FArray1D<Real> &z_slice) {
  // z_slice: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _z_slice_desc;
  _z_slice_desc.rank = 1;
  _z_slice_desc.data_ptr = z_slice.data();
  _z_slice_desc.dims[0] = z_slice.size();
  fortran_bbi_slice_calc(/* void* */ ele.get_fortran_ptr(),
                         /* int& */ n_slice,
                         /* Bmad::array_descriptor_t& */ _z_slice_desc);
}
FixedArray2D<Real, 6, 6> Bmad::beam_envelope_ibs(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species
) {
  // sigma_mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  // ibs_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _ibs_mat_desc;
  _ibs_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> ibs_mat;
  double _ibs_mat_vec[6 * 6];
  _ibs_mat_desc.data_ptr = _ibs_mat_vec;
  _ibs_mat_desc.dims[0] = 6;
  _ibs_mat_desc.dims[1] = 6;
  fortran_beam_envelope_ibs(
      /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
      /* Bmad::array_descriptor_t& */ _ibs_mat_desc,
      /* bool& */ tail_cut,
      /* double& */ tau,
      /* double& */ energy,
      /* double& */ n_part,
      /* int& */ species
  );
  vec_to_matrix(_ibs_mat_vec, ibs_mat);
  return ibs_mat;
}
void Bmad::beam_equal_beam(BeamStruct &beam1, BeamStruct &beam2) {
  fortran_beam_equal_beam(/* void* */ beam1.get_fortran_ptr(), /* void* */ beam2.get_fortran_ptr());
}
Bmad::BeamInitSetup Bmad::beam_init_setup(
    BeamInitStruct &beam_init_in,
    EleStruct &ele,
    int species,
    optional_ref<NormalModesStruct> modes
) {
  auto *_modes = modes.has_value() ? modes->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err_flag{};
  BeamInitStruct _beam_init_set;
  fortran_beam_init_setup(/* void* */ beam_init_in.get_fortran_ptr(),
                          /* void* */ ele.get_fortran_ptr(),
                          /* int& */ species,
                          /* void* */ _modes,
                          /* bool& */ _err_flag,
                          /* void* */ _beam_init_set.get_fortran_ptr());
  return BeamInitSetup{_err_flag, std::move(_beam_init_set)};
}
Bmad::BeamTilts Bmad::beam_tilts(FixedArray2D<Real, 6, 6> S) {
  // S: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _S_desc;
  _S_desc.rank = 2;
  double _S_vec[6 * 6];
  _S_desc.data_ptr = _S_vec;
  _S_desc.dims[0] = 6;
  _S_desc.dims[1] = 6;
  matrix_to_vec(S, _S_vec);
  double _angle_xy{};
  double _angle_xz{};
  double _angle_yz{};
  double _angle_xpz{};
  double _angle_ypz{};
  fortran_beam_tilts(
      /* Bmad::array_descriptor_t& */ _S_desc,
      /* double& */ _angle_xy,
      /* double& */ _angle_xz,
      /* double& */ _angle_yz,
      /* double& */ _angle_xpz,
      /* double& */ _angle_ypz
  );
  return BeamTilts{_angle_xy, _angle_xz, _angle_yz, _angle_xpz, _angle_ypz};
}
Fibre Bmad::beambeam_fibre_setup(EleStruct &ele) {
  Fibre _ptc_fibre;
  fortran_beambeam_fibre_setup(/* void* */ ele.get_fortran_ptr(),
                               /* void* */ _ptc_fibre.get_fortran_ptr());
  return std::move(_ptc_fibre);
}
void Bmad::bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<bool> track_spin
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  bool track_spin_lvalue;
  auto *_track_spin{&track_spin_lvalue};
  if (track_spin.has_value()) {
    track_spin_lvalue = track_spin.value();
  } else {
    _track_spin = nullptr;
  }
  fortran_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ param.get_fortran_ptr(),
                         /* int& */ particle_at,
                         /* void* */ orb.get_fortran_ptr(),
                         /* Bmad::array_descriptor_t& */ _mat6_desc,
                         /* bool* */ _make_matrix,
                         /* bool* */ _track_spin);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
EmFieldStruct Bmad::bend_exact_multipole_field(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential
) {
  EmFieldStruct _field;
  bool calc_dfield_lvalue;
  auto *_calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto *_calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_bend_exact_multipole_field(/* void* */ ele.get_fortran_ptr(),
                                     /* void* */ param.get_fortran_ptr(),
                                     /* void* */ orbit.get_fortran_ptr(),
                                     /* bool& */ local_ref_frame,
                                     /* void* */ _field.get_fortran_ptr(),
                                     /* bool* */ _calc_dfield,
                                     /* bool* */ _calc_potential);
  return std::move(_field);
}
bool Bmad::bend_length_has_been_set(EleStruct &ele) {
  bool _is_set{};
  fortran_bend_length_has_been_set(/* void* */ ele.get_fortran_ptr(), /* bool& */ _is_set);
  return _is_set;
}
double Bmad::bend_photon_e_rel_init(std::optional<double> r_in) {
  double r_in_lvalue;
  auto *_r_in{&r_in_lvalue};
  if (r_in.has_value()) {
    r_in_lvalue = r_in.value();
  } else {
    _r_in = nullptr;
  }
  double _E_rel{};
  fortran_bend_photon_e_rel_init(/* double* */ _r_in, /* double& */ _E_rel);
  return _E_rel;
}
double Bmad::bend_photon_energy_integ_prob(double E_photon, double g_bend, double gamma) {
  double _integ_prob{};
  fortran_bend_photon_energy_integ_prob(
      /* double& */ E_photon,
      /* double& */ g_bend,
      /* double& */ gamma,
      /* double& */ _integ_prob
  );
  return _integ_prob;
}
double Bmad::bend_photon_energy_normalized_probability(double E_rel) {
  double _prob{};
  fortran_bend_photon_energy_normalized_probability(/* double& */ E_rel, /* double& */ _prob);
  return _prob;
}
CoordStruct Bmad::bend_photon_init(
    double g_bend_x,
    double g_bend_y,
    double gamma,
    std::optional<double> E_min,
    std::optional<double> E_max,
    std::optional<double> E_integ_prob,
    std::optional<double> vert_angle_min,
    std::optional<double> vert_angle_max,
    std::optional<bool> vert_angle_symmetric,
    std::optional<double> emit_probability
) {
  CoordStruct _orbit;
  double E_min_lvalue;
  auto *_E_min{&E_min_lvalue};
  if (E_min.has_value()) {
    E_min_lvalue = E_min.value();
  } else {
    _E_min = nullptr;
  }
  double E_max_lvalue;
  auto *_E_max{&E_max_lvalue};
  if (E_max.has_value()) {
    E_max_lvalue = E_max.value();
  } else {
    _E_max = nullptr;
  }
  double E_integ_prob_lvalue;
  auto *_E_integ_prob{&E_integ_prob_lvalue};
  if (E_integ_prob.has_value()) {
    E_integ_prob_lvalue = E_integ_prob.value();
  } else {
    _E_integ_prob = nullptr;
  }
  double vert_angle_min_lvalue;
  auto *_vert_angle_min{&vert_angle_min_lvalue};
  if (vert_angle_min.has_value()) {
    vert_angle_min_lvalue = vert_angle_min.value();
  } else {
    _vert_angle_min = nullptr;
  }
  double vert_angle_max_lvalue;
  auto *_vert_angle_max{&vert_angle_max_lvalue};
  if (vert_angle_max.has_value()) {
    vert_angle_max_lvalue = vert_angle_max.value();
  } else {
    _vert_angle_max = nullptr;
  }
  bool vert_angle_symmetric_lvalue;
  auto *_vert_angle_symmetric{&vert_angle_symmetric_lvalue};
  if (vert_angle_symmetric.has_value()) {
    vert_angle_symmetric_lvalue = vert_angle_symmetric.value();
  } else {
    _vert_angle_symmetric = nullptr;
  }
  double emit_probability_lvalue;
  auto *_emit_probability{&emit_probability_lvalue};
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
      /* double* */ _emit_probability
  );
  return std::move(_orbit);
}
CoordStruct Bmad::bend_photon_polarization_init(
    double g_bend_x,
    double g_bend_y,
    double E_rel,
    double gamma_phi
) {
  CoordStruct _orbit;
  fortran_bend_photon_polarization_init(
      /* double& */ g_bend_x,
      /* double& */ g_bend_y,
      /* double& */ E_rel,
      /* double& */ gamma_phi,
      /* void* */ _orbit.get_fortran_ptr()
  );
  return std::move(_orbit);
}
double Bmad::bend_photon_vert_angle_init(
    double E_rel,
    double gamma,
    std::optional<double> r_in,
    std::optional<bool> invert
) {
  double r_in_lvalue;
  auto *_r_in{&r_in_lvalue};
  if (r_in.has_value()) {
    r_in_lvalue = r_in.value();
  } else {
    _r_in = nullptr;
  }
  bool invert_lvalue;
  auto *_invert{&invert_lvalue};
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
      /* double& */ _phi
  );
  return _phi;
}
Bmad::BendShift Bmad::bend_shift(
    FloorPositionStruct &position1,
    double g,
    double delta_s,
    std::optional<double> ref_tilt
) {
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  double ref_tilt_lvalue;
  auto *_ref_tilt{&ref_tilt_lvalue};
  if (ref_tilt.has_value()) {
    ref_tilt_lvalue = ref_tilt.value();
  } else {
    _ref_tilt = nullptr;
  }
  FloorPositionStruct _position2;
  fortran_bend_shift(/* void* */ position1.get_fortran_ptr(),
                     /* double& */ g,
                     /* double& */ delta_s,
                     /* Bmad::array_descriptor_t& */ _w_mat_desc,
                     /* double* */ _ref_tilt,
                     /* void* */ _position2.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return BendShift{w_mat, std::move(_position2)};
}
double Bmad::bend_vert_angle_integ_prob(double vert_angle, double E_rel, double gamma) {
  double _integ_prob{};
  fortran_bend_vert_angle_integ_prob(
      /* double& */ vert_angle,
      /* double& */ E_rel,
      /* double& */ gamma,
      /* double& */ _integ_prob
  );
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
    double L
) {
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
      /* double& */ _sigma_z
  );
  return _sigma_z;
}
Bmad::BmadParser Bmad::bmad_parser(
    std::string lat_file,
    std::optional<bool> make_mats6,
    std::optional<std::string> use_line
) {
  auto _lat_file = lat_file.c_str();
  LatStruct _lat;
  bool make_mats6_lvalue;
  auto *_make_mats6{&make_mats6_lvalue};
  if (make_mats6.has_value()) {
    make_mats6_lvalue = make_mats6.value();
  } else {
    _make_mats6 = nullptr;
  }
  bool _digested_read_ok{};
  const char *_use_line = use_line.has_value() ? use_line->c_str() : nullptr;
  bool _err_flag{};
  LatStruct _parse_lat;
  fortran_bmad_parser(
      /* const char* */ _lat_file,
      /* void* */ _lat.get_fortran_ptr(),
      /* bool* */ _make_mats6,
      /* bool& */ _digested_read_ok,
      /* const char* */ _use_line,
      /* bool& */ _err_flag,
      /* void* */ _parse_lat.get_fortran_ptr()
  );
  return BmadParser{std::move(_lat), _digested_read_ok, _err_flag, std::move(_parse_lat)};
}
void Bmad::bmad_parser2(
    std::string lat_file,
    LatStruct &lat,
    std::optional<CoordStructArray1D> orbit,
    std::optional<bool> make_mats6,
    std::optional<bool> err_flag,
    optional_ref<LatStruct> parse_lat
) {
  auto _lat_file = lat_file.c_str();
  // orbit: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  if (orbit) {
    _orbit_desc.data_ptr = orbit->data();
    _orbit_desc.dims[0] = orbit->size();
  } else {
    _orbit_desc.data_ptr = nullptr;
    _orbit_desc.dims[0] = 0;
  }
  _orbit_desc.strides[0] = 1;
  bool make_mats6_lvalue;
  auto *_make_mats6{&make_mats6_lvalue};
  if (make_mats6.has_value()) {
    make_mats6_lvalue = make_mats6.value();
  } else {
    _make_mats6 = nullptr;
  }
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  auto *_parse_lat =
      parse_lat.has_value() ? parse_lat->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_bmad_parser2(
      /* const char* */ _lat_file,
      /* void* */ lat.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _orbit_desc,
      /* bool* */ _make_mats6,
      /* bool* */ _err_flag,
      /* void* */ _parse_lat
  );
}
void Bmad::bmad_patch_parameters_to_ptc(FixedArray1D<Real, 3> ang, FixedArray2D<Real, 3, 3> exi) {
  // ang: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _ang_desc;
  _ang_desc.rank = 1;
  _ang_desc.data_ptr = ang.data();
  _ang_desc.dims[0] = ang.size();
  // exi: inout NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _exi_desc;
  _exi_desc.rank = 2;
  double _exi_vec[3 * 3];
  _exi_desc.data_ptr = _exi_vec;
  _exi_desc.dims[0] = 3;
  _exi_desc.dims[1] = 3;
  matrix_to_vec(exi, _exi_vec);
  fortran_bmad_patch_parameters_to_ptc(
      /* Bmad::array_descriptor_t& */ _ang_desc,
      /* Bmad::array_descriptor_t& */ _exi_desc
  );
  vec_to_matrix(_exi_vec, exi);
}
void Bmad::bp_set_ran_status() { fortran_bp_set_ran_status(); }
void Bmad::branch_equal_branch(BranchStruct &branch1, BranchStruct &branch2) {
  fortran_branch_equal_branch(/* void* */ branch1.get_fortran_ptr(),
                              /* void* */ branch2.get_fortran_ptr());
}
std::string Bmad::branch_name(BranchStruct &branch) {
  char _name[4096];
  fortran_branch_name(/* void* */ branch.get_fortran_ptr(), /* const char* */ _name);
  return _name;
}
void Bmad::branch_to_ptc_m_u(BranchStruct &branch) {
  fortran_branch_to_ptc_m_u(/* void* */ branch.get_fortran_ptr());
}
void Bmad::bunch_equal_bunch(BunchStruct &bunch1, BunchStruct &bunch2) {
  fortran_bunch_equal_bunch(/* void* */ bunch1.get_fortran_ptr(),
                            /* void* */ bunch2.get_fortran_ptr());
}
FixedArray2D<Real, 2, 2> Bmad::c_to_cbar(EleStruct &ele) {
  // cbar_mat: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _cbar_mat_desc;
  _cbar_mat_desc.rank = 2;
  FixedArray2D<Real, 2, 2> cbar_mat;
  double _cbar_mat_vec[2 * 2];
  _cbar_mat_desc.data_ptr = _cbar_mat_vec;
  _cbar_mat_desc.dims[0] = 2;
  _cbar_mat_desc.dims[1] = 2;
  fortran_c_to_cbar(/* void* */ ele.get_fortran_ptr(),
                    /* Bmad::array_descriptor_t& */ _cbar_mat_desc);
  vec_to_matrix(_cbar_mat_vec, cbar_mat);
  return cbar_mat;
}
Bmad::CalcBunchParams Bmad::calc_bunch_params(
    BunchStruct &bunch,
    std::optional<bool> print_err,
    std::optional<bool> is_time_coords,
    optional_ref<EleStruct> ele
) {
  BunchParamsStruct _bunch_params;
  bool _error{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  // n_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _n_mat_desc;
  _n_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> n_mat;
  double _n_mat_vec[6 * 6];
  _n_mat_desc.data_ptr = _n_mat_vec;
  _n_mat_desc.dims[0] = 6;
  _n_mat_desc.dims[1] = 6;
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_calc_bunch_params(/* void* */ bunch.get_fortran_ptr(),
                            /* void* */ _bunch_params.get_fortran_ptr(),
                            /* bool& */ _error,
                            /* bool* */ _print_err,
                            /* Bmad::array_descriptor_t& */ _n_mat_desc,
                            /* bool* */ _is_time_coords,
                            /* void* */ _ele);
  vec_to_matrix(_n_mat_vec, n_mat);
  return CalcBunchParams{std::move(_bunch_params), _error, n_mat};
}
bool Bmad::calc_bunch_params_slice(
    BunchStruct &bunch,
    BunchParamsStruct &bunch_params,
    int plane,
    double slice_center,
    double slice_spread,
    std::optional<bool> print_err,
    std::optional<bool> is_time_coords,
    optional_ref<EleStruct> ele
) {
  bool _err{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_calc_bunch_params_slice(/* void* */ bunch.get_fortran_ptr(),
                                  /* void* */ bunch_params.get_fortran_ptr(),
                                  /* int& */ plane,
                                  /* double& */ slice_center,
                                  /* double& */ slice_spread,
                                  /* bool& */ _err,
                                  /* bool* */ _print_err,
                                  /* bool* */ _is_time_coords,
                                  /* void* */ _ele);
  return _err;
}
bool Bmad::calc_bunch_params_z_slice(
    BunchStruct &bunch,
    BunchParamsStruct &bunch_params,
    FixedArray1D<Real, 2> slice_bounds,
    std::optional<bool> print_err,
    std::optional<bool> is_time_coords,
    optional_ref<EleStruct> ele
) {
  // slice_bounds: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _slice_bounds_desc;
  _slice_bounds_desc.rank = 1;
  _slice_bounds_desc.data_ptr = slice_bounds.data();
  _slice_bounds_desc.dims[0] = slice_bounds.size();
  bool _err{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_calc_bunch_params_z_slice(/* void* */ bunch.get_fortran_ptr(),
                                    /* void* */ bunch_params.get_fortran_ptr(),
                                    /* Bmad::array_descriptor_t& */ _slice_bounds_desc,
                                    /* bool& */ _err,
                                    /* bool* */ _print_err,
                                    /* bool* */ _is_time_coords,
                                    /* void* */ _ele);
  return _err;
}
BunchParamsStruct Bmad::calc_bunch_sigma_matrix_etc(
    CoordStructArray1D particle,
    FArray1D<Real> &charge,
    std::optional<bool> is_time_coords,
    optional_ref<EleStruct> ele
) {
  // particle: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _particle_desc;
  _particle_desc.rank = 1;
  _particle_desc.data_ptr = particle.data();
  _particle_desc.dims[0] = particle.size();
  _particle_desc.strides[0] = 1;
  // charge: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _charge_desc;
  _charge_desc.rank = 1;
  _charge_desc.data_ptr = charge.data();
  _charge_desc.dims[0] = charge.size();
  BunchParamsStruct _bunch_params;
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_calc_bunch_sigma_matrix_etc(
      /* Bmad::array_descriptor_t& */ _particle_desc,
      /* Bmad::array_descriptor_t& */ _charge_desc,
      /* void* */ _bunch_params.get_fortran_ptr(),
      /* bool* */ _is_time_coords,
      /* void* */ _ele
  );
  return std::move(_bunch_params);
}
Bmad::CalcEmittancesAndTwissFromSigmaMatrix Bmad::calc_emittances_and_twiss_from_sigma_matrix(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<bool> print_err
) {
  // sigma_mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  BunchParamsStruct _bunch_params;
  bool _error{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  // n_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _n_mat_desc;
  _n_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> n_mat;
  double _n_mat_vec[6 * 6];
  _n_mat_desc.data_ptr = _n_mat_vec;
  _n_mat_desc.dims[0] = 6;
  _n_mat_desc.dims[1] = 6;
  fortran_calc_emittances_and_twiss_from_sigma_matrix(
      /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
      /* void* */ _bunch_params.get_fortran_ptr(),
      /* bool& */ _error,
      /* bool* */ _print_err,
      /* Bmad::array_descriptor_t& */ _n_mat_desc
  );
  vec_to_matrix(_n_mat_vec, n_mat);
  return CalcEmittancesAndTwissFromSigmaMatrix{std::move(_bunch_params), _error, n_mat};
}
BunchParamsStruct Bmad::calc_spin_params(BunchStruct &bunch) {
  BunchParamsStruct _bunch_params;
  fortran_calc_spin_params(/* void* */ bunch.get_fortran_ptr(),
                           /* void* */ _bunch_params.get_fortran_ptr());
  return std::move(_bunch_params);
}
EleStruct Bmad::calc_super_slave_key(
    EleStruct &lord1,
    EleStruct &lord2,
    std::optional<bool> create_jumbo_slave
) {
  EleStruct _slave;
  bool create_jumbo_slave_lvalue;
  auto *_create_jumbo_slave{&create_jumbo_slave_lvalue};
  if (create_jumbo_slave.has_value()) {
    create_jumbo_slave_lvalue = create_jumbo_slave.value();
  } else {
    _create_jumbo_slave = nullptr;
  }
  fortran_calc_super_slave_key(/* void* */ lord1.get_fortran_ptr(),
                               /* void* */ lord2.get_fortran_ptr(),
                               /* void* */ _slave.get_fortran_ptr(),
                               /* bool* */ _create_jumbo_slave);
  return std::move(_slave);
}
Bmad::CalcWallRadius
Bmad::calc_wall_radius(Wall3dVertexStructArray1D v, double cos_ang, double sin_ang) {
  // v: Wall3dVertexStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _v_desc;
  _v_desc.rank = 1;
  _v_desc.data_ptr = v.data();
  _v_desc.dims[0] = v.size();
  _v_desc.strides[0] = 1;
  double _r_wall{};
  double _dr_dtheta{};
  int _ix_vertex{};
  fortran_calc_wall_radius(
      /* Bmad::array_descriptor_t& */ _v_desc,
      /* double& */ cos_ang,
      /* double& */ sin_ang,
      /* double& */ _r_wall,
      /* double& */ _dr_dtheta,
      /* int& */ _ix_vertex
  );
  return CalcWallRadius{_r_wall, _dr_dtheta, _ix_vertex};
}
void Bmad::calc_z_tune(BranchStruct &branch) {
  fortran_calc_z_tune(/* void* */ branch.get_fortran_ptr());
}
void Bmad::canonical_to_angle_coords(CoordStruct &orbit, std::optional<std::string> coord_type) {
  const char *_coord_type = coord_type.has_value() ? coord_type->c_str() : nullptr;
  fortran_canonical_to_angle_coords(/* void* */ orbit.get_fortran_ptr(),
                                    /* const char* */ _coord_type);
}
FixedArray2D<Real, 2, 2>
Bmad::cbar_to_c(FixedArray2D<Real, 2, 2> cbar_mat, TwissStruct &a, TwissStruct &b) {
  // cbar_mat: in NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _cbar_mat_desc;
  _cbar_mat_desc.rank = 2;
  double _cbar_mat_vec[2 * 2];
  _cbar_mat_desc.data_ptr = _cbar_mat_vec;
  _cbar_mat_desc.dims[0] = 2;
  _cbar_mat_desc.dims[1] = 2;
  matrix_to_vec(cbar_mat, _cbar_mat_vec);
  // c_mat: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _c_mat_desc;
  _c_mat_desc.rank = 2;
  FixedArray2D<Real, 2, 2> c_mat;
  double _c_mat_vec[2 * 2];
  _c_mat_desc.data_ptr = _c_mat_vec;
  _c_mat_desc.dims[0] = 2;
  _c_mat_desc.dims[1] = 2;
  fortran_cbar_to_c(
      /* Bmad::array_descriptor_t& */ _cbar_mat_desc,
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _c_mat_desc
  );
  vec_to_matrix(_c_mat_vec, c_mat);
  return c_mat;
}
void Bmad::check_aperture_limit(
    CoordStruct &orb,
    EleStruct &ele,
    int particle_at,
    LatParamStruct &param,
    optional_ref<CoordStruct> old_orb,
    std::optional<bool> check_momentum
) {
  auto *_old_orb =
      old_orb.has_value() ? old_orb->get().get_fortran_ptr() : nullptr; // input, optional
  bool check_momentum_lvalue;
  auto *_check_momentum{&check_momentum_lvalue};
  if (check_momentum.has_value()) {
    check_momentum_lvalue = check_momentum.value();
  } else {
    _check_momentum = nullptr;
  }
  fortran_check_aperture_limit(/* void* */ orb.get_fortran_ptr(),
                               /* void* */ ele.get_fortran_ptr(),
                               /* int& */ particle_at,
                               /* void* */ param.get_fortran_ptr(),
                               /* void* */ _old_orb,
                               /* bool* */ _check_momentum);
}
bool Bmad::check_controller_controls(int ele_key, ControlStructArray1D contrl, std::string name) {
  // contrl: ControlStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _contrl_desc;
  _contrl_desc.rank = 1;
  _contrl_desc.data_ptr = contrl.data();
  _contrl_desc.dims[0] = contrl.size();
  _contrl_desc.strides[0] = 1;
  auto _name = name.c_str();
  bool _err{};
  fortran_check_controller_controls(
      /* int& */ ele_key,
      /* Bmad::array_descriptor_t& */ _contrl_desc,
      /* const char* */ _name,
      /* bool& */ _err
  );
  return _err;
}
void Bmad::check_for_superimpose_problem(
    BranchStruct &branch,
    EleStruct &super_ele,
    bool err_flag,
    bool wrap,
    optional_ref<EleStruct> ref_ele
) {
  auto *_ref_ele =
      ref_ele.has_value() ? ref_ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_check_for_superimpose_problem(/* void* */ branch.get_fortran_ptr(),
                                        /* void* */ super_ele.get_fortran_ptr(),
                                        /* bool& */ err_flag,
                                        /* void* */ _ref_ele,
                                        /* bool& */ wrap);
}
Bmad::CheckIfSInBounds
Bmad::check_if_s_in_bounds(BranchStruct &branch, double s, std::optional<bool> print_err) {
  bool _err_flag{};
  double _translated_s{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_check_if_s_in_bounds(/* void* */ branch.get_fortran_ptr(),
                               /* double& */ s,
                               /* bool& */ _err_flag,
                               /* double& */ _translated_s,
                               /* bool* */ _print_err);
  return CheckIfSInBounds{_err_flag, _translated_s};
}
Bmad::ChooseQuadsForSetTune
Bmad::choose_quads_for_set_tune(BranchStruct &branch, std::optional<std::string> mask) {
  // intent=out allocatable general array
  auto dk1{RealAlloc1D()};
  // intent=out allocatable type array
  auto eles{ElePointerStructAlloc1D()};
  const char *_mask = mask.has_value() ? mask->c_str() : nullptr;
  bool _err_flag{};
  fortran_choose_quads_for_set_tune(/* void* */ branch.get_fortran_ptr(),
                                    /* void* */ dk1.get_fortran_ptr(),
                                    /* void* */ eles.get_fortran_ptr(),
                                    /* const char* */ _mask,
                                    /* bool& */ _err_flag);
  return ChooseQuadsForSetTune{std::move(dk1), std::move(eles), _err_flag};
}
Bmad::ChromCalc Bmad::chrom_calc(
    LatStruct &lat,
    double &delta_e,
    std::optional<double> pz,
    std::optional<int> ix_branch,
    optional_ref<CoordStruct> orb0
) {
  double _chrom_a{};
  double _chrom_b{};
  bool _err_flag{};
  double pz_lvalue;
  auto *_pz{&pz_lvalue};
  if (pz.has_value()) {
    pz_lvalue = pz.value();
  } else {
    _pz = nullptr;
  }
  LatStruct _low_E_lat;
  LatStruct _high_E_lat;
  // intent=out allocatable type array
  auto low_E_orb{CoordStructAlloc1D()};
  // intent=out allocatable type array
  auto high_E_orb{CoordStructAlloc1D()};
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  auto *_orb0 = orb0.has_value() ? orb0->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_chrom_calc(/* void* */ lat.get_fortran_ptr(),
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
      std::move(high_E_orb)
  };
}
bool Bmad::chrom_tune(
    LatStruct &lat,
    double &delta_e,
    double target_x,
    double target_y,
    double err_tol
) {
  bool _err_flag{};
  fortran_chrom_tune(/* void* */ lat.get_fortran_ptr(),
                     /* double& */ delta_e,
                     /* double& */ target_x,
                     /* double& */ target_y,
                     /* double& */ err_tol,
                     /* bool& */ _err_flag);
  return _err_flag;
}
double Bmad::classical_radius(int species) {
  double _radius{};
  fortran_classical_radius(/* int& */ species, /* double& */ _radius);
  return _radius;
}
LatStruct Bmad::clear_lat_1turn_mats() {
  LatStruct _lat;
  fortran_clear_lat_1turn_mats(/* void* */ _lat.get_fortran_ptr());
  return std::move(_lat);
}
void Bmad::clear_taylor_maps_from_elements(LatStruct &lat) {
  fortran_clear_taylor_maps_from_elements(/* void* */ lat.get_fortran_ptr());
}
bool Bmad::closed_orbit_calc(
    LatStruct &lat,
    CoordStructAlloc1D closed_orb,
    std::optional<int> i_dim,
    std::optional<int> direction,
    std::optional<int> ix_branch,
    std::optional<bool> print_err
) {
  // intent=inout allocatable type array
  int i_dim_lvalue;
  auto *_i_dim{&i_dim_lvalue};
  if (i_dim.has_value()) {
    i_dim_lvalue = i_dim.value();
  } else {
    _i_dim = nullptr;
  }
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_closed_orbit_calc(/* void* */ lat.get_fortran_ptr(),
                            /* void* */ closed_orb.get_fortran_ptr(),
                            /* int* */ _i_dim,
                            /* int* */ _direction,
                            /* int* */ _ix_branch,
                            /* bool& */ _err_flag,
                            /* bool* */ _print_err);
  return _err_flag;
}
Bmad::ClosedOrbitFromTracking Bmad::closed_orbit_from_tracking(
    LatStruct &lat,
    int i_dim,
    std::optional<FArray1D<Real>> eps_rel,
    std::optional<FArray1D<Real>> eps_abs,
    optional_ref<CoordStruct> init_guess
) {
  // intent=out allocatable type array
  auto closed_orb{CoordStructAlloc1D()};
  // eps_rel: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _eps_rel_desc;
  _eps_rel_desc.rank = 1;
  if (eps_rel.has_value()) {
    _eps_rel_desc.data_ptr = eps_rel->data();
    _eps_rel_desc.dims[0] = eps_rel->size();
  } else {
    _eps_rel_desc.data_ptr = nullptr;
    _eps_rel_desc.dims[0] = 0;
  }
  // eps_abs: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _eps_abs_desc;
  _eps_abs_desc.rank = 1;
  if (eps_abs.has_value()) {
    _eps_abs_desc.data_ptr = eps_abs->data();
    _eps_abs_desc.dims[0] = eps_abs->size();
  } else {
    _eps_abs_desc.data_ptr = nullptr;
    _eps_abs_desc.dims[0] = 0;
  }
  auto *_init_guess =
      init_guess.has_value() ? init_guess->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err_flag{};
  fortran_closed_orbit_from_tracking(/* void* */ lat.get_fortran_ptr(),
                                     /* void* */ closed_orb.get_fortran_ptr(),
                                     /* int& */ i_dim,
                                     /* Bmad::array_descriptor_t& */ _eps_rel_desc,
                                     /* Bmad::array_descriptor_t& */ _eps_abs_desc,
                                     /* void* */ _init_guess,
                                     /* bool& */ _err_flag);
  return ClosedOrbitFromTracking{std::move(closed_orb), _err_flag};
}
void Bmad::cmplx_re_str(std::complex<double> cmp, std::string str_out) {
  auto _str_out = str_out.c_str();
  fortran_cmplx_re_str(/* std::complex<double>& */ cmp, /* const char* */ _str_out);
}
bool Bmad::combine_consecutive_elements(LatStruct &lat) {
  bool _error{};
  fortran_combine_consecutive_elements(/* void* */ lat.get_fortran_ptr(), /* bool& */ _error);
  return _error;
}
void Bmad::complex_taylor_clean(ComplexTaylorStruct &complex_taylor) {
  fortran_complex_taylor_clean(/* void* */ complex_taylor.get_fortran_ptr());
}
void Bmad::complex_taylor_coef(
    ComplexTaylorStruct &complex_taylor,
    FArray1D<Int> &exp,
    std::complex<double> coef
) {
  // exp: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _exp_desc;
  _exp_desc.rank = 1;
  _exp_desc.data_ptr = exp.data();
  _exp_desc.dims[0] = exp.size();
  fortran_complex_taylor_coef1(/* void* */ complex_taylor.get_fortran_ptr(),
                               /* Bmad::array_descriptor_t& */ _exp_desc,
                               /* std::complex<double>& */ coef);
}
void Bmad::complex_taylor_coef(
    ComplexTaylorStruct &complex_taylor,
    std::complex<double> coef,
    std::optional<int> i1,
    std::optional<int> i2,
    std::optional<int> i3,
    std::optional<int> i4,
    std::optional<int> i5,
    std::optional<int> i6,
    std::optional<int> i7,
    std::optional<int> i8,
    std::optional<int> i9
) {
  int i1_lvalue;
  auto *_i1{&i1_lvalue};
  if (i1.has_value()) {
    i1_lvalue = i1.value();
  } else {
    _i1 = nullptr;
  }
  int i2_lvalue;
  auto *_i2{&i2_lvalue};
  if (i2.has_value()) {
    i2_lvalue = i2.value();
  } else {
    _i2 = nullptr;
  }
  int i3_lvalue;
  auto *_i3{&i3_lvalue};
  if (i3.has_value()) {
    i3_lvalue = i3.value();
  } else {
    _i3 = nullptr;
  }
  int i4_lvalue;
  auto *_i4{&i4_lvalue};
  if (i4.has_value()) {
    i4_lvalue = i4.value();
  } else {
    _i4 = nullptr;
  }
  int i5_lvalue;
  auto *_i5{&i5_lvalue};
  if (i5.has_value()) {
    i5_lvalue = i5.value();
  } else {
    _i5 = nullptr;
  }
  int i6_lvalue;
  auto *_i6{&i6_lvalue};
  if (i6.has_value()) {
    i6_lvalue = i6.value();
  } else {
    _i6 = nullptr;
  }
  int i7_lvalue;
  auto *_i7{&i7_lvalue};
  if (i7.has_value()) {
    i7_lvalue = i7.value();
  } else {
    _i7 = nullptr;
  }
  int i8_lvalue;
  auto *_i8{&i8_lvalue};
  if (i8.has_value()) {
    i8_lvalue = i8.value();
  } else {
    _i8 = nullptr;
  }
  int i9_lvalue;
  auto *_i9{&i9_lvalue};
  if (i9.has_value()) {
    i9_lvalue = i9.value();
  } else {
    _i9 = nullptr;
  }
  fortran_complex_taylor_coef2(/* void* */ complex_taylor.get_fortran_ptr(),
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
    ComplexTaylorStruct &complex_taylor1,
    ComplexTaylorStruct &complex_taylor2
) {
  fortran_complex_taylor_equal_complex_taylor(/* void* */ complex_taylor1.get_fortran_ptr(),
                                              /* void* */ complex_taylor2.get_fortran_ptr());
}
int Bmad::complex_taylor_exponent_index(FixedArray1D<Int, 6> expn) {
  // expn: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _expn_desc;
  _expn_desc.rank = 1;
  _expn_desc.data_ptr = expn.data();
  _expn_desc.dims[0] = expn.size();
  int _index{};
  fortran_complex_taylor_exponent_index(
      /* Bmad::array_descriptor_t& */ _expn_desc,
      /* int& */ _index
  );
  return _index;
}
void Bmad::complex_taylor_make_unit(ComplexTaylorStructArray1D complex_taylor) {
  // complex_taylor: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_desc;
  _complex_taylor_desc.rank = 1;
  _complex_taylor_desc.data_ptr = complex_taylor.data();
  _complex_taylor_desc.dims[0] = complex_taylor.size();
  _complex_taylor_desc.strides[0] = 1;
  fortran_complex_taylor_make_unit(/* Bmad::array_descriptor_t& */ _complex_taylor_desc);
}
Bmad::ComplexTaylorToMat6 Bmad::complex_taylor_to_mat6(
    ComplexTaylorStructArray1D a_complex_taylor,
    FArray1D<Complex> &r_in,
    std::optional<FArray1D<Complex>> r_out
) {
  // a_complex_taylor: ComplexTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _a_complex_taylor_desc;
  _a_complex_taylor_desc.rank = 1;
  _a_complex_taylor_desc.data_ptr = a_complex_taylor.data();
  _a_complex_taylor_desc.dims[0] = a_complex_taylor.size();
  _a_complex_taylor_desc.strides[0] = 1;
  // r_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r_in_desc;
  _r_in_desc.rank = 1;
  _r_in_desc.data_ptr = r_in.data();
  _r_in_desc.dims[0] = r_in.size();
  // vec0: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec0_desc;
  _vec0_desc.rank = 1;
  FixedArray1D<Complex, 6> _vec0;
  _vec0_desc.data_ptr = _vec0.data();
  _vec0_desc.dims[0] = _vec0.size();
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> mat6;
  std::complex<double> _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  // r_out: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r_out_desc;
  _r_out_desc.rank = 1;
  if (r_out.has_value()) {
    _r_out_desc.data_ptr = r_out->data();
    _r_out_desc.dims[0] = r_out->size();
  } else {
    _r_out_desc.data_ptr = nullptr;
    _r_out_desc.dims[0] = 0;
  }
  fortran_complex_taylor_to_mat6(
      /* Bmad::array_descriptor_t& */ _a_complex_taylor_desc,
      /* Bmad::array_descriptor_t& */ _r_in_desc,
      /* Bmad::array_descriptor_t& */ _vec0_desc,
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* Bmad::array_descriptor_t& */ _r_out_desc
  );
  vec_to_matrix(_mat6_vec, mat6);
  return ComplexTaylorToMat6{_vec0, mat6};
}
void Bmad::complex_taylors_equal_complex_taylors(
    ComplexTaylorStructArray1D complex_taylor1,
    ComplexTaylorStructArray1D complex_taylor2
) {
  // complex_taylor1: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor1_desc;
  _complex_taylor1_desc.rank = 1;
  _complex_taylor1_desc.data_ptr = complex_taylor1.data();
  _complex_taylor1_desc.dims[0] = complex_taylor1.size();
  _complex_taylor1_desc.strides[0] = 1;
  // complex_taylor2: ComplexTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor2_desc;
  _complex_taylor2_desc.rank = 1;
  _complex_taylor2_desc.data_ptr = complex_taylor2.data();
  _complex_taylor2_desc.dims[0] = complex_taylor2.size();
  _complex_taylor2_desc.strides[0] = 1;
  fortran_complex_taylors_equal_complex_taylors(
      /* Bmad::array_descriptor_t& */ _complex_taylor1_desc,
      /* Bmad::array_descriptor_t& */ _complex_taylor2_desc
  );
}
void Bmad::compute_slave_coupler(EleStruct &slave) {
  fortran_compute_slave_coupler(/* void* */ slave.get_fortran_ptr());
}
bool Bmad::concat_ele_taylor(
    TaylorStructArray1D orb_taylor,
    EleStruct &ele,
    std::optional<TaylorStructArray1D> spin_taylor
) {
  // orb_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb_taylor_desc;
  _orb_taylor_desc.rank = 1;
  _orb_taylor_desc.data_ptr = orb_taylor.data();
  _orb_taylor_desc.dims[0] = orb_taylor.size();
  _orb_taylor_desc.strides[0] = 1;
  bool _err_flag{};
  // spin_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_taylor_desc;
  _spin_taylor_desc.rank = 1;
  if (spin_taylor) {
    _spin_taylor_desc.data_ptr = spin_taylor->data();
    _spin_taylor_desc.dims[0] = spin_taylor->size();
  } else {
    _spin_taylor_desc.data_ptr = nullptr;
    _spin_taylor_desc.dims[0] = 0;
  }
  _spin_taylor_desc.strides[0] = 1;
  fortran_concat_ele_taylor(
      /* Bmad::array_descriptor_t& */ _orb_taylor_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* Bmad::array_descriptor_t& */ _spin_taylor_desc
  );
  return _err_flag;
}
void Bmad::concat_taylor(
    TaylorStructArray1D taylor1,
    TaylorStructArray1D taylor2,
    TaylorStructArray1D taylor3
) {
  // taylor1: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor1_desc;
  _taylor1_desc.rank = 1;
  _taylor1_desc.data_ptr = taylor1.data();
  _taylor1_desc.dims[0] = taylor1.size();
  _taylor1_desc.strides[0] = 1;
  // taylor2: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor2_desc;
  _taylor2_desc.rank = 1;
  _taylor2_desc.data_ptr = taylor2.data();
  _taylor2_desc.dims[0] = taylor2.size();
  _taylor2_desc.strides[0] = 1;
  // taylor3: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor3_desc;
  _taylor3_desc.rank = 1;
  _taylor3_desc.data_ptr = taylor3.data();
  _taylor3_desc.dims[0] = taylor3.size();
  _taylor3_desc.strides[0] = 1;
  fortran_concat_taylor(
      /* Bmad::array_descriptor_t& */ _taylor1_desc,
      /* Bmad::array_descriptor_t& */ _taylor2_desc,
      /* Bmad::array_descriptor_t& */ _taylor3_desc
  );
}
Bmad::ConcatTransferMat Bmad::concat_transfer_mat(
    FixedArray2D<Real, 6, 6> mat_1,
    FixedArray1D<Real, 6> vec_1,
    FixedArray2D<Real, 6, 6> mat_0,
    FixedArray1D<Real, 6> vec_0
) {
  // mat_1: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_1_desc;
  _mat_1_desc.rank = 2;
  double _mat_1_vec[6 * 6];
  _mat_1_desc.data_ptr = _mat_1_vec;
  _mat_1_desc.dims[0] = 6;
  _mat_1_desc.dims[1] = 6;
  matrix_to_vec(mat_1, _mat_1_vec);
  // vec_1: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec_1_desc;
  _vec_1_desc.rank = 1;
  _vec_1_desc.data_ptr = vec_1.data();
  _vec_1_desc.dims[0] = vec_1.size();
  // mat_0: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_0_desc;
  _mat_0_desc.rank = 2;
  double _mat_0_vec[6 * 6];
  _mat_0_desc.data_ptr = _mat_0_vec;
  _mat_0_desc.dims[0] = 6;
  _mat_0_desc.dims[1] = 6;
  matrix_to_vec(mat_0, _mat_0_vec);
  // vec_0: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec_0_desc;
  _vec_0_desc.rank = 1;
  _vec_0_desc.data_ptr = vec_0.data();
  _vec_0_desc.dims[0] = vec_0.size();
  // mat_out: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_out_desc;
  _mat_out_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat_out;
  double _mat_out_vec[6 * 6];
  _mat_out_desc.data_ptr = _mat_out_vec;
  _mat_out_desc.dims[0] = 6;
  _mat_out_desc.dims[1] = 6;
  // vec_out: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec_out_desc;
  _vec_out_desc.rank = 1;
  FixedArray1D<Real, 6> _vec_out;
  _vec_out_desc.data_ptr = _vec_out.data();
  _vec_out_desc.dims[0] = _vec_out.size();
  fortran_concat_transfer_mat(
      /* Bmad::array_descriptor_t& */ _mat_1_desc,
      /* Bmad::array_descriptor_t& */ _vec_1_desc,
      /* Bmad::array_descriptor_t& */ _mat_0_desc,
      /* Bmad::array_descriptor_t& */ _vec_0_desc,
      /* Bmad::array_descriptor_t& */ _mat_out_desc,
      /* Bmad::array_descriptor_t& */ _vec_out_desc
  );
  vec_to_matrix(_mat_out_vec, mat_out);
  return ConcatTransferMat{mat_out, _vec_out};
}
void Bmad::control_bookkeeper(
    LatStruct &lat,
    optional_ref<EleStruct> ele,
    std::optional<bool> err_flag
) {
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_control_bookkeeper(/* void* */ lat.get_fortran_ptr(),
                             /* void* */ _ele,
                             /* bool* */ _err_flag);
}
void Bmad::convert_bend_exact_multipole(
    double g,
    int out_type,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> an,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> bn
) {
  // an: inout NOT (CppWrapperGeneralArgumentArray) (['0:n_pole_maxx'])
  Bmad::array_descriptor_t _an_desc;
  _an_desc.rank = 1;
  _an_desc.data_ptr = an.data();
  _an_desc.dims[0] = an.size();
  // bn: inout NOT (CppWrapperGeneralArgumentArray) (['0:n_pole_maxx'])
  Bmad::array_descriptor_t _bn_desc;
  _bn_desc.rank = 1;
  _bn_desc.data_ptr = bn.data();
  _bn_desc.dims[0] = bn.size();
  fortran_convert_bend_exact_multipole(
      /* double& */ g,
      /* int& */ out_type,
      /* Bmad::array_descriptor_t& */ _an_desc,
      /* Bmad::array_descriptor_t& */ _bn_desc
  );
}
Bmad::ConvertCoords
Bmad::convert_coords(std::string in_type_str, CoordStruct &coord_in, EleStruct &ele) {
  auto _in_type_str = in_type_str.c_str();
  char _out_type_str[4096];
  CoordStruct _coord_out;
  bool _err_flag{};
  fortran_convert_coords(
      /* const char* */ _in_type_str,
      /* void* */ coord_in.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _out_type_str,
      /* void* */ _coord_out.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return ConvertCoords{_out_type_str, std::move(_coord_out), _err_flag};
}
EmFieldStruct Bmad::convert_field_ele_to_lab(
    EleStruct &ele,
    double s_here,
    bool forward_transform,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential
) {
  EmFieldStruct _field;
  bool calc_dfield_lvalue;
  auto *_calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto *_calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_convert_field_ele_to_lab(/* void* */ ele.get_fortran_ptr(),
                                   /* double& */ s_here,
                                   /* bool& */ forward_transform,
                                   /* void* */ _field.get_fortran_ptr(),
                                   /* bool* */ _calc_dfield,
                                   /* bool* */ _calc_potential);
  return std::move(_field);
}
void Bmad::convert_local_cartesian_to_local_curvilinear(
    double x,
    double z,
    double g,
    double xout,
    double sout
) {
  fortran_convert_local_cartesian_to_local_curvilinear(
      /* double& */ x,
      /* double& */ z,
      /* double& */ g,
      /* double& */ xout,
      /* double& */ sout
  );
}
void Bmad::convert_local_curvilinear_to_local_cartesian(
    double x,
    double s,
    double g,
    double xout,
    double zout
) {
  fortran_convert_local_curvilinear_to_local_cartesian(
      /* double& */ x,
      /* double& */ s,
      /* double& */ g,
      /* double& */ xout,
      /* double& */ zout
  );
}
void Bmad::convert_particle_coordinates_s_to_t(
    CoordStruct &particle,
    double s_body,
    int orientation
) {
  fortran_convert_particle_coordinates_s_to_t(/* void* */ particle.get_fortran_ptr(),
                                              /* double& */ s_body,
                                              /* int& */ orientation);
}
double Bmad::convert_particle_coordinates_t_to_s(
    CoordStruct &particle,
    EleStruct &ele,
    std::optional<bool> use_downstream_p0c
) {
  double _s_body{};
  bool use_downstream_p0c_lvalue;
  auto *_use_downstream_p0c{&use_downstream_p0c_lvalue};
  if (use_downstream_p0c.has_value()) {
    use_downstream_p0c_lvalue = use_downstream_p0c.value();
  } else {
    _use_downstream_p0c = nullptr;
  }
  fortran_convert_particle_coordinates_t_to_s(/* void* */ particle.get_fortran_ptr(),
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
      /* bool& */ _err_flag
  );
  return ConvertPcTo{_E_tot, _gamma, _kinetic, _beta, _brho, _beta1, _err_flag};
}
Bmad::ConvertTotalEnergyTo
Bmad::convert_total_energy_to(double E_tot, int particle, std::optional<bool> print_err) {
  double _gamma{};
  double _kinetic{};
  double _beta{};
  double _pc{};
  double _brho{};
  double _beta1{};
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
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
      /* bool* */ _print_err
  );
  return ConvertTotalEnergyTo{_gamma, _kinetic, _beta, _pc, _brho, _beta1, _err_flag};
}
Bmad::ConverterDistributionParser Bmad::converter_distribution_parser(EleStruct &ele) {
  char _delim[4096];
  bool _delim_found{};
  bool _err_flag{};
  fortran_converter_distribution_parser(/* void* */ ele.get_fortran_ptr(),
                                        /* const char* */ _delim,
                                        /* bool& */ _delim_found,
                                        /* bool& */ _err_flag);
  return ConverterDistributionParser{_delim, _delim_found, _err_flag};
}
CoordStruct Bmad::coord_equal_coord(CoordStruct &coord2) {
  CoordStruct _coord1;
  fortran_coord_equal_coord(/* void* */ _coord1.get_fortran_ptr(),
                            /* void* */ coord2.get_fortran_ptr());
  return std::move(_coord1);
}
std::string Bmad::coord_state_name(int coord_state, std::optional<bool> one_word) {
  bool one_word_lvalue;
  auto *_one_word{&one_word_lvalue};
  if (one_word.has_value()) {
    one_word_lvalue = one_word.value();
  } else {
    _one_word = nullptr;
  }
  char _state_str[4096];
  fortran_coord_state_name(
      /* int& */ coord_state,
      /* bool* */ _one_word,
      /* const char* */ _state_str
  );
  return _state_str;
}
Bmad::CoordsBodyToLocal Bmad::coords_body_to_local(
    FloorPositionStruct &body_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles
) {
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  bool calculate_angles_lvalue;
  auto *_calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  FloorPositionStruct _local_position;
  fortran_coords_body_to_local(/* void* */ body_position.get_fortran_ptr(),
                               /* void* */ ele.get_fortran_ptr(),
                               /* Bmad::array_descriptor_t& */ _w_mat_desc,
                               /* bool* */ _calculate_angles,
                               /* void* */ _local_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsBodyToLocal{w_mat, std::move(_local_position)};
}
Bmad::CoordsBodyToRelExit Bmad::coords_body_to_rel_exit(
    FloorPositionStruct &body_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles
) {
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  bool calculate_angles_lvalue;
  auto *_calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  FloorPositionStruct _rel_exit;
  fortran_coords_body_to_rel_exit(/* void* */ body_position.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* Bmad::array_descriptor_t& */ _w_mat_desc,
                                  /* bool* */ _calculate_angles,
                                  /* void* */ _rel_exit.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsBodyToRelExit{w_mat, std::move(_rel_exit)};
}
Bmad::CoordsCurvilinearToFloor
Bmad::coords_curvilinear_to_floor(FixedArray1D<Real, 3> xys, BranchStruct &branch) {
  // xys: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _xys_desc;
  _xys_desc.rank = 1;
  _xys_desc.data_ptr = xys.data();
  _xys_desc.dims[0] = xys.size();
  bool _err_flag{};
  FloorPositionStruct _global;
  fortran_coords_curvilinear_to_floor(
      /* Bmad::array_descriptor_t& */ _xys_desc,
      /* void* */ branch.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _global.get_fortran_ptr()
  );
  return CoordsCurvilinearToFloor{_err_flag, std::move(_global)};
}
Bmad::CoordsFloorToCurvilinear
Bmad::coords_floor_to_curvilinear(FloorPositionStruct &floor_coords, EleStruct &ele0) {
  void *_ele1;
  int _status{};
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  FloorPositionStruct _local_coords;
  fortran_coords_floor_to_curvilinear(/* void* */ floor_coords.get_fortran_ptr(),
                                      /* void* */ ele0.get_fortran_ptr(),
                                      /* void* */ &_ele1,
                                      /* int& */ _status,
                                      /* Bmad::array_descriptor_t& */ _w_mat_desc,
                                      /* void* */ _local_coords.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsFloorToCurvilinear{
      std::move((_ele1 ? std::make_optional<EleStruct>(_ele1) : std::nullopt)),
      _status,
      w_mat,
      std::move(_local_coords)
  };
}
Bmad::CoordsFloorToLocalCurvilinear Bmad::coords_floor_to_local_curvilinear(
    FloorPositionStruct &global_position,
    EleStruct &ele,
    std::optional<int> relative_to
) {
  int _status{};
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  int relative_to_lvalue;
  auto *_relative_to{&relative_to_lvalue};
  if (relative_to.has_value()) {
    relative_to_lvalue = relative_to.value();
  } else {
    _relative_to = nullptr;
  }
  FloorPositionStruct _local_position;
  fortran_coords_floor_to_local_curvilinear(/* void* */ global_position.get_fortran_ptr(),
                                            /* void* */ ele.get_fortran_ptr(),
                                            /* int& */ _status,
                                            /* Bmad::array_descriptor_t& */ _w_mat_desc,
                                            /* int* */ _relative_to,
                                            /* void* */ _local_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsFloorToLocalCurvilinear{_status, w_mat, std::move(_local_position)};
}
FloorPositionStruct Bmad::coords_floor_to_relative(
    FloorPositionStruct &floor0,
    FloorPositionStruct &global_position,
    std::optional<bool> calculate_angles,
    std::optional<bool> is_delta_position
) {
  bool calculate_angles_lvalue;
  auto *_calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  bool is_delta_position_lvalue;
  auto *_is_delta_position{&is_delta_position_lvalue};
  if (is_delta_position.has_value()) {
    is_delta_position_lvalue = is_delta_position.value();
  } else {
    _is_delta_position = nullptr;
  }
  FloorPositionStruct _local_position;
  fortran_coords_floor_to_relative(/* void* */ floor0.get_fortran_ptr(),
                                   /* void* */ global_position.get_fortran_ptr(),
                                   /* bool* */ _calculate_angles,
                                   /* bool* */ _is_delta_position,
                                   /* void* */ _local_position.get_fortran_ptr());
  return std::move(_local_position);
}
Bmad::CoordsLocalCurvilinearToBody Bmad::coords_local_curvilinear_to_body(
    FloorPositionStruct &local_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles
) {
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  bool calculate_angles_lvalue;
  auto *_calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  FloorPositionStruct _body_position;
  fortran_coords_local_curvilinear_to_body(/* void* */ local_position.get_fortran_ptr(),
                                           /* void* */ ele.get_fortran_ptr(),
                                           /* Bmad::array_descriptor_t& */ _w_mat_desc,
                                           /* bool* */ _calculate_angles,
                                           /* void* */ _body_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsLocalCurvilinearToBody{w_mat, std::move(_body_position)};
}
Bmad::CoordsLocalCurvilinearToFloor Bmad::coords_local_curvilinear_to_floor(
    FloorPositionStruct &local_position,
    EleStruct &ele,
    std::optional<bool> in_body_frame,
    std::optional<bool> calculate_angles,
    std::optional<int> end_origin,
    std::optional<bool> downstream_dir_ref
) {
  bool in_body_frame_lvalue;
  auto *_in_body_frame{&in_body_frame_lvalue};
  if (in_body_frame.has_value()) {
    in_body_frame_lvalue = in_body_frame.value();
  } else {
    _in_body_frame = nullptr;
  }
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  bool calculate_angles_lvalue;
  auto *_calculate_angles{&calculate_angles_lvalue};
  if (calculate_angles.has_value()) {
    calculate_angles_lvalue = calculate_angles.value();
  } else {
    _calculate_angles = nullptr;
  }
  int end_origin_lvalue;
  auto *_end_origin{&end_origin_lvalue};
  if (end_origin.has_value()) {
    end_origin_lvalue = end_origin.value();
  } else {
    _end_origin = nullptr;
  }
  bool downstream_dir_ref_lvalue;
  auto *_downstream_dir_ref{&downstream_dir_ref_lvalue};
  if (downstream_dir_ref.has_value()) {
    downstream_dir_ref_lvalue = downstream_dir_ref.value();
  } else {
    _downstream_dir_ref = nullptr;
  }
  FloorPositionStruct _global_position;
  fortran_coords_local_curvilinear_to_floor(/* void* */ local_position.get_fortran_ptr(),
                                            /* void* */ ele.get_fortran_ptr(),
                                            /* bool* */ _in_body_frame,
                                            /* Bmad::array_descriptor_t& */ _w_mat_desc,
                                            /* bool* */ _calculate_angles,
                                            /* int* */ _end_origin,
                                            /* bool* */ _downstream_dir_ref,
                                            /* void* */ _global_position.get_fortran_ptr());
  vec_to_matrix(_w_mat_vec, w_mat);
  return CoordsLocalCurvilinearToFloor{w_mat, std::move(_global_position)};
}
FloorPositionStruct Bmad::coords_relative_to_floor(
    FloorPositionStruct &floor0,
    FixedArray1D<Real, 3> dr,
    std::optional<double> theta,
    std::optional<double> phi,
    std::optional<double> psi
) {
  // dr: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _dr_desc;
  _dr_desc.rank = 1;
  _dr_desc.data_ptr = dr.data();
  _dr_desc.dims[0] = dr.size();
  double theta_lvalue;
  auto *_theta{&theta_lvalue};
  if (theta.has_value()) {
    theta_lvalue = theta.value();
  } else {
    _theta = nullptr;
  }
  double phi_lvalue;
  auto *_phi{&phi_lvalue};
  if (phi.has_value()) {
    phi_lvalue = phi.value();
  } else {
    _phi = nullptr;
  }
  double psi_lvalue;
  auto *_psi{&psi_lvalue};
  if (psi.has_value()) {
    psi_lvalue = psi.value();
  } else {
    _psi = nullptr;
  }
  FloorPositionStruct _floor1;
  fortran_coords_relative_to_floor(/* void* */ floor0.get_fortran_ptr(),
                                   /* Bmad::array_descriptor_t& */ _dr_desc,
                                   /* double* */ _theta,
                                   /* double* */ _phi,
                                   /* double* */ _psi,
                                   /* void* */ _floor1.get_fortran_ptr());
  return std::move(_floor1);
}
void Bmad::coulombfun(double u, double v, double w, double gam, double res) {
  fortran_coulombfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ res
  );
}
void Bmad::create_concatenated_wall3d(LatStruct &lat, bool err) {
  fortran_create_concatenated_wall3d(/* void* */ lat.get_fortran_ptr(), /* bool& */ err);
}
Bmad::CreateElementSlice Bmad::create_element_slice(
    EleStruct &ele_in,
    double l_slice,
    double offset,
    LatParamStruct &param,
    bool include_upstream_end,
    bool include_downstream_end,
    optional_ref<EleStruct> old_slice,
    optional_ref<CoordStruct> orb_in
) {
  EleStruct _sliced_ele;
  bool _err_flag{};
  auto *_old_slice =
      old_slice.has_value() ? old_slice->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_orb_in = orb_in.has_value() ? orb_in->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_create_element_slice(/* void* */ _sliced_ele.get_fortran_ptr(),
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
void Bmad::create_feedback(
    EleStruct &lord,
    CharacterAlloc1D &input,
    CharacterAlloc1D &output,
    bool err_flag
) {
  // intent=in character array container
  // intent=in character array container
  fortran_create_feedback(/* void* */ lord.get_fortran_ptr(),
                          /* void* */ input.get_fortran_ptr(),
                          /* void* */ output.get_fortran_ptr(),
                          /* bool& */ err_flag);
}
bool Bmad::create_field_overlap(LatStruct &lat, std::string lord_name, std::string slave_name) {
  auto _lord_name = lord_name.c_str();
  auto _slave_name = slave_name.c_str();
  bool _err_flag{};
  fortran_create_field_overlap(/* void* */ lat.get_fortran_ptr(),
                               /* const char* */ _lord_name,
                               /* const char* */ _slave_name,
                               /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::create_girder(
    LatStruct &lat,
    int ix_girder,
    ControlStructArray1D contrl,
    EleStruct &girder_info,
    bool err_flag
) {
  // contrl: ControlStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _contrl_desc;
  _contrl_desc.rank = 1;
  _contrl_desc.data_ptr = contrl.data();
  _contrl_desc.dims[0] = contrl.size();
  _contrl_desc.strides[0] = 1;
  fortran_create_girder(/* void* */ lat.get_fortran_ptr(),
                        /* int& */ ix_girder,
                        /* Bmad::array_descriptor_t& */ _contrl_desc,
                        /* void* */ girder_info.get_fortran_ptr(),
                        /* bool& */ err_flag);
}
void Bmad::create_group(EleStruct &lord, ControlStructArray1D contrl, bool err) {
  // contrl: ControlStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _contrl_desc;
  _contrl_desc.rank = 1;
  _contrl_desc.data_ptr = contrl.data();
  _contrl_desc.dims[0] = contrl.size();
  _contrl_desc.strides[0] = 1;
  fortran_create_group(/* void* */ lord.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _contrl_desc,
                       /* bool& */ err);
}
NametableStruct Bmad::create_lat_ele_nametable(LatStruct &lat) {
  NametableStruct _nametable;
  fortran_create_lat_ele_nametable(/* void* */ lat.get_fortran_ptr(),
                                   /* void* */ _nametable.get_fortran_ptr());
  return std::move(_nametable);
}
void Bmad::create_overlay(EleStruct &lord, ControlStructArray1D contrl, bool err) {
  // contrl: ControlStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _contrl_desc;
  _contrl_desc.rank = 1;
  _contrl_desc.data_ptr = contrl.data();
  _contrl_desc.dims[0] = contrl.size();
  _contrl_desc.strides[0] = 1;
  fortran_create_overlay(/* void* */ lord.get_fortran_ptr(),
                         /* Bmad::array_descriptor_t& */ _contrl_desc,
                         /* bool& */ err);
}
Bmad::CreatePlanarWigglerModel
Bmad::create_planar_wiggler_model(EleStruct &wiggler_in, std::optional<bool> print_err) {
  LatStruct _lat;
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_create_planar_wiggler_model(/* void* */ wiggler_in.get_fortran_ptr(),
                                      /* void* */ _lat.get_fortran_ptr(),
                                      /* bool& */ _err_flag,
                                      /* bool* */ _print_err);
  return CreatePlanarWigglerModel{std::move(_lat), _err_flag};
}
void Bmad::create_ramper(EleStruct &lord, ControlStructArray1D contrl, bool err) {
  // contrl: ControlStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _contrl_desc;
  _contrl_desc.rank = 1;
  _contrl_desc.data_ptr = contrl.data();
  _contrl_desc.dims[0] = contrl.size();
  _contrl_desc.strides[0] = 1;
  fortran_create_ramper(/* void* */ lord.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _contrl_desc,
                        /* bool& */ err);
}
void Bmad::create_sol_quad_model(EleStruct &sol_quad, LatStruct &lat) {
  fortran_create_sol_quad_model(/* void* */ sol_quad.get_fortran_ptr(),
                                /* void* */ lat.get_fortran_ptr());
}
void Bmad::create_unique_ele_names(LatStruct &lat, int key, std::string suffix) {
  auto _suffix = suffix.c_str();
  fortran_create_unique_ele_names(/* void* */ lat.get_fortran_ptr(),
                                  /* int& */ key,
                                  /* const char* */ _suffix);
}
CartesianMapStruct Bmad::create_wiggler_cartesian_map(EleStruct &ele) {
  CartesianMapStruct _cart_map;
  fortran_create_wiggler_cartesian_map(/* void* */ ele.get_fortran_ptr(),
                                       /* void* */ _cart_map.get_fortran_ptr());
  return std::move(_cart_map);
}
void Bmad::crystal_attribute_bookkeeper(EleStruct &ele) {
  fortran_crystal_attribute_bookkeeper(/* void* */ ele.get_fortran_ptr());
}
void Bmad::crystal_h_misalign(EleStruct &ele, CoordStruct &orbit, FixedArray1D<Real, 3> h_vec) {
  // h_vec: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _h_vec_desc;
  _h_vec_desc.rank = 1;
  _h_vec_desc.data_ptr = h_vec.data();
  _h_vec_desc.dims[0] = h_vec.size();
  fortran_crystal_h_misalign(/* void* */ ele.get_fortran_ptr(),
                             /* void* */ orbit.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _h_vec_desc);
}
bool Bmad::crystal_type_to_crystal_params(EleStruct &ele) {
  bool _err_flag{};
  fortran_crystal_type_to_crystal_params(/* void* */ ele.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
int Bmad::custom_attribute_ubound_index(int ele_class) {
  int _ix_ubound{};
  fortran_custom_attribute_ubound_index(/* int& */ ele_class, /* int& */ _ix_ubound);
  return _ix_ubound;
}
Bmad::CustomEleAttribNameList Bmad::custom_ele_attrib_name_list() {
  // intent=out allocatable general array
  auto index_list{IntAlloc1D()};
  // intent=out character array container
  auto name_list{CharacterAlloc1D()};
  fortran_custom_ele_attrib_name_list(/* void* */ index_list.get_fortran_ptr(),
                                      /* void* */ name_list.get_fortran_ptr());
  return CustomEleAttribNameList{std::move(index_list), std::move(name_list)};
}
void Bmad::damping_matrix_d(
    double gamma,
    double g_tot,
    double B0,
    double B1,
    double delta,
    int species,
    FixedArray2D<Real, 6, 6> mat
) {
  // mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_desc;
  _mat_desc.rank = 2;
  double _mat_vec[6 * 6];
  _mat_desc.data_ptr = _mat_vec;
  _mat_desc.dims[0] = 6;
  _mat_desc.dims[1] = 6;
  matrix_to_vec(mat, _mat_vec);
  fortran_damping_matrix_d(
      /* double& */ gamma,
      /* double& */ g_tot,
      /* double& */ B0,
      /* double& */ B1,
      /* double& */ delta,
      /* int& */ species,
      /* Bmad::array_descriptor_t& */ _mat_desc
  );
  vec_to_matrix(_mat_vec, mat);
}
void Bmad::deallocate_ele_pointers(
    EleStruct &ele,
    std::optional<bool> nullify_only,
    std::optional<bool> nullify_branch,
    std::optional<bool> dealloc_poles
) {
  bool nullify_only_lvalue;
  auto *_nullify_only{&nullify_only_lvalue};
  if (nullify_only.has_value()) {
    nullify_only_lvalue = nullify_only.value();
  } else {
    _nullify_only = nullptr;
  }
  bool nullify_branch_lvalue;
  auto *_nullify_branch{&nullify_branch_lvalue};
  if (nullify_branch.has_value()) {
    nullify_branch_lvalue = nullify_branch.value();
  } else {
    _nullify_branch = nullptr;
  }
  bool dealloc_poles_lvalue;
  auto *_dealloc_poles{&dealloc_poles_lvalue};
  if (dealloc_poles.has_value()) {
    dealloc_poles_lvalue = dealloc_poles.value();
  } else {
    _dealloc_poles = nullptr;
  }
  fortran_deallocate_ele_pointers(/* void* */ ele.get_fortran_ptr(),
                                  /* bool* */ _nullify_only,
                                  /* bool* */ _nullify_branch,
                                  /* bool* */ _dealloc_poles);
}
void Bmad::deallocate_expression_tree(ExpressionTreeStruct &tree) {
  fortran_deallocate_expression_tree(/* void* */ tree.get_fortran_ptr());
}
void Bmad::deallocate_lat_pointers(LatStruct &lat) {
  fortran_deallocate_lat_pointers(/* void* */ lat.get_fortran_ptr());
}
int Bmad::default_tracking_species(LatParamStruct &param) {
  int _species{};
  fortran_default_tracking_species(/* void* */ param.get_fortran_ptr(), /* int& */ _species);
  return _species;
}
FixedArray1D<Int, 2> Bmad::detector_pixel_pt(CoordStruct &orbit, EleStruct &ele) {
  // ix_pix: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _ix_pix_desc;
  _ix_pix_desc.rank = 1;
  FixedArray1D<Int, 2> _ix_pix;
  _ix_pix_desc.data_ptr = _ix_pix.data();
  _ix_pix_desc.dims[0] = _ix_pix.size();
  fortran_detector_pixel_pt(/* void* */ orbit.get_fortran_ptr(),
                            /* void* */ ele.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _ix_pix_desc);
  return _ix_pix;
}
int Bmad::diffraction_plate_or_mask_hit_spot(EleStruct &ele, CoordStruct &orbit) {
  int _ix_section{};
  fortran_diffraction_plate_or_mask_hit_spot(/* void* */ ele.get_fortran_ptr(),
                                             /* void* */ orbit.get_fortran_ptr(),
                                             /* int& */ _ix_section);
  return _ix_section;
}
void Bmad::diffusion_matrix_b(
    double gamma,
    double g_tot,
    int species,
    FixedArray2D<Real, 6, 6> mat
) {
  // mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_desc;
  _mat_desc.rank = 2;
  double _mat_vec[6 * 6];
  _mat_desc.data_ptr = _mat_vec;
  _mat_desc.dims[0] = 6;
  _mat_desc.dims[1] = 6;
  matrix_to_vec(mat, _mat_vec);
  fortran_diffusion_matrix_b(
      /* double& */ gamma,
      /* double& */ g_tot,
      /* int& */ species,
      /* Bmad::array_descriptor_t& */ _mat_desc
  );
  vec_to_matrix(_mat_vec, mat);
}
Bmad::DistanceToAperture
Bmad::distance_to_aperture(CoordStruct &orbit, int particle_at, EleStruct &ele) {
  bool _no_aperture_here{};
  double _dist{};
  fortran_distance_to_aperture(/* void* */ orbit.get_fortran_ptr(),
                               /* int& */ particle_at,
                               /* void* */ ele.get_fortran_ptr(),
                               /* bool& */ _no_aperture_here,
                               /* double& */ _dist);
  return DistanceToAperture{_no_aperture_here, _dist};
}
bool Bmad::do_mode_flip(EleStruct &ele) {
  bool _err_flag{};
  fortran_do_mode_flip(/* void* */ ele.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::dpc_given_de(double pc_old, double mass, double dE, double dpc) {
  fortran_dpc_given_de(
      /* double& */ pc_old,
      /* double& */ mass,
      /* double& */ dE,
      /* double& */ dpc
  );
}
void Bmad::drift_and_pipe_track_methods_adjustment(LatStruct &lat) {
  fortran_drift_and_pipe_track_methods_adjustment(/* void* */ lat.get_fortran_ptr());
}
void Bmad::drift_multipass_name_correction(LatStruct &lat) {
  fortran_drift_multipass_name_correction(/* void* */ lat.get_fortran_ptr());
}
void Bmad::drift_orbit_time(
    CoordStruct &orbit,
    double beta0,
    std::optional<double> delta_s,
    std::optional<double> delta_t
) {
  double delta_s_lvalue;
  auto *_delta_s{&delta_s_lvalue};
  if (delta_s.has_value()) {
    delta_s_lvalue = delta_s.value();
  } else {
    _delta_s = nullptr;
  }
  double delta_t_lvalue;
  auto *_delta_t{&delta_t_lvalue};
  if (delta_t.has_value()) {
    delta_t_lvalue = delta_t.value();
  } else {
    _delta_t = nullptr;
  }
  fortran_drift_orbit_time(/* void* */ orbit.get_fortran_ptr(),
                           /* double& */ beta0,
                           /* double* */ _delta_s,
                           /* double* */ _delta_t);
}
void Bmad::drift_particle_to_s(CoordStruct &p, double s, BranchStruct &branch) {
  fortran_drift_particle_to_s(/* void* */ p.get_fortran_ptr(),
                              /* double& */ s,
                              /* void* */ branch.get_fortran_ptr());
}
void Bmad::drift_particle_to_t(CoordStruct &p, double t, BranchStruct &branch) {
  fortran_drift_particle_to_t(/* void* */ p.get_fortran_ptr(),
                              /* double& */ t,
                              /* void* */ branch.get_fortran_ptr());
}
double Bmad::dspline_len(
    double s_chord0,
    double s_chord1,
    SplineStruct &spline,
    std::optional<double> dtheta_ref
) {
  double dtheta_ref_lvalue;
  auto *_dtheta_ref{&dtheta_ref_lvalue};
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
      /* double& */ _dlen
  );
  return _dlen;
}
AperturePointStruct Bmad::dynamic_aperture_point(
    BranchStruct &branch,
    EleStruct &ele0,
    CoordStruct &orb0,
    double theta_xy,
    ApertureParamStruct &ap_param,
    std::optional<bool> check_xy_init
) {
  AperturePointStruct _ap_point;
  bool check_xy_init_lvalue;
  auto *_check_xy_init{&check_xy_init_lvalue};
  if (check_xy_init.has_value()) {
    check_xy_init_lvalue = check_xy_init.value();
  } else {
    _check_xy_init = nullptr;
  }
  fortran_dynamic_aperture_point(/* void* */ branch.get_fortran_ptr(),
                                 /* void* */ ele0.get_fortran_ptr(),
                                 /* void* */ orb0.get_fortran_ptr(),
                                 /* double& */ theta_xy,
                                 /* void* */ ap_param.get_fortran_ptr(),
                                 /* void* */ _ap_point.get_fortran_ptr(),
                                 /* bool* */ _check_xy_init);
  return std::move(_ap_point);
}
ApertureScanStructAlloc1D Bmad::dynamic_aperture_scan(
    ApertureParamStruct &aperture_param,
    FArray1D<Real> &pz_start,
    LatStruct &lat,
    std::optional<bool> print_timing
) {
  // intent=out allocatable type array
  auto aperture_scan{ApertureScanStructAlloc1D()};
  // pz_start: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _pz_start_desc;
  _pz_start_desc.rank = 1;
  _pz_start_desc.data_ptr = pz_start.data();
  _pz_start_desc.dims[0] = pz_start.size();
  bool print_timing_lvalue;
  auto *_print_timing{&print_timing_lvalue};
  if (print_timing.has_value()) {
    print_timing_lvalue = print_timing.value();
  } else {
    _print_timing = nullptr;
  }
  fortran_dynamic_aperture_scan(/* void* */ aperture_scan.get_fortran_ptr(),
                                /* void* */ aperture_param.get_fortran_ptr(),
                                /* Bmad::array_descriptor_t& */ _pz_start_desc,
                                /* void* */ lat.get_fortran_ptr(),
                                /* bool* */ _print_timing);
  return std::move(aperture_scan);
}
double Bmad::e_accel_field(
    EleStruct &ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking
) {
  bool bmad_standard_tracking_lvalue;
  auto *_bmad_standard_tracking{&bmad_standard_tracking_lvalue};
  if (bmad_standard_tracking.has_value()) {
    bmad_standard_tracking_lvalue = bmad_standard_tracking.value();
  } else {
    _bmad_standard_tracking = nullptr;
  }
  double _field{};
  fortran_e_accel_field(/* void* */ ele.get_fortran_ptr(),
                        /* int& */ voltage_or_gradient,
                        /* bool* */ _bmad_standard_tracking,
                        /* double& */ _field);
  return _field;
}
double Bmad::e_crit_photon(double gamma, double g_bend) {
  double _E_crit{};
  fortran_e_crit_photon(/* double& */ gamma, /* double& */ g_bend, /* double& */ _E_crit);
  return _E_crit;
}
Bmad::EigenDecomp6mat Bmad::eigen_decomp_6mat(FixedArray2D<Real, 6, 6> mat) {
  // mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat_desc;
  _mat_desc.rank = 2;
  double _mat_vec[6 * 6];
  _mat_desc.data_ptr = _mat_vec;
  _mat_desc.dims[0] = 6;
  _mat_desc.dims[1] = 6;
  matrix_to_vec(mat, _mat_vec);
  // eval: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _eval_desc;
  _eval_desc.rank = 1;
  FixedArray1D<Complex, 6> _eval;
  _eval_desc.data_ptr = _eval.data();
  _eval_desc.dims[0] = _eval.size();
  // evec: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _evec_desc;
  _evec_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> evec;
  std::complex<double> _evec_vec[6 * 6];
  _evec_desc.data_ptr = _evec_vec;
  _evec_desc.dims[0] = 6;
  _evec_desc.dims[1] = 6;
  bool _err_flag{};
  // tunes: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _tunes_desc;
  _tunes_desc.rank = 1;
  FixedArray1D<Real, 3> _tunes;
  _tunes_desc.data_ptr = _tunes.data();
  _tunes_desc.dims[0] = _tunes.size();
  fortran_eigen_decomp_6mat(
      /* Bmad::array_descriptor_t& */ _mat_desc,
      /* Bmad::array_descriptor_t& */ _eval_desc,
      /* Bmad::array_descriptor_t& */ _evec_desc,
      /* bool& */ _err_flag,
      /* Bmad::array_descriptor_t& */ _tunes_desc
  );
  vec_to_matrix(_evec_vec, evec);
  return EigenDecomp6mat{_eval, evec, _err_flag, _tunes};
}
void Bmad::ele_compute_ref_energy_and_time(
    EleStruct &ele0,
    EleStruct &ele,
    LatParamStruct &param,
    bool err_flag
) {
  fortran_ele_compute_ref_energy_and_time(/* void* */ ele0.get_fortran_ptr(),
                                          /* void* */ ele.get_fortran_ptr(),
                                          /* void* */ param.get_fortran_ptr(),
                                          /* bool& */ err_flag);
}
void Bmad::ele_equal_ele(EleStruct &ele_out, EleStruct &ele_in) {
  fortran_ele_equal_ele(/* void* */ ele_out.get_fortran_ptr(),
                        /* void* */ ele_in.get_fortran_ptr());
}
EleStruct Bmad::ele_equals_ele(EleStruct &ele_in, bool update_nametable) {
  EleStruct _ele_out;
  fortran_ele_equals_ele(/* void* */ _ele_out.get_fortran_ptr(),
                         /* void* */ ele_in.get_fortran_ptr(),
                         /* bool& */ update_nametable);
  return std::move(_ele_out);
}
void Bmad::ele_finalizer(EleStruct &ele) {
  fortran_ele_finalizer(/* void* */ ele.get_fortran_ptr());
}
std::string Bmad::ele_full_name(EleStruct &ele, std::optional<std::string> template_) {
  const char *_template_ = template_.has_value() ? template_->c_str() : nullptr;
  char _str[4096];
  fortran_ele_full_name(/* void* */ ele.get_fortran_ptr(),
                        /* const char* */ _template_,
                        /* const char* */ _str);
  return _str;
}
FloorPositionStruct Bmad::ele_geometry(
    FloorPositionStruct &floor_start,
    EleStruct &ele,
    std::optional<double> len_scale,
    std::optional<bool> ignore_patch_err
) {
  FloorPositionStruct _floor_end;
  double len_scale_lvalue;
  auto *_len_scale{&len_scale_lvalue};
  if (len_scale.has_value()) {
    len_scale_lvalue = len_scale.value();
  } else {
    _len_scale = nullptr;
  }
  bool ignore_patch_err_lvalue;
  auto *_ignore_patch_err{&ignore_patch_err_lvalue};
  if (ignore_patch_err.has_value()) {
    ignore_patch_err_lvalue = ignore_patch_err.value();
  } else {
    _ignore_patch_err = nullptr;
  }
  fortran_ele_geometry(/* void* */ floor_start.get_fortran_ptr(),
                       /* void* */ ele.get_fortran_ptr(),
                       /* void* */ _floor_end.get_fortran_ptr(),
                       /* double* */ _len_scale,
                       /* bool* */ _ignore_patch_err);
  return std::move(_floor_end);
}
FloorPositionStruct
Bmad::ele_geometry_with_misalignments(EleStruct &ele, std::optional<double> len_scale) {
  double len_scale_lvalue;
  auto *_len_scale{&len_scale_lvalue};
  if (len_scale.has_value()) {
    len_scale_lvalue = len_scale.value();
  } else {
    _len_scale = nullptr;
  }
  FloorPositionStruct _floor;
  fortran_ele_geometry_with_misalignments(/* void* */ ele.get_fortran_ptr(),
                                          /* double* */ _len_scale,
                                          /* void* */ _floor.get_fortran_ptr());
  return std::move(_floor);
}
bool Bmad::ele_has_constant_ds_dt_ref(EleStruct &ele) {
  bool _is_const{};
  fortran_ele_has_constant_ds_dt_ref(/* void* */ ele.get_fortran_ptr(), /* bool& */ _is_const);
  return _is_const;
}
void Bmad::ele_has_nonzero_kick(EleStruct &ele, bool has_kick) {
  fortran_ele_has_nonzero_kick(/* void* */ ele.get_fortran_ptr(), /* bool& */ has_kick);
}
bool Bmad::ele_has_nonzero_offset(EleStruct &ele) {
  bool _has_offset{};
  fortran_ele_has_nonzero_offset(/* void* */ ele.get_fortran_ptr(), /* bool& */ _has_offset);
  return _has_offset;
}
bool Bmad::ele_is_monitor(EleStruct &ele, std::optional<bool> print_warning) {
  bool print_warning_lvalue;
  auto *_print_warning{&print_warning_lvalue};
  if (print_warning.has_value()) {
    print_warning_lvalue = print_warning.value();
  } else {
    _print_warning = nullptr;
  }
  bool _is_monitor{};
  fortran_ele_is_monitor(/* void* */ ele.get_fortran_ptr(),
                         /* bool* */ _print_warning,
                         /* bool& */ _is_monitor);
  return _is_monitor;
}
LatEleLocStruct Bmad::ele_loc(EleStruct &ele) {
  LatEleLocStruct _loc;
  fortran_ele_loc(/* void* */ ele.get_fortran_ptr(), /* void* */ _loc.get_fortran_ptr());
  return std::move(_loc);
}
std::string Bmad::ele_loc_name(
    EleStruct &ele,
    std::optional<bool> show_branch0,
    std::optional<std::string> parens
) {
  bool show_branch0_lvalue;
  auto *_show_branch0{&show_branch0_lvalue};
  if (show_branch0.has_value()) {
    show_branch0_lvalue = show_branch0.value();
  } else {
    _show_branch0 = nullptr;
  }
  const char *_parens = parens.has_value() ? parens->c_str() : nullptr;
  char _str[4096];
  fortran_ele_loc_name(/* void* */ ele.get_fortran_ptr(),
                       /* bool* */ _show_branch0,
                       /* const char* */ _parens,
                       /* const char* */ _str);
  return _str;
}
Bmad::EleMisalignmentLSCalc Bmad::ele_misalignment_l_s_calc(EleStruct &ele) {
  // L_mis: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _L_mis_desc;
  _L_mis_desc.rank = 1;
  FixedArray1D<Real, 3> _L_mis;
  _L_mis_desc.data_ptr = _L_mis.data();
  _L_mis_desc.dims[0] = _L_mis.size();
  // S_mis: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _S_mis_desc;
  _S_mis_desc.rank = 2;
  FixedArray2D<Real, 3, 3> S_mis;
  double _S_mis_vec[3 * 3];
  _S_mis_desc.data_ptr = _S_mis_vec;
  _S_mis_desc.dims[0] = 3;
  _S_mis_desc.dims[1] = 3;
  fortran_ele_misalignment_l_s_calc(/* void* */ ele.get_fortran_ptr(),
                                    /* Bmad::array_descriptor_t& */ _L_mis_desc,
                                    /* Bmad::array_descriptor_t& */ _S_mis_desc);
  vec_to_matrix(_S_mis_vec, S_mis);
  return EleMisalignmentLSCalc{_L_mis, S_mis};
}
int Bmad::ele_nametable_index(EleStruct &ele) {
  int _ix_nt{};
  fortran_ele_nametable_index(/* void* */ ele.get_fortran_ptr(), /* int& */ _ix_nt);
  return _ix_nt;
}
LatEleOrderStruct Bmad::ele_order_calc(LatStruct &lat) {
  LatEleOrderStruct _order;
  fortran_ele_order_calc(/* void* */ lat.get_fortran_ptr(), /* void* */ _order.get_fortran_ptr());
  return std::move(_order);
}
void Bmad::ele_reference_energy_correction(
    EleStruct &ele,
    CoordStruct &orbit,
    int particle_at,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_ele_reference_energy_correction(/* void* */ ele.get_fortran_ptr(),
                                          /* void* */ orbit.get_fortran_ptr(),
                                          /* int& */ particle_at,
                                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                                          /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
int Bmad::ele_rf_step_index(double E_ref, double s_rel, EleStruct &ele) {
  int _ix_step{};
  fortran_ele_rf_step_index(
      /* double& */ E_ref,
      /* double& */ s_rel,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ _ix_step
  );
  return _ix_step;
}
Bmad::EleToFibre Bmad::ele_to_fibre(
    EleStruct &ele,
    bool use_offsets,
    std::optional<int> integ_order,
    std::optional<int> steps,
    std::optional<bool> for_layout,
    optional_ref<CoordStruct> ref_in
) {
  void *_ptc_fibre;
  bool _err_flag{};
  int integ_order_lvalue;
  auto *_integ_order{&integ_order_lvalue};
  if (integ_order.has_value()) {
    integ_order_lvalue = integ_order.value();
  } else {
    _integ_order = nullptr;
  }
  int steps_lvalue;
  auto *_steps{&steps_lvalue};
  if (steps.has_value()) {
    steps_lvalue = steps.value();
  } else {
    _steps = nullptr;
  }
  bool for_layout_lvalue;
  auto *_for_layout{&for_layout_lvalue};
  if (for_layout.has_value()) {
    for_layout_lvalue = for_layout.value();
  } else {
    _for_layout = nullptr;
  }
  auto *_ref_in = ref_in.has_value() ? ref_in->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_ele_to_fibre(/* void* */ ele.get_fortran_ptr(),
                       /* void* */ &_ptc_fibre,
                       /* bool& */ use_offsets,
                       /* bool& */ _err_flag,
                       /* int* */ _integ_order,
                       /* int* */ _steps,
                       /* bool* */ _for_layout,
                       /* void* */ _ref_in);
  return EleToFibre{
      std::move((_ptc_fibre ? std::make_optional<Fibre>(_ptc_fibre) : std::nullopt)),
      _err_flag
  };
}
int Bmad::ele_to_ptc_magnetic_bn_an(EleStruct &ele, FArray1D<Real> &bn, FArray1D<Real> &an) {
  // bn: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _bn_desc;
  _bn_desc.rank = 1;
  _bn_desc.data_ptr = bn.data();
  _bn_desc.dims[0] = bn.size();
  // an: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _an_desc;
  _an_desc.rank = 1;
  _an_desc.data_ptr = an.data();
  _an_desc.dims[0] = an.size();
  int _n_max{};
  fortran_ele_to_ptc_magnetic_bn_an(/* void* */ ele.get_fortran_ptr(),
                                    /* Bmad::array_descriptor_t& */ _bn_desc,
                                    /* Bmad::array_descriptor_t& */ _an_desc,
                                    /* int& */ _n_max);
  return _n_max;
}
void Bmad::ele_to_spin_taylor(EleStruct &ele, LatParamStruct &param, CoordStruct &orb0) {
  fortran_ele_to_spin_taylor(/* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* void* */ orb0.get_fortran_ptr());
}
Bmad::EleToTaylor Bmad::ele_to_taylor(
    EleStruct &ele,
    optional_ref<CoordStruct> orb0,
    std::optional<bool> taylor_map_includes_offsets,
    std::optional<bool> include_damping
) {
  auto *_orb0 = orb0.has_value() ? orb0->get().get_fortran_ptr() : nullptr; // input, optional
  bool taylor_map_includes_offsets_lvalue;
  auto *_taylor_map_includes_offsets{&taylor_map_includes_offsets_lvalue};
  if (taylor_map_includes_offsets.has_value()) {
    taylor_map_includes_offsets_lvalue = taylor_map_includes_offsets.value();
  } else {
    _taylor_map_includes_offsets = nullptr;
  }
  bool include_damping_lvalue;
  auto *_include_damping{&include_damping_lvalue};
  if (include_damping.has_value()) {
    include_damping_lvalue = include_damping.value();
  } else {
    _include_damping = nullptr;
  }
  // orbital_taylor: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbital_taylor_desc;
  _orbital_taylor_desc.rank = 1;
  // Output-only type array
  auto orbital_taylor = TaylorStructAlloc1D(6);
  _orbital_taylor_desc.data_ptr = orbital_taylor.get_fortran_ptr();
  _orbital_taylor_desc.dims[0] = orbital_taylor.size();

  _orbital_taylor_desc.strides[0] = 1;
  // spin_taylor: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_taylor_desc;
  _spin_taylor_desc.rank = 1;
  // Output-only type array
  auto spin_taylor = TaylorStructAlloc1D(4);
  _spin_taylor_desc.data_ptr = spin_taylor.get_fortran_ptr();
  _spin_taylor_desc.dims[0] = spin_taylor.size();

  _spin_taylor_desc.strides[0] = 1;
  fortran_ele_to_taylor(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ _orb0,
                        /* bool* */ _taylor_map_includes_offsets,
                        /* bool* */ _include_damping,
                        /* Bmad::array_descriptor_t& */ _orbital_taylor_desc,
                        /* Bmad::array_descriptor_t& */ _spin_taylor_desc);
  return EleToTaylor{std::move(std::move(orbital_taylor)), std::move(std::move(spin_taylor))};
}
std::string Bmad::ele_unique_name(EleStruct &ele, LatEleOrderStruct &order) {
  char _unique_name[4096];
  fortran_ele_unique_name(/* void* */ ele.get_fortran_ptr(),
                          /* void* */ order.get_fortran_ptr(),
                          /* const char* */ _unique_name);
  return _unique_name;
}
bool Bmad::ele_value_has_changed(
    EleStruct &ele,
    FArray1D<Int> &list,
    FArray1D<Real> &abs_tol,
    bool set_old
) {
  // list: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _list_desc;
  _list_desc.rank = 1;
  _list_desc.data_ptr = list.data();
  _list_desc.dims[0] = list.size();
  // abs_tol: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _abs_tol_desc;
  _abs_tol_desc.rank = 1;
  _abs_tol_desc.data_ptr = abs_tol.data();
  _abs_tol_desc.dims[0] = abs_tol.size();
  bool _has_changed{};
  fortran_ele_value_has_changed(/* void* */ ele.get_fortran_ptr(),
                                /* Bmad::array_descriptor_t& */ _list_desc,
                                /* Bmad::array_descriptor_t& */ _abs_tol_desc,
                                /* bool& */ set_old,
                                /* bool& */ _has_changed);
  return _has_changed;
}
void Bmad::ele_vec_equal_ele_vec(EleStructArray1D ele1, EleStructArray1D ele2) {
  // ele1: EleStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele1_desc;
  _ele1_desc.rank = 1;
  _ele1_desc.data_ptr = ele1.data();
  _ele1_desc.dims[0] = ele1.size();
  _ele1_desc.strides[0] = 1;
  // ele2: EleStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele2_desc;
  _ele2_desc.rank = 1;
  _ele2_desc.data_ptr = ele2.data();
  _ele2_desc.dims[0] = ele2.size();
  _ele2_desc.strides[0] = 1;
  fortran_ele_vec_equal_ele_vec(
      /* Bmad::array_descriptor_t& */ _ele1_desc,
      /* Bmad::array_descriptor_t& */ _ele2_desc
  );
}
Bmad::ElecMultipoleField Bmad::elec_multipole_field(double a, double b, int n, CoordStruct &coord) {
  double _Ex{};
  double _Ey{};
  // dE: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _dE_desc;
  _dE_desc.rank = 2;
  FixedArray2D<Real, 2, 2> dE;
  double _dE_vec[2 * 2];
  _dE_desc.data_ptr = _dE_vec;
  _dE_desc.dims[0] = 2;
  _dE_desc.dims[1] = 2;
  bool _compute_dE{};
  fortran_elec_multipole_field(
      /* double& */ a,
      /* double& */ b,
      /* int& */ n,
      /* void* */ coord.get_fortran_ptr(),
      /* double& */ _Ex,
      /* double& */ _Ey,
      /* Bmad::array_descriptor_t& */ _dE_desc,
      /* bool& */ _compute_dE
  );
  vec_to_matrix(_dE_vec, dE);
  return ElecMultipoleField{_Ex, _Ey, dE, _compute_dE};
}
Bmad::ElementAtSBranch
Bmad::element_at_s(BranchStruct &branch, double s, bool choose_max, std::optional<bool> print_err) {
  bool _err_flag{};
  double _s_eff{};
  CoordStruct _position;
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _ix_ele{};
  fortran_element_at_s_branch(/* void* */ branch.get_fortran_ptr(),
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
    LatStruct &lat,
    double s,
    bool choose_max,
    std::optional<int> ix_branch,
    std::optional<bool> print_err
) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  double _s_eff{};
  CoordStruct _position;
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _ix_ele{};
  fortran_element_at_s_lat(/* void* */ lat.get_fortran_ptr(),
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
    EleStruct &ele,
    LatParamStruct &param,
    int i_slice,
    int n_slice_tot,
    EleStruct &sliced_ele,
    std::optional<double> s_start,
    std::optional<double> s_end
) {
  double s_start_lvalue;
  auto *_s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto *_s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  fortran_element_slice_iterator(/* void* */ ele.get_fortran_ptr(),
                                 /* void* */ param.get_fortran_ptr(),
                                 /* int& */ i_slice,
                                 /* int& */ n_slice_tot,
                                 /* void* */ sliced_ele.get_fortran_ptr(),
                                 /* double* */ _s_start,
                                 /* double* */ _s_end);
}
void Bmad::ellipinc_test() { fortran_ellipinc_test(); }
Bmad::EmFieldCalc Bmad::em_field_calc(
    EleStruct &ele,
    LatParamStruct &param,
    double s_pos,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential,
    std::optional<bool> use_overlap,
    std::optional<bool> grid_allow_s_out_of_bounds,
    std::optional<double> rf_time,
    std::optional<ElePointerStructAlloc1D> used_eles,
    std::optional<bool> print_err,
    optional_ref<EleStruct> original_ele
) {
  EmFieldStruct _field;
  bool calc_dfield_lvalue;
  auto *_calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool _err_flag{};
  bool calc_potential_lvalue;
  auto *_calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  bool use_overlap_lvalue;
  auto *_use_overlap{&use_overlap_lvalue};
  if (use_overlap.has_value()) {
    use_overlap_lvalue = use_overlap.value();
  } else {
    _use_overlap = nullptr;
  }
  bool grid_allow_s_out_of_bounds_lvalue;
  auto *_grid_allow_s_out_of_bounds{&grid_allow_s_out_of_bounds_lvalue};
  if (grid_allow_s_out_of_bounds.has_value()) {
    grid_allow_s_out_of_bounds_lvalue = grid_allow_s_out_of_bounds.value();
  } else {
    _grid_allow_s_out_of_bounds = nullptr;
  }
  double rf_time_lvalue;
  auto *_rf_time{&rf_time_lvalue};
  if (rf_time.has_value()) {
    rf_time_lvalue = rf_time.value();
  } else {
    _rf_time = nullptr;
  }
  // intent=in allocatable type array
  auto *_used_eles =
      used_eles.has_value() ? used_eles->get_fortran_ptr() : nullptr; // input, optional
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  auto *_original_ele =
      original_ele.has_value() ? original_ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_em_field_calc(/* void* */ ele.get_fortran_ptr(),
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
EmFieldStruct Bmad::em_field_derivatives(
    EleStruct &ele,
    LatParamStruct &param,
    double s_pos,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> grid_allow_s_out_of_bounds,
    std::optional<double> rf_time
) {
  EmFieldStruct _dfield;
  bool grid_allow_s_out_of_bounds_lvalue;
  auto *_grid_allow_s_out_of_bounds{&grid_allow_s_out_of_bounds_lvalue};
  if (grid_allow_s_out_of_bounds.has_value()) {
    grid_allow_s_out_of_bounds_lvalue = grid_allow_s_out_of_bounds.value();
  } else {
    _grid_allow_s_out_of_bounds = nullptr;
  }
  double rf_time_lvalue;
  auto *_rf_time{&rf_time_lvalue};
  if (rf_time.has_value()) {
    rf_time_lvalue = rf_time.value();
  } else {
    _rf_time = nullptr;
  }
  fortran_em_field_derivatives(/* void* */ ele.get_fortran_ptr(),
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
    EleStruct &ele,
    LatParamStruct &param,
    double rf_time,
    CoordStruct &orbit,
    bool err_flag,
    std::optional<bool> print_err,
    optional_ref<EmFieldStruct> extra_field
) {
  // dvec_dt: out NOT (CppWrapperGeneralArgumentArray) (['10'])
  Bmad::array_descriptor_t _dvec_dt_desc;
  _dvec_dt_desc.rank = 1;
  FixedArray1D<Real, 10> _dvec_dt;
  _dvec_dt_desc.data_ptr = _dvec_dt.data();
  _dvec_dt_desc.dims[0] = _dvec_dt.size();
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  auto *_extra_field =
      extra_field.has_value() ? extra_field->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_em_field_kick_vector_time(/* void* */ ele.get_fortran_ptr(),
                                    /* void* */ param.get_fortran_ptr(),
                                    /* double& */ rf_time,
                                    /* void* */ orbit.get_fortran_ptr(),
                                    /* Bmad::array_descriptor_t& */ _dvec_dt_desc,
                                    /* bool& */ err_flag,
                                    /* bool* */ _print_err,
                                    /* void* */ _extra_field);
  return _dvec_dt;
}
void Bmad::em_field_plus_em_field(
    EmFieldStruct &field1,
    EmFieldStruct &field2,
    EmFieldStruct &field_tot
) {
  fortran_em_field_plus_em_field(/* void* */ field1.get_fortran_ptr(),
                                 /* void* */ field2.get_fortran_ptr(),
                                 /* void* */ field_tot.get_fortran_ptr());
}
Bmad::Emit6d Bmad::emit_6d(
    EleStruct &ele_ref,
    bool include_opening_angle,
    std::optional<CoordStructArray1D> closed_orbit
) {
  NormalModesStruct _mode;
  // sigma_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> sigma_mat;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  // closed_orbit: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _closed_orbit_desc;
  _closed_orbit_desc.rank = 1;
  if (closed_orbit) {
    _closed_orbit_desc.data_ptr = closed_orbit->data();
    _closed_orbit_desc.dims[0] = closed_orbit->size();
  } else {
    _closed_orbit_desc.data_ptr = nullptr;
    _closed_orbit_desc.dims[0] = 0;
  }
  _closed_orbit_desc.strides[0] = 1;
  RadIntAllEleStruct _rad_int_by_ele;
  fortran_emit_6d(/* void* */ ele_ref.get_fortran_ptr(),
                  /* bool& */ include_opening_angle,
                  /* void* */ _mode.get_fortran_ptr(),
                  /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
                  /* Bmad::array_descriptor_t& */ _closed_orbit_desc,
                  /* void* */ _rad_int_by_ele.get_fortran_ptr());
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  return Emit6d{std::move(_mode), sigma_mat, std::move(_rad_int_by_ele)};
}
bool Bmad::entering_element(CoordStruct &orbit, int particle_at) {
  bool _is_entering{};
  fortran_entering_element(/* void* */ orbit.get_fortran_ptr(),
                           /* int& */ particle_at,
                           /* bool& */ _is_entering);
  return _is_entering;
}
void Bmad::envelope_radints(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    FixedArray1D<Real, 3> alpha,
    FixedArray1D<Real, 3> emit
) {
  // Lambda: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Lambda_desc;
  _Lambda_desc.rank = 2;
  std::complex<double> _Lambda_vec[6 * 6];
  _Lambda_desc.data_ptr = _Lambda_vec;
  _Lambda_desc.dims[0] = 6;
  _Lambda_desc.dims[1] = 6;
  matrix_to_vec(Lambda, _Lambda_vec);
  // Theta: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Theta_desc;
  _Theta_desc.rank = 2;
  std::complex<double> _Theta_vec[6 * 6];
  _Theta_desc.data_ptr = _Theta_vec;
  _Theta_desc.dims[0] = 6;
  _Theta_desc.dims[1] = 6;
  matrix_to_vec(Theta, _Theta_vec);
  // Iota: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Iota_desc;
  _Iota_desc.rank = 2;
  std::complex<double> _Iota_vec[6 * 6];
  _Iota_desc.data_ptr = _Iota_vec;
  _Iota_desc.dims[0] = 6;
  _Iota_desc.dims[1] = 6;
  matrix_to_vec(Iota, _Iota_vec);
  // alpha: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _alpha_desc;
  _alpha_desc.rank = 1;
  _alpha_desc.data_ptr = alpha.data();
  _alpha_desc.dims[0] = alpha.size();
  // emit: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _emit_desc;
  _emit_desc.rank = 1;
  _emit_desc.data_ptr = emit.data();
  _emit_desc.dims[0] = emit.size();
  fortran_envelope_radints(
      /* Bmad::array_descriptor_t& */ _Lambda_desc,
      /* Bmad::array_descriptor_t& */ _Theta_desc,
      /* Bmad::array_descriptor_t& */ _Iota_desc,
      /* Bmad::array_descriptor_t& */ _alpha_desc,
      /* Bmad::array_descriptor_t& */ _emit_desc
  );
  vec_to_matrix(_Lambda_vec, Lambda);
  vec_to_matrix(_Theta_vec, Theta);
  vec_to_matrix(_Iota_vec, Iota);
}
Bmad::EnvelopeRadintsIbs Bmad::envelope_radints_ibs(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    EleStructArray1D eles,
    NormalModesStruct &mode,
    bool tail_cut,
    double npart,
    int species
) {
  // Lambda: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Lambda_desc;
  _Lambda_desc.rank = 2;
  std::complex<double> _Lambda_vec[6 * 6];
  _Lambda_desc.data_ptr = _Lambda_vec;
  _Lambda_desc.dims[0] = 6;
  _Lambda_desc.dims[1] = 6;
  matrix_to_vec(Lambda, _Lambda_vec);
  // Theta: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Theta_desc;
  _Theta_desc.rank = 2;
  std::complex<double> _Theta_vec[6 * 6];
  _Theta_desc.data_ptr = _Theta_vec;
  _Theta_desc.dims[0] = 6;
  _Theta_desc.dims[1] = 6;
  matrix_to_vec(Theta, _Theta_vec);
  // Iota: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Iota_desc;
  _Iota_desc.rank = 2;
  std::complex<double> _Iota_vec[6 * 6];
  _Iota_desc.data_ptr = _Iota_vec;
  _Iota_desc.dims[0] = 6;
  _Iota_desc.dims[1] = 6;
  matrix_to_vec(Iota, _Iota_vec);
  // eles: EleStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _eles_desc;
  _eles_desc.rank = 1;
  _eles_desc.data_ptr = eles.data();
  _eles_desc.dims[0] = eles.size();
  _eles_desc.strides[0] = 1;
  // alpha: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _alpha_desc;
  _alpha_desc.rank = 1;
  FixedArray1D<Real, 3> _alpha;
  _alpha_desc.data_ptr = _alpha.data();
  _alpha_desc.dims[0] = _alpha.size();
  // emit: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _emit_desc;
  _emit_desc.rank = 1;
  FixedArray1D<Real, 3> _emit;
  _emit_desc.data_ptr = _emit.data();
  _emit_desc.dims[0] = _emit.size();
  fortran_envelope_radints_ibs(
      /* Bmad::array_descriptor_t& */ _Lambda_desc,
      /* Bmad::array_descriptor_t& */ _Theta_desc,
      /* Bmad::array_descriptor_t& */ _Iota_desc,
      /* Bmad::array_descriptor_t& */ _eles_desc,
      /* Bmad::array_descriptor_t& */ _alpha_desc,
      /* Bmad::array_descriptor_t& */ _emit_desc,
      /* void* */ mode.get_fortran_ptr(),
      /* bool& */ tail_cut,
      /* double& */ npart,
      /* int& */ species
  );
  return EnvelopeRadintsIbs{_alpha, _emit};
}
void Bmad::eq_ac_kicker(AcKickerStruct &f1, AcKickerStruct &f2, bool is_eq) {
  fortran_eq_ac_kicker(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_ac_kicker_freq(AcKickerFreqStruct &f1, AcKickerFreqStruct &f2, bool is_eq) {
  fortran_eq_ac_kicker_freq(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_ac_kicker_time(AcKickerTimeStruct &f1, AcKickerTimeStruct &f2, bool is_eq) {
  fortran_eq_ac_kicker_time(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_anormal_mode(AnormalModeStruct &f1, AnormalModeStruct &f2, bool is_eq) {
  fortran_eq_anormal_mode(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_aperture_param(ApertureParamStruct &f1, ApertureParamStruct &f2, bool is_eq) {
  fortran_eq_aperture_param(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_aperture_point(AperturePointStruct &f1, AperturePointStruct &f2, bool is_eq) {
  fortran_eq_aperture_point(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_aperture_scan(ApertureScanStruct &f1, ApertureScanStruct &f2, bool is_eq) {
  fortran_eq_aperture_scan(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_beam(BeamStruct &f1, BeamStruct &f2, bool is_eq) {
  fortran_eq_beam(/* void* */ f1.get_fortran_ptr(),
                  /* void* */ f2.get_fortran_ptr(),
                  /* bool& */ is_eq);
}
void Bmad::eq_beam_init(BeamInitStruct &f1, BeamInitStruct &f2, bool is_eq) {
  fortran_eq_beam_init(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_bmad_common(BmadCommonStruct &f1, BmadCommonStruct &f2, bool is_eq) {
  fortran_eq_bmad_common(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_bookkeeping_state(
    BookkeepingStateStruct &f1,
    BookkeepingStateStruct &f2,
    bool is_eq
) {
  fortran_eq_bookkeeping_state(/* void* */ f1.get_fortran_ptr(),
                               /* void* */ f2.get_fortran_ptr(),
                               /* bool& */ is_eq);
}
void Bmad::eq_bpm_phase_coupling(
    BpmPhaseCouplingStruct &f1,
    BpmPhaseCouplingStruct &f2,
    bool is_eq
) {
  fortran_eq_bpm_phase_coupling(/* void* */ f1.get_fortran_ptr(),
                                /* void* */ f2.get_fortran_ptr(),
                                /* bool& */ is_eq);
}
void Bmad::eq_branch(BranchStruct &f1, BranchStruct &f2, bool is_eq) {
  fortran_eq_branch(/* void* */ f1.get_fortran_ptr(),
                    /* void* */ f2.get_fortran_ptr(),
                    /* bool& */ is_eq);
}
void Bmad::eq_bunch(BunchStruct &f1, BunchStruct &f2, bool is_eq) {
  fortran_eq_bunch(/* void* */ f1.get_fortran_ptr(),
                   /* void* */ f2.get_fortran_ptr(),
                   /* bool& */ is_eq);
}
void Bmad::eq_bunch_params(BunchParamsStruct &f1, BunchParamsStruct &f2, bool is_eq) {
  fortran_eq_bunch_params(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map(CartesianMapStruct &f1, CartesianMapStruct &f2, bool is_eq) {
  fortran_eq_cartesian_map(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map_term(
    CartesianMapTermStruct &f1,
    CartesianMapTermStruct &f2,
    bool is_eq
) {
  fortran_eq_cartesian_map_term(/* void* */ f1.get_fortran_ptr(),
                                /* void* */ f2.get_fortran_ptr(),
                                /* bool& */ is_eq);
}
void Bmad::eq_cartesian_map_term1(
    CartesianMapTerm1Struct &f1,
    CartesianMapTerm1Struct &f2,
    bool is_eq
) {
  fortran_eq_cartesian_map_term1(/* void* */ f1.get_fortran_ptr(),
                                 /* void* */ f2.get_fortran_ptr(),
                                 /* bool& */ is_eq);
}
void Bmad::eq_complex_taylor(ComplexTaylorStruct &f1, ComplexTaylorStruct &f2, bool is_eq) {
  fortran_eq_complex_taylor(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_complex_taylor_term(
    ComplexTaylorTermStruct &f1,
    ComplexTaylorTermStruct &f2,
    bool is_eq
) {
  fortran_eq_complex_taylor_term(/* void* */ f1.get_fortran_ptr(),
                                 /* void* */ f2.get_fortran_ptr(),
                                 /* bool& */ is_eq);
}
void Bmad::eq_control(ControlStruct &f1, ControlStruct &f2, bool is_eq) {
  fortran_eq_control(/* void* */ f1.get_fortran_ptr(),
                     /* void* */ f2.get_fortran_ptr(),
                     /* bool& */ is_eq);
}
void Bmad::eq_control_ramp1(ControlRamp1Struct &f1, ControlRamp1Struct &f2, bool is_eq) {
  fortran_eq_control_ramp1(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_control_var1(ControlVar1Struct &f1, ControlVar1Struct &f2, bool is_eq) {
  fortran_eq_control_var1(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_controller(ControllerStruct &f1, ControllerStruct &f2, bool is_eq) {
  fortran_eq_controller(/* void* */ f1.get_fortran_ptr(),
                        /* void* */ f2.get_fortran_ptr(),
                        /* bool& */ is_eq);
}
void Bmad::eq_coord(CoordStruct &f1, CoordStruct &f2, bool is_eq) {
  fortran_eq_coord(/* void* */ f1.get_fortran_ptr(),
                   /* void* */ f2.get_fortran_ptr(),
                   /* bool& */ is_eq);
}
void Bmad::eq_coord_array(CoordArrayStruct &f1, CoordArrayStruct &f2, bool is_eq) {
  fortran_eq_coord_array(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map(CylindricalMapStruct &f1, CylindricalMapStruct &f2, bool is_eq) {
  fortran_eq_cylindrical_map(/* void* */ f1.get_fortran_ptr(),
                             /* void* */ f2.get_fortran_ptr(),
                             /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map_term(
    CylindricalMapTermStruct &f1,
    CylindricalMapTermStruct &f2,
    bool is_eq
) {
  fortran_eq_cylindrical_map_term(/* void* */ f1.get_fortran_ptr(),
                                  /* void* */ f2.get_fortran_ptr(),
                                  /* bool& */ is_eq);
}
void Bmad::eq_cylindrical_map_term1(
    CylindricalMapTerm1Struct &f1,
    CylindricalMapTerm1Struct &f2,
    bool is_eq
) {
  fortran_eq_cylindrical_map_term1(/* void* */ f1.get_fortran_ptr(),
                                   /* void* */ f2.get_fortran_ptr(),
                                   /* bool& */ is_eq);
}
void Bmad::eq_ele(EleStruct &f1, EleStruct &f2, bool is_eq) {
  fortran_eq_ele(/* void* */ f1.get_fortran_ptr(),
                 /* void* */ f2.get_fortran_ptr(),
                 /* bool& */ is_eq);
}
void Bmad::eq_ellipse_beam_init(EllipseBeamInitStruct &f1, EllipseBeamInitStruct &f2, bool is_eq) {
  fortran_eq_ellipse_beam_init(/* void* */ f1.get_fortran_ptr(),
                               /* void* */ f2.get_fortran_ptr(),
                               /* bool& */ is_eq);
}
void Bmad::eq_em_field(EmFieldStruct &f1, EmFieldStruct &f2, bool is_eq) {
  fortran_eq_em_field(/* void* */ f1.get_fortran_ptr(),
                      /* void* */ f2.get_fortran_ptr(),
                      /* bool& */ is_eq);
}
void Bmad::eq_expression_atom(ExpressionAtomStruct &f1, ExpressionAtomStruct &f2, bool is_eq) {
  fortran_eq_expression_atom(/* void* */ f1.get_fortran_ptr(),
                             /* void* */ f2.get_fortran_ptr(),
                             /* bool& */ is_eq);
}
void Bmad::eq_floor_position(FloorPositionStruct &f1, FloorPositionStruct &f2, bool is_eq) {
  fortran_eq_floor_position(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_gen_grad1(GenGrad1Struct &f1, GenGrad1Struct &f2, bool is_eq) {
  fortran_eq_gen_grad1(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_gen_grad_map(GenGradMapStruct &f1, GenGradMapStruct &f2, bool is_eq) {
  fortran_eq_gen_grad_map(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_gg_taylor(GgTaylorStruct &f1, GgTaylorStruct &f2, bool is_eq) {
  fortran_eq_gg_taylor(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_gg_taylor_term(GgTaylorTermStruct &f1, GgTaylorTermStruct &f2, bool is_eq) {
  fortran_eq_gg_taylor_term(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_grid_beam_init(GridBeamInitStruct &f1, GridBeamInitStruct &f2, bool is_eq) {
  fortran_eq_grid_beam_init(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_grid_field(GridFieldStruct &f1, GridFieldStruct &f2, bool is_eq) {
  fortran_eq_grid_field(/* void* */ f1.get_fortran_ptr(),
                        /* void* */ f2.get_fortran_ptr(),
                        /* bool& */ is_eq);
}
void Bmad::eq_grid_field_pt(GridFieldPtStruct &f1, GridFieldPtStruct &f2, bool is_eq) {
  fortran_eq_grid_field_pt(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_grid_field_pt1(GridFieldPt1Struct &f1, GridFieldPt1Struct &f2, bool is_eq) {
  fortran_eq_grid_field_pt1(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_high_energy_space_charge(
    HighEnergySpaceChargeStruct &f1,
    HighEnergySpaceChargeStruct &f2,
    bool is_eq
) {
  fortran_eq_high_energy_space_charge(/* void* */ f1.get_fortran_ptr(),
                                      /* void* */ f2.get_fortran_ptr(),
                                      /* bool& */ is_eq);
}
void Bmad::eq_interval1_coef(Interval1CoefStruct &f1, Interval1CoefStruct &f2, bool is_eq) {
  fortran_eq_interval1_coef(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_kv_beam_init(KvBeamInitStruct &f1, KvBeamInitStruct &f2, bool is_eq) {
  fortran_eq_kv_beam_init(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_lat(LatStruct &f1, LatStruct &f2, bool is_eq) {
  fortran_eq_lat(/* void* */ f1.get_fortran_ptr(),
                 /* void* */ f2.get_fortran_ptr(),
                 /* bool& */ is_eq);
}
void Bmad::eq_lat_ele_loc(LatEleLocStruct &f1, LatEleLocStruct &f2, bool is_eq) {
  fortran_eq_lat_ele_loc(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_lat_param(LatParamStruct &f1, LatParamStruct &f2, bool is_eq) {
  fortran_eq_lat_param(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_linac_normal_mode(LinacNormalModeStruct &f1, LinacNormalModeStruct &f2, bool is_eq) {
  fortran_eq_linac_normal_mode(/* void* */ f1.get_fortran_ptr(),
                               /* void* */ f2.get_fortran_ptr(),
                               /* bool& */ is_eq);
}
void Bmad::eq_mode3(Mode3Struct &f1, Mode3Struct &f2, bool is_eq) {
  fortran_eq_mode3(/* void* */ f1.get_fortran_ptr(),
                   /* void* */ f2.get_fortran_ptr(),
                   /* bool& */ is_eq);
}
void Bmad::eq_mode_info(ModeInfoStruct &f1, ModeInfoStruct &f2, bool is_eq) {
  fortran_eq_mode_info(/* void* */ f1.get_fortran_ptr(),
                       /* void* */ f2.get_fortran_ptr(),
                       /* bool& */ is_eq);
}
void Bmad::eq_normal_modes(NormalModesStruct &f1, NormalModesStruct &f2, bool is_eq) {
  fortran_eq_normal_modes(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_photon_element(PhotonElementStruct &f1, PhotonElementStruct &f2, bool is_eq) {
  fortran_eq_photon_element(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_photon_material(PhotonMaterialStruct &f1, PhotonMaterialStruct &f2, bool is_eq) {
  fortran_eq_photon_material(/* void* */ f1.get_fortran_ptr(),
                             /* void* */ f2.get_fortran_ptr(),
                             /* bool& */ is_eq);
}
void Bmad::eq_photon_reflect_surface(
    PhotonReflectSurfaceStruct &f1,
    PhotonReflectSurfaceStruct &f2,
    bool is_eq
) {
  fortran_eq_photon_reflect_surface(/* void* */ f1.get_fortran_ptr(),
                                    /* void* */ f2.get_fortran_ptr(),
                                    /* bool& */ is_eq);
}
void Bmad::eq_photon_reflect_table(
    PhotonReflectTableStruct &f1,
    PhotonReflectTableStruct &f2,
    bool is_eq
) {
  fortran_eq_photon_reflect_table(/* void* */ f1.get_fortran_ptr(),
                                  /* void* */ f2.get_fortran_ptr(),
                                  /* bool& */ is_eq);
}
void Bmad::eq_photon_target(PhotonTargetStruct &f1, PhotonTargetStruct &f2, bool is_eq) {
  fortran_eq_photon_target(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_pixel_detec(PixelDetecStruct &f1, PixelDetecStruct &f2, bool is_eq) {
  fortran_eq_pixel_detec(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_pixel_pt(PixelPtStruct &f1, PixelPtStruct &f2, bool is_eq) {
  fortran_eq_pixel_pt(/* void* */ f1.get_fortran_ptr(),
                      /* void* */ f2.get_fortran_ptr(),
                      /* bool& */ is_eq);
}
void Bmad::eq_pre_tracker(PreTrackerStruct &f1, PreTrackerStruct &f2, bool is_eq) {
  fortran_eq_pre_tracker(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_rad_int1(RadInt1Struct &f1, RadInt1Struct &f2, bool is_eq) {
  fortran_eq_rad_int1(/* void* */ f1.get_fortran_ptr(),
                      /* void* */ f2.get_fortran_ptr(),
                      /* bool& */ is_eq);
}
void Bmad::eq_rad_int_all_ele(RadIntAllEleStruct &f1, RadIntAllEleStruct &f2, bool is_eq) {
  fortran_eq_rad_int_all_ele(/* void* */ f1.get_fortran_ptr(),
                             /* void* */ f2.get_fortran_ptr(),
                             /* bool& */ is_eq);
}
void Bmad::eq_rad_int_branch(RadIntBranchStruct &f1, RadIntBranchStruct &f2, bool is_eq) {
  fortran_eq_rad_int_branch(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_rad_map(RadMapStruct &f1, RadMapStruct &f2, bool is_eq) {
  fortran_eq_rad_map(/* void* */ f1.get_fortran_ptr(),
                     /* void* */ f2.get_fortran_ptr(),
                     /* bool& */ is_eq);
}
void Bmad::eq_rad_map_ele(RadMapEleStruct &f1, RadMapEleStruct &f2, bool is_eq) {
  fortran_eq_rad_map_ele(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_ramper_lord(RamperLordStruct &f1, RamperLordStruct &f2, bool is_eq) {
  fortran_eq_ramper_lord(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_space_charge_common(
    SpaceChargeCommonStruct &f1,
    SpaceChargeCommonStruct &f2,
    bool is_eq
) {
  fortran_eq_space_charge_common(/* void* */ f1.get_fortran_ptr(),
                                 /* void* */ f2.get_fortran_ptr(),
                                 /* bool& */ is_eq);
}
void Bmad::eq_spin_polar(SpinPolarStruct &f1, SpinPolarStruct &f2, bool is_eq) {
  fortran_eq_spin_polar(/* void* */ f1.get_fortran_ptr(),
                        /* void* */ f2.get_fortran_ptr(),
                        /* bool& */ is_eq);
}
void Bmad::eq_spline(SplineStruct &f1, SplineStruct &f2, bool is_eq) {
  fortran_eq_spline(/* void* */ f1.get_fortran_ptr(),
                    /* void* */ f2.get_fortran_ptr(),
                    /* bool& */ is_eq);
}
void Bmad::eq_strong_beam(StrongBeamStruct &f1, StrongBeamStruct &f2, bool is_eq) {
  fortran_eq_strong_beam(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_surface_curvature(
    SurfaceCurvatureStruct &f1,
    SurfaceCurvatureStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_curvature(/* void* */ f1.get_fortran_ptr(),
                               /* void* */ f2.get_fortran_ptr(),
                               /* bool& */ is_eq);
}
void Bmad::eq_surface_displacement(
    SurfaceDisplacementStruct &f1,
    SurfaceDisplacementStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_displacement(/* void* */ f1.get_fortran_ptr(),
                                  /* void* */ f2.get_fortran_ptr(),
                                  /* bool& */ is_eq);
}
void Bmad::eq_surface_displacement_pt(
    SurfaceDisplacementPtStruct &f1,
    SurfaceDisplacementPtStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_displacement_pt(/* void* */ f1.get_fortran_ptr(),
                                     /* void* */ f2.get_fortran_ptr(),
                                     /* bool& */ is_eq);
}
void Bmad::eq_surface_h_misalign(
    SurfaceHMisalignStruct &f1,
    SurfaceHMisalignStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_h_misalign(/* void* */ f1.get_fortran_ptr(),
                                /* void* */ f2.get_fortran_ptr(),
                                /* bool& */ is_eq);
}
void Bmad::eq_surface_h_misalign_pt(
    SurfaceHMisalignPtStruct &f1,
    SurfaceHMisalignPtStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_h_misalign_pt(/* void* */ f1.get_fortran_ptr(),
                                   /* void* */ f2.get_fortran_ptr(),
                                   /* bool& */ is_eq);
}
void Bmad::eq_surface_segmented(
    SurfaceSegmentedStruct &f1,
    SurfaceSegmentedStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_segmented(/* void* */ f1.get_fortran_ptr(),
                               /* void* */ f2.get_fortran_ptr(),
                               /* bool& */ is_eq);
}
void Bmad::eq_surface_segmented_pt(
    SurfaceSegmentedPtStruct &f1,
    SurfaceSegmentedPtStruct &f2,
    bool is_eq
) {
  fortran_eq_surface_segmented_pt(/* void* */ f1.get_fortran_ptr(),
                                  /* void* */ f2.get_fortran_ptr(),
                                  /* bool& */ is_eq);
}
void Bmad::eq_target_point(TargetPointStruct &f1, TargetPointStruct &f2, bool is_eq) {
  fortran_eq_target_point(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_taylor(TaylorStruct &f1, TaylorStruct &f2, bool is_eq) {
  fortran_eq_taylor(/* void* */ f1.get_fortran_ptr(),
                    /* void* */ f2.get_fortran_ptr(),
                    /* bool& */ is_eq);
}
void Bmad::eq_taylor_term(TaylorTermStruct &f1, TaylorTermStruct &f2, bool is_eq) {
  fortran_eq_taylor_term(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_track(TrackStruct &f1, TrackStruct &f2, bool is_eq) {
  fortran_eq_track(/* void* */ f1.get_fortran_ptr(),
                   /* void* */ f2.get_fortran_ptr(),
                   /* bool& */ is_eq);
}
void Bmad::eq_track_point(TrackPointStruct &f1, TrackPointStruct &f2, bool is_eq) {
  fortran_eq_track_point(/* void* */ f1.get_fortran_ptr(),
                         /* void* */ f2.get_fortran_ptr(),
                         /* bool& */ is_eq);
}
void Bmad::eq_twiss(TwissStruct &f1, TwissStruct &f2, bool is_eq) {
  fortran_eq_twiss(/* void* */ f1.get_fortran_ptr(),
                   /* void* */ f2.get_fortran_ptr(),
                   /* bool& */ is_eq);
}
void Bmad::eq_wake(WakeStruct &f1, WakeStruct &f2, bool is_eq) {
  fortran_eq_wake(/* void* */ f1.get_fortran_ptr(),
                  /* void* */ f2.get_fortran_ptr(),
                  /* bool& */ is_eq);
}
void Bmad::eq_wake_lr(WakeLrStruct &f1, WakeLrStruct &f2, bool is_eq) {
  fortran_eq_wake_lr(/* void* */ f1.get_fortran_ptr(),
                     /* void* */ f2.get_fortran_ptr(),
                     /* bool& */ is_eq);
}
void Bmad::eq_wake_lr_mode(WakeLrModeStruct &f1, WakeLrModeStruct &f2, bool is_eq) {
  fortran_eq_wake_lr_mode(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_wake_sr(WakeSrStruct &f1, WakeSrStruct &f2, bool is_eq) {
  fortran_eq_wake_sr(/* void* */ f1.get_fortran_ptr(),
                     /* void* */ f2.get_fortran_ptr(),
                     /* bool& */ is_eq);
}
void Bmad::eq_wake_sr_mode(WakeSrModeStruct &f1, WakeSrModeStruct &f2, bool is_eq) {
  fortran_eq_wake_sr_mode(/* void* */ f1.get_fortran_ptr(),
                          /* void* */ f2.get_fortran_ptr(),
                          /* bool& */ is_eq);
}
void Bmad::eq_wake_sr_z_long(WakeSrZLongStruct &f1, WakeSrZLongStruct &f2, bool is_eq) {
  fortran_eq_wake_sr_z_long(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_wall3d(Wall3dStruct &f1, Wall3dStruct &f2, bool is_eq) {
  fortran_eq_wall3d(/* void* */ f1.get_fortran_ptr(),
                    /* void* */ f2.get_fortran_ptr(),
                    /* bool& */ is_eq);
}
void Bmad::eq_wall3d_section(Wall3dSectionStruct &f1, Wall3dSectionStruct &f2, bool is_eq) {
  fortran_eq_wall3d_section(/* void* */ f1.get_fortran_ptr(),
                            /* void* */ f2.get_fortran_ptr(),
                            /* bool& */ is_eq);
}
void Bmad::eq_wall3d_vertex(Wall3dVertexStruct &f1, Wall3dVertexStruct &f2, bool is_eq) {
  fortran_eq_wall3d_vertex(/* void* */ f1.get_fortran_ptr(),
                           /* void* */ f2.get_fortran_ptr(),
                           /* bool& */ is_eq);
}
void Bmad::eq_xy_disp(XyDispStruct &f1, XyDispStruct &f2, bool is_eq) {
  fortran_eq_xy_disp(/* void* */ f1.get_fortran_ptr(),
                     /* void* */ f2.get_fortran_ptr(),
                     /* bool& */ is_eq);
}
void Bmad::equal_sign_here(EleStruct &ele, std::string delim, bool is_here) {
  auto _delim = delim.c_str();
  fortran_equal_sign_here(/* void* */ ele.get_fortran_ptr(),
                          /* const char* */ _delim,
                          /* bool& */ is_here);
}
bool Bmad::equivalent_taylor_attributes(EleStruct &ele_taylor, EleStruct &ele2) {
  bool _equiv{};
  fortran_equivalent_taylor_attributes(/* void* */ ele_taylor.get_fortran_ptr(),
                                       /* void* */ ele2.get_fortran_ptr(),
                                       /* bool& */ _equiv);
  return _equiv;
}
void Bmad::etdiv(double A, double B, double C, double D, double E, double F) {
  fortran_etdiv(
      /* double& */ A,
      /* double& */ B,
      /* double& */ C,
      /* double& */ D,
      /* double& */ E,
      /* double& */ F
  );
}
Bmad::EvaluateArrayIndex
Bmad::evaluate_array_index(std::string delim_list1, std::string delim_list2) {
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
      /* int& */ _this_index
  );
  return EvaluateArrayIndex{_err_flag, _word2, _delim2, _this_index};
}
Bmad::EvaluateLogical Bmad::evaluate_logical(std::string word) {
  auto _word = word.c_str();
  int _iostat{};
  bool _this_logic{};
  fortran_evaluate_logical(/* const char* */ _word, /* int& */ _iostat, /* bool& */ _this_logic);
  return EvaluateLogical{_iostat, _this_logic};
}
void Bmad::exact_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_exact_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                               /* void* */ param.get_fortran_ptr(),
                               /* int& */ particle_at,
                               /* void* */ orb.get_fortran_ptr(),
                               /* Bmad::array_descriptor_t& */ _mat6_desc,
                               /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::exp_bessi0(double t, double B1, double B2, double func_retval__) {
  fortran_exp_bessi0(
      /* double& */ t,
      /* double& */ B1,
      /* double& */ B2,
      /* double& */ func_retval__
  );
}
void Bmad::expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string &delim,
    bool delim_found,
    bool is_ok
) {
  auto _delim_list = delim_list.c_str();
  auto _ele_name = ele_name.c_str();
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_expect_one_of(
      /* const char* */ _delim_list,
      /* bool& */ check_input_delim,
      /* const char* */ _ele_name,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ is_ok
  );
}
Bmad::ExpectThis Bmad::expect_this(
    std::string expecting,
    bool check_delim,
    bool call_check,
    std::string err_str,
    EleStruct &ele,
    bool is_ok
) {
  auto _expecting = expecting.c_str();
  auto _err_str = err_str.c_str();
  char _delim[4096];
  bool _delim_found{};
  fortran_expect_this(
      /* const char* */ _expecting,
      /* bool& */ check_delim,
      /* bool& */ call_check,
      /* const char* */ _err_str,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* bool& */ is_ok
  );
  return ExpectThis{_delim, _delim_found};
}
std::string
Bmad::expression_stack_to_string(ExpressionAtomStructArray1D stack, std::optional<bool> polish) {
  // stack: ExpressionAtomStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _stack_desc;
  _stack_desc.rank = 1;
  _stack_desc.data_ptr = stack.data();
  _stack_desc.dims[0] = stack.size();
  _stack_desc.strides[0] = 1;
  bool polish_lvalue;
  auto *_polish{&polish_lvalue};
  if (polish.has_value()) {
    polish_lvalue = polish.value();
  } else {
    _polish = nullptr;
  }
  char _str[4096];
  fortran_expression_stack_to_string(
      /* Bmad::array_descriptor_t& */ _stack_desc,
      /* bool* */ _polish,
      /* const char* */ _str
  );
  return _str;
}
Bmad::ExpressionStackValue Bmad::expression_stack_value(
    ExpressionAtomStructArray1D stack,
    std::optional<ControlVar1StructArray1D> var,
    std::optional<bool> use_old
) {
  // stack: ExpressionAtomStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _stack_desc;
  _stack_desc.rank = 1;
  _stack_desc.data_ptr = stack.data();
  _stack_desc.dims[0] = stack.size();
  _stack_desc.strides[0] = 1;
  bool _err_flag{};
  char _err_str[4096];
  // var: ControlVar1Struct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _var_desc;
  _var_desc.rank = 1;
  if (var) {
    _var_desc.data_ptr = var->data();
    _var_desc.dims[0] = var->size();
  } else {
    _var_desc.data_ptr = nullptr;
    _var_desc.dims[0] = 0;
  }
  _var_desc.strides[0] = 1;
  bool use_old_lvalue;
  auto *_use_old{&use_old_lvalue};
  if (use_old.has_value()) {
    use_old_lvalue = use_old.value();
  } else {
    _use_old = nullptr;
  }
  double _value{};
  fortran_expression_stack_value(
      /* Bmad::array_descriptor_t& */ _stack_desc,
      /* bool& */ _err_flag,
      /* const char* */ _err_str,
      /* Bmad::array_descriptor_t& */ _var_desc,
      /* bool* */ _use_old,
      /* double& */ _value
  );
  return ExpressionStackValue{_err_flag, _err_str, _value};
}
Bmad::ExpressionStringToStack Bmad::expression_string_to_stack(std::string string) {
  auto _string = string.c_str();
  // intent=out allocatable type array
  auto stack{ExpressionAtomStructAlloc1D()};
  int _n_stack{};
  bool _err_flag{};
  char _err_str[4096];
  fortran_expression_string_to_stack(
      /* const char* */ _string,
      /* void* */ stack.get_fortran_ptr(),
      /* int& */ _n_stack,
      /* bool& */ _err_flag,
      /* const char* */ _err_str
  );
  return ExpressionStringToStack{std::move(stack), _n_stack, _err_flag, _err_str};
}
Bmad::ExpressionStringToTree
Bmad::expression_string_to_tree(std::string string, ExpressionTreeStruct &root_tree) {
  auto _string = string.c_str();
  bool _err_flag{};
  char _err_str[4096];
  fortran_expression_string_to_tree(
      /* const char* */ _string,
      /* void* */ root_tree.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* const char* */ _err_str
  );
  return ExpressionStringToTree{_err_flag, _err_str};
}
std::string Bmad::expression_tree_to_string(
    ExpressionTreeStruct &tree,
    std::optional<bool> include_root,
    std::optional<int> n_node,
    optional_ref<ExpressionTreeStruct> parent
) {
  bool include_root_lvalue;
  auto *_include_root{&include_root_lvalue};
  if (include_root.has_value()) {
    include_root_lvalue = include_root.value();
  } else {
    _include_root = nullptr;
  }
  int n_node_lvalue;
  auto *_n_node{&n_node_lvalue};
  if (n_node.has_value()) {
    n_node_lvalue = n_node.value();
  } else {
    _n_node = nullptr;
  }
  auto *_parent = parent.has_value() ? parent->get().get_fortran_ptr() : nullptr; // input, optional
  char _str_out[4096];
  fortran_expression_tree_to_string(/* void* */ tree.get_fortran_ptr(),
                                    /* bool* */ _include_root,
                                    /* int* */ _n_node,
                                    /* void* */ _parent,
                                    /* const char* */ _str_out);
  return _str_out;
}
Bmad::ExpressionValue Bmad::expression_value(
    std::string expression,
    std::optional<ControlVar1StructArray1D> var,
    std::optional<bool> use_old
) {
  auto _expression = expression.c_str();
  bool _err_flag{};
  char _err_str[4096];
  // var: ControlVar1Struct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _var_desc;
  _var_desc.rank = 1;
  if (var) {
    _var_desc.data_ptr = var->data();
    _var_desc.dims[0] = var->size();
  } else {
    _var_desc.data_ptr = nullptr;
    _var_desc.dims[0] = 0;
  }
  _var_desc.strides[0] = 1;
  bool use_old_lvalue;
  auto *_use_old{&use_old_lvalue};
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
      /* Bmad::array_descriptor_t& */ _var_desc,
      /* bool* */ _use_old,
      /* double& */ _value
  );
  return ExpressionValue{_err_flag, _err_str, _value};
}
int Bmad::fft1(FArray1D<Real> &a, FArray1D<Real> &b, int n, int isn) {
  // a: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _a_desc;
  _a_desc.rank = 1;
  _a_desc.data_ptr = a.data();
  _a_desc.dims[0] = a.size();
  // b: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  _b_desc.data_ptr = b.data();
  _b_desc.dims[0] = b.size();
  int _ierr{};
  fortran_fft1(
      /* Bmad::array_descriptor_t& */ _a_desc,
      /* Bmad::array_descriptor_t& */ _b_desc,
      /* int& */ n,
      /* int& */ isn,
      /* int& */ _ierr
  );
  return _ierr;
}
bool Bmad::fibre_to_ele(
    Fibre &ptc_fibre,
    BranchStruct &branch,
    int &ix_ele,
    std::optional<bool> from_mad
) {
  bool _err_flag{};
  bool from_mad_lvalue;
  auto *_from_mad{&from_mad_lvalue};
  if (from_mad.has_value()) {
    from_mad_lvalue = from_mad.value();
  } else {
    _from_mad = nullptr;
  }
  fortran_fibre_to_ele(/* void* */ ptc_fibre.get_fortran_ptr(),
                       /* void* */ branch.get_fortran_ptr(),
                       /* int& */ ix_ele,
                       /* bool& */ _err_flag,
                       /* bool* */ _from_mad);
  return _err_flag;
}
bool Bmad::field_attribute_free(EleStruct &ele, std::string attrib_name) {
  auto _attrib_name = attrib_name.c_str();
  bool _free{};
  fortran_field_attribute_free(/* void* */ ele.get_fortran_ptr(),
                               /* const char* */ _attrib_name,
                               /* bool& */ _free);
  return _free;
}
void Bmad::finalize_reflectivity_table(PhotonReflectTableStruct &table, bool in_degrees) {
  fortran_finalize_reflectivity_table(/* void* */ table.get_fortran_ptr(), /* bool& */ in_degrees);
}
Bmad::FindElementEnds Bmad::find_element_ends(EleStruct &ele, std::optional<int> ix_multipass) {
  void *_ele1;
  void *_ele2;
  int ix_multipass_lvalue;
  auto *_ix_multipass{&ix_multipass_lvalue};
  if (ix_multipass.has_value()) {
    ix_multipass_lvalue = ix_multipass.value();
  } else {
    _ix_multipass = nullptr;
  }
  fortran_find_element_ends(/* void* */ ele.get_fortran_ptr(),
                            /* void* */ &_ele1,
                            /* void* */ &_ele2,
                            /* int* */ _ix_multipass);
  return FindElementEnds{
      std::move((_ele1 ? std::make_optional<EleStruct>(_ele1) : std::nullopt)),
      std::move((_ele2 ? std::make_optional<EleStruct>(_ele2) : std::nullopt))
  };
}
double Bmad::find_fwhm(double bound, FixedArray1D<Real, 8> args) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  double _fwhm{};
  fortran_find_fwhm(
      /* double& */ bound,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ _fwhm
  );
  return _fwhm;
}
Bmad::FindMatchingFieldmap Bmad::find_matching_fieldmap(
    std::string file_name,
    EleStruct &ele,
    int fm_type,
    std::optional<bool> ignore_slaves
) {
  auto _file_name = file_name.c_str();
  void *_match_ele;
  int _ix_field{};
  bool ignore_slaves_lvalue;
  auto *_ignore_slaves{&ignore_slaves_lvalue};
  if (ignore_slaves.has_value()) {
    ignore_slaves_lvalue = ignore_slaves.value();
  } else {
    _ignore_slaves = nullptr;
  }
  fortran_find_matching_fieldmap(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ fm_type,
      /* void* */ &_match_ele,
      /* int& */ _ix_field,
      /* bool* */ _ignore_slaves
  );
  return FindMatchingFieldmap{
      std::move((_match_ele ? std::make_optional<EleStruct>(_match_ele) : std::nullopt)),
      _ix_field
  };
}
double Bmad::find_normalization(double bound, double p0, FixedArray1D<Real, 8> args) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  double _pnrml{};
  fortran_find_normalization(
      /* double& */ bound,
      /* double& */ p0,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ _pnrml
  );
  return _pnrml;
}
Bmad::FloorAnglesToWMat Bmad::floor_angles_to_w_mat(double theta, double phi, double psi) {
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  // w_mat_inv: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_inv_desc;
  _w_mat_inv_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat_inv;
  double _w_mat_inv_vec[3 * 3];
  _w_mat_inv_desc.data_ptr = _w_mat_inv_vec;
  _w_mat_inv_desc.dims[0] = 3;
  _w_mat_inv_desc.dims[1] = 3;
  fortran_floor_angles_to_w_mat(
      /* double& */ theta,
      /* double& */ phi,
      /* double& */ psi,
      /* Bmad::array_descriptor_t& */ _w_mat_desc,
      /* Bmad::array_descriptor_t& */ _w_mat_inv_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  vec_to_matrix(_w_mat_inv_vec, w_mat_inv);
  return FloorAnglesToWMat{w_mat, w_mat_inv};
}
Bmad::FloorWMatToAngles Bmad::floor_w_mat_to_angles(
    FixedArray2D<Real, 3, 3> w_mat,
    optional_ref<FloorPositionStruct> floor0
) {
  // w_mat: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  matrix_to_vec(w_mat, _w_mat_vec);
  double _theta{};
  double _phi{};
  double _psi{};
  auto *_floor0 = floor0.has_value() ? floor0->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_floor_w_mat_to_angles(
      /* Bmad::array_descriptor_t& */ _w_mat_desc,
      /* double& */ _theta,
      /* double& */ _phi,
      /* double& */ _psi,
      /* void* */ _floor0
  );
  return FloorWMatToAngles{_theta, _phi, _psi};
}
ComplexTaylorStruct Bmad::form_complex_taylor(TaylorStruct &re_taylor, TaylorStruct &im_taylor) {
  ComplexTaylorStruct _complex_taylor;
  fortran_form_complex_taylor(/* void* */ re_taylor.get_fortran_ptr(),
                              /* void* */ im_taylor.get_fortran_ptr(),
                              /* void* */ _complex_taylor.get_fortran_ptr());
  return std::move(_complex_taylor);
}
Bmad::FormDigestedBmadFileName
Bmad::form_digested_bmad_file_name(std::string lat_file, std::optional<std::string> use_line) {
  auto _lat_file = lat_file.c_str();
  char _digested_file[4096];
  char _full_lat_file[4096];
  const char *_use_line = use_line.has_value() ? use_line->c_str() : nullptr;
  fortran_form_digested_bmad_file_name(
      /* const char* */ _lat_file,
      /* const char* */ _digested_file,
      /* const char* */ _full_lat_file,
      /* const char* */ _use_line
  );
  return FormDigestedBmadFileName{_digested_file, _full_lat_file};
}
bool Bmad::fringe_here(EleStruct &ele, CoordStruct &orbit, int particle_at) {
  bool _is_here{};
  fortran_fringe_here(/* void* */ ele.get_fortran_ptr(),
                      /* void* */ orbit.get_fortran_ptr(),
                      /* int& */ particle_at,
                      /* bool& */ _is_here);
  return _is_here;
}
FixedArray1D<Real, 3>
Bmad::g_bend_from_em_field(FixedArray1D<Real, 3> b, FixedArray1D<Real, 3> e, CoordStruct &orbit) {
  // b: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  _b_desc.data_ptr = b.data();
  _b_desc.dims[0] = b.size();
  // e: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _e_desc;
  _e_desc.rank = 1;
  _e_desc.data_ptr = e.data();
  _e_desc.dims[0] = e.size();
  // g_bend: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _g_bend_desc;
  _g_bend_desc.rank = 1;
  FixedArray1D<Real, 3> _g_bend;
  _g_bend_desc.data_ptr = _g_bend.data();
  _g_bend_desc.dims[0] = _g_bend.size();
  fortran_g_bend_from_em_field(
      /* Bmad::array_descriptor_t& */ _b_desc,
      /* Bmad::array_descriptor_t& */ _e_desc,
      /* void* */ orbit.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _g_bend_desc
  );
  return _g_bend;
}
Bmad::GBendingStrengthFromEmField Bmad::g_bending_strength_from_em_field(
    EleStruct &ele,
    LatParamStruct &param,
    double s_rel,
    CoordStruct &orbit,
    bool local_ref_frame
) {
  // g: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _g_desc;
  _g_desc.rank = 1;
  FixedArray1D<Real, 3> _g;
  _g_desc.data_ptr = _g.data();
  _g_desc.dims[0] = _g.size();
  // dg: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _dg_desc;
  _dg_desc.rank = 2;
  FixedArray2D<Real, 3, 3> dg;
  double _dg_vec[3 * 3];
  _dg_desc.data_ptr = _dg_vec;
  _dg_desc.dims[0] = 3;
  _dg_desc.dims[1] = 3;
  fortran_g_bending_strength_from_em_field(/* void* */ ele.get_fortran_ptr(),
                                           /* void* */ param.get_fortran_ptr(),
                                           /* double& */ s_rel,
                                           /* void* */ orbit.get_fortran_ptr(),
                                           /* bool& */ local_ref_frame,
                                           /* Bmad::array_descriptor_t& */ _g_desc,
                                           /* Bmad::array_descriptor_t& */ _dg_desc);
  vec_to_matrix(_dg_vec, dg);
  return GBendingStrengthFromEmField{_g, dg};
}
void Bmad::g_integrals_calc(LatStruct &lat) {
  fortran_g_integrals_calc(/* void* */ lat.get_fortran_ptr());
}
double Bmad::gamma_ref(EleStruct &ele) {
  double _gamma{};
  fortran_gamma_ref(/* void* */ ele.get_fortran_ptr(), /* double& */ _gamma);
  return _gamma;
}
GgTaylorStructArray1D
Bmad::gen_grad1_to_gg_taylor(EleStruct &ele, GenGradMapStruct &gen_grad, int iz) {
  // gg_taylor: GgTaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _gg_taylor_desc;
  _gg_taylor_desc.rank = 1;
  // Output-only type array
  auto gg_taylor = GgTaylorStructAlloc1D(3);
  _gg_taylor_desc.data_ptr = gg_taylor.get_fortran_ptr();
  _gg_taylor_desc.dims[0] = gg_taylor.size();

  _gg_taylor_desc.strides[0] = 1;
  fortran_gen_grad1_to_gg_taylor(/* void* */ ele.get_fortran_ptr(),
                                 /* void* */ gen_grad.get_fortran_ptr(),
                                 /* int& */ iz,
                                 /* Bmad::array_descriptor_t& */ _gg_taylor_desc);
  return std::move(std::move(gg_taylor));
}
GgTaylorStructArray1D
Bmad::gen_grad_at_s_to_gg_taylor(EleStruct &ele, GenGradMapStruct &gen_grad, double s_pos) {
  // gg_taylor: GgTaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _gg_taylor_desc;
  _gg_taylor_desc.rank = 1;
  // Output-only type array
  auto gg_taylor = GgTaylorStructAlloc1D(3);
  _gg_taylor_desc.data_ptr = gg_taylor.get_fortran_ptr();
  _gg_taylor_desc.dims[0] = gg_taylor.size();

  _gg_taylor_desc.strides[0] = 1;
  fortran_gen_grad_at_s_to_gg_taylor(/* void* */ ele.get_fortran_ptr(),
                                     /* void* */ gen_grad.get_fortran_ptr(),
                                     /* double& */ s_pos,
                                     /* Bmad::array_descriptor_t& */ _gg_taylor_desc);
  return std::move(std::move(gg_taylor));
}
void Bmad::gen_grad_field(
    FArray1D<Real> &deriv,
    GenGrad1Struct &gg,
    double rho,
    double theta,
    FixedArray1D<Real, 3> field
) {
  // deriv: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _deriv_desc;
  _deriv_desc.rank = 1;
  _deriv_desc.data_ptr = deriv.data();
  _deriv_desc.dims[0] = deriv.size();
  // field: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _field_desc;
  _field_desc.rank = 1;
  _field_desc.data_ptr = field.data();
  _field_desc.dims[0] = field.size();
  fortran_gen_grad_field(
      /* Bmad::array_descriptor_t& */ _deriv_desc,
      /* void* */ gg.get_fortran_ptr(),
      /* double& */ rho,
      /* double& */ theta,
      /* Bmad::array_descriptor_t& */ _field_desc
  );
}
double Bmad::get_bl_from_fwhm(double bound, FixedArray1D<Real, 8> args) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  double _sigma{};
  fortran_get_bl_from_fwhm(
      /* double& */ bound,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ _sigma
  );
  return _sigma;
}
void Bmad::get_called_file(std::string delim, std::string call_file, bool err) {
  auto _delim = delim.c_str();
  auto _call_file = call_file.c_str();
  fortran_get_called_file(/* const char* */ _delim, /* const char* */ _call_file, /* bool& */ err);
}
Bmad::GetEmitFromSigmaMat Bmad::get_emit_from_sigma_mat(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<FixedArray2D<Real, 6, 6>> Nmat
) {
  // sigma_mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  // normal: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _normal_desc;
  _normal_desc.rank = 1;
  FixedArray1D<Real, 3> _normal;
  _normal_desc.data_ptr = _normal.data();
  _normal_desc.dims[0] = _normal.size();
  // Nmat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Nmat_desc;
  _Nmat_desc.rank = 2;
  double _Nmat_vec[6 * 6];
  _Nmat_desc.data_ptr = _Nmat_vec;
  if (Nmat.has_value()) {
    matrix_to_vec(Nmat.value(), _Nmat_vec);
    _Nmat_desc.dims[0] = 6;
    _Nmat_desc.dims[1] = 6;
  } else {
    _Nmat_desc.data_ptr = nullptr;
  }
  bool _err_flag{};
  fortran_get_emit_from_sigma_mat(
      /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
      /* Bmad::array_descriptor_t& */ _normal_desc,
      /* Bmad::array_descriptor_t& */ _Nmat_desc,
      /* bool& */ _err_flag
  );
  return GetEmitFromSigmaMat{_normal, _err_flag};
}
void Bmad::get_list_of_names(
    EleStruct &ele,
    std::string err_str,
    CharacterAlloc1D &name_list,
    std::string delim,
    bool delim_found,
    bool err_flag
) {
  auto _err_str = err_str.c_str();
  // intent=inout character array container
  auto _delim = delim.c_str();
  fortran_get_list_of_names(/* void* */ ele.get_fortran_ptr(),
                            /* const char* */ _err_str,
                            /* void* */ name_list.get_fortran_ptr(),
                            /* const char* */ _delim,
                            /* bool& */ delim_found,
                            /* bool& */ err_flag);
}
Bmad::GetNextWord Bmad::get_next_word(
    std::string word,
    std::string delim_list,
    std::optional<bool> upper_case_word,
    std::optional<bool> call_check
) {
  auto _word = word.c_str();
  int _ix_word{};
  auto _delim_list = delim_list.c_str();
  char _delim[4096];
  bool _delim_found{};
  bool upper_case_word_lvalue;
  auto *_upper_case_word{&upper_case_word_lvalue};
  if (upper_case_word.has_value()) {
    upper_case_word_lvalue = upper_case_word.value();
  } else {
    _upper_case_word = nullptr;
  }
  bool call_check_lvalue;
  auto *_call_check{&call_check_lvalue};
  if (call_check.has_value()) {
    call_check_lvalue = call_check.value();
  } else {
    _call_check = nullptr;
  }
  bool _err_flag{};
  fortran_get_next_word(
      /* const char* */ _word,
      /* int& */ _ix_word,
      /* const char* */ _delim_list,
      /* const char* */ _delim,
      /* bool& */ _delim_found,
      /* bool* */ _upper_case_word,
      /* bool* */ _call_check,
      /* bool& */ _err_flag
  );
  return GetNextWord{_ix_word, _delim, _delim_found, _err_flag};
}
void Bmad::get_sequence_args(
    std::string seq_name,
    CharacterAlloc1D &arg_list,
    std::string delim,
    bool err_flag
) {
  auto _seq_name = seq_name.c_str();
  // intent=inout character array container
  auto _delim = delim.c_str();
  fortran_get_sequence_args(
      /* const char* */ _seq_name,
      /* void* */ arg_list.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ err_flag
  );
}
Bmad::GetSlaveList Bmad::get_slave_list(EleStruct &lord) {
  // intent=out allocatable type array
  auto slaves{ElePointerStructAlloc1D()};
  int _n_slave{};
  fortran_get_slave_list(/* void* */ lord.get_fortran_ptr(),
                         /* void* */ slaves.get_fortran_ptr(),
                         /* int& */ _n_slave);
  return GetSlaveList{std::move(slaves), _n_slave};
}
void Bmad::get_switch(
    std::string name,
    CharacterAlloc1D &name_list,
    int switch_,
    bool err,
    EleStruct &ele,
    std::string delim,
    bool delim_found
) {
  auto _name = name.c_str();
  // intent=inout character array container
  auto _delim = delim.c_str();
  fortran_get_switch(
      /* const char* */ _name,
      /* void* */ name_list.get_fortran_ptr(),
      /* int& */ switch_,
      /* bool& */ err,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* bool& */ delim_found
  );
}
void Bmad::gg_taylor_equal_gg_taylor(GgTaylorStruct &gg_taylor1, GgTaylorStruct &gg_taylor2) {
  fortran_gg_taylor_equal_gg_taylor(/* void* */ gg_taylor1.get_fortran_ptr(),
                                    /* void* */ gg_taylor2.get_fortran_ptr());
}
void Bmad::gg_taylors_equal_gg_taylors(
    GgTaylorStructArray1D gg_taylor1,
    GgTaylorStructArray1D gg_taylor2
) {
  // gg_taylor1: GgTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _gg_taylor1_desc;
  _gg_taylor1_desc.rank = 1;
  _gg_taylor1_desc.data_ptr = gg_taylor1.data();
  _gg_taylor1_desc.dims[0] = gg_taylor1.size();
  _gg_taylor1_desc.strides[0] = 1;
  // gg_taylor2: GgTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _gg_taylor2_desc;
  _gg_taylor2_desc.rank = 1;
  _gg_taylor2_desc.data_ptr = gg_taylor2.data();
  _gg_taylor2_desc.dims[0] = gg_taylor2.size();
  _gg_taylor2_desc.strides[0] = 1;
  fortran_gg_taylors_equal_gg_taylors(
      /* Bmad::array_descriptor_t& */ _gg_taylor1_desc,
      /* Bmad::array_descriptor_t& */ _gg_taylor2_desc
  );
}
void Bmad::gpt_field_grid_scaling(
    EleStruct &ele,
    int dimensions,
    double field_scale,
    double ref_time
) {
  fortran_gpt_field_grid_scaling(/* void* */ ele.get_fortran_ptr(),
                                 /* int& */ dimensions,
                                 /* double& */ field_scale,
                                 /* double& */ ref_time);
}
void Bmad::gpt_max_field_reference(GridFieldPt1Struct &pt0, EleStruct &ele, double field_value) {
  fortran_gpt_max_field_reference(/* void* */ pt0.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* double& */ field_value);
}
Bmad::GptToParticleBunch Bmad::gpt_to_particle_bunch(std::string gpt_file, EleStruct &ele) {
  auto _gpt_file = gpt_file.c_str();
  BunchStruct _bunch;
  bool _err_flag{};
  fortran_gpt_to_particle_bunch(
      /* const char* */ _gpt_file,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _bunch.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return GptToParticleBunch{std::move(_bunch), _err_flag};
}
double Bmad::gradient_shift_sr_wake(EleStruct &ele, LatParamStruct &param) {
  double _grad_shift{};
  fortran_gradient_shift_sr_wake(/* void* */ ele.get_fortran_ptr(),
                                 /* void* */ param.get_fortran_ptr(),
                                 /* double& */ _grad_shift);
  return _grad_shift;
}
GridFieldPt1Struct Bmad::grid_field_interpolate(
    EleStruct &ele,
    CoordStruct &orbit,
    GridFieldStruct &grid,
    bool err_flag,
    double x1,
    std::optional<double> x2,
    std::optional<double> x3,
    std::optional<bool> allow_s_out_of_bounds,
    std::optional<bool> print_err
) {
  GridFieldPt1Struct _g_field;
  double x2_lvalue;
  auto *_x2{&x2_lvalue};
  if (x2.has_value()) {
    x2_lvalue = x2.value();
  } else {
    _x2 = nullptr;
  }
  double x3_lvalue;
  auto *_x3{&x3_lvalue};
  if (x3.has_value()) {
    x3_lvalue = x3.value();
  } else {
    _x3 = nullptr;
  }
  bool allow_s_out_of_bounds_lvalue;
  auto *_allow_s_out_of_bounds{&allow_s_out_of_bounds_lvalue};
  if (allow_s_out_of_bounds.has_value()) {
    allow_s_out_of_bounds_lvalue = allow_s_out_of_bounds.value();
  } else {
    _allow_s_out_of_bounds = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_grid_field_interpolate(/* void* */ ele.get_fortran_ptr(),
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
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_hard_multipole_edge_kick(/* void* */ ele.get_fortran_ptr(),
                                   /* void* */ param.get_fortran_ptr(),
                                   /* int& */ particle_at,
                                   /* void* */ orbit.get_fortran_ptr(),
                                   /* Bmad::array_descriptor_t& */ _mat6_desc,
                                   /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::has_attribute(EleStruct &ele, std::string attrib, bool has_it) {
  auto _attrib = attrib.c_str();
  fortran_has_attribute(/* void* */ ele.get_fortran_ptr(),
                        /* const char* */ _attrib,
                        /* bool& */ has_it);
}
bool Bmad::has_curvature(PhotonElementStruct &phot_ele) {
  bool _curved{};
  fortran_has_curvature(/* void* */ phot_ele.get_fortran_ptr(), /* bool& */ _curved);
  return _curved;
}
bool Bmad::has_orientation_attributes(EleStruct &ele) {
  bool _has_attribs{};
  fortran_has_orientation_attributes(/* void* */ ele.get_fortran_ptr(), /* bool& */ _has_attribs);
  return _has_attribs;
}
void Bmad::hdf5_write_beam(
    std::string file_name,
    BunchStructArray1D bunches,
    bool append,
    bool error,
    optional_ref<LatStruct> lat,
    std::optional<bool> alive_only
) {
  auto _file_name = file_name.c_str();
  // bunches: BunchStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _bunches_desc;
  _bunches_desc.rank = 1;
  _bunches_desc.data_ptr = bunches.data();
  _bunches_desc.dims[0] = bunches.size();
  _bunches_desc.strides[0] = 1;
  auto *_lat = lat.has_value() ? lat->get().get_fortran_ptr() : nullptr; // input, optional
  bool alive_only_lvalue;
  auto *_alive_only{&alive_only_lvalue};
  if (alive_only.has_value()) {
    alive_only_lvalue = alive_only.value();
  } else {
    _alive_only = nullptr;
  }
  fortran_hdf5_write_beam(
      /* const char* */ _file_name,
      /* Bmad::array_descriptor_t& */ _bunches_desc,
      /* bool& */ append,
      /* bool& */ error,
      /* void* */ _lat,
      /* bool* */ _alive_only
  );
}
void Bmad::hdf5_write_grid_field(
    std::string file_name,
    EleStruct &ele,
    GridFieldStructArray1D g_field,
    bool err_flag
) {
  auto _file_name = file_name.c_str();
  // g_field: GridFieldStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _g_field_desc;
  _g_field_desc.rank = 1;
  _g_field_desc.data_ptr = g_field.data();
  _g_field_desc.dims[0] = g_field.size();
  _g_field_desc.strides[0] = 1;
  fortran_hdf5_write_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _g_field_desc,
      /* bool& */ err_flag
  );
}
void Bmad::hwang_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_hwang_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                               /* void* */ param.get_fortran_ptr(),
                               /* int& */ particle_at,
                               /* void* */ orb.get_fortran_ptr(),
                               /* Bmad::array_descriptor_t& */ _mat6_desc,
                               /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species,
    FixedArray2D<Real, 6, 6> ibs_mat
) {
  // sigma_mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  // ibs_mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _ibs_mat_desc;
  _ibs_mat_desc.rank = 2;
  double _ibs_mat_vec[6 * 6];
  _ibs_mat_desc.data_ptr = _ibs_mat_vec;
  _ibs_mat_desc.dims[0] = 6;
  _ibs_mat_desc.dims[1] = 6;
  matrix_to_vec(ibs_mat, _ibs_mat_vec);
  fortran_ibs_matrix_c(
      /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
      /* bool& */ tail_cut,
      /* double& */ tau,
      /* double& */ energy,
      /* double& */ n_part,
      /* int& */ species,
      /* Bmad::array_descriptor_t& */ _ibs_mat_desc
  );
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  vec_to_matrix(_ibs_mat_vec, ibs_mat);
}
void Bmad::igfcoulombfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
) {
  fortran_igfcoulombfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res
  );
}
void Bmad::igfexfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
) {
  fortran_igfexfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res
  );
}
void Bmad::igfeyfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
) {
  fortran_igfeyfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res
  );
}
void Bmad::igfezfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
) {
  fortran_igfezfun(
      /* double& */ u,
      /* double& */ v,
      /* double& */ w,
      /* double& */ gam,
      /* double& */ dx,
      /* double& */ dy,
      /* double& */ dz,
      /* double& */ res
  );
}
void Bmad::init_attribute_name1(
    bool &is_ok,
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state,
    std::optional<bool> override
) {
  auto _name = name.c_str();
  int attrib_state_lvalue;
  auto *_attrib_state{&attrib_state_lvalue};
  if (attrib_state.has_value()) {
    attrib_state_lvalue = attrib_state.value();
  } else {
    _attrib_state = nullptr;
  }
  bool override_lvalue;
  auto *_override{&override_lvalue};
  if (override.has_value()) {
    override_lvalue = override.value();
  } else {
    _override = nullptr;
  }
  fortran_init_attribute_name1(
      /* bool& */ is_ok,
      /* int& */ ix_key,
      /* int& */ ix_attrib,
      /* const char* */ _name,
      /* int* */ _attrib_state,
      /* bool* */ _override
  );
}
void Bmad::init_attribute_name_array() { fortran_init_attribute_name_array(); }
Bmad::InitBeamDistribution Bmad::init_beam_distribution(
    EleStruct &ele,
    LatParamStruct &param,
    BeamInitStruct &beam_init,
    optional_ref<NormalModesStruct> modes,
    std::optional<bool> print_p0c_shift_warning,
    std::optional<bool> conserve_momentum
) {
  BeamStruct _beam;
  bool _err_flag{};
  auto *_modes = modes.has_value() ? modes->get().get_fortran_ptr() : nullptr; // input, optional
  BeamInitStruct _beam_init_set;
  bool print_p0c_shift_warning_lvalue;
  auto *_print_p0c_shift_warning{&print_p0c_shift_warning_lvalue};
  if (print_p0c_shift_warning.has_value()) {
    print_p0c_shift_warning_lvalue = print_p0c_shift_warning.value();
  } else {
    _print_p0c_shift_warning = nullptr;
  }
  bool conserve_momentum_lvalue;
  auto *_conserve_momentum{&conserve_momentum_lvalue};
  if (conserve_momentum.has_value()) {
    conserve_momentum_lvalue = conserve_momentum.value();
  } else {
    _conserve_momentum = nullptr;
  }
  fortran_init_beam_distribution(/* void* */ ele.get_fortran_ptr(),
                                 /* void* */ param.get_fortran_ptr(),
                                 /* void* */ beam_init.get_fortran_ptr(),
                                 /* void* */ _beam.get_fortran_ptr(),
                                 /* bool& */ _err_flag,
                                 /* void* */ _modes,
                                 /* void* */ _beam_init_set.get_fortran_ptr(),
                                 /* bool* */ _print_p0c_shift_warning,
                                 /* bool* */ _conserve_momentum);
  return InitBeamDistribution{std::move(_beam), _err_flag, std::move(_beam_init_set)};
}
void Bmad::init_bmad() { fortran_init_bmad(); }
void Bmad::init_bmad_parser_common(optional_ref<LatStruct> lat) {
  auto *_lat = lat.has_value() ? lat->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_init_bmad_parser_common(/* void* */ _lat);
}
Bmad::InitBunchDistribution Bmad::init_bunch_distribution(
    EleStruct &ele,
    LatParamStruct &param,
    BeamInitStruct &beam_init,
    int ix_bunch,
    optional_ref<NormalModesStruct> modes,
    std::optional<bool> print_p0c_shift_warning,
    std::optional<bool> conserve_momentum
) {
  BunchStruct _bunch;
  bool _err_flag{};
  auto *_modes = modes.has_value() ? modes->get().get_fortran_ptr() : nullptr; // input, optional
  BeamInitStruct _beam_init_used;
  bool print_p0c_shift_warning_lvalue;
  auto *_print_p0c_shift_warning{&print_p0c_shift_warning_lvalue};
  if (print_p0c_shift_warning.has_value()) {
    print_p0c_shift_warning_lvalue = print_p0c_shift_warning.value();
  } else {
    _print_p0c_shift_warning = nullptr;
  }
  bool conserve_momentum_lvalue;
  auto *_conserve_momentum{&conserve_momentum_lvalue};
  if (conserve_momentum.has_value()) {
    conserve_momentum_lvalue = conserve_momentum.value();
  } else {
    _conserve_momentum = nullptr;
  }
  fortran_init_bunch_distribution(/* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* void* */ beam_init.get_fortran_ptr(),
                                  /* int& */ ix_bunch,
                                  /* void* */ _bunch.get_fortran_ptr(),
                                  /* bool& */ _err_flag,
                                  /* void* */ _modes,
                                  /* void* */ _beam_init_used.get_fortran_ptr(),
                                  /* bool* */ _print_p0c_shift_warning,
                                  /* bool* */ _conserve_momentum);
  return InitBunchDistribution{std::move(_bunch), _err_flag, std::move(_beam_init_used)};
}
void Bmad::init_complex_taylor_series(
    ComplexTaylorStruct &complex_taylor,
    int n_term,
    std::optional<bool> save
) {
  bool save_lvalue;
  auto *_save{&save_lvalue};
  if (save.has_value()) {
    save_lvalue = save.value();
  } else {
    _save = nullptr;
  }
  fortran_init_complex_taylor_series(/* void* */ complex_taylor.get_fortran_ptr(),
                                     /* int& */ n_term,
                                     /* bool* */ _save);
}
void Bmad::init_coord(
    CoordStruct &orb,
    FixedArray1D<Real, 6> vec,
    optional_ref<EleStruct> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin,
    std::optional<double> s_pos,
    std::optional<bool> random_on
) {
  // vec: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  int element_end_lvalue;
  auto *_element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto *_particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto *_E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto *_t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto *_shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  // spin: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _spin_desc;
  _spin_desc.rank = 1;
  if (spin.has_value()) {
    _spin_desc.data_ptr = spin->data();
    _spin_desc.dims[0] = spin->size();
  } else {
    _spin_desc.data_ptr = nullptr;
    _spin_desc.dims[0] = 0;
  }
  double s_pos_lvalue;
  auto *_s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  bool random_on_lvalue;
  auto *_random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_coord1(/* void* */ orb.get_fortran_ptr(),
                      /* Bmad::array_descriptor_t& */ _vec_desc,
                      /* void* */ _ele,
                      /* int* */ _element_end,
                      /* int* */ _particle,
                      /* int* */ _direction,
                      /* double* */ _E_photon,
                      /* double* */ _t_offset,
                      /* bool* */ _shift_vec6,
                      /* Bmad::array_descriptor_t& */ _spin_desc,
                      /* double* */ _s_pos,
                      /* bool* */ _random_on);
}
CoordStruct Bmad::init_coord(
    CoordStruct &orb_in,
    optional_ref<EleStruct> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin,
    std::optional<double> s_pos,
    std::optional<bool> random_on
) {
  CoordStruct _orb_out;
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  int element_end_lvalue;
  auto *_element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto *_particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto *_E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto *_t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto *_shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  // spin: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _spin_desc;
  _spin_desc.rank = 1;
  if (spin.has_value()) {
    _spin_desc.data_ptr = spin->data();
    _spin_desc.dims[0] = spin->size();
  } else {
    _spin_desc.data_ptr = nullptr;
    _spin_desc.dims[0] = 0;
  }
  double s_pos_lvalue;
  auto *_s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  bool random_on_lvalue;
  auto *_random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_coord2(/* void* */ _orb_out.get_fortran_ptr(),
                      /* void* */ orb_in.get_fortran_ptr(),
                      /* void* */ _ele,
                      /* int* */ _element_end,
                      /* int* */ _particle,
                      /* int* */ _direction,
                      /* double* */ _E_photon,
                      /* double* */ _t_offset,
                      /* bool* */ _shift_vec6,
                      /* Bmad::array_descriptor_t& */ _spin_desc,
                      /* double* */ _s_pos,
                      /* bool* */ _random_on);
  return std::move(_orb_out);
}
void Bmad::init_coord(
    CoordStruct &orb,
    optional_ref<EleStruct> ele,
    std::optional<int> element_end,
    std::optional<int> particle,
    std::optional<int> direction,
    std::optional<double> E_photon,
    std::optional<double> t_offset,
    std::optional<bool> shift_vec6,
    std::optional<FixedArray1D<Real, 3>> spin
) {
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  int element_end_lvalue;
  auto *_element_end{&element_end_lvalue};
  if (element_end.has_value()) {
    element_end_lvalue = element_end.value();
  } else {
    _element_end = nullptr;
  }
  int particle_lvalue;
  auto *_particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  double E_photon_lvalue;
  auto *_E_photon{&E_photon_lvalue};
  if (E_photon.has_value()) {
    E_photon_lvalue = E_photon.value();
  } else {
    _E_photon = nullptr;
  }
  double t_offset_lvalue;
  auto *_t_offset{&t_offset_lvalue};
  if (t_offset.has_value()) {
    t_offset_lvalue = t_offset.value();
  } else {
    _t_offset = nullptr;
  }
  bool shift_vec6_lvalue;
  auto *_shift_vec6{&shift_vec6_lvalue};
  if (shift_vec6.has_value()) {
    shift_vec6_lvalue = shift_vec6.value();
  } else {
    _shift_vec6 = nullptr;
  }
  // spin: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _spin_desc;
  _spin_desc.rank = 1;
  if (spin.has_value()) {
    _spin_desc.data_ptr = spin->data();
    _spin_desc.dims[0] = spin->size();
  } else {
    _spin_desc.data_ptr = nullptr;
    _spin_desc.dims[0] = 0;
  }
  fortran_init_coord3(/* void* */ orb.get_fortran_ptr(),
                      /* void* */ _ele,
                      /* int* */ _element_end,
                      /* int* */ _particle,
                      /* int* */ _direction,
                      /* double* */ _E_photon,
                      /* double* */ _t_offset,
                      /* bool* */ _shift_vec6,
                      /* Bmad::array_descriptor_t& */ _spin_desc);
}
void Bmad::init_custom(LatStruct &lat) { fortran_init_custom(/* void* */ lat.get_fortran_ptr()); }
EleStruct Bmad::init_ele(
    std::optional<int> key,
    std::optional<int> sub_key,
    std::optional<int> ix_ele,
    optional_ref<BranchStruct> branch
) {
  EleStruct _ele;
  int key_lvalue;
  auto *_key{&key_lvalue};
  if (key.has_value()) {
    key_lvalue = key.value();
  } else {
    _key = nullptr;
  }
  int sub_key_lvalue;
  auto *_sub_key{&sub_key_lvalue};
  if (sub_key.has_value()) {
    sub_key_lvalue = sub_key.value();
  } else {
    _sub_key = nullptr;
  }
  int ix_ele_lvalue;
  auto *_ix_ele{&ix_ele_lvalue};
  if (ix_ele.has_value()) {
    ix_ele_lvalue = ix_ele.value();
  } else {
    _ix_ele = nullptr;
  }
  auto *_branch = branch.has_value() ? branch->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_init_ele(/* void* */ _ele.get_fortran_ptr(),
                   /* int* */ _key,
                   /* int* */ _sub_key,
                   /* int* */ _ix_ele,
                   /* void* */ _branch);
  return std::move(_ele);
}
void Bmad::init_gg_taylor_series(
    GgTaylorStruct &gg_taylor,
    int n_term,
    std::optional<bool> save_old
) {
  bool save_old_lvalue;
  auto *_save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  fortran_init_gg_taylor_series(/* void* */ gg_taylor.get_fortran_ptr(),
                                /* int& */ n_term,
                                /* bool* */ _save_old);
}
LatStruct Bmad::init_lat(std::optional<int> n, std::optional<bool> init_beginning_ele) {
  LatStruct _lat;
  int n_lvalue;
  auto *_n{&n_lvalue};
  if (n.has_value()) {
    n_lvalue = n.value();
  } else {
    _n = nullptr;
  }
  bool init_beginning_ele_lvalue;
  auto *_init_beginning_ele{&init_beginning_ele_lvalue};
  if (init_beginning_ele.has_value()) {
    init_beginning_ele_lvalue = init_beginning_ele.value();
  } else {
    _init_beginning_ele = nullptr;
  }
  fortran_init_lat(/* void* */ _lat.get_fortran_ptr(),
                   /* int* */ _n,
                   /* bool* */ _init_beginning_ele);
  return std::move(_lat);
}
void Bmad::init_multipole_cache(EleStruct &ele) {
  fortran_init_multipole_cache(/* void* */ ele.get_fortran_ptr());
}
CoordStruct Bmad::init_photon_from_a_photon_init_ele(
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> random_on
) {
  CoordStruct _orbit;
  bool random_on_lvalue;
  auto *_random_on{&random_on_lvalue};
  if (random_on.has_value()) {
    random_on_lvalue = random_on.value();
  } else {
    _random_on = nullptr;
  }
  fortran_init_photon_from_a_photon_init_ele(/* void* */ ele.get_fortran_ptr(),
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
    std::optional<double> energy_integ_prob
) {
  double vert_angle_min_lvalue;
  auto *_vert_angle_min{&vert_angle_min_lvalue};
  if (vert_angle_min.has_value()) {
    vert_angle_min_lvalue = vert_angle_min.value();
  } else {
    _vert_angle_min = nullptr;
  }
  double vert_angle_max_lvalue;
  auto *_vert_angle_max{&vert_angle_max_lvalue};
  if (vert_angle_max.has_value()) {
    vert_angle_max_lvalue = vert_angle_max.value();
  } else {
    _vert_angle_max = nullptr;
  }
  bool vert_angle_symmetric_lvalue;
  auto *_vert_angle_symmetric{&vert_angle_symmetric_lvalue};
  if (vert_angle_symmetric.has_value()) {
    vert_angle_symmetric_lvalue = vert_angle_symmetric.value();
  } else {
    _vert_angle_symmetric = nullptr;
  }
  double energy_integ_prob_lvalue;
  auto *_energy_integ_prob{&energy_integ_prob_lvalue};
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
      /* double& */ _integ_prob
  );
  return InitPhotonIntegProb{_E_photon, _integ_prob};
}
BunchStruct Bmad::init_spin_distribution(BeamInitStruct &beam_init, EleStruct &ele) {
  BunchStruct _bunch;
  fortran_init_spin_distribution(/* void* */ beam_init.get_fortran_ptr(),
                                 /* void* */ _bunch.get_fortran_ptr(),
                                 /* void* */ ele.get_fortran_ptr());
  return std::move(_bunch);
}
void Bmad::init_surface_segment(PhotonElementStruct &phot, int ix_pt, int iy_pt) {
  fortran_init_surface_segment(/* void* */ phot.get_fortran_ptr(),
                               /* int& */ ix_pt,
                               /* int& */ iy_pt);
}
void Bmad::init_taylor_series(TaylorStruct &bmad_taylor, int n_term, std::optional<bool> save_old) {
  bool save_old_lvalue;
  auto *_save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  fortran_init_taylor_series(/* void* */ bmad_taylor.get_fortran_ptr(),
                             /* int& */ n_term,
                             /* bool* */ _save_old);
}
std::optional<WakeStruct> Bmad::init_wake(
    int n_sr_long,
    int n_sr_trans,
    int n_sr_z,
    int n_lr_mode,
    std::optional<bool> always_allocate
) {
  void *_wake;
  bool always_allocate_lvalue;
  auto *_always_allocate{&always_allocate_lvalue};
  if (always_allocate.has_value()) {
    always_allocate_lvalue = always_allocate.value();
  } else {
    _always_allocate = nullptr;
  }
  fortran_init_wake(/* void* */ &_wake,
                    /* int& */ n_sr_long,
                    /* int& */ n_sr_trans,
                    /* int& */ n_sr_z,
                    /* int& */ n_lr_mode,
                    /* bool* */ _always_allocate);
  return std::move((_wake ? std::make_optional<WakeStruct>(_wake) : std::nullopt));
}
void Bmad::insert_element(
    LatStruct &lat,
    EleStruct &insert_ele,
    int ix_ele,
    std::optional<int> ix_branch,
    std::optional<CoordStructAlloc1D> orbit
) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  // intent=inout allocatable type array
  auto *_orbit = orbit.has_value() ? orbit->get_fortran_ptr() : nullptr; // input, optional
  fortran_insert_element(/* void* */ lat.get_fortran_ptr(),
                         /* void* */ insert_ele.get_fortran_ptr(),
                         /* int& */ ix_ele,
                         /* int* */ _ix_branch,
                         /* void* */ _orbit);
}
void Bmad::integrand_base(double t, FArray1D<Real> &args, double func_retval__) {
  // args: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  fortran_integrand_base(
      /* double& */ t,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ func_retval__
  );
}
double Bmad::integrate_psi(double bound, double p0, FixedArray1D<Real, 8> args) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  double _result{};
  fortran_integrate_psi(
      /* double& */ bound,
      /* double& */ p0,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ _result
  );
  return _result;
}
void Bmad::integrated_mats(
    EleStructArray1D eles,
    CoordStructArray1D coos,
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    NormalModesStruct &mode
) {
  // eles: EleStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _eles_desc;
  _eles_desc.rank = 1;
  _eles_desc.data_ptr = eles.data();
  _eles_desc.dims[0] = eles.size();
  _eles_desc.strides[0] = 1;
  // coos: CoordStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _coos_desc;
  _coos_desc.rank = 1;
  _coos_desc.data_ptr = coos.data();
  _coos_desc.dims[0] = coos.size();
  _coos_desc.strides[0] = 1;
  // Lambda: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Lambda_desc;
  _Lambda_desc.rank = 2;
  std::complex<double> _Lambda_vec[6 * 6];
  _Lambda_desc.data_ptr = _Lambda_vec;
  _Lambda_desc.dims[0] = 6;
  _Lambda_desc.dims[1] = 6;
  matrix_to_vec(Lambda, _Lambda_vec);
  // Theta: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Theta_desc;
  _Theta_desc.rank = 2;
  std::complex<double> _Theta_vec[6 * 6];
  _Theta_desc.data_ptr = _Theta_vec;
  _Theta_desc.dims[0] = 6;
  _Theta_desc.dims[1] = 6;
  matrix_to_vec(Theta, _Theta_vec);
  // Iota: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Iota_desc;
  _Iota_desc.rank = 2;
  std::complex<double> _Iota_vec[6 * 6];
  _Iota_desc.data_ptr = _Iota_vec;
  _Iota_desc.dims[0] = 6;
  _Iota_desc.dims[1] = 6;
  matrix_to_vec(Iota, _Iota_vec);
  fortran_integrated_mats(
      /* Bmad::array_descriptor_t& */ _eles_desc,
      /* Bmad::array_descriptor_t& */ _coos_desc,
      /* Bmad::array_descriptor_t& */ _Lambda_desc,
      /* Bmad::array_descriptor_t& */ _Theta_desc,
      /* Bmad::array_descriptor_t& */ _Iota_desc,
      /* void* */ mode.get_fortran_ptr()
  );
  vec_to_matrix(_Lambda_vec, Lambda);
  vec_to_matrix(_Theta_vec, Theta);
  vec_to_matrix(_Iota_vec, Iota);
}
void Bmad::integration_timer(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &start,
    CoordStruct &orb_max,
    double tol
) {
  fortran_integration_timer_ele(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ param.get_fortran_ptr(),
                                /* void* */ start.get_fortran_ptr(),
                                /* void* */ orb_max.get_fortran_ptr(),
                                /* double& */ tol);
}
void Bmad::integration_timer(
    Fibre &a_fibre,
    FixedArray1D<Real, 6> orbit,
    FixedArray1D<Real, 6> orbit_max,
    double tol_dp
) {
  // orbit: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  _orbit_desc.data_ptr = orbit.data();
  _orbit_desc.dims[0] = orbit.size();
  // orbit_max: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orbit_max_desc;
  _orbit_max_desc.rank = 1;
  _orbit_max_desc.data_ptr = orbit_max.data();
  _orbit_max_desc.dims[0] = orbit_max.size();
  fortran_integration_timer_fibre(/* void* */ a_fibre.get_fortran_ptr(),
                                  /* Bmad::array_descriptor_t& */ _orbit_desc,
                                  /* Bmad::array_descriptor_t& */ _orbit_max_desc,
                                  /* double& */ tol_dp);
}
FixedArray1D<Real, 3> Bmad::ion_kick(
    CoordStruct &orbit,
    FixedArray1D<Real, 2> r_beam,
    double n_beam_part,
    TwissStruct &a_twiss,
    TwissStruct &b_twiss,
    double sig_ee
) {
  // r_beam: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _r_beam_desc;
  _r_beam_desc.rank = 1;
  _r_beam_desc.data_ptr = r_beam.data();
  _r_beam_desc.dims[0] = r_beam.size();
  // kick: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _kick_desc;
  _kick_desc.rank = 1;
  FixedArray1D<Real, 3> _kick;
  _kick_desc.data_ptr = _kick.data();
  _kick_desc.dims[0] = _kick.size();
  fortran_ion_kick(/* void* */ orbit.get_fortran_ptr(),
                   /* Bmad::array_descriptor_t& */ _r_beam_desc,
                   /* double& */ n_beam_part,
                   /* void* */ a_twiss.get_fortran_ptr(),
                   /* void* */ b_twiss.get_fortran_ptr(),
                   /* double& */ sig_ee,
                   /* Bmad::array_descriptor_t& */ _kick_desc);
  return _kick;
}
bool Bmad::is_attribute(int ix_attrib, int which) {
  bool _is_attrib{};
  fortran_is_attribute(/* int& */ ix_attrib, /* int& */ which, /* bool& */ _is_attrib);
  return _is_attrib;
}
int Bmad::key_name_to_key_index(std::string key_str, std::optional<bool> abbrev_allowed) {
  auto _key_str = key_str.c_str();
  bool abbrev_allowed_lvalue;
  auto *_abbrev_allowed{&abbrev_allowed_lvalue};
  if (abbrev_allowed.has_value()) {
    abbrev_allowed_lvalue = abbrev_allowed.value();
  } else {
    _abbrev_allowed = nullptr;
  }
  int _key_index{};
  fortran_key_name_to_key_index(
      /* const char* */ _key_str,
      /* bool* */ _abbrev_allowed,
      /* int& */ _key_index
  );
  return _key_index;
}
Bmad::KickVectorCalc Bmad::kick_vector_calc(
    EleStruct &ele,
    LatParamStruct &param,
    double s_body,
    CoordStruct &orbit,
    std::optional<bool> print_err
) {
  // dr_ds: out NOT (CppWrapperGeneralArgumentArray) (['11'])
  Bmad::array_descriptor_t _dr_ds_desc;
  _dr_ds_desc.rank = 1;
  FixedArray1D<Real, 11> _dr_ds;
  _dr_ds_desc.data_ptr = _dr_ds.data();
  _dr_ds_desc.dims[0] = _dr_ds.size();
  bool _err{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_kick_vector_calc(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* double& */ s_body,
                           /* void* */ orbit.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _dr_ds_desc,
                           /* bool& */ _err,
                           /* bool* */ _print_err);
  return KickVectorCalc{_dr_ds, _err};
}
void Bmad::kill_complex_taylor(ComplexTaylorStructArray1D complex_taylor) {
  // complex_taylor: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_desc;
  _complex_taylor_desc.rank = 1;
  _complex_taylor_desc.data_ptr = complex_taylor.data();
  _complex_taylor_desc.dims[0] = complex_taylor.size();
  _complex_taylor_desc.strides[0] = 1;
  fortran_kill_complex_taylor(/* Bmad::array_descriptor_t& */ _complex_taylor_desc);
}
void Bmad::kill_ptc_layouts(LatStruct &lat) {
  fortran_kill_ptc_layouts(/* void* */ lat.get_fortran_ptr());
}
void Bmad::kill_taylor(TaylorStructArray1D bmad_taylor) {
  // bmad_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _bmad_taylor_desc;
  _bmad_taylor_desc.rank = 1;
  _bmad_taylor_desc.data_ptr = bmad_taylor.data();
  _bmad_taylor_desc.dims[0] = bmad_taylor.size();
  _bmad_taylor_desc.strides[0] = 1;
  fortran_kill_taylor(/* Bmad::array_descriptor_t& */ _bmad_taylor_desc);
}
std::string Bmad::kind_name(int this_kind) {
  char _kind_str[4096];
  fortran_kind_name(/* int* */ &this_kind, /* const char* */ _kind_str);
  return _kind_str;
}
Bmad::KnotInterpolate Bmad::knot_interpolate(
    FArray1D<Real> &x_knot,
    FArray1D<Real> &y_knot,
    double x_pt,
    int interpolation
) {
  // x_knot: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_knot_desc;
  _x_knot_desc.rank = 1;
  _x_knot_desc.data_ptr = x_knot.data();
  _x_knot_desc.dims[0] = x_knot.size();
  // y_knot: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _y_knot_desc;
  _y_knot_desc.rank = 1;
  _y_knot_desc.data_ptr = y_knot.data();
  _y_knot_desc.dims[0] = y_knot.size();
  bool _err_flag{};
  double _y_pt{};
  fortran_knot_interpolate(
      /* Bmad::array_descriptor_t& */ _x_knot_desc,
      /* Bmad::array_descriptor_t& */ _y_knot_desc,
      /* double& */ x_pt,
      /* int& */ interpolation,
      /* bool& */ _err_flag,
      /* double& */ _y_pt
  );
  return KnotInterpolate{_err_flag, _y_pt};
}
void Bmad::knots_to_string(FArray1D<Real> &x_knot, FArray1D<Real> &y_knot, std::string str) {
  // x_knot: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_knot_desc;
  _x_knot_desc.rank = 1;
  _x_knot_desc.data_ptr = x_knot.data();
  _x_knot_desc.dims[0] = x_knot.size();
  // y_knot: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _y_knot_desc;
  _y_knot_desc.rank = 1;
  _y_knot_desc.data_ptr = y_knot.data();
  _y_knot_desc.dims[0] = y_knot.size();
  auto _str = str.c_str();
  fortran_knots_to_string(
      /* Bmad::array_descriptor_t& */ _x_knot_desc,
      /* Bmad::array_descriptor_t& */ _y_knot_desc,
      /* const char* */ _str
  );
}
void Bmad::lafun(double x, double y, double z, double res) {
  fortran_lafun(/* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
bool Bmad::lat_compute_ref_energy_and_time(LatStruct &lat) {
  bool _err_flag{};
  fortran_lat_compute_ref_energy_and_time(/* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::lat_ele_locator(
    std::string loc_str,
    LatStruct &lat,
    ElePointerStructAlloc1D eles,
    int &n_loc,
    std::optional<bool> above_ubound_is_err,
    std::optional<int> ix_dflt_branch,
    std::optional<bool> order_by_index,
    std::optional<bool> append_eles
) {
  auto _loc_str = loc_str.c_str();
  // intent=inout allocatable type array
  bool _err{};
  bool above_ubound_is_err_lvalue;
  auto *_above_ubound_is_err{&above_ubound_is_err_lvalue};
  if (above_ubound_is_err.has_value()) {
    above_ubound_is_err_lvalue = above_ubound_is_err.value();
  } else {
    _above_ubound_is_err = nullptr;
  }
  int ix_dflt_branch_lvalue;
  auto *_ix_dflt_branch{&ix_dflt_branch_lvalue};
  if (ix_dflt_branch.has_value()) {
    ix_dflt_branch_lvalue = ix_dflt_branch.value();
  } else {
    _ix_dflt_branch = nullptr;
  }
  bool order_by_index_lvalue;
  auto *_order_by_index{&order_by_index_lvalue};
  if (order_by_index.has_value()) {
    order_by_index_lvalue = order_by_index.value();
  } else {
    _order_by_index = nullptr;
  }
  bool append_eles_lvalue;
  auto *_append_eles{&append_eles_lvalue};
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
      /* bool* */ _append_eles
  );
  return _err;
}
void Bmad::lat_equal_lat(LatStruct &lat_out, LatStruct &lat_in) {
  fortran_lat_equal_lat(/* void* */ lat_out.get_fortran_ptr(),
                        /* void* */ lat_in.get_fortran_ptr());
}
void Bmad::lat_geometry(LatStruct &lat) { fortran_lat_geometry(/* void* */ lat.get_fortran_ptr()); }
bool Bmad::lat_make_mat6(
    LatStruct &lat,
    std::optional<int> ix_ele,
    std::optional<CoordStructArray1D> ref_orb,
    std::optional<int> ix_branch
) {
  int ix_ele_lvalue;
  auto *_ix_ele{&ix_ele_lvalue};
  if (ix_ele.has_value()) {
    ix_ele_lvalue = ix_ele.value();
  } else {
    _ix_ele = nullptr;
  }
  // ref_orb: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ref_orb_desc;
  _ref_orb_desc.rank = 1;
  if (ref_orb) {
    _ref_orb_desc.data_ptr = ref_orb->data();
    _ref_orb_desc.dims[0] = ref_orb->size();
  } else {
    _ref_orb_desc.data_ptr = nullptr;
    _ref_orb_desc.dims[0] = 0;
  }
  _ref_orb_desc.strides[0] = 1;
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  fortran_lat_make_mat6(/* void* */ lat.get_fortran_ptr(),
                        /* int* */ _ix_ele,
                        /* Bmad::array_descriptor_t& */ _ref_orb_desc,
                        /* int* */ _ix_branch,
                        /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::lat_sanity_check(LatStruct &lat) {
  bool _err_flag{};
  fortran_lat_sanity_check(/* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::lat_to_ptc_layout(LatStruct &lat) {
  fortran_lat_to_ptc_layout(/* void* */ lat.get_fortran_ptr());
}
void Bmad::lat_vec_equal_lat_vec(LatStructArray1D lat1, LatStructArray1D lat2) {
  // lat1: LatStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _lat1_desc;
  _lat1_desc.rank = 1;
  _lat1_desc.data_ptr = lat1.data();
  _lat1_desc.dims[0] = lat1.size();
  _lat1_desc.strides[0] = 1;
  // lat2: LatStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _lat2_desc;
  _lat2_desc.rank = 1;
  _lat2_desc.data_ptr = lat2.data();
  _lat2_desc.dims[0] = lat2.size();
  _lat2_desc.strides[0] = 1;
  fortran_lat_vec_equal_lat_vec(
      /* Bmad::array_descriptor_t& */ _lat1_desc,
      /* Bmad::array_descriptor_t& */ _lat2_desc
  );
}
bool Bmad::lattice_bookkeeper(LatStruct &lat) {
  bool _err_flag{};
  fortran_lattice_bookkeeper(/* void* */ lat.get_fortran_ptr(), /* bool& */ _err_flag);
  return _err_flag;
}
void Bmad::lcavity_rf_step_setup(EleStruct &ele) {
  fortran_lcavity_rf_step_setup(/* void* */ ele.get_fortran_ptr());
}
void Bmad::linear_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_linear_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ param.get_fortran_ptr(),
                                /* int& */ particle_at,
                                /* void* */ orb.get_fortran_ptr(),
                                /* Bmad::array_descriptor_t& */ _mat6_desc,
                                /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
Bmad::LinearCoef Bmad::linear_coef(ExpressionAtomStructArray1D stack) {
  // stack: ExpressionAtomStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _stack_desc;
  _stack_desc.rank = 1;
  _stack_desc.data_ptr = stack.data();
  _stack_desc.dims[0] = stack.size();
  _stack_desc.strides[0] = 1;
  bool _err_flag{};
  double _coef{};
  fortran_linear_coef(
      /* Bmad::array_descriptor_t& */ _stack_desc,
      /* bool& */ _err_flag,
      /* double& */ _coef
  );
  return LinearCoef{_err_flag, _coef};
}
TaylorStructArray1D Bmad::linear_to_spin_taylor(FixedArray2D<Real, 4, 7> q_map) {
  // q_map: in NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _q_map_desc;
  _q_map_desc.rank = 2;
  double _q_map_vec[4 * 7];
  _q_map_desc.data_ptr = _q_map_vec;
  _q_map_desc.dims[0] = 4;
  _q_map_desc.dims[1] = 7;
  matrix_to_vec(q_map, _q_map_vec);
  // spin_taylor: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_taylor_desc;
  _spin_taylor_desc.rank = 1;
  // Output-only type array
  auto spin_taylor = TaylorStructAlloc1D(4);
  _spin_taylor_desc.data_ptr = spin_taylor.get_fortran_ptr();
  _spin_taylor_desc.dims[0] = spin_taylor.size();

  _spin_taylor_desc.strides[0] = 1;
  fortran_linear_to_spin_taylor(
      /* Bmad::array_descriptor_t& */ _q_map_desc,
      /* Bmad::array_descriptor_t& */ _spin_taylor_desc
  );
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
      /* bool& */ _err_flag
  );
  return LoadParseLine{_end_of_file, _err_flag};
}
bool Bmad::lord_edge_aligned(EleStruct &slave, int slave_edge, EleStruct &lord) {
  bool _is_aligned{};
  fortran_lord_edge_aligned(/* void* */ slave.get_fortran_ptr(),
                            /* int& */ slave_edge,
                            /* void* */ lord.get_fortran_ptr(),
                            /* bool& */ _is_aligned);
  return _is_aligned;
}
double Bmad::low_energy_z_correction(
    CoordStruct &orbit,
    EleStruct &ele,
    double ds,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  double _dz{};
  fortran_low_energy_z_correction(/* void* */ orbit.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* double& */ ds,
                                  /* Bmad::array_descriptor_t& */ _mat6_desc,
                                  /* bool* */ _make_matrix,
                                  /* double& */ _dz);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return _dz;
}
MadMapStruct Bmad::mad_add_offsets_and_multipoles(EleStruct &ele) {
  MadMapStruct _map;
  fortran_mad_add_offsets_and_multipoles(/* void* */ ele.get_fortran_ptr(),
                                         /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_concat_map2(MadMapStruct &map1, MadMapStruct &map2) {
  MadMapStruct _map3;
  fortran_mad_concat_map2(/* void* */ map1.get_fortran_ptr(),
                          /* void* */ map2.get_fortran_ptr(),
                          /* void* */ _map3.get_fortran_ptr());
  return std::move(_map3);
}
MadMapStruct Bmad::mad_drift(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_drift(/* void* */ ele.get_fortran_ptr(),
                    /* void* */ energy.get_fortran_ptr(),
                    /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_elsep(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_elsep(/* void* */ ele.get_fortran_ptr(),
                    /* void* */ energy.get_fortran_ptr(),
                    /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
void Bmad::mad_map_to_taylor(
    MadMapStruct &map,
    MadEnergyStruct &energy,
    TaylorStructArray1D taylor
) {
  // taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_desc;
  _taylor_desc.rank = 1;
  _taylor_desc.data_ptr = taylor.data();
  _taylor_desc.dims[0] = taylor.size();
  _taylor_desc.strides[0] = 1;
  fortran_mad_map_to_taylor(/* void* */ map.get_fortran_ptr(),
                            /* void* */ energy.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _taylor_desc);
}
MadMapStruct Bmad::mad_quadrupole(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_quadrupole(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ energy.get_fortran_ptr(),
                         /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_rfcavity(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_rfcavity(/* void* */ ele.get_fortran_ptr(),
                       /* void* */ energy.get_fortran_ptr(),
                       /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_sbend(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_sbend(/* void* */ ele.get_fortran_ptr(),
                    /* void* */ energy.get_fortran_ptr(),
                    /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_sbend_body(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_sbend_body(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ energy.get_fortran_ptr(),
                         /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_sbend_fringe(EleStruct &ele, MadEnergyStruct &energy, bool into) {
  MadMapStruct _map;
  fortran_mad_sbend_fringe(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ energy.get_fortran_ptr(),
                           /* bool& */ into,
                           /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_sextupole(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_sextupole(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ energy.get_fortran_ptr(),
                        /* void* */ _map.get_fortran_ptr());
  return std::move(_map);
}
MadMapStruct Bmad::mad_solenoid(EleStruct &ele, MadEnergyStruct &energy) {
  MadMapStruct _map;
  fortran_mad_solenoid(/* void* */ ele.get_fortran_ptr(),
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
      /* double& */ _f
  );
  return MadTmfoc{_c, _s, _d, _f};
}
void Bmad::mad_tmsymm(FixedArray3D<Real, 6, 6, 6> te) {
  // te: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6', '6'])
  Bmad::array_descriptor_t _te_desc;
  _te_desc.rank = 3;
  double _te_vec[6 * 6 * 6];
  _te_desc.data_ptr = _te_vec;
  _te_desc.dims[0] = 6;
  _te_desc.dims[1] = 6;
  _te_desc.dims[2] = 6;
  tensor_to_vec(te, _te_vec);
  fortran_mad_tmsymm(/* Bmad::array_descriptor_t& */ _te_desc);
  vec_to_tensor(_te_vec, te);
}
void Bmad::mad_tmtilt(MadMapStruct &map, double tilt) {
  fortran_mad_tmtilt(/* void* */ map.get_fortran_ptr(), /* double& */ tilt);
}
CoordStruct Bmad::mad_track1(CoordStruct &c0, MadMapStruct &map) {
  CoordStruct _c1;
  fortran_mad_track1(/* void* */ c0.get_fortran_ptr(),
                     /* void* */ map.get_fortran_ptr(),
                     /* void* */ _c1.get_fortran_ptr());
  return std::move(_c1);
}
void Bmad::make_g2_mats(
    TwissStruct &twiss,
    FixedArray2D<Real, 2, 2> g2_mat,
    FixedArray2D<Real, 2, 2> g2_inv_mat
) {
  // g2_mat: inout NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _g2_mat_desc;
  _g2_mat_desc.rank = 2;
  double _g2_mat_vec[2 * 2];
  _g2_mat_desc.data_ptr = _g2_mat_vec;
  _g2_mat_desc.dims[0] = 2;
  _g2_mat_desc.dims[1] = 2;
  matrix_to_vec(g2_mat, _g2_mat_vec);
  // g2_inv_mat: inout NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _g2_inv_mat_desc;
  _g2_inv_mat_desc.rank = 2;
  double _g2_inv_mat_vec[2 * 2];
  _g2_inv_mat_desc.data_ptr = _g2_inv_mat_vec;
  _g2_inv_mat_desc.dims[0] = 2;
  _g2_inv_mat_desc.dims[1] = 2;
  matrix_to_vec(g2_inv_mat, _g2_inv_mat_vec);
  fortran_make_g2_mats(/* void* */ twiss.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _g2_mat_desc,
                       /* Bmad::array_descriptor_t& */ _g2_inv_mat_desc);
  vec_to_matrix(_g2_mat_vec, g2_mat);
  vec_to_matrix(_g2_inv_mat_vec, g2_inv_mat);
}
Bmad::MakeGMats Bmad::make_g_mats(EleStruct &ele) {
  // g_mat: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _g_mat_desc;
  _g_mat_desc.rank = 2;
  FixedArray2D<Real, 4, 4> g_mat;
  double _g_mat_vec[4 * 4];
  _g_mat_desc.data_ptr = _g_mat_vec;
  _g_mat_desc.dims[0] = 4;
  _g_mat_desc.dims[1] = 4;
  // g_inv_mat: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _g_inv_mat_desc;
  _g_inv_mat_desc.rank = 2;
  FixedArray2D<Real, 4, 4> g_inv_mat;
  double _g_inv_mat_vec[4 * 4];
  _g_inv_mat_desc.data_ptr = _g_inv_mat_vec;
  _g_inv_mat_desc.dims[0] = 4;
  _g_inv_mat_desc.dims[1] = 4;
  fortran_make_g_mats(/* void* */ ele.get_fortran_ptr(),
                      /* Bmad::array_descriptor_t& */ _g_mat_desc,
                      /* Bmad::array_descriptor_t& */ _g_inv_mat_desc);
  vec_to_matrix(_g_mat_vec, g_mat);
  vec_to_matrix(_g_inv_mat_vec, g_inv_mat);
  return MakeGMats{g_mat, g_inv_mat};
}
Bmad::MakeHvbp Bmad::make_hvbp(FixedArray2D<Real, 6, 6> N) {
  // N: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _N_desc;
  _N_desc.rank = 2;
  double _N_vec[6 * 6];
  _N_desc.data_ptr = _N_vec;
  _N_desc.dims[0] = 6;
  _N_desc.dims[1] = 6;
  matrix_to_vec(N, _N_vec);
  // B: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _B_desc;
  _B_desc.rank = 2;
  FixedArray2D<Real, 6, 6> B;
  double _B_vec[6 * 6];
  _B_desc.data_ptr = _B_vec;
  _B_desc.dims[0] = 6;
  _B_desc.dims[1] = 6;
  // V: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _V_desc;
  _V_desc.rank = 2;
  FixedArray2D<Real, 6, 6> V;
  double _V_vec[6 * 6];
  _V_desc.data_ptr = _V_vec;
  _V_desc.dims[0] = 6;
  _V_desc.dims[1] = 6;
  // H: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _H_desc;
  _H_desc.rank = 2;
  FixedArray2D<Real, 6, 6> H;
  double _H_vec[6 * 6];
  _H_desc.data_ptr = _H_vec;
  _H_desc.dims[0] = 6;
  _H_desc.dims[1] = 6;
  // Vbar: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Vbar_desc;
  _Vbar_desc.rank = 2;
  FixedArray2D<Real, 6, 6> Vbar;
  double _Vbar_vec[6 * 6];
  _Vbar_desc.data_ptr = _Vbar_vec;
  _Vbar_desc.dims[0] = 6;
  _Vbar_desc.dims[1] = 6;
  // Hbar: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Hbar_desc;
  _Hbar_desc.rank = 2;
  FixedArray2D<Real, 6, 6> Hbar;
  double _Hbar_vec[6 * 6];
  _Hbar_desc.data_ptr = _Hbar_vec;
  _Hbar_desc.dims[0] = 6;
  _Hbar_desc.dims[1] = 6;
  fortran_make_hvbp(
      /* Bmad::array_descriptor_t& */ _N_desc,
      /* Bmad::array_descriptor_t& */ _B_desc,
      /* Bmad::array_descriptor_t& */ _V_desc,
      /* Bmad::array_descriptor_t& */ _H_desc,
      /* Bmad::array_descriptor_t& */ _Vbar_desc,
      /* Bmad::array_descriptor_t& */ _Hbar_desc
  );
  vec_to_matrix(_B_vec, B);
  vec_to_matrix(_V_vec, V);
  vec_to_matrix(_H_vec, H);
  vec_to_matrix(_Vbar_vec, Vbar);
  vec_to_matrix(_Hbar_vec, Hbar);
  return MakeHvbp{B, V, H, Vbar, Hbar};
}
LatStruct Bmad::make_hybrid_lat(
    LatStruct &lat_in,
    std::optional<bool> use_taylor,
    std::optional<CoordArrayStructArray1D> orb0_arr
) {
  LatStruct _lat_out;
  bool use_taylor_lvalue;
  auto *_use_taylor{&use_taylor_lvalue};
  if (use_taylor.has_value()) {
    use_taylor_lvalue = use_taylor.value();
  } else {
    _use_taylor = nullptr;
  }
  // orb0_arr: CoordArrayStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb0_arr_desc;
  _orb0_arr_desc.rank = 1;
  if (orb0_arr) {
    _orb0_arr_desc.data_ptr = orb0_arr->data();
    _orb0_arr_desc.dims[0] = orb0_arr->size();
  } else {
    _orb0_arr_desc.data_ptr = nullptr;
    _orb0_arr_desc.dims[0] = 0;
  }
  _orb0_arr_desc.strides[0] = 1;
  fortran_make_hybrid_lat(/* void* */ lat_in.get_fortran_ptr(),
                          /* void* */ _lat_out.get_fortran_ptr(),
                          /* bool* */ _use_taylor,
                          /* Bmad::array_descriptor_t& */ _orb0_arr_desc);
  return std::move(_lat_out);
}
Bmad::MakeMadMap Bmad::make_mad_map(EleStruct &ele, LatParamStruct &param) {
  MadEnergyStruct _energy;
  MadMapStruct _map;
  fortran_make_mad_map(/* void* */ ele.get_fortran_ptr(),
                       /* void* */ param.get_fortran_ptr(),
                       /* void* */ _energy.get_fortran_ptr(),
                       /* void* */ _map.get_fortran_ptr());
  return MakeMadMap{std::move(_energy), std::move(_map)};
}
Bmad::MakeMat6
Bmad::make_mat6(EleStruct &ele, LatParamStruct &param, optional_ref<CoordStruct> start_orb) {
  auto *_start_orb =
      start_orb.has_value() ? start_orb->get().get_fortran_ptr() : nullptr; // input, optional
  CoordStruct _end_orb;
  bool _err_flag{};
  fortran_make_mat6(/* void* */ ele.get_fortran_ptr(),
                    /* void* */ param.get_fortran_ptr(),
                    /* void* */ _start_orb,
                    /* void* */ _end_orb.get_fortran_ptr(),
                    /* bool& */ _err_flag);
  return MakeMat6{std::move(_end_orb), _err_flag};
}
Bmad::MakeMat6Bmad
Bmad::make_mat6_bmad(EleStruct &ele, LatParamStruct &param, CoordStruct &start_orb) {
  CoordStruct _end_orb;
  bool _err{};
  fortran_make_mat6_bmad(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ param.get_fortran_ptr(),
                         /* void* */ start_orb.get_fortran_ptr(),
                         /* void* */ _end_orb.get_fortran_ptr(),
                         /* bool& */ _err);
  return MakeMat6Bmad{std::move(_end_orb), _err};
}
Bmad::MakeMat6BmadPhoton
Bmad::make_mat6_bmad_photon(EleStruct &ele, LatParamStruct &param, CoordStruct &start_orb) {
  CoordStruct _end_orb;
  bool _err{};
  fortran_make_mat6_bmad_photon(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ param.get_fortran_ptr(),
                                /* void* */ start_orb.get_fortran_ptr(),
                                /* void* */ _end_orb.get_fortran_ptr(),
                                /* bool& */ _err);
  return MakeMat6BmadPhoton{std::move(_end_orb), _err};
}
void Bmad::make_mat6_high_energy_space_charge(EleStruct &ele, LatParamStruct &param) {
  fortran_make_mat6_high_energy_space_charge(/* void* */ ele.get_fortran_ptr(),
                                             /* void* */ param.get_fortran_ptr());
}
CoordStruct Bmad::make_mat6_mad(EleStruct &ele, LatParamStruct &param, CoordStruct &c0) {
  CoordStruct _c1;
  fortran_make_mat6_mad(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* void* */ c0.get_fortran_ptr(),
                        /* void* */ _c1.get_fortran_ptr());
  return std::move(_c1);
}
CoordStruct Bmad::make_mat6_symp_lie_ptc(EleStruct &ele, CoordStruct &start_orb) {
  CoordStruct _end_orb;
  fortran_make_mat6_symp_lie_ptc(/* void* */ ele.get_fortran_ptr(),
                                 /* void* */ start_orb.get_fortran_ptr(),
                                 /* void* */ _end_orb.get_fortran_ptr());
  return std::move(_end_orb);
}
CoordStruct
Bmad::make_mat6_taylor(EleStruct &ele, CoordStruct &start_orb, std::optional<bool> err_flag) {
  CoordStruct _end_orb;
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_make_mat6_taylor(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ start_orb.get_fortran_ptr(),
                           /* void* */ _end_orb.get_fortran_ptr(),
                           /* bool* */ _err_flag);
  return std::move(_end_orb);
}
Bmad::MakeMat6Tracking Bmad::make_mat6_tracking(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &start_orb,
    std::optional<bool> spin_only
) {
  CoordStruct _end_orb;
  bool _err_flag{};
  bool spin_only_lvalue;
  auto *_spin_only{&spin_only_lvalue};
  if (spin_only.has_value()) {
    spin_only_lvalue = spin_only.value();
  } else {
    _spin_only = nullptr;
  }
  fortran_make_mat6_tracking(/* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* void* */ start_orb.get_fortran_ptr(),
                             /* void* */ _end_orb.get_fortran_ptr(),
                             /* bool& */ _err_flag,
                             /* bool* */ _spin_only);
  return MakeMat6Tracking{std::move(_end_orb), _err_flag};
}
Bmad::MakeN
Bmad::make_n(FixedArray2D<Real, 6, 6> t6, std::optional<FixedArray1D<Real, 3>> abz_tunes) {
  // t6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_desc;
  _t6_desc.rank = 2;
  double _t6_vec[6 * 6];
  _t6_desc.data_ptr = _t6_vec;
  _t6_desc.dims[0] = 6;
  _t6_desc.dims[1] = 6;
  matrix_to_vec(t6, _t6_vec);
  // N: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _N_desc;
  _N_desc.rank = 2;
  FixedArray2D<Real, 6, 6> N;
  double _N_vec[6 * 6];
  _N_desc.data_ptr = _N_vec;
  _N_desc.dims[0] = 6;
  _N_desc.dims[1] = 6;
  bool _err_flag{};
  // abz_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  if (abz_tunes.has_value()) {
    _abz_tunes_desc.data_ptr = abz_tunes->data();
    _abz_tunes_desc.dims[0] = abz_tunes->size();
  } else {
    _abz_tunes_desc.data_ptr = nullptr;
    _abz_tunes_desc.dims[0] = 0;
  }
  // tunes_out: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _tunes_out_desc;
  _tunes_out_desc.rank = 1;
  FixedArray1D<Real, 3> _tunes_out;
  _tunes_out_desc.data_ptr = _tunes_out.data();
  _tunes_out_desc.dims[0] = _tunes_out.size();
  // U: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _U_desc;
  _U_desc.rank = 2;
  FixedArray2D<Real, 6, 6> U;
  double _U_vec[6 * 6];
  _U_desc.data_ptr = _U_vec;
  _U_desc.dims[0] = 6;
  _U_desc.dims[1] = 6;
  fortran_make_n(
      /* Bmad::array_descriptor_t& */ _t6_desc,
      /* Bmad::array_descriptor_t& */ _N_desc,
      /* bool& */ _err_flag,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc,
      /* Bmad::array_descriptor_t& */ _tunes_out_desc,
      /* Bmad::array_descriptor_t& */ _U_desc
  );
  vec_to_matrix(_N_vec, N);
  vec_to_matrix(_U_vec, U);
  return MakeN{N, _err_flag, _tunes_out, U};
}
Bmad::MakePbrh Bmad::make_pbrh(FixedArray2D<Real, 6, 6> M, FixedArray1D<Real, 3> abz_tunes) {
  // M: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _M_desc;
  _M_desc.rank = 2;
  double _M_vec[6 * 6];
  _M_desc.data_ptr = _M_vec;
  _M_desc.dims[0] = 6;
  _M_desc.dims[1] = 6;
  matrix_to_vec(M, _M_vec);
  // P: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _P_desc;
  _P_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> P;
  std::complex<double> _P_vec[6 * 6];
  _P_desc.data_ptr = _P_vec;
  _P_desc.dims[0] = 6;
  _P_desc.dims[1] = 6;
  // Bp: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Bp_desc;
  _Bp_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> Bp;
  std::complex<double> _Bp_vec[6 * 6];
  _Bp_desc.data_ptr = _Bp_vec;
  _Bp_desc.dims[0] = 6;
  _Bp_desc.dims[1] = 6;
  // R: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _R_desc;
  _R_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> R;
  std::complex<double> _R_vec[6 * 6];
  _R_desc.data_ptr = _R_vec;
  _R_desc.dims[0] = 6;
  _R_desc.dims[1] = 6;
  // H: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _H_desc;
  _H_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> H;
  std::complex<double> _H_vec[6 * 6];
  _H_desc.data_ptr = _H_vec;
  _H_desc.dims[0] = 6;
  _H_desc.dims[1] = 6;
  // abz_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  _abz_tunes_desc.data_ptr = abz_tunes.data();
  _abz_tunes_desc.dims[0] = abz_tunes.size();
  fortran_make_pbrh(
      /* Bmad::array_descriptor_t& */ _M_desc,
      /* Bmad::array_descriptor_t& */ _P_desc,
      /* Bmad::array_descriptor_t& */ _Bp_desc,
      /* Bmad::array_descriptor_t& */ _R_desc,
      /* Bmad::array_descriptor_t& */ _H_desc,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc
  );
  vec_to_matrix(_P_vec, P);
  vec_to_matrix(_Bp_vec, Bp);
  vec_to_matrix(_R_vec, R);
  vec_to_matrix(_H_vec, H);
  return MakePbrh{P, Bp, R, H};
}
Bmad::MakeSmatFromAbc
Bmad::make_smat_from_abc(FixedArray2D<Real, 6, 6> t6, NormalModesStruct &mode) {
  // t6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_desc;
  _t6_desc.rank = 2;
  double _t6_vec[6 * 6];
  _t6_desc.data_ptr = _t6_vec;
  _t6_desc.dims[0] = 6;
  _t6_desc.dims[1] = 6;
  matrix_to_vec(t6, _t6_vec);
  // sigma_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> sigma_mat;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  bool _err_flag{};
  // Nout: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Nout_desc;
  _Nout_desc.rank = 2;
  FixedArray2D<Real, 6, 6> Nout;
  double _Nout_vec[6 * 6];
  _Nout_desc.data_ptr = _Nout_vec;
  _Nout_desc.dims[0] = 6;
  _Nout_desc.dims[1] = 6;
  fortran_make_smat_from_abc(
      /* Bmad::array_descriptor_t& */ _t6_desc,
      /* void* */ mode.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
      /* bool& */ _err_flag,
      /* Bmad::array_descriptor_t& */ _Nout_desc
  );
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  vec_to_matrix(_Nout_vec, Nout);
  return MakeSmatFromAbc{sigma_mat, _err_flag, Nout};
}
void Bmad::make_unit_mad_map(MadMapStruct &map) {
  fortran_make_unit_mad_map(/* void* */ map.get_fortran_ptr());
}
void Bmad::make_v(
    FixedArray2D<Real, 6, 6> M,
    FixedArray2D<Complex, 6, 6> V,
    FixedArray1D<Real, 3> abz_tunes
) {
  // M: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _M_desc;
  _M_desc.rank = 2;
  double _M_vec[6 * 6];
  _M_desc.data_ptr = _M_vec;
  _M_desc.dims[0] = 6;
  _M_desc.dims[1] = 6;
  matrix_to_vec(M, _M_vec);
  // V: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _V_desc;
  _V_desc.rank = 2;
  std::complex<double> _V_vec[6 * 6];
  _V_desc.data_ptr = _V_vec;
  _V_desc.dims[0] = 6;
  _V_desc.dims[1] = 6;
  matrix_to_vec(V, _V_vec);
  // abz_tunes: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  _abz_tunes_desc.data_ptr = abz_tunes.data();
  _abz_tunes_desc.dims[0] = abz_tunes.size();
  fortran_make_v(
      /* Bmad::array_descriptor_t& */ _M_desc,
      /* Bmad::array_descriptor_t& */ _V_desc,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc
  );
  vec_to_matrix(_M_vec, M);
  vec_to_matrix(_V_vec, V);
}
Bmad::MakeVMats Bmad::make_v_mats(EleStruct &ele) {
  // v_mat: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _v_mat_desc;
  _v_mat_desc.rank = 2;
  FixedArray2D<Real, 4, 4> v_mat;
  double _v_mat_vec[4 * 4];
  _v_mat_desc.data_ptr = _v_mat_vec;
  _v_mat_desc.dims[0] = 4;
  _v_mat_desc.dims[1] = 4;
  // v_inv_mat: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _v_inv_mat_desc;
  _v_inv_mat_desc.rank = 2;
  FixedArray2D<Real, 4, 4> v_inv_mat;
  double _v_inv_mat_vec[4 * 4];
  _v_inv_mat_desc.data_ptr = _v_inv_mat_vec;
  _v_inv_mat_desc.dims[0] = 4;
  _v_inv_mat_desc.dims[1] = 4;
  fortran_make_v_mats(/* void* */ ele.get_fortran_ptr(),
                      /* Bmad::array_descriptor_t& */ _v_mat_desc,
                      /* Bmad::array_descriptor_t& */ _v_inv_mat_desc);
  vec_to_matrix(_v_mat_vec, v_mat);
  vec_to_matrix(_v_inv_mat_vec, v_inv_mat);
  return MakeVMats{v_mat, v_inv_mat};
}
void Bmad::makeup_control_slave(LatStruct &lat, EleStruct &slave, bool err_flag) {
  fortran_makeup_control_slave(/* void* */ lat.get_fortran_ptr(),
                               /* void* */ slave.get_fortran_ptr(),
                               /* bool& */ err_flag);
}
void Bmad::makeup_group_lord(LatStruct &lat, EleStruct &lord, bool err_flag) {
  fortran_makeup_group_lord(/* void* */ lat.get_fortran_ptr(),
                            /* void* */ lord.get_fortran_ptr(),
                            /* bool& */ err_flag);
}
void Bmad::makeup_multipass_slave(LatStruct &lat, EleStruct &slave, bool err_flag) {
  fortran_makeup_multipass_slave(/* void* */ lat.get_fortran_ptr(),
                                 /* void* */ slave.get_fortran_ptr(),
                                 /* bool& */ err_flag);
}
void Bmad::makeup_super_slave(LatStruct &lat, EleStruct &slave, bool err_flag) {
  fortran_makeup_super_slave(/* void* */ lat.get_fortran_ptr(),
                             /* void* */ slave.get_fortran_ptr(),
                             /* bool& */ err_flag);
}
bool Bmad::makeup_super_slave1(
    EleStruct &slave,
    EleStruct &lord,
    double offset,
    LatParamStruct &param,
    bool include_upstream_end,
    bool include_downstream_end
) {
  bool _err_flag{};
  fortran_makeup_super_slave1(/* void* */ slave.get_fortran_ptr(),
                              /* void* */ lord.get_fortran_ptr(),
                              /* double& */ offset,
                              /* void* */ param.get_fortran_ptr(),
                              /* bool& */ include_upstream_end,
                              /* bool& */ include_downstream_end,
                              /* bool& */ _err_flag);
  return _err_flag;
}
SpinOrbitMap1Struct Bmad::map1_inverse(SpinOrbitMap1Struct &map1) {
  SpinOrbitMap1Struct _inv_map1;
  fortran_map1_inverse(/* void* */ map1.get_fortran_ptr(), /* void* */ _inv_map1.get_fortran_ptr());
  return std::move(_inv_map1);
}
SpinOrbitMap1Struct Bmad::map1_make_unit() {
  SpinOrbitMap1Struct _map1;
  fortran_map1_make_unit(/* void* */ _map1.get_fortran_ptr());
  return std::move(_map1);
}
void Bmad::map1_times_map1(
    SpinOrbitMap1Struct &map2,
    SpinOrbitMap1Struct &map1,
    SpinOrbitMap1Struct &map_out
) {
  fortran_map1_times_map1(/* void* */ map2.get_fortran_ptr(),
                          /* void* */ map1.get_fortran_ptr(),
                          /* void* */ map_out.get_fortran_ptr());
}
TaylorStructArray1D Bmad::map_to_angle_coords(TaylorStructArray1D t_canon) {
  // t_canon: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _t_canon_desc;
  _t_canon_desc.rank = 1;
  _t_canon_desc.data_ptr = t_canon.data();
  _t_canon_desc.dims[0] = t_canon.size();
  _t_canon_desc.strides[0] = 1;
  // t_angle: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _t_angle_desc;
  _t_angle_desc.rank = 1;
  // Output-only type array
  auto t_angle = TaylorStructAlloc1D(6);
  _t_angle_desc.data_ptr = t_angle.get_fortran_ptr();
  _t_angle_desc.dims[0] = t_angle.size();

  _t_angle_desc.strides[0] = 1;
  fortran_map_to_angle_coords(
      /* Bmad::array_descriptor_t& */ _t_canon_desc,
      /* Bmad::array_descriptor_t& */ _t_angle_desc
  );
  return std::move(std::move(t_angle));
}
void Bmad::mark_patch_regions(BranchStruct &branch) {
  fortran_mark_patch_regions(/* void* */ branch.get_fortran_ptr());
}
double Bmad::master_parameter_value(int master_parameter, EleStruct &ele) {
  double _value{};
  fortran_master_parameter_value(
      /* int& */ master_parameter,
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _value
  );
  return _value;
}
FixedArray2D<Real, 4, 4> Bmad::mat4_multipole(double knl, double tilt, int n, CoordStruct &orbit) {
  // kick_mat: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _kick_mat_desc;
  _kick_mat_desc.rank = 2;
  FixedArray2D<Real, 4, 4> kick_mat;
  double _kick_mat_vec[4 * 4];
  _kick_mat_desc.data_ptr = _kick_mat_vec;
  _kick_mat_desc.dims[0] = 4;
  _kick_mat_desc.dims[1] = 4;
  fortran_mat4_multipole(
      /* double& */ knl,
      /* double& */ tilt,
      /* int& */ n,
      /* void* */ orbit.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _kick_mat_desc
  );
  vec_to_matrix(_kick_mat_vec, kick_mat);
  return kick_mat;
}
void Bmad::mat6_add_offsets(EleStruct &ele, LatParamStruct &param) {
  fortran_mat6_add_offsets(/* void* */ ele.get_fortran_ptr(), /* void* */ param.get_fortran_ptr());
}
void Bmad::mat6_add_pitch(
    double x_pitch_tot,
    double y_pitch_tot,
    int orientation,
    FixedArray2D<Real, 6, 6> mat6
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  matrix_to_vec(mat6, _mat6_vec);
  fortran_mat6_add_pitch(
      /* double& */ x_pitch_tot,
      /* double& */ y_pitch_tot,
      /* int& */ orientation,
      /* Bmad::array_descriptor_t& */ _mat6_desc
  );
  vec_to_matrix(_mat6_vec, mat6);
}
ComplexTaylorStructArray1D
Bmad::mat6_to_complex_taylor(FixedArray1D<Complex, 6> vec0, FixedArray2D<Complex, 6, 6> mat6) {
  // vec0: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec0_desc;
  _vec0_desc.rank = 1;
  _vec0_desc.data_ptr = vec0.data();
  _vec0_desc.dims[0] = vec0.size();
  // mat6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  std::complex<double> _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  matrix_to_vec(mat6, _mat6_vec);
  // complex_taylor: ComplexTaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_desc;
  _complex_taylor_desc.rank = 1;
  // Output-only type array
  auto complex_taylor = ComplexTaylorStructAlloc1D(6);
  _complex_taylor_desc.data_ptr = complex_taylor.get_fortran_ptr();
  _complex_taylor_desc.dims[0] = complex_taylor.size();

  _complex_taylor_desc.strides[0] = 1;
  fortran_mat6_to_complex_taylor(
      /* Bmad::array_descriptor_t& */ _vec0_desc,
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* Bmad::array_descriptor_t& */ _complex_taylor_desc
  );
  return std::move(std::move(complex_taylor));
}
Bmad::MatSympDecouple Bmad::mat_symp_decouple(FixedArray2D<Real, 4, 4> t0, bool type_out) {
  // t0: in NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _t0_desc;
  _t0_desc.rank = 2;
  double _t0_vec[4 * 4];
  _t0_desc.data_ptr = _t0_vec;
  _t0_desc.dims[0] = 4;
  _t0_desc.dims[1] = 4;
  matrix_to_vec(t0, _t0_vec);
  int _stat{};
  // U: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _U_desc;
  _U_desc.rank = 2;
  FixedArray2D<Real, 4, 4> U;
  double _U_vec[4 * 4];
  _U_desc.data_ptr = _U_vec;
  _U_desc.dims[0] = 4;
  _U_desc.dims[1] = 4;
  // V: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _V_desc;
  _V_desc.rank = 2;
  FixedArray2D<Real, 4, 4> V;
  double _V_vec[4 * 4];
  _V_desc.data_ptr = _V_vec;
  _V_desc.dims[0] = 4;
  _V_desc.dims[1] = 4;
  // Ubar: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _Ubar_desc;
  _Ubar_desc.rank = 2;
  FixedArray2D<Real, 4, 4> Ubar;
  double _Ubar_vec[4 * 4];
  _Ubar_desc.data_ptr = _Ubar_vec;
  _Ubar_desc.dims[0] = 4;
  _Ubar_desc.dims[1] = 4;
  // Vbar: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _Vbar_desc;
  _Vbar_desc.rank = 2;
  FixedArray2D<Real, 4, 4> Vbar;
  double _Vbar_vec[4 * 4];
  _Vbar_desc.data_ptr = _Vbar_vec;
  _Vbar_desc.dims[0] = 4;
  _Vbar_desc.dims[1] = 4;
  // G: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _G_desc;
  _G_desc.rank = 2;
  FixedArray2D<Real, 4, 4> G;
  double _G_vec[4 * 4];
  _G_desc.data_ptr = _G_vec;
  _G_desc.dims[0] = 4;
  _G_desc.dims[1] = 4;
  TwissStruct _twiss1;
  TwissStruct _twiss2;
  double _gamma{};
  fortran_mat_symp_decouple(
      /* Bmad::array_descriptor_t& */ _t0_desc,
      /* int& */ _stat,
      /* Bmad::array_descriptor_t& */ _U_desc,
      /* Bmad::array_descriptor_t& */ _V_desc,
      /* Bmad::array_descriptor_t& */ _Ubar_desc,
      /* Bmad::array_descriptor_t& */ _Vbar_desc,
      /* Bmad::array_descriptor_t& */ _G_desc,
      /* void* */ _twiss1.get_fortran_ptr(),
      /* void* */ _twiss2.get_fortran_ptr(),
      /* double& */ _gamma,
      /* bool& */ type_out
  );
  vec_to_matrix(_U_vec, U);
  vec_to_matrix(_V_vec, V);
  vec_to_matrix(_Ubar_vec, Ubar);
  vec_to_matrix(_Vbar_vec, Vbar);
  vec_to_matrix(_G_vec, G);
  return MatSympDecouple{
      _stat,
      U,
      V,
      Ubar,
      Vbar,
      G,
      std::move(_twiss1),
      std::move(_twiss2),
      _gamma
  };
}
Bmad::MatchEleToMat6 Bmad::match_ele_to_mat6(
    EleStruct &ele,
    CoordStruct &start_orb,
    std::optional<bool> include_delta_time,
    std::optional<bool> set_trombone
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  // vec0: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _vec0_desc;
  _vec0_desc.rank = 1;
  FixedArray1D<Real, 6> _vec0;
  _vec0_desc.data_ptr = _vec0.data();
  _vec0_desc.dims[0] = _vec0.size();
  bool _err_flag{};
  bool include_delta_time_lvalue;
  auto *_include_delta_time{&include_delta_time_lvalue};
  if (include_delta_time.has_value()) {
    include_delta_time_lvalue = include_delta_time.value();
  } else {
    _include_delta_time = nullptr;
  }
  bool set_trombone_lvalue;
  auto *_set_trombone{&set_trombone_lvalue};
  if (set_trombone.has_value()) {
    set_trombone_lvalue = set_trombone.value();
  } else {
    _set_trombone = nullptr;
  }
  fortran_match_ele_to_mat6(/* void* */ ele.get_fortran_ptr(),
                            /* void* */ start_orb.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _mat6_desc,
                            /* Bmad::array_descriptor_t& */ _vec0_desc,
                            /* bool& */ _err_flag,
                            /* bool* */ _include_delta_time,
                            /* bool* */ _set_trombone);
  vec_to_matrix(_mat6_vec, mat6);
  return MatchEleToMat6{mat6, _vec0, _err_flag};
}
double Bmad::mexp(double x, int m) {
  double _this_exp{};
  fortran_mexp(/* double& */ x, /* int& */ m, /* double& */ _this_exp);
  return _this_exp;
}
int Bmad::mfft1(FArray1D<Real> &a, FArray1D<Real> &b, FArray1D<Int> &n, int ndim, int isn) {
  // a: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _a_desc;
  _a_desc.rank = 1;
  _a_desc.data_ptr = a.data();
  _a_desc.dims[0] = a.size();
  // b: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  _b_desc.data_ptr = b.data();
  _b_desc.dims[0] = b.size();
  // n: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _n_desc;
  _n_desc.rank = 1;
  _n_desc.data_ptr = n.data();
  _n_desc.dims[0] = n.size();
  int _ierr{};
  fortran_mfft1(
      /* Bmad::array_descriptor_t& */ _a_desc,
      /* Bmad::array_descriptor_t& */ _b_desc,
      /* Bmad::array_descriptor_t& */ _n_desc,
      /* int& */ ndim,
      /* int& */ isn,
      /* int& */ _ierr
  );
  return _ierr;
}
std::optional<Fibre> Bmad::misalign_ptc_fibre(EleStruct &ele, bool use_offsets, bool for_layout) {
  void *_ptc_fibre;
  fortran_misalign_ptc_fibre(/* void* */ ele.get_fortran_ptr(),
                             /* bool& */ use_offsets,
                             /* void* */ &_ptc_fibre,
                             /* bool& */ for_layout);
  return std::move((_ptc_fibre ? std::make_optional<Fibre>(_ptc_fibre) : std::nullopt));
}
double Bmad::momentum_compaction(BranchStruct &branch) {
  double _mom_comp{};
  fortran_momentum_compaction(/* void* */ branch.get_fortran_ptr(), /* double& */ _mom_comp);
  return _mom_comp;
}
Bmad::MultiTurnTrackingAnalysis
Bmad::multi_turn_tracking_analysis(CoordStructArray1D track, int i_dim) {
  // track: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _track_desc;
  _track_desc.rank = 1;
  _track_desc.data_ptr = track.data();
  _track_desc.dims[0] = track.size();
  _track_desc.strides[0] = 1;
  CoordStruct _track0;
  EleStruct _ele;
  bool _stable{};
  double _growth_rate{};
  double _chi{};
  bool _err_flag{};
  fortran_multi_turn_tracking_analysis(
      /* Bmad::array_descriptor_t& */ _track_desc,
      /* int& */ i_dim,
      /* void* */ _track0.get_fortran_ptr(),
      /* void* */ _ele.get_fortran_ptr(),
      /* bool& */ _stable,
      /* double& */ _growth_rate,
      /* double& */ _chi,
      /* bool& */ _err_flag
  );
  return MultiTurnTrackingAnalysis{
      std::move(_track0),
      std::move(_ele),
      _stable,
      _growth_rate,
      _chi,
      _err_flag
  };
}
bool Bmad::multilayer_type_to_multilayer_params(EleStruct &ele) {
  bool _err_flag{};
  fortran_multilayer_type_to_multilayer_params(/* void* */ ele.get_fortran_ptr(),
                                               /* bool& */ _err_flag);
  return _err_flag;
}
Bmad::MultipassChain Bmad::multipass_chain(EleStruct &ele, std::optional<bool> use_super_lord) {
  int _ix_pass{};
  int _n_links{};
  // intent=out allocatable type array
  auto chain_ele{ElePointerStructAlloc1D()};
  bool use_super_lord_lvalue;
  auto *_use_super_lord{&use_super_lord_lvalue};
  if (use_super_lord.has_value()) {
    use_super_lord_lvalue = use_super_lord.value();
  } else {
    _use_super_lord = nullptr;
  }
  fortran_multipass_chain(/* void* */ ele.get_fortran_ptr(),
                          /* int& */ _ix_pass,
                          /* int& */ _n_links,
                          /* void* */ chain_ele.get_fortran_ptr(),
                          /* bool* */ _use_super_lord);
  return MultipassChain{_ix_pass, _n_links, std::move(chain_ele)};
}
Bmad::Multipole1AbToKt Bmad::multipole1_ab_to_kt(double an, double bn, int n) {
  double _knl{};
  double _tn{};
  fortran_multipole1_ab_to_kt(
      /* double& */ an,
      /* double& */ bn,
      /* int& */ n,
      /* double& */ _knl,
      /* double& */ _tn
  );
  return Multipole1AbToKt{_knl, _tn};
}
Bmad::Multipole1KtToAb Bmad::multipole1_kt_to_ab(double knl, double knsl, double tn, int n) {
  double _an{};
  double _bn{};
  fortran_multipole1_kt_to_ab(
      /* double& */ knl,
      /* double& */ knsl,
      /* double& */ tn,
      /* int& */ n,
      /* double& */ _an,
      /* double& */ _bn
  );
  return Multipole1KtToAb{_an, _bn};
}
void Bmad::multipole_ab_to_kt(
    FArray1D<Real> &an,
    FArray1D<Real> &bn,
    FArray1D<Real> &knl,
    FArray1D<Real> &tn
) {
  // an: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _an_desc;
  _an_desc.rank = 1;
  _an_desc.data_ptr = an.data();
  _an_desc.dims[0] = an.size();
  // bn: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _bn_desc;
  _bn_desc.rank = 1;
  _bn_desc.data_ptr = bn.data();
  _bn_desc.dims[0] = bn.size();
  // knl: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knl_desc;
  _knl_desc.rank = 1;
  _knl_desc.data_ptr = knl.data();
  _knl_desc.dims[0] = knl.size();
  // tn: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _tn_desc;
  _tn_desc.rank = 1;
  _tn_desc.data_ptr = tn.data();
  _tn_desc.dims[0] = tn.size();
  fortran_multipole_ab_to_kt(
      /* Bmad::array_descriptor_t& */ _an_desc,
      /* Bmad::array_descriptor_t& */ _bn_desc,
      /* Bmad::array_descriptor_t& */ _knl_desc,
      /* Bmad::array_descriptor_t& */ _tn_desc
  );
}
Bmad::MultipoleEleToAb Bmad::multipole_ele_to_ab(
    EleStruct &ele,
    bool use_ele_tilt,
    std::optional<int> pole_type,
    std::optional<int> include_kicks,
    std::optional<bool> original
) {
  int _ix_pole_max{};
  // a: out NOT (CppWrapperGeneralArgumentArray) (['0:n_pole_maxx'])
  Bmad::array_descriptor_t _a_desc;
  _a_desc.rank = 1;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> _a;
  _a_desc.data_ptr = _a.data();
  _a_desc.dims[0] = _a.size();
  // b: out NOT (CppWrapperGeneralArgumentArray) (['0:n_pole_maxx'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> _b;
  _b_desc.data_ptr = _b.data();
  _b_desc.dims[0] = _b.size();
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  int include_kicks_lvalue;
  auto *_include_kicks{&include_kicks_lvalue};
  if (include_kicks.has_value()) {
    include_kicks_lvalue = include_kicks.value();
  } else {
    _include_kicks = nullptr;
  }
  double _b1{};
  bool original_lvalue;
  auto *_original{&original_lvalue};
  if (original.has_value()) {
    original_lvalue = original.value();
  } else {
    _original = nullptr;
  }
  fortran_multipole_ele_to_ab(/* void* */ ele.get_fortran_ptr(),
                              /* bool& */ use_ele_tilt,
                              /* int& */ _ix_pole_max,
                              /* Bmad::array_descriptor_t& */ _a_desc,
                              /* Bmad::array_descriptor_t& */ _b_desc,
                              /* int* */ _pole_type,
                              /* int* */ _include_kicks,
                              /* double& */ _b1,
                              /* bool* */ _original);
  return MultipoleEleToAb{_ix_pole_max, _a, _b, _b1};
}
int Bmad::multipole_ele_to_kt(
    EleStruct &ele,
    bool use_ele_tilt,
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    std::optional<int> pole_type,
    std::optional<int> include_kicks
) {
  int _ix_pole_max{};
  // knl: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knl_desc;
  _knl_desc.rank = 1;
  _knl_desc.data_ptr = knl.data();
  _knl_desc.dims[0] = knl.size();
  // tilt: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _tilt_desc;
  _tilt_desc.rank = 1;
  _tilt_desc.data_ptr = tilt.data();
  _tilt_desc.dims[0] = tilt.size();
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  int include_kicks_lvalue;
  auto *_include_kicks{&include_kicks_lvalue};
  if (include_kicks.has_value()) {
    include_kicks_lvalue = include_kicks.value();
  } else {
    _include_kicks = nullptr;
  }
  fortran_multipole_ele_to_kt(/* void* */ ele.get_fortran_ptr(),
                              /* bool& */ use_ele_tilt,
                              /* int& */ _ix_pole_max,
                              /* Bmad::array_descriptor_t& */ _knl_desc,
                              /* Bmad::array_descriptor_t& */ _tilt_desc,
                              /* int* */ _pole_type,
                              /* int* */ _include_kicks);
  return _ix_pole_max;
}
EleStruct Bmad::multipole_init(int who, std::optional<bool> zero) {
  EleStruct _ele;
  bool zero_lvalue;
  auto *_zero{&zero_lvalue};
  if (zero.has_value()) {
    zero_lvalue = zero.value();
  } else {
    _zero = nullptr;
  }
  fortran_multipole_init(/* void* */ _ele.get_fortran_ptr(), /* int& */ who, /* bool* */ _zero);
  return std::move(_ele);
}
void Bmad::multipole_kick(
    double knl,
    double tilt,
    int n,
    int ref_species,
    int ele_orientation,
    CoordStruct &coord,
    std::optional<int> pole_type,
    std::optional<bool> ref_orb_offset
) {
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  bool ref_orb_offset_lvalue;
  auto *_ref_orb_offset{&ref_orb_offset_lvalue};
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
      /* bool* */ _ref_orb_offset
  );
}
FixedArray2D<Real, 6, 6> Bmad::multipole_kick_mat(
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    int ref_species,
    EleStruct &ele,
    CoordStruct &orbit,
    double factor
) {
  // knl: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knl_desc;
  _knl_desc.rank = 1;
  _knl_desc.data_ptr = knl.data();
  _knl_desc.dims[0] = knl.size();
  // tilt: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _tilt_desc;
  _tilt_desc.rank = 1;
  _tilt_desc.data_ptr = tilt.data();
  _tilt_desc.dims[0] = tilt.size();
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  fortran_multipole_kick_mat(
      /* Bmad::array_descriptor_t& */ _knl_desc,
      /* Bmad::array_descriptor_t& */ _tilt_desc,
      /* int& */ ref_species,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* double& */ factor,
      /* Bmad::array_descriptor_t& */ _mat6_desc
  );
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::multipole_kicks(
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<int> pole_type,
    std::optional<bool> ref_orb_offset
) {
  // knl: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knl_desc;
  _knl_desc.rank = 1;
  _knl_desc.data_ptr = knl.data();
  _knl_desc.dims[0] = knl.size();
  // tilt: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _tilt_desc;
  _tilt_desc.rank = 1;
  _tilt_desc.data_ptr = tilt.data();
  _tilt_desc.dims[0] = tilt.size();
  int pole_type_lvalue;
  auto *_pole_type{&pole_type_lvalue};
  if (pole_type.has_value()) {
    pole_type_lvalue = pole_type.value();
  } else {
    _pole_type = nullptr;
  }
  bool ref_orb_offset_lvalue;
  auto *_ref_orb_offset{&ref_orb_offset_lvalue};
  if (ref_orb_offset.has_value()) {
    ref_orb_offset_lvalue = ref_orb_offset.value();
  } else {
    _ref_orb_offset = nullptr;
  }
  fortran_multipole_kicks(
      /* Bmad::array_descriptor_t& */ _knl_desc,
      /* Bmad::array_descriptor_t& */ _tilt_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* int* */ _pole_type,
      /* bool* */ _ref_orb_offset
  );
}
void Bmad::multipole_kt_to_ab(
    FArray1D<Real> &knl,
    FArray1D<Real> &knsl,
    FArray1D<Real> &tn,
    FArray1D<Real> &an,
    FArray1D<Real> &bn
) {
  // knl: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knl_desc;
  _knl_desc.rank = 1;
  _knl_desc.data_ptr = knl.data();
  _knl_desc.dims[0] = knl.size();
  // knsl: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _knsl_desc;
  _knsl_desc.rank = 1;
  _knsl_desc.data_ptr = knsl.data();
  _knsl_desc.dims[0] = knsl.size();
  // tn: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _tn_desc;
  _tn_desc.rank = 1;
  _tn_desc.data_ptr = tn.data();
  _tn_desc.dims[0] = tn.size();
  // an: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _an_desc;
  _an_desc.rank = 1;
  _an_desc.data_ptr = an.data();
  _an_desc.dims[0] = an.size();
  // bn: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _bn_desc;
  _bn_desc.rank = 1;
  _bn_desc.data_ptr = bn.data();
  _bn_desc.dims[0] = bn.size();
  fortran_multipole_kt_to_ab(
      /* Bmad::array_descriptor_t& */ _knl_desc,
      /* Bmad::array_descriptor_t& */ _knsl_desc,
      /* Bmad::array_descriptor_t& */ _tn_desc,
      /* Bmad::array_descriptor_t& */ _an_desc,
      /* Bmad::array_descriptor_t& */ _bn_desc
  );
}
void Bmad::multipole_spin_tracking(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_multipole_spin_tracking(/* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* void* */ orbit.get_fortran_ptr());
}
void Bmad::mytan(double y, double x, double arg) {
  fortran_mytan(/* double& */ y, /* double& */ x, /* double& */ arg);
}
int Bmad::n_attrib_string_max_len() {
  int _max_len{};
  fortran_n_attrib_string_max_len(/* int& */ _max_len);
  return _max_len;
}
int Bmad::new_control(LatStruct &lat, std::optional<std::string> ele_name) {
  int _ix_ele{};
  const char *_ele_name = ele_name.has_value() ? ele_name->c_str() : nullptr;
  fortran_new_control(/* void* */ lat.get_fortran_ptr(),
                      /* int& */ _ix_ele,
                      /* const char* */ _ele_name);
  return _ix_ele;
}
int Bmad::nint_chk(double re_val) {
  int _int_val{};
  fortran_nint_chk(/* double& */ re_val, /* int& */ _int_val);
  return _int_val;
}
void Bmad::normal_form_complex_taylors(
    TaylorStructArray1D one_turn_taylor,
    bool rf_on,
    std::optional<ComplexTaylorStructArray1D> F,
    std::optional<ComplexTaylorStructArray1D> L,
    std::optional<TaylorStructArray1D> A,
    std::optional<TaylorStructArray1D> A_inverse,
    std::optional<int> order
) {
  // one_turn_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _one_turn_taylor_desc;
  _one_turn_taylor_desc.rank = 1;
  _one_turn_taylor_desc.data_ptr = one_turn_taylor.data();
  _one_turn_taylor_desc.dims[0] = one_turn_taylor.size();
  _one_turn_taylor_desc.strides[0] = 1;
  // F: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _F_desc;
  _F_desc.rank = 1;
  if (F) {
    _F_desc.data_ptr = F->data();
    _F_desc.dims[0] = F->size();
  } else {
    _F_desc.data_ptr = nullptr;
    _F_desc.dims[0] = 0;
  }
  _F_desc.strides[0] = 1;
  // L: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _L_desc;
  _L_desc.rank = 1;
  if (L) {
    _L_desc.data_ptr = L->data();
    _L_desc.dims[0] = L->size();
  } else {
    _L_desc.data_ptr = nullptr;
    _L_desc.dims[0] = 0;
  }
  _L_desc.strides[0] = 1;
  // A: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _A_desc;
  _A_desc.rank = 1;
  if (A) {
    _A_desc.data_ptr = A->data();
    _A_desc.dims[0] = A->size();
  } else {
    _A_desc.data_ptr = nullptr;
    _A_desc.dims[0] = 0;
  }
  _A_desc.strides[0] = 1;
  // A_inverse: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _A_inverse_desc;
  _A_inverse_desc.rank = 1;
  if (A_inverse) {
    _A_inverse_desc.data_ptr = A_inverse->data();
    _A_inverse_desc.dims[0] = A_inverse->size();
  } else {
    _A_inverse_desc.data_ptr = nullptr;
    _A_inverse_desc.dims[0] = 0;
  }
  _A_inverse_desc.strides[0] = 1;
  int order_lvalue;
  auto *_order{&order_lvalue};
  if (order.has_value()) {
    order_lvalue = order.value();
  } else {
    _order = nullptr;
  }
  fortran_normal_form_complex_taylors(
      /* Bmad::array_descriptor_t& */ _one_turn_taylor_desc,
      /* bool& */ rf_on,
      /* Bmad::array_descriptor_t& */ _F_desc,
      /* Bmad::array_descriptor_t& */ _L_desc,
      /* Bmad::array_descriptor_t& */ _A_desc,
      /* Bmad::array_descriptor_t& */ _A_inverse_desc,
      /* int* */ _order
  );
}
Bmad::NormalFormTaylors Bmad::normal_form_taylors(TaylorStructArray1D one_turn_taylor, bool rf_on) {
  // one_turn_taylor: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _one_turn_taylor_desc;
  _one_turn_taylor_desc.rank = 1;
  _one_turn_taylor_desc.data_ptr = one_turn_taylor.data();
  _one_turn_taylor_desc.dims[0] = one_turn_taylor.size();
  _one_turn_taylor_desc.strides[0] = 1;
  // dhdj: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _dhdj_desc;
  _dhdj_desc.rank = 1;
  // Output-only type array
  auto dhdj = TaylorStructAlloc1D(6);
  _dhdj_desc.data_ptr = dhdj.get_fortran_ptr();
  _dhdj_desc.dims[0] = dhdj.size();

  _dhdj_desc.strides[0] = 1;
  // A: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _A_desc;
  _A_desc.rank = 1;
  // Output-only type array
  auto A = TaylorStructAlloc1D(6);
  _A_desc.data_ptr = A.get_fortran_ptr();
  _A_desc.dims[0] = A.size();

  _A_desc.strides[0] = 1;
  // A_inverse: TaylorStruct out (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _A_inverse_desc;
  _A_inverse_desc.rank = 1;
  // Output-only type array
  auto A_inverse = TaylorStructAlloc1D(6);
  _A_inverse_desc.data_ptr = A_inverse.get_fortran_ptr();
  _A_inverse_desc.dims[0] = A_inverse.size();

  _A_inverse_desc.strides[0] = 1;
  fortran_normal_form_taylors(
      /* Bmad::array_descriptor_t& */ _one_turn_taylor_desc,
      /* bool& */ rf_on,
      /* Bmad::array_descriptor_t& */ _dhdj_desc,
      /* Bmad::array_descriptor_t& */ _A_desc,
      /* Bmad::array_descriptor_t& */ _A_inverse_desc
  );
  return NormalFormTaylors{
      std::move(std::move(dhdj)),
      std::move(std::move(A)),
      std::move(std::move(A_inverse))
  };
}
Bmad::NormalMode3Calc Bmad::normal_mode3_calc(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<bool> above_transition,
    std::optional<FixedArray1D<Real, 3>> abz_tunes
) {
  // t6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_desc;
  _t6_desc.rank = 2;
  double _t6_vec[6 * 6];
  _t6_desc.data_ptr = _t6_vec;
  _t6_desc.dims[0] = 6;
  _t6_desc.dims[1] = 6;
  matrix_to_vec(t6, _t6_vec);
  // tune: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _tune_desc;
  _tune_desc.rank = 1;
  FixedArray1D<Real, 3> _tune;
  _tune_desc.data_ptr = _tune.data();
  _tune_desc.dims[0] = _tune.size();
  // B: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _B_desc;
  _B_desc.rank = 2;
  FixedArray2D<Real, 6, 6> B;
  double _B_vec[6 * 6];
  _B_desc.data_ptr = _B_vec;
  _B_desc.dims[0] = 6;
  _B_desc.dims[1] = 6;
  // HV: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _HV_desc;
  _HV_desc.rank = 2;
  FixedArray2D<Real, 6, 6> HV;
  double _HV_vec[6 * 6];
  _HV_desc.data_ptr = _HV_vec;
  _HV_desc.dims[0] = 6;
  _HV_desc.dims[1] = 6;
  bool above_transition_lvalue;
  auto *_above_transition{&above_transition_lvalue};
  if (above_transition.has_value()) {
    above_transition_lvalue = above_transition.value();
  } else {
    _above_transition = nullptr;
  }
  // abz_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  if (abz_tunes.has_value()) {
    _abz_tunes_desc.data_ptr = abz_tunes->data();
    _abz_tunes_desc.dims[0] = abz_tunes->size();
  } else {
    _abz_tunes_desc.data_ptr = nullptr;
    _abz_tunes_desc.dims[0] = 0;
  }
  fortran_normal_mode3_calc(
      /* Bmad::array_descriptor_t& */ _t6_desc,
      /* Bmad::array_descriptor_t& */ _tune_desc,
      /* Bmad::array_descriptor_t& */ _B_desc,
      /* Bmad::array_descriptor_t& */ _HV_desc,
      /* bool* */ _above_transition,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc
  );
  vec_to_matrix(_t6_vec, t6);
  vec_to_matrix(_B_vec, B);
  vec_to_matrix(_HV_vec, HV);
  return NormalMode3Calc{_tune, B, HV};
}
void Bmad::normal_mode_dispersion(EleStruct &ele, std::optional<bool> reverse) {
  bool reverse_lvalue;
  auto *_reverse{&reverse_lvalue};
  if (reverse.has_value()) {
    reverse_lvalue = reverse.value();
  } else {
    _reverse = nullptr;
  }
  fortran_normal_mode_dispersion(/* void* */ ele.get_fortran_ptr(), /* bool* */ _reverse);
}
bool Bmad::normalize_evecs(FixedArray2D<Complex, 6, 6> evec) {
  // evec: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _evec_desc;
  _evec_desc.rank = 2;
  std::complex<double> _evec_vec[6 * 6];
  _evec_desc.data_ptr = _evec_vec;
  _evec_desc.dims[0] = 6;
  _evec_desc.dims[1] = 6;
  matrix_to_vec(evec, _evec_vec);
  bool _err_flag{};
  fortran_normalize_evecs(/* Bmad::array_descriptor_t& */ _evec_desc, /* bool& */ _err_flag);
  vec_to_matrix(_evec_vec, evec);
  return _err_flag;
}
int Bmad::num_field_eles(EleStruct &ele) {
  int _n_field_ele{};
  fortran_num_field_eles(/* void* */ ele.get_fortran_ptr(), /* int& */ _n_field_ele);
  return _n_field_ele;
}
int Bmad::num_lords(EleStruct &slave, int lord_type) {
  int _num{};
  fortran_num_lords(/* void* */ slave.get_fortran_ptr(), /* int& */ lord_type, /* int& */ _num);
  return _num;
}
Bmad::OdeintBmad Bmad::odeint_bmad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    double s1_body,
    double s2_body,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  bool _err_flag{};
  TrackStruct _track;
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_odeint_bmad(/* void* */ orbit.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* void* */ param.get_fortran_ptr(),
                      /* double& */ s1_body,
                      /* double& */ s2_body,
                      /* bool& */ _err_flag,
                      /* void* */ _track.get_fortran_ptr(),
                      /* Bmad::array_descriptor_t& */ _mat6_desc,
                      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return OdeintBmad{_err_flag, std::move(_track)};
}
Bmad::OdeintBmadTime Bmad::odeint_bmad_time(
    CoordStruct &orb,
    EleStruct &ele,
    LatParamStruct &param,
    int t_dir,
    double &rf_time,
    optional_ref<TrackStruct> track,
    std::optional<double> t_end,
    optional_ref<EmFieldStruct> extra_field
) {
  bool _err_flag{};
  auto *_track = track.has_value() ? track->get().get_fortran_ptr() : nullptr; // input, optional
  double t_end_lvalue;
  auto *_t_end{&t_end_lvalue};
  if (t_end.has_value()) {
    t_end_lvalue = t_end.value();
  } else {
    _t_end = nullptr;
  }
  double _dt_step{};
  auto *_extra_field =
      extra_field.has_value() ? extra_field->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_odeint_bmad_time(/* void* */ orb.get_fortran_ptr(),
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
    EleStruct &ele,
    bool set,
    CoordStruct &orbit,
    std::optional<bool> set_tilt,
    std::optional<bool> set_hvkicks,
    std::optional<int> drift_to_edge,
    std::optional<double> s_pos,
    std::optional<bool> set_spin,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    optional_ref<double> time
) {
  bool set_tilt_lvalue;
  auto *_set_tilt{&set_tilt_lvalue};
  if (set_tilt.has_value()) {
    set_tilt_lvalue = set_tilt.value();
  } else {
    _set_tilt = nullptr;
  }
  bool set_hvkicks_lvalue;
  auto *_set_hvkicks{&set_hvkicks_lvalue};
  if (set_hvkicks.has_value()) {
    set_hvkicks_lvalue = set_hvkicks.value();
  } else {
    _set_hvkicks = nullptr;
  }
  int drift_to_edge_lvalue;
  auto *_drift_to_edge{&drift_to_edge_lvalue};
  if (drift_to_edge.has_value()) {
    drift_to_edge_lvalue = drift_to_edge.value();
  } else {
    _drift_to_edge = nullptr;
  }
  double s_pos_lvalue;
  auto *_s_pos{&s_pos_lvalue};
  if (s_pos.has_value()) {
    s_pos_lvalue = s_pos.value();
  } else {
    _s_pos = nullptr;
  }
  double _s_out{};
  bool set_spin_lvalue;
  auto *_set_spin{&set_spin_lvalue};
  if (set_spin.has_value()) {
    set_spin_lvalue = set_spin.value();
  } else {
    _set_spin = nullptr;
  }
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  // spin_qrot: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _spin_qrot_desc;
  _spin_qrot_desc.rank = 1;
  FixedArray1D<Real, 4> _spin_qrot;
  _spin_qrot_desc.data_ptr = _spin_qrot.data();
  _spin_qrot_desc.dims[0] = _spin_qrot.size();
  auto *_time = time.has_value() ? &time->get() : nullptr; // inout, optional
  fortran_offset_particle(/* void* */ ele.get_fortran_ptr(),
                          /* bool& */ set,
                          /* void* */ orbit.get_fortran_ptr(),
                          /* bool* */ _set_tilt,
                          /* bool* */ _set_hvkicks,
                          /* int* */ _drift_to_edge,
                          /* double* */ _s_pos,
                          /* double& */ _s_out,
                          /* bool* */ _set_spin,
                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                          /* bool* */ _make_matrix,
                          /* Bmad::array_descriptor_t& */ _spin_qrot_desc,
                          /* double* */ _time);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return OffsetParticle{_s_out, _spin_qrot};
}
void Bmad::offset_photon(
    EleStruct &ele,
    CoordStruct &orbit,
    bool set,
    std::optional<bool> offset_position_only,
    std::optional<FixedArray2D<Real, 3, 3>> rot_mat
) {
  bool offset_position_only_lvalue;
  auto *_offset_position_only{&offset_position_only_lvalue};
  if (offset_position_only.has_value()) {
    offset_position_only_lvalue = offset_position_only.value();
  } else {
    _offset_position_only = nullptr;
  }
  // rot_mat: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _rot_mat_desc;
  _rot_mat_desc.rank = 2;
  double _rot_mat_vec[3 * 3];
  _rot_mat_desc.data_ptr = _rot_mat_vec;
  if (rot_mat.has_value()) {
    matrix_to_vec(rot_mat.value(), _rot_mat_vec);
    _rot_mat_desc.dims[0] = 3;
    _rot_mat_desc.dims[1] = 3;
  } else {
    _rot_mat_desc.data_ptr = nullptr;
  }
  fortran_offset_photon(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr(),
                        /* bool& */ set,
                        /* bool* */ _offset_position_only,
                        /* Bmad::array_descriptor_t& */ _rot_mat_desc);
}
FixedArray2D<Real, 4, 4> Bmad::one_turn_mat_at_ele(EleStruct &ele, double phi_a, double phi_b) {
  // mat4: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _mat4_desc;
  _mat4_desc.rank = 2;
  FixedArray2D<Real, 4, 4> mat4;
  double _mat4_vec[4 * 4];
  _mat4_desc.data_ptr = _mat4_vec;
  _mat4_desc.dims[0] = 4;
  _mat4_desc.dims[1] = 4;
  fortran_one_turn_mat_at_ele(/* void* */ ele.get_fortran_ptr(),
                              /* double& */ phi_a,
                              /* double& */ phi_b,
                              /* Bmad::array_descriptor_t& */ _mat4_desc);
  vec_to_matrix(_mat4_vec, mat4);
  return mat4;
}
Bmad::OpenBinaryFile
Bmad::open_binary_file(std::string file_name, std::string action, std::string r_name) {
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
      /* bool& */ _is_ok
  );
  return OpenBinaryFile{_iu, _iver, _is_ok};
}
Bmad::OrbitAmplitudeCalc Bmad::orbit_amplitude_calc(EleStruct &ele, CoordStruct &orb) {
  double _amp_a{};
  double _amp_b{};
  double _amp_na{};
  double _amp_nb{};
  fortran_orbit_amplitude_calc(/* void* */ ele.get_fortran_ptr(),
                               /* void* */ orb.get_fortran_ptr(),
                               /* double& */ _amp_a,
                               /* double& */ _amp_b,
                               /* double& */ _amp_na,
                               /* double& */ _amp_nb);
  return OrbitAmplitudeCalc{_amp_a, _amp_b, _amp_na, _amp_nb};
}
void Bmad::orbit_reference_energy_correction(
    CoordStruct &orbit,
    double p0c_new,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_orbit_reference_energy_correction(/* void* */ orbit.get_fortran_ptr(),
                                            /* double& */ p0c_new,
                                            /* Bmad::array_descriptor_t& */ _mat6_desc,
                                            /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
FixedArray1D<Real, 6> Bmad::orbit_to_floor_phase_space(CoordStruct &orbit, EleStruct &ele) {
  // floor_phase_space: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _floor_phase_space_desc;
  _floor_phase_space_desc.rank = 1;
  FixedArray1D<Real, 6> _floor_phase_space;
  _floor_phase_space_desc.data_ptr = _floor_phase_space.data();
  _floor_phase_space_desc.dims[0] = _floor_phase_space.size();
  fortran_orbit_to_floor_phase_space(/* void* */ orbit.get_fortran_ptr(),
                                     /* void* */ ele.get_fortran_ptr(),
                                     /* Bmad::array_descriptor_t& */ _floor_phase_space_desc);
  return _floor_phase_space;
}
FloorPositionStruct Bmad::orbit_to_local_curvilinear(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<int> z_direction,
    std::optional<int> relative_to
) {
  int z_direction_lvalue;
  auto *_z_direction{&z_direction_lvalue};
  if (z_direction.has_value()) {
    z_direction_lvalue = z_direction.value();
  } else {
    _z_direction = nullptr;
  }
  int relative_to_lvalue;
  auto *_relative_to{&relative_to_lvalue};
  if (relative_to.has_value()) {
    relative_to_lvalue = relative_to.value();
  } else {
    _relative_to = nullptr;
  }
  FloorPositionStruct _local_position;
  fortran_orbit_to_local_curvilinear(/* void* */ orbit.get_fortran_ptr(),
                                     /* void* */ ele.get_fortran_ptr(),
                                     /* int* */ _z_direction,
                                     /* int* */ _relative_to,
                                     /* void* */ _local_position.get_fortran_ptr());
  return std::move(_local_position);
}
Bmad::OrbitTooLarge Bmad::orbit_too_large(CoordStruct &orbit, std::optional<bool> check_momentum) {
  LatParamStruct _param;
  bool check_momentum_lvalue;
  auto *_check_momentum{&check_momentum_lvalue};
  if (check_momentum.has_value()) {
    check_momentum_lvalue = check_momentum.value();
  } else {
    _check_momentum = nullptr;
  }
  bool _is_too_large{};
  fortran_orbit_too_large(/* void* */ orbit.get_fortran_ptr(),
                          /* void* */ _param.get_fortran_ptr(),
                          /* bool* */ _check_momentum,
                          /* bool& */ _is_too_large);
  return OrbitTooLarge{std::move(_param), _is_too_large};
}
Bmad::OrderEvecsByNSimilarity Bmad::order_evecs_by_n_similarity(
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray2D<Real, 6, 6> Nmat
) {
  // evec: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _evec_desc;
  _evec_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> evec;
  std::complex<double> _evec_vec[6 * 6];
  _evec_desc.data_ptr = _evec_vec;
  _evec_desc.dims[0] = 6;
  _evec_desc.dims[1] = 6;
  // eval: inout NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _eval_desc;
  _eval_desc.rank = 1;
  _eval_desc.data_ptr = eval.data();
  _eval_desc.dims[0] = eval.size();
  // mat_tunes: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _mat_tunes_desc;
  _mat_tunes_desc.rank = 1;
  _mat_tunes_desc.data_ptr = mat_tunes.data();
  _mat_tunes_desc.dims[0] = mat_tunes.size();
  // Nmat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _Nmat_desc;
  _Nmat_desc.rank = 2;
  double _Nmat_vec[6 * 6];
  _Nmat_desc.data_ptr = _Nmat_vec;
  _Nmat_desc.dims[0] = 6;
  _Nmat_desc.dims[1] = 6;
  matrix_to_vec(Nmat, _Nmat_vec);
  bool _err_flag{};
  fortran_order_evecs_by_n_similarity(
      /* Bmad::array_descriptor_t& */ _evec_desc,
      /* Bmad::array_descriptor_t& */ _eval_desc,
      /* Bmad::array_descriptor_t& */ _mat_tunes_desc,
      /* Bmad::array_descriptor_t& */ _Nmat_desc,
      /* bool& */ _err_flag
  );
  vec_to_matrix(_evec_vec, evec);
  return OrderEvecsByNSimilarity{evec, _err_flag};
}
void Bmad::order_evecs_by_plane_dominance(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    std::optional<FixedArray1D<Real, 3>> mat_tunes
) {
  // evec: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _evec_desc;
  _evec_desc.rank = 2;
  std::complex<double> _evec_vec[6 * 6];
  _evec_desc.data_ptr = _evec_vec;
  _evec_desc.dims[0] = 6;
  _evec_desc.dims[1] = 6;
  matrix_to_vec(evec, _evec_vec);
  // eval: inout NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _eval_desc;
  _eval_desc.rank = 1;
  _eval_desc.data_ptr = eval.data();
  _eval_desc.dims[0] = eval.size();
  // mat_tunes: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _mat_tunes_desc;
  _mat_tunes_desc.rank = 1;
  if (mat_tunes.has_value()) {
    _mat_tunes_desc.data_ptr = mat_tunes->data();
    _mat_tunes_desc.dims[0] = mat_tunes->size();
  } else {
    _mat_tunes_desc.data_ptr = nullptr;
    _mat_tunes_desc.dims[0] = 0;
  }
  fortran_order_evecs_by_plane_dominance(
      /* Bmad::array_descriptor_t& */ _evec_desc,
      /* Bmad::array_descriptor_t& */ _eval_desc,
      /* Bmad::array_descriptor_t& */ _mat_tunes_desc
  );
  vec_to_matrix(_evec_vec, evec);
}
bool Bmad::order_evecs_by_tune(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray1D<Real, 3> abz_tunes
) {
  // evec: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _evec_desc;
  _evec_desc.rank = 2;
  std::complex<double> _evec_vec[6 * 6];
  _evec_desc.data_ptr = _evec_vec;
  _evec_desc.dims[0] = 6;
  _evec_desc.dims[1] = 6;
  matrix_to_vec(evec, _evec_vec);
  // eval: inout NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _eval_desc;
  _eval_desc.rank = 1;
  _eval_desc.data_ptr = eval.data();
  _eval_desc.dims[0] = eval.size();
  // mat_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _mat_tunes_desc;
  _mat_tunes_desc.rank = 1;
  _mat_tunes_desc.data_ptr = mat_tunes.data();
  _mat_tunes_desc.dims[0] = mat_tunes.size();
  // abz_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  _abz_tunes_desc.data_ptr = abz_tunes.data();
  _abz_tunes_desc.dims[0] = abz_tunes.size();
  bool _err_flag{};
  fortran_order_evecs_by_tune(
      /* Bmad::array_descriptor_t& */ _evec_desc,
      /* Bmad::array_descriptor_t& */ _eval_desc,
      /* Bmad::array_descriptor_t& */ _mat_tunes_desc,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc,
      /* bool& */ _err_flag
  );
  vec_to_matrix(_evec_vec, evec);
  return _err_flag;
}
void Bmad::order_particles_in_z(BunchStruct &bunch) {
  fortran_order_particles_in_z(/* void* */ bunch.get_fortran_ptr());
}
void Bmad::order_super_lord_slaves(LatStruct &lat, int ix_lord) {
  fortran_order_super_lord_slaves(/* void* */ lat.get_fortran_ptr(), /* int& */ ix_lord);
}
void Bmad::osc_alloc_freespace_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
) {
  // nlo: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nlo_desc;
  _nlo_desc.rank = 1;
  _nlo_desc.data_ptr = nlo.data();
  _nlo_desc.dims[0] = nlo.size();
  // nhi: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nhi_desc;
  _nhi_desc.rank = 1;
  _nhi_desc.data_ptr = nhi.data();
  _nhi_desc.dims[0] = nhi.size();
  // npad: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _npad_desc;
  _npad_desc.rank = 1;
  _npad_desc.data_ptr = npad.data();
  _npad_desc.dims[0] = npad.size();
  fortran_osc_alloc_freespace_array(
      /* Bmad::array_descriptor_t& */ _nlo_desc,
      /* Bmad::array_descriptor_t& */ _nhi_desc,
      /* Bmad::array_descriptor_t& */ _npad_desc
  );
}
void Bmad::osc_alloc_image_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
) {
  // nlo: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nlo_desc;
  _nlo_desc.rank = 1;
  _nlo_desc.data_ptr = nlo.data();
  _nlo_desc.dims[0] = nlo.size();
  // nhi: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nhi_desc;
  _nhi_desc.rank = 1;
  _nhi_desc.data_ptr = nhi.data();
  _nhi_desc.dims[0] = nhi.size();
  // npad: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _npad_desc;
  _npad_desc.rank = 1;
  _npad_desc.data_ptr = npad.data();
  _npad_desc.dims[0] = npad.size();
  fortran_osc_alloc_image_array(
      /* Bmad::array_descriptor_t& */ _nlo_desc,
      /* Bmad::array_descriptor_t& */ _nhi_desc,
      /* Bmad::array_descriptor_t& */ _npad_desc
  );
}
void Bmad::osc_alloc_rectpipe_arrays(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
) {
  // nlo: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nlo_desc;
  _nlo_desc.rank = 1;
  _nlo_desc.data_ptr = nlo.data();
  _nlo_desc.dims[0] = nlo.size();
  // nhi: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nhi_desc;
  _nhi_desc.rank = 1;
  _nhi_desc.data_ptr = nhi.data();
  _nhi_desc.dims[0] = nhi.size();
  // npad: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _npad_desc;
  _npad_desc.rank = 1;
  _npad_desc.data_ptr = npad.data();
  _npad_desc.dims[0] = npad.size();
  fortran_osc_alloc_rectpipe_arrays(
      /* Bmad::array_descriptor_t& */ _nlo_desc,
      /* Bmad::array_descriptor_t& */ _nhi_desc,
      /* Bmad::array_descriptor_t& */ _npad_desc
  );
}
void Bmad::osc_getgrnpipe(
    double gam,
    double a,
    double b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad
) {
  // delta: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _delta_desc;
  _delta_desc.rank = 1;
  _delta_desc.data_ptr = delta.data();
  _delta_desc.dims[0] = delta.size();
  // umin: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _umin_desc;
  _umin_desc.rank = 1;
  _umin_desc.data_ptr = umin.data();
  _umin_desc.dims[0] = umin.size();
  // npad: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _npad_desc;
  _npad_desc.rank = 1;
  _npad_desc.data_ptr = npad.data();
  _npad_desc.dims[0] = npad.size();
  fortran_osc_getgrnpipe(
      /* double& */ gam,
      /* double& */ a,
      /* double& */ b,
      /* Bmad::array_descriptor_t& */ _delta_desc,
      /* Bmad::array_descriptor_t& */ _umin_desc,
      /* Bmad::array_descriptor_t& */ _npad_desc
  );
}
void Bmad::osc_read_rectpipe_grn() { fortran_osc_read_rectpipe_grn(); }
void Bmad::osc_write_rectpipe_grn(
    double apipe,
    double bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double gamma
) {
  // delta: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _delta_desc;
  _delta_desc.rank = 1;
  _delta_desc.data_ptr = delta.data();
  _delta_desc.dims[0] = delta.size();
  // umin: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _umin_desc;
  _umin_desc.rank = 1;
  _umin_desc.data_ptr = umin.data();
  _umin_desc.dims[0] = umin.size();
  // umax: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _umax_desc;
  _umax_desc.rank = 1;
  _umax_desc.data_ptr = umax.data();
  _umax_desc.dims[0] = umax.size();
  // nlo: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nlo_desc;
  _nlo_desc.rank = 1;
  _nlo_desc.data_ptr = nlo.data();
  _nlo_desc.dims[0] = nlo.size();
  // nhi: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _nhi_desc;
  _nhi_desc.rank = 1;
  _nhi_desc.data_ptr = nhi.data();
  _nhi_desc.dims[0] = nhi.size();
  fortran_osc_write_rectpipe_grn(
      /* double& */ apipe,
      /* double& */ bpipe,
      /* Bmad::array_descriptor_t& */ _delta_desc,
      /* Bmad::array_descriptor_t& */ _umin_desc,
      /* Bmad::array_descriptor_t& */ _umax_desc,
      /* Bmad::array_descriptor_t& */ _nlo_desc,
      /* Bmad::array_descriptor_t& */ _nhi_desc,
      /* double& */ gamma
  );
}
void Bmad::parse_cartesian_map(
    CartesianMapStruct &ct_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
) {
  auto _delim = delim.c_str();
  fortran_parse_cartesian_map(/* void* */ ct_map.get_fortran_ptr(),
                              /* void* */ ele.get_fortran_ptr(),
                              /* void* */ lat.get_fortran_ptr(),
                              /* const char* */ _delim,
                              /* bool& */ delim_found,
                              /* bool& */ err_flag);
}
void Bmad::parse_cylindrical_map(
    CylindricalMapStruct &cl_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
) {
  auto _cl_map = &cl_map; // input, required, pointer
  auto _delim = delim.c_str();
  fortran_parse_cylindrical_map(/* void* */ &cl_map,
                                /* void* */ ele.get_fortran_ptr(),
                                /* void* */ lat.get_fortran_ptr(),
                                /* const char* */ _delim,
                                /* bool& */ delim_found,
                                /* bool& */ err_flag);
}
void Bmad::parse_gen_grad_map(
    GenGradMapStruct &gg_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
) {
  auto _gg_map = &gg_map; // input, required, pointer
  auto _delim = delim.c_str();
  fortran_parse_gen_grad_map(/* void* */ &gg_map,
                             /* void* */ ele.get_fortran_ptr(),
                             /* void* */ lat.get_fortran_ptr(),
                             /* const char* */ _delim,
                             /* bool& */ delim_found,
                             /* bool& */ err_flag);
}
void Bmad::parse_grid_field(
    GridFieldStruct &g_field,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
) {
  auto _g_field = &g_field; // input, required, pointer
  auto _delim = delim.c_str();
  fortran_parse_grid_field(/* void* */ &g_field,
                           /* void* */ ele.get_fortran_ptr(),
                           /* void* */ lat.get_fortran_ptr(),
                           /* const char* */ _delim,
                           /* bool& */ delim_found,
                           /* bool& */ err_flag);
}
void Bmad::parse_integer_list(
    std::string err_str,
    LatStruct &lat,
    FArray1D<Int> &int_array,
    bool exact_size,
    std::string delim,
    bool delim_found,
    bool is_ok,
    std::optional<std::string> open_delim,
    std::optional<std::string> separator,
    std::optional<std::string> close_delim,
    std::optional<int> default_value
) {
  auto _err_str = err_str.c_str();
  // int_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _int_array_desc;
  _int_array_desc.rank = 1;
  _int_array_desc.data_ptr = int_array.data();
  _int_array_desc.dims[0] = int_array.size();
  auto _delim = delim.c_str();
  const char *_open_delim = open_delim.has_value() ? open_delim->c_str() : nullptr;
  const char *_separator = separator.has_value() ? separator->c_str() : nullptr;
  const char *_close_delim = close_delim.has_value() ? close_delim->c_str() : nullptr;
  int default_value_lvalue;
  auto *_default_value{&default_value_lvalue};
  if (default_value.has_value()) {
    default_value_lvalue = default_value.value();
  } else {
    _default_value = nullptr;
  }
  fortran_parse_integer_list(
      /* const char* */ _err_str,
      /* void* */ lat.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _int_array_desc,
      /* bool& */ exact_size,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* const char* */ _open_delim,
      /* const char* */ _separator,
      /* const char* */ _close_delim,
      /* int* */ _default_value,
      /* bool& */ is_ok
  );
}
Bmad::ParseIntegerList2 Bmad::parse_integer_list2(
    std::string err_str,
    LatStruct &lat,
    IntAlloc1D &int_array,
    std::optional<int> num_expected,
    std::optional<std::string> open_delim,
    std::optional<std::string> separator,
    std::optional<std::string> close_delim,
    std::optional<int> default_value
) {
  auto _err_str = err_str.c_str();
  // intent=inout allocatable general array
  int _num_found{};
  char _delim[4096];
  bool _delim_found{};
  int num_expected_lvalue;
  auto *_num_expected{&num_expected_lvalue};
  if (num_expected.has_value()) {
    num_expected_lvalue = num_expected.value();
  } else {
    _num_expected = nullptr;
  }
  const char *_open_delim = open_delim.has_value() ? open_delim->c_str() : nullptr;
  const char *_separator = separator.has_value() ? separator->c_str() : nullptr;
  const char *_close_delim = close_delim.has_value() ? close_delim->c_str() : nullptr;
  int default_value_lvalue;
  auto *_default_value{&default_value_lvalue};
  if (default_value.has_value()) {
    default_value_lvalue = default_value.value();
  } else {
    _default_value = nullptr;
  }
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
      /* bool& */ _is_ok
  );
  return ParseIntegerList2{_num_found, _delim, _delim_found, _is_ok};
}
Bmad::ParseRealList Bmad::parse_real_list(
    LatStruct &lat,
    std::string err_str,
    FArray1D<Real> &real_array,
    bool exact_size,
    bool is_ok,
    std::optional<std::string> open_delim,
    std::optional<std::string> separator,
    std::optional<std::string> close_delim,
    std::optional<double> default_value
) {
  auto _err_str = err_str.c_str();
  // real_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _real_array_desc;
  _real_array_desc.rank = 1;
  _real_array_desc.data_ptr = real_array.data();
  _real_array_desc.dims[0] = real_array.size();
  char _delim[4096];
  bool _delim_found{};
  const char *_open_delim = open_delim.has_value() ? open_delim->c_str() : nullptr;
  const char *_separator = separator.has_value() ? separator->c_str() : nullptr;
  const char *_close_delim = close_delim.has_value() ? close_delim->c_str() : nullptr;
  double default_value_lvalue;
  auto *_default_value{&default_value_lvalue};
  if (default_value.has_value()) {
    default_value_lvalue = default_value.value();
  } else {
    _default_value = nullptr;
  }
  int _num_found{};
  fortran_parse_real_list(/* void* */ lat.get_fortran_ptr(),
                          /* const char* */ _err_str,
                          /* Bmad::array_descriptor_t& */ _real_array_desc,
                          /* bool& */ exact_size,
                          /* const char* */ _delim,
                          /* bool& */ _delim_found,
                          /* const char* */ _open_delim,
                          /* const char* */ _separator,
                          /* const char* */ _close_delim,
                          /* double* */ _default_value,
                          /* int& */ _num_found,
                          /* bool& */ is_ok);
  return ParseRealList{_delim, _delim_found, _num_found};
}
Bmad::ParseRealList2 Bmad::parse_real_list2(
    LatStruct &lat,
    std::string err_str,
    RealAlloc1D &real_array,
    std::optional<int> num_expected,
    std::optional<std::string> open_brace,
    std::optional<std::string> separator,
    std::optional<std::string> close_brace,
    std::optional<double> default_value,
    std::optional<bool> single_value
) {
  auto _err_str = err_str.c_str();
  // intent=inout allocatable general array
  int _num_found{};
  char _delim[4096];
  bool _delim_found{};
  int num_expected_lvalue;
  auto *_num_expected{&num_expected_lvalue};
  if (num_expected.has_value()) {
    num_expected_lvalue = num_expected.value();
  } else {
    _num_expected = nullptr;
  }
  const char *_open_brace = open_brace.has_value() ? open_brace->c_str() : nullptr;
  const char *_separator = separator.has_value() ? separator->c_str() : nullptr;
  const char *_close_brace = close_brace.has_value() ? close_brace->c_str() : nullptr;
  double default_value_lvalue;
  auto *_default_value{&default_value_lvalue};
  if (default_value.has_value()) {
    default_value_lvalue = default_value.value();
  } else {
    _default_value = nullptr;
  }
  bool single_value_lvalue;
  auto *_single_value{&single_value_lvalue};
  if (single_value.has_value()) {
    single_value_lvalue = single_value.value();
  } else {
    _single_value = nullptr;
  }
  bool _is_ok{};
  fortran_parse_real_list2(/* void* */ lat.get_fortran_ptr(),
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
void Bmad::parser_add_constant(std::string word, LatStruct &lat, bool redef_is_error) {
  auto _word = word.c_str();
  fortran_parser_add_constant(
      /* const char* */ _word,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ redef_is_error
  );
}
void Bmad::parser_call_check(
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool call_found,
    std::optional<bool> err_flag
) {
  auto _word = word.c_str();
  auto _delim = delim.c_str();
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_parser_call_check(
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ call_found,
      /* bool* */ _err_flag
  );
}
Bmad::ParserFastComplexRead
Bmad::parser_fast_complex_read(FArray1D<Complex> &cmplx_vec, EleStruct &ele, std::string err_str) {
  // cmplx_vec: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _cmplx_vec_desc;
  _cmplx_vec_desc.rank = 1;
  _cmplx_vec_desc.data_ptr = cmplx_vec.data();
  _cmplx_vec_desc.dims[0] = cmplx_vec.size();
  char _delim[4096];
  auto _err_str = err_str.c_str();
  bool _is_ok{};
  fortran_parser_fast_complex_read(
      /* Bmad::array_descriptor_t& */ _cmplx_vec_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim,
      /* const char* */ _err_str,
      /* bool& */ _is_ok
  );
  return ParserFastComplexRead{_delim, _is_ok};
}
void Bmad::parser_fast_integer_read(
    FArray1D<Int> &int_vec,
    EleStruct &ele,
    std::string delim_wanted,
    std::string err_str,
    bool is_ok
) {
  // int_vec: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _int_vec_desc;
  _int_vec_desc.rank = 1;
  _int_vec_desc.data_ptr = int_vec.data();
  _int_vec_desc.dims[0] = int_vec.size();
  auto _delim_wanted = delim_wanted.c_str();
  auto _err_str = err_str.c_str();
  fortran_parser_fast_integer_read(
      /* Bmad::array_descriptor_t& */ _int_vec_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _delim_wanted,
      /* const char* */ _err_str,
      /* bool& */ is_ok
  );
}
Bmad::ParserFastRealRead Bmad::parser_fast_real_read(
    FArray1D<Real> &real_vec,
    EleStruct &ele,
    std::string end_delims,
    std::string err_str,
    std::optional<bool> exact_size
) {
  // real_vec: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _real_vec_desc;
  _real_vec_desc.rank = 1;
  _real_vec_desc.data_ptr = real_vec.data();
  _real_vec_desc.dims[0] = real_vec.size();
  auto _end_delims = end_delims.c_str();
  char _delim[4096];
  auto _err_str = err_str.c_str();
  bool exact_size_lvalue;
  auto *_exact_size{&exact_size_lvalue};
  if (exact_size.has_value()) {
    exact_size_lvalue = exact_size.value();
  } else {
    _exact_size = nullptr;
  }
  int _n_real{};
  bool _is_ok{};
  fortran_parser_fast_real_read(
      /* Bmad::array_descriptor_t& */ _real_vec_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _end_delims,
      /* const char* */ _delim,
      /* const char* */ _err_str,
      /* bool* */ _exact_size,
      /* int& */ _n_real,
      /* bool& */ _is_ok
  );
  return ParserFastRealRead{_delim, _n_real, _is_ok};
}
void Bmad::parser_file_stack(
    std::string how,
    std::optional<std::string> file_name_in,
    std::optional<bool> finished,
    std::optional<bool> err,
    std::optional<bool> open_file,
    std::optional<bool> abort_on_open_error
) {
  auto _how = how.c_str();
  const char *_file_name_in = file_name_in.has_value() ? file_name_in->c_str() : nullptr;
  bool finished_lvalue;
  auto *_finished{&finished_lvalue};
  if (finished.has_value()) {
    finished_lvalue = finished.value();
  } else {
    _finished = nullptr;
  }
  bool err_lvalue;
  auto *_err{&err_lvalue};
  if (err.has_value()) {
    err_lvalue = err.value();
  } else {
    _err = nullptr;
  }
  bool open_file_lvalue;
  auto *_open_file{&open_file_lvalue};
  if (open_file.has_value()) {
    open_file_lvalue = open_file.value();
  } else {
    _open_file = nullptr;
  }
  bool abort_on_open_error_lvalue;
  auto *_abort_on_open_error{&abort_on_open_error_lvalue};
  if (abort_on_open_error.has_value()) {
    abort_on_open_error_lvalue = abort_on_open_error.value();
  } else {
    _abort_on_open_error = nullptr;
  }
  fortran_parser_file_stack(
      /* const char* */ _how,
      /* const char* */ _file_name_in,
      /* bool* */ _finished,
      /* bool* */ _err,
      /* bool* */ _open_file,
      /* bool* */ _abort_on_open_error
  );
}
void Bmad::parser_get_integer(
    int int_val,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool err,
    std::optional<std::string> str1,
    std::optional<std::string> str2
) {
  auto _word = word.c_str();
  auto _delim = delim.c_str();
  const char *_str1 = str1.has_value() ? str1->c_str() : nullptr;
  const char *_str2 = str2.has_value() ? str2->c_str() : nullptr;
  fortran_parser_get_integer(
      /* int& */ int_val,
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err,
      /* const char* */ _str1,
      /* const char* */ _str2
  );
}
void Bmad::parser_get_logical(
    std::string attrib_name,
    bool this_logic,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool err
) {
  auto _attrib_name = attrib_name.c_str();
  auto _ele_name = ele_name.c_str();
  auto _delim = delim.c_str();
  fortran_parser_get_logical(
      /* const char* */ _attrib_name,
      /* bool& */ this_logic,
      /* const char* */ _ele_name,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* bool& */ err
  );
}
void Bmad::parser_identify_fork_to_element(LatStruct &lat) {
  fortran_parser_identify_fork_to_element(/* void* */ lat.get_fortran_ptr());
}
void Bmad::parser_init_custom_elements(LatStruct &lat) {
  fortran_parser_init_custom_elements(/* void* */ lat.get_fortran_ptr());
}
void Bmad::parser_print_line(LatStruct &lat, bool end_of_file) {
  fortran_parser_print_line(/* void* */ lat.get_fortran_ptr(), /* bool& */ end_of_file);
}
void Bmad::parser_read_lr_wake(EleStruct &ele, std::string delim, bool delim_found, bool err_flag) {
  auto _delim = delim.c_str();
  fortran_parser_read_lr_wake(/* void* */ ele.get_fortran_ptr(),
                              /* const char* */ _delim,
                              /* bool& */ delim_found,
                              /* bool& */ err_flag);
}
void Bmad::parser_read_old_format_lr_wake(EleStruct &ele, std::string lr_file_name) {
  auto _lr_file_name = lr_file_name.c_str();
  fortran_parser_read_old_format_lr_wake(/* void* */ ele.get_fortran_ptr(),
                                         /* const char* */ _lr_file_name);
}
void Bmad::parser_read_old_format_sr_wake(EleStruct &ele, std::string sr_file_name) {
  auto _sr_file_name = sr_file_name.c_str();
  fortran_parser_read_old_format_sr_wake(/* void* */ ele.get_fortran_ptr(),
                                         /* const char* */ _sr_file_name);
}
void Bmad::parser_read_sr_wake(EleStruct &ele, std::string delim, bool delim_found, bool err_flag) {
  auto _delim = delim.c_str();
  fortran_parser_read_sr_wake(/* void* */ ele.get_fortran_ptr(),
                              /* const char* */ _delim,
                              /* bool& */ delim_found,
                              /* bool& */ err_flag);
}
ControlStruct
Bmad::parser_transfer_control_struct(ControlStruct &con_in, EleStruct &lord, int ix_var) {
  ControlStruct _con_out;
  fortran_parser_transfer_control_struct(/* void* */ con_in.get_fortran_ptr(),
                                         /* void* */ _con_out.get_fortran_ptr(),
                                         /* void* */ lord.get_fortran_ptr(),
                                         /* int& */ ix_var);
  return std::move(_con_out);
}
CoordStruct Bmad::particle_in_global_frame(
    CoordStruct &orb,
    BranchStruct &branch,
    std::optional<bool> in_time_coordinates,
    std::optional<bool> in_body_frame,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat_out
) {
  bool in_time_coordinates_lvalue;
  auto *_in_time_coordinates{&in_time_coordinates_lvalue};
  if (in_time_coordinates.has_value()) {
    in_time_coordinates_lvalue = in_time_coordinates.value();
  } else {
    _in_time_coordinates = nullptr;
  }
  bool in_body_frame_lvalue;
  auto *_in_body_frame{&in_body_frame_lvalue};
  if (in_body_frame.has_value()) {
    in_body_frame_lvalue = in_body_frame.value();
  } else {
    _in_body_frame = nullptr;
  }
  // w_mat_out: inout NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_out_desc;
  _w_mat_out_desc.rank = 2;
  double _w_mat_out_vec[3 * 3];
  _w_mat_out_desc.data_ptr = _w_mat_out_vec;
  if (w_mat_out.has_value()) {
    matrix_to_vec(w_mat_out.value(), _w_mat_out_vec);
    _w_mat_out_desc.dims[0] = 3;
    _w_mat_out_desc.dims[1] = 3;
  } else {
    _w_mat_out_desc.data_ptr = nullptr;
  }
  CoordStruct _particle;
  fortran_particle_in_global_frame(/* void* */ orb.get_fortran_ptr(),
                                   /* void* */ branch.get_fortran_ptr(),
                                   /* bool* */ _in_time_coordinates,
                                   /* bool* */ _in_body_frame,
                                   /* Bmad::array_descriptor_t& */ _w_mat_out_desc,
                                   /* void* */ _particle.get_fortran_ptr());
  if (w_mat_out.has_value())
    vec_to_matrix(_w_mat_out_vec, w_mat_out.value());
  return std::move(_particle);
}
bool Bmad::particle_is_moving_backwards(CoordStruct &orbit) {
  bool _is_moving_backwards{};
  fortran_particle_is_moving_backwards(/* void* */ orbit.get_fortran_ptr(),
                                       /* bool& */ _is_moving_backwards);
  return _is_moving_backwards;
}
bool Bmad::particle_is_moving_forward(CoordStruct &orbit, std::optional<int> dir) {
  int dir_lvalue;
  auto *_dir{&dir_lvalue};
  if (dir.has_value()) {
    dir_lvalue = dir.value();
  } else {
    _dir = nullptr;
  }
  bool _is_moving_forward{};
  fortran_particle_is_moving_forward(/* void* */ orbit.get_fortran_ptr(),
                                     /* int* */ _dir,
                                     /* bool& */ _is_moving_forward);
  return _is_moving_forward;
}
long double Bmad::particle_rf_time(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<bool> reference_active_edge,
    std::optional<double> s_rel,
    std::optional<bool> time_coords,
    std::optional<double> rf_freq,
    std::optional<bool> abs_time
) {
  bool reference_active_edge_lvalue;
  auto *_reference_active_edge{&reference_active_edge_lvalue};
  if (reference_active_edge.has_value()) {
    reference_active_edge_lvalue = reference_active_edge.value();
  } else {
    _reference_active_edge = nullptr;
  }
  double s_rel_lvalue;
  auto *_s_rel{&s_rel_lvalue};
  if (s_rel.has_value()) {
    s_rel_lvalue = s_rel.value();
  } else {
    _s_rel = nullptr;
  }
  bool time_coords_lvalue;
  auto *_time_coords{&time_coords_lvalue};
  if (time_coords.has_value()) {
    time_coords_lvalue = time_coords.value();
  } else {
    _time_coords = nullptr;
  }
  double rf_freq_lvalue;
  auto *_rf_freq{&rf_freq_lvalue};
  if (rf_freq.has_value()) {
    rf_freq_lvalue = rf_freq.value();
  } else {
    _rf_freq = nullptr;
  }
  bool abs_time_lvalue;
  auto *_abs_time{&abs_time_lvalue};
  if (abs_time.has_value()) {
    abs_time_lvalue = abs_time.value();
  } else {
    _abs_time = nullptr;
  }
  long double _time{};
  fortran_particle_rf_time(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* bool* */ _reference_active_edge,
                           /* double* */ _s_rel,
                           /* bool* */ _time_coords,
                           /* double* */ _rf_freq,
                           /* bool* */ _abs_time,
                           /* long double& */ _time);
  return _time;
}
bool Bmad::patch_flips_propagation_direction(double x_pitch, double y_pitch) {
  bool _is_flip{};
  fortran_patch_flips_propagation_direction(
      /* double& */ x_pitch,
      /* double& */ y_pitch,
      /* bool& */ _is_flip
  );
  return _is_flip;
}
double Bmad::patch_length(EleStruct &patch, std::optional<int> ref_coords) {
  int ref_coords_lvalue;
  auto *_ref_coords{&ref_coords_lvalue};
  if (ref_coords.has_value()) {
    ref_coords_lvalue = ref_coords.value();
  } else {
    _ref_coords = nullptr;
  }
  double _length{};
  fortran_patch_length(/* void* */ patch.get_fortran_ptr(),
                       /* int* */ _ref_coords,
                       /* double& */ _length);
  return _length;
}
Bmad::PhotonAbsorptionAndPhaseShift
Bmad::photon_absorption_and_phase_shift(std::string material, double Energy) {
  auto _material = material.c_str();
  double _absorption{};
  double _phase_shift{};
  bool _err_flag{};
  fortran_photon_absorption_and_phase_shift(
      /* const char* */ _material,
      /* double& */ Energy,
      /* double& */ _absorption,
      /* double& */ _phase_shift,
      /* bool& */ _err_flag
  );
  return PhotonAbsorptionAndPhaseShift{_absorption, _phase_shift, _err_flag};
}
Bmad::PhotonAddToDetectorStatistics Bmad::photon_add_to_detector_statistics(
    CoordStruct &orbit0,
    CoordStruct &orbit,
    EleStruct &ele,
    optional_ref<PixelPtStruct> pixel_pt
) {
  int _ix_pt{};
  int _iy_pt{};
  auto *_pixel_pt =
      pixel_pt.has_value() ? pixel_pt->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_photon_add_to_detector_statistics(/* void* */ orbit0.get_fortran_ptr(),
                                            /* void* */ orbit.get_fortran_ptr(),
                                            /* void* */ ele.get_fortran_ptr(),
                                            /* int& */ _ix_pt,
                                            /* int& */ _iy_pt,
                                            /* void* */ _pixel_pt);
  return PhotonAddToDetectorStatistics{_ix_pt, _iy_pt};
}
Bmad::PhotonReflection
Bmad::photon_reflection(double graze_angle_in, double energy, PhotonReflectSurfaceStruct &surface) {
  double _graze_angle_out{};
  double _phi_out{};
  fortran_photon_reflection(
      /* double& */ graze_angle_in,
      /* double& */ energy,
      /* void* */ surface.get_fortran_ptr(),
      /* double& */ _graze_angle_out,
      /* double& */ _phi_out
  );
  return PhotonReflection{_graze_angle_out, _phi_out};
}
PhotonReflectSurfaceStruct Bmad::photon_reflection_std_surface_init() {
  PhotonReflectSurfaceStruct _surface;
  fortran_photon_reflection_std_surface_init(/* void* */ _surface.get_fortran_ptr());
  return std::move(_surface);
}
Bmad::PhotonReflectivity
Bmad::photon_reflectivity(double angle, double energy, PhotonReflectSurfaceStruct &surface) {
  double _p_reflect{};
  double _rel_p_specular{};
  fortran_photon_reflectivity(
      /* double& */ angle,
      /* double& */ energy,
      /* void* */ surface.get_fortran_ptr(),
      /* double& */ _p_reflect,
      /* double& */ _rel_p_specular
  );
  return PhotonReflectivity{_p_reflect, _rel_p_specular};
}
TargetPointStruct Bmad::photon_target_corner_calc(
    EleStruct &aperture_ele,
    double x_lim,
    double y_lim,
    double z_lim,
    EleStruct &source_ele
) {
  TargetPointStruct _corner;
  fortran_photon_target_corner_calc(/* void* */ aperture_ele.get_fortran_ptr(),
                                    /* double& */ x_lim,
                                    /* double& */ y_lim,
                                    /* double& */ z_lim,
                                    /* void* */ source_ele.get_fortran_ptr(),
                                    /* void* */ _corner.get_fortran_ptr());
  return std::move(_corner);
}
void Bmad::photon_target_setup(EleStruct &ele) {
  fortran_photon_target_setup(/* void* */ ele.get_fortran_ptr());
}
int Bmad::photon_type(EleStruct &ele) {
  int _e_type{};
  fortran_photon_type(/* void* */ ele.get_fortran_ptr(), /* int& */ _e_type);
  return _e_type;
}
int Bmad::physical_ele_end(
    int track_end,
    CoordStruct &orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end
) {
  bool return_stream_end_lvalue;
  auto *_return_stream_end{&return_stream_end_lvalue};
  if (return_stream_end.has_value()) {
    return_stream_end_lvalue = return_stream_end.value();
  } else {
    _return_stream_end = nullptr;
  }
  int _physical_end{};
  fortran_physical_ele_end(
      /* int& */ track_end,
      /* void* */ orbit.get_fortran_ptr(),
      /* int& */ ele_orientation,
      /* bool* */ _return_stream_end,
      /* int& */ _physical_end
  );
  return _physical_end;
}
void Bmad::point_photon_emission(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    int direction,
    double max_target_area,
    std::optional<FixedArray2D<Real, 3, 3>> w_to_surface
) {
  // w_to_surface: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_to_surface_desc;
  _w_to_surface_desc.rank = 2;
  double _w_to_surface_vec[3 * 3];
  _w_to_surface_desc.data_ptr = _w_to_surface_vec;
  if (w_to_surface.has_value()) {
    matrix_to_vec(w_to_surface.value(), _w_to_surface_vec);
    _w_to_surface_desc.dims[0] = 3;
    _w_to_surface_desc.dims[1] = 3;
  } else {
    _w_to_surface_desc.data_ptr = nullptr;
  }
  fortran_point_photon_emission(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ param.get_fortran_ptr(),
                                /* void* */ orbit.get_fortran_ptr(),
                                /* int& */ direction,
                                /* double& */ max_target_area,
                                /* Bmad::array_descriptor_t& */ _w_to_surface_desc);
}
std::optional<BranchStruct> Bmad::pointer_to_branch(EleStruct &ele) {
  void *_branch_ptr;
  fortran_pointer_to_branch_given_ele(/* void* */ ele.get_fortran_ptr(), /* void* */ &_branch_ptr);
  return std::move((_branch_ptr ? std::make_optional<BranchStruct>(_branch_ptr) : std::nullopt));
}
std::optional<BranchStruct> Bmad::pointer_to_branch(
    std::string branch_name,
    LatStruct &lat,
    std::optional<bool> parameter_is_branch0,
    std::optional<int> blank_branch
) {
  auto _branch_name = branch_name.c_str();
  bool parameter_is_branch0_lvalue;
  auto *_parameter_is_branch0{&parameter_is_branch0_lvalue};
  if (parameter_is_branch0.has_value()) {
    parameter_is_branch0_lvalue = parameter_is_branch0.value();
  } else {
    _parameter_is_branch0 = nullptr;
  }
  int blank_branch_lvalue;
  auto *_blank_branch{&blank_branch_lvalue};
  if (blank_branch.has_value()) {
    blank_branch_lvalue = blank_branch.value();
  } else {
    _blank_branch = nullptr;
  }
  void *_branch_ptr;
  fortran_pointer_to_branch_given_name(
      /* const char* */ _branch_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _parameter_is_branch0,
      /* int* */ _blank_branch,
      /* void* */ &_branch_ptr
  );
  return std::move((_branch_ptr ? std::make_optional<BranchStruct>(_branch_ptr) : std::nullopt));
}
std::optional<EleStruct>
Bmad::pointer_to_ele(LatStruct &lat, int ix_ele, std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  void *_ele_ptr;
  fortran_pointer_to_ele1(/* void* */ lat.get_fortran_ptr(),
                          /* int& */ ix_ele,
                          /* int* */ _ix_branch,
                          /* void* */ &_ele_ptr);
  return std::move((_ele_ptr ? std::make_optional<EleStruct>(_ele_ptr) : std::nullopt));
}
std::optional<EleStruct> Bmad::pointer_to_ele(LatStruct &lat, LatEleLocStruct &ele_loc) {
  void *_ele_ptr;
  fortran_pointer_to_ele2(/* void* */ lat.get_fortran_ptr(),
                          /* void* */ ele_loc.get_fortran_ptr(),
                          /* void* */ &_ele_ptr);
  return std::move((_ele_ptr ? std::make_optional<EleStruct>(_ele_ptr) : std::nullopt));
}
std::optional<EleStruct> Bmad::pointer_to_ele(LatStruct &lat, std::string ele_name) {
  auto _ele_name = ele_name.c_str();
  void *_ele_ptr;
  fortran_pointer_to_ele3(/* void* */ lat.get_fortran_ptr(),
                          /* const char* */ _ele_name,
                          /* void* */ &_ele_ptr);
  return std::move((_ele_ptr ? std::make_optional<EleStruct>(_ele_ptr) : std::nullopt));
}
std::optional<EleStruct> Bmad::pointer_to_ele(LatStruct &lat, EleStruct &foreign_ele) {
  void *_ele_ptr;
  fortran_pointer_to_ele4(/* void* */ lat.get_fortran_ptr(),
                          /* void* */ foreign_ele.get_fortran_ptr(),
                          /* void* */ &_ele_ptr);
  return std::move((_ele_ptr ? std::make_optional<EleStruct>(_ele_ptr) : std::nullopt));
}
Bmad::PointerToElementAtS Bmad::pointer_to_element_at_s(
    BranchStruct &branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err
) {
  bool _err_flag{};
  double _s_eff{};
  CoordStruct _position;
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  void *_ele;
  fortran_pointer_to_element_at_s(/* void* */ branch.get_fortran_ptr(),
                                  /* double& */ s,
                                  /* bool& */ choose_max,
                                  /* bool& */ _err_flag,
                                  /* double& */ _s_eff,
                                  /* void* */ _position.get_fortran_ptr(),
                                  /* bool* */ _print_err,
                                  /* void* */ &_ele);
  return PointerToElementAtS{
      _err_flag,
      _s_eff,
      std::move(_position),
      std::move((_ele ? std::make_optional<EleStruct>(_ele) : std::nullopt))
  };
}
std::optional<Fibre> Bmad::pointer_to_fibre(EleStruct &ele) {
  void *_assoc_fibre;
  fortran_pointer_to_fibre(/* void* */ ele.get_fortran_ptr(), /* void* */ &_assoc_fibre);
  return std::move((_assoc_fibre ? std::make_optional<Fibre>(_assoc_fibre) : std::nullopt));
}
Bmad::PointerToFieldEle Bmad::pointer_to_field_ele(EleStruct &ele, int ix_field_ele) {
  double _dz_offset{};
  void *_field_ele;
  fortran_pointer_to_field_ele(/* void* */ ele.get_fortran_ptr(),
                               /* int& */ ix_field_ele,
                               /* double& */ _dz_offset,
                               /* void* */ &_field_ele);
  return PointerToFieldEle{
      _dz_offset,
      std::move((_field_ele ? std::make_optional<EleStruct>(_field_ele) : std::nullopt))
  };
}
Bmad::PointerToGirder Bmad::pointer_to_girder(EleStruct &ele) {
  int _ix_slave_back{};
  void *_girder;
  fortran_pointer_to_girder(/* void* */ ele.get_fortran_ptr(),
                            /* int& */ _ix_slave_back,
                            /* void* */ &_girder);
  return PointerToGirder{
      _ix_slave_back,
      std::move((_girder ? std::make_optional<EleStruct>(_girder) : std::nullopt))
  };
}
Bmad::PointerToLord
Bmad::pointer_to_lord(EleStruct &slave, int ix_lord, std::optional<int> lord_type) {
  void *_control;
  int _ix_slave_back{};
  int lord_type_lvalue;
  auto *_lord_type{&lord_type_lvalue};
  if (lord_type.has_value()) {
    lord_type_lvalue = lord_type.value();
  } else {
    _lord_type = nullptr;
  }
  int _ix_control{};
  int _ix_ic{};
  void *_lord_ptr;
  fortran_pointer_to_lord(/* void* */ slave.get_fortran_ptr(),
                          /* int& */ ix_lord,
                          /* void* */ &_control,
                          /* int& */ _ix_slave_back,
                          /* int* */ _lord_type,
                          /* int& */ _ix_control,
                          /* int& */ _ix_ic,
                          /* void* */ &_lord_ptr);
  return PointerToLord{
      std::move((_control ? std::make_optional<ControlStruct>(_control) : std::nullopt)),
      _ix_slave_back,
      _ix_control,
      _ix_ic,
      std::move((_lord_ptr ? std::make_optional<EleStruct>(_lord_ptr) : std::nullopt))
  };
}
Bmad::PointerToMultipassLord Bmad::pointer_to_multipass_lord(EleStruct &ele) {
  int _ix_pass{};
  void *_super_lord;
  void *_multi_lord;
  fortran_pointer_to_multipass_lord(/* void* */ ele.get_fortran_ptr(),
                                    /* int& */ _ix_pass,
                                    /* void* */ &_super_lord,
                                    /* void* */ &_multi_lord);
  return PointerToMultipassLord{
      _ix_pass,
      std::move((_super_lord ? std::make_optional<EleStruct>(_super_lord) : std::nullopt)),
      std::move((_multi_lord ? std::make_optional<EleStruct>(_multi_lord) : std::nullopt))
  };
}
void Bmad::pointer_to_next_ele(
    EleStruct &this_ele,
    EleStruct &next_ele,
    std::optional<int> offset,
    std::optional<bool> skip_beginning,
    std::optional<bool> follow_fork
) {
  int offset_lvalue;
  auto *_offset{&offset_lvalue};
  if (offset.has_value()) {
    offset_lvalue = offset.value();
  } else {
    _offset = nullptr;
  }
  bool skip_beginning_lvalue;
  auto *_skip_beginning{&skip_beginning_lvalue};
  if (skip_beginning.has_value()) {
    skip_beginning_lvalue = skip_beginning.value();
  } else {
    _skip_beginning = nullptr;
  }
  bool follow_fork_lvalue;
  auto *_follow_fork{&follow_fork_lvalue};
  if (follow_fork.has_value()) {
    follow_fork_lvalue = follow_fork.value();
  } else {
    _follow_fork = nullptr;
  }
  auto _next_ele = &next_ele; // input, required, pointer
  fortran_pointer_to_next_ele(/* void* */ this_ele.get_fortran_ptr(),
                              /* int* */ _offset,
                              /* bool* */ _skip_beginning,
                              /* bool* */ _follow_fork,
                              /* void* */ &next_ele);
}
Bmad::PointerToSlave
Bmad::pointer_to_slave(EleStruct &lord, int ix_slave, std::optional<int> slave_type) {
  void *_control;
  int slave_type_lvalue;
  auto *_slave_type{&slave_type_lvalue};
  if (slave_type.has_value()) {
    slave_type_lvalue = slave_type.value();
  } else {
    _slave_type = nullptr;
  }
  int _ix_lord_back{};
  int _ix_control{};
  int _ix_ic{};
  void *_slave_ptr;
  fortran_pointer_to_slave(/* void* */ lord.get_fortran_ptr(),
                           /* int& */ ix_slave,
                           /* void* */ &_control,
                           /* int* */ _slave_type,
                           /* int& */ _ix_lord_back,
                           /* int& */ _ix_control,
                           /* int& */ _ix_ic,
                           /* void* */ &_slave_ptr);
  return PointerToSlave{
      std::move((_control ? std::make_optional<ControlStruct>(_control) : std::nullopt)),
      _ix_lord_back,
      _ix_control,
      _ix_ic,
      std::move((_slave_ptr ? std::make_optional<EleStruct>(_slave_ptr) : std::nullopt))
  };
}
Bmad::PointerToSuperLord
Bmad::pointer_to_super_lord(EleStruct &slave, std::optional<int> lord_type) {
  void *_control;
  int _ix_slave_back{};
  int _ix_control{};
  int _ix_ic{};
  int lord_type_lvalue;
  auto *_lord_type{&lord_type_lvalue};
  if (lord_type.has_value()) {
    lord_type_lvalue = lord_type.value();
  } else {
    _lord_type = nullptr;
  }
  void *_lord_ptr;
  fortran_pointer_to_super_lord(/* void* */ slave.get_fortran_ptr(),
                                /* void* */ &_control,
                                /* int& */ _ix_slave_back,
                                /* int& */ _ix_control,
                                /* int& */ _ix_ic,
                                /* int* */ _lord_type,
                                /* void* */ &_lord_ptr);
  return PointerToSuperLord{
      std::move((_control ? std::make_optional<ControlStruct>(_control) : std::nullopt)),
      _ix_slave_back,
      _ix_control,
      _ix_ic,
      std::move((_lord_ptr ? std::make_optional<EleStruct>(_lord_ptr) : std::nullopt))
  };
}
Bmad::PointerToSurfaceDisplacementPt Bmad::pointer_to_surface_displacement_pt(
    EleStruct &ele,
    bool nearest,
    double x,
    double y,
    std::optional<bool> extend_grid
) {
  int _ix{};
  int _iy{};
  bool extend_grid_lvalue;
  auto *_extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  double _xx{};
  double _yy{};
  void *_pt;
  fortran_pointer_to_surface_displacement_pt(/* void* */ ele.get_fortran_ptr(),
                                             /* bool& */ nearest,
                                             /* double& */ x,
                                             /* double& */ y,
                                             /* int& */ _ix,
                                             /* int& */ _iy,
                                             /* bool* */ _extend_grid,
                                             /* double& */ _xx,
                                             /* double& */ _yy,
                                             /* void* */ &_pt);
  return PointerToSurfaceDisplacementPt{
      _ix,
      _iy,
      _xx,
      _yy,
      std::move((_pt ? std::make_optional<SurfaceDisplacementPtStruct>(_pt) : std::nullopt))
  };
}
Bmad::PointerToSurfaceSegmentedPt Bmad::pointer_to_surface_segmented_pt(
    EleStruct &ele,
    bool nearest,
    double x,
    double y,
    std::optional<bool> extend_grid
) {
  int _ix{};
  int _iy{};
  bool extend_grid_lvalue;
  auto *_extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  double _xx{};
  double _yy{};
  void *_pt;
  fortran_pointer_to_surface_segmented_pt(/* void* */ ele.get_fortran_ptr(),
                                          /* bool& */ nearest,
                                          /* double& */ x,
                                          /* double& */ y,
                                          /* int& */ _ix,
                                          /* int& */ _iy,
                                          /* bool* */ _extend_grid,
                                          /* double& */ _xx,
                                          /* double& */ _yy,
                                          /* void* */ &_pt);
  return PointerToSurfaceSegmentedPt{
      _ix,
      _iy,
      _xx,
      _yy,
      std::move((_pt ? std::make_optional<SurfaceSegmentedPtStruct>(_pt) : std::nullopt))
  };
}
Bmad::PointerToWakeEle Bmad::pointer_to_wake_ele(EleStruct &ele) {
  double _delta_s{};
  void *_wake_ele;
  fortran_pointer_to_wake_ele(/* void* */ ele.get_fortran_ptr(),
                              /* double& */ _delta_s,
                              /* void* */ &_wake_ele);
  return PointerToWakeEle{
      _delta_s,
      std::move((_wake_ele ? std::make_optional<EleStruct>(_wake_ele) : std::nullopt))
  };
}
Bmad::PointerToWall3d Bmad::pointer_to_wall3d(EleStruct &ele, std::optional<int> ix_wall) {
  int ix_wall_lvalue;
  auto *_ix_wall{&ix_wall_lvalue};
  if (ix_wall.has_value()) {
    ix_wall_lvalue = ix_wall.value();
  } else {
    _ix_wall = nullptr;
  }
  double _ds_offset{};
  bool _is_branch_wall{};
  void *_wall3d;
  fortran_pointer_to_wall3d(/* void* */ ele.get_fortran_ptr(),
                            /* int* */ _ix_wall,
                            /* double& */ _ds_offset,
                            /* bool& */ _is_branch_wall,
                            /* void* */ &_wall3d);
  return PointerToWall3d{
      _ds_offset,
      _is_branch_wall,
      std::move((_wall3d ? std::make_optional<Wall3dStruct>(_wall3d) : std::nullopt))
  };
}
FixedArray1D<Complex, 2> Bmad::polar_to_spinor(SpinPolarStruct &polar) {
  // spinor: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _spinor_desc;
  _spinor_desc.rank = 1;
  FixedArray1D<Complex, 2> _spinor;
  _spinor_desc.data_ptr = _spinor.data();
  _spinor_desc.dims[0] = _spinor.size();
  fortran_polar_to_spinor(/* void* */ polar.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _spinor_desc);
  return _spinor;
}
FixedArray1D<Real, 3> Bmad::polar_to_vec(SpinPolarStruct &polar) {
  // vec: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  FixedArray1D<Real, 3> _vec;
  _vec_desc.data_ptr = _vec.data();
  _vec_desc.dims[0] = _vec.size();
  fortran_polar_to_vec(/* void* */ polar.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _vec_desc);
  return _vec;
}
Bmad::ProjectEmitToXyz Bmad::project_emit_to_xyz(LatStruct &ring, int ix, NormalModesStruct &mode) {
  double _sigma_x{};
  double _sigma_y{};
  double _sigma_z{};
  fortran_project_emit_to_xyz(/* void* */ ring.get_fortran_ptr(),
                              /* int& */ ix,
                              /* void* */ mode.get_fortran_ptr(),
                              /* double& */ _sigma_x,
                              /* double& */ _sigma_y,
                              /* double& */ _sigma_z);
  return ProjectEmitToXyz{_sigma_x, _sigma_y, _sigma_z};
}
double Bmad::psi_prime_sca(double t, double p, FixedArray1D<Real, 8> args) {
  double _dpdt{};
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  fortran_psi_prime_sca(
      /* double& */ t,
      /* double& */ p,
      /* double& */ _dpdt,
      /* Bmad::array_descriptor_t& */ _args_desc
  );
  return _dpdt;
}
void Bmad::ptc_bookkeeper(LatStruct &lat) {
  fortran_ptc_bookkeeper(/* void* */ lat.get_fortran_ptr());
}
void Bmad::ptc_calculate_tracking_step_size(
    Layout &ptc_layout,
    double kl_max,
    std::optional<double> ds_max,
    optional_ref<BoolAlloc1D> even_steps,
    std::optional<double> r_typical,
    std::optional<double> dx_tol_bend,
    std::optional<bool> use_2nd_order,
    std::optional<FixedArray1D<Int, 2>> crossover,
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler
) {
  double ds_max_lvalue;
  auto *_ds_max{&ds_max_lvalue};
  if (ds_max.has_value()) {
    ds_max_lvalue = ds_max.value();
  } else {
    _ds_max = nullptr;
  }
  // intent=in allocatable general array
  auto *_even_steps =
      even_steps.has_value() ? even_steps->get().get_fortran_ptr() : nullptr; // input, optional
  double r_typical_lvalue;
  auto *_r_typical{&r_typical_lvalue};
  if (r_typical.has_value()) {
    r_typical_lvalue = r_typical.value();
  } else {
    _r_typical = nullptr;
  }
  double dx_tol_bend_lvalue;
  auto *_dx_tol_bend{&dx_tol_bend_lvalue};
  if (dx_tol_bend.has_value()) {
    dx_tol_bend_lvalue = dx_tol_bend.value();
  } else {
    _dx_tol_bend = nullptr;
  }
  bool use_2nd_order_lvalue;
  auto *_use_2nd_order{&use_2nd_order_lvalue};
  if (use_2nd_order.has_value()) {
    use_2nd_order_lvalue = use_2nd_order.value();
  } else {
    _use_2nd_order = nullptr;
  }
  // crossover: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _crossover_desc;
  _crossover_desc.rank = 1;
  if (crossover.has_value()) {
    _crossover_desc.data_ptr = crossover->data();
    _crossover_desc.dims[0] = crossover->size();
  } else {
    _crossover_desc.data_ptr = nullptr;
    _crossover_desc.dims[0] = 0;
  }
  // crossover_wiggler: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _crossover_wiggler_desc;
  _crossover_wiggler_desc.rank = 1;
  if (crossover_wiggler.has_value()) {
    _crossover_wiggler_desc.data_ptr = crossover_wiggler->data();
    _crossover_wiggler_desc.dims[0] = crossover_wiggler->size();
  } else {
    _crossover_wiggler_desc.data_ptr = nullptr;
    _crossover_wiggler_desc.dims[0] = 0;
  }
  fortran_ptc_calculate_tracking_step_size(/* void* */ ptc_layout.get_fortran_ptr(),
                                           /* double& */ kl_max,
                                           /* double* */ _ds_max,
                                           /* void* */ _even_steps,
                                           /* double* */ _r_typical,
                                           /* double* */ _dx_tol_bend,
                                           /* bool* */ _use_2nd_order,
                                           /* Bmad::array_descriptor_t& */ _crossover_desc,
                                           /* Bmad::array_descriptor_t& */ _crossover_wiggler_desc);
}
Bmad::PtcCheckForLostParticle Bmad::ptc_check_for_lost_particle(bool do_reset) {
  int _state{};
  void *_ptc_fibre;
  fortran_ptc_check_for_lost_particle(
      /* int& */ _state,
      /* void* */ &_ptc_fibre,
      /* bool& */ do_reset
  );
  return PtcCheckForLostParticle{
      _state,
      std::move((_ptc_fibre ? std::make_optional<Fibre>(_ptc_fibre) : std::nullopt))
  };
}
CoordStructAlloc1D
Bmad::ptc_closed_orbit_calc(BranchStruct &branch, std::optional<bool> radiation_damping_on) {
  // intent=out allocatable type array
  auto closed_orbit{CoordStructAlloc1D()};
  bool radiation_damping_on_lvalue;
  auto *_radiation_damping_on{&radiation_damping_on_lvalue};
  if (radiation_damping_on.has_value()) {
    radiation_damping_on_lvalue = radiation_damping_on.value();
  } else {
    _radiation_damping_on = nullptr;
  }
  fortran_ptc_closed_orbit_calc(/* void* */ branch.get_fortran_ptr(),
                                /* void* */ closed_orbit.get_fortran_ptr(),
                                /* bool* */ _radiation_damping_on);
  return std::move(closed_orbit);
}
Bmad::PtcEmitCalc Bmad::ptc_emit_calc(EleStruct &ele, FixedArray2D<Real, 6, 6> sigma_mat) {
  NormalModesStruct _norm_mode;
  // sigma_mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  CoordStruct _closed_orb;
  fortran_ptc_emit_calc(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ _norm_mode.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
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
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler
) {
  bool even_lvalue;
  auto *_even{&even_lvalue};
  if (even.has_value()) {
    even_lvalue = even.value();
  } else {
    _even = nullptr;
  }
  // crossover: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _crossover_desc;
  _crossover_desc.rank = 1;
  if (crossover.has_value()) {
    _crossover_desc.data_ptr = crossover->data();
    _crossover_desc.dims[0] = crossover->size();
  } else {
    _crossover_desc.data_ptr = nullptr;
    _crossover_desc.dims[0] = 0;
  }
  // crossover_wiggler: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _crossover_wiggler_desc;
  _crossover_wiggler_desc.rank = 1;
  if (crossover_wiggler.has_value()) {
    _crossover_wiggler_desc.data_ptr = crossover_wiggler->data();
    _crossover_wiggler_desc.dims[0] = crossover_wiggler->size();
  } else {
    _crossover_wiggler_desc.data_ptr = nullptr;
    _crossover_wiggler_desc.dims[0] = 0;
  }
  fortran_ptc_layouts_resplit(
      /* double& */ dKL_max,
      /* double& */ l_max,
      /* bool& */ l_max_drift_only,
      /* double& */ bend_dorb,
      /* double& */ sex_dx,
      /* bool* */ _even,
      /* Bmad::array_descriptor_t& */ _crossover_desc,
      /* Bmad::array_descriptor_t& */ _crossover_wiggler_desc
  );
}
void Bmad::ptc_one_turn_mat_and_closed_orbit_calc(BranchStruct &branch, std::optional<double> pz) {
  double pz_lvalue;
  auto *_pz{&pz_lvalue};
  if (pz.has_value()) {
    pz_lvalue = pz.value();
  } else {
    _pz = nullptr;
  }
  fortran_ptc_one_turn_mat_and_closed_orbit_calc(/* void* */ branch.get_fortran_ptr(),
                                                 /* double* */ _pz);
}
void Bmad::ptc_ran_seed_put(int iseed) { fortran_ptc_ran_seed_put(/* int& */ iseed); }
Bmad::PtcReadFlatFile Bmad::ptc_read_flat_file(
    CharacterAlloc1D &flat_file,
    std::optional<bool> create_end_marker,
    std::optional<bool> from_mad
) {
  // intent=in character array container
  bool _err_flag{};
  LatStruct _lat;
  bool create_end_marker_lvalue;
  auto *_create_end_marker{&create_end_marker_lvalue};
  if (create_end_marker.has_value()) {
    create_end_marker_lvalue = create_end_marker.value();
  } else {
    _create_end_marker = nullptr;
  }
  bool from_mad_lvalue;
  auto *_from_mad{&from_mad_lvalue};
  if (from_mad.has_value()) {
    from_mad_lvalue = from_mad.value();
  } else {
    _from_mad = nullptr;
  }
  fortran_ptc_read_flat_file(/* void* */ flat_file.get_fortran_ptr(),
                             /* bool& */ _err_flag,
                             /* void* */ _lat.get_fortran_ptr(),
                             /* bool* */ _create_end_marker,
                             /* bool* */ _from_mad);
  return PtcReadFlatFile{_err_flag, std::move(_lat)};
}
void Bmad::ptc_set_rf_state_for_c_normal(bool nocavity) {
  fortran_ptc_set_rf_state_for_c_normal(/* bool& */ nocavity);
}
void Bmad::ptc_set_taylor_order_if_needed() { fortran_ptc_set_taylor_order_if_needed(); }
Bmad::PtcSpinCalc Bmad::ptc_spin_calc(EleStruct &ele, FixedArray2D<Real, 6, 6> sigma_mat) {
  NormalModesStruct _norm_mode;
  // sigma_mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_desc;
  _sigma_mat_desc.rank = 2;
  double _sigma_mat_vec[6 * 6];
  _sigma_mat_desc.data_ptr = _sigma_mat_vec;
  _sigma_mat_desc.dims[0] = 6;
  _sigma_mat_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat, _sigma_mat_vec);
  CoordStruct _closed_orb;
  fortran_ptc_spin_calc(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ _norm_mode.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _sigma_mat_desc,
                        /* void* */ _closed_orb.get_fortran_ptr());
  vec_to_matrix(_sigma_mat_vec, sigma_mat);
  return PtcSpinCalc{std::move(_norm_mode), std::move(_closed_orb)};
}
Bmad::PtcTrackAll Bmad::ptc_track_all(BranchStruct &branch, CoordStructAlloc1D orbit) {
  // intent=inout allocatable type array
  int _track_state{};
  bool _err_flag{};
  fortran_ptc_track_all(/* void* */ branch.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr(),
                        /* int& */ _track_state,
                        /* bool& */ _err_flag);
  return PtcTrackAll{_track_state, _err_flag};
}
bool Bmad::ptc_transfer_map_with_spin(
    BranchStruct &branch,
    TaylorStructArray1D t_map,
    TaylorStructArray1D s_map,
    CoordStruct &orb0,
    std::optional<int> ix1,
    std::optional<int> ix2,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start
) {
  // t_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _t_map_desc;
  _t_map_desc.rank = 1;
  _t_map_desc.data_ptr = t_map.data();
  _t_map_desc.dims[0] = t_map.size();
  _t_map_desc.strides[0] = 1;
  // s_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _s_map_desc;
  _s_map_desc.rank = 1;
  _s_map_desc.data_ptr = s_map.data();
  _s_map_desc.dims[0] = s_map.size();
  _s_map_desc.strides[0] = 1;
  bool _err_flag{};
  int ix1_lvalue;
  auto *_ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto *_ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  bool one_turn_lvalue;
  auto *_one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto *_unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  fortran_ptc_transfer_map_with_spin(/* void* */ branch.get_fortran_ptr(),
                                     /* Bmad::array_descriptor_t& */ _t_map_desc,
                                     /* Bmad::array_descriptor_t& */ _s_map_desc,
                                     /* void* */ orb0.get_fortran_ptr(),
                                     /* bool& */ _err_flag,
                                     /* int* */ _ix1,
                                     /* int* */ _ix2,
                                     /* bool* */ _one_turn,
                                     /* bool* */ _unit_start);
  return _err_flag;
}
FixedArray2D<Real, 6, 6>
Bmad::pwd_mat(LatStruct &lat, FixedArray2D<Real, 6, 6> t6, double inductance, double sig_z) {
  // t6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_desc;
  _t6_desc.rank = 2;
  double _t6_vec[6 * 6];
  _t6_desc.data_ptr = _t6_vec;
  _t6_desc.dims[0] = 6;
  _t6_desc.dims[1] = 6;
  matrix_to_vec(t6, _t6_vec);
  // t6_pwd: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_pwd_desc;
  _t6_pwd_desc.rank = 2;
  FixedArray2D<Real, 6, 6> t6_pwd;
  double _t6_pwd_vec[6 * 6];
  _t6_pwd_desc.data_ptr = _t6_pwd_vec;
  _t6_pwd_desc.dims[0] = 6;
  _t6_pwd_desc.dims[1] = 6;
  fortran_pwd_mat(/* void* */ lat.get_fortran_ptr(),
                  /* Bmad::array_descriptor_t& */ _t6_desc,
                  /* double& */ inductance,
                  /* double& */ sig_z,
                  /* Bmad::array_descriptor_t& */ _t6_pwd_desc);
  vec_to_matrix(_t6_pwd_vec, t6_pwd);
  return t6_pwd;
}
Bmad::Rad1DampAndStocMats Bmad::rad1_damp_and_stoc_mats(
    EleStruct &ele,
    bool include_opening_angle,
    CoordStruct &orb_in,
    CoordStruct &orb_out,
    double g2_tol,
    double g3_tol,
    optional_ref<EleStruct> ele0
) {
  RadMapStruct _rad_map;
  bool _err_flag{};
  auto *_ele0 = ele0.has_value() ? ele0->get().get_fortran_ptr() : nullptr; // input, optional
  RadInt1Struct _rad_int1;
  fortran_rad1_damp_and_stoc_mats(/* void* */ ele.get_fortran_ptr(),
                                  /* bool& */ include_opening_angle,
                                  /* void* */ orb_in.get_fortran_ptr(),
                                  /* void* */ orb_out.get_fortran_ptr(),
                                  /* void* */ _rad_map.get_fortran_ptr(),
                                  /* double& */ g2_tol,
                                  /* double& */ g3_tol,
                                  /* bool& */ _err_flag,
                                  /* void* */ _ele0,
                                  /* void* */ _rad_int1.get_fortran_ptr());
  return Rad1DampAndStocMats{std::move(_rad_map), _err_flag, std::move(_rad_int1)};
}
Bmad::RadDampAndStocMats Bmad::rad_damp_and_stoc_mats(
    EleStruct &ele1,
    EleStruct &ele2,
    bool include_opening_angle,
    std::optional<CoordStructArray1D> closed_orbit
) {
  RadMapStruct _rmap;
  NormalModesStruct _mode;
  // xfer_nodamp_mat: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _xfer_nodamp_mat_desc;
  _xfer_nodamp_mat_desc.rank = 2;
  FixedArray2D<Real, 6, 6> xfer_nodamp_mat;
  double _xfer_nodamp_mat_vec[6 * 6];
  _xfer_nodamp_mat_desc.data_ptr = _xfer_nodamp_mat_vec;
  _xfer_nodamp_mat_desc.dims[0] = 6;
  _xfer_nodamp_mat_desc.dims[1] = 6;
  bool _err_flag{};
  // closed_orbit: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _closed_orbit_desc;
  _closed_orbit_desc.rank = 1;
  if (closed_orbit) {
    _closed_orbit_desc.data_ptr = closed_orbit->data();
    _closed_orbit_desc.dims[0] = closed_orbit->size();
  } else {
    _closed_orbit_desc.data_ptr = nullptr;
    _closed_orbit_desc.dims[0] = 0;
  }
  _closed_orbit_desc.strides[0] = 1;
  RadIntBranchStruct _rad_int_branch;
  fortran_rad_damp_and_stoc_mats(/* void* */ ele1.get_fortran_ptr(),
                                 /* void* */ ele2.get_fortran_ptr(),
                                 /* bool& */ include_opening_angle,
                                 /* void* */ _rmap.get_fortran_ptr(),
                                 /* void* */ _mode.get_fortran_ptr(),
                                 /* Bmad::array_descriptor_t& */ _xfer_nodamp_mat_desc,
                                 /* bool& */ _err_flag,
                                 /* Bmad::array_descriptor_t& */ _closed_orbit_desc,
                                 /* void* */ _rad_int_branch.get_fortran_ptr());
  vec_to_matrix(_xfer_nodamp_mat_vec, xfer_nodamp_mat);
  return RadDampAndStocMats{
      std::move(_rmap),
      std::move(_mode),
      xfer_nodamp_mat,
      _err_flag,
      std::move(_rad_int_branch)
  };
}
Bmad::RadGIntegrals Bmad::rad_g_integrals(
    EleStruct &ele,
    int where,
    CoordStruct &orb_in,
    CoordStruct &orb_out,
    double int_g2,
    double g_tol,
    double g2_tol,
    double g3_tol
) {
  // int_g: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _int_g_desc;
  _int_g_desc.rank = 1;
  FixedArray1D<Real, 2> _int_g;
  _int_g_desc.data_ptr = _int_g.data();
  _int_g_desc.dims[0] = _int_g.size();
  double _int_g3{};
  fortran_rad_g_integrals(/* void* */ ele.get_fortran_ptr(),
                          /* int& */ where,
                          /* void* */ orb_in.get_fortran_ptr(),
                          /* void* */ orb_out.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _int_g_desc,
                          /* double& */ int_g2,
                          /* double& */ _int_g3,
                          /* double& */ g_tol,
                          /* double& */ g2_tol,
                          /* double& */ g3_tol);
  return RadGIntegrals{_int_g, _int_g3};
}
Bmad::RadiationIntegrals Bmad::radiation_integrals(
    LatStruct &lat,
    CoordStructArray1D orbit,
    optional_ref<int> ix_cache,
    std::optional<int> ix_branch
) {
  // orbit: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  _orbit_desc.data_ptr = orbit.data();
  _orbit_desc.dims[0] = orbit.size();
  _orbit_desc.strides[0] = 1;
  NormalModesStruct _mode;
  auto *_ix_cache = ix_cache.has_value() ? &ix_cache->get() : nullptr; // inout, optional
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  RadIntAllEleStruct _rad_int_by_ele;
  fortran_radiation_integrals(/* void* */ lat.get_fortran_ptr(),
                              /* Bmad::array_descriptor_t& */ _orbit_desc,
                              /* void* */ _mode.get_fortran_ptr(),
                              /* int* */ _ix_cache,
                              /* int* */ _ix_branch,
                              /* void* */ _rad_int_by_ele.get_fortran_ptr());
  return RadiationIntegrals{std::move(_mode), std::move(_rad_int_by_ele)};
}
bool Bmad::radiation_map_setup(EleStruct &ele, optional_ref<CoordStruct> ref_orbit_in) {
  bool _err_flag{};
  auto *_ref_orbit_in =
      ref_orbit_in.has_value() ? ref_orbit_in->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_radiation_map_setup(/* void* */ ele.get_fortran_ptr(),
                              /* bool& */ _err_flag,
                              /* void* */ _ref_orbit_in);
  return _err_flag;
}
void Bmad::ramper_slave_setup(LatStruct &lat, std::optional<bool> force_setup) {
  bool force_setup_lvalue;
  auto *_force_setup{&force_setup_lvalue};
  if (force_setup.has_value()) {
    force_setup_lvalue = force_setup.value();
  } else {
    _force_setup = nullptr;
  }
  fortran_ramper_slave_setup(/* void* */ lat.get_fortran_ptr(), /* bool* */ _force_setup);
}
Bmad::RamperValue Bmad::ramper_value(EleStruct &ramper, ControlRamp1Struct &r1) {
  bool _err_flag{};
  double _value{};
  fortran_ramper_value(/* void* */ ramper.get_fortran_ptr(),
                       /* void* */ r1.get_fortran_ptr(),
                       /* bool& */ _err_flag,
                       /* double& */ _value);
  return RamperValue{_err_flag, _value};
}
bool Bmad::randomize_lr_wake_frequencies(EleStruct &ele) {
  bool _set_done{};
  fortran_randomize_lr_wake_frequencies(/* void* */ ele.get_fortran_ptr(), /* bool& */ _set_done);
  return _set_done;
}
void Bmad::rchomp(double rel, int plc, std::string out) {
  auto _out = out.c_str();
  fortran_rchomp(/* double& */ rel, /* int& */ plc, /* const char* */ _out);
}
void Bmad::re_allocate_eles(
    ElePointerStructAlloc1D eles,
    int n,
    std::optional<bool> save_old,
    std::optional<bool> exact
) {
  // intent=inout allocatable type array
  bool save_old_lvalue;
  auto *_save_old{&save_old_lvalue};
  if (save_old.has_value()) {
    save_old_lvalue = save_old.value();
  } else {
    _save_old = nullptr;
  }
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_allocate_eles(/* void* */ eles.get_fortran_ptr(),
                           /* int& */ n,
                           /* bool* */ _save_old,
                           /* bool* */ _exact);
}
void Bmad::re_allocate(Wall3dSectionStructAlloc1D section, int n, std::optional<bool> exact) {
  // intent=inout allocatable type array
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_allocate_wall3d_section_array(/* void* */ section.get_fortran_ptr(),
                                           /* int& */ n,
                                           /* bool* */ _exact);
}
void Bmad::re_allocate(Wall3dVertexStructAlloc1D v, int n, std::optional<bool> exact) {
  // intent=inout allocatable type array
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_allocate_wall3d_vertex_array(/* void* */ v.get_fortran_ptr(),
                                          /* int& */ n,
                                          /* bool* */ _exact);
}
void Bmad::re_associate_node_array(ExpressionTreeStruct &tree, int n, std::optional<bool> exact) {
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_re_associate_node_array(/* void* */ tree.get_fortran_ptr(),
                                  /* int& */ n,
                                  /* bool* */ _exact);
}
void Bmad::re_str(long double rel, std::string str_out) {
  auto _str_out = str_out.c_str();
  fortran_re_str_qp(/* long double& */ rel, /* const char* */ _str_out);
}
void Bmad::re_str(double rel, std::string str_out) {
  auto _str_out = str_out.c_str();
  fortran_re_str_rp(/* double& */ rel, /* const char* */ _str_out);
}
Bmad::ReadBeamAscii Bmad::read_beam_ascii(std::string file_name, BeamInitStruct &beam_init) {
  auto _file_name = file_name.c_str();
  BeamStruct _beam;
  bool _err_flag{};
  fortran_read_beam_ascii(
      /* const char* */ _file_name,
      /* void* */ _beam.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return ReadBeamAscii{std::move(_beam), _err_flag};
}
Bmad::ReadBeamFile Bmad::read_beam_file(
    std::string file_name,
    BeamInitStruct &beam_init,
    optional_ref<EleStruct> ele,
    std::optional<bool> print_mom_shift_warning,
    std::optional<bool> conserve_momentum
) {
  auto _file_name = file_name.c_str();
  BeamStruct _beam;
  bool _err_flag{};
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  bool print_mom_shift_warning_lvalue;
  auto *_print_mom_shift_warning{&print_mom_shift_warning_lvalue};
  if (print_mom_shift_warning.has_value()) {
    print_mom_shift_warning_lvalue = print_mom_shift_warning.value();
  } else {
    _print_mom_shift_warning = nullptr;
  }
  bool conserve_momentum_lvalue;
  auto *_conserve_momentum{&conserve_momentum_lvalue};
  if (conserve_momentum.has_value()) {
    conserve_momentum_lvalue = conserve_momentum.value();
  } else {
    _conserve_momentum = nullptr;
  }
  fortran_read_beam_file(
      /* const char* */ _file_name,
      /* void* */ _beam.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _ele,
      /* bool* */ _print_mom_shift_warning,
      /* bool* */ _conserve_momentum
  );
  return ReadBeamFile{std::move(_beam), _err_flag};
}
Bmad::ReadBinaryCartesianMap
Bmad::read_binary_cartesian_map(std::string file_name, EleStruct &ele) {
  auto _file_name = file_name.c_str();
  CartesianMapStruct _cart_map;
  bool _err_flag{};
  fortran_read_binary_cartesian_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _cart_map.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return ReadBinaryCartesianMap{std::move(_cart_map), _err_flag};
}
Bmad::ReadBinaryCylindricalMap
Bmad::read_binary_cylindrical_map(std::string file_name, EleStruct &ele) {
  auto _file_name = file_name.c_str();
  CylindricalMapStruct _cl_map;
  bool _err_flag{};
  fortran_read_binary_cylindrical_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _cl_map.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return ReadBinaryCylindricalMap{std::move(_cl_map), _err_flag};
}
Bmad::ReadBinaryGridField Bmad::read_binary_grid_field(std::string file_name, EleStruct &ele) {
  auto _file_name = file_name.c_str();
  GridFieldStruct _g_field;
  bool _err_flag{};
  fortran_read_binary_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _g_field.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return ReadBinaryGridField{std::move(_g_field), _err_flag};
}
Bmad::ReadDigestedBmadFile Bmad::read_digested_bmad_file(std::string digested_file) {
  auto _digested_file = digested_file.c_str();
  LatStruct _lat;
  int _inc_version{};
  bool _err_flag{};
  bool _parser_calling{};
  // intent=out character array container
  auto lat_files{CharacterAlloc1D()};
  fortran_read_digested_bmad_file(
      /* const char* */ _digested_file,
      /* void* */ _lat.get_fortran_ptr(),
      /* int& */ _inc_version,
      /* bool& */ _err_flag,
      /* bool& */ _parser_calling,
      /* void* */ lat_files.get_fortran_ptr()
  );
  return ReadDigestedBmadFile{
      std::move(_lat),
      _inc_version,
      _err_flag,
      _parser_calling,
      std::move(lat_files)
  };
}
PhotonReflectSurfaceStruct Bmad::read_surface_reflection_file(std::string file_name) {
  auto _file_name = file_name.c_str();
  PhotonReflectSurfaceStruct _surface;
  fortran_read_surface_reflection_file(
      /* const char* */ _file_name,
      /* void* */ _surface.get_fortran_ptr()
  );
  return std::move(_surface);
}
void Bmad::reallocate_beam(
    BeamStruct &beam,
    int n_bunch,
    std::optional<int> n_particle,
    std::optional<bool> extend
) {
  int n_particle_lvalue;
  auto *_n_particle{&n_particle_lvalue};
  if (n_particle.has_value()) {
    n_particle_lvalue = n_particle.value();
  } else {
    _n_particle = nullptr;
  }
  bool extend_lvalue;
  auto *_extend{&extend_lvalue};
  if (extend.has_value()) {
    extend_lvalue = extend.value();
  } else {
    _extend = nullptr;
  }
  fortran_reallocate_beam(/* void* */ beam.get_fortran_ptr(),
                          /* int& */ n_bunch,
                          /* int* */ _n_particle,
                          /* bool* */ _extend);
}
void Bmad::reallocate_bp_com_const() { fortran_reallocate_bp_com_const(); }
BunchStruct Bmad::reallocate_bunch(int n_particle, std::optional<bool> save) {
  BunchStruct _bunch;
  bool save_lvalue;
  auto *_save{&save_lvalue};
  if (save.has_value()) {
    save_lvalue = save.value();
  } else {
    _save = nullptr;
  }
  fortran_reallocate_bunch(/* void* */ _bunch.get_fortran_ptr(),
                           /* int& */ n_particle,
                           /* bool* */ _save);
  return std::move(_bunch);
}
void Bmad::reallocate_control(LatStruct &lat, int n) {
  fortran_reallocate_control(/* void* */ lat.get_fortran_ptr(), /* int& */ n);
}
void Bmad::reallocate_coord(CoordArrayStructAlloc1D coord_array, LatStruct &lat) {
  // intent=inout allocatable type array
  fortran_reallocate_coord_array(/* void* */ coord_array.get_fortran_ptr(),
                                 /* void* */ lat.get_fortran_ptr());
}
void Bmad::reallocate_coord(
    CoordStructAlloc1D coord,
    LatStruct &lat,
    std::optional<int> ix_branch
) {
  // intent=inout allocatable type array
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  fortran_reallocate_coord_lat(/* void* */ coord.get_fortran_ptr(),
                               /* void* */ lat.get_fortran_ptr(),
                               /* int* */ _ix_branch);
}
void Bmad::reallocate_coord(CoordStructAlloc1D coord, int n_coord) {
  // intent=inout allocatable type array
  fortran_reallocate_coord_n(/* void* */ coord.get_fortran_ptr(), /* int& */ n_coord);
}
void Bmad::reallocate_expression_stack(
    ExpressionAtomStructAlloc1D stack,
    int n,
    std::optional<bool> exact
) {
  // intent=inout allocatable type array
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_reallocate_expression_stack(/* void* */ stack.get_fortran_ptr(),
                                      /* int& */ n,
                                      /* bool* */ _exact);
}
double Bmad::rel_tracking_charge_to_mass(CoordStruct &orbit, int ref_species) {
  double _rel_charge{};
  fortran_rel_tracking_charge_to_mass(/* void* */ orbit.get_fortran_ptr(),
                                      /* int& */ ref_species,
                                      /* double& */ _rel_charge);
  return _rel_charge;
}
void Bmad::relative_mode_flip(EleStruct &ele1, EleStruct &ele2, bool func_retval__) {
  fortran_relative_mode_flip(/* void* */ ele1.get_fortran_ptr(),
                             /* void* */ ele2.get_fortran_ptr(),
                             /* bool& */ func_retval__);
}
void Bmad::release_rad_int_cache(int &ix_cache) {
  fortran_release_rad_int_cache(/* int& */ ix_cache);
}
void Bmad::remove_constant_taylor(
    TaylorStructArray1D taylor_in,
    TaylorStructArray1D taylor_out,
    FArray1D<Real> &c0,
    bool remove_higher_order_terms
) {
  // taylor_in: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_in_desc;
  _taylor_in_desc.rank = 1;
  _taylor_in_desc.data_ptr = taylor_in.data();
  _taylor_in_desc.dims[0] = taylor_in.size();
  _taylor_in_desc.strides[0] = 1;
  // taylor_out: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_out_desc;
  _taylor_out_desc.rank = 1;
  _taylor_out_desc.data_ptr = taylor_out.data();
  _taylor_out_desc.dims[0] = taylor_out.size();
  _taylor_out_desc.strides[0] = 1;
  // c0: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _c0_desc;
  _c0_desc.rank = 1;
  _c0_desc.data_ptr = c0.data();
  _c0_desc.dims[0] = c0.size();
  fortran_remove_constant_taylor(
      /* Bmad::array_descriptor_t& */ _taylor_in_desc,
      /* Bmad::array_descriptor_t& */ _taylor_out_desc,
      /* Bmad::array_descriptor_t& */ _c0_desc,
      /* bool& */ remove_higher_order_terms
  );
}
BunchStruct Bmad::remove_dead_from_bunch(BunchStruct &bunch_in) {
  BunchStruct _bunch_out;
  fortran_remove_dead_from_bunch(/* void* */ bunch_in.get_fortran_ptr(),
                                 /* void* */ _bunch_out.get_fortran_ptr());
  return std::move(_bunch_out);
}
void Bmad::remove_eles_from_lat(LatStruct &lat, std::optional<bool> check_sanity) {
  bool check_sanity_lvalue;
  auto *_check_sanity{&check_sanity_lvalue};
  if (check_sanity.has_value()) {
    check_sanity_lvalue = check_sanity.value();
  } else {
    _check_sanity = nullptr;
  }
  fortran_remove_eles_from_lat(/* void* */ lat.get_fortran_ptr(), /* bool* */ _check_sanity);
}
void Bmad::remove_lord_slave_link(EleStruct &lord, EleStruct &slave) {
  fortran_remove_lord_slave_link(/* void* */ lord.get_fortran_ptr(),
                                 /* void* */ slave.get_fortran_ptr());
}
LatStruct Bmad::reverse_lat(LatStruct &lat_in, std::optional<bool> track_antiparticle) {
  LatStruct _lat_rev;
  bool track_antiparticle_lvalue;
  auto *_track_antiparticle{&track_antiparticle_lvalue};
  if (track_antiparticle.has_value()) {
    track_antiparticle_lvalue = track_antiparticle.value();
  } else {
    _track_antiparticle = nullptr;
  }
  fortran_reverse_lat(/* void* */ lat_in.get_fortran_ptr(),
                      /* void* */ _lat_rev.get_fortran_ptr(),
                      /* bool* */ _track_antiparticle);
  return std::move(_lat_rev);
}
void Bmad::rf_coupler_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    double phase,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_rf_coupler_kick(/* void* */ ele.get_fortran_ptr(),
                          /* void* */ param.get_fortran_ptr(),
                          /* int& */ particle_at,
                          /* double& */ phase,
                          /* void* */ orbit.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                          /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
bool Bmad::rf_is_on(BranchStruct &branch, std::optional<int> ix_ele1, std::optional<int> ix_ele2) {
  int ix_ele1_lvalue;
  auto *_ix_ele1{&ix_ele1_lvalue};
  if (ix_ele1.has_value()) {
    ix_ele1_lvalue = ix_ele1.value();
  } else {
    _ix_ele1 = nullptr;
  }
  int ix_ele2_lvalue;
  auto *_ix_ele2{&ix_ele2_lvalue};
  if (ix_ele2.has_value()) {
    ix_ele2_lvalue = ix_ele2.value();
  } else {
    _ix_ele2 = nullptr;
  }
  bool _is_on{};
  fortran_rf_is_on(/* void* */ branch.get_fortran_ptr(),
                   /* int* */ _ix_ele1,
                   /* int* */ _ix_ele2,
                   /* bool& */ _is_on);
  return _is_on;
}
double Bmad::rf_ref_time_offset(EleStruct &ele, std::optional<double> ds) {
  double ds_lvalue;
  auto *_ds{&ds_lvalue};
  if (ds.has_value()) {
    ds_lvalue = ds.value();
  } else {
    _ds = nullptr;
  }
  double _time{};
  fortran_rf_ref_time_offset(/* void* */ ele.get_fortran_ptr(),
                             /* double* */ _ds,
                             /* double& */ _time);
  return _time;
}
void Bmad::rfun(
    double u,
    double v,
    double w,
    double gam,
    double a,
    double b,
    double hz,
    int i,
    int j,
    double res
) {
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
      /* double& */ res
  );
}
void Bmad::rk_adaptive_time_step(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orb,
    int t_dir,
    double rf_time,
    double dt_try,
    double dt_did,
    double dt_next,
    bool err_flag,
    optional_ref<EmFieldStruct> extra_field
) {
  auto *_extra_field =
      extra_field.has_value() ? extra_field->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_rk_adaptive_time_step(/* void* */ ele.get_fortran_ptr(),
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
    EleStruct &ele,
    LatParamStruct &param,
    double rf_time,
    CoordStruct &orb,
    double dt,
    CoordStruct &new_orb,
    bool err_flag,
    std::optional<FixedArray1D<Real, 10>> dr_dt,
    std::optional<bool> print_err,
    optional_ref<EmFieldStruct> extra_field
) {
  // r_err: out NOT (CppWrapperGeneralArgumentArray) (['10'])
  Bmad::array_descriptor_t _r_err_desc;
  _r_err_desc.rank = 1;
  FixedArray1D<Real, 10> _r_err;
  _r_err_desc.data_ptr = _r_err.data();
  _r_err_desc.dims[0] = _r_err.size();
  // dr_dt: in NOT (CppWrapperGeneralArgumentArray) (['10'])
  Bmad::array_descriptor_t _dr_dt_desc;
  _dr_dt_desc.rank = 1;
  if (dr_dt.has_value()) {
    _dr_dt_desc.data_ptr = dr_dt->data();
    _dr_dt_desc.dims[0] = dr_dt->size();
  } else {
    _dr_dt_desc.data_ptr = nullptr;
    _dr_dt_desc.dims[0] = 0;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  auto *_extra_field =
      extra_field.has_value() ? extra_field->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_rk_time_step1(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* double& */ rf_time,
                        /* void* */ orb.get_fortran_ptr(),
                        /* double& */ dt,
                        /* void* */ new_orb.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _r_err_desc,
                        /* Bmad::array_descriptor_t& */ _dr_dt_desc,
                        /* bool& */ err_flag,
                        /* bool* */ _print_err,
                        /* void* */ _extra_field);
  return _r_err;
}
void Bmad::rotate3(FixedArray1D<Real, 3> vec, double angle, FixedArray1D<Real, 3> rvec) {
  // vec: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  // rvec: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _rvec_desc;
  _rvec_desc.rank = 1;
  _rvec_desc.data_ptr = rvec.data();
  _rvec_desc.dims[0] = rvec.size();
  fortran_rotate3(
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* double& */ angle,
      /* Bmad::array_descriptor_t& */ _rvec_desc
  );
}
void Bmad::rotate_em_field(
    EmFieldStruct &field,
    FixedArray2D<Real, 3, 3> w_mat,
    FixedArray2D<Real, 3, 3> w_inv,
    std::optional<bool> calc_dfield,
    std::optional<bool> calc_potential
) {
  // w_mat: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  matrix_to_vec(w_mat, _w_mat_vec);
  // w_inv: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_inv_desc;
  _w_inv_desc.rank = 2;
  double _w_inv_vec[3 * 3];
  _w_inv_desc.data_ptr = _w_inv_vec;
  _w_inv_desc.dims[0] = 3;
  _w_inv_desc.dims[1] = 3;
  matrix_to_vec(w_inv, _w_inv_vec);
  bool calc_dfield_lvalue;
  auto *_calc_dfield{&calc_dfield_lvalue};
  if (calc_dfield.has_value()) {
    calc_dfield_lvalue = calc_dfield.value();
  } else {
    _calc_dfield = nullptr;
  }
  bool calc_potential_lvalue;
  auto *_calc_potential{&calc_potential_lvalue};
  if (calc_potential.has_value()) {
    calc_potential_lvalue = calc_potential.value();
  } else {
    _calc_potential = nullptr;
  }
  fortran_rotate_em_field(/* void* */ field.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _w_mat_desc,
                          /* Bmad::array_descriptor_t& */ _w_inv_desc,
                          /* bool* */ _calc_dfield,
                          /* bool* */ _calc_potential);
}
void Bmad::rotate_field_zx(EmFieldStruct &field, double theta) {
  fortran_rotate_field_zx(/* void* */ field.get_fortran_ptr(), /* double& */ theta);
}
void Bmad::rotate_for_curved_surface(
    EleStruct &ele,
    CoordStruct &orbit,
    bool set,
    FixedArray2D<Real, 3, 3> rot_mat
) {
  // rot_mat: inout NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _rot_mat_desc;
  _rot_mat_desc.rank = 2;
  double _rot_mat_vec[3 * 3];
  _rot_mat_desc.data_ptr = _rot_mat_vec;
  _rot_mat_desc.dims[0] = 3;
  _rot_mat_desc.dims[1] = 3;
  matrix_to_vec(rot_mat, _rot_mat_vec);
  fortran_rotate_for_curved_surface(/* void* */ ele.get_fortran_ptr(),
                                    /* void* */ orbit.get_fortran_ptr(),
                                    /* bool& */ set,
                                    /* Bmad::array_descriptor_t& */ _rot_mat_desc);
  vec_to_matrix(_rot_mat_vec, rot_mat);
}
FixedArray1D<Real, 4> Bmad::rotate_spin(FixedArray1D<Real, 3> rot_vec, FixedArray1D<Real, 3> spin) {
  // rot_vec: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _rot_vec_desc;
  _rot_vec_desc.rank = 1;
  _rot_vec_desc.data_ptr = rot_vec.data();
  _rot_vec_desc.dims[0] = rot_vec.size();
  // spin: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _spin_desc;
  _spin_desc.rank = 1;
  _spin_desc.data_ptr = spin.data();
  _spin_desc.dims[0] = spin.size();
  // qrot: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _qrot_desc;
  _qrot_desc.rank = 1;
  FixedArray1D<Real, 4> _qrot;
  _qrot_desc.data_ptr = _qrot.data();
  _qrot_desc.dims[0] = _qrot.size();
  fortran_rotate_spin(
      /* Bmad::array_descriptor_t& */ _rot_vec_desc,
      /* Bmad::array_descriptor_t& */ _spin_desc,
      /* Bmad::array_descriptor_t& */ _qrot_desc
  );
  return _qrot;
}
void Bmad::rotate_spin_a_step(CoordStruct &orbit, EmFieldStruct &field, EleStruct &ele, double ds) {
  fortran_rotate_spin_a_step(/* void* */ orbit.get_fortran_ptr(),
                             /* void* */ field.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* double& */ ds);
}
void Bmad::rotate_spin_given_field(
    CoordStruct &orbit,
    int sign_z_vel,
    std::optional<FixedArray1D<Real, 3>> BL,
    std::optional<FixedArray1D<Real, 3>> EL,
    std::optional<FixedArray1D<Real, 4>> qrot
) {
  // BL: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _BL_desc;
  _BL_desc.rank = 1;
  if (BL.has_value()) {
    _BL_desc.data_ptr = BL->data();
    _BL_desc.dims[0] = BL->size();
  } else {
    _BL_desc.data_ptr = nullptr;
    _BL_desc.dims[0] = 0;
  }
  // EL: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _EL_desc;
  _EL_desc.rank = 1;
  if (EL.has_value()) {
    _EL_desc.data_ptr = EL->data();
    _EL_desc.dims[0] = EL->size();
  } else {
    _EL_desc.data_ptr = nullptr;
    _EL_desc.dims[0] = 0;
  }
  // qrot: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _qrot_desc;
  _qrot_desc.rank = 1;
  if (qrot.has_value()) {
    _qrot_desc.data_ptr = qrot->data();
    _qrot_desc.dims[0] = qrot->size();
  } else {
    _qrot_desc.data_ptr = nullptr;
    _qrot_desc.dims[0] = 0;
  }
  fortran_rotate_spin_given_field(/* void* */ orbit.get_fortran_ptr(),
                                  /* int& */ sign_z_vel,
                                  /* Bmad::array_descriptor_t& */ _BL_desc,
                                  /* Bmad::array_descriptor_t& */ _EL_desc,
                                  /* Bmad::array_descriptor_t& */ _qrot_desc);
}
double Bmad::s_body_calc(CoordStruct &orbit, EleStruct &ele) {
  double _s_body{};
  fortran_s_body_calc(/* void* */ orbit.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* double& */ _s_body);
  return _s_body;
}
void Bmad::s_calc(LatStruct &lat) { fortran_s_calc(/* void* */ lat.get_fortran_ptr()); }
void Bmad::sad_mult_hard_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sad_mult_hard_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                                       /* void* */ param.get_fortran_ptr(),
                                       /* int& */ particle_at,
                                       /* void* */ orbit.get_fortran_ptr(),
                                       /* Bmad::array_descriptor_t& */ _mat6_desc,
                                       /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::sad_soft_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sad_soft_bend_edge_kick(/* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* int& */ particle_at,
                                  /* void* */ orb.get_fortran_ptr(),
                                  /* Bmad::array_descriptor_t& */ _mat6_desc,
                                  /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::save_a_beam_step(
    EleStruct &ele,
    BeamStruct &beam,
    std::optional<BunchTrackStructArray1D> bunch_tracks,
    std::optional<double> s_body,
    std::optional<bool> is_time_coords
) {
  // bunch_tracks: BunchTrackStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _bunch_tracks_desc;
  _bunch_tracks_desc.rank = 1;
  if (bunch_tracks) {
    _bunch_tracks_desc.data_ptr = bunch_tracks->data();
    _bunch_tracks_desc.dims[0] = bunch_tracks->size();
  } else {
    _bunch_tracks_desc.data_ptr = nullptr;
    _bunch_tracks_desc.dims[0] = 0;
  }
  _bunch_tracks_desc.strides[0] = 1;
  double s_body_lvalue;
  auto *_s_body{&s_body_lvalue};
  if (s_body.has_value()) {
    s_body_lvalue = s_body.value();
  } else {
    _s_body = nullptr;
  }
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  fortran_save_a_beam_step(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ beam.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _bunch_tracks_desc,
                           /* double* */ _s_body,
                           /* bool* */ _is_time_coords);
}
void Bmad::save_a_bunch_step(
    EleStruct &ele,
    BunchStruct &bunch,
    optional_ref<BunchTrackStruct> bunch_track,
    std::optional<double> s_body,
    std::optional<bool> is_time_coords
) {
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  double s_body_lvalue;
  auto *_s_body{&s_body_lvalue};
  if (s_body.has_value()) {
    s_body_lvalue = s_body.value();
  } else {
    _s_body = nullptr;
  }
  bool is_time_coords_lvalue;
  auto *_is_time_coords{&is_time_coords_lvalue};
  if (is_time_coords.has_value()) {
    is_time_coords_lvalue = is_time_coords.value();
  } else {
    _is_time_coords = nullptr;
  }
  fortran_save_a_bunch_step(/* void* */ ele.get_fortran_ptr(),
                            /* void* */ bunch.get_fortran_ptr(),
                            /* void* */ _bunch_track,
                            /* double* */ _s_body,
                            /* bool* */ _is_time_coords);
}
void Bmad::save_a_step(
    TrackStruct &track,
    EleStruct &ele,
    LatParamStruct &param,
    bool local_ref_frame,
    CoordStruct &orb,
    double s_rel,
    std::optional<bool> save_field,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<double> rf_time,
    optional_ref<StrongBeamStruct> strong_beam
) {
  bool save_field_lvalue;
  auto *_save_field{&save_field_lvalue};
  if (save_field.has_value()) {
    save_field_lvalue = save_field.value();
  } else {
    _save_field = nullptr;
  }
  // mat6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  double rf_time_lvalue;
  auto *_rf_time{&rf_time_lvalue};
  if (rf_time.has_value()) {
    rf_time_lvalue = rf_time.value();
  } else {
    _rf_time = nullptr;
  }
  auto *_strong_beam =
      strong_beam.has_value() ? strong_beam->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_save_a_step(/* void* */ track.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* void* */ param.get_fortran_ptr(),
                      /* bool& */ local_ref_frame,
                      /* void* */ orb.get_fortran_ptr(),
                      /* double& */ s_rel,
                      /* bool* */ _save_field,
                      /* Bmad::array_descriptor_t& */ _mat6_desc,
                      /* bool* */ _make_matrix,
                      /* double* */ _rf_time,
                      /* void* */ _strong_beam);
}
void Bmad::sbend_body_with_k1_map(
    EleStruct &ele,
    double dg,
    double b1,
    LatParamStruct &param,
    int n_step,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_sbend_body_with_k1_map(/* void* */ ele.get_fortran_ptr(),
                                 /* double& */ dg,
                                 /* double& */ b1,
                                 /* void* */ param.get_fortran_ptr(),
                                 /* int& */ n_step,
                                 /* void* */ orbit.get_fortran_ptr(),
                                 /* Bmad::array_descriptor_t& */ _mat6_desc,
                                 /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
double Bmad::sc_adaptive_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool &include_image,
    double t_now,
    double &dt_step,
    EmFieldStructArray1D sc_field
) {
  double _dt_next{};
  // sc_field: EmFieldStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _sc_field_desc;
  _sc_field_desc.rank = 1;
  _sc_field_desc.data_ptr = sc_field.data();
  _sc_field_desc.dims[0] = sc_field.size();
  _sc_field_desc.strides[0] = 1;
  fortran_sc_adaptive_step(/* void* */ bunch.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* bool& */ include_image,
                           /* double& */ t_now,
                           /* double& */ dt_step,
                           /* double& */ _dt_next,
                           /* Bmad::array_descriptor_t& */ _sc_field_desc);
  return _dt_next;
}
int Bmad::sc_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool &include_image,
    double t_end,
    EmFieldStructArray1D sc_field
) {
  // sc_field: EmFieldStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _sc_field_desc;
  _sc_field_desc.rank = 1;
  _sc_field_desc.data_ptr = sc_field.data();
  _sc_field_desc.dims[0] = sc_field.size();
  _sc_field_desc.strides[0] = 1;
  int _n_emit{};
  fortran_sc_step(/* void* */ bunch.get_fortran_ptr(),
                  /* void* */ ele.get_fortran_ptr(),
                  /* bool& */ include_image,
                  /* double& */ t_end,
                  /* Bmad::array_descriptor_t& */ _sc_field_desc,
                  /* int& */ _n_emit);
  return _n_emit;
}
CoordStruct Bmad::set_active_fixer(EleStruct &fixer, std::optional<bool> turn_on) {
  bool turn_on_lvalue;
  auto *_turn_on{&turn_on_lvalue};
  if (turn_on.has_value()) {
    turn_on_lvalue = turn_on.value();
  } else {
    _turn_on = nullptr;
  }
  CoordStruct _orbit;
  fortran_set_active_fixer(/* void* */ fixer.get_fortran_ptr(),
                           /* bool* */ _turn_on,
                           /* void* */ _orbit.get_fortran_ptr());
  return std::move(_orbit);
}
bool Bmad::set_custom_attribute_name(std::string custom_name, std::optional<int> custom_index) {
  auto _custom_name = custom_name.c_str();
  bool _err_flag{};
  int custom_index_lvalue;
  auto *_custom_index{&custom_index_lvalue};
  if (custom_index.has_value()) {
    custom_index_lvalue = custom_index.value();
  } else {
    _custom_index = nullptr;
  }
  fortran_set_custom_attribute_name(
      /* const char* */ _custom_name,
      /* bool& */ _err_flag,
      /* int* */ _custom_index
  );
  return _err_flag;
}
Bmad::SetEleAttribute Bmad::set_ele_attribute(
    EleStruct &ele,
    std::string set_string,
    std::optional<bool> err_print_flag,
    std::optional<bool> set_lords
) {
  auto _set_string = set_string.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  bool set_lords_lvalue;
  auto *_set_lords{&set_lords_lvalue};
  if (set_lords.has_value()) {
    set_lords_lvalue = set_lords.value();
  } else {
    _set_lords = nullptr;
  }
  int _err_id{};
  fortran_set_ele_attribute(/* void* */ ele.get_fortran_ptr(),
                            /* const char* */ _set_string,
                            /* bool& */ _err_flag,
                            /* bool* */ _err_print_flag,
                            /* bool* */ _set_lords,
                            /* int& */ _err_id);
  return SetEleAttribute{_err_flag, _err_id};
}
void Bmad::set_ele_defaults(EleStruct &ele, std::optional<bool> do_allocate) {
  bool do_allocate_lvalue;
  auto *_do_allocate{&do_allocate_lvalue};
  if (do_allocate.has_value()) {
    do_allocate_lvalue = do_allocate.value();
  } else {
    _do_allocate = nullptr;
  }
  fortran_set_ele_defaults(/* void* */ ele.get_fortran_ptr(), /* bool* */ _do_allocate);
}
void Bmad::set_ele_name(EleStruct &ele, std::string name) {
  auto _name = name.c_str();
  fortran_set_ele_name(/* void* */ ele.get_fortran_ptr(), /* const char* */ _name);
}
bool Bmad::set_ele_real_attribute(
    EleStruct &ele,
    std::string attrib_name,
    double value,
    std::optional<bool> err_print_flag
) {
  auto _attrib_name = attrib_name.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  fortran_set_ele_real_attribute(/* void* */ ele.get_fortran_ptr(),
                                 /* const char* */ _attrib_name,
                                 /* double& */ value,
                                 /* bool& */ _err_flag,
                                 /* bool* */ _err_print_flag);
  return _err_flag;
}
void Bmad::set_ele_status_stale(
    EleStruct &ele,
    int status_group,
    std::optional<bool> set_slaves,
    std::optional<ElePointerStructAlloc1D> old_eles
) {
  bool set_slaves_lvalue;
  auto *_set_slaves{&set_slaves_lvalue};
  if (set_slaves.has_value()) {
    set_slaves_lvalue = set_slaves.value();
  } else {
    _set_slaves = nullptr;
  }
  // intent=in allocatable type array
  auto *_old_eles = old_eles.has_value() ? old_eles->get_fortran_ptr() : nullptr; // input, optional
  fortran_set_ele_status_stale(/* void* */ ele.get_fortran_ptr(),
                               /* int& */ status_group,
                               /* bool* */ _set_slaves,
                               /* void* */ _old_eles);
}
void Bmad::set_flags_for_changed_attribute(
    EleStruct &ele,
    int attrib,
    std::optional<bool> set_dependent
) {
  bool set_dependent_lvalue;
  auto *_set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_integer_attribute(/* void* */ ele.get_fortran_ptr(),
                                                  /* int& */ attrib,
                                                  /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(LatStruct &lat, std::optional<bool> set_dependent) {
  bool set_dependent_lvalue;
  auto *_set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_lat_attribute(/* void* */ lat.get_fortran_ptr(),
                                              /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(
    EleStruct &ele,
    bool attrib,
    std::optional<bool> set_dependent
) {
  bool set_dependent_lvalue;
  auto *_set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_logical_attribute(/* void* */ ele.get_fortran_ptr(),
                                                  /* bool& */ attrib,
                                                  /* bool* */ _set_dependent);
}
void Bmad::set_flags_for_changed_attribute(
    EleStruct &ele,
    std::optional<double> attrib,
    std::optional<bool> set_dependent
) {
  double attrib_lvalue;
  auto *_attrib{&attrib_lvalue};
  if (attrib.has_value()) {
    attrib_lvalue = attrib.value();
  } else {
    _attrib = nullptr;
  }
  bool set_dependent_lvalue;
  auto *_set_dependent{&set_dependent_lvalue};
  if (set_dependent.has_value()) {
    set_dependent_lvalue = set_dependent.value();
  } else {
    _set_dependent = nullptr;
  }
  fortran_set_flags_for_changed_real_attribute(/* void* */ ele.get_fortran_ptr(),
                                               /* double* */ _attrib,
                                               /* bool* */ _set_dependent);
}
void Bmad::set_fringe_on_off(double &fringe_at, int ele_end, int on_or_off) {
  fortran_set_fringe_on_off(/* double& */ fringe_at, /* int& */ ele_end, /* int& */ on_or_off);
}
void Bmad::set_lords_status_stale(
    EleStruct &ele,
    int stat_group,
    std::optional<bool> control_bookkeeping,
    std::optional<int> flag
) {
  bool control_bookkeeping_lvalue;
  auto *_control_bookkeeping{&control_bookkeeping_lvalue};
  if (control_bookkeeping.has_value()) {
    control_bookkeeping_lvalue = control_bookkeeping.value();
  } else {
    _control_bookkeeping = nullptr;
  }
  int flag_lvalue;
  auto *_flag{&flag_lvalue};
  if (flag.has_value()) {
    flag_lvalue = flag.value();
  } else {
    _flag = nullptr;
  }
  fortran_set_lords_status_stale(/* void* */ ele.get_fortran_ptr(),
                                 /* int& */ stat_group,
                                 /* bool* */ _control_bookkeeping,
                                 /* int* */ _flag);
}
void Bmad::set_on_off(
    int key,
    LatStruct &lat,
    int switch_,
    std::optional<CoordStructArray1D> orb,
    std::optional<bool> use_ref_orb,
    std::optional<int> ix_branch,
    optional_ref<RealAlloc1D> saved_values,
    std::optional<std::string> attribute,
    std::optional<int> set_val
) {
  // orb: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb_desc;
  _orb_desc.rank = 1;
  if (orb) {
    _orb_desc.data_ptr = orb->data();
    _orb_desc.dims[0] = orb->size();
  } else {
    _orb_desc.data_ptr = nullptr;
    _orb_desc.dims[0] = 0;
  }
  _orb_desc.strides[0] = 1;
  bool use_ref_orb_lvalue;
  auto *_use_ref_orb{&use_ref_orb_lvalue};
  if (use_ref_orb.has_value()) {
    use_ref_orb_lvalue = use_ref_orb.value();
  } else {
    _use_ref_orb = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  // intent=inout allocatable general array
  auto *_saved_values =
      saved_values.has_value() ? saved_values->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_attribute = attribute.has_value() ? attribute->c_str() : nullptr;
  int set_val_lvalue;
  auto *_set_val{&set_val_lvalue};
  if (set_val.has_value()) {
    set_val_lvalue = set_val.value();
  } else {
    _set_val = nullptr;
  }
  fortran_set_on_off(
      /* int& */ key,
      /* void* */ lat.get_fortran_ptr(),
      /* int& */ switch_,
      /* Bmad::array_descriptor_t& */ _orb_desc,
      /* bool* */ _use_ref_orb,
      /* int* */ _ix_branch,
      /* void* */ _saved_values,
      /* const char* */ _attribute,
      /* int* */ _set_val
  );
}
void Bmad::set_orbit_to_zero(
    CoordStructArray1D orbit,
    int n1,
    int n2,
    std::optional<int> ix_noset
) {
  // orbit: CoordStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  _orbit_desc.data_ptr = orbit.data();
  _orbit_desc.dims[0] = orbit.size();
  _orbit_desc.strides[0] = 1;
  int ix_noset_lvalue;
  auto *_ix_noset{&ix_noset_lvalue};
  if (ix_noset.has_value()) {
    ix_noset_lvalue = ix_noset.value();
  } else {
    _ix_noset = nullptr;
  }
  fortran_set_orbit_to_zero(
      /* Bmad::array_descriptor_t& */ _orbit_desc,
      /* int& */ n1,
      /* int& */ n2,
      /* int* */ _ix_noset
  );
}
void Bmad::set_ptc(
    std::optional<double> e_tot,
    std::optional<int> particle,
    std::optional<int> taylor_order,
    std::optional<int> integ_order,
    std::optional<int> n_step,
    std::optional<bool> no_cavity,
    std::optional<bool> force_init
) {
  double e_tot_lvalue;
  auto *_e_tot{&e_tot_lvalue};
  if (e_tot.has_value()) {
    e_tot_lvalue = e_tot.value();
  } else {
    _e_tot = nullptr;
  }
  int particle_lvalue;
  auto *_particle{&particle_lvalue};
  if (particle.has_value()) {
    particle_lvalue = particle.value();
  } else {
    _particle = nullptr;
  }
  int taylor_order_lvalue;
  auto *_taylor_order{&taylor_order_lvalue};
  if (taylor_order.has_value()) {
    taylor_order_lvalue = taylor_order.value();
  } else {
    _taylor_order = nullptr;
  }
  int integ_order_lvalue;
  auto *_integ_order{&integ_order_lvalue};
  if (integ_order.has_value()) {
    integ_order_lvalue = integ_order.value();
  } else {
    _integ_order = nullptr;
  }
  int n_step_lvalue;
  auto *_n_step{&n_step_lvalue};
  if (n_step.has_value()) {
    n_step_lvalue = n_step.value();
  } else {
    _n_step = nullptr;
  }
  bool no_cavity_lvalue;
  auto *_no_cavity{&no_cavity_lvalue};
  if (no_cavity.has_value()) {
    no_cavity_lvalue = no_cavity.value();
  } else {
    _no_cavity = nullptr;
  }
  bool force_init_lvalue;
  auto *_force_init{&force_init_lvalue};
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
      /* bool* */ _force_init
  );
}
bool Bmad::set_ptc_base_state(std::string component, bool set_val) {
  auto _component = component.c_str();
  bool _old_val{};
  fortran_set_ptc_base_state(
      /* const char* */ _component,
      /* bool& */ set_val,
      /* bool& */ _old_val
  );
  return _old_val;
}
void Bmad::set_ptc_com_pointers() { fortran_set_ptc_com_pointers(); }
void Bmad::set_ptc_quiet(int channel, bool set, int &old_val) {
  fortran_set_ptc_quiet(/* int& */ channel, /* bool& */ set, /* int& */ old_val);
}
void Bmad::set_ptc_verbose(bool on) { fortran_set_ptc_verbose(/* bool& */ on); }
void Bmad::set_pwd_ele(LatStruct &lat, NormalModesStruct &mode0, double inductance) {
  fortran_set_pwd_ele(/* void* */ lat.get_fortran_ptr(),
                      /* void* */ mode0.get_fortran_ptr(),
                      /* double& */ inductance);
}
BookkeepingStateStruct Bmad::set_status_flags(int stat) {
  BookkeepingStateStruct _bookkeeping_state;
  fortran_set_status_flags(/* void* */ _bookkeeping_state.get_fortran_ptr(), /* int& */ stat);
  return std::move(_bookkeeping_state);
}
bool Bmad::set_tune(
    double phi_a_set,
    double phi_b_set,
    FArray1D<Real> &dk1,
    ElePointerStructArray1D eles,
    BranchStruct &branch,
    CoordStructAlloc1D orb,
    std::optional<bool> print_err
) {
  // dk1: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _dk1_desc;
  _dk1_desc.rank = 1;
  _dk1_desc.data_ptr = dk1.data();
  _dk1_desc.dims[0] = dk1.size();
  // eles: ElePointerStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _eles_desc;
  _eles_desc.rank = 1;
  _eles_desc.data_ptr = eles.data();
  _eles_desc.dims[0] = eles.size();
  _eles_desc.strides[0] = 1;
  // intent=inout allocatable type array
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool _ok{};
  fortran_set_tune(
      /* double& */ phi_a_set,
      /* double& */ phi_b_set,
      /* Bmad::array_descriptor_t& */ _dk1_desc,
      /* Bmad::array_descriptor_t& */ _eles_desc,
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ orb.get_fortran_ptr(),
      /* bool* */ _print_err,
      /* bool& */ _ok
  );
  return _ok;
}
void Bmad::set_twiss(
    BranchStruct &branch,
    EleStruct &twiss_ele,
    int ix_ele,
    bool match_deta_ds,
    bool err_flag,
    std::optional<bool> print_err
) {
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_twiss(/* void* */ branch.get_fortran_ptr(),
                    /* void* */ twiss_ele.get_fortran_ptr(),
                    /* int& */ ix_ele,
                    /* bool& */ match_deta_ds,
                    /* bool& */ err_flag,
                    /* bool* */ _print_err);
}
bool Bmad::set_z_tune(BranchStruct &branch, double z_tune, std::optional<bool> print_err) {
  bool _ok{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_z_tune(/* void* */ branch.get_fortran_ptr(),
                     /* double& */ z_tune,
                     /* bool& */ _ok,
                     /* bool* */ _print_err);
  return _ok;
}
void Bmad::settable_dep_var_bookkeeping(EleStruct &ele) {
  fortran_settable_dep_var_bookkeeping(/* void* */ ele.get_fortran_ptr());
}
void Bmad::setup_high_energy_space_charge_calc(
    bool calc_on,
    BranchStruct &branch,
    double n_part,
    NormalModesStruct &mode,
    optional_ref<BeamInitStruct> beam_init,
    std::optional<CoordStructArray1D> closed_orb
) {
  auto *_beam_init =
      beam_init.has_value() ? beam_init->get().get_fortran_ptr() : nullptr; // input, optional
  // closed_orb: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _closed_orb_desc;
  _closed_orb_desc.rank = 1;
  if (closed_orb) {
    _closed_orb_desc.data_ptr = closed_orb->data();
    _closed_orb_desc.dims[0] = closed_orb->size();
  } else {
    _closed_orb_desc.data_ptr = nullptr;
    _closed_orb_desc.dims[0] = 0;
  }
  _closed_orb_desc.strides[0] = 1;
  fortran_setup_high_energy_space_charge_calc(
      /* bool& */ calc_on,
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ n_part,
      /* void* */ mode.get_fortran_ptr(),
      /* void* */ _beam_init,
      /* Bmad::array_descriptor_t& */ _closed_orb_desc
  );
}
FixedArray2D<Real, 6, 6>
Bmad::sigma_mat_ptc_to_bmad(FixedArray2D<Real, 6, 6> sigma_mat_ptc, double beta0) {
  // sigma_mat_ptc: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_ptc_desc;
  _sigma_mat_ptc_desc.rank = 2;
  double _sigma_mat_ptc_vec[6 * 6];
  _sigma_mat_ptc_desc.data_ptr = _sigma_mat_ptc_vec;
  _sigma_mat_ptc_desc.dims[0] = 6;
  _sigma_mat_ptc_desc.dims[1] = 6;
  matrix_to_vec(sigma_mat_ptc, _sigma_mat_ptc_vec);
  // sigma_mat_bmad: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _sigma_mat_bmad_desc;
  _sigma_mat_bmad_desc.rank = 2;
  FixedArray2D<Real, 6, 6> sigma_mat_bmad;
  double _sigma_mat_bmad_vec[6 * 6];
  _sigma_mat_bmad_desc.data_ptr = _sigma_mat_bmad_vec;
  _sigma_mat_bmad_desc.dims[0] = 6;
  _sigma_mat_bmad_desc.dims[1] = 6;
  fortran_sigma_mat_ptc_to_bmad(
      /* Bmad::array_descriptor_t& */ _sigma_mat_ptc_desc,
      /* double& */ beta0,
      /* Bmad::array_descriptor_t& */ _sigma_mat_bmad_desc
  );
  vec_to_matrix(_sigma_mat_bmad_vec, sigma_mat_bmad);
  return sigma_mat_bmad;
}
bool Bmad::significant_difference(
    double value1,
    double value2,
    std::optional<double> abs_tol,
    std::optional<double> rel_tol
) {
  double abs_tol_lvalue;
  auto *_abs_tol{&abs_tol_lvalue};
  if (abs_tol.has_value()) {
    abs_tol_lvalue = abs_tol.value();
  } else {
    _abs_tol = nullptr;
  }
  double rel_tol_lvalue;
  auto *_rel_tol{&rel_tol_lvalue};
  if (rel_tol.has_value()) {
    rel_tol_lvalue = rel_tol.value();
  } else {
    _rel_tol = nullptr;
  }
  bool _is_different{};
  fortran_significant_difference(
      /* double& */ value1,
      /* double& */ value2,
      /* double* */ _abs_tol,
      /* double* */ _rel_tol,
      /* bool& */ _is_different
  );
  return _is_different;
}
void Bmad::skip_ele_blender(EleStruct &ele, bool skip) {
  fortran_skip_ele_blender(/* void* */ ele.get_fortran_ptr(), /* bool& */ skip);
}
bool Bmad::slice_lattice(
    LatStruct &lat,
    std::string ele_list,
    std::optional<bool> do_bookkeeping,
    std::optional<bool> set_phase_zero
) {
  auto _ele_list = ele_list.c_str();
  bool _error{};
  bool do_bookkeeping_lvalue;
  auto *_do_bookkeeping{&do_bookkeeping_lvalue};
  if (do_bookkeeping.has_value()) {
    do_bookkeeping_lvalue = do_bookkeeping.value();
  } else {
    _do_bookkeeping = nullptr;
  }
  bool set_phase_zero_lvalue;
  auto *_set_phase_zero{&set_phase_zero_lvalue};
  if (set_phase_zero.has_value()) {
    set_phase_zero_lvalue = set_phase_zero.value();
  } else {
    _set_phase_zero = nullptr;
  }
  fortran_slice_lattice(/* void* */ lat.get_fortran_ptr(),
                        /* const char* */ _ele_list,
                        /* bool& */ _error,
                        /* bool* */ _do_bookkeeping,
                        /* bool* */ _set_phase_zero);
  return _error;
}
void Bmad::soft_quadrupole_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_soft_quadrupole_edge_kick(/* void* */ ele.get_fortran_ptr(),
                                    /* void* */ param.get_fortran_ptr(),
                                    /* int& */ particle_at,
                                    /* void* */ orbit.get_fortran_ptr(),
                                    /* Bmad::array_descriptor_t& */ _mat6_desc,
                                    /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::sol_quad_mat6_calc(
    double ks_in,
    double k1_in,
    double tilt,
    double length,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
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
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* bool* */ _make_matrix
  );
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
double Bmad::solve_psi_adaptive(double t0, double t1, double p0, FixedArray1D<Real, 8> args) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  double _p1{};
  fortran_solve_psi_adaptive(
      /* double& */ t0,
      /* double& */ t1,
      /* double& */ p0,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* double& */ _p1
  );
  return _p1;
}
void Bmad::solve_psi_fixed_steps(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args,
    FArray1D<Real> &t,
    FArray1D<Real> &p
) {
  // args: in NOT (CppWrapperGeneralArgumentArray) (['1:8'])
  Bmad::array_descriptor_t _args_desc;
  _args_desc.rank = 1;
  _args_desc.data_ptr = args.data();
  _args_desc.dims[0] = args.size();
  // t: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _t_desc;
  _t_desc.rank = 1;
  _t_desc.data_ptr = t.data();
  _t_desc.dims[0] = t.size();
  // p: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _p_desc;
  _p_desc.rank = 1;
  _p_desc.data_ptr = p.data();
  _p_desc.dims[0] = p.size();
  fortran_solve_psi_fixed_steps(
      /* double& */ t0,
      /* double& */ t1,
      /* double& */ p0,
      /* Bmad::array_descriptor_t& */ _args_desc,
      /* Bmad::array_descriptor_t& */ _t_desc,
      /* Bmad::array_descriptor_t& */ _p_desc
  );
}
ComplexTaylorStruct Bmad::sort_complex_taylor_terms(ComplexTaylorStruct &complex_taylor_in) {
  ComplexTaylorStruct _complex_taylor_sorted;
  fortran_sort_complex_taylor_terms(/* void* */ complex_taylor_in.get_fortran_ptr(),
                                    /* void* */ _complex_taylor_sorted.get_fortran_ptr());
  return std::move(_complex_taylor_sorted);
}
Bmad::SpinDnDpzFromMat8 Bmad::spin_dn_dpz_from_mat8(
    FixedArray2D<Real, 8, 8> mat_1turn,
    std::optional<FixedArray2D<Real, 3, 3>> dn_dpz_partial
) {
  // mat_1turn: in NOT (CppWrapperGeneralArgumentArray) (['8', '8'])
  Bmad::array_descriptor_t _mat_1turn_desc;
  _mat_1turn_desc.rank = 2;
  double _mat_1turn_vec[8 * 8];
  _mat_1turn_desc.data_ptr = _mat_1turn_vec;
  _mat_1turn_desc.dims[0] = 8;
  _mat_1turn_desc.dims[1] = 8;
  matrix_to_vec(mat_1turn, _mat_1turn_vec);
  // dn_dpz_partial: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _dn_dpz_partial_desc;
  _dn_dpz_partial_desc.rank = 2;
  double _dn_dpz_partial_vec[3 * 3];
  _dn_dpz_partial_desc.data_ptr = _dn_dpz_partial_vec;
  if (dn_dpz_partial.has_value()) {
    matrix_to_vec(dn_dpz_partial.value(), _dn_dpz_partial_vec);
    _dn_dpz_partial_desc.dims[0] = 3;
    _dn_dpz_partial_desc.dims[1] = 3;
  } else {
    _dn_dpz_partial_desc.data_ptr = nullptr;
  }
  bool _error{};
  // dn_dpz: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _dn_dpz_desc;
  _dn_dpz_desc.rank = 1;
  FixedArray1D<Real, 3> _dn_dpz;
  _dn_dpz_desc.data_ptr = _dn_dpz.data();
  _dn_dpz_desc.dims[0] = _dn_dpz.size();
  fortran_spin_dn_dpz_from_mat8(
      /* Bmad::array_descriptor_t& */ _mat_1turn_desc,
      /* Bmad::array_descriptor_t& */ _dn_dpz_partial_desc,
      /* bool& */ _error,
      /* Bmad::array_descriptor_t& */ _dn_dpz_desc
  );
  return SpinDnDpzFromMat8{_error, _dn_dpz};
}
Bmad::SpinDnDpzFromQmap Bmad::spin_dn_dpz_from_qmap(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> q_map,
    FixedArray2D<Real, 3, 3> dn_dpz_partial,
    FixedArray2D<Real, 3, 3> dn_dpz_partial2,
    std::optional<FixedArray1D<Real, 3>> n0
) {
  // orb_mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _orb_mat_desc;
  _orb_mat_desc.rank = 2;
  double _orb_mat_vec[6 * 6];
  _orb_mat_desc.data_ptr = _orb_mat_vec;
  _orb_mat_desc.dims[0] = 6;
  _orb_mat_desc.dims[1] = 6;
  matrix_to_vec(orb_mat, _orb_mat_vec);
  // q_map: in NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _q_map_desc;
  _q_map_desc.rank = 2;
  double _q_map_vec[4 * 7];
  _q_map_desc.data_ptr = _q_map_vec;
  _q_map_desc.dims[0] = 4;
  _q_map_desc.dims[1] = 7;
  matrix_to_vec(q_map, _q_map_vec);
  // dn_dpz_partial: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _dn_dpz_partial_desc;
  _dn_dpz_partial_desc.rank = 2;
  double _dn_dpz_partial_vec[3 * 3];
  _dn_dpz_partial_desc.data_ptr = _dn_dpz_partial_vec;
  _dn_dpz_partial_desc.dims[0] = 3;
  _dn_dpz_partial_desc.dims[1] = 3;
  matrix_to_vec(dn_dpz_partial, _dn_dpz_partial_vec);
  // dn_dpz_partial2: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _dn_dpz_partial2_desc;
  _dn_dpz_partial2_desc.rank = 2;
  double _dn_dpz_partial2_vec[3 * 3];
  _dn_dpz_partial2_desc.data_ptr = _dn_dpz_partial2_vec;
  _dn_dpz_partial2_desc.dims[0] = 3;
  _dn_dpz_partial2_desc.dims[1] = 3;
  matrix_to_vec(dn_dpz_partial2, _dn_dpz_partial2_vec);
  bool _error{};
  // n0: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _n0_desc;
  _n0_desc.rank = 1;
  if (n0.has_value()) {
    _n0_desc.data_ptr = n0->data();
    _n0_desc.dims[0] = n0->size();
  } else {
    _n0_desc.data_ptr = nullptr;
    _n0_desc.dims[0] = 0;
  }
  // dn_dpz: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _dn_dpz_desc;
  _dn_dpz_desc.rank = 1;
  FixedArray1D<Real, 3> _dn_dpz;
  _dn_dpz_desc.data_ptr = _dn_dpz.data();
  _dn_dpz_desc.dims[0] = _dn_dpz.size();
  fortran_spin_dn_dpz_from_qmap(
      /* Bmad::array_descriptor_t& */ _orb_mat_desc,
      /* Bmad::array_descriptor_t& */ _q_map_desc,
      /* Bmad::array_descriptor_t& */ _dn_dpz_partial_desc,
      /* Bmad::array_descriptor_t& */ _dn_dpz_partial2_desc,
      /* bool& */ _error,
      /* Bmad::array_descriptor_t& */ _n0_desc,
      /* Bmad::array_descriptor_t& */ _dn_dpz_desc
  );
  return SpinDnDpzFromQmap{_error, _dn_dpz};
}
void Bmad::spin_map1_normalize(FixedArray2D<Real, 4, 7> spin1) {
  // spin1: inout NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _spin1_desc;
  _spin1_desc.rank = 2;
  double _spin1_vec[4 * 7];
  _spin1_desc.data_ptr = _spin1_vec;
  _spin1_desc.dims[0] = 4;
  _spin1_desc.dims[1] = 7;
  matrix_to_vec(spin1, _spin1_vec);
  fortran_spin_map1_normalize(/* Bmad::array_descriptor_t& */ _spin1_desc);
  vec_to_matrix(_spin1_vec, spin1);
}
Bmad::SpinMat8ResonanceStrengths Bmad::spin_mat8_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 6, 6> mat8
) {
  // orb_evec: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orb_evec_desc;
  _orb_evec_desc.rank = 1;
  _orb_evec_desc.data_ptr = orb_evec.data();
  _orb_evec_desc.dims[0] = orb_evec.size();
  // mat8: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat8_desc;
  _mat8_desc.rank = 2;
  double _mat8_vec[6 * 6];
  _mat8_desc.data_ptr = _mat8_vec;
  _mat8_desc.dims[0] = 6;
  _mat8_desc.dims[1] = 6;
  matrix_to_vec(mat8, _mat8_vec);
  double _xi_sum{};
  double _xi_diff{};
  fortran_spin_mat8_resonance_strengths(
      /* Bmad::array_descriptor_t& */ _orb_evec_desc,
      /* Bmad::array_descriptor_t& */ _mat8_desc,
      /* double& */ _xi_sum,
      /* double& */ _xi_diff
  );
  return SpinMat8ResonanceStrengths{_xi_sum, _xi_diff};
}
Bmad::SpinMatToEigen
Bmad::spin_mat_to_eigen(FixedArray2D<Real, 6, 6> orb_mat, FixedArray2D<Real, 4, 7> spin_map) {
  // orb_mat: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _orb_mat_desc;
  _orb_mat_desc.rank = 2;
  double _orb_mat_vec[6 * 6];
  _orb_mat_desc.data_ptr = _orb_mat_vec;
  _orb_mat_desc.dims[0] = 6;
  _orb_mat_desc.dims[1] = 6;
  matrix_to_vec(orb_mat, _orb_mat_vec);
  // spin_map: in NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _spin_map_desc;
  _spin_map_desc.rank = 2;
  double _spin_map_vec[4 * 7];
  _spin_map_desc.data_ptr = _spin_map_vec;
  _spin_map_desc.dims[0] = 4;
  _spin_map_desc.dims[1] = 7;
  matrix_to_vec(spin_map, _spin_map_vec);
  // orb_eval: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orb_eval_desc;
  _orb_eval_desc.rank = 1;
  FixedArray1D<Complex, 6> _orb_eval;
  _orb_eval_desc.data_ptr = _orb_eval.data();
  _orb_eval_desc.dims[0] = _orb_eval.size();
  // orb_evec: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _orb_evec_desc;
  _orb_evec_desc.rank = 2;
  FixedArray2D<Complex, 6, 6> orb_evec;
  std::complex<double> _orb_evec_vec[6 * 6];
  _orb_evec_desc.data_ptr = _orb_evec_vec;
  _orb_evec_desc.dims[0] = 6;
  _orb_evec_desc.dims[1] = 6;
  // n0: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _n0_desc;
  _n0_desc.rank = 1;
  FixedArray1D<Real, 3> _n0;
  _n0_desc.data_ptr = _n0.data();
  _n0_desc.dims[0] = _n0.size();
  // spin_evec: out NOT (CppWrapperGeneralArgumentArray) (['6', '3'])
  Bmad::array_descriptor_t _spin_evec_desc;
  _spin_evec_desc.rank = 2;
  FixedArray2D<Complex, 6, 3> spin_evec;
  std::complex<double> _spin_evec_vec[6 * 3];
  _spin_evec_desc.data_ptr = _spin_evec_vec;
  _spin_evec_desc.dims[0] = 6;
  _spin_evec_desc.dims[1] = 3;
  bool _error{};
  fortran_spin_mat_to_eigen(
      /* Bmad::array_descriptor_t& */ _orb_mat_desc,
      /* Bmad::array_descriptor_t& */ _spin_map_desc,
      /* Bmad::array_descriptor_t& */ _orb_eval_desc,
      /* Bmad::array_descriptor_t& */ _orb_evec_desc,
      /* Bmad::array_descriptor_t& */ _n0_desc,
      /* Bmad::array_descriptor_t& */ _spin_evec_desc,
      /* bool& */ _error
  );
  vec_to_matrix(_orb_evec_vec, orb_evec);
  vec_to_matrix(_spin_evec_vec, spin_evec);
  return SpinMatToEigen{_orb_eval, orb_evec, _n0, spin_evec, _error};
}
void Bmad::spin_omega(
    EmFieldStruct &field,
    CoordStruct &orbit,
    int sign_z_vel,
    FixedArray1D<Real, 3> omega,
    std::optional<bool> phase_space_coords
) {
  bool phase_space_coords_lvalue;
  auto *_phase_space_coords{&phase_space_coords_lvalue};
  if (phase_space_coords.has_value()) {
    phase_space_coords_lvalue = phase_space_coords.value();
  } else {
    _phase_space_coords = nullptr;
  }
  // omega: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _omega_desc;
  _omega_desc.rank = 1;
  _omega_desc.data_ptr = omega.data();
  _omega_desc.dims[0] = omega.size();
  fortran_spin_omega(/* void* */ field.get_fortran_ptr(),
                     /* void* */ orbit.get_fortran_ptr(),
                     /* int& */ sign_z_vel,
                     /* bool* */ _phase_space_coords,
                     /* Bmad::array_descriptor_t& */ _omega_desc);
}
Bmad::SpinQuatResonanceStrengths Bmad::spin_quat_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 4, 7> spin_q
) {
  // orb_evec: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orb_evec_desc;
  _orb_evec_desc.rank = 1;
  _orb_evec_desc.data_ptr = orb_evec.data();
  _orb_evec_desc.dims[0] = orb_evec.size();
  // spin_q: in NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _spin_q_desc;
  _spin_q_desc.rank = 2;
  double _spin_q_vec[4 * 7];
  _spin_q_desc.data_ptr = _spin_q_vec;
  _spin_q_desc.dims[0] = 4;
  _spin_q_desc.dims[1] = 7;
  matrix_to_vec(spin_q, _spin_q_vec);
  double _xi_sum{};
  double _xi_diff{};
  fortran_spin_quat_resonance_strengths(
      /* Bmad::array_descriptor_t& */ _orb_evec_desc,
      /* Bmad::array_descriptor_t& */ _spin_q_desc,
      /* double& */ _xi_sum,
      /* double& */ _xi_diff
  );
  return SpinQuatResonanceStrengths{_xi_sum, _xi_diff};
}
FixedArray2D<Real, 4, 7> Bmad::spin_taylor_to_linear(
    TaylorStructArray1D spin_taylor,
    bool normalize,
    FixedArray1D<Real, 6> dref_orb,
    bool is_on
) {
  // spin_taylor: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_taylor_desc;
  _spin_taylor_desc.rank = 1;
  _spin_taylor_desc.data_ptr = spin_taylor.data();
  _spin_taylor_desc.dims[0] = spin_taylor.size();
  _spin_taylor_desc.strides[0] = 1;
  // dref_orb: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _dref_orb_desc;
  _dref_orb_desc.rank = 1;
  _dref_orb_desc.data_ptr = dref_orb.data();
  _dref_orb_desc.dims[0] = dref_orb.size();
  // spin_map1: out NOT (CppWrapperGeneralArgumentArray) (['0:3', '0:6'])
  Bmad::array_descriptor_t _spin_map1_desc;
  _spin_map1_desc.rank = 2;
  FixedArray2D<Real, 4, 7> spin_map1;
  double _spin_map1_vec[4 * 7];
  _spin_map1_desc.data_ptr = _spin_map1_vec;
  _spin_map1_desc.dims[0] = 4;
  _spin_map1_desc.dims[1] = 7;
  fortran_spin_taylor_to_linear(
      /* Bmad::array_descriptor_t& */ _spin_taylor_desc,
      /* bool& */ normalize,
      /* Bmad::array_descriptor_t& */ _dref_orb_desc,
      /* bool& */ is_on,
      /* Bmad::array_descriptor_t& */ _spin_map1_desc
  );
  vec_to_matrix(_spin_map1_vec, spin_map1);
  return spin_map1;
}
SpinPolarStruct Bmad::spinor_to_polar(FixedArray1D<Complex, 2> spinor) {
  // spinor: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _spinor_desc;
  _spinor_desc.rank = 1;
  _spinor_desc.data_ptr = spinor.data();
  _spinor_desc.dims[0] = spinor.size();
  SpinPolarStruct _polar;
  fortran_spinor_to_polar(
      /* Bmad::array_descriptor_t& */ _spinor_desc,
      /* void* */ _polar.get_fortran_ptr()
  );
  return std::move(_polar);
}
FixedArray1D<Real, 3> Bmad::spinor_to_vec(FixedArray1D<Complex, 2> spinor) {
  // spinor: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _spinor_desc;
  _spinor_desc.rank = 1;
  _spinor_desc.data_ptr = spinor.data();
  _spinor_desc.dims[0] = spinor.size();
  // vec: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  FixedArray1D<Real, 3> _vec;
  _vec_desc.data_ptr = _vec.data();
  _vec_desc.dims[0] = _vec.size();
  fortran_spinor_to_vec(
      /* Bmad::array_descriptor_t& */ _spinor_desc,
      /* Bmad::array_descriptor_t& */ _vec_desc
  );
  return _vec;
}
void Bmad::spline_fit_orbit(
    CoordStruct &start_orb,
    CoordStruct &end_orb,
    FixedArray1D<Real, 4> spline_x,
    FixedArray1D<Real, 4> spline_y
) {
  // spline_x: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _spline_x_desc;
  _spline_x_desc.rank = 1;
  _spline_x_desc.data_ptr = spline_x.data();
  _spline_x_desc.dims[0] = spline_x.size();
  // spline_y: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _spline_y_desc;
  _spline_y_desc.rank = 1;
  _spline_y_desc.data_ptr = spline_y.data();
  _spline_y_desc.dims[0] = spline_y.size();
  fortran_spline_fit_orbit(/* void* */ start_orb.get_fortran_ptr(),
                           /* void* */ end_orb.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _spline_x_desc,
                           /* Bmad::array_descriptor_t& */ _spline_y_desc);
}
CharacterAlloc1D Bmad::split_expression_string(
    std::string expr,
    int width,
    int indent,
    std::optional<std::string> break_str
) {
  auto _expr = expr.c_str();
  // intent=out character array container
  auto lines{CharacterAlloc1D()};
  const char *_break_str = break_str.has_value() ? break_str->c_str() : nullptr;
  fortran_split_expression_string(
      /* const char* */ _expr,
      /* int& */ width,
      /* int& */ indent,
      /* void* */ lines.get_fortran_ptr(),
      /* const char* */ _break_str
  );
  return std::move(lines);
}
Bmad::SplitLat Bmad::split_lat(
    LatStruct &lat,
    double s_split,
    int ix_branch,
    std::optional<bool> add_suffix,
    std::optional<bool> check_sanity,
    std::optional<bool> save_null_drift,
    std::optional<bool> choose_max,
    std::optional<int> ix_insert
) {
  int _ix_split{};
  bool _split_done{};
  bool add_suffix_lvalue;
  auto *_add_suffix{&add_suffix_lvalue};
  if (add_suffix.has_value()) {
    add_suffix_lvalue = add_suffix.value();
  } else {
    _add_suffix = nullptr;
  }
  bool check_sanity_lvalue;
  auto *_check_sanity{&check_sanity_lvalue};
  if (check_sanity.has_value()) {
    check_sanity_lvalue = check_sanity.value();
  } else {
    _check_sanity = nullptr;
  }
  bool save_null_drift_lvalue;
  auto *_save_null_drift{&save_null_drift_lvalue};
  if (save_null_drift.has_value()) {
    save_null_drift_lvalue = save_null_drift.value();
  } else {
    _save_null_drift = nullptr;
  }
  bool _err_flag{};
  bool choose_max_lvalue;
  auto *_choose_max{&choose_max_lvalue};
  if (choose_max.has_value()) {
    choose_max_lvalue = choose_max.value();
  } else {
    _choose_max = nullptr;
  }
  int ix_insert_lvalue;
  auto *_ix_insert{&ix_insert_lvalue};
  if (ix_insert.has_value()) {
    ix_insert_lvalue = ix_insert.value();
  } else {
    _ix_insert = nullptr;
  }
  fortran_split_lat(/* void* */ lat.get_fortran_ptr(),
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
    EleStruct &ele,
    std::optional<FixedArray1D<Real, 6>> start_orbit
) {
  // start_orbit: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _start_orbit_desc;
  _start_orbit_desc.rank = 1;
  if (start_orbit.has_value()) {
    _start_orbit_desc.data_ptr = start_orbit->data();
    _start_orbit_desc.dims[0] = start_orbit->size();
  } else {
    _start_orbit_desc.data_ptr = nullptr;
    _start_orbit_desc.dims[0] = 0;
  }
  fortran_sprint_spin_taylor_map(/* void* */ ele.get_fortran_ptr(),
                                 /* Bmad::array_descriptor_t& */ _start_orbit_desc);
}
void Bmad::sr_longitudinal_wake_particle(EleStruct &ele, CoordStruct &orbit) {
  fortran_sr_longitudinal_wake_particle(/* void* */ ele.get_fortran_ptr(),
                                        /* void* */ orbit.get_fortran_ptr());
}
void Bmad::sr_transverse_wake_particle(EleStruct &ele, CoordStruct &orbit) {
  fortran_sr_transverse_wake_particle(/* void* */ ele.get_fortran_ptr(),
                                      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::sr_z_long_wake(EleStruct &ele, BunchStruct &bunch, double z_ave) {
  fortran_sr_z_long_wake(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ bunch.get_fortran_ptr(),
                         /* double& */ z_ave);
}
SummationRdtStruct Bmad::srdt_calc(
    LatStruct &lat,
    int order,
    std::optional<int> n_slices_gen_opt,
    std::optional<int> n_slices_sxt_opt,
    std::optional<SummationRdtStructAlloc1D> per_ele_out
) {
  SummationRdtStruct _srdt_sums;
  int n_slices_gen_opt_lvalue;
  auto *_n_slices_gen_opt{&n_slices_gen_opt_lvalue};
  if (n_slices_gen_opt.has_value()) {
    n_slices_gen_opt_lvalue = n_slices_gen_opt.value();
  } else {
    _n_slices_gen_opt = nullptr;
  }
  int n_slices_sxt_opt_lvalue;
  auto *_n_slices_sxt_opt{&n_slices_sxt_opt_lvalue};
  if (n_slices_sxt_opt.has_value()) {
    n_slices_sxt_opt_lvalue = n_slices_sxt_opt.value();
  } else {
    _n_slices_sxt_opt = nullptr;
  }
  // intent=inout allocatable type array
  auto *_per_ele_out =
      per_ele_out.has_value() ? per_ele_out->get_fortran_ptr() : nullptr; // input, optional
  fortran_srdt_calc(/* void* */ lat.get_fortran_ptr(),
                    /* void* */ _srdt_sums.get_fortran_ptr(),
                    /* int& */ order,
                    /* int* */ _n_slices_gen_opt,
                    /* int* */ _n_slices_sxt_opt,
                    /* void* */ _per_ele_out);
  return std::move(_srdt_sums);
}
RealAlloc1D Bmad::srdt_lsq_solution(
    LatStruct &lat,
    FArray1D<Int> &var_indexes,
    std::optional<int> n_slices_gen_opt,
    std::optional<int> n_slices_sxt_opt,
    std::optional<double> chrom_set_x_opt,
    std::optional<double> chrom_set_y_opt,
    std::optional<FixedArray1D<Real, 10>> weight_in
) {
  // var_indexes: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _var_indexes_desc;
  _var_indexes_desc.rank = 1;
  _var_indexes_desc.data_ptr = var_indexes.data();
  _var_indexes_desc.dims[0] = var_indexes.size();
  // intent=out allocatable general array
  auto ls_soln{RealAlloc1D()};
  int n_slices_gen_opt_lvalue;
  auto *_n_slices_gen_opt{&n_slices_gen_opt_lvalue};
  if (n_slices_gen_opt.has_value()) {
    n_slices_gen_opt_lvalue = n_slices_gen_opt.value();
  } else {
    _n_slices_gen_opt = nullptr;
  }
  int n_slices_sxt_opt_lvalue;
  auto *_n_slices_sxt_opt{&n_slices_sxt_opt_lvalue};
  if (n_slices_sxt_opt.has_value()) {
    n_slices_sxt_opt_lvalue = n_slices_sxt_opt.value();
  } else {
    _n_slices_sxt_opt = nullptr;
  }
  double chrom_set_x_opt_lvalue;
  auto *_chrom_set_x_opt{&chrom_set_x_opt_lvalue};
  if (chrom_set_x_opt.has_value()) {
    chrom_set_x_opt_lvalue = chrom_set_x_opt.value();
  } else {
    _chrom_set_x_opt = nullptr;
  }
  double chrom_set_y_opt_lvalue;
  auto *_chrom_set_y_opt{&chrom_set_y_opt_lvalue};
  if (chrom_set_y_opt.has_value()) {
    chrom_set_y_opt_lvalue = chrom_set_y_opt.value();
  } else {
    _chrom_set_y_opt = nullptr;
  }
  // weight_in: in NOT (CppWrapperGeneralArgumentArray) (['10'])
  Bmad::array_descriptor_t _weight_in_desc;
  _weight_in_desc.rank = 1;
  if (weight_in.has_value()) {
    _weight_in_desc.data_ptr = weight_in->data();
    _weight_in_desc.dims[0] = weight_in->size();
  } else {
    _weight_in_desc.data_ptr = nullptr;
    _weight_in_desc.dims[0] = 0;
  }
  fortran_srdt_lsq_solution(/* void* */ lat.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _var_indexes_desc,
                            /* void* */ ls_soln.get_fortran_ptr(),
                            /* int* */ _n_slices_gen_opt,
                            /* int* */ _n_slices_sxt_opt,
                            /* double* */ _chrom_set_x_opt,
                            /* double* */ _chrom_set_y_opt,
                            /* Bmad::array_descriptor_t& */ _weight_in_desc);
  return std::move(ls_soln);
}
bool Bmad::start_branch_at(LatStruct &lat, std::string ele_start, bool move_end_marker) {
  auto _ele_start = ele_start.c_str();
  bool _error{};
  fortran_start_branch_at(/* void* */ lat.get_fortran_ptr(),
                          /* const char* */ _ele_start,
                          /* bool& */ move_end_marker,
                          /* bool& */ _error);
  return _error;
}
int Bmad::stream_ele_end(int physical_end, int ele_orientation) {
  int _stream_end{};
  fortran_stream_ele_end(
      /* int& */ physical_end,
      /* int& */ ele_orientation,
      /* int& */ _stream_end
  );
  return _stream_end;
}
std::string Bmad::string_attrib(std::string attrib_name, EleStruct &ele) {
  auto _attrib_name = attrib_name.c_str();
  char _attrib_value[4096];
  fortran_string_attrib(
      /* const char* */ _attrib_name,
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _attrib_value
  );
  return _attrib_value;
}
Bmad::StrongBeamSigmaCalc Bmad::strong_beam_sigma_calc(EleStruct &ele, double s_pos) {
  // sigma: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _sigma_desc;
  _sigma_desc.rank = 1;
  FixedArray1D<Real, 2> _sigma;
  _sigma_desc.data_ptr = _sigma.data();
  _sigma_desc.dims[0] = _sigma.size();
  double _bbi_const{};
  // dsigma_ds: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _dsigma_ds_desc;
  _dsigma_ds_desc.rank = 1;
  FixedArray1D<Real, 2> _dsigma_ds;
  _dsigma_ds_desc.data_ptr = _dsigma_ds.data();
  _dsigma_ds_desc.dims[0] = _dsigma_ds.size();
  fortran_strong_beam_sigma_calc(/* void* */ ele.get_fortran_ptr(),
                                 /* double& */ s_pos,
                                 /* Bmad::array_descriptor_t& */ _sigma_desc,
                                 /* double& */ _bbi_const,
                                 /* Bmad::array_descriptor_t& */ _dsigma_ds_desc);
  return StrongBeamSigmaCalc{_sigma, _bbi_const, _dsigma_ds};
}
double Bmad::strong_beam_strength(EleStruct &ele) {
  double _strength{};
  fortran_strong_beam_strength(/* void* */ ele.get_fortran_ptr(), /* double& */ _strength);
  return _strength;
}
Bmad::SurfaceGridDisplacement Bmad::surface_grid_displacement(
    EleStruct &ele,
    double x,
    double y,
    std::optional<bool> extend_grid
) {
  bool _err_flag{};
  double _z{};
  // dz_dxy: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _dz_dxy_desc;
  _dz_dxy_desc.rank = 1;
  FixedArray1D<Real, 2> _dz_dxy;
  _dz_dxy_desc.data_ptr = _dz_dxy.data();
  _dz_dxy_desc.dims[0] = _dz_dxy.size();
  bool extend_grid_lvalue;
  auto *_extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  fortran_surface_grid_displacement(/* void* */ ele.get_fortran_ptr(),
                                    /* double& */ x,
                                    /* double& */ y,
                                    /* bool& */ _err_flag,
                                    /* double& */ _z,
                                    /* Bmad::array_descriptor_t& */ _dz_dxy_desc,
                                    /* bool* */ _extend_grid);
  return SurfaceGridDisplacement{_err_flag, _z, _dz_dxy};
}
Bmad::SwitchAttribValueName
Bmad::switch_attrib_value_name(std::string attrib_name, double attrib_value, EleStruct &ele) {
  auto _attrib_name = attrib_name.c_str();
  bool _is_default{};
  // intent=out character array container
  auto name_list{CharacterAlloc1D()};
  char _attrib_val_name[4096];
  fortran_switch_attrib_value_name(
      /* const char* */ _attrib_name,
      /* double& */ attrib_value,
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _is_default,
      /* void* */ name_list.get_fortran_ptr(),
      /* const char* */ _attrib_val_name
  );
  return SwitchAttribValueName{_is_default, std::move(name_list), _attrib_val_name};
}
TrackStruct Bmad::symp_lie_bmad(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<bool> offset_ele
) {
  TrackStruct _track;
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  bool offset_ele_lvalue;
  auto *_offset_ele{&offset_ele_lvalue};
  if (offset_ele.has_value()) {
    offset_ele_lvalue = offset_ele.value();
  } else {
    _offset_ele = nullptr;
  }
  fortran_symp_lie_bmad(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr(),
                        /* void* */ _track.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _mat6_desc,
                        /* bool* */ _make_matrix,
                        /* bool* */ _offset_ele);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return std::move(_track);
}
Bmad::T6ToB123 Bmad::t6_to_b123(FixedArray2D<Real, 6, 6> t6, FixedArray1D<Real, 3> abz_tunes) {
  // t6: in NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _t6_desc;
  _t6_desc.rank = 2;
  double _t6_vec[6 * 6];
  _t6_desc.data_ptr = _t6_vec;
  _t6_desc.dims[0] = 6;
  _t6_desc.dims[1] = 6;
  matrix_to_vec(t6, _t6_vec);
  // abz_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _abz_tunes_desc;
  _abz_tunes_desc.rank = 1;
  _abz_tunes_desc.data_ptr = abz_tunes.data();
  _abz_tunes_desc.dims[0] = abz_tunes.size();
  // B1: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _B1_desc;
  _B1_desc.rank = 2;
  FixedArray2D<Real, 6, 6> B1;
  double _B1_vec[6 * 6];
  _B1_desc.data_ptr = _B1_vec;
  _B1_desc.dims[0] = 6;
  _B1_desc.dims[1] = 6;
  // B2: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _B2_desc;
  _B2_desc.rank = 2;
  FixedArray2D<Real, 6, 6> B2;
  double _B2_vec[6 * 6];
  _B2_desc.data_ptr = _B2_vec;
  _B2_desc.dims[0] = 6;
  _B2_desc.dims[1] = 6;
  // B3: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _B3_desc;
  _B3_desc.rank = 2;
  FixedArray2D<Real, 6, 6> B3;
  double _B3_vec[6 * 6];
  _B3_desc.data_ptr = _B3_vec;
  _B3_desc.dims[0] = 6;
  _B3_desc.dims[1] = 6;
  bool _err_flag{};
  fortran_t6_to_b123(
      /* Bmad::array_descriptor_t& */ _t6_desc,
      /* Bmad::array_descriptor_t& */ _abz_tunes_desc,
      /* Bmad::array_descriptor_t& */ _B1_desc,
      /* Bmad::array_descriptor_t& */ _B2_desc,
      /* Bmad::array_descriptor_t& */ _B3_desc,
      /* bool& */ _err_flag
  );
  vec_to_matrix(_B1_vec, B1);
  vec_to_matrix(_B2_vec, B2);
  vec_to_matrix(_B3_vec, B3);
  return T6ToB123{B1, B2, B3, _err_flag};
}
void Bmad::taper_mag_strengths(
    LatStruct &lat,
    optional_ref<LatStruct> ref_lat,
    std::optional<std::string> except,
    std::optional<bool> err_flag
) {
  auto *_ref_lat =
      ref_lat.has_value() ? ref_lat->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_except = except.has_value() ? except->c_str() : nullptr;
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  fortran_taper_mag_strengths(/* void* */ lat.get_fortran_ptr(),
                              /* void* */ _ref_lat,
                              /* const char* */ _except,
                              /* bool* */ _err_flag);
}
void Bmad::target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double &y_min,
    double &y_max,
    double &phi_min,
    double &phi_max,
    std::optional<bool> initial
) {
  // r_corner1: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r_corner1_desc;
  _r_corner1_desc.rank = 1;
  _r_corner1_desc.data_ptr = r_corner1.data();
  _r_corner1_desc.dims[0] = r_corner1.size();
  // r_corner2: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r_corner2_desc;
  _r_corner2_desc.rank = 1;
  _r_corner2_desc.data_ptr = r_corner2.data();
  _r_corner2_desc.dims[0] = r_corner2.size();
  bool initial_lvalue;
  auto *_initial{&initial_lvalue};
  if (initial.has_value()) {
    initial_lvalue = initial.value();
  } else {
    _initial = nullptr;
  }
  fortran_target_min_max_calc(
      /* Bmad::array_descriptor_t& */ _r_corner1_desc,
      /* Bmad::array_descriptor_t& */ _r_corner2_desc,
      /* double& */ y_min,
      /* double& */ y_max,
      /* double& */ phi_min,
      /* double& */ phi_max,
      /* bool* */ _initial
  );
}
Bmad::TargetRotMats Bmad::target_rot_mats(FixedArray1D<Real, 3> r_center) {
  // r_center: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r_center_desc;
  _r_center_desc.rank = 1;
  _r_center_desc.data_ptr = r_center.data();
  _r_center_desc.dims[0] = r_center.size();
  // w_to_target: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_to_target_desc;
  _w_to_target_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_to_target;
  double _w_to_target_vec[3 * 3];
  _w_to_target_desc.data_ptr = _w_to_target_vec;
  _w_to_target_desc.dims[0] = 3;
  _w_to_target_desc.dims[1] = 3;
  // w_to_ele: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_to_ele_desc;
  _w_to_ele_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_to_ele;
  double _w_to_ele_vec[3 * 3];
  _w_to_ele_desc.data_ptr = _w_to_ele_vec;
  _w_to_ele_desc.dims[0] = 3;
  _w_to_ele_desc.dims[1] = 3;
  fortran_target_rot_mats(
      /* Bmad::array_descriptor_t& */ _r_center_desc,
      /* Bmad::array_descriptor_t& */ _w_to_target_desc,
      /* Bmad::array_descriptor_t& */ _w_to_ele_desc
  );
  vec_to_matrix(_w_to_target_vec, w_to_target);
  vec_to_matrix(_w_to_ele_vec, w_to_ele);
  return TargetRotMats{w_to_target, w_to_ele};
}
void Bmad::taylor_equal_taylor(TaylorStruct &taylor1, TaylorStruct &taylor2) {
  fortran_taylor_equal_taylor(/* void* */ taylor1.get_fortran_ptr(),
                              /* void* */ taylor2.get_fortran_ptr());
}
bool Bmad::taylor_inverse(TaylorStructArray1D taylor_in, TaylorStructArray1D taylor_inv) {
  // taylor_in: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_in_desc;
  _taylor_in_desc.rank = 1;
  _taylor_in_desc.data_ptr = taylor_in.data();
  _taylor_in_desc.dims[0] = taylor_in.size();
  _taylor_in_desc.strides[0] = 1;
  // taylor_inv: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_inv_desc;
  _taylor_inv_desc.rank = 1;
  _taylor_inv_desc.data_ptr = taylor_inv.data();
  _taylor_inv_desc.dims[0] = taylor_inv.size();
  _taylor_inv_desc.strides[0] = 1;
  bool _err{};
  fortran_taylor_inverse(
      /* Bmad::array_descriptor_t& */ _taylor_in_desc,
      /* Bmad::array_descriptor_t& */ _taylor_inv_desc,
      /* bool& */ _err
  );
  return _err;
}
bool Bmad::taylor_propagate1(
    TaylorStructArray1D orb_taylor,
    EleStruct &ele,
    LatParamStruct &param,
    optional_ref<CoordStruct> ref_in,
    std::optional<TaylorStructArray1D> spin_taylor
) {
  // orb_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb_taylor_desc;
  _orb_taylor_desc.rank = 1;
  _orb_taylor_desc.data_ptr = orb_taylor.data();
  _orb_taylor_desc.dims[0] = orb_taylor.size();
  _orb_taylor_desc.strides[0] = 1;
  bool _err_flag{};
  auto *_ref_in = ref_in.has_value() ? ref_in->get().get_fortran_ptr() : nullptr; // input, optional
  // spin_taylor: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_taylor_desc;
  _spin_taylor_desc.rank = 1;
  if (spin_taylor) {
    _spin_taylor_desc.data_ptr = spin_taylor->data();
    _spin_taylor_desc.dims[0] = spin_taylor->size();
  } else {
    _spin_taylor_desc.data_ptr = nullptr;
    _spin_taylor_desc.dims[0] = 0;
  }
  _spin_taylor_desc.strides[0] = 1;
  fortran_taylor_propagate1(
      /* Bmad::array_descriptor_t& */ _orb_taylor_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* void* */ _ref_in,
      /* Bmad::array_descriptor_t& */ _spin_taylor_desc
  );
  return _err_flag;
}
MadMapStruct Bmad::taylor_to_mad_map(TaylorStructArray1D taylor, MadEnergyStruct &energy) {
  // taylor: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_desc;
  _taylor_desc.rank = 1;
  _taylor_desc.data_ptr = taylor.data();
  _taylor_desc.dims[0] = taylor.size();
  _taylor_desc.strides[0] = 1;
  MadMapStruct _map;
  fortran_taylor_to_mad_map(
      /* Bmad::array_descriptor_t& */ _taylor_desc,
      /* void* */ energy.get_fortran_ptr(),
      /* void* */ _map.get_fortran_ptr()
  );
  return std::move(_map);
}
void Bmad::taylors_equal_taylors(TaylorStructArray1D taylor1, TaylorStructArray1D taylor2) {
  // taylor1: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor1_desc;
  _taylor1_desc.rank = 1;
  _taylor1_desc.data_ptr = taylor1.data();
  _taylor1_desc.dims[0] = taylor1.size();
  _taylor1_desc.strides[0] = 1;
  // taylor2: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor2_desc;
  _taylor2_desc.rank = 1;
  _taylor2_desc.data_ptr = taylor2.data();
  _taylor2_desc.dims[0] = taylor2.size();
  _taylor2_desc.strides[0] = 1;
  fortran_taylors_equal_taylors(
      /* Bmad::array_descriptor_t& */ _taylor1_desc,
      /* Bmad::array_descriptor_t& */ _taylor2_desc
  );
}
void Bmad::tilt_coords(
    double tilt_val,
    FArray1D<Real> &coord,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // coord: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _coord_desc;
  _coord_desc.rank = 1;
  _coord_desc.data_ptr = coord.data();
  _coord_desc.dims[0] = coord.size();
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_tilt_coords(
      /* double& */ tilt_val,
      /* Bmad::array_descriptor_t& */ _coord_desc,
      /* Bmad::array_descriptor_t& */ _mat6_desc,
      /* bool* */ _make_matrix
  );
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::tilt_coords_photon(
    double tilt_val,
    FArray1D<Real> &coord,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat
) {
  // coord: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _coord_desc;
  _coord_desc.rank = 1;
  _coord_desc.data_ptr = coord.data();
  _coord_desc.dims[0] = coord.size();
  // w_mat: inout NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  if (w_mat.has_value()) {
    matrix_to_vec(w_mat.value(), _w_mat_vec);
    _w_mat_desc.dims[0] = 3;
    _w_mat_desc.dims[1] = 3;
  } else {
    _w_mat_desc.data_ptr = nullptr;
  }
  fortran_tilt_coords_photon(
      /* double& */ tilt_val,
      /* Bmad::array_descriptor_t& */ _coord_desc,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  if (w_mat.has_value())
    vec_to_matrix(_w_mat_vec, w_mat.value());
}
void Bmad::tilt_mat6(FixedArray2D<Real, 6, 6> mat6, double tilt) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  matrix_to_vec(mat6, _mat6_vec);
  fortran_tilt_mat6(/* Bmad::array_descriptor_t& */ _mat6_desc, /* double& */ tilt);
  vec_to_matrix(_mat6_vec, mat6);
}
Bmad::ToEtaReading
Bmad::to_eta_reading(FArray1D<Real> &eta_actual, EleStruct &ele, int axis, bool add_noise) {
  // eta_actual: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _eta_actual_desc;
  _eta_actual_desc.rank = 1;
  _eta_actual_desc.data_ptr = eta_actual.data();
  _eta_actual_desc.dims[0] = eta_actual.size();
  double _reading{};
  bool _err{};
  fortran_to_eta_reading(
      /* Bmad::array_descriptor_t& */ _eta_actual_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ axis,
      /* bool& */ add_noise,
      /* double& */ _reading,
      /* bool& */ _err
  );
  return ToEtaReading{_reading, _err};
}
Bmad::ToFieldmapCoords Bmad::to_fieldmap_coords(
    EleStruct &ele,
    CoordStruct &local_orb,
    double s_body,
    int ele_anchor_pt,
    FixedArray1D<Real, 3> r0,
    bool curved_ref_frame
) {
  // r0: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r0_desc;
  _r0_desc.rank = 1;
  _r0_desc.data_ptr = r0.data();
  _r0_desc.dims[0] = r0.size();
  double _x{};
  double _y{};
  double _z{};
  double _cos_ang{};
  double _sin_ang{};
  bool _err_flag{};
  fortran_to_fieldmap_coords(/* void* */ ele.get_fortran_ptr(),
                             /* void* */ local_orb.get_fortran_ptr(),
                             /* double& */ s_body,
                             /* int& */ ele_anchor_pt,
                             /* Bmad::array_descriptor_t& */ _r0_desc,
                             /* bool& */ curved_ref_frame,
                             /* double& */ _x,
                             /* double& */ _y,
                             /* double& */ _z,
                             /* double& */ _cos_ang,
                             /* double& */ _sin_ang,
                             /* bool& */ _err_flag);
  return ToFieldmapCoords{_x, _y, _z, _cos_ang, _sin_ang, _err_flag};
}
Bmad::ToOrbitReading
Bmad::to_orbit_reading(CoordStruct &orb, EleStruct &ele, int axis, bool add_noise) {
  double _reading{};
  bool _err{};
  fortran_to_orbit_reading(/* void* */ orb.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* int& */ axis,
                           /* bool& */ add_noise,
                           /* double& */ _reading,
                           /* bool& */ _err);
  return ToOrbitReading{_reading, _err};
}
Bmad::ToPhaseAndCouplingReading
Bmad::to_phase_and_coupling_reading(EleStruct &ele, bool add_noise) {
  BpmPhaseCouplingStruct _reading;
  bool _err{};
  fortran_to_phase_and_coupling_reading(/* void* */ ele.get_fortran_ptr(),
                                        /* bool& */ add_noise,
                                        /* void* */ _reading.get_fortran_ptr(),
                                        /* bool& */ _err);
  return ToPhaseAndCouplingReading{std::move(_reading), _err};
}
CoordStruct Bmad::to_photon_angle_coords(CoordStruct &orb_in, EleStruct &ele) {
  CoordStruct _orb_out;
  fortran_to_photon_angle_coords(/* void* */ orb_in.get_fortran_ptr(),
                                 /* void* */ ele.get_fortran_ptr(),
                                 /* void* */ _orb_out.get_fortran_ptr());
  return std::move(_orb_out);
}
CoordStruct Bmad::to_surface_coords(CoordStruct &lab_orbit, EleStruct &ele) {
  CoordStruct _surface_orbit;
  fortran_to_surface_coords(/* void* */ lab_orbit.get_fortran_ptr(),
                            /* void* */ ele.get_fortran_ptr(),
                            /* void* */ _surface_orbit.get_fortran_ptr());
  return std::move(_surface_orbit);
}
double Bmad::touschek_lifetime(NormalModesStruct &mode, LatStruct &lat) {
  double _Tl{};
  fortran_touschek_lifetime(/* void* */ mode.get_fortran_ptr(),
                            /* double& */ _Tl,
                            /* void* */ lat.get_fortran_ptr());
  return _Tl;
}
double Bmad::touschek_rate1(
    NormalModesStruct &mode,
    LatStruct &lat,
    std::optional<int> ix,
    std::optional<double> s
) {
  double _rate{};
  int ix_lvalue;
  auto *_ix{&ix_lvalue};
  if (ix.has_value()) {
    ix_lvalue = ix.value();
  } else {
    _ix = nullptr;
  }
  double s_lvalue;
  auto *_s{&s_lvalue};
  if (s.has_value()) {
    s_lvalue = s.value();
  } else {
    _s = nullptr;
  }
  fortran_touschek_rate1(/* void* */ mode.get_fortran_ptr(),
                         /* double& */ _rate,
                         /* void* */ lat.get_fortran_ptr(),
                         /* int* */ _ix,
                         /* double* */ _s);
  return _rate;
}
void Bmad::touschek_rate1_zap(
    NormalModesStruct &mode,
    double rate,
    LatStruct &lat,
    std::optional<int> ix,
    std::optional<double> s
) {
  int ix_lvalue;
  auto *_ix{&ix_lvalue};
  if (ix.has_value()) {
    ix_lvalue = ix.value();
  } else {
    _ix = nullptr;
  }
  double s_lvalue;
  auto *_s{&s_lvalue};
  if (s.has_value()) {
    s_lvalue = s.value();
  } else {
    _s = nullptr;
  }
  fortran_touschek_rate1_zap(/* void* */ mode.get_fortran_ptr(),
                             /* double& */ rate,
                             /* void* */ lat.get_fortran_ptr(),
                             /* int* */ _ix,
                             /* double* */ _s);
}
Bmad::Track1 Bmad::track1(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    optional_ref<TrackStruct> track,
    std::optional<bool> ignore_radiation,
    std::optional<bool> make_map1,
    std::optional<bool> init_to_edge
) {
  CoordStruct _end_orb;
  auto *_track = track.has_value() ? track->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err_flag{};
  bool ignore_radiation_lvalue;
  auto *_ignore_radiation{&ignore_radiation_lvalue};
  if (ignore_radiation.has_value()) {
    ignore_radiation_lvalue = ignore_radiation.value();
  } else {
    _ignore_radiation = nullptr;
  }
  bool make_map1_lvalue;
  auto *_make_map1{&make_map1_lvalue};
  if (make_map1.has_value()) {
    make_map1_lvalue = make_map1.value();
  } else {
    _make_map1 = nullptr;
  }
  bool init_to_edge_lvalue;
  auto *_init_to_edge{&init_to_edge_lvalue};
  if (init_to_edge.has_value()) {
    init_to_edge_lvalue = init_to_edge.value();
  } else {
    _init_to_edge = nullptr;
  }
  fortran_track1(/* void* */ start_orb.get_fortran_ptr(),
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
    BeamStruct &beam,
    EleStruct &ele,
    std::optional<CoordStructArray1D> centroid,
    std::optional<int> direction
) {
  bool _err{};
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  if (centroid) {
    _centroid_desc.data_ptr = centroid->data();
    _centroid_desc.dims[0] = centroid->size();
  } else {
    _centroid_desc.data_ptr = nullptr;
    _centroid_desc.dims[0] = 0;
  }
  _centroid_desc.strides[0] = 1;
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  fortran_track1_beam(/* void* */ beam.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* bool& */ _err,
                      /* Bmad::array_descriptor_t& */ _centroid_desc,
                      /* int* */ _direction);
  return _err;
}
Bmad::Track1Bmad Bmad::track1_bmad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  bool _err_flag{};
  TrackStruct _track;
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_bmad(/* void* */ orbit.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* void* */ param.get_fortran_ptr(),
                      /* bool& */ _err_flag,
                      /* void* */ _track.get_fortran_ptr(),
                      /* Bmad::array_descriptor_t& */ _mat6_desc,
                      /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return Track1Bmad{_err_flag, std::move(_track)};
}
bool Bmad::track1_bmad_photon(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param) {
  bool _err_flag{};
  fortran_track1_bmad_photon(/* void* */ orbit.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* bool& */ _err_flag);
  return _err_flag;
}
bool Bmad::track1_bunch(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<CoordStructArray1D> centroid,
    std::optional<int> direction,
    optional_ref<BunchTrackStruct> bunch_track
) {
  bool _err{};
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  if (centroid) {
    _centroid_desc.data_ptr = centroid->data();
    _centroid_desc.dims[0] = centroid->size();
  } else {
    _centroid_desc.data_ptr = nullptr;
    _centroid_desc.dims[0] = 0;
  }
  _centroid_desc.strides[0] = 1;
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track1_bunch(/* void* */ bunch.get_fortran_ptr(),
                       /* void* */ ele.get_fortran_ptr(),
                       /* bool& */ _err,
                       /* Bmad::array_descriptor_t& */ _centroid_desc,
                       /* int* */ _direction,
                       /* void* */ _bunch_track);
  return _err;
}
bool Bmad::track1_bunch_csr(
    BunchStruct &bunch,
    EleStruct &ele,
    CoordStructArray1D centroid,
    std::optional<double> s_start,
    std::optional<double> s_end,
    optional_ref<BunchTrackStruct> bunch_track
) {
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  _centroid_desc.data_ptr = centroid.data();
  _centroid_desc.dims[0] = centroid.size();
  _centroid_desc.strides[0] = 1;
  bool _err{};
  double s_start_lvalue;
  auto *_s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto *_s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track1_bunch_csr(/* void* */ bunch.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _centroid_desc,
                           /* bool& */ _err,
                           /* double* */ _s_start,
                           /* double* */ _s_end,
                           /* void* */ _bunch_track);
  return _err;
}
bool Bmad::track1_bunch_csr3d(
    BunchStruct &bunch,
    EleStruct &ele,
    CoordStructArray1D centroid,
    std::optional<double> s_start,
    std::optional<double> s_end,
    optional_ref<BunchTrackStruct> bunch_track
) {
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  _centroid_desc.data_ptr = centroid.data();
  _centroid_desc.dims[0] = centroid.size();
  _centroid_desc.strides[0] = 1;
  bool _err{};
  double s_start_lvalue;
  auto *_s_start{&s_start_lvalue};
  if (s_start.has_value()) {
    s_start_lvalue = s_start.value();
  } else {
    _s_start = nullptr;
  }
  double s_end_lvalue;
  auto *_s_end{&s_end_lvalue};
  if (s_end.has_value()) {
    s_end_lvalue = s_end.value();
  } else {
    _s_end = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track1_bunch_csr3d(/* void* */ bunch.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _centroid_desc,
                             /* bool& */ _err,
                             /* double* */ _s_start,
                             /* double* */ _s_end,
                             /* void* */ _bunch_track);
  return _err;
}
void Bmad::track1_bunch_hom(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<int> direction,
    optional_ref<BunchTrackStruct> bunch_track
) {
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track1_bunch_hom(/* void* */ bunch.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* int* */ _direction,
                           /* void* */ _bunch_track);
}
bool Bmad::track1_bunch_space_charge(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<bool> track_to_same_s,
    optional_ref<BunchTrackStruct> bunch_track
) {
  bool _err{};
  bool track_to_same_s_lvalue;
  auto *_track_to_same_s{&track_to_same_s_lvalue};
  if (track_to_same_s.has_value()) {
    track_to_same_s_lvalue = track_to_same_s.value();
  } else {
    _track_to_same_s = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track1_bunch_space_charge(/* void* */ bunch.get_fortran_ptr(),
                                    /* void* */ ele.get_fortran_ptr(),
                                    /* bool& */ _err,
                                    /* bool* */ _track_to_same_s,
                                    /* void* */ _bunch_track);
  return _err;
}
void Bmad::track1_crystal(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_crystal(/* void* */ ele.get_fortran_ptr(),
                         /* void* */ param.get_fortran_ptr(),
                         /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_diffraction_plate_or_mask(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit
) {
  fortran_track1_diffraction_plate_or_mask(/* void* */ ele.get_fortran_ptr(),
                                           /* void* */ param.get_fortran_ptr(),
                                           /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_high_energy_space_charge(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit
) {
  fortran_track1_high_energy_space_charge(/* void* */ ele.get_fortran_ptr(),
                                          /* void* */ param.get_fortran_ptr(),
                                          /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_lens(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_lens(/* void* */ ele.get_fortran_ptr(),
                      /* void* */ param.get_fortran_ptr(),
                      /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_linear(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param) {
  fortran_track1_linear(/* void* */ orbit.get_fortran_ptr(),
                        /* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr());
}
void Bmad::track1_lr_wake(BunchStruct &bunch, EleStruct &ele) {
  fortran_track1_lr_wake(/* void* */ bunch.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
void Bmad::track1_mad(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param) {
  fortran_track1_mad(/* void* */ orbit.get_fortran_ptr(),
                     /* void* */ ele.get_fortran_ptr(),
                     /* void* */ param.get_fortran_ptr());
}
void Bmad::track1_mirror(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_mirror(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_mosaic_crystal(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_mosaic_crystal(/* void* */ ele.get_fortran_ptr(),
                                /* void* */ param.get_fortran_ptr(),
                                /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_multilayer_mirror(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_multilayer_mirror(/* void* */ ele.get_fortran_ptr(),
                                   /* void* */ param.get_fortran_ptr(),
                                   /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_radiation(CoordStruct &orbit, EleStruct &ele, int edge) {
  fortran_track1_radiation(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* int& */ edge);
}
void Bmad::track1_radiation_center(
    CoordStruct &orbit,
    EleStruct &ele1,
    EleStruct &ele2,
    std::optional<bool> rad_damp,
    std::optional<bool> rad_fluct
) {
  bool rad_damp_lvalue;
  auto *_rad_damp{&rad_damp_lvalue};
  if (rad_damp.has_value()) {
    rad_damp_lvalue = rad_damp.value();
  } else {
    _rad_damp = nullptr;
  }
  bool rad_fluct_lvalue;
  auto *_rad_fluct{&rad_fluct_lvalue};
  if (rad_fluct.has_value()) {
    rad_fluct_lvalue = rad_fluct.value();
  } else {
    _rad_fluct = nullptr;
  }
  fortran_track1_radiation_center(/* void* */ orbit.get_fortran_ptr(),
                                  /* void* */ ele1.get_fortran_ptr(),
                                  /* void* */ ele2.get_fortran_ptr(),
                                  /* bool* */ _rad_damp,
                                  /* bool* */ _rad_fluct);
}
Bmad::Track1RungeKutta Bmad::track1_runge_kutta(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  bool _err_flag{};
  TrackStruct _track;
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_runge_kutta(/* void* */ orbit.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* bool& */ _err_flag,
                             /* void* */ _track.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _mat6_desc,
                             /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
  return Track1RungeKutta{_err_flag, std::move(_track)};
}
void Bmad::track1_sample(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit) {
  fortran_track1_sample(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr());
}
void Bmad::track1_spin(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &end_orb,
    std::optional<bool> make_quaternion
) {
  bool make_quaternion_lvalue;
  auto *_make_quaternion{&make_quaternion_lvalue};
  if (make_quaternion.has_value()) {
    make_quaternion_lvalue = make_quaternion.value();
  } else {
    _make_quaternion = nullptr;
  }
  fortran_track1_spin(/* void* */ start_orb.get_fortran_ptr(),
                      /* void* */ ele.get_fortran_ptr(),
                      /* void* */ param.get_fortran_ptr(),
                      /* void* */ end_orb.get_fortran_ptr(),
                      /* bool* */ _make_quaternion);
}
void Bmad::track1_spin_integration(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &end_orb
) {
  fortran_track1_spin_integration(/* void* */ start_orb.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* void* */ end_orb.get_fortran_ptr());
}
CoordStruct
Bmad::track1_spin_taylor(CoordStruct &start_orb, EleStruct &ele, LatParamStruct &param) {
  CoordStruct _end_orb;
  fortran_track1_spin_taylor(/* void* */ start_orb.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* void* */ _end_orb.get_fortran_ptr());
  return std::move(_end_orb);
}
void Bmad::track1_sr_wake(BunchStruct &bunch, EleStruct &ele) {
  fortran_track1_sr_wake(/* void* */ bunch.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
TrackStruct Bmad::track1_symp_lie_ptc(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param) {
  TrackStruct _track;
  fortran_track1_symp_lie_ptc(/* void* */ orbit.get_fortran_ptr(),
                              /* void* */ ele.get_fortran_ptr(),
                              /* void* */ param.get_fortran_ptr(),
                              /* void* */ _track.get_fortran_ptr());
  return std::move(_track);
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track1_taylor(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<TaylorStructArray1D> taylor,
    std::optional<bool> make_matrix
) {
  // taylor: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _taylor_desc;
  _taylor_desc.rank = 1;
  if (taylor) {
    _taylor_desc.data_ptr = taylor->data();
    _taylor_desc.dims[0] = taylor->size();
  } else {
    _taylor_desc.data_ptr = nullptr;
    _taylor_desc.dims[0] = 0;
  }
  _taylor_desc.strides[0] = 1;
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track1_taylor(/* void* */ orbit.get_fortran_ptr(),
                        /* void* */ ele.get_fortran_ptr(),
                        /* Bmad::array_descriptor_t& */ _taylor_desc,
                        /* Bmad::array_descriptor_t& */ _mat6_desc,
                        /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::Track1TimeRungeKutta Bmad::track1_time_runge_kutta(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<double> t_end,
    optional_ref<double> dt_step
) {
  bool _err_flag{};
  TrackStruct _track;
  double t_end_lvalue;
  auto *_t_end{&t_end_lvalue};
  if (t_end.has_value()) {
    t_end_lvalue = t_end.value();
  } else {
    _t_end = nullptr;
  }
  auto *_dt_step = dt_step.has_value() ? &dt_step->get() : nullptr; // inout, optional
  fortran_track1_time_runge_kutta(/* void* */ orbit.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* bool& */ _err_flag,
                                  /* void* */ _track.get_fortran_ptr(),
                                  /* double* */ _t_end,
                                  /* double* */ _dt_step);
  return Track1TimeRungeKutta{_err_flag, std::move(_track)};
}
Bmad::TrackABeambeam Bmad::track_a_beambeam(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  TrackStruct _track;
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_beambeam(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* void* */ _track.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _mat6_desc,
                           /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return TrackABeambeam{std::move(_track), mat6};
}
void Bmad::track_a_bend(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_bend(/* void* */ orbit.get_fortran_ptr(),
                       /* void* */ ele.get_fortran_ptr(),
                       /* void* */ param.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _mat6_desc,
                       /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_bend_photon(CoordStruct &orb, EleStruct &ele, double length) {
  fortran_track_a_bend_photon(/* void* */ orb.get_fortran_ptr(),
                              /* void* */ ele.get_fortran_ptr(),
                              /* double& */ length);
}
void Bmad::track_a_capillary(CoordStruct &orb, EleStruct &ele) {
  fortran_track_a_capillary(/* void* */ orb.get_fortran_ptr(), /* void* */ ele.get_fortran_ptr());
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_converter(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_converter(/* void* */ orbit.get_fortran_ptr(),
                            /* void* */ ele.get_fortran_ptr(),
                            /* void* */ param.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _mat6_desc,
                            /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_crab_cavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_crab_cavity(/* void* */ orbit.get_fortran_ptr(),
                              /* void* */ ele.get_fortran_ptr(),
                              /* void* */ param.get_fortran_ptr(),
                              /* Bmad::array_descriptor_t& */ _mat6_desc,
                              /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_drift(
    CoordStruct &orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    std::optional<int> ele_orientation,
    std::optional<bool> include_ref_motion,
    optional_ref<double> time
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  int ele_orientation_lvalue;
  auto *_ele_orientation{&ele_orientation_lvalue};
  if (ele_orientation.has_value()) {
    ele_orientation_lvalue = ele_orientation.value();
  } else {
    _ele_orientation = nullptr;
  }
  bool include_ref_motion_lvalue;
  auto *_include_ref_motion{&include_ref_motion_lvalue};
  if (include_ref_motion.has_value()) {
    include_ref_motion_lvalue = include_ref_motion.value();
  } else {
    _include_ref_motion = nullptr;
  }
  auto *_time = time.has_value() ? &time->get() : nullptr; // inout, optional
  fortran_track_a_drift(/* void* */ orb.get_fortran_ptr(),
                        /* double& */ length,
                        /* Bmad::array_descriptor_t& */ _mat6_desc,
                        /* bool* */ _make_matrix,
                        /* int* */ _ele_orientation,
                        /* bool* */ _include_ref_motion,
                        /* double* */ _time);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_drift_photon(CoordStruct &orb, double length, bool phase_relative_to_ref) {
  fortran_track_a_drift_photon(/* void* */ orb.get_fortran_ptr(),
                               /* double& */ length,
                               /* bool& */ phase_relative_to_ref);
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_foil(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_foil(/* void* */ orbit.get_fortran_ptr(),
                       /* void* */ ele.get_fortran_ptr(),
                       /* void* */ param.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _mat6_desc,
                       /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_gkicker(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_gkicker(/* void* */ orbit.get_fortran_ptr(),
                          /* void* */ ele.get_fortran_ptr(),
                          /* void* */ param.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                          /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_lcavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_lcavity(/* void* */ orbit.get_fortran_ptr(),
                          /* void* */ ele.get_fortran_ptr(),
                          /* void* */ param.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                          /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
void Bmad::track_a_lcavity_old(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_lcavity_old(/* void* */ orbit.get_fortran_ptr(),
                              /* void* */ ele.get_fortran_ptr(),
                              /* void* */ param.get_fortran_ptr(),
                              /* Bmad::array_descriptor_t& */ _mat6_desc,
                              /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_mask(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_mask(/* void* */ orbit.get_fortran_ptr(),
                       /* void* */ ele.get_fortran_ptr(),
                       /* void* */ param.get_fortran_ptr(),
                       /* Bmad::array_descriptor_t& */ _mat6_desc,
                       /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_match(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> err_flag,
    std::optional<bool> make_matrix
) {
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_match(/* void* */ orbit.get_fortran_ptr(),
                        /* void* */ ele.get_fortran_ptr(),
                        /* void* */ param.get_fortran_ptr(),
                        /* bool* */ _err_flag,
                        /* Bmad::array_descriptor_t& */ _mat6_desc,
                        /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::TrackAPatch Bmad::track_a_patch(
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<bool> drift_to_exit,
    std::optional<bool> track_spin,
    std::optional<bool> make_matrix
) {
  bool drift_to_exit_lvalue;
  auto *_drift_to_exit{&drift_to_exit_lvalue};
  if (drift_to_exit.has_value()) {
    drift_to_exit_lvalue = drift_to_exit.value();
  } else {
    _drift_to_exit = nullptr;
  }
  double _s_ent{};
  double _ds_ref{};
  bool track_spin_lvalue;
  auto *_track_spin{&track_spin_lvalue};
  if (track_spin.has_value()) {
    track_spin_lvalue = track_spin.value();
  } else {
    _track_spin = nullptr;
  }
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_patch(/* void* */ ele.get_fortran_ptr(),
                        /* void* */ orbit.get_fortran_ptr(),
                        /* bool* */ _drift_to_exit,
                        /* double& */ _s_ent,
                        /* double& */ _ds_ref,
                        /* bool* */ _track_spin,
                        /* Bmad::array_descriptor_t& */ _mat6_desc,
                        /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return TrackAPatch{_s_ent, _ds_ref, mat6};
}
void Bmad::track_a_patch_photon(
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<bool> drift_to_exit,
    std::optional<bool> use_z_pos
) {
  bool drift_to_exit_lvalue;
  auto *_drift_to_exit{&drift_to_exit_lvalue};
  if (drift_to_exit.has_value()) {
    drift_to_exit_lvalue = drift_to_exit.value();
  } else {
    _drift_to_exit = nullptr;
  }
  bool use_z_pos_lvalue;
  auto *_use_z_pos{&use_z_pos_lvalue};
  if (use_z_pos.has_value()) {
    use_z_pos_lvalue = use_z_pos.value();
  } else {
    _use_z_pos = nullptr;
  }
  fortran_track_a_patch_photon(/* void* */ ele.get_fortran_ptr(),
                               /* void* */ orbit.get_fortran_ptr(),
                               /* bool* */ _drift_to_exit,
                               /* bool* */ _use_z_pos);
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_pickup(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> err_flag,
    std::optional<bool> make_matrix
) {
  bool err_flag_lvalue;
  auto *_err_flag{&err_flag_lvalue};
  if (err_flag.has_value()) {
    err_flag_lvalue = err_flag.value();
  } else {
    _err_flag = nullptr;
  }
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_pickup(/* void* */ orbit.get_fortran_ptr(),
                         /* void* */ ele.get_fortran_ptr(),
                         /* void* */ param.get_fortran_ptr(),
                         /* bool* */ _err_flag,
                         /* Bmad::array_descriptor_t& */ _mat6_desc,
                         /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_quadrupole(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_quadrupole(/* void* */ orbit.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* void* */ param.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _mat6_desc,
                             /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_rfcavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_rfcavity(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _mat6_desc,
                           /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_sad_mult(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_sad_mult(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _mat6_desc,
                           /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_sol_quad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_sol_quad(/* void* */ orbit.get_fortran_ptr(),
                           /* void* */ ele.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _mat6_desc,
                           /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
void Bmad::track_a_thick_multipole(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix
) {
  // mat6: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  if (mat6.has_value()) {
    matrix_to_vec(mat6.value(), _mat6_vec);
    _mat6_desc.dims[0] = 6;
    _mat6_desc.dims[1] = 6;
  } else {
    _mat6_desc.data_ptr = nullptr;
  }
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_thick_multipole(/* void* */ orbit.get_fortran_ptr(),
                                  /* void* */ ele.get_fortran_ptr(),
                                  /* void* */ param.get_fortran_ptr(),
                                  /* Bmad::array_descriptor_t& */ _mat6_desc,
                                  /* bool* */ _make_matrix);
  if (mat6.has_value())
    vec_to_matrix(_mat6_vec, mat6.value());
}
std::optional<FixedArray2D<Real, 6, 6>> Bmad::track_a_wiggler(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix
) {
  // mat6: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _mat6_desc;
  _mat6_desc.rank = 2;
  FixedArray2D<Real, 6, 6> mat6;
  double _mat6_vec[6 * 6];
  _mat6_desc.data_ptr = _mat6_vec;
  _mat6_desc.dims[0] = 6;
  _mat6_desc.dims[1] = 6;
  bool make_matrix_lvalue;
  auto *_make_matrix{&make_matrix_lvalue};
  if (make_matrix.has_value()) {
    make_matrix_lvalue = make_matrix.value();
  } else {
    _make_matrix = nullptr;
  }
  fortran_track_a_wiggler(/* void* */ orbit.get_fortran_ptr(),
                          /* void* */ ele.get_fortran_ptr(),
                          /* void* */ param.get_fortran_ptr(),
                          /* Bmad::array_descriptor_t& */ _mat6_desc,
                          /* bool* */ _make_matrix);
  vec_to_matrix(_mat6_vec, mat6);
  return mat6;
}
Bmad::TrackAZeroLengthElement
Bmad::track_a_zero_length_element(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param) {
  bool _err_flag{};
  TrackStruct _track;
  fortran_track_a_zero_length_element(/* void* */ orbit.get_fortran_ptr(),
                                      /* void* */ ele.get_fortran_ptr(),
                                      /* void* */ param.get_fortran_ptr(),
                                      /* bool& */ _err_flag,
                                      /* void* */ _track.get_fortran_ptr());
  return TrackAZeroLengthElement{_err_flag, std::move(_track)};
}
Bmad::TrackAll Bmad::track_all(
    LatStruct &lat,
    CoordStructAlloc1D orbit,
    std::optional<int> ix_branch,
    std::optional<bool> init_lost
) {
  // intent=inout allocatable type array
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  bool _err_flag{};
  // intent=out allocatable type array
  auto orbit0{CoordStructAlloc1D()};
  bool init_lost_lvalue;
  auto *_init_lost{&init_lost_lvalue};
  if (init_lost.has_value()) {
    init_lost_lvalue = init_lost.value();
  } else {
    _init_lost = nullptr;
  }
  fortran_track_all(/* void* */ lat.get_fortran_ptr(),
                    /* void* */ orbit.get_fortran_ptr(),
                    /* int* */ _ix_branch,
                    /* int& */ _track_state,
                    /* bool& */ _err_flag,
                    /* void* */ orbit0.get_fortran_ptr(),
                    /* bool* */ _init_lost);
  return TrackAll{_track_state, _err_flag, std::move(orbit0)};
}
bool Bmad::track_beam(
    LatStruct &lat,
    BeamStruct &beam,
    optional_ref<EleStruct> ele1,
    optional_ref<EleStruct> ele2,
    std::optional<CoordStructArray1D> centroid,
    std::optional<int> direction,
    std::optional<BunchTrackStructArray1D> bunch_tracks
) {
  auto *_ele1 = ele1.has_value() ? ele1->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_ele2 = ele2.has_value() ? ele2->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err{};
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  if (centroid) {
    _centroid_desc.data_ptr = centroid->data();
    _centroid_desc.dims[0] = centroid->size();
  } else {
    _centroid_desc.data_ptr = nullptr;
    _centroid_desc.dims[0] = 0;
  }
  _centroid_desc.strides[0] = 1;
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  // bunch_tracks: BunchTrackStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _bunch_tracks_desc;
  _bunch_tracks_desc.rank = 1;
  if (bunch_tracks) {
    _bunch_tracks_desc.data_ptr = bunch_tracks->data();
    _bunch_tracks_desc.dims[0] = bunch_tracks->size();
  } else {
    _bunch_tracks_desc.data_ptr = nullptr;
    _bunch_tracks_desc.dims[0] = 0;
  }
  _bunch_tracks_desc.strides[0] = 1;
  fortran_track_beam(/* void* */ lat.get_fortran_ptr(),
                     /* void* */ beam.get_fortran_ptr(),
                     /* void* */ _ele1,
                     /* void* */ _ele2,
                     /* bool& */ _err,
                     /* Bmad::array_descriptor_t& */ _centroid_desc,
                     /* int* */ _direction,
                     /* Bmad::array_descriptor_t& */ _bunch_tracks_desc);
  return _err;
}
bool Bmad::track_bunch(
    LatStruct &lat,
    BunchStruct &bunch,
    optional_ref<EleStruct> ele1,
    optional_ref<EleStruct> ele2,
    std::optional<CoordStructArray1D> centroid,
    std::optional<int> direction,
    optional_ref<BunchTrackStruct> bunch_track
) {
  auto *_ele1 = ele1.has_value() ? ele1->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_ele2 = ele2.has_value() ? ele2->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err{};
  // centroid: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _centroid_desc;
  _centroid_desc.rank = 1;
  if (centroid) {
    _centroid_desc.data_ptr = centroid->data();
    _centroid_desc.dims[0] = centroid->size();
  } else {
    _centroid_desc.data_ptr = nullptr;
    _centroid_desc.dims[0] = 0;
  }
  _centroid_desc.strides[0] = 1;
  int direction_lvalue;
  auto *_direction{&direction_lvalue};
  if (direction.has_value()) {
    direction_lvalue = direction.value();
  } else {
    _direction = nullptr;
  }
  auto *_bunch_track =
      bunch_track.has_value() ? bunch_track->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_track_bunch(/* void* */ lat.get_fortran_ptr(),
                      /* void* */ bunch.get_fortran_ptr(),
                      /* void* */ _ele1,
                      /* void* */ _ele2,
                      /* bool& */ _err,
                      /* Bmad::array_descriptor_t& */ _centroid_desc,
                      /* int* */ _direction,
                      /* void* */ _bunch_track);
  return _err;
}
void Bmad::track_bunch_time(
    BunchStruct &bunch,
    BranchStruct &branch,
    double t_end,
    double s_end,
    std::optional<FArray1D<Real>> dt_step,
    std::optional<EmFieldStructArray1D> extra_field
) {
  // dt_step: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _dt_step_desc;
  _dt_step_desc.rank = 1;
  if (dt_step.has_value()) {
    _dt_step_desc.data_ptr = dt_step->data();
    _dt_step_desc.dims[0] = dt_step->size();
  } else {
    _dt_step_desc.data_ptr = nullptr;
    _dt_step_desc.dims[0] = 0;
  }
  // extra_field: EmFieldStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _extra_field_desc;
  _extra_field_desc.rank = 1;
  if (extra_field) {
    _extra_field_desc.data_ptr = extra_field->data();
    _extra_field_desc.dims[0] = extra_field->size();
  } else {
    _extra_field_desc.data_ptr = nullptr;
    _extra_field_desc.dims[0] = 0;
  }
  _extra_field_desc.strides[0] = 1;
  fortran_track_bunch_time(/* void* */ bunch.get_fortran_ptr(),
                           /* void* */ branch.get_fortran_ptr(),
                           /* double& */ t_end,
                           /* double& */ s_end,
                           /* Bmad::array_descriptor_t& */ _dt_step_desc,
                           /* Bmad::array_descriptor_t& */ _extra_field_desc);
}
void Bmad::track_bunch_to_s(BunchStruct &bunch, double s, BranchStruct &branch) {
  fortran_track_bunch_to_s(/* void* */ bunch.get_fortran_ptr(),
                           /* double& */ s,
                           /* void* */ branch.get_fortran_ptr());
}
void Bmad::track_bunch_to_t(BunchStruct &bunch, double t_target, BranchStruct &branch) {
  fortran_track_bunch_to_t(/* void* */ bunch.get_fortran_ptr(),
                           /* double& */ t_target,
                           /* void* */ branch.get_fortran_ptr());
}
void Bmad::track_complex_taylor(
    FArray1D<Complex> &start_orb,
    ComplexTaylorStructArray1D complex_taylor,
    FArray1D<Complex> &end_orb
) {
  // start_orb: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _start_orb_desc;
  _start_orb_desc.rank = 1;
  _start_orb_desc.data_ptr = start_orb.data();
  _start_orb_desc.dims[0] = start_orb.size();
  // complex_taylor: ComplexTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_desc;
  _complex_taylor_desc.rank = 1;
  _complex_taylor_desc.data_ptr = complex_taylor.data();
  _complex_taylor_desc.dims[0] = complex_taylor.size();
  _complex_taylor_desc.strides[0] = 1;
  // end_orb: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _end_orb_desc;
  _end_orb_desc.rank = 1;
  _end_orb_desc.data_ptr = end_orb.data();
  _end_orb_desc.dims[0] = end_orb.size();
  fortran_track_complex_taylor(
      /* Bmad::array_descriptor_t& */ _start_orb_desc,
      /* Bmad::array_descriptor_t& */ _complex_taylor_desc,
      /* Bmad::array_descriptor_t& */ _end_orb_desc
  );
}
Bmad::TrackFromSToS Bmad::track_from_s_to_s(
    LatStruct &lat,
    double s_start,
    double s_end,
    CoordStruct &orbit_start,
    std::optional<int> ix_branch,
    std::optional<int> ix_ele_end
) {
  CoordStruct _orbit_end;
  // intent=out allocatable type array
  auto all_orb{CoordStructAlloc1D()};
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  int ix_ele_end_lvalue;
  auto *_ix_ele_end{&ix_ele_end_lvalue};
  if (ix_ele_end.has_value()) {
    ix_ele_end_lvalue = ix_ele_end.value();
  } else {
    _ix_ele_end = nullptr;
  }
  fortran_track_from_s_to_s(/* void* */ lat.get_fortran_ptr(),
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
    LatStruct &lat,
    CoordStructArray1D orbit,
    int ix_start,
    int ix_end,
    int direction,
    std::optional<int> ix_branch
) {
  // orbit: CoordStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  _orbit_desc.data_ptr = orbit.data();
  _orbit_desc.dims[0] = orbit.size();
  _orbit_desc.strides[0] = 1;
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  int _track_state{};
  fortran_track_many(/* void* */ lat.get_fortran_ptr(),
                     /* Bmad::array_descriptor_t& */ _orbit_desc,
                     /* int& */ ix_start,
                     /* int& */ ix_end,
                     /* int& */ direction,
                     /* int* */ _ix_branch,
                     /* int& */ _track_state);
  return _track_state;
}
FixedArray2D<Real, 3, 3>
Bmad::track_to_surface(EleStruct &ele, CoordStruct &orbit, LatParamStruct &param) {
  // w_surface: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_surface_desc;
  _w_surface_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_surface;
  double _w_surface_vec[3 * 3];
  _w_surface_desc.data_ptr = _w_surface_vec;
  _w_surface_desc.dims[0] = 3;
  _w_surface_desc.dims[1] = 3;
  fortran_track_to_surface(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ orbit.get_fortran_ptr(),
                           /* void* */ param.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _w_surface_desc);
  vec_to_matrix(_w_surface_vec, w_surface);
  return w_surface;
}
Bmad::TrackUntilDead Bmad::track_until_dead(CoordStruct &start_orb, LatStruct &lat) {
  CoordStruct _end_orb;
  TrackStruct _track;
  fortran_track_until_dead(/* void* */ start_orb.get_fortran_ptr(),
                           /* void* */ lat.get_fortran_ptr(),
                           /* void* */ _end_orb.get_fortran_ptr(),
                           /* void* */ _track.get_fortran_ptr());
  return TrackUntilDead{std::move(_end_orb), std::move(_track)};
}
Bmad::TrackingRadMapSetup
Bmad::tracking_rad_map_setup(EleStruct &ele, double tollerance, int ref_edge) {
  RadMapStruct _rad_map;
  bool _err_flag{};
  fortran_tracking_rad_map_setup(/* void* */ ele.get_fortran_ptr(),
                                 /* double& */ tollerance,
                                 /* int& */ ref_edge,
                                 /* void* */ _rad_map.get_fortran_ptr(),
                                 /* bool& */ _err_flag);
  return TrackingRadMapSetup{std::move(_rad_map), _err_flag};
}
std::optional<AcKickerStruct> Bmad::transfer_ac_kick(AcKickerStruct &ac_in) {
  auto _ac_in = &ac_in; // input, required, pointer
  void *_ac_out;
  fortran_transfer_ac_kick(/* void* */ &ac_in, /* void* */ &_ac_out);
  return std::move((_ac_out ? std::make_optional<AcKickerStruct>(_ac_out) : std::nullopt));
}
BranchStruct Bmad::transfer_branch(BranchStruct &branch1) {
  BranchStruct _branch2;
  fortran_transfer_branch(/* void* */ branch1.get_fortran_ptr(),
                          /* void* */ _branch2.get_fortran_ptr());
  return std::move(_branch2);
}
BranchStruct Bmad::transfer_branch_parameters(BranchStruct &branch_in) {
  BranchStruct _branch_out;
  fortran_transfer_branch_parameters(/* void* */ branch_in.get_fortran_ptr(),
                                     /* void* */ _branch_out.get_fortran_ptr());
  return std::move(_branch_out);
}
void Bmad::transfer_branches(BranchStructArray1D branch1, BranchStructArray1D branch2) {
  // branch1: BranchStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _branch1_desc;
  _branch1_desc.rank = 1;
  _branch1_desc.data_ptr = branch1.data();
  _branch1_desc.dims[0] = branch1.size();
  _branch1_desc.strides[0] = 1;
  // branch2: BranchStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _branch2_desc;
  _branch2_desc.rank = 1;
  _branch2_desc.data_ptr = branch2.data();
  _branch2_desc.dims[0] = branch2.size();
  _branch2_desc.strides[0] = 1;
  fortran_transfer_branches(
      /* Bmad::array_descriptor_t& */ _branch1_desc,
      /* Bmad::array_descriptor_t& */ _branch2_desc
  );
}
EleStruct Bmad::transfer_ele(EleStruct &ele1, std::optional<bool> nullify_pointers) {
  EleStruct _ele2;
  bool nullify_pointers_lvalue;
  auto *_nullify_pointers{&nullify_pointers_lvalue};
  if (nullify_pointers.has_value()) {
    nullify_pointers_lvalue = nullify_pointers.value();
  } else {
    _nullify_pointers = nullptr;
  }
  fortran_transfer_ele(/* void* */ ele1.get_fortran_ptr(),
                       /* void* */ _ele2.get_fortran_ptr(),
                       /* bool* */ _nullify_pointers);
  return std::move(_ele2);
}
EleStruct Bmad::transfer_ele_taylor(EleStruct &ele_in, std::optional<int> taylor_order) {
  EleStruct _ele_out;
  int taylor_order_lvalue;
  auto *_taylor_order{&taylor_order_lvalue};
  if (taylor_order.has_value()) {
    taylor_order_lvalue = taylor_order.value();
  } else {
    _taylor_order = nullptr;
  }
  fortran_transfer_ele_taylor(/* void* */ ele_in.get_fortran_ptr(),
                              /* void* */ _ele_out.get_fortran_ptr(),
                              /* int* */ _taylor_order);
  return std::move(_ele_out);
}
void Bmad::transfer_eles(EleStructArray1D ele1, EleStructArray1D ele2) {
  // ele1: EleStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele1_desc;
  _ele1_desc.rank = 1;
  _ele1_desc.data_ptr = ele1.data();
  _ele1_desc.dims[0] = ele1.size();
  _ele1_desc.strides[0] = 1;
  // ele2: EleStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele2_desc;
  _ele2_desc.rank = 1;
  _ele2_desc.data_ptr = ele2.data();
  _ele2_desc.dims[0] = ele2.size();
  _ele2_desc.strides[0] = 1;
  fortran_transfer_eles(
      /* Bmad::array_descriptor_t& */ _ele1_desc,
      /* Bmad::array_descriptor_t& */ _ele2_desc
  );
}
EleStruct Bmad::transfer_fieldmap(EleStruct &ele_in, int who) {
  EleStruct _ele_out;
  fortran_transfer_fieldmap(/* void* */ ele_in.get_fortran_ptr(),
                            /* void* */ _ele_out.get_fortran_ptr(),
                            /* int& */ who);
  return std::move(_ele_out);
}
bool Bmad::transfer_fixer_params(
    EleStruct &fixer,
    bool to_stored,
    optional_ref<CoordStruct> orbit,
    std::optional<std::string> who
) {
  auto *_orbit = orbit.has_value() ? orbit->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_who = who.has_value() ? who->c_str() : nullptr;
  bool _is_ok{};
  fortran_transfer_fixer_params(/* void* */ fixer.get_fortran_ptr(),
                                /* bool& */ to_stored,
                                /* void* */ _orbit,
                                /* const char* */ _who,
                                /* bool& */ _is_ok);
  return _is_ok;
}
LatStruct Bmad::transfer_lat(LatStruct &lat1) {
  LatStruct _lat2;
  fortran_transfer_lat(/* void* */ lat1.get_fortran_ptr(), /* void* */ _lat2.get_fortran_ptr());
  return std::move(_lat2);
}
LatStruct Bmad::transfer_lat_parameters(LatStruct &lat_in) {
  LatStruct _lat_out;
  fortran_transfer_lat_parameters(/* void* */ lat_in.get_fortran_ptr(),
                                  /* void* */ _lat_out.get_fortran_ptr());
  return std::move(_lat_out);
}
bool Bmad::transfer_map_calc(
    LatStruct &lat,
    TaylorStructArray1D orb_map,
    std::optional<int> ix1,
    std::optional<int> ix2,
    optional_ref<CoordStruct> ref_orb,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start,
    std::optional<bool> concat_if_possible,
    std::optional<TaylorStructArray1D> spin_map
) {
  // orb_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb_map_desc;
  _orb_map_desc.rank = 1;
  _orb_map_desc.data_ptr = orb_map.data();
  _orb_map_desc.dims[0] = orb_map.size();
  _orb_map_desc.strides[0] = 1;
  bool _err_flag{};
  int ix1_lvalue;
  auto *_ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto *_ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  auto *_ref_orb =
      ref_orb.has_value() ? ref_orb->get().get_fortran_ptr() : nullptr; // input, optional
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto *_one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto *_unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  bool concat_if_possible_lvalue;
  auto *_concat_if_possible{&concat_if_possible_lvalue};
  if (concat_if_possible.has_value()) {
    concat_if_possible_lvalue = concat_if_possible.value();
  } else {
    _concat_if_possible = nullptr;
  }
  // spin_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_map_desc;
  _spin_map_desc.rank = 1;
  if (spin_map) {
    _spin_map_desc.data_ptr = spin_map->data();
    _spin_map_desc.dims[0] = spin_map->size();
  } else {
    _spin_map_desc.data_ptr = nullptr;
    _spin_map_desc.dims[0] = 0;
  }
  _spin_map_desc.strides[0] = 1;
  fortran_transfer_map_calc(/* void* */ lat.get_fortran_ptr(),
                            /* Bmad::array_descriptor_t& */ _orb_map_desc,
                            /* bool& */ _err_flag,
                            /* int* */ _ix1,
                            /* int* */ _ix2,
                            /* void* */ _ref_orb,
                            /* int* */ _ix_branch,
                            /* bool* */ _one_turn,
                            /* bool* */ _unit_start,
                            /* bool* */ _concat_if_possible,
                            /* Bmad::array_descriptor_t& */ _spin_map_desc);
  return _err_flag;
}
Bmad::TransferMapFromSToS Bmad::transfer_map_from_s_to_s(
    LatStruct &lat,
    TaylorStructArray1D t_map,
    std::optional<double> s1,
    std::optional<double> s2,
    optional_ref<CoordStruct> ref_orb_in,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn,
    std::optional<bool> unit_start,
    std::optional<bool> concat_if_possible,
    std::optional<TaylorStructArray1D> spin_map
) {
  // t_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _t_map_desc;
  _t_map_desc.rank = 1;
  _t_map_desc.data_ptr = t_map.data();
  _t_map_desc.dims[0] = t_map.size();
  _t_map_desc.strides[0] = 1;
  double s1_lvalue;
  auto *_s1{&s1_lvalue};
  if (s1.has_value()) {
    s1_lvalue = s1.value();
  } else {
    _s1 = nullptr;
  }
  double s2_lvalue;
  auto *_s2{&s2_lvalue};
  if (s2.has_value()) {
    s2_lvalue = s2.value();
  } else {
    _s2 = nullptr;
  }
  auto *_ref_orb_in =
      ref_orb_in.has_value() ? ref_orb_in->get().get_fortran_ptr() : nullptr; // input, optional
  CoordStruct _ref_orb_out;
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto *_one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  bool unit_start_lvalue;
  auto *_unit_start{&unit_start_lvalue};
  if (unit_start.has_value()) {
    unit_start_lvalue = unit_start.value();
  } else {
    _unit_start = nullptr;
  }
  bool _err_flag{};
  bool concat_if_possible_lvalue;
  auto *_concat_if_possible{&concat_if_possible_lvalue};
  if (concat_if_possible.has_value()) {
    concat_if_possible_lvalue = concat_if_possible.value();
  } else {
    _concat_if_possible = nullptr;
  }
  // spin_map: TaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spin_map_desc;
  _spin_map_desc.rank = 1;
  if (spin_map) {
    _spin_map_desc.data_ptr = spin_map->data();
    _spin_map_desc.dims[0] = spin_map->size();
  } else {
    _spin_map_desc.data_ptr = nullptr;
    _spin_map_desc.dims[0] = 0;
  }
  _spin_map_desc.strides[0] = 1;
  fortran_transfer_map_from_s_to_s(/* void* */ lat.get_fortran_ptr(),
                                   /* Bmad::array_descriptor_t& */ _t_map_desc,
                                   /* double* */ _s1,
                                   /* double* */ _s2,
                                   /* void* */ _ref_orb_in,
                                   /* void* */ _ref_orb_out.get_fortran_ptr(),
                                   /* int* */ _ix_branch,
                                   /* bool* */ _one_turn,
                                   /* bool* */ _unit_start,
                                   /* bool& */ _err_flag,
                                   /* bool* */ _concat_if_possible,
                                   /* Bmad::array_descriptor_t& */ _spin_map_desc);
  return TransferMapFromSToS{std::move(_ref_orb_out), _err_flag};
}
FixedArray2D<Real, 2, 2> Bmad::transfer_mat2_from_twiss(TwissStruct &twiss1, TwissStruct &twiss2) {
  // mat: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _mat_desc;
  _mat_desc.rank = 2;
  FixedArray2D<Real, 2, 2> mat;
  double _mat_vec[2 * 2];
  _mat_desc.data_ptr = _mat_vec;
  _mat_desc.dims[0] = 2;
  _mat_desc.dims[1] = 2;
  fortran_transfer_mat2_from_twiss(/* void* */ twiss1.get_fortran_ptr(),
                                   /* void* */ twiss2.get_fortran_ptr(),
                                   /* Bmad::array_descriptor_t& */ _mat_desc);
  vec_to_matrix(_mat_vec, mat);
  return mat;
}
FixedArray2D<Real, 6, 6> Bmad::transfer_mat_from_twiss(
    EleStruct &ele1,
    EleStruct &ele2,
    FixedArray1D<Real, 6> orb1,
    FixedArray1D<Real, 6> orb2
) {
  // orb1: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orb1_desc;
  _orb1_desc.rank = 1;
  _orb1_desc.data_ptr = orb1.data();
  _orb1_desc.dims[0] = orb1.size();
  // orb2: in NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _orb2_desc;
  _orb2_desc.rank = 1;
  _orb2_desc.data_ptr = orb2.data();
  _orb2_desc.dims[0] = orb2.size();
  // m: out NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _m_desc;
  _m_desc.rank = 2;
  FixedArray2D<Real, 6, 6> m;
  double _m_vec[6 * 6];
  _m_desc.data_ptr = _m_vec;
  _m_desc.dims[0] = 6;
  _m_desc.dims[1] = 6;
  fortran_transfer_mat_from_twiss(/* void* */ ele1.get_fortran_ptr(),
                                  /* void* */ ele2.get_fortran_ptr(),
                                  /* Bmad::array_descriptor_t& */ _orb1_desc,
                                  /* Bmad::array_descriptor_t& */ _orb2_desc,
                                  /* Bmad::array_descriptor_t& */ _m_desc);
  vec_to_matrix(_m_vec, m);
  return m;
}
void Bmad::transfer_matrix_calc(
    LatStruct &lat,
    FixedArray2D<Real, 6, 6> xfer_mat,
    std::optional<FixedArray1D<Real, 6>> xfer_vec,
    std::optional<int> ix1,
    std::optional<int> ix2,
    std::optional<int> ix_branch,
    std::optional<bool> one_turn
) {
  // xfer_mat: inout NOT (CppWrapperGeneralArgumentArray) (['6', '6'])
  Bmad::array_descriptor_t _xfer_mat_desc;
  _xfer_mat_desc.rank = 2;
  double _xfer_mat_vec[6 * 6];
  _xfer_mat_desc.data_ptr = _xfer_mat_vec;
  _xfer_mat_desc.dims[0] = 6;
  _xfer_mat_desc.dims[1] = 6;
  matrix_to_vec(xfer_mat, _xfer_mat_vec);
  // xfer_vec: inout NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _xfer_vec_desc;
  _xfer_vec_desc.rank = 1;
  if (xfer_vec.has_value()) {
    _xfer_vec_desc.data_ptr = xfer_vec->data();
    _xfer_vec_desc.dims[0] = xfer_vec->size();
  } else {
    _xfer_vec_desc.data_ptr = nullptr;
    _xfer_vec_desc.dims[0] = 0;
  }
  int ix1_lvalue;
  auto *_ix1{&ix1_lvalue};
  if (ix1.has_value()) {
    ix1_lvalue = ix1.value();
  } else {
    _ix1 = nullptr;
  }
  int ix2_lvalue;
  auto *_ix2{&ix2_lvalue};
  if (ix2.has_value()) {
    ix2_lvalue = ix2.value();
  } else {
    _ix2 = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool one_turn_lvalue;
  auto *_one_turn{&one_turn_lvalue};
  if (one_turn.has_value()) {
    one_turn_lvalue = one_turn.value();
  } else {
    _one_turn = nullptr;
  }
  fortran_transfer_matrix_calc(/* void* */ lat.get_fortran_ptr(),
                               /* Bmad::array_descriptor_t& */ _xfer_mat_desc,
                               /* Bmad::array_descriptor_t& */ _xfer_vec_desc,
                               /* int* */ _ix1,
                               /* int* */ _ix2,
                               /* int* */ _ix_branch,
                               /* bool* */ _one_turn);
  vec_to_matrix(_xfer_mat_vec, xfer_mat);
}
EleStruct Bmad::transfer_twiss(EleStruct &ele_in, std::optional<bool> reverse) {
  EleStruct _ele_out;
  bool reverse_lvalue;
  auto *_reverse{&reverse_lvalue};
  if (reverse.has_value()) {
    reverse_lvalue = reverse.value();
  } else {
    _reverse = nullptr;
  }
  fortran_transfer_twiss(/* void* */ ele_in.get_fortran_ptr(),
                         /* void* */ _ele_out.get_fortran_ptr(),
                         /* bool* */ _reverse);
  return std::move(_ele_out);
}
std::optional<WakeStruct> Bmad::transfer_wake(WakeStruct &wake_in) {
  auto _wake_in = &wake_in; // input, required, pointer
  void *_wake_out;
  fortran_transfer_wake(/* void* */ &wake_in, /* void* */ &_wake_out);
  return std::move((_wake_out ? std::make_optional<WakeStruct>(_wake_out) : std::nullopt));
}
void Bmad::truncate_complex_taylor_to_order(
    ComplexTaylorStructArray1D complex_taylor_in,
    int order,
    ComplexTaylorStructArray1D complex_taylor_out
) {
  // complex_taylor_in: ComplexTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_in_desc;
  _complex_taylor_in_desc.rank = 1;
  _complex_taylor_in_desc.data_ptr = complex_taylor_in.data();
  _complex_taylor_in_desc.dims[0] = complex_taylor_in.size();
  _complex_taylor_in_desc.strides[0] = 1;
  // complex_taylor_out: ComplexTaylorStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_out_desc;
  _complex_taylor_out_desc.rank = 1;
  _complex_taylor_out_desc.data_ptr = complex_taylor_out.data();
  _complex_taylor_out_desc.dims[0] = complex_taylor_out.size();
  _complex_taylor_out_desc.strides[0] = 1;
  fortran_truncate_complex_taylor_to_order(
      /* Bmad::array_descriptor_t& */ _complex_taylor_in_desc,
      /* int& */ order,
      /* Bmad::array_descriptor_t& */ _complex_taylor_out_desc
  );
}
Bmad::Twiss1Propagate Bmad::twiss1_propagate(
    TwissStruct &twiss1,
    FixedArray2D<Real, 2, 2> mat2,
    int ele_key,
    double length
) {
  // mat2: in NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _mat2_desc;
  _mat2_desc.rank = 2;
  double _mat2_vec[2 * 2];
  _mat2_desc.data_ptr = _mat2_vec;
  _mat2_desc.dims[0] = 2;
  _mat2_desc.dims[1] = 2;
  matrix_to_vec(mat2, _mat2_vec);
  TwissStruct _twiss2;
  bool _err{};
  fortran_twiss1_propagate(/* void* */ twiss1.get_fortran_ptr(),
                           /* Bmad::array_descriptor_t& */ _mat2_desc,
                           /* int& */ ele_key,
                           /* double& */ length,
                           /* void* */ _twiss2.get_fortran_ptr(),
                           /* bool& */ _err);
  return Twiss1Propagate{std::move(_twiss2), _err};
}
FixedArray1D<Real, 3>
Bmad::twiss3_at_start(LatStruct &lat, bool err_flag, std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  // tune3: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _tune3_desc;
  _tune3_desc.rank = 1;
  FixedArray1D<Real, 3> _tune3;
  _tune3_desc.data_ptr = _tune3.data();
  _tune3_desc.dims[0] = _tune3.size();
  fortran_twiss3_at_start(/* void* */ lat.get_fortran_ptr(),
                          /* bool& */ err_flag,
                          /* int* */ _ix_branch,
                          /* Bmad::array_descriptor_t& */ _tune3_desc);
  return _tune3;
}
void Bmad::twiss3_from_twiss2(EleStruct &ele) {
  fortran_twiss3_from_twiss2(/* void* */ ele.get_fortran_ptr());
}
void Bmad::twiss3_propagate1(EleStruct &ele1, EleStruct &ele2, bool err_flag) {
  fortran_twiss3_propagate1(/* void* */ ele1.get_fortran_ptr(),
                            /* void* */ ele2.get_fortran_ptr(),
                            /* bool& */ err_flag);
}
void Bmad::twiss3_propagate_all(LatStruct &lat, std::optional<int> ix_branch) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  fortran_twiss3_propagate_all(/* void* */ lat.get_fortran_ptr(), /* int* */ _ix_branch);
}
int Bmad::twiss_and_track(
    LatStruct &lat,
    CoordArrayStructAlloc1D orb_array,
    std::optional<bool> print_err,
    std::optional<bool> calc_chrom,
    std::optional<bool> use_particle_start
) {
  // intent=inout allocatable type array
  int _status{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool calc_chrom_lvalue;
  auto *_calc_chrom{&calc_chrom_lvalue};
  if (calc_chrom.has_value()) {
    calc_chrom_lvalue = calc_chrom.value();
  } else {
    _calc_chrom = nullptr;
  }
  bool use_particle_start_lvalue;
  auto *_use_particle_start{&use_particle_start_lvalue};
  if (use_particle_start.has_value()) {
    use_particle_start_lvalue = use_particle_start.value();
  } else {
    _use_particle_start = nullptr;
  }
  fortran_twiss_and_track_all(/* void* */ lat.get_fortran_ptr(),
                              /* void* */ orb_array.get_fortran_ptr(),
                              /* int& */ _status,
                              /* bool* */ _print_err,
                              /* bool* */ _calc_chrom,
                              /* bool* */ _use_particle_start);
  return _status;
}
bool Bmad::twiss_and_track_at_s(
    LatStruct &lat,
    double s,
    optional_ref<EleStruct> ele_at_s,
    std::optional<CoordStructArray1D> orb,
    optional_ref<CoordStruct> orb_at_s,
    std::optional<int> ix_branch,
    std::optional<bool> use_last,
    std::optional<bool> compute_floor_coords
) {
  auto *_ele_at_s =
      ele_at_s.has_value() ? ele_at_s->get().get_fortran_ptr() : nullptr; // input, optional
  // orb: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orb_desc;
  _orb_desc.rank = 1;
  if (orb) {
    _orb_desc.data_ptr = orb->data();
    _orb_desc.dims[0] = orb->size();
  } else {
    _orb_desc.data_ptr = nullptr;
    _orb_desc.dims[0] = 0;
  }
  _orb_desc.strides[0] = 1;
  auto *_orb_at_s =
      orb_at_s.has_value() ? orb_at_s->get().get_fortran_ptr() : nullptr; // input, optional
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err{};
  bool use_last_lvalue;
  auto *_use_last{&use_last_lvalue};
  if (use_last.has_value()) {
    use_last_lvalue = use_last.value();
  } else {
    _use_last = nullptr;
  }
  bool compute_floor_coords_lvalue;
  auto *_compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  fortran_twiss_and_track_at_s(/* void* */ lat.get_fortran_ptr(),
                               /* double& */ s,
                               /* void* */ _ele_at_s,
                               /* Bmad::array_descriptor_t& */ _orb_desc,
                               /* void* */ _orb_at_s,
                               /* int* */ _ix_branch,
                               /* bool& */ _err,
                               /* bool* */ _use_last,
                               /* bool* */ _compute_floor_coords);
  return _err;
}
int Bmad::twiss_and_track(
    LatStruct &lat,
    CoordStructAlloc1D orb,
    std::optional<int> ix_branch,
    std::optional<bool> print_err,
    std::optional<bool> calc_chrom,
    optional_ref<CoordStruct> orb_start,
    std::optional<bool> use_particle_start
) {
  // intent=inout allocatable type array
  int _status{};
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool calc_chrom_lvalue;
  auto *_calc_chrom{&calc_chrom_lvalue};
  if (calc_chrom.has_value()) {
    calc_chrom_lvalue = calc_chrom.value();
  } else {
    _calc_chrom = nullptr;
  }
  auto *_orb_start =
      orb_start.has_value() ? orb_start->get().get_fortran_ptr() : nullptr; // input, optional
  bool use_particle_start_lvalue;
  auto *_use_particle_start{&use_particle_start_lvalue};
  if (use_particle_start.has_value()) {
    use_particle_start_lvalue = use_particle_start.value();
  } else {
    _use_particle_start = nullptr;
  }
  fortran_twiss_and_track_branch(/* void* */ lat.get_fortran_ptr(),
                                 /* void* */ orb.get_fortran_ptr(),
                                 /* int& */ _status,
                                 /* int* */ _ix_branch,
                                 /* bool* */ _print_err,
                                 /* bool* */ _calc_chrom,
                                 /* void* */ _orb_start,
                                 /* bool* */ _use_particle_start);
  return _status;
}
Bmad::TwissAndTrackFromSToS Bmad::twiss_and_track_from_s_to_s(
    BranchStruct &branch,
    CoordStruct &orbit_start,
    double s_end,
    optional_ref<EleStruct> ele_start,
    std::optional<bool> compute_floor_coords,
    std::optional<bool> compute_twiss
) {
  CoordStruct _orbit_end;
  auto *_ele_start =
      ele_start.has_value() ? ele_start->get().get_fortran_ptr() : nullptr; // input, optional
  EleStruct _ele_end;
  bool _err{};
  bool compute_floor_coords_lvalue;
  auto *_compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  bool compute_twiss_lvalue;
  auto *_compute_twiss{&compute_twiss_lvalue};
  if (compute_twiss.has_value()) {
    compute_twiss_lvalue = compute_twiss.value();
  } else {
    _compute_twiss = nullptr;
  }
  fortran_twiss_and_track_from_s_to_s(/* void* */ branch.get_fortran_ptr(),
                                      /* void* */ orbit_start.get_fortran_ptr(),
                                      /* double& */ s_end,
                                      /* void* */ _orbit_end.get_fortran_ptr(),
                                      /* void* */ _ele_start,
                                      /* void* */ _ele_end.get_fortran_ptr(),
                                      /* bool& */ _err,
                                      /* bool* */ _compute_floor_coords,
                                      /* bool* */ _compute_twiss);
  return TwissAndTrackFromSToS{std::move(_orbit_end), std::move(_ele_end), _err};
}
Bmad::TwissAndTrackIntraEle Bmad::twiss_and_track_intra_ele(
    EleStruct &ele,
    LatParamStruct &param,
    double l_start,
    double l_end,
    bool track_upstream_end,
    bool track_downstream_end,
    optional_ref<CoordStruct> orbit_start,
    optional_ref<EleStruct> ele_start,
    optional_ref<EleStruct> ele_end,
    std::optional<bool> compute_floor_coords,
    std::optional<bool> compute_twiss,
    std::optional<bool> reuse_ele_end
) {
  auto *_orbit_start =
      orbit_start.has_value() ? orbit_start->get().get_fortran_ptr() : nullptr; // input, optional
  CoordStruct _orbit_end;
  auto *_ele_start =
      ele_start.has_value() ? ele_start->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_ele_end =
      ele_end.has_value() ? ele_end->get().get_fortran_ptr() : nullptr; // input, optional
  bool _err{};
  bool compute_floor_coords_lvalue;
  auto *_compute_floor_coords{&compute_floor_coords_lvalue};
  if (compute_floor_coords.has_value()) {
    compute_floor_coords_lvalue = compute_floor_coords.value();
  } else {
    _compute_floor_coords = nullptr;
  }
  bool compute_twiss_lvalue;
  auto *_compute_twiss{&compute_twiss_lvalue};
  if (compute_twiss.has_value()) {
    compute_twiss_lvalue = compute_twiss.value();
  } else {
    _compute_twiss = nullptr;
  }
  bool reuse_ele_end_lvalue;
  auto *_reuse_ele_end{&reuse_ele_end_lvalue};
  if (reuse_ele_end.has_value()) {
    reuse_ele_end_lvalue = reuse_ele_end.value();
  } else {
    _reuse_ele_end = nullptr;
  }
  fortran_twiss_and_track_intra_ele(/* void* */ ele.get_fortran_ptr(),
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
Bmad::TwissAtElement Bmad::twiss_at_element(EleStruct &ele) {
  EleStruct _start;
  EleStruct _end;
  EleStruct _average;
  fortran_twiss_at_element(/* void* */ ele.get_fortran_ptr(),
                           /* void* */ _start.get_fortran_ptr(),
                           /* void* */ _end.get_fortran_ptr(),
                           /* void* */ _average.get_fortran_ptr());
  return TwissAtElement{std::move(_start), std::move(_end), std::move(_average)};
}
int Bmad::twiss_at_start(
    LatStruct &lat,
    std::optional<int> ix_branch,
    std::optional<bool> type_out
) {
  int _status{};
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool type_out_lvalue;
  auto *_type_out{&type_out_lvalue};
  if (type_out.has_value()) {
    type_out_lvalue = type_out.value();
  } else {
    _type_out = nullptr;
  }
  fortran_twiss_at_start(/* void* */ lat.get_fortran_ptr(),
                         /* int& */ _status,
                         /* int* */ _ix_branch,
                         /* bool* */ _type_out);
  return _status;
}
Bmad::TwissFromTracking Bmad::twiss_from_tracking(
    LatStruct &lat,
    CoordStruct &ref_orb0,
    std::optional<FArray1D<Real>> d_orb
) {
  double _symp_err{};
  bool _err_flag{};
  // d_orb: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _d_orb_desc;
  _d_orb_desc.rank = 1;
  if (d_orb.has_value()) {
    _d_orb_desc.data_ptr = d_orb->data();
    _d_orb_desc.dims[0] = d_orb->size();
  } else {
    _d_orb_desc.data_ptr = nullptr;
    _d_orb_desc.dims[0] = 0;
  }
  fortran_twiss_from_tracking(/* void* */ lat.get_fortran_ptr(),
                              /* void* */ ref_orb0.get_fortran_ptr(),
                              /* double& */ _symp_err,
                              /* bool& */ _err_flag,
                              /* Bmad::array_descriptor_t& */ _d_orb_desc);
  return TwissFromTracking{_symp_err, _err_flag};
}
bool Bmad::twiss_propagate1(EleStruct &ele1, EleStruct &ele2, std::optional<bool> forward) {
  bool _err_flag{};
  bool forward_lvalue;
  auto *_forward{&forward_lvalue};
  if (forward.has_value()) {
    forward_lvalue = forward.value();
  } else {
    _forward = nullptr;
  }
  fortran_twiss_propagate1(/* void* */ ele1.get_fortran_ptr(),
                           /* void* */ ele2.get_fortran_ptr(),
                           /* bool& */ _err_flag,
                           /* bool* */ _forward);
  return _err_flag;
}
bool Bmad::twiss_propagate_all(
    LatStruct &lat,
    std::optional<int> ix_branch,
    std::optional<int> ie_start,
    std::optional<int> ie_end
) {
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool _err_flag{};
  int ie_start_lvalue;
  auto *_ie_start{&ie_start_lvalue};
  if (ie_start.has_value()) {
    ie_start_lvalue = ie_start.value();
  } else {
    _ie_start = nullptr;
  }
  int ie_end_lvalue;
  auto *_ie_end{&ie_end_lvalue};
  if (ie_end.has_value()) {
    ie_end_lvalue = ie_end.value();
  } else {
    _ie_end = nullptr;
  }
  fortran_twiss_propagate_all(/* void* */ lat.get_fortran_ptr(),
                              /* int* */ _ix_branch,
                              /* bool& */ _err_flag,
                              /* int* */ _ie_start,
                              /* int* */ _ie_end);
  return _err_flag;
}
FixedArray2D<Real, 2, 2> Bmad::twiss_to_1_turn_mat(TwissStruct &twiss, double phi) {
  // mat2: out NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _mat2_desc;
  _mat2_desc.rank = 2;
  FixedArray2D<Real, 2, 2> mat2;
  double _mat2_vec[2 * 2];
  _mat2_desc.data_ptr = _mat2_vec;
  _mat2_desc.dims[0] = 2;
  _mat2_desc.dims[1] = 2;
  fortran_twiss_to_1_turn_mat(/* void* */ twiss.get_fortran_ptr(),
                              /* double& */ phi,
                              /* Bmad::array_descriptor_t& */ _mat2_desc);
  vec_to_matrix(_mat2_vec, mat2);
  return mat2;
}
Bmad::TypeComplexTaylors Bmad::type_complex_taylors(
    ComplexTaylorStructArray1D complex_taylor,
    std::optional<int> max_order,
    std::optional<int> file_id,
    std::optional<std::string> out_type,
    std::optional<bool> clean
) {
  // complex_taylor: ComplexTaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _complex_taylor_desc;
  _complex_taylor_desc.rank = 1;
  _complex_taylor_desc.data_ptr = complex_taylor.data();
  _complex_taylor_desc.dims[0] = complex_taylor.size();
  _complex_taylor_desc.strides[0] = 1;
  int max_order_lvalue;
  auto *_max_order{&max_order_lvalue};
  if (max_order.has_value()) {
    max_order_lvalue = max_order.value();
  } else {
    _max_order = nullptr;
  }
  // intent=out character array container
  auto lines{CharacterAlloc1D()};
  int _n_lines{};
  int file_id_lvalue;
  auto *_file_id{&file_id_lvalue};
  if (file_id.has_value()) {
    file_id_lvalue = file_id.value();
  } else {
    _file_id = nullptr;
  }
  const char *_out_type = out_type.has_value() ? out_type->c_str() : nullptr;
  bool clean_lvalue;
  auto *_clean{&clean_lvalue};
  if (clean.has_value()) {
    clean_lvalue = clean.value();
  } else {
    _clean = nullptr;
  }
  fortran_type_complex_taylors(
      /* Bmad::array_descriptor_t& */ _complex_taylor_desc,
      /* int* */ _max_order,
      /* void* */ lines.get_fortran_ptr(),
      /* int& */ _n_lines,
      /* int* */ _file_id,
      /* const char* */ _out_type,
      /* bool* */ _clean
  );
  return TypeComplexTaylors{std::move(lines), _n_lines};
}
void Bmad::type_coord(CoordStruct &coord) {
  fortran_type_coord(/* void* */ coord.get_fortran_ptr());
}
Bmad::TypeEle Bmad::type_ele(
    EleStruct &ele,
    std::optional<bool> type_zero_attrib,
    std::optional<int> type_mat6,
    std::optional<bool> type_taylor,
    std::optional<int> twiss_out,
    std::optional<int> type_control,
    std::optional<bool> type_wake,
    std::optional<bool> type_floor_coords,
    std::optional<int> type_field,
    std::optional<bool> type_wall,
    std::optional<bool> type_rad_kick,
    std::optional<bool> type_internal
) {
  bool type_zero_attrib_lvalue;
  auto *_type_zero_attrib{&type_zero_attrib_lvalue};
  if (type_zero_attrib.has_value()) {
    type_zero_attrib_lvalue = type_zero_attrib.value();
  } else {
    _type_zero_attrib = nullptr;
  }
  int type_mat6_lvalue;
  auto *_type_mat6{&type_mat6_lvalue};
  if (type_mat6.has_value()) {
    type_mat6_lvalue = type_mat6.value();
  } else {
    _type_mat6 = nullptr;
  }
  bool type_taylor_lvalue;
  auto *_type_taylor{&type_taylor_lvalue};
  if (type_taylor.has_value()) {
    type_taylor_lvalue = type_taylor.value();
  } else {
    _type_taylor = nullptr;
  }
  int twiss_out_lvalue;
  auto *_twiss_out{&twiss_out_lvalue};
  if (twiss_out.has_value()) {
    twiss_out_lvalue = twiss_out.value();
  } else {
    _twiss_out = nullptr;
  }
  int type_control_lvalue;
  auto *_type_control{&type_control_lvalue};
  if (type_control.has_value()) {
    type_control_lvalue = type_control.value();
  } else {
    _type_control = nullptr;
  }
  bool type_wake_lvalue;
  auto *_type_wake{&type_wake_lvalue};
  if (type_wake.has_value()) {
    type_wake_lvalue = type_wake.value();
  } else {
    _type_wake = nullptr;
  }
  bool type_floor_coords_lvalue;
  auto *_type_floor_coords{&type_floor_coords_lvalue};
  if (type_floor_coords.has_value()) {
    type_floor_coords_lvalue = type_floor_coords.value();
  } else {
    _type_floor_coords = nullptr;
  }
  int type_field_lvalue;
  auto *_type_field{&type_field_lvalue};
  if (type_field.has_value()) {
    type_field_lvalue = type_field.value();
  } else {
    _type_field = nullptr;
  }
  bool type_wall_lvalue;
  auto *_type_wall{&type_wall_lvalue};
  if (type_wall.has_value()) {
    type_wall_lvalue = type_wall.value();
  } else {
    _type_wall = nullptr;
  }
  bool type_rad_kick_lvalue;
  auto *_type_rad_kick{&type_rad_kick_lvalue};
  if (type_rad_kick.has_value()) {
    type_rad_kick_lvalue = type_rad_kick.value();
  } else {
    _type_rad_kick = nullptr;
  }
  bool type_internal_lvalue;
  auto *_type_internal{&type_internal_lvalue};
  if (type_internal.has_value()) {
    type_internal_lvalue = type_internal.value();
  } else {
    _type_internal = nullptr;
  }
  // intent=out character array container
  auto lines{CharacterAlloc1D()};
  int _n_lines{};
  fortran_type_ele(/* void* */ ele.get_fortran_ptr(),
                   /* bool* */ _type_zero_attrib,
                   /* int* */ _type_mat6,
                   /* bool* */ _type_taylor,
                   /* int* */ _twiss_out,
                   /* int* */ _type_control,
                   /* bool* */ _type_wake,
                   /* bool* */ _type_floor_coords,
                   /* int* */ _type_field,
                   /* bool* */ _type_wall,
                   /* bool* */ _type_rad_kick,
                   /* bool* */ _type_internal,
                   /* void* */ lines.get_fortran_ptr(),
                   /* int& */ _n_lines);
  return TypeEle{std::move(lines), _n_lines};
}
void Bmad::type_end_stuff(
    CharacterAlloc1D &li,
    int nl,
    optional_ref<CharacterAlloc1D> lines,
    std::optional<int> n_lines
) {
  // intent=inout character array container
  // intent=inout character array container
  auto *_lines = lines.has_value() ? lines->get().get_fortran_ptr() : nullptr; // input, optional
  int n_lines_lvalue;
  auto *_n_lines{&n_lines_lvalue};
  if (n_lines.has_value()) {
    n_lines_lvalue = n_lines.value();
  } else {
    _n_lines = nullptr;
  }
  fortran_type_end_stuff(/* void* */ li.get_fortran_ptr(),
                         /* int& */ nl,
                         /* void* */ _lines,
                         /* int* */ _n_lines);
}
void Bmad::type_expression_tree(ExpressionTreeStruct &tree, std::optional<int> indent) {
  int indent_lvalue;
  auto *_indent{&indent_lvalue};
  if (indent.has_value()) {
    indent_lvalue = indent.value();
  } else {
    _indent = nullptr;
  }
  fortran_type_expression_tree(/* void* */ tree.get_fortran_ptr(), /* int* */ _indent);
}
Bmad::TypePtcFibre Bmad::type_ptc_fibre(Fibre &ptc_fibre, std::optional<bool> print_coords) {
  auto _ptc_fibre = &ptc_fibre; // input, required, pointer
  bool print_coords_lvalue;
  auto *_print_coords{&print_coords_lvalue};
  if (print_coords.has_value()) {
    print_coords_lvalue = print_coords.value();
  } else {
    _print_coords = nullptr;
  }
  // intent=out character array container
  auto lines{CharacterAlloc1D()};
  int _n_lines{};
  fortran_type_ptc_fibre(/* void* */ &ptc_fibre,
                         /* bool* */ _print_coords,
                         /* void* */ lines.get_fortran_ptr(),
                         /* int& */ _n_lines);
  return TypePtcFibre{std::move(lines), _n_lines};
}
void Bmad::type_ptc_layout(Layout &lay) {
  fortran_type_ptc_layout(/* void* */ lay.get_fortran_ptr());
}
void Bmad::type_taylors(
    TaylorStructArray1D bmad_taylor,
    std::optional<int> max_order,
    optional_ref<CharacterAlloc1D> lines,
    optional_ref<int> n_lines,
    std::optional<int> file_id,
    std::optional<std::string> out_style,
    std::optional<bool> clean,
    std::optional<std::string> out_var_suffix,
    std::optional<bool> append
) {
  // bmad_taylor: TaylorStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _bmad_taylor_desc;
  _bmad_taylor_desc.rank = 1;
  _bmad_taylor_desc.data_ptr = bmad_taylor.data();
  _bmad_taylor_desc.dims[0] = bmad_taylor.size();
  _bmad_taylor_desc.strides[0] = 1;
  int max_order_lvalue;
  auto *_max_order{&max_order_lvalue};
  if (max_order.has_value()) {
    max_order_lvalue = max_order.value();
  } else {
    _max_order = nullptr;
  }
  // intent=inout character array container
  auto *_lines = lines.has_value() ? lines->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_n_lines = n_lines.has_value() ? &n_lines->get() : nullptr; // inout, optional
  int file_id_lvalue;
  auto *_file_id{&file_id_lvalue};
  if (file_id.has_value()) {
    file_id_lvalue = file_id.value();
  } else {
    _file_id = nullptr;
  }
  const char *_out_style = out_style.has_value() ? out_style->c_str() : nullptr;
  bool clean_lvalue;
  auto *_clean{&clean_lvalue};
  if (clean.has_value()) {
    clean_lvalue = clean.value();
  } else {
    _clean = nullptr;
  }
  const char *_out_var_suffix = out_var_suffix.has_value() ? out_var_suffix->c_str() : nullptr;
  bool append_lvalue;
  auto *_append{&append_lvalue};
  if (append.has_value()) {
    append_lvalue = append.value();
  } else {
    _append = nullptr;
  }
  fortran_type_taylors(
      /* Bmad::array_descriptor_t& */ _bmad_taylor_desc,
      /* int* */ _max_order,
      /* void* */ _lines,
      /* int* */ _n_lines,
      /* int* */ _file_id,
      /* const char* */ _out_style,
      /* bool* */ _clean,
      /* const char* */ _out_var_suffix,
      /* bool* */ _append
  );
}
int Bmad::type_twiss(
    EleStruct &ele,
    std::optional<int> frequency_units,
    std::optional<bool> compact_format,
    optional_ref<CharacterAlloc1D> lines
) {
  int frequency_units_lvalue;
  auto *_frequency_units{&frequency_units_lvalue};
  if (frequency_units.has_value()) {
    frequency_units_lvalue = frequency_units.value();
  } else {
    _frequency_units = nullptr;
  }
  bool compact_format_lvalue;
  auto *_compact_format{&compact_format_lvalue};
  if (compact_format.has_value()) {
    compact_format_lvalue = compact_format.value();
  } else {
    _compact_format = nullptr;
  }
  // intent=inout character array container
  auto *_lines = lines.has_value() ? lines->get().get_fortran_ptr() : nullptr; // input, optional
  int _n_lines{};
  fortran_type_twiss(/* void* */ ele.get_fortran_ptr(),
                     /* int* */ _frequency_units,
                     /* bool* */ _compact_format,
                     /* void* */ _lines,
                     /* int& */ _n_lines);
  return _n_lines;
}
void Bmad::update_ele_from_fibre(EleStruct &ele) {
  fortran_update_ele_from_fibre(/* void* */ ele.get_fortran_ptr());
}
bool Bmad::update_fibre_from_ele(EleStruct &ele) {
  bool _survey_needed{};
  fortran_update_fibre_from_ele(/* void* */ ele.get_fortran_ptr(), /* bool& */ _survey_needed);
  return _survey_needed;
}
void Bmad::update_floor_angles(
    FloorPositionStruct &floor,
    optional_ref<FloorPositionStruct> floor0
) {
  auto *_floor0 = floor0.has_value() ? floor0->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_update_floor_angles(/* void* */ floor.get_fortran_ptr(), /* void* */ _floor0);
}
bool Bmad::valid_field_calc(EleStruct &ele, int field_calc) {
  bool _is_valid{};
  fortran_valid_field_calc(/* void* */ ele.get_fortran_ptr(),
                           /* int& */ field_calc,
                           /* bool& */ _is_valid);
  return _is_valid;
}
bool Bmad::valid_fringe_type(EleStruct &ele, int fringe_type) {
  bool _is_valid{};
  fortran_valid_fringe_type(/* void* */ ele.get_fortran_ptr(),
                            /* int& */ fringe_type,
                            /* bool& */ _is_valid);
  return _is_valid;
}
bool Bmad::valid_mat6_calc_method(EleStruct &ele, int species, int mat6_calc_method) {
  bool _is_valid{};
  fortran_valid_mat6_calc_method(/* void* */ ele.get_fortran_ptr(),
                                 /* int& */ species,
                                 /* int& */ mat6_calc_method,
                                 /* bool& */ _is_valid);
  return _is_valid;
}
bool Bmad::valid_spin_tracking_method(EleStruct &ele, int spin_tracking_method) {
  bool _is_valid{};
  fortran_valid_spin_tracking_method(/* void* */ ele.get_fortran_ptr(),
                                     /* int& */ spin_tracking_method,
                                     /* bool& */ _is_valid);
  return _is_valid;
}
bool Bmad::valid_tracking_method(EleStruct &ele, int species, int tracking_method) {
  bool _is_valid{};
  fortran_valid_tracking_method(/* void* */ ele.get_fortran_ptr(),
                                /* int& */ species,
                                /* int& */ tracking_method,
                                /* bool& */ _is_valid);
  return _is_valid;
}
Bmad::ValueOfAttribute Bmad::value_of_attribute(
    EleStruct &ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<double> err_value
) {
  auto _attrib_name = attrib_name.c_str();
  bool _err_flag{};
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  double err_value_lvalue;
  auto *_err_value{&err_value_lvalue};
  if (err_value.has_value()) {
    err_value_lvalue = err_value.value();
  } else {
    _err_value = nullptr;
  }
  double _value{};
  fortran_value_of_attribute(/* void* */ ele.get_fortran_ptr(),
                             /* const char* */ _attrib_name,
                             /* bool& */ _err_flag,
                             /* bool* */ _err_print_flag,
                             /* double* */ _err_value,
                             /* double& */ _value);
  return ValueOfAttribute{_err_flag, _value};
}
void Bmad::value_to_line(
    std::string line,
    double value,
    std::string str,
    std::string typ,
    std::optional<bool> ignore_if_zero,
    std::optional<bool> use_comma
) {
  auto _line = line.c_str();
  auto _str = str.c_str();
  auto _typ = typ.c_str();
  bool ignore_if_zero_lvalue;
  auto *_ignore_if_zero{&ignore_if_zero_lvalue};
  if (ignore_if_zero.has_value()) {
    ignore_if_zero_lvalue = ignore_if_zero.value();
  } else {
    _ignore_if_zero = nullptr;
  }
  bool use_comma_lvalue;
  auto *_use_comma{&use_comma_lvalue};
  if (use_comma.has_value()) {
    use_comma_lvalue = use_comma.value();
  } else {
    _use_comma = nullptr;
  }
  fortran_value_to_line(
      /* const char* */ _line,
      /* double& */ value,
      /* const char* */ _str,
      /* const char* */ _typ,
      /* bool* */ _ignore_if_zero,
      /* bool* */ _use_comma
  );
}
SpinPolarStruct Bmad::vec_to_polar(FixedArray1D<Real, 3> vec, std::optional<double> phase) {
  // vec: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  double phase_lvalue;
  auto *_phase{&phase_lvalue};
  if (phase.has_value()) {
    phase_lvalue = phase.value();
  } else {
    _phase = nullptr;
  }
  SpinPolarStruct _polar;
  fortran_vec_to_polar(
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* double* */ _phase,
      /* void* */ _polar.get_fortran_ptr()
  );
  return std::move(_polar);
}
FixedArray1D<Complex, 2>
Bmad::vec_to_spinor(FixedArray1D<Real, 3> vec, std::optional<double> phase) {
  // vec: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  double phase_lvalue;
  auto *_phase{&phase_lvalue};
  if (phase.has_value()) {
    phase_lvalue = phase.value();
  } else {
    _phase = nullptr;
  }
  // spinor: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _spinor_desc;
  _spinor_desc.rank = 1;
  FixedArray1D<Complex, 2> _spinor;
  _spinor_desc.data_ptr = _spinor.data();
  _spinor_desc.dims[0] = _spinor.size();
  fortran_vec_to_spinor(
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* double* */ _phase,
      /* Bmad::array_descriptor_t& */ _spinor_desc
  );
  return _spinor;
}
bool Bmad::verify_valid_name(
    std::string name,
    int ix_name,
    std::optional<bool> pure_name,
    std::optional<bool> include_wild
) {
  auto _name = name.c_str();
  bool pure_name_lvalue;
  auto *_pure_name{&pure_name_lvalue};
  if (pure_name.has_value()) {
    pure_name_lvalue = pure_name.value();
  } else {
    _pure_name = nullptr;
  }
  bool include_wild_lvalue;
  auto *_include_wild{&include_wild_lvalue};
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
      /* bool& */ _is_valid
  );
  return _is_valid;
}
FixedArray2D<Real, 3, 3> Bmad::w_mat_for_bend_angle(
    double angle,
    double ref_tilt,
    std::optional<FixedArray1D<Real, 3>> r_vec
) {
  // r_vec: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r_vec_desc;
  _r_vec_desc.rank = 1;
  if (r_vec.has_value()) {
    _r_vec_desc.data_ptr = r_vec->data();
    _r_vec_desc.dims[0] = r_vec->size();
  } else {
    _r_vec_desc.data_ptr = nullptr;
    _r_vec_desc.dims[0] = 0;
  }
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_w_mat_for_bend_angle(
      /* double& */ angle,
      /* double& */ ref_tilt,
      /* Bmad::array_descriptor_t& */ _r_vec_desc,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray2D<Real, 3, 3> Bmad::w_mat_for_tilt(double tilt, std::optional<bool> return_inverse) {
  bool return_inverse_lvalue;
  auto *_return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_w_mat_for_tilt(
      /* double& */ tilt,
      /* bool* */ _return_inverse,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray2D<Real, 3, 3>
Bmad::w_mat_for_x_pitch(double x_pitch, std::optional<bool> return_inverse) {
  bool return_inverse_lvalue;
  auto *_return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_w_mat_for_x_pitch(
      /* double& */ x_pitch,
      /* bool* */ _return_inverse,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray2D<Real, 3, 3>
Bmad::w_mat_for_y_pitch(double y_pitch, std::optional<bool> return_inverse) {
  bool return_inverse_lvalue;
  auto *_return_inverse{&return_inverse_lvalue};
  if (return_inverse.has_value()) {
    return_inverse_lvalue = return_inverse.value();
  } else {
    _return_inverse = nullptr;
  }
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_w_mat_for_y_pitch(
      /* double& */ y_pitch,
      /* bool* */ _return_inverse,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
Bmad::Wall3dDRadius
Bmad::wall3d_d_radius(FArray1D<Real> &position, EleStruct &ele, std::optional<int> ix_wall) {
  // position: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _position_desc;
  _position_desc.rank = 1;
  _position_desc.data_ptr = position.data();
  _position_desc.dims[0] = position.size();
  int ix_wall_lvalue;
  auto *_ix_wall{&ix_wall_lvalue};
  if (ix_wall.has_value()) {
    ix_wall_lvalue = ix_wall.value();
  } else {
    _ix_wall = nullptr;
  }
  // perp: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _perp_desc;
  _perp_desc.rank = 1;
  FixedArray1D<Real, 3> _perp;
  _perp_desc.data_ptr = _perp.data();
  _perp_desc.dims[0] = _perp.size();
  int _ix_section{};
  bool _no_wall_here{};
  // origin: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _origin_desc;
  _origin_desc.rank = 1;
  FixedArray1D<Real, 3> _origin;
  _origin_desc.data_ptr = _origin.data();
  _origin_desc.dims[0] = _origin.size();
  double _radius_wall{};
  bool _err_flag{};
  double _d_radius{};
  fortran_wall3d_d_radius(
      /* Bmad::array_descriptor_t& */ _position_desc,
      /* void* */ ele.get_fortran_ptr(),
      /* int* */ _ix_wall,
      /* Bmad::array_descriptor_t& */ _perp_desc,
      /* int& */ _ix_section,
      /* bool& */ _no_wall_here,
      /* Bmad::array_descriptor_t& */ _origin_desc,
      /* double& */ _radius_wall,
      /* bool& */ _err_flag,
      /* double& */ _d_radius
  );
  return Wall3dDRadius{
      _perp,
      _ix_section,
      _no_wall_here,
      _origin,
      _radius_wall,
      _err_flag,
      _d_radius
  };
}
bool Bmad::wall3d_initializer(Wall3dStruct &wall3d) {
  bool _err{};
  fortran_wall3d_initializer(/* void* */ wall3d.get_fortran_ptr(), /* bool& */ _err);
  return _err;
}
bool Bmad::wall3d_section_initializer(Wall3dSectionStruct &section) {
  bool _err{};
  fortran_wall3d_section_initializer(/* void* */ section.get_fortran_ptr(), /* bool& */ _err);
  return _err;
}
FixedArray1D<Real, 6> Bmad::wall3d_to_position(CoordStruct &orbit, EleStruct &ele) {
  // position: out NOT (CppWrapperGeneralArgumentArray) (['6'])
  Bmad::array_descriptor_t _position_desc;
  _position_desc.rank = 1;
  FixedArray1D<Real, 6> _position;
  _position_desc.data_ptr = _position.data();
  _position_desc.dims[0] = _position.size();
  fortran_wall3d_to_position(/* void* */ orbit.get_fortran_ptr(),
                             /* void* */ ele.get_fortran_ptr(),
                             /* Bmad::array_descriptor_t& */ _position_desc);
  return _position;
}
void Bmad::word_to_value(
    std::string word,
    LatStruct &lat,
    double value,
    bool err_flag,
    optional_ref<EleStruct> ele
) {
  auto _word = word.c_str();
  auto *_ele = ele.has_value() ? ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_word_to_value(
      /* const char* */ _word,
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ value,
      /* bool& */ err_flag,
      /* void* */ _ele
  );
}
void Bmad::write_ascii_beam_file(
    std::string file_name,
    BeamStruct &beam,
    std::optional<bool> new_file,
    std::optional<bool> alive_only
) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto *_new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  bool alive_only_lvalue;
  auto *_alive_only{&alive_only_lvalue};
  if (alive_only.has_value()) {
    alive_only_lvalue = alive_only.value();
  } else {
    _alive_only = nullptr;
  }
  fortran_write_ascii_beam_file(
      /* const char* */ _file_name,
      /* void* */ beam.get_fortran_ptr(),
      /* bool* */ _new_file,
      /* bool* */ _alive_only
  );
}
void Bmad::write_astra_bend(
    int iu,
    double strength,
    int id,
    FixedArray1D<Real, 2> d1,
    FixedArray1D<Real, 2> d2,
    FixedArray1D<Real, 2> d3,
    FixedArray1D<Real, 2> d4
) {
  // d1: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _d1_desc;
  _d1_desc.rank = 1;
  _d1_desc.data_ptr = d1.data();
  _d1_desc.dims[0] = d1.size();
  // d2: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _d2_desc;
  _d2_desc.rank = 1;
  _d2_desc.data_ptr = d2.data();
  _d2_desc.dims[0] = d2.size();
  // d3: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _d3_desc;
  _d3_desc.rank = 1;
  _d3_desc.data_ptr = d3.data();
  _d3_desc.dims[0] = d3.size();
  // d4: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _d4_desc;
  _d4_desc.rank = 1;
  _d4_desc.data_ptr = d4.data();
  _d4_desc.dims[0] = d4.size();
  fortran_write_astra_bend(
      /* int& */ iu,
      /* double& */ strength,
      /* int& */ id,
      /* Bmad::array_descriptor_t& */ _d1_desc,
      /* Bmad::array_descriptor_t& */ _d2_desc,
      /* Bmad::array_descriptor_t& */ _d3_desc,
      /* Bmad::array_descriptor_t& */ _d4_desc
  );
}
Bmad::WriteAstraFieldGridFile
Bmad::write_astra_field_grid_file(int astra_file_unit, EleStruct &ele, std::optional<double> dz) {
  double _maxfield{};
  double dz_lvalue;
  auto *_dz{&dz_lvalue};
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
      /* bool& */ _err
  );
  return WriteAstraFieldGridFile{_maxfield, _err};
}
Bmad::WriteAstraFieldGridFile3d Bmad::write_astra_field_grid_file_3d(
    std::string base_filename,
    EleStruct &ele,
    std::optional<double> dz
) {
  auto _base_filename = base_filename.c_str();
  double _maxfield{};
  double dz_lvalue;
  auto *_dz{&dz_lvalue};
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
      /* bool& */ _err
  );
  return WriteAstraFieldGridFile3d{_maxfield, _err};
}
void Bmad::write_beam_file(
    std::string file_name,
    BeamStruct &beam,
    std::optional<bool> new_file,
    std::optional<int> file_format,
    optional_ref<LatStruct> lat,
    std::optional<bool> alive_only
) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto *_new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  int file_format_lvalue;
  auto *_file_format{&file_format_lvalue};
  if (file_format.has_value()) {
    file_format_lvalue = file_format.value();
  } else {
    _file_format = nullptr;
  }
  auto *_lat = lat.has_value() ? lat->get().get_fortran_ptr() : nullptr; // input, optional
  bool alive_only_lvalue;
  auto *_alive_only{&alive_only_lvalue};
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
      /* bool* */ _alive_only
  );
}
void Bmad::write_beam_floor_positions(
    std::string file_name,
    BeamStruct &beam,
    EleStruct &ele,
    std::optional<bool> new_file
) {
  auto _file_name = file_name.c_str();
  bool new_file_lvalue;
  auto *_new_file{&new_file_lvalue};
  if (new_file.has_value()) {
    new_file_lvalue = new_file.value();
  } else {
    _new_file = nullptr;
  }
  fortran_write_beam_floor_positions(
      /* const char* */ _file_name,
      /* void* */ beam.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _new_file
  );
}
bool Bmad::write_binary_cartesian_map(
    std::string file_name,
    EleStruct &ele,
    CartesianMapStruct &cart_map
) {
  auto _file_name = file_name.c_str();
  bool _err_flag{};
  fortran_write_binary_cartesian_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cart_map.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return _err_flag;
}
bool Bmad::write_binary_cylindrical_map(
    std::string file_name,
    EleStruct &ele,
    CylindricalMapStruct &cl_map
) {
  auto _file_name = file_name.c_str();
  bool _err_flag{};
  fortran_write_binary_cylindrical_map(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ cl_map.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return _err_flag;
}
bool Bmad::write_binary_grid_field(
    std::string file_name,
    EleStruct &ele,
    GridFieldStruct &g_field
) {
  auto _file_name = file_name.c_str();
  bool _err_flag{};
  fortran_write_binary_grid_field(
      /* const char* */ _file_name,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ g_field.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return _err_flag;
}
void Bmad::write_blender_ele(int iu, EleStruct &ele, std::optional<bool> old_format) {
  bool old_format_lvalue;
  auto *_old_format{&old_format_lvalue};
  if (old_format.has_value()) {
    old_format_lvalue = old_format.value();
  } else {
    _old_format = nullptr;
  }
  fortran_write_blender_ele(
      /* int& */ iu,
      /* void* */ ele.get_fortran_ptr(),
      /* bool* */ _old_format
  );
}
void Bmad::write_blender_lat_layout(std::string file_name, LatStruct &lat) {
  auto _file_name = file_name.c_str();
  fortran_write_blender_lat_layout(/* const char* */ _file_name, /* void* */ lat.get_fortran_ptr());
}
bool Bmad::write_bmad_lattice_file(
    std::string bmad_file,
    LatStruct &lat,
    std::optional<int> output_form,
    optional_ref<CoordStruct> orbit0
) {
  auto _bmad_file = bmad_file.c_str();
  bool _err{};
  int output_form_lvalue;
  auto *_output_form{&output_form_lvalue};
  if (output_form.has_value()) {
    output_form_lvalue = output_form.value();
  } else {
    _output_form = nullptr;
  }
  auto *_orbit0 = orbit0.has_value() ? orbit0->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_write_bmad_lattice_file(
      /* const char* */ _bmad_file,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err,
      /* int* */ _output_form,
      /* void* */ _orbit0
  );
  return _err;
}
Bmad::WriteGptFieldGridFile1d
Bmad::write_gpt_field_grid_file_1d(int gpt_file_unit, EleStruct &ele, std::optional<double> dz) {
  double _maxfield{};
  double _ref_time{};
  double dz_lvalue;
  auto *_dz{&dz_lvalue};
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
      /* bool& */ _err
  );
  return WriteGptFieldGridFile1d{_maxfield, _ref_time, _err};
}
Bmad::WriteGptFieldGridFile2d Bmad::write_gpt_field_grid_file_2d(
    int gpt_file_unit,
    EleStruct &ele,
    std::optional<double> dr,
    std::optional<double> dz,
    std::optional<double> r_max
) {
  double _maxfield{};
  double _ref_time{};
  double dr_lvalue;
  auto *_dr{&dr_lvalue};
  if (dr.has_value()) {
    dr_lvalue = dr.value();
  } else {
    _dr = nullptr;
  }
  double dz_lvalue;
  auto *_dz{&dz_lvalue};
  if (dz.has_value()) {
    dz_lvalue = dz.value();
  } else {
    _dz = nullptr;
  }
  double r_max_lvalue;
  auto *_r_max{&r_max_lvalue};
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
      /* bool& */ _err
  );
  return WriteGptFieldGridFile2d{_maxfield, _ref_time, _err};
}
Bmad::WriteGptFieldGridFile3d Bmad::write_gpt_field_grid_file_3d(
    std::string base_filename,
    EleStruct &ele,
    std::optional<double> dz
) {
  auto _base_filename = base_filename.c_str();
  double _maxfield{};
  double _ref_time{};
  double dz_lvalue;
  auto *_dz{&dz_lvalue};
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
      /* bool& */ _err
  );
  return WriteGptFieldGridFile3d{_maxfield, _ref_time, _err};
}
void Bmad::write_lat_line(
    std::string &line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split,
    std::optional<bool> ampersand_at_ends
) {
  auto _line = line.c_str(); // ptr, inout, required
  bool do_split_lvalue;
  auto *_do_split{&do_split_lvalue};
  if (do_split.has_value()) {
    do_split_lvalue = do_split.value();
  } else {
    _do_split = nullptr;
  }
  bool ampersand_at_ends_lvalue;
  auto *_ampersand_at_ends{&ampersand_at_ends_lvalue};
  if (ampersand_at_ends.has_value()) {
    ampersand_at_ends_lvalue = ampersand_at_ends.value();
  } else {
    _ampersand_at_ends = nullptr;
  }
  fortran_write_lat_line(
      /* const char* */ _line,
      /* int& */ iu,
      /* bool& */ end_is_neigh,
      /* bool* */ _do_split,
      /* bool* */ _ampersand_at_ends
  );
}
bool Bmad::write_lattice_in_elegant_format(
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch
) {
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto *_ref_orbit =
      ref_orbit.has_value() ? ref_orbit->get_fortran_ptr() : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto *_use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto *_include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto *_dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
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
      /* bool& */ _err
  );
  return _err;
}
bool Bmad::write_lattice_in_foreign_format(
    std::string out_type,
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch
) {
  auto _out_type = out_type.c_str();
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto *_ref_orbit =
      ref_orbit.has_value() ? ref_orbit->get_fortran_ptr() : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto *_use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto *_include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto *_dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
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
      /* bool& */ _err
  );
  return _err;
}
bool Bmad::write_lattice_in_mad_format(
    std::string out_type,
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit,
    std::optional<bool> use_matrix_model,
    std::optional<bool> include_apertures,
    std::optional<double> dr12_drift_max,
    std::optional<int> ix_branch
) {
  auto _out_type = out_type.c_str();
  auto _out_file_name = out_file_name.c_str();
  // intent=in allocatable type array
  auto *_ref_orbit =
      ref_orbit.has_value() ? ref_orbit->get_fortran_ptr() : nullptr; // input, optional
  bool use_matrix_model_lvalue;
  auto *_use_matrix_model{&use_matrix_model_lvalue};
  if (use_matrix_model.has_value()) {
    use_matrix_model_lvalue = use_matrix_model.value();
  } else {
    _use_matrix_model = nullptr;
  }
  bool include_apertures_lvalue;
  auto *_include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  double dr12_drift_max_lvalue;
  auto *_dr12_drift_max{&dr12_drift_max_lvalue};
  if (dr12_drift_max.has_value()) {
    dr12_drift_max_lvalue = dr12_drift_max.value();
  } else {
    _dr12_drift_max = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
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
      /* bool& */ _err
  );
  return _err;
}
Bmad::WriteLatticeInPals Bmad::write_lattice_in_pals(LatStruct &lat) {
  char _pals_file[4096];
  bool _err_flag{};
  fortran_write_lattice_in_pals(
      /* const char* */ _pals_file,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return WriteLatticeInPals{_pals_file, _err_flag};
}
void Bmad::write_lattice_in_sad_format(
    std::string out_file_name,
    LatStruct &lat,
    std::optional<bool> include_apertures,
    std::optional<int> ix_branch,
    std::optional<bool> err
) {
  auto _out_file_name = out_file_name.c_str();
  bool include_apertures_lvalue;
  auto *_include_apertures{&include_apertures_lvalue};
  if (include_apertures.has_value()) {
    include_apertures_lvalue = include_apertures.value();
  } else {
    _include_apertures = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool err_lvalue;
  auto *_err{&err_lvalue};
  if (err.has_value()) {
    err_lvalue = err.value();
  } else {
    _err = nullptr;
  }
  fortran_write_lattice_in_sad_format(
      /* const char* */ _out_file_name,
      /* void* */ lat.get_fortran_ptr(),
      /* bool* */ _include_apertures,
      /* int* */ _ix_branch,
      /* bool* */ _err
  );
}
Bmad::WriteLatticeInScibmad Bmad::write_lattice_in_scibmad(LatStruct &lat) {
  char _scibmad_file[4096];
  bool _err_flag{};
  fortran_write_lattice_in_scibmad(
      /* const char* */ _scibmad_file,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return WriteLatticeInScibmad{_scibmad_file, _err_flag};
}
void Bmad::write_line_element(std::string line, int iu, EleStruct &ele, LatStruct &lat) {
  auto _line = line.c_str();
  fortran_write_line_element(
      /* const char* */ _line,
      /* int& */ iu,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ lat.get_fortran_ptr()
  );
}
Bmad::WriteOpalFieldGridFile
Bmad::write_opal_field_grid_file(int opal_file_unit, EleStruct &ele, LatParamStruct &param) {
  double _maxfield{};
  bool _err{};
  fortran_write_opal_field_grid_file(
      /* int& */ opal_file_unit,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ param.get_fortran_ptr(),
      /* double& */ _maxfield,
      /* bool& */ _err
  );
  return WriteOpalFieldGridFile{_maxfield, _err};
}
bool Bmad::write_opal_lattice_file(int opal_file_unit, LatStruct &lat) {
  bool _err{};
  fortran_write_opal_lattice_file(
      /* int& */ opal_file_unit,
      /* void* */ lat.get_fortran_ptr(),
      /* bool& */ _err
  );
  return _err;
}
bool Bmad::write_time_particle_distribution(
    int time_file_unit,
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<std::string> style,
    optional_ref<BranchStruct> branch,
    std::optional<std::string> format
) {
  const char *_style = style.has_value() ? style->c_str() : nullptr;
  auto *_branch = branch.has_value() ? branch->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_format = format.has_value() ? format->c_str() : nullptr;
  bool _err{};
  fortran_write_time_particle_distribution(
      /* int& */ time_file_unit,
      /* void* */ bunch.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* const char* */ _style,
      /* void* */ _branch,
      /* const char* */ _format,
      /* bool& */ _err
  );
  return _err;
}
void Bmad::xlafun(double x, double y, double z, double res) {
  fortran_xlafun(/* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
int Bmad::xraylib_nist_compound(std::string name) {
  auto _name = name.c_str();
  int _indx{};
  fortran_xraylib_nist_compound(/* const char* */ _name, /* int& */ _indx);
  return _indx;
}
void Bmad::ylafun(double x, double y, double z, double res) {
  fortran_ylafun(/* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}
Bmad::ZAtSurface
Bmad::z_at_surface(EleStruct &ele, double x, double y, std::optional<bool> extend_grid) {
  bool _err_flag{};
  bool extend_grid_lvalue;
  auto *_extend_grid{&extend_grid_lvalue};
  if (extend_grid.has_value()) {
    extend_grid_lvalue = extend_grid.value();
  } else {
    _extend_grid = nullptr;
  }
  // dz_dxy: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _dz_dxy_desc;
  _dz_dxy_desc.rank = 1;
  FixedArray1D<Real, 2> _dz_dxy;
  _dz_dxy_desc.data_ptr = _dz_dxy.data();
  _dz_dxy_desc.dims[0] = _dz_dxy.size();
  double _z{};
  fortran_z_at_surface(/* void* */ ele.get_fortran_ptr(),
                       /* double& */ x,
                       /* double& */ y,
                       /* bool& */ _err_flag,
                       /* bool* */ _extend_grid,
                       /* Bmad::array_descriptor_t& */ _dz_dxy_desc,
                       /* double& */ _z);
  return ZAtSurface{_err_flag, _dz_dxy, _z};
}
void Bmad::zero_ele_kicks(EleStruct &ele) {
  fortran_zero_ele_kicks(/* void* */ ele.get_fortran_ptr());
}
void Bmad::zero_ele_offsets(EleStruct &ele) {
  fortran_zero_ele_offsets(/* void* */ ele.get_fortran_ptr());
}
void Bmad::zero_lr_wakes_in_lat(LatStruct &lat) {
  fortran_zero_lr_wakes_in_lat(/* void* */ lat.get_fortran_ptr());
}
void Bmad::zlafun(double x, double y, double z, double res) {
  fortran_zlafun(/* double& */ x, /* double& */ y, /* double& */ z, /* double& */ res);
}