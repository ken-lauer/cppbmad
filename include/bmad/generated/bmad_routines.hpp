#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace Bmad {
extern "C" void fortran_ab_multipole_kick(
    double &a /* 0D_NOT_real in */,
    double &b /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    int &ref_species /* 0D_NOT_integer in */,
    int &ele_orientation /* 0D_NOT_integer in */,
    void *coord /* 0D_NOT_type in */,
    double &kx /* 0D_NOT_real out */,
    double &ky /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &dk /* 2D_NOT_real out */,
    int *pole_type /* 0D_NOT_integer in */,
    double *scale /* 0D_NOT_real in */
);
struct AbMultipoleKick {
  double kx;
  double ky;
  std::optional<FixedArray2D<Real, 2, 2>> dk;
};
Bmad::AbMultipoleKick ab_multipole_kick(
    double a,
    double b,
    int n,
    int ref_species,
    int ele_orientation,
    CoordStruct &coord,
    std::optional<int> pole_type = std::nullopt,
    std::optional<double> scale = std::nullopt
);
extern "C" void fortran_ab_multipole_kicks(
    Bmad::array_descriptor_t &an /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &bn /* 1D_NOT_real in */,
    int &ix_pole_max /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    int *pole_type /* 0D_NOT_integer in */,
    double *scale /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void ab_multipole_kicks(
    FArray1D<Real> &an,
    FArray1D<Real> &bn,
    int ix_pole_max,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<int> pole_type = std::nullopt,
    std::optional<double> scale = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_absolute_photon_position(
    void *e_orb /* 0D_NOT_type in */,
    void *photon_orb /* 0D_NOT_type inout */
);
void absolute_photon_position(CoordStruct &e_orb, CoordStruct &photon_orb);
extern "C" bool fortran_absolute_time_tracking(
    void *ele /* 0D_NOT_type in */,
    bool &is_abs_time /* 0D_NOT_logical out */
);
bool absolute_time_tracking(EleStruct &ele);
extern "C" bool fortran_ac_kicker_amp(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    double *true_time /* 0D_NOT_real in */,
    double &ac_amp /* 0D_NOT_real out */
);
double
ac_kicker_amp(EleStruct &ele, CoordStruct &orbit, std::optional<double> true_time = std::nullopt);
extern "C" void fortran_action_to_xyz(
    void *ring /* 0D_NOT_type in */,
    int &ix /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &J /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &X /* 1D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ActionToXyz {
  FixedArray1D<Real, 6> X;
  bool err_flag;
};
Bmad::ActionToXyz action_to_xyz(LatStruct &ring, int ix, FixedArray1D<Real, 6> J);
extern "C" void fortran_add_lattice_control_structs(
    void *ele /* 0D_NOT_type in */,
    int *n_add_slave /* 0D_NOT_integer in */,
    int *n_add_lord /* 0D_NOT_integer in */,
    int *n_add_slave_field /* 0D_NOT_integer in */,
    int *n_add_lord_field /* 0D_NOT_integer in */,
    bool *add_at_end /* 0D_NOT_logical in */
);
void add_lattice_control_structs(
    EleStruct &ele,
    std::optional<int> n_add_slave = std::nullopt,
    std::optional<int> n_add_lord = std::nullopt,
    std::optional<int> n_add_slave_field = std::nullopt,
    std::optional<int> n_add_lord_field = std::nullopt,
    std::optional<bool> add_at_end = std::nullopt
);

// Skipped unusable routine add_ptc_layout_to_list:
// - Untranslated type: ptc_branch1_struct (0D)
extern "C" void fortran_add_superimpose(
    void *lat /* 0D_NOT_type inout */,
    void *super_ele_in /* 0D_NOT_type in */,
    int &ix_branch /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *super_ele_out /* 0D_PTR_type out */,
    bool *save_null_drift /* 0D_NOT_logical in */,
    bool *create_jumbo_slave /* 0D_NOT_logical in */,
    int *ix_insert /* 0D_NOT_integer in */,
    bool *mangle_slave_names /* 0D_NOT_logical in */,
    bool *wrap /* 0D_NOT_logical in */
);
struct AddSuperimpose {
  bool err_flag;
  std::optional<EleStruct> super_ele_out;
};
Bmad::AddSuperimpose add_superimpose(
    LatStruct &lat,
    EleStruct &super_ele_in,
    int ix_branch,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> create_jumbo_slave = std::nullopt,
    std::optional<int> ix_insert = std::nullopt,
    std::optional<bool> mangle_slave_names = std::nullopt,
    std::optional<bool> wrap = std::nullopt
);
extern "C" void fortran_add_this_multipass(
    void *lat /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &m_slaves /* 1D_NOT_type inout */,
    void *lord_in /* 0D_NOT_type inout */
);
void add_this_multipass(
    LatStruct &lat,
    LatEleLocStructArray1D m_slaves,
    optional_ref<EleStruct> lord_in = std::nullopt
);

// Skipped unusable routine add_this_name_to_list:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_add_this_taylor_term(
    void *ele /* 0D_NOT_type inout */,
    int &i_out /* 0D_NOT_integer in */,
    double &coef /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &expn /* 1D_NOT_integer inout */
);
void add_this_taylor_term(EleStruct &ele, int i_out, double coef, FixedArray1D<Int, 6> expn);
extern "C" void fortran_adjust_super_slave_names(
    void *lat /* 0D_NOT_type inout */,
    int &ix1_lord /* 0D_NOT_integer in */,
    int &ix2_lord /* 0D_NOT_integer in */,
    bool *first_time /* 0D_NOT_logical in */
);
void adjust_super_slave_names(
    LatStruct &lat,
    int ix1_lord,
    int ix2_lord,
    std::optional<bool> first_time = std::nullopt
);
extern "C" void fortran_allocate_branch_array(
    void *lat /* 0D_NOT_type inout */,
    int &upper_bound /* 0D_NOT_integer in */
);
void allocate_branch_array(LatStruct &lat, int upper_bound);

// Skipped unusable routine allocate_element_array:
// - Routine in configuration skip list
extern "C" void fortran_allocate_lat_ele_array(
    void *lat /* 0D_NOT_type inout */,
    int *upper_bound /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *do_ramper_slave_setup /* 0D_NOT_logical in */
);
void allocate_lat_ele_array(
    LatStruct &lat,
    std::optional<int> upper_bound = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> do_ramper_slave_setup = std::nullopt
);

// Skipped unusable routine allocate_plat:
// - Untranslated type: parser_lat_struct (0D)

// Skipped unusable routine aml_parser:
// - Routine in configuration skip list
extern "C" bool fortran_angle_between_polars(
    void *polar1 /* 0D_NOT_type in */,
    void *polar2 /* 0D_NOT_type in */,
    double &angle /* 0D_NOT_real out */
);
double angle_between_polars(SpinPolarStruct &polar1, SpinPolarStruct &polar2);
extern "C" void fortran_angle_to_canonical_coords(
    void *orbit /* 0D_NOT_type inout */,
    const char *coord_type /* 0D_NOT_character in */
);
void angle_to_canonical_coords(
    CoordStruct &orbit,
    std::optional<std::string> coord_type = std::nullopt
);
extern "C" void fortran_aperture_bookkeeper(void *ele /* 0D_NOT_type inout */);
void aperture_bookkeeper(EleStruct &ele);
extern "C" void fortran_apply_all_rampers(
    void *lat /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool apply_all_rampers(LatStruct &lat);

// Skipped unusable routine apply_element_edge_kick:
// - Untranslated type: fringe_field_info_struct (0D)

// Skipped unusable routine apply_element_edge_kick_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_apply_energy_kick(
    double &dE /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &ddE_dr /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void apply_energy_kick(
    double dE,
    CoordStruct &orbit,
    FixedArray1D<Real, 2> ddE_dr,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);

// Skipped unusable routine apply_fft_3d_kicks:
// - Untranslated type: csr_struct (0D)
extern "C" void fortran_apply_patch_to_ptc_fibre(void *ele /* 0D_NOT_type in */);
void apply_patch_to_ptc_fibre(EleStruct &ele);
extern "C" void fortran_apply_rampers_to_slave(
    void *slave /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool apply_rampers_to_slave(EleStruct &slave);
extern "C" bool fortran_array_re_str(
    Bmad::array_descriptor_t &arr /* 1D_NOT_real inout */,
    const char *parens_in /* 0D_NOT_character in */,
    const char *str_out /* 0D_NOT_character in */
);
void array_re_str(
    FArray1D<Real> &arr,
    std::string str_out,
    std::optional<std::string> parens_in = std::nullopt
);
extern "C" bool fortran_astra_max_field_reference(
    void *pt0 /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type inout */,
    double &field_value /* 0D_NOT_real in */
);
void astra_max_field_reference(GridFieldPt1Struct &pt0, EleStruct &ele, double field_value);
extern "C" bool fortran_at_this_ele_end(
    int &now_at /* 0D_NOT_integer in */,
    int &where_at /* 0D_NOT_integer in */,
    bool &is_at_this_end /* 0D_NOT_logical out */
);
bool at_this_ele_end(int now_at, int where_at);
extern "C" void fortran_attribute_bookkeeper(
    void *ele /* 0D_NOT_type inout */,
    bool *force_bookkeeping /* 0D_NOT_logical in */
);
void attribute_bookkeeper(EleStruct &ele, std::optional<bool> force_bookkeeping = std::nullopt);
extern "C" bool fortran_attribute_free1(
    int &ix_ele /* 0D_NOT_integer in */,
    const char *attrib_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    bool *except_overlay /* 0D_NOT_logical in */,
    bool *dependent_attribs_free /* 0D_NOT_logical in */,
    int &why_not_free /* 0D_NOT_integer out */,
    bool &free /* 0D_NOT_logical out */
);
struct AttributeFree1 {
  int why_not_free;
  bool free;
};
Bmad::AttributeFree1 attribute_free(
    int ix_ele,
    std::string attrib_name,
    LatStruct &lat,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> except_overlay = std::nullopt,
    std::optional<bool> dependent_attribs_free = std::nullopt
);
extern "C" bool fortran_attribute_free2(
    void *ele /* 0D_NOT_type in */,
    const char *attrib_name /* 0D_NOT_character in */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    bool *except_overlay /* 0D_NOT_logical in */,
    bool *dependent_attribs_free /* 0D_NOT_logical in */,
    int &why_not_free /* 0D_NOT_integer out */,
    bool &free /* 0D_NOT_logical out */
);
struct AttributeFree2 {
  int why_not_free;
  bool free;
};
Bmad::AttributeFree2 attribute_free(
    EleStruct &ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> except_overlay = std::nullopt,
    std::optional<bool> dependent_attribs_free = std::nullopt
);
extern "C" bool fortran_attribute_free3(
    int &ix_ele /* 0D_NOT_integer in */,
    int &ix_branch /* 0D_NOT_integer in */,
    const char *attrib_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    bool *except_overlay /* 0D_NOT_logical in */,
    bool *dependent_attribs_free /* 0D_NOT_logical in */,
    int &why_not_free /* 0D_NOT_integer out */,
    bool &free /* 0D_NOT_logical out */
);
struct AttributeFree3 {
  int why_not_free;
  bool free;
};
Bmad::AttributeFree3 attribute_free(
    int ix_ele,
    int ix_branch,
    std::string attrib_name,
    LatStruct &lat,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> except_overlay = std::nullopt,
    std::optional<bool> dependent_attribs_free = std::nullopt
);
extern "C" bool fortran_attribute_index1(
    void *ele /* 0D_NOT_type in */,
    const char *name /* 0D_NOT_character in */,
    const char *full_name /* 0D_NOT_character out */,
    bool *can_abbreviate /* 0D_NOT_logical in */,
    bool *print_error /* 0D_NOT_logical in */,
    int &attrib_index /* 0D_NOT_integer out */
);
struct AttributeIndex1 {
  std::string full_name;
  int attrib_index;
};
Bmad::AttributeIndex1 attribute_index(
    EleStruct &ele,
    std::string name,
    std::optional<bool> can_abbreviate = std::nullopt,
    std::optional<bool> print_error = std::nullopt
);
extern "C" bool fortran_attribute_index2(
    int &key /* 0D_NOT_integer in */,
    const char *name /* 0D_NOT_character in */,
    const char *full_name /* 0D_NOT_character out */,
    bool *can_abbreviate /* 0D_NOT_logical in */,
    bool *print_error /* 0D_NOT_logical in */,
    int &attrib_index /* 0D_NOT_integer out */
);
struct AttributeIndex2 {
  std::string full_name;
  int attrib_index;
};
Bmad::AttributeIndex2 attribute_index(
    int key,
    std::string name,
    std::optional<bool> can_abbreviate = std::nullopt,
    std::optional<bool> print_error = std::nullopt
);

// Skipped unusable routine attribute_info:
// - Untranslated type: ele_attribute_struct (0D)
extern "C" bool fortran_attribute_name1(
    int &key /* 0D_NOT_integer in */,
    int &ix_att /* 0D_NOT_integer in */,
    bool *show_private /* 0D_NOT_logical in */,
    const char *attrib_name /* 0D_NOT_character out */
);
std::string attribute_name(int key, int ix_att, std::optional<bool> show_private = std::nullopt);
extern "C" bool fortran_attribute_name2(
    void *ele /* 0D_NOT_type in */,
    int &ix_att /* 0D_NOT_integer in */,
    bool *show_private /* 0D_NOT_logical in */,
    const char *attrib_name /* 0D_NOT_character out */
);
std::string
attribute_name(EleStruct &ele, int ix_att, std::optional<bool> show_private = std::nullopt);

// Skipped unusable routine attribute_set_bookkeeping:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_attribute_type(
    const char *attrib_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    int &attrib_type /* 0D_NOT_integer out */
);
int attribute_type(std::string attrib_name, optional_ref<EleStruct> ele = std::nullopt);
extern "C" bool fortran_attribute_units(
    const char *attrib_name /* 0D_NOT_character in */,
    const char *unrecognized_units /* 0D_NOT_character in */,
    const char *attrib_units /* 0D_NOT_character out */
);
std::string attribute_units(
    std::string attrib_name,
    std::optional<std::string> unrecognized_units = std::nullopt
);
extern "C" void fortran_autoscale_phase_and_amp(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *scale_phase /* 0D_NOT_logical in */,
    bool *scale_amp /* 0D_NOT_logical in */,
    bool *call_bookkeeper /* 0D_NOT_logical in */
);
bool autoscale_phase_and_amp(
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> scale_phase = std::nullopt,
    std::optional<bool> scale_amp = std::nullopt,
    std::optional<bool> call_bookkeeper = std::nullopt
);
extern "C" bool fortran_average_twiss(
    double &frac1 /* 0D_NOT_real in */,
    void *twiss1 /* 0D_NOT_type in */,
    void *twiss2 /* 0D_NOT_type inout */,
    void *ave_twiss /* 0D_NOT_type out */
);
TwissStruct average_twiss(double frac1, TwissStruct &twiss1, TwissStruct &twiss2);

// Skipped unusable routine bane1:
// - Untranslated type: ibs_struct (0D)
extern "C" void fortran_bbi_kick(
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &sigma /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &nk /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &dnk /* 2D_NOT_real out */,
    bool *linear_kick /* 0D_NOT_logical in */
);
struct BbiKick {
  FixedArray1D<Real, 2> nk;
  FixedArray2D<Real, 2, 2> dnk;
};
Bmad::BbiKick bbi_kick(
    double x,
    double y,
    FixedArray1D<Real, 2> sigma,
    std::optional<bool> linear_kick = std::nullopt
);
extern "C" void fortran_bbi_slice_calc(
    void *ele /* 0D_NOT_type in */,
    int &n_slice /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &z_slice /* 1D_NOT_real inout */
);
void bbi_slice_calc(EleStruct &ele, int n_slice, FArray1D<Real> &z_slice);
extern "C" void fortran_beam_envelope_ibs(
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &ibs_mat /* 2D_NOT_real out */,
    bool &tail_cut /* 0D_NOT_logical in */,
    double &tau /* 0D_NOT_real in */,
    double &energy /* 0D_NOT_real in */,
    double &n_part /* 0D_NOT_real in */,
    int &species /* 0D_NOT_integer in */
);
FixedArray2D<Real, 6, 6> beam_envelope_ibs(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species
);
extern "C" void
fortran_beam_equal_beam(void *beam1 /* 0D_NOT_type inout */, void *beam2 /* 0D_NOT_type in */);
void beam_equal_beam(BeamStruct &beam1, BeamStruct &beam2);
extern "C" bool fortran_beam_init_setup(
    void *beam_init_in /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    int &species /* 0D_NOT_integer in */,
    void *modes /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *beam_init_set /* 0D_NOT_type out */
);
struct BeamInitSetup {
  bool err_flag;
  BeamInitStruct beam_init_set;
};
Bmad::BeamInitSetup beam_init_setup(
    BeamInitStruct &beam_init_in,
    EleStruct &ele,
    int species,
    optional_ref<NormalModesStruct> modes = std::nullopt
);
extern "C" void fortran_beam_tilts(
    Bmad::array_descriptor_t &S /* 2D_NOT_real in */,
    double &angle_xy /* 0D_NOT_real out */,
    double &angle_xz /* 0D_NOT_real out */,
    double &angle_yz /* 0D_NOT_real out */,
    double &angle_xpz /* 0D_NOT_real out */,
    double &angle_ypz /* 0D_NOT_real out */
);
struct BeamTilts {
  double angle_xy;
  double angle_xz;
  double angle_yz;
  double angle_xpz;
  double angle_ypz;
};
Bmad::BeamTilts beam_tilts(FixedArray2D<Real, 6, 6> S);
extern "C" void
fortran_beambeam_fibre_setup(void *ele /* 0D_NOT_type in */, void *ptc_fibre /* 0D_NOT_type out */);
Fibre beambeam_fibre_setup(EleStruct &ele);
extern "C" void fortran_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */,
    bool *track_spin /* 0D_NOT_logical in */
);
void bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> track_spin = std::nullopt
);
extern "C" void fortran_bend_exact_multipole_field(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    bool &local_ref_frame /* 0D_NOT_logical in */,
    void *field /* 0D_NOT_type out */,
    bool *calc_dfield /* 0D_NOT_logical in */,
    bool *calc_potential /* 0D_NOT_logical in */
);
EmFieldStruct bend_exact_multipole_field(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt
);
extern "C" bool fortran_bend_length_has_been_set(
    void *ele /* 0D_NOT_type in */,
    bool &is_set /* 0D_NOT_logical out */
);
bool bend_length_has_been_set(EleStruct &ele);
extern "C" bool fortran_bend_photon_e_rel_init(
    double *r_in /* 0D_NOT_real in */,
    double &E_rel /* 0D_NOT_real out */
);
double bend_photon_e_rel_init(std::optional<double> r_in = std::nullopt);
extern "C" bool fortran_bend_photon_energy_integ_prob(
    double &E_photon /* 0D_NOT_real in */,
    double &g_bend /* 0D_NOT_real in */,
    double &gamma /* 0D_NOT_real in */,
    double &integ_prob /* 0D_NOT_real out */
);
double bend_photon_energy_integ_prob(double E_photon, double g_bend, double gamma);
extern "C" bool fortran_bend_photon_energy_normalized_probability(
    double &E_rel /* 0D_NOT_real in */,
    double &prob /* 0D_NOT_real out */
);
double bend_photon_energy_normalized_probability(double E_rel);
extern "C" void fortran_bend_photon_init(
    double &g_bend_x /* 0D_NOT_real in */,
    double &g_bend_y /* 0D_NOT_real in */,
    double &gamma /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type out */,
    double *E_min /* 0D_NOT_real in */,
    double *E_max /* 0D_NOT_real in */,
    double *E_integ_prob /* 0D_NOT_real in */,
    double *vert_angle_min /* 0D_NOT_real in */,
    double *vert_angle_max /* 0D_NOT_real in */,
    bool *vert_angle_symmetric /* 0D_NOT_logical in */,
    double *emit_probability /* 0D_NOT_real in */
);
CoordStruct bend_photon_init(
    double g_bend_x,
    double g_bend_y,
    double gamma,
    std::optional<double> E_min = std::nullopt,
    std::optional<double> E_max = std::nullopt,
    std::optional<double> E_integ_prob = std::nullopt,
    std::optional<double> vert_angle_min = std::nullopt,
    std::optional<double> vert_angle_max = std::nullopt,
    std::optional<bool> vert_angle_symmetric = std::nullopt,
    std::optional<double> emit_probability = std::nullopt
);
extern "C" void fortran_bend_photon_polarization_init(
    double &g_bend_x /* 0D_NOT_real in */,
    double &g_bend_y /* 0D_NOT_real in */,
    double &E_rel /* 0D_NOT_real in */,
    double &gamma_phi /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type out */
);
CoordStruct
bend_photon_polarization_init(double g_bend_x, double g_bend_y, double E_rel, double gamma_phi);
extern "C" bool fortran_bend_photon_vert_angle_init(
    double &E_rel /* 0D_NOT_real in */,
    double &gamma /* 0D_NOT_real in */,
    double *r_in /* 0D_NOT_real in */,
    bool *invert /* 0D_NOT_logical in */,
    double &phi /* 0D_NOT_real out */
);
double bend_photon_vert_angle_init(
    double E_rel,
    double gamma,
    std::optional<double> r_in = std::nullopt,
    std::optional<bool> invert = std::nullopt
);
extern "C" bool fortran_bend_shift(
    void *position1 /* 0D_NOT_type in */,
    double &g /* 0D_NOT_real in */,
    double &delta_s /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    double *ref_tilt /* 0D_NOT_real in */,
    void *position2 /* 0D_NOT_type out */
);
struct BendShift {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct position2;
};
Bmad::BendShift bend_shift(
    FloorPositionStruct &position1,
    double g,
    double delta_s,
    std::optional<double> ref_tilt = std::nullopt
);
extern "C" bool fortran_bend_vert_angle_integ_prob(
    double &vert_angle /* 0D_NOT_real in */,
    double &E_rel /* 0D_NOT_real in */,
    double &gamma /* 0D_NOT_real in */,
    double &integ_prob /* 0D_NOT_real out */
);
double bend_vert_angle_integ_prob(double vert_angle, double E_rel, double gamma);

// Skipped unusable routine bjmt1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine bjmt_integrand:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine bl_via_mat:
// - Untranslated type: ibs_sim_param_struct (0D)
extern "C" void fortran_bl_via_vlassov(
    double &current /* 0D_NOT_real in */,
    double &alpha /* 0D_NOT_real in */,
    double &Energy /* 0D_NOT_real in */,
    double &sigma_p /* 0D_NOT_real in */,
    double &Vrf /* 0D_NOT_real in */,
    double &omega /* 0D_NOT_real in */,
    double &U0 /* 0D_NOT_real in */,
    double &circ /* 0D_NOT_real in */,
    double &R /* 0D_NOT_real in */,
    double &L /* 0D_NOT_real in */,
    double &sigma_z /* 0D_NOT_real out */
);
double bl_via_vlassov(
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
);

// Skipped unusable routine bmad_and_xsif_parser:
// - Routine in configuration skip list
extern "C" void fortran_bmad_parser(
    const char *lat_file /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type out */,
    bool *make_mats6 /* 0D_NOT_logical in */,
    bool &digested_read_ok /* 0D_NOT_logical out */,
    const char *use_line /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *parse_lat /* 0D_NOT_type out */
);
struct BmadParser {
  LatStruct lat;
  bool digested_read_ok;
  bool err_flag;
  LatStruct parse_lat;
};
Bmad::BmadParser bmad_parser(
    std::string lat_file,
    std::optional<bool> make_mats6 = std::nullopt,
    std::optional<std::string> use_line = std::nullopt
);
extern "C" void fortran_bmad_parser2(
    const char *lat_file /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &orbit /* 1D_NOT_type in */,
    bool *make_mats6 /* 0D_NOT_logical in */,
    bool *err_flag /* 0D_NOT_logical in */,
    void *parse_lat /* 0D_NOT_type in */
);
void bmad_parser2(
    std::string lat_file,
    LatStruct &lat,
    std::optional<CoordStructArray1D> orbit = std::nullopt,
    std::optional<bool> make_mats6 = std::nullopt,
    std::optional<bool> err_flag = std::nullopt,
    optional_ref<LatStruct> parse_lat = std::nullopt
);

// Skipped unusable routine bmad_parser_string_attribute_set:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_bmad_patch_parameters_to_ptc(
    Bmad::array_descriptor_t &ang /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &exi /* 2D_NOT_real inout */
);
void bmad_patch_parameters_to_ptc(FixedArray1D<Real, 3> ang, FixedArray2D<Real, 3, 3> exi);

// Skipped unusable routine bmad_taylor_equal_damap:
// - Untranslated type: damap (0D)

// Skipped unusable routine bmad_taylors_equal_ptc_taylors:
// - Untranslated type: taylor (1D)

// Skipped unusable routine bmad_taylors_equal_reals_8:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_bp_set_ran_status();
void bp_set_ran_status();
extern "C" void fortran_branch_equal_branch(
    void *branch1 /* 0D_NOT_type inout */,
    void *branch2 /* 0D_NOT_type in */
);
void branch_equal_branch(BranchStruct &branch1, BranchStruct &branch2);
extern "C" bool
fortran_branch_name(void *branch /* 0D_NOT_type in */, const char *name /* 0D_NOT_character out */);
std::string branch_name(BranchStruct &branch);
extern "C" void fortran_branch_to_ptc_m_u(void *branch /* 0D_NOT_type inout */);
void branch_to_ptc_m_u(BranchStruct &branch);
extern "C" void
fortran_bunch_equal_bunch(void *bunch1 /* 0D_NOT_type inout */, void *bunch2 /* 0D_NOT_type in */);
void bunch_equal_bunch(BunchStruct &bunch1, BunchStruct &bunch2);

// Skipped unusable routine c_multi:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_c_to_cbar(
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &cbar_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 2, 2> c_to_cbar(EleStruct &ele);
extern "C" void fortran_calc_bunch_params(
    void *bunch /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type out */,
    bool &error /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &n_mat /* 2D_NOT_real out */,
    bool *is_time_coords /* 0D_NOT_logical in */,
    void *ele /* 0D_NOT_type in */
);
struct CalcBunchParams {
  BunchParamsStruct bunch_params;
  bool error;
  std::optional<FixedArray2D<Real, 6, 6>> n_mat;
};
Bmad::CalcBunchParams calc_bunch_params(
    BunchStruct &bunch,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleStruct> ele = std::nullopt
);
extern "C" void fortran_calc_bunch_params_slice(
    void *bunch /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type inout */,
    int &plane /* 0D_NOT_integer in */,
    double &slice_center /* 0D_NOT_real in */,
    double &slice_spread /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    bool *is_time_coords /* 0D_NOT_logical in */,
    void *ele /* 0D_NOT_type in */
);
bool calc_bunch_params_slice(
    BunchStruct &bunch,
    BunchParamsStruct &bunch_params,
    int plane,
    double slice_center,
    double slice_spread,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleStruct> ele = std::nullopt
);
extern "C" void fortran_calc_bunch_params_z_slice(
    void *bunch /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &slice_bounds /* 1D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    bool *is_time_coords /* 0D_NOT_logical in */,
    void *ele /* 0D_NOT_type in */
);
bool calc_bunch_params_z_slice(
    BunchStruct &bunch,
    BunchParamsStruct &bunch_params,
    FixedArray1D<Real, 2> slice_bounds,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleStruct> ele = std::nullopt
);
extern "C" void fortran_calc_bunch_sigma_matrix_etc(
    Bmad::array_descriptor_t &particle /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &charge /* 1D_NOT_real in */,
    void *bunch_params /* 0D_NOT_type out */,
    bool *is_time_coords /* 0D_NOT_logical in */,
    void *ele /* 0D_NOT_type inout */
);
BunchParamsStruct calc_bunch_sigma_matrix_etc(
    CoordStructArray1D particle,
    FArray1D<Real> &charge,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleStruct> ele = std::nullopt
);

// Skipped unusable routine calc_density_derivative_complex:
// - Variable in sized array: 3D_NOT_real
// - Variable inout sized array: 3D_NOT_complex
extern "C" void fortran_calc_emittances_and_twiss_from_sigma_matrix(
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real in */,
    void *bunch_params /* 0D_NOT_type out */,
    bool &error /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &n_mat /* 2D_NOT_real out */
);
struct CalcEmittancesAndTwissFromSigmaMatrix {
  BunchParamsStruct bunch_params;
  bool error;
  std::optional<FixedArray2D<Real, 6, 6>> n_mat;
};
Bmad::CalcEmittancesAndTwissFromSigmaMatrix calc_emittances_and_twiss_from_sigma_matrix(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<bool> print_err = std::nullopt
);

// Skipped unusable routine calc_next_fringe_edge:
// - Untranslated type: fringe_field_info_struct (0D)
extern "C" void fortran_calc_spin_params(
    void *bunch /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type out */
);
BunchParamsStruct calc_spin_params(BunchStruct &bunch);
extern "C" void fortran_calc_super_slave_key(
    void *lord1 /* 0D_NOT_type in */,
    void *lord2 /* 0D_NOT_type in */,
    void *slave /* 0D_NOT_type out */,
    bool *create_jumbo_slave /* 0D_NOT_logical in */
);
EleStruct calc_super_slave_key(
    EleStruct &lord1,
    EleStruct &lord2,
    std::optional<bool> create_jumbo_slave = std::nullopt
);
extern "C" void fortran_calc_wall_radius(
    Bmad::array_descriptor_t &v /* 1D_NOT_type in */,
    double &cos_ang /* 0D_NOT_real in */,
    double &sin_ang /* 0D_NOT_real in */,
    double &r_wall /* 0D_NOT_real out */,
    double &dr_dtheta /* 0D_NOT_real out */,
    int &ix_vertex /* 0D_NOT_integer out */
);
struct CalcWallRadius {
  double r_wall;
  double dr_dtheta;
  int ix_vertex;
};
Bmad::CalcWallRadius calc_wall_radius(Wall3dVertexStructArray1D v, double cos_ang, double sin_ang);

// Skipped unusable routine calc_wiggler_g_params:
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)
extern "C" void fortran_calc_z_tune(void *branch /* 0D_NOT_type inout */);
void calc_z_tune(BranchStruct &branch);
extern "C" void fortran_canonical_to_angle_coords(
    void *orbit /* 0D_NOT_type inout */,
    const char *coord_type /* 0D_NOT_character in */
);
void canonical_to_angle_coords(
    CoordStruct &orbit,
    std::optional<std::string> coord_type = std::nullopt
);

// Skipped unusable routine capillary_photon_hit_spot_calc:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_propagate_photon_a_step:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_reflect_photon:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_track_photon_to_wall:
// - Untranslated type: photon_track_struct (0D)
extern "C" void fortran_cbar_to_c(
    Bmad::array_descriptor_t &cbar_mat /* 2D_NOT_real in */,
    void *a /* 0D_NOT_type in */,
    void *b /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &c_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 2, 2>
cbar_to_c(FixedArray2D<Real, 2, 2> cbar_mat, TwissStruct &a, TwissStruct &b);

// Skipped unusable routine ccfft3d:
// - Routine module (fft_interface_mod) in configuration skip list

// Skipped unusable routine ccfftam:
// - Routine in configuration skip list
extern "C" void fortran_check_aperture_limit(
    void *orb /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *param /* 0D_NOT_type inout */,
    void *old_orb /* 0D_NOT_type in */,
    bool *check_momentum /* 0D_NOT_logical in */
);
void check_aperture_limit(
    CoordStruct &orb,
    EleStruct &ele,
    int particle_at,
    LatParamStruct &param,
    optional_ref<CoordStruct> old_orb = std::nullopt,
    std::optional<bool> check_momentum = std::nullopt
);

// Skipped unusable routine check_aperture_limit_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_check_controller_controls(
    int &ele_key /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &contrl /* 1D_NOT_type in */,
    const char *name /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical out */
);
bool check_controller_controls(int ele_key, ControlStructArray1D contrl, std::string name);
extern "C" void fortran_check_for_superimpose_problem(
    void *branch /* 0D_NOT_type inout */,
    void *super_ele /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */,
    void *ref_ele /* 0D_NOT_type inout */,
    bool &wrap /* 0D_NOT_logical in */
);
void check_for_superimpose_problem(
    BranchStruct &branch,
    EleStruct &super_ele,
    bool err_flag,
    bool wrap,
    optional_ref<EleStruct> ref_ele = std::nullopt
);
extern "C" void fortran_check_if_s_in_bounds(
    void *branch /* 0D_NOT_type in */,
    double &s /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &translated_s /* 0D_NOT_real out */,
    bool *print_err /* 0D_NOT_logical in */
);
struct CheckIfSInBounds {
  bool err_flag;
  double translated_s;
};
Bmad::CheckIfSInBounds
check_if_s_in_bounds(BranchStruct &branch, double s, std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_choose_quads_for_set_tune(
    void *branch /* 0D_NOT_type in */,
    void *dk1 /* 1D_ALLOC_real out */,
    void *eles /* 1D_ALLOC_type out */,
    const char *mask /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ChooseQuadsForSetTune {
  RealAlloc1D dk1;
  ElePointerStructAlloc1D eles;
  bool err_flag;
};
Bmad::ChooseQuadsForSetTune
choose_quads_for_set_tune(BranchStruct &branch, std::optional<std::string> mask = std::nullopt);
extern "C" void fortran_chrom_calc(
    void *lat /* 0D_NOT_type in */,
    double &delta_e /* 0D_NOT_real inout */,
    double &chrom_a /* 0D_NOT_real out */,
    double &chrom_b /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    double *pz /* 0D_NOT_real in */,
    void *low_E_lat /* 0D_NOT_type out */,
    void *high_E_lat /* 0D_NOT_type out */,
    void *low_E_orb /* 1D_ALLOC_type out */,
    void *high_E_orb /* 1D_ALLOC_type out */,
    int *ix_branch /* 0D_NOT_integer in */,
    void *orb0 /* 0D_NOT_type in */
);
struct ChromCalc {
  double chrom_a;
  double chrom_b;
  bool err_flag;
  LatStruct low_E_lat;
  LatStruct high_E_lat;
  CoordStructAlloc1D low_E_orb;
  CoordStructAlloc1D high_E_orb;
};
Bmad::ChromCalc chrom_calc(
    LatStruct &lat,
    double &delta_e,
    std::optional<double> pz = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<CoordStruct> orb0 = std::nullopt
);
extern "C" void fortran_chrom_tune(
    void *lat /* 0D_NOT_type inout */,
    double &delta_e /* 0D_NOT_real inout */,
    double &target_x /* 0D_NOT_real in */,
    double &target_y /* 0D_NOT_real in */,
    double &err_tol /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool chrom_tune(LatStruct &lat, double &delta_e, double target_x, double target_y, double err_tol);

// Skipped unusable routine cimp1:
// - Untranslated type: ibs_struct (0D)
extern "C" bool fortran_classical_radius(
    int &species /* 0D_NOT_integer in */,
    double &radius /* 0D_NOT_real out */
);
double classical_radius(int species);
extern "C" void fortran_clear_lat_1turn_mats(void *lat /* 0D_NOT_type out */);
LatStruct clear_lat_1turn_mats();
extern "C" void fortran_clear_taylor_maps_from_elements(void *lat /* 0D_NOT_type inout */);
void clear_taylor_maps_from_elements(LatStruct &lat);
extern "C" void fortran_closed_orbit_calc(
    void *lat /* 0D_NOT_type inout */,
    void *closed_orb /* 1D_ALLOC_type inout */,
    int *i_dim /* 0D_NOT_integer in */,
    int *direction /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
bool closed_orbit_calc(
    LatStruct &lat,
    CoordStructAlloc1D closed_orb,
    std::optional<int> i_dim = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_closed_orbit_from_tracking(
    void *lat /* 0D_NOT_type in */,
    void *closed_orb /* 1D_ALLOC_type out */,
    int &i_dim /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &eps_rel /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &eps_abs /* 1D_NOT_real in */,
    void *init_guess /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ClosedOrbitFromTracking {
  CoordStructAlloc1D closed_orb;
  bool err_flag;
};
Bmad::ClosedOrbitFromTracking closed_orbit_from_tracking(
    LatStruct &lat,
    int i_dim,
    std::optional<FArray1D<Real>> eps_rel = std::nullopt,
    std::optional<FArray1D<Real>> eps_abs = std::nullopt,
    optional_ref<CoordStruct> init_guess = std::nullopt
);
extern "C" bool fortran_cmplx_re_str(
    std::complex<double> &cmp /* 0D_NOT_complex in */,
    const char *str_out /* 0D_NOT_character in */
);
void cmplx_re_str(std::complex<double> cmp, std::string str_out);
extern "C" void fortran_combine_consecutive_elements(
    void *lat /* 0D_NOT_type inout */,
    bool &error /* 0D_NOT_logical out */
);
bool combine_consecutive_elements(LatStruct &lat);
extern "C" void fortran_complex_taylor_clean(void *complex_taylor /* 0D_NOT_type inout */);
void complex_taylor_clean(ComplexTaylorStruct &complex_taylor);
extern "C" bool fortran_complex_taylor_coef1(
    void *complex_taylor /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &exp /* 1D_NOT_integer in */,
    std::complex<double> &coef /* 0D_NOT_complex in */
);
void complex_taylor_coef(
    ComplexTaylorStruct &complex_taylor,
    FArray1D<Int> &exp,
    std::complex<double> coef
);
extern "C" bool fortran_complex_taylor_coef2(
    void *complex_taylor /* 0D_NOT_type in */,
    int *i1 /* 0D_NOT_integer in */,
    int *i2 /* 0D_NOT_integer in */,
    int *i3 /* 0D_NOT_integer in */,
    int *i4 /* 0D_NOT_integer in */,
    int *i5 /* 0D_NOT_integer in */,
    int *i6 /* 0D_NOT_integer in */,
    int *i7 /* 0D_NOT_integer in */,
    int *i8 /* 0D_NOT_integer in */,
    int *i9 /* 0D_NOT_integer in */,
    std::complex<double> &coef /* 0D_NOT_complex in */
);
void complex_taylor_coef(
    ComplexTaylorStruct &complex_taylor,
    std::complex<double> coef,
    std::optional<int> i1 = std::nullopt,
    std::optional<int> i2 = std::nullopt,
    std::optional<int> i3 = std::nullopt,
    std::optional<int> i4 = std::nullopt,
    std::optional<int> i5 = std::nullopt,
    std::optional<int> i6 = std::nullopt,
    std::optional<int> i7 = std::nullopt,
    std::optional<int> i8 = std::nullopt,
    std::optional<int> i9 = std::nullopt
);

// Skipped unusable routine complex_taylor_equal_c_taylor:
// - Untranslated type: c_taylor (0D)
extern "C" void fortran_complex_taylor_equal_complex_taylor(
    void *complex_taylor1 /* 0D_NOT_type inout */,
    void *complex_taylor2 /* 0D_NOT_type in */
);
void complex_taylor_equal_complex_taylor(
    ComplexTaylorStruct &complex_taylor1,
    ComplexTaylorStruct &complex_taylor2
);
extern "C" bool fortran_complex_taylor_exponent_index(
    Bmad::array_descriptor_t &expn /* 1D_NOT_integer in */,
    int &index /* 0D_NOT_integer out */
);
int complex_taylor_exponent_index(FixedArray1D<Int, 6> expn);
extern "C" void
fortran_complex_taylor_make_unit(Bmad::array_descriptor_t &complex_taylor /* 1D_NOT_type inout */);
void complex_taylor_make_unit(ComplexTaylorStructArray1D complex_taylor);
extern "C" void fortran_complex_taylor_to_mat6(
    Bmad::array_descriptor_t &a_complex_taylor /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &r_in /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &vec0 /* 1D_NOT_complex out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &r_out /* 1D_NOT_complex inout */
);
struct ComplexTaylorToMat6 {
  FixedArray1D<Complex, 6> vec0;
  FixedArray2D<Complex, 6, 6> mat6;
};
Bmad::ComplexTaylorToMat6 complex_taylor_to_mat6(
    ComplexTaylorStructArray1D a_complex_taylor,
    FArray1D<Complex> &r_in,
    std::optional<FArray1D<Complex>> r_out = std::nullopt
);

// Skipped unusable routine complex_taylors_equal_c_taylors:
// - Untranslated type: c_taylor (1D)
extern "C" void fortran_complex_taylors_equal_complex_taylors(
    Bmad::array_descriptor_t &complex_taylor1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &complex_taylor2 /* 1D_NOT_type in */
);
void complex_taylors_equal_complex_taylors(
    ComplexTaylorStructArray1D complex_taylor1,
    ComplexTaylorStructArray1D complex_taylor2
);
extern "C" void fortran_compute_slave_coupler(void *slave /* 0D_NOT_type inout */);
void compute_slave_coupler(EleStruct &slave);

// Skipped unusable routine compute_super_lord_s:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_concat_ele_taylor(
    Bmad::array_descriptor_t &orb_taylor /* 1D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &spin_taylor /* 1D_NOT_type inout */
);
bool concat_ele_taylor(
    TaylorStructArray1D orb_taylor,
    EleStruct &ele,
    std::optional<TaylorStructArray1D> spin_taylor = std::nullopt
);

// Skipped unusable routine concat_real_8:
// - Untranslated type: real_8 (1D)
// - Untranslated type: real_8 (1D)
// - Untranslated type: real_8 (1D)
extern "C" void fortran_concat_taylor(
    Bmad::array_descriptor_t &taylor1 /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &taylor2 /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &taylor3 /* 1D_NOT_type inout */
);
void concat_taylor(
    TaylorStructArray1D taylor1,
    TaylorStructArray1D taylor2,
    TaylorStructArray1D taylor3
);
extern "C" void fortran_concat_transfer_mat(
    Bmad::array_descriptor_t &mat_1 /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &vec_1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &mat_0 /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &vec_0 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &mat_out /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &vec_out /* 1D_NOT_real out */
);
struct ConcatTransferMat {
  FixedArray2D<Real, 6, 6> mat_out;
  FixedArray1D<Real, 6> vec_out;
};
Bmad::ConcatTransferMat concat_transfer_mat(
    FixedArray2D<Real, 6, 6> mat_1,
    FixedArray1D<Real, 6> vec_1,
    FixedArray2D<Real, 6, 6> mat_0,
    FixedArray1D<Real, 6> vec_0
);
extern "C" void fortran_control_bookkeeper(
    void *lat /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    bool *err_flag /* 0D_NOT_logical in */
);
void control_bookkeeper(
    LatStruct &lat,
    optional_ref<EleStruct> ele = std::nullopt,
    std::optional<bool> err_flag = std::nullopt
);

// Skipped unusable routine conv3d:
// - Array bounds handling: "Enum 'RILO' found in bounds 'rilo' but not in provided map."
// - Array bounds handling: "Enum 'G1ILO' found in bounds 'g1ilo' but not in provided map."
// - Array bounds handling: "Enum 'CILO' found in bounds 'cilo' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_convert_bend_exact_multipole(
    double &g /* 0D_NOT_real in */,
    int &out_type /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &an /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &bn /* 1D_NOT_real inout */
);
void convert_bend_exact_multipole(
    double g,
    int out_type,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> an,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> bn
);
extern "C" void fortran_convert_coords(
    const char *in_type_str /* 0D_NOT_character in */,
    void *coord_in /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    const char *out_type_str /* 0D_NOT_character out */,
    void *coord_out /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ConvertCoords {
  std::string out_type_str;
  CoordStruct coord_out;
  bool err_flag;
};
Bmad::ConvertCoords convert_coords(std::string in_type_str, CoordStruct &coord_in, EleStruct &ele);
extern "C" void fortran_convert_field_ele_to_lab(
    void *ele /* 0D_NOT_type in */,
    double &s_here /* 0D_NOT_real in */,
    bool &forward_transform /* 0D_NOT_logical in */,
    void *field /* 0D_NOT_type out */,
    bool *calc_dfield /* 0D_NOT_logical in */,
    bool *calc_potential /* 0D_NOT_logical in */
);
EmFieldStruct convert_field_ele_to_lab(
    EleStruct &ele,
    double s_here,
    bool forward_transform,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt
);
extern "C" void fortran_convert_local_cartesian_to_local_curvilinear(
    double &x /* 0D_NOT_real in */,
    double &z /* 0D_NOT_real in */,
    double &g /* 0D_NOT_real in */,
    double &xout /* 0D_NOT_real in */,
    double &sout /* 0D_NOT_real in */
);
void convert_local_cartesian_to_local_curvilinear(
    double x,
    double z,
    double g,
    double xout,
    double sout
);
extern "C" void fortran_convert_local_curvilinear_to_local_cartesian(
    double &x /* 0D_NOT_real in */,
    double &s /* 0D_NOT_real in */,
    double &g /* 0D_NOT_real in */,
    double &xout /* 0D_NOT_real in */,
    double &zout /* 0D_NOT_real in */
);
void convert_local_curvilinear_to_local_cartesian(
    double x,
    double s,
    double g,
    double xout,
    double zout
);
extern "C" void fortran_convert_particle_coordinates_s_to_t(
    void *particle /* 0D_NOT_type inout */,
    double &s_body /* 0D_NOT_real in */,
    int &orientation /* 0D_NOT_integer in */
);
void convert_particle_coordinates_s_to_t(CoordStruct &particle, double s_body, int orientation);
extern "C" void fortran_convert_particle_coordinates_t_to_s(
    void *particle /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    double &s_body /* 0D_NOT_real out */,
    bool *use_downstream_p0c /* 0D_NOT_logical in */
);
double convert_particle_coordinates_t_to_s(
    CoordStruct &particle,
    EleStruct &ele,
    std::optional<bool> use_downstream_p0c = std::nullopt
);
extern "C" void fortran_convert_pc_to(
    double &pc /* 0D_NOT_real in */,
    int &particle /* 0D_NOT_integer in */,
    double &E_tot /* 0D_NOT_real out */,
    double &gamma /* 0D_NOT_real out */,
    double &kinetic /* 0D_NOT_real out */,
    double &beta /* 0D_NOT_real out */,
    double &brho /* 0D_NOT_real out */,
    double &beta1 /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ConvertPcTo {
  double E_tot;
  double gamma;
  double kinetic;
  double beta;
  double brho;
  double beta1;
  bool err_flag;
};
Bmad::ConvertPcTo convert_pc_to(double pc, int particle);
extern "C" void fortran_convert_total_energy_to(
    double &E_tot /* 0D_NOT_real in */,
    int &particle /* 0D_NOT_integer in */,
    double &gamma /* 0D_NOT_real out */,
    double &kinetic /* 0D_NOT_real out */,
    double &beta /* 0D_NOT_real out */,
    double &pc /* 0D_NOT_real out */,
    double &brho /* 0D_NOT_real out */,
    double &beta1 /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
struct ConvertTotalEnergyTo {
  double gamma;
  double kinetic;
  double beta;
  double pc;
  double brho;
  double beta1;
  bool err_flag;
};
Bmad::ConvertTotalEnergyTo
convert_total_energy_to(double E_tot, int particle, std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_converter_distribution_parser(
    void *ele /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ConverterDistributionParser {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
Bmad::ConverterDistributionParser converter_distribution_parser(EleStruct &ele);
extern "C" void
fortran_coord_equal_coord(void *coord1 /* 0D_NOT_type out */, void *coord2 /* 0D_NOT_type in */);
CoordStruct coord_equal_coord(CoordStruct &coord2);
extern "C" bool fortran_coord_state_name(
    int &coord_state /* 0D_NOT_integer in */,
    bool *one_word /* 0D_NOT_logical in */,
    const char *state_str /* 0D_NOT_character out */
);
std::string coord_state_name(int coord_state, std::optional<bool> one_word = std::nullopt);
extern "C" bool fortran_coords_body_to_local(
    void *body_position /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    bool *calculate_angles /* 0D_NOT_logical in */,
    void *local_position /* 0D_NOT_type out */
);
struct CoordsBodyToLocal {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct local_position;
};
Bmad::CoordsBodyToLocal coords_body_to_local(
    FloorPositionStruct &body_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles = std::nullopt
);
extern "C" bool fortran_coords_body_to_rel_exit(
    void *body_position /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    bool *calculate_angles /* 0D_NOT_logical in */,
    void *rel_exit /* 0D_NOT_type out */
);
struct CoordsBodyToRelExit {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct rel_exit;
};
Bmad::CoordsBodyToRelExit coords_body_to_rel_exit(
    FloorPositionStruct &body_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles = std::nullopt
);
extern "C" bool fortran_coords_curvilinear_to_floor(
    Bmad::array_descriptor_t &xys /* 1D_NOT_real in */,
    void *branch /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *global /* 0D_NOT_type out */
);
struct CoordsCurvilinearToFloor {
  bool err_flag;
  FloorPositionStruct global;
};
Bmad::CoordsCurvilinearToFloor
coords_curvilinear_to_floor(FixedArray1D<Real, 3> xys, BranchStruct &branch);
extern "C" bool fortran_coords_floor_to_curvilinear(
    void *floor_coords /* 0D_NOT_type in */,
    void *ele0 /* 0D_NOT_type in */,
    void *ele1 /* 0D_PTR_type out */,
    int &status /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    void *local_coords /* 0D_NOT_type out */
);
struct CoordsFloorToCurvilinear {
  std::optional<EleStruct> ele1;
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct local_coords;
};
Bmad::CoordsFloorToCurvilinear
coords_floor_to_curvilinear(FloorPositionStruct &floor_coords, EleStruct &ele0);
extern "C" bool fortran_coords_floor_to_local_curvilinear(
    void *global_position /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    int &status /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    int *relative_to /* 0D_NOT_integer in */,
    void *local_position /* 0D_NOT_type out */
);
struct CoordsFloorToLocalCurvilinear {
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct local_position;
};
Bmad::CoordsFloorToLocalCurvilinear coords_floor_to_local_curvilinear(
    FloorPositionStruct &global_position,
    EleStruct &ele,
    std::optional<int> relative_to = std::nullopt
);
extern "C" bool fortran_coords_floor_to_relative(
    void *floor0 /* 0D_NOT_type in */,
    void *global_position /* 0D_NOT_type in */,
    bool *calculate_angles /* 0D_NOT_logical in */,
    bool *is_delta_position /* 0D_NOT_logical in */,
    void *local_position /* 0D_NOT_type out */
);
FloorPositionStruct coords_floor_to_relative(
    FloorPositionStruct &floor0,
    FloorPositionStruct &global_position,
    std::optional<bool> calculate_angles = std::nullopt,
    std::optional<bool> is_delta_position = std::nullopt
);
extern "C" bool fortran_coords_local_curvilinear_to_body(
    void *local_position /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    bool *calculate_angles /* 0D_NOT_logical in */,
    void *body_position /* 0D_NOT_type out */
);
struct CoordsLocalCurvilinearToBody {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct body_position;
};
Bmad::CoordsLocalCurvilinearToBody coords_local_curvilinear_to_body(
    FloorPositionStruct &local_position,
    EleStruct &ele,
    std::optional<bool> calculate_angles = std::nullopt
);
extern "C" bool fortran_coords_local_curvilinear_to_floor(
    void *local_position /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    bool *in_body_frame /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    bool *calculate_angles /* 0D_NOT_logical in */,
    int *end_origin /* 0D_NOT_integer in */,
    bool *downstream_dir_ref /* 0D_NOT_logical in */,
    void *global_position /* 0D_NOT_type out */
);
struct CoordsLocalCurvilinearToFloor {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  FloorPositionStruct global_position;
};
Bmad::CoordsLocalCurvilinearToFloor coords_local_curvilinear_to_floor(
    FloorPositionStruct &local_position,
    EleStruct &ele,
    std::optional<bool> in_body_frame = std::nullopt,
    std::optional<bool> calculate_angles = std::nullopt,
    std::optional<int> end_origin = std::nullopt,
    std::optional<bool> downstream_dir_ref = std::nullopt
);
extern "C" bool fortran_coords_relative_to_floor(
    void *floor0 /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &dr /* 1D_NOT_real in */,
    double *theta /* 0D_NOT_real in */,
    double *phi /* 0D_NOT_real in */,
    double *psi /* 0D_NOT_real in */,
    void *floor1 /* 0D_NOT_type out */
);
FloorPositionStruct coords_relative_to_floor(
    FloorPositionStruct &floor0,
    FixedArray1D<Real, 3> dr,
    std::optional<double> theta = std::nullopt,
    std::optional<double> phi = std::nullopt,
    std::optional<double> psi = std::nullopt
);

// Skipped unusable routine cos_phi:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" bool fortran_coulombfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void coulombfun(double u, double v, double w, double gam, double res);
extern "C" void fortran_create_concatenated_wall3d(
    void *lat /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical in */
);
void create_concatenated_wall3d(LatStruct &lat, bool err);
extern "C" void fortran_create_element_slice(
    void *sliced_ele /* 0D_NOT_type out */,
    void *ele_in /* 0D_NOT_type in */,
    double &l_slice /* 0D_NOT_real in */,
    double &offset /* 0D_NOT_real in */,
    void *param /* 0D_NOT_type in */,
    bool &include_upstream_end /* 0D_NOT_logical in */,
    bool &include_downstream_end /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *old_slice /* 0D_NOT_type in */,
    void *orb_in /* 0D_NOT_type in */
);
struct CreateElementSlice {
  EleStruct sliced_ele;
  bool err_flag;
};
Bmad::CreateElementSlice create_element_slice(
    EleStruct &ele_in,
    double l_slice,
    double offset,
    LatParamStruct &param,
    bool include_upstream_end,
    bool include_downstream_end,
    optional_ref<EleStruct> old_slice = std::nullopt,
    optional_ref<CoordStruct> orb_in = std::nullopt
);

// Skipped unusable routine create_feedback:
// - Variable-sized in character array: 1D_NOT_character
// - Variable-sized in character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_create_field_overlap(
    void *lat /* 0D_NOT_type inout */,
    const char *lord_name /* 0D_NOT_character in */,
    const char *slave_name /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool create_field_overlap(LatStruct &lat, std::string lord_name, std::string slave_name);
extern "C" void fortran_create_girder(
    void *lat /* 0D_NOT_type inout */,
    int &ix_girder /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &contrl /* 1D_NOT_type in */,
    void *girder_info /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void create_girder(
    LatStruct &lat,
    int ix_girder,
    ControlStructArray1D contrl,
    EleStruct &girder_info,
    bool err_flag
);
extern "C" void fortran_create_group(
    void *lord /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &contrl /* 1D_NOT_type in */,
    bool &err /* 0D_NOT_logical in */
);
void create_group(EleStruct &lord, ControlStructArray1D contrl, bool err);
extern "C" void fortran_create_lat_ele_nametable(
    void *lat /* 0D_NOT_type in */,
    void *nametable /* 0D_NOT_type out */
);
NametableStruct create_lat_ele_nametable(LatStruct &lat);
extern "C" void fortran_create_overlay(
    void *lord /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &contrl /* 1D_NOT_type in */,
    bool &err /* 0D_NOT_logical in */
);
void create_overlay(EleStruct &lord, ControlStructArray1D contrl, bool err);
extern "C" void fortran_create_planar_wiggler_model(
    void *wiggler_in /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
struct CreatePlanarWigglerModel {
  LatStruct lat;
  bool err_flag;
};
Bmad::CreatePlanarWigglerModel
create_planar_wiggler_model(EleStruct &wiggler_in, std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_create_ramper(
    void *lord /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &contrl /* 1D_NOT_type in */,
    bool &err /* 0D_NOT_logical in */
);
void create_ramper(EleStruct &lord, ControlStructArray1D contrl, bool err);
extern "C" void fortran_create_sol_quad_model(
    void *sol_quad /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */
);
void create_sol_quad_model(EleStruct &sol_quad, LatStruct &lat);
extern "C" void fortran_create_unique_ele_names(
    void *lat /* 0D_NOT_type inout */,
    int &key /* 0D_NOT_integer in */,
    const char *suffix /* 0D_NOT_character in */
);
void create_unique_ele_names(LatStruct &lat, int key, std::string suffix);
extern "C" void fortran_create_wiggler_cartesian_map(
    void *ele /* 0D_NOT_type in */,
    void *cart_map /* 0D_NOT_type out */
);
CartesianMapStruct create_wiggler_cartesian_map(EleStruct &ele);
extern "C" void fortran_crystal_attribute_bookkeeper(void *ele /* 0D_NOT_type inout */);
void crystal_attribute_bookkeeper(EleStruct &ele);

// Skipped unusable routine crystal_diffraction_field_calc:
// - Untranslated type: crystal_param_struct (0D)
extern "C" void fortran_crystal_h_misalign(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &h_vec /* 1D_NOT_real inout */
);
void crystal_h_misalign(EleStruct &ele, CoordStruct &orbit, FixedArray1D<Real, 3> h_vec);
extern "C" void fortran_crystal_type_to_crystal_params(
    void *ele /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool crystal_type_to_crystal_params(EleStruct &ele);

// Skipped unusable routine csr3d_steady_state_solver:
// - Variable in sized array: 3D_NOT_real
// - Variable inout sized array: 4D_NOT_real

// Skipped unusable routine csr_and_sc_apply_kicks:
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine csr_bin_kicks:
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine csr_bin_particles:
// - Untranslated type: csr_struct (0D)
extern "C" bool fortran_custom_attribute_ubound_index(
    int &ele_class /* 0D_NOT_integer in */,
    int &ix_ubound /* 0D_NOT_integer out */
);
int custom_attribute_ubound_index(int ele_class);

// Skipped unusable routine custom_ele_attrib_name_list:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine damap_equal_bmad_taylor:
// - Untranslated type: damap (0D)
extern "C" bool fortran_damping_matrix_d(
    double &gamma /* 0D_NOT_real in */,
    double &g_tot /* 0D_NOT_real in */,
    double &B0 /* 0D_NOT_real in */,
    double &B1 /* 0D_NOT_real in */,
    double &delta /* 0D_NOT_real in */,
    int &species /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &mat /* 2D_NOT_real inout */
);
void damping_matrix_d(
    double gamma,
    double g_tot,
    double B0,
    double B1,
    double delta,
    int species,
    FixedArray2D<Real, 6, 6> mat
);

// Skipped unusable routine deallocate_ele_array_pointers:
// - Routine in configuration skip list
extern "C" void fortran_deallocate_ele_pointers(
    void *ele /* 0D_NOT_type inout */,
    bool *nullify_only /* 0D_NOT_logical in */,
    bool *nullify_branch /* 0D_NOT_logical in */,
    bool *dealloc_poles /* 0D_NOT_logical in */
);
void deallocate_ele_pointers(
    EleStruct &ele,
    std::optional<bool> nullify_only = std::nullopt,
    std::optional<bool> nullify_branch = std::nullopt,
    std::optional<bool> dealloc_poles = std::nullopt
);
extern "C" void fortran_deallocate_expression_tree(void *tree /* 0D_NOT_type inout */);
void deallocate_expression_tree(ExpressionTreeStruct &tree);
extern "C" void fortran_deallocate_lat_pointers(void *lat /* 0D_NOT_type inout */);
void deallocate_lat_pointers(LatStruct &lat);
extern "C" bool fortran_default_tracking_species(
    void *param /* 0D_NOT_type in */,
    int &species /* 0D_NOT_integer out */
);
int default_tracking_species(LatParamStruct &param);

// Skipped unusable routine deposit_particles:
// - Untranslated type: mesh3d_struct (0D)
extern "C" bool fortran_detector_pixel_pt(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &ix_pix /* 1D_NOT_integer out */
);
FixedArray1D<Int, 2> detector_pixel_pt(CoordStruct &orbit, EleStruct &ele);
extern "C" bool fortran_diffraction_plate_or_mask_hit_spot(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    int &ix_section /* 0D_NOT_integer out */
);
int diffraction_plate_or_mask_hit_spot(EleStruct &ele, CoordStruct &orbit);
extern "C" bool fortran_diffusion_matrix_b(
    double &gamma /* 0D_NOT_real in */,
    double &g_tot /* 0D_NOT_real in */,
    int &species /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &mat /* 2D_NOT_real inout */
);
void diffusion_matrix_b(double gamma, double g_tot, int species, FixedArray2D<Real, 6, 6> mat);
extern "C" bool fortran_distance_to_aperture(
    void *orbit /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    bool &no_aperture_here /* 0D_NOT_logical out */,
    double &dist /* 0D_NOT_real out */
);
struct DistanceToAperture {
  bool no_aperture_here;
  double dist;
};
Bmad::DistanceToAperture distance_to_aperture(CoordStruct &orbit, int particle_at, EleStruct &ele);

// Skipped unusable routine distance_to_aperture_custom_def:
// - Routine in configuration skip list
extern "C" void
fortran_do_mode_flip(void *ele /* 0D_NOT_type inout */, bool &err_flag /* 0D_NOT_logical out */);
bool do_mode_flip(EleStruct &ele);
extern "C" bool fortran_dpc_given_de(
    double &pc_old /* 0D_NOT_real in */,
    double &mass /* 0D_NOT_real in */,
    double &dE /* 0D_NOT_real in */,
    double &dpc /* 0D_NOT_real in */
);
void dpc_given_de(double pc_old, double mass, double dE, double dpc);
extern "C" void fortran_drift_and_pipe_track_methods_adjustment(void *lat /* 0D_NOT_type inout */);
void drift_and_pipe_track_methods_adjustment(LatStruct &lat);
extern "C" void fortran_drift_multipass_name_correction(void *lat /* 0D_NOT_type inout */);
void drift_multipass_name_correction(LatStruct &lat);
extern "C" void fortran_drift_orbit_time(
    void *orbit /* 0D_NOT_type inout */,
    double &beta0 /* 0D_NOT_real in */,
    double *delta_s /* 0D_NOT_real in */,
    double *delta_t /* 0D_NOT_real in */
);
void drift_orbit_time(
    CoordStruct &orbit,
    double beta0,
    std::optional<double> delta_s = std::nullopt,
    std::optional<double> delta_t = std::nullopt
);
extern "C" void fortran_drift_particle_to_s(
    void *p /* 0D_NOT_type inout */,
    double &s /* 0D_NOT_real in */,
    void *branch /* 0D_NOT_type in */
);
void drift_particle_to_s(CoordStruct &p, double s, BranchStruct &branch);
extern "C" void fortran_drift_particle_to_t(
    void *p /* 0D_NOT_type inout */,
    double &t /* 0D_NOT_real in */,
    void *branch /* 0D_NOT_type in */
);
void drift_particle_to_t(CoordStruct &p, double t, BranchStruct &branch);
extern "C" bool fortran_dspline_len(
    double &s_chord0 /* 0D_NOT_real in */,
    double &s_chord1 /* 0D_NOT_real in */,
    void *spline /* 0D_NOT_type in */,
    double *dtheta_ref /* 0D_NOT_real in */,
    double &dlen /* 0D_NOT_real out */
);
double dspline_len(
    double s_chord0,
    double s_chord1,
    SplineStruct &spline,
    std::optional<double> dtheta_ref = std::nullopt
);
extern "C" void fortran_dynamic_aperture_point(
    void *branch /* 0D_NOT_type in */,
    void *ele0 /* 0D_NOT_type in */,
    void *orb0 /* 0D_NOT_type in */,
    double &theta_xy /* 0D_NOT_real in */,
    void *ap_param /* 0D_NOT_type in */,
    void *ap_point /* 0D_NOT_type out */,
    bool *check_xy_init /* 0D_NOT_logical in */
);
AperturePointStruct dynamic_aperture_point(
    BranchStruct &branch,
    EleStruct &ele0,
    CoordStruct &orb0,
    double theta_xy,
    ApertureParamStruct &ap_param,
    std::optional<bool> check_xy_init = std::nullopt
);
extern "C" void fortran_dynamic_aperture_scan(
    void *aperture_scan /* 1D_ALLOC_type out */,
    void *aperture_param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &pz_start /* 1D_NOT_real in */,
    void *lat /* 0D_NOT_type in */,
    bool *print_timing /* 0D_NOT_logical in */
);
ApertureScanStructAlloc1D dynamic_aperture_scan(
    ApertureParamStruct &aperture_param,
    FArray1D<Real> &pz_start,
    LatStruct &lat,
    std::optional<bool> print_timing = std::nullopt
);
extern "C" bool fortran_e_accel_field(
    void *ele /* 0D_NOT_type in */,
    int &voltage_or_gradient /* 0D_NOT_integer in */,
    bool *bmad_standard_tracking /* 0D_NOT_logical in */,
    double &field /* 0D_NOT_real out */
);
double e_accel_field(
    EleStruct &ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking = std::nullopt
);
extern "C" bool fortran_e_crit_photon(
    double &gamma /* 0D_NOT_real in */,
    double &g_bend /* 0D_NOT_real in */,
    double &E_crit /* 0D_NOT_real out */
);
double e_crit_photon(double gamma, double g_bend);
extern "C" void fortran_eigen_decomp_6mat(
    Bmad::array_descriptor_t &mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &eval /* 1D_NOT_complex out */,
    Bmad::array_descriptor_t &evec /* 2D_NOT_complex out */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &tunes /* 1D_NOT_real out */
);
struct EigenDecomp6mat {
  FixedArray1D<Complex, 6> eval;
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
  FixedArray1D<Real, 3> tunes;
};
Bmad::EigenDecomp6mat eigen_decomp_6mat(FixedArray2D<Real, 6, 6> mat);

// Skipped unusable routine eigensys:
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_ele_compute_ref_energy_and_time(
    void *ele0 /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void ele_compute_ref_energy_and_time(
    EleStruct &ele0,
    EleStruct &ele,
    LatParamStruct &param,
    bool err_flag
);
extern "C" void
fortran_ele_equal_ele(void *ele_out /* 0D_NOT_type inout */, void *ele_in /* 0D_NOT_type in */);
void ele_equal_ele(EleStruct &ele_out, EleStruct &ele_in);
extern "C" void fortran_ele_equals_ele(
    void *ele_out /* 0D_NOT_type out */,
    void *ele_in /* 0D_NOT_type in */,
    bool &update_nametable /* 0D_NOT_logical in */
);
EleStruct ele_equals_ele(EleStruct &ele_in, bool update_nametable);
extern "C" void fortran_ele_finalizer(void *ele /* 0D_NOT_type inout */);
void ele_finalizer(EleStruct &ele);
extern "C" bool fortran_ele_full_name(
    void *ele /* 0D_NOT_type in */,
    const char *template_ /* 0D_NOT_character in */,
    const char *str /* 0D_ALLOC_character out */
);
std::string ele_full_name(EleStruct &ele, std::optional<std::string> template_ = std::nullopt);
extern "C" void fortran_ele_geometry(
    void *floor_start /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *floor_end /* 0D_NOT_type out */,
    double *len_scale /* 0D_NOT_real in */,
    bool *ignore_patch_err /* 0D_NOT_logical in */
);
FloorPositionStruct ele_geometry(
    FloorPositionStruct &floor_start,
    EleStruct &ele,
    std::optional<double> len_scale = std::nullopt,
    std::optional<bool> ignore_patch_err = std::nullopt
);

// Skipped unusable routine ele_geometry_hook_def:
// - Routine in configuration skip list
extern "C" bool fortran_ele_geometry_with_misalignments(
    void *ele /* 0D_NOT_type in */,
    double *len_scale /* 0D_NOT_real in */,
    void *floor /* 0D_NOT_type out */
);
FloorPositionStruct
ele_geometry_with_misalignments(EleStruct &ele, std::optional<double> len_scale = std::nullopt);
extern "C" bool fortran_ele_has_constant_ds_dt_ref(
    void *ele /* 0D_NOT_type in */,
    bool &is_const /* 0D_NOT_logical out */
);
bool ele_has_constant_ds_dt_ref(EleStruct &ele);
extern "C" bool fortran_ele_has_nonzero_kick(
    void *ele /* 0D_NOT_type inout */,
    bool &has_kick /* 0D_NOT_logical in */
);
void ele_has_nonzero_kick(EleStruct &ele, bool has_kick);
extern "C" bool fortran_ele_has_nonzero_offset(
    void *ele /* 0D_NOT_type in */,
    bool &has_offset /* 0D_NOT_logical out */
);
bool ele_has_nonzero_offset(EleStruct &ele);
extern "C" bool fortran_ele_is_monitor(
    void *ele /* 0D_NOT_type in */,
    bool *print_warning /* 0D_NOT_logical in */,
    bool &is_monitor /* 0D_NOT_logical out */
);
bool ele_is_monitor(EleStruct &ele, std::optional<bool> print_warning = std::nullopt);
extern "C" bool fortran_ele_loc(void *ele /* 0D_NOT_type in */, void *loc /* 0D_NOT_type out */);
LatEleLocStruct ele_loc(EleStruct &ele);
extern "C" bool fortran_ele_loc_name(
    void *ele /* 0D_NOT_type in */,
    bool *show_branch0 /* 0D_NOT_logical in */,
    const char *parens /* 0D_NOT_character in */,
    const char *str /* 0D_NOT_character out */
);
std::string ele_loc_name(
    EleStruct &ele,
    std::optional<bool> show_branch0 = std::nullopt,
    std::optional<std::string> parens = std::nullopt
);
extern "C" void fortran_ele_misalignment_l_s_calc(
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &L_mis /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &S_mis /* 2D_NOT_real out */
);
struct EleMisalignmentLSCalc {
  FixedArray1D<Real, 3> L_mis;
  FixedArray2D<Real, 3, 3> S_mis;
};
Bmad::EleMisalignmentLSCalc ele_misalignment_l_s_calc(EleStruct &ele);
extern "C" bool
fortran_ele_nametable_index(void *ele /* 0D_NOT_type in */, int &ix_nt /* 0D_NOT_integer out */);
int ele_nametable_index(EleStruct &ele);
extern "C" void
fortran_ele_order_calc(void *lat /* 0D_NOT_type in */, void *order /* 0D_NOT_type out */);
LatEleOrderStruct ele_order_calc(LatStruct &lat);
extern "C" void fortran_ele_reference_energy_correction(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    int &particle_at /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void ele_reference_energy_correction(
    EleStruct &ele,
    CoordStruct &orbit,
    int particle_at,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_ele_rf_step_index(
    double &E_ref /* 0D_NOT_real in */,
    double &s_rel /* 0D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    int &ix_step /* 0D_NOT_integer out */
);
int ele_rf_step_index(double E_ref, double s_rel, EleStruct &ele);
extern "C" void fortran_ele_to_fibre(
    void *ele /* 0D_NOT_type in */,
    void *ptc_fibre /* 0D_PTR_type out */,
    bool &use_offsets /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    int *integ_order /* 0D_NOT_integer in */,
    int *steps /* 0D_NOT_integer in */,
    bool *for_layout /* 0D_NOT_logical in */,
    void *ref_in /* 0D_NOT_type in */
);
struct EleToFibre {
  std::optional<Fibre> ptc_fibre;
  bool err_flag;
};
Bmad::EleToFibre ele_to_fibre(
    EleStruct &ele,
    bool use_offsets,
    std::optional<int> integ_order = std::nullopt,
    std::optional<int> steps = std::nullopt,
    std::optional<bool> for_layout = std::nullopt,
    optional_ref<CoordStruct> ref_in = std::nullopt
);

// Skipped unusable routine ele_to_fibre_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_ele_to_ptc_magnetic_bn_an(
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &bn /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &an /* 1D_NOT_real inout */,
    int &n_max /* 0D_NOT_integer out */
);
int ele_to_ptc_magnetic_bn_an(EleStruct &ele, FArray1D<Real> &bn, FArray1D<Real> &an);
extern "C" void fortran_ele_to_spin_taylor(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *orb0 /* 0D_NOT_type in */
);
void ele_to_spin_taylor(EleStruct &ele, LatParamStruct &param, CoordStruct &orb0);
extern "C" void fortran_ele_to_taylor(
    void *ele /* 0D_NOT_type in */,
    void *orb0 /* 0D_NOT_type in */,
    bool *taylor_map_includes_offsets /* 0D_NOT_logical in */,
    bool *include_damping /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &orbital_taylor /* 1D_NOT_type out */,
    Bmad::array_descriptor_t &spin_taylor /* 1D_NOT_type out */
);
struct EleToTaylor {
  TaylorStructArray1D orbital_taylor;
  TaylorStructArray1D spin_taylor;
};
Bmad::EleToTaylor ele_to_taylor(
    EleStruct &ele,
    optional_ref<CoordStruct> orb0 = std::nullopt,
    std::optional<bool> taylor_map_includes_offsets = std::nullopt,
    std::optional<bool> include_damping = std::nullopt
);
extern "C" bool fortran_ele_unique_name(
    void *ele /* 0D_NOT_type in */,
    void *order /* 0D_NOT_type in */,
    const char *unique_name /* 0D_NOT_character out */
);
std::string ele_unique_name(EleStruct &ele, LatEleOrderStruct &order);
extern "C" bool fortran_ele_value_has_changed(
    void *ele /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &list /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &abs_tol /* 1D_NOT_real in */,
    bool &set_old /* 0D_NOT_logical in */,
    bool &has_changed /* 0D_NOT_logical out */
);
bool ele_value_has_changed(
    EleStruct &ele,
    FArray1D<Int> &list,
    FArray1D<Real> &abs_tol,
    bool set_old
);
extern "C" void fortran_ele_vec_equal_ele_vec(
    Bmad::array_descriptor_t &ele1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &ele2 /* 1D_NOT_type in */
);
void ele_vec_equal_ele_vec(EleStructArray1D ele1, EleStructArray1D ele2);
extern "C" void fortran_elec_multipole_field(
    double &a /* 0D_NOT_real in */,
    double &b /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    void *coord /* 0D_NOT_type in */,
    double &Ex /* 0D_NOT_real out */,
    double &Ey /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &dE /* 2D_NOT_real out */,
    bool &compute_dE /* 0D_NOT_logical out */
);
struct ElecMultipoleField {
  double Ex;
  double Ey;
  std::optional<FixedArray2D<Real, 2, 2>> dE;
  bool compute_dE;
};
Bmad::ElecMultipoleField elec_multipole_field(double a, double b, int n, CoordStruct &coord);
extern "C" bool fortran_element_at_s_branch(
    void *branch /* 0D_NOT_type in */,
    double &s /* 0D_NOT_real in */,
    bool &choose_max /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &s_eff /* 0D_NOT_real out */,
    void *position /* 0D_NOT_type out */,
    bool *print_err /* 0D_NOT_logical in */,
    int &ix_ele /* 0D_NOT_integer out */
);
struct ElementAtSBranch {
  bool err_flag;
  double s_eff;
  CoordStruct position;
  int ix_ele;
};
Bmad::ElementAtSBranch element_at_s(
    BranchStruct &branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool fortran_element_at_s_lat(
    void *lat /* 0D_NOT_type in */,
    double &s /* 0D_NOT_real in */,
    bool &choose_max /* 0D_NOT_logical in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &s_eff /* 0D_NOT_real out */,
    void *position /* 0D_NOT_type out */,
    bool *print_err /* 0D_NOT_logical in */,
    int &ix_ele /* 0D_NOT_integer out */
);
struct ElementAtSLat {
  bool err_flag;
  double s_eff;
  CoordStruct position;
  int ix_ele;
};
Bmad::ElementAtSLat element_at_s(
    LatStruct &lat,
    double s,
    bool choose_max,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_element_slice_iterator(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &i_slice /* 0D_NOT_integer in */,
    int &n_slice_tot /* 0D_NOT_integer in */,
    void *sliced_ele /* 0D_NOT_type inout */,
    double *s_start /* 0D_NOT_real in */,
    double *s_end /* 0D_NOT_real in */
);
void element_slice_iterator(
    EleStruct &ele,
    LatParamStruct &param,
    int i_slice,
    int n_slice_tot,
    EleStruct &sliced_ele,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt
);
extern "C" void fortran_ellipinc_test();
void ellipinc_test();
extern "C" void fortran_em_field_calc(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &s_pos /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    bool &local_ref_frame /* 0D_NOT_logical in */,
    void *field /* 0D_NOT_type out */,
    bool *calc_dfield /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *calc_potential /* 0D_NOT_logical in */,
    bool *use_overlap /* 0D_NOT_logical in */,
    bool *grid_allow_s_out_of_bounds /* 0D_NOT_logical in */,
    double *rf_time /* 0D_NOT_real in */,
    void *used_eles /* 1D_ALLOC_type in */,
    bool *print_err /* 0D_NOT_logical in */,
    void *original_ele /* 0D_NOT_type in */
);
struct EmFieldCalc {
  EmFieldStruct field;
  bool err_flag;
};
Bmad::EmFieldCalc em_field_calc(
    EleStruct &ele,
    LatParamStruct &param,
    double s_pos,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt,
    std::optional<bool> use_overlap = std::nullopt,
    std::optional<bool> grid_allow_s_out_of_bounds = std::nullopt,
    std::optional<double> rf_time = std::nullopt,
    std::optional<ElePointerStructAlloc1D> used_eles = std::nullopt,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EleStruct> original_ele = std::nullopt
);

// Skipped unusable routine em_field_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_em_field_derivatives(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &s_pos /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    bool &local_ref_frame /* 0D_NOT_logical in */,
    void *dfield /* 0D_NOT_type out */,
    bool *grid_allow_s_out_of_bounds /* 0D_NOT_logical in */,
    double *rf_time /* 0D_NOT_real in */
);
EmFieldStruct em_field_derivatives(
    EleStruct &ele,
    LatParamStruct &param,
    double s_pos,
    CoordStruct &orbit,
    bool local_ref_frame,
    std::optional<bool> grid_allow_s_out_of_bounds = std::nullopt,
    std::optional<double> rf_time = std::nullopt
);
extern "C" void fortran_em_field_kick_vector_time(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &rf_time /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &dvec_dt /* 1D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */,
    void *extra_field /* 0D_NOT_type in */
);
FixedArray1D<Real, 10> em_field_kick_vector_time(
    EleStruct &ele,
    LatParamStruct &param,
    double rf_time,
    CoordStruct &orbit,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EmFieldStruct> extra_field = std::nullopt
);
extern "C" bool fortran_em_field_plus_em_field(
    void *field1 /* 0D_NOT_type in */,
    void *field2 /* 0D_NOT_type in */,
    void *field_tot /* 0D_NOT_type inout */
);
void em_field_plus_em_field(EmFieldStruct &field1, EmFieldStruct &field2, EmFieldStruct &field_tot);
extern "C" void fortran_emit_6d(
    void *ele_ref /* 0D_NOT_type in */,
    bool &include_opening_angle /* 0D_NOT_logical in */,
    void *mode /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &closed_orbit /* 1D_NOT_type in */,
    void *rad_int_by_ele /* 0D_NOT_type out */
);
struct Emit6d {
  NormalModesStruct mode;
  FixedArray2D<Real, 6, 6> sigma_mat;
  RadIntAllEleStruct rad_int_by_ele;
};
Bmad::Emit6d emit_6d(
    EleStruct &ele_ref,
    bool include_opening_angle,
    std::optional<CoordStructArray1D> closed_orbit = std::nullopt
);
extern "C" bool fortran_entering_element(
    void *orbit /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    bool &is_entering /* 0D_NOT_logical out */
);
bool entering_element(CoordStruct &orbit, int particle_at);
extern "C" void fortran_envelope_radints(
    Bmad::array_descriptor_t &Lambda /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &Theta /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &Iota /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &alpha /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &emit /* 1D_NOT_real inout */
);
void envelope_radints(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    FixedArray1D<Real, 3> alpha,
    FixedArray1D<Real, 3> emit
);
extern "C" void fortran_envelope_radints_ibs(
    Bmad::array_descriptor_t &Lambda /* 2D_NOT_complex in */,
    Bmad::array_descriptor_t &Theta /* 2D_NOT_complex in */,
    Bmad::array_descriptor_t &Iota /* 2D_NOT_complex in */,
    Bmad::array_descriptor_t &eles /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &alpha /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &emit /* 1D_NOT_real out */,
    void *mode /* 0D_NOT_type in */,
    bool &tail_cut /* 0D_NOT_logical in */,
    double &npart /* 0D_NOT_real in */,
    int &species /* 0D_NOT_integer in */
);
struct EnvelopeRadintsIbs {
  FixedArray1D<Real, 3> alpha;
  FixedArray1D<Real, 3> emit;
};
Bmad::EnvelopeRadintsIbs envelope_radints_ibs(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    EleStructArray1D eles,
    NormalModesStruct &mode,
    bool tail_cut,
    double npart,
    int species
);
extern "C" bool fortran_eq_ac_kicker(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ac_kicker(AcKickerStruct &f1, AcKickerStruct &f2, bool is_eq);
extern "C" bool fortran_eq_ac_kicker_freq(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ac_kicker_freq(AcKickerFreqStruct &f1, AcKickerFreqStruct &f2, bool is_eq);
extern "C" bool fortran_eq_ac_kicker_time(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ac_kicker_time(AcKickerTimeStruct &f1, AcKickerTimeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_anormal_mode(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_anormal_mode(AnormalModeStruct &f1, AnormalModeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_aperture_param(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_aperture_param(ApertureParamStruct &f1, ApertureParamStruct &f2, bool is_eq);
extern "C" bool fortran_eq_aperture_point(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_aperture_point(AperturePointStruct &f1, AperturePointStruct &f2, bool is_eq);
extern "C" bool fortran_eq_aperture_scan(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_aperture_scan(ApertureScanStruct &f1, ApertureScanStruct &f2, bool is_eq);
extern "C" bool fortran_eq_beam(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_beam(BeamStruct &f1, BeamStruct &f2, bool is_eq);
extern "C" bool fortran_eq_beam_init(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_beam_init(BeamInitStruct &f1, BeamInitStruct &f2, bool is_eq);
extern "C" bool fortran_eq_bmad_common(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_bmad_common(BmadCommonStruct &f1, BmadCommonStruct &f2, bool is_eq);
extern "C" bool fortran_eq_bookkeeping_state(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_bookkeeping_state(BookkeepingStateStruct &f1, BookkeepingStateStruct &f2, bool is_eq);
extern "C" bool fortran_eq_bpm_phase_coupling(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_bpm_phase_coupling(BpmPhaseCouplingStruct &f1, BpmPhaseCouplingStruct &f2, bool is_eq);
extern "C" bool fortran_eq_branch(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_branch(BranchStruct &f1, BranchStruct &f2, bool is_eq);
extern "C" bool fortran_eq_bunch(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_bunch(BunchStruct &f1, BunchStruct &f2, bool is_eq);
extern "C" bool fortran_eq_bunch_params(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_bunch_params(BunchParamsStruct &f1, BunchParamsStruct &f2, bool is_eq);
extern "C" bool fortran_eq_cartesian_map(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cartesian_map(CartesianMapStruct &f1, CartesianMapStruct &f2, bool is_eq);
extern "C" bool fortran_eq_cartesian_map_term(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cartesian_map_term(CartesianMapTermStruct &f1, CartesianMapTermStruct &f2, bool is_eq);
extern "C" bool fortran_eq_cartesian_map_term1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cartesian_map_term1(CartesianMapTerm1Struct &f1, CartesianMapTerm1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_complex_taylor(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_complex_taylor(ComplexTaylorStruct &f1, ComplexTaylorStruct &f2, bool is_eq);
extern "C" bool fortran_eq_complex_taylor_term(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_complex_taylor_term(ComplexTaylorTermStruct &f1, ComplexTaylorTermStruct &f2, bool is_eq);
extern "C" bool fortran_eq_control(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_control(ControlStruct &f1, ControlStruct &f2, bool is_eq);
extern "C" bool fortran_eq_control_ramp1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_control_ramp1(ControlRamp1Struct &f1, ControlRamp1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_control_var1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_control_var1(ControlVar1Struct &f1, ControlVar1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_controller(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_controller(ControllerStruct &f1, ControllerStruct &f2, bool is_eq);
extern "C" bool fortran_eq_coord(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_coord(CoordStruct &f1, CoordStruct &f2, bool is_eq);
extern "C" bool fortran_eq_coord_array(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_coord_array(CoordArrayStruct &f1, CoordArrayStruct &f2, bool is_eq);
extern "C" bool fortran_eq_cylindrical_map(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cylindrical_map(CylindricalMapStruct &f1, CylindricalMapStruct &f2, bool is_eq);
extern "C" bool fortran_eq_cylindrical_map_term(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cylindrical_map_term(
    CylindricalMapTermStruct &f1,
    CylindricalMapTermStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_cylindrical_map_term1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_cylindrical_map_term1(
    CylindricalMapTerm1Struct &f1,
    CylindricalMapTerm1Struct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_ele(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ele(EleStruct &f1, EleStruct &f2, bool is_eq);
extern "C" bool fortran_eq_ellipse_beam_init(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ellipse_beam_init(EllipseBeamInitStruct &f1, EllipseBeamInitStruct &f2, bool is_eq);
extern "C" bool fortran_eq_em_field(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_em_field(EmFieldStruct &f1, EmFieldStruct &f2, bool is_eq);
extern "C" bool fortran_eq_expression_atom(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_expression_atom(ExpressionAtomStruct &f1, ExpressionAtomStruct &f2, bool is_eq);
extern "C" bool fortran_eq_floor_position(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_floor_position(FloorPositionStruct &f1, FloorPositionStruct &f2, bool is_eq);
extern "C" bool fortran_eq_gen_grad1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_gen_grad1(GenGrad1Struct &f1, GenGrad1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_gen_grad_map(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_gen_grad_map(GenGradMapStruct &f1, GenGradMapStruct &f2, bool is_eq);
extern "C" bool fortran_eq_gg_taylor(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_gg_taylor(GgTaylorStruct &f1, GgTaylorStruct &f2, bool is_eq);
extern "C" bool fortran_eq_gg_taylor_term(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_gg_taylor_term(GgTaylorTermStruct &f1, GgTaylorTermStruct &f2, bool is_eq);
extern "C" bool fortran_eq_grid_beam_init(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_grid_beam_init(GridBeamInitStruct &f1, GridBeamInitStruct &f2, bool is_eq);
extern "C" bool fortran_eq_grid_field(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_grid_field(GridFieldStruct &f1, GridFieldStruct &f2, bool is_eq);
extern "C" bool fortran_eq_grid_field_pt(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_grid_field_pt(GridFieldPtStruct &f1, GridFieldPtStruct &f2, bool is_eq);
extern "C" bool fortran_eq_grid_field_pt1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_grid_field_pt1(GridFieldPt1Struct &f1, GridFieldPt1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_high_energy_space_charge(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_high_energy_space_charge(
    HighEnergySpaceChargeStruct &f1,
    HighEnergySpaceChargeStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_interval1_coef(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_interval1_coef(Interval1CoefStruct &f1, Interval1CoefStruct &f2, bool is_eq);
extern "C" bool fortran_eq_kv_beam_init(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_kv_beam_init(KvBeamInitStruct &f1, KvBeamInitStruct &f2, bool is_eq);
extern "C" bool fortran_eq_lat(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_lat(LatStruct &f1, LatStruct &f2, bool is_eq);
extern "C" bool fortran_eq_lat_ele_loc(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_lat_ele_loc(LatEleLocStruct &f1, LatEleLocStruct &f2, bool is_eq);
extern "C" bool fortran_eq_lat_param(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_lat_param(LatParamStruct &f1, LatParamStruct &f2, bool is_eq);
extern "C" bool fortran_eq_linac_normal_mode(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_linac_normal_mode(LinacNormalModeStruct &f1, LinacNormalModeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_mode3(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_mode3(Mode3Struct &f1, Mode3Struct &f2, bool is_eq);
extern "C" bool fortran_eq_mode_info(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_mode_info(ModeInfoStruct &f1, ModeInfoStruct &f2, bool is_eq);
extern "C" bool fortran_eq_normal_modes(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_normal_modes(NormalModesStruct &f1, NormalModesStruct &f2, bool is_eq);
extern "C" bool fortran_eq_photon_element(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_photon_element(PhotonElementStruct &f1, PhotonElementStruct &f2, bool is_eq);
extern "C" bool fortran_eq_photon_material(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_photon_material(PhotonMaterialStruct &f1, PhotonMaterialStruct &f2, bool is_eq);
extern "C" bool fortran_eq_photon_reflect_surface(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_photon_reflect_surface(
    PhotonReflectSurfaceStruct &f1,
    PhotonReflectSurfaceStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_photon_reflect_table(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_photon_reflect_table(
    PhotonReflectTableStruct &f1,
    PhotonReflectTableStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_photon_target(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_photon_target(PhotonTargetStruct &f1, PhotonTargetStruct &f2, bool is_eq);
extern "C" bool fortran_eq_pixel_detec(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_pixel_detec(PixelDetecStruct &f1, PixelDetecStruct &f2, bool is_eq);
extern "C" bool fortran_eq_pixel_pt(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_pixel_pt(PixelPtStruct &f1, PixelPtStruct &f2, bool is_eq);
extern "C" bool fortran_eq_pre_tracker(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_pre_tracker(PreTrackerStruct &f1, PreTrackerStruct &f2, bool is_eq);
extern "C" bool fortran_eq_rad_int1(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_rad_int1(RadInt1Struct &f1, RadInt1Struct &f2, bool is_eq);
extern "C" bool fortran_eq_rad_int_all_ele(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_rad_int_all_ele(RadIntAllEleStruct &f1, RadIntAllEleStruct &f2, bool is_eq);
extern "C" bool fortran_eq_rad_int_branch(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_rad_int_branch(RadIntBranchStruct &f1, RadIntBranchStruct &f2, bool is_eq);
extern "C" bool fortran_eq_rad_map(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_rad_map(RadMapStruct &f1, RadMapStruct &f2, bool is_eq);
extern "C" bool fortran_eq_rad_map_ele(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_rad_map_ele(RadMapEleStruct &f1, RadMapEleStruct &f2, bool is_eq);
extern "C" bool fortran_eq_ramper_lord(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_ramper_lord(RamperLordStruct &f1, RamperLordStruct &f2, bool is_eq);
extern "C" bool fortran_eq_space_charge_common(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_space_charge_common(SpaceChargeCommonStruct &f1, SpaceChargeCommonStruct &f2, bool is_eq);
extern "C" bool fortran_eq_spin_polar(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_spin_polar(SpinPolarStruct &f1, SpinPolarStruct &f2, bool is_eq);
extern "C" bool fortran_eq_spline(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_spline(SplineStruct &f1, SplineStruct &f2, bool is_eq);
extern "C" bool fortran_eq_strong_beam(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_strong_beam(StrongBeamStruct &f1, StrongBeamStruct &f2, bool is_eq);
extern "C" bool fortran_eq_surface_curvature(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_curvature(SurfaceCurvatureStruct &f1, SurfaceCurvatureStruct &f2, bool is_eq);
extern "C" bool fortran_eq_surface_displacement(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_displacement(
    SurfaceDisplacementStruct &f1,
    SurfaceDisplacementStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_surface_displacement_pt(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_displacement_pt(
    SurfaceDisplacementPtStruct &f1,
    SurfaceDisplacementPtStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_surface_h_misalign(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_h_misalign(SurfaceHMisalignStruct &f1, SurfaceHMisalignStruct &f2, bool is_eq);
extern "C" bool fortran_eq_surface_h_misalign_pt(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_h_misalign_pt(
    SurfaceHMisalignPtStruct &f1,
    SurfaceHMisalignPtStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_surface_segmented(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_segmented(SurfaceSegmentedStruct &f1, SurfaceSegmentedStruct &f2, bool is_eq);
extern "C" bool fortran_eq_surface_segmented_pt(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_surface_segmented_pt(
    SurfaceSegmentedPtStruct &f1,
    SurfaceSegmentedPtStruct &f2,
    bool is_eq
);
extern "C" bool fortran_eq_target_point(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_target_point(TargetPointStruct &f1, TargetPointStruct &f2, bool is_eq);
extern "C" bool fortran_eq_taylor(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_taylor(TaylorStruct &f1, TaylorStruct &f2, bool is_eq);
extern "C" bool fortran_eq_taylor_term(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_taylor_term(TaylorTermStruct &f1, TaylorTermStruct &f2, bool is_eq);
extern "C" bool fortran_eq_track(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_track(TrackStruct &f1, TrackStruct &f2, bool is_eq);
extern "C" bool fortran_eq_track_point(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_track_point(TrackPointStruct &f1, TrackPointStruct &f2, bool is_eq);
extern "C" bool fortran_eq_twiss(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_twiss(TwissStruct &f1, TwissStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake(WakeStruct &f1, WakeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake_lr(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake_lr(WakeLrStruct &f1, WakeLrStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake_lr_mode(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake_lr_mode(WakeLrModeStruct &f1, WakeLrModeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake_sr(WakeSrStruct &f1, WakeSrStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr_mode(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake_sr_mode(WakeSrModeStruct &f1, WakeSrModeStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr_z_long(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wake_sr_z_long(WakeSrZLongStruct &f1, WakeSrZLongStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wall3d(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wall3d(Wall3dStruct &f1, Wall3dStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wall3d_section(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wall3d_section(Wall3dSectionStruct &f1, Wall3dSectionStruct &f2, bool is_eq);
extern "C" bool fortran_eq_wall3d_vertex(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_wall3d_vertex(Wall3dVertexStruct &f1, Wall3dVertexStruct &f2, bool is_eq);
extern "C" bool fortran_eq_xy_disp(
    void *f1 /* 0D_NOT_type in */,
    void *f2 /* 0D_NOT_type in */,
    bool &is_eq /* 0D_NOT_logical in */
);
void eq_xy_disp(XyDispStruct &f1, XyDispStruct &f2, bool is_eq);
extern "C" bool fortran_equal_sign_here(
    void *ele /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &is_here /* 0D_NOT_logical in */
);
void equal_sign_here(EleStruct &ele, std::string delim, bool is_here);
extern "C" bool fortran_equivalent_taylor_attributes(
    void *ele_taylor /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool &equiv /* 0D_NOT_logical out */
);
bool equivalent_taylor_attributes(EleStruct &ele_taylor, EleStruct &ele2);
extern "C" void fortran_etdiv(
    double &A /* 0D_NOT_real in */,
    double &B /* 0D_NOT_real in */,
    double &C /* 0D_NOT_real in */,
    double &D /* 0D_NOT_real in */,
    double &E /* 0D_NOT_real in */,
    double &F /* 0D_NOT_real in */
);
void etdiv(double A, double B, double C, double D, double E, double F);

// Skipped unusable routine ety:
// - Array bounds handling: "Enum 'NM' found in bounds 'NM' but not in provided map."
// - Array bounds handling: "Enum 'IGH' found in bounds 'IGH' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine ety2:
// - Array bounds handling: "Enum 'NM' found in bounds 'NM' but not in provided map."
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Array bounds handling: "Enum 'N' found in bounds 'N' but not in provided map."
// - Array bounds handling: "Enum 'NM' found in bounds 'NM' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine etyt:
// - Array bounds handling: "Enum 'NM' found in bounds 'NM' but not in provided map."
// - Array bounds handling: "Enum 'IGH' found in bounds 'IGH' but not in provided map."
// - Array bounds handling: "Enum 'NM' found in bounds 'NM' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_evaluate_array_index(
    bool &err_flag /* 0D_NOT_logical out */,
    const char *delim_list1 /* 0D_NOT_character in */,
    const char *word2 /* 0D_NOT_character out */,
    const char *delim_list2 /* 0D_NOT_character in */,
    const char *delim2 /* 0D_NOT_character out */,
    int &this_index /* 0D_NOT_integer out */
);
struct EvaluateArrayIndex {
  bool err_flag;
  std::string word2;
  std::string delim2;
  int this_index;
};
Bmad::EvaluateArrayIndex evaluate_array_index(std::string delim_list1, std::string delim_list2);
extern "C" bool fortran_evaluate_logical(
    const char *word /* 0D_NOT_character in */,
    int &iostat /* 0D_NOT_integer out */,
    bool &this_logic /* 0D_NOT_logical out */
);
struct EvaluateLogical {
  int iostat;
  bool this_logic;
};
Bmad::EvaluateLogical evaluate_logical(std::string word);
extern "C" void fortran_exact_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void exact_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_exp_bessi0(
    double &t /* 0D_NOT_real in */,
    double &B1 /* 0D_NOT_real in */,
    double &B2 /* 0D_NOT_real in */,
    double &func_retval__ /* 0D_NOT_real in */
);
void exp_bessi0(double t, double B1, double B2, double func_retval__);
extern "C" bool fortran_expect_one_of(
    const char *delim_list /* 0D_NOT_character in */,
    bool &check_input_delim /* 0D_NOT_logical in */,
    const char *ele_name /* 0D_NOT_character in */,
    const char *delim /* 0D_NOT_character inout */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &is_ok /* 0D_NOT_logical in */
);
void expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string &delim,
    bool delim_found,
    bool is_ok
);
extern "C" bool fortran_expect_this(
    const char *expecting /* 0D_NOT_character in */,
    bool &check_delim /* 0D_NOT_logical in */,
    bool &call_check /* 0D_NOT_logical in */,
    const char *err_str /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    bool &is_ok /* 0D_NOT_logical in */
);
struct ExpectThis {
  std::string delim;
  bool delim_found;
};
Bmad::ExpectThis expect_this(
    std::string expecting,
    bool check_delim,
    bool call_check,
    std::string err_str,
    EleStruct &ele,
    bool is_ok
);
extern "C" bool fortran_expression_stack_to_string(
    Bmad::array_descriptor_t &stack /* 1D_NOT_type in */,
    bool *polish /* 0D_NOT_logical in */,
    const char *str /* 0D_ALLOC_character out */
);
std::string expression_stack_to_string(
    ExpressionAtomStructArray1D stack,
    std::optional<bool> polish = std::nullopt
);
extern "C" bool fortran_expression_stack_value(
    Bmad::array_descriptor_t &stack /* 1D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    const char *err_str /* 0D_NOT_character out */,
    Bmad::array_descriptor_t &var /* 1D_NOT_type in */,
    bool *use_old /* 0D_NOT_logical in */,
    double &value /* 0D_NOT_real out */
);
struct ExpressionStackValue {
  bool err_flag;
  std::string err_str;
  double value;
};
Bmad::ExpressionStackValue expression_stack_value(
    ExpressionAtomStructArray1D stack,
    std::optional<ControlVar1StructArray1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt
);
extern "C" void fortran_expression_string_to_stack(
    const char *string /* 0D_NOT_character in */,
    void *stack /* 1D_ALLOC_type out */,
    int &n_stack /* 0D_NOT_integer out */,
    bool &err_flag /* 0D_NOT_logical out */,
    const char *err_str /* 0D_NOT_character out */
);
struct ExpressionStringToStack {
  ExpressionAtomStructAlloc1D stack;
  int n_stack;
  bool err_flag;
  std::string err_str;
};
Bmad::ExpressionStringToStack expression_string_to_stack(std::string string);
extern "C" void fortran_expression_string_to_tree(
    const char *string /* 0D_NOT_character in */,
    void *root_tree /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    const char *err_str /* 0D_NOT_character out */
);
struct ExpressionStringToTree {
  bool err_flag;
  std::string err_str;
};
Bmad::ExpressionStringToTree
expression_string_to_tree(std::string string, ExpressionTreeStruct &root_tree);
extern "C" bool fortran_expression_tree_to_string(
    void *tree /* 0D_NOT_type in */,
    bool *include_root /* 0D_NOT_logical in */,
    int *n_node /* 0D_NOT_integer in */,
    void *parent /* 0D_NOT_type in */,
    const char *str_out /* 0D_ALLOC_character out */
);
std::string expression_tree_to_string(
    ExpressionTreeStruct &tree,
    std::optional<bool> include_root = std::nullopt,
    std::optional<int> n_node = std::nullopt,
    optional_ref<ExpressionTreeStruct> parent = std::nullopt
);
extern "C" bool fortran_expression_value(
    const char *expression /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */,
    const char *err_str /* 0D_NOT_character out */,
    Bmad::array_descriptor_t &var /* 1D_NOT_type in */,
    bool *use_old /* 0D_NOT_logical in */,
    double &value /* 0D_NOT_real out */
);
struct ExpressionValue {
  bool err_flag;
  std::string err_str;
  double value;
};
Bmad::ExpressionValue expression_value(
    std::string expression,
    std::optional<ControlVar1StructArray1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt
);
extern "C" void fortran_fft1(
    Bmad::array_descriptor_t &a /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &b /* 1D_NOT_real inout */,
    int &n /* 0D_NOT_integer in */,
    int &isn /* 0D_NOT_integer in */,
    int &ierr /* 0D_NOT_integer out */
);
int fft1(FArray1D<Real> &a, FArray1D<Real> &b, int n, int isn);

// Skipped unusable routine fftconvcorr3d:
// - Array bounds handling: "Enum 'RILO' found in bounds 'rilo' but not in provided map."
// - Array bounds handling: "Enum 'G1ILO' found in bounds 'g1ilo' but not in provided map."
// - Array bounds handling: "Enum 'G2ILO' found in bounds 'g2ilo' but not in provided map."
// - Array bounds handling: "Enum 'G3ILO' found in bounds 'g3ilo' but not in provided map."
// - Array bounds handling: "Enum 'G4ILO' found in bounds 'g4ilo' but not in provided map."
// - Array bounds handling: "Enum 'CILO' found in bounds 'cilo' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_fibre_to_ele(
    void *ptc_fibre /* 0D_NOT_type in */,
    void *branch /* 0D_NOT_type inout */,
    int &ix_ele /* 0D_NOT_integer inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *from_mad /* 0D_NOT_logical in */
);
bool fibre_to_ele(
    Fibre &ptc_fibre,
    BranchStruct &branch,
    int &ix_ele,
    std::optional<bool> from_mad = std::nullopt
);
extern "C" bool fortran_field_attribute_free(
    void *ele /* 0D_NOT_type in */,
    const char *attrib_name /* 0D_NOT_character in */,
    bool &free /* 0D_NOT_logical out */
);
bool field_attribute_free(EleStruct &ele, std::string attrib_name);

// Skipped unusable routine field_interpolate_3d:
// - Variable in sized array: 3D_NOT_real
extern "C" void fortran_finalize_reflectivity_table(
    void *table /* 0D_NOT_type inout */,
    bool &in_degrees /* 0D_NOT_logical in */
);
void finalize_reflectivity_table(PhotonReflectTableStruct &table, bool in_degrees);
extern "C" void fortran_find_element_ends(
    void *ele /* 0D_NOT_type in */,
    void *ele1 /* 0D_PTR_type out */,
    void *ele2 /* 0D_PTR_type out */,
    int *ix_multipass /* 0D_NOT_integer in */
);
struct FindElementEnds {
  std::optional<EleStruct> ele1;
  std::optional<EleStruct> ele2;
};
Bmad::FindElementEnds
find_element_ends(EleStruct &ele, std::optional<int> ix_multipass = std::nullopt);
extern "C" void fortran_find_fwhm(
    double &bound /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    double &fwhm /* 0D_NOT_real out */
);
double find_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_find_matching_fieldmap(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    int &fm_type /* 0D_NOT_integer in */,
    void *match_ele /* 0D_PTR_type out */,
    int &ix_field /* 0D_NOT_integer out */,
    bool *ignore_slaves /* 0D_NOT_logical in */
);
struct FindMatchingFieldmap {
  std::optional<EleStruct> match_ele;
  int ix_field;
};
Bmad::FindMatchingFieldmap find_matching_fieldmap(
    std::string file_name,
    EleStruct &ele,
    int fm_type,
    std::optional<bool> ignore_slaves = std::nullopt
);
extern "C" void fortran_find_normalization(
    double &bound /* 0D_NOT_real in */,
    double &p0 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    double &pnrml /* 0D_NOT_real out */
);
double find_normalization(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_floor_angles_to_w_mat(
    double &theta /* 0D_NOT_real in */,
    double &phi /* 0D_NOT_real in */,
    double &psi /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &w_mat_inv /* 2D_NOT_real out */
);
struct FloorAnglesToWMat {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat_inv;
};
Bmad::FloorAnglesToWMat floor_angles_to_w_mat(double theta, double phi, double psi);
extern "C" void fortran_floor_w_mat_to_angles(
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real in */,
    double &theta /* 0D_NOT_real out */,
    double &phi /* 0D_NOT_real out */,
    double &psi /* 0D_NOT_real out */,
    void *floor0 /* 0D_NOT_type in */
);
struct FloorWMatToAngles {
  double theta;
  double phi;
  double psi;
};
Bmad::FloorWMatToAngles floor_w_mat_to_angles(
    FixedArray2D<Real, 3, 3> w_mat,
    optional_ref<FloorPositionStruct> floor0 = std::nullopt
);
extern "C" void fortran_form_complex_taylor(
    void *re_taylor /* 0D_NOT_type in */,
    void *im_taylor /* 0D_NOT_type in */,
    void *complex_taylor /* 0D_NOT_type out */
);
ComplexTaylorStruct form_complex_taylor(TaylorStruct &re_taylor, TaylorStruct &im_taylor);
extern "C" void fortran_form_digested_bmad_file_name(
    const char *lat_file /* 0D_NOT_character in */,
    const char *digested_file /* 0D_NOT_character out */,
    const char *full_lat_file /* 0D_NOT_character out */,
    const char *use_line /* 0D_NOT_character in */
);
struct FormDigestedBmadFileName {
  std::string digested_file;
  std::string full_lat_file;
};
Bmad::FormDigestedBmadFileName form_digested_bmad_file_name(
    std::string lat_file,
    std::optional<std::string> use_line = std::nullopt
);
extern "C" bool fortran_fringe_here(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    bool &is_here /* 0D_NOT_logical out */
);
bool fringe_here(EleStruct &ele, CoordStruct &orbit, int particle_at);
extern "C" bool fortran_g_bend_from_em_field(
    Bmad::array_descriptor_t &b /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &e /* 1D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &g_bend /* 1D_NOT_real out */
);
FixedArray1D<Real, 3>
g_bend_from_em_field(FixedArray1D<Real, 3> b, FixedArray1D<Real, 3> e, CoordStruct &orbit);
extern "C" void fortran_g_bending_strength_from_em_field(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &s_rel /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    bool &local_ref_frame /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &g /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &dg /* 2D_NOT_real out */
);
struct GBendingStrengthFromEmField {
  FixedArray1D<Real, 3> g;
  std::optional<FixedArray2D<Real, 3, 3>> dg;
};
Bmad::GBendingStrengthFromEmField g_bending_strength_from_em_field(
    EleStruct &ele,
    LatParamStruct &param,
    double s_rel,
    CoordStruct &orbit,
    bool local_ref_frame
);
extern "C" void fortran_g_integrals_calc(void *lat /* 0D_NOT_type inout */);
void g_integrals_calc(LatStruct &lat);
extern "C" bool
fortran_gamma_ref(void *ele /* 0D_NOT_type in */, double &gamma /* 0D_NOT_real out */);
double gamma_ref(EleStruct &ele);
extern "C" void fortran_gen_grad1_to_gg_taylor(
    void *ele /* 0D_NOT_type in */,
    void *gen_grad /* 0D_NOT_type in */,
    int &iz /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &gg_taylor /* 1D_NOT_type out */
);
GgTaylorStructArray1D gen_grad1_to_gg_taylor(EleStruct &ele, GenGradMapStruct &gen_grad, int iz);
extern "C" void fortran_gen_grad_at_s_to_gg_taylor(
    void *ele /* 0D_NOT_type in */,
    void *gen_grad /* 0D_NOT_type in */,
    double &s_pos /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &gg_taylor /* 1D_NOT_type out */
);
GgTaylorStructArray1D
gen_grad_at_s_to_gg_taylor(EleStruct &ele, GenGradMapStruct &gen_grad, double s_pos);
extern "C" bool fortran_gen_grad_field(
    Bmad::array_descriptor_t &deriv /* 1D_NOT_real inout */,
    void *gg /* 0D_NOT_type inout */,
    double &rho /* 0D_NOT_real in */,
    double &theta /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &field /* 1D_NOT_real inout */
);
void gen_grad_field(
    FArray1D<Real> &deriv,
    GenGrad1Struct &gg,
    double rho,
    double theta,
    FixedArray1D<Real, 3> field
);

// Skipped unusable routine get_astra_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_get_bl_from_fwhm(
    double &bound /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    double &sigma /* 0D_NOT_real out */
);
double get_bl_from_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_get_called_file(
    const char *delim /* 0D_NOT_character in */,
    const char *call_file /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical in */
);
void get_called_file(std::string delim, std::string call_file, bool err);

// Skipped unusable routine get_cgrn_csr3d:
// - Variable inout sized array: 3D_NOT_complex
extern "C" void fortran_get_emit_from_sigma_mat(
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &normal /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &Nmat /* 2D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct GetEmitFromSigmaMat {
  FixedArray1D<Real, 3> normal;
  bool err_flag;
};
Bmad::GetEmitFromSigmaMat get_emit_from_sigma_mat(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<FixedArray2D<Real, 6, 6>> Nmat = std::nullopt
);

// Skipped unusable routine get_gpt_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)

// Skipped unusable routine get_list_of_names:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_get_next_word(
    const char *word /* 0D_NOT_character in */,
    int &ix_word /* 0D_NOT_integer out */,
    const char *delim_list /* 0D_NOT_character in */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    bool *upper_case_word /* 0D_NOT_logical in */,
    bool *call_check /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct GetNextWord {
  int ix_word;
  std::string delim;
  bool delim_found;
  bool err_flag;
};
Bmad::GetNextWord get_next_word(
    std::string word,
    std::string delim_list,
    std::optional<bool> upper_case_word = std::nullopt,
    std::optional<bool> call_check = std::nullopt
);

// Skipped unusable routine get_opal_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)

// Skipped unusable routine get_overlay_group_names:
// - Untranslated type: parser_ele_struct (0D)
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine get_sequence_args:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_get_slave_list(
    void *lord /* 0D_NOT_type in */,
    void *slaves /* 1D_ALLOC_type out */,
    int &n_slave /* 0D_NOT_integer out */
);
struct GetSlaveList {
  ElePointerStructAlloc1D slaves;
  int n_slave;
};
Bmad::GetSlaveList get_slave_list(EleStruct &lord);

// Skipped unusable routine get_switch:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine getrhotilde:
// - Array bounds handling: "Enum 'ILO' found in bounds 'ilo' but not in provided map."
// - Array bounds handling: "Enum 'ILO' found in bounds 'ilo' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_gg_taylor_equal_gg_taylor(
    void *gg_taylor1 /* 0D_NOT_type inout */,
    void *gg_taylor2 /* 0D_NOT_type in */
);
void gg_taylor_equal_gg_taylor(GgTaylorStruct &gg_taylor1, GgTaylorStruct &gg_taylor2);
extern "C" void fortran_gg_taylors_equal_gg_taylors(
    Bmad::array_descriptor_t &gg_taylor1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &gg_taylor2 /* 1D_NOT_type in */
);
void gg_taylors_equal_gg_taylors(
    GgTaylorStructArray1D gg_taylor1,
    GgTaylorStructArray1D gg_taylor2
);
extern "C" void fortran_gpt_field_grid_scaling(
    void *ele /* 0D_NOT_type inout */,
    int &dimensions /* 0D_NOT_integer in */,
    double &field_scale /* 0D_NOT_real in */,
    double &ref_time /* 0D_NOT_real in */
);
void gpt_field_grid_scaling(EleStruct &ele, int dimensions, double field_scale, double ref_time);
extern "C" bool fortran_gpt_max_field_reference(
    void *pt0 /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type inout */,
    double &field_value /* 0D_NOT_real in */
);
void gpt_max_field_reference(GridFieldPt1Struct &pt0, EleStruct &ele, double field_value);
extern "C" void fortran_gpt_to_particle_bunch(
    const char *gpt_file /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *bunch /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct GptToParticleBunch {
  BunchStruct bunch;
  bool err_flag;
};
Bmad::GptToParticleBunch gpt_to_particle_bunch(std::string gpt_file, EleStruct &ele);
extern "C" bool fortran_gradient_shift_sr_wake(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &grad_shift /* 0D_NOT_real out */
);
double gradient_shift_sr_wake(EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_grid_field_interpolate(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    void *grid /* 0D_NOT_type in */,
    void *g_field /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical in */,
    double &x1 /* 0D_NOT_real in */,
    double *x2 /* 0D_NOT_real in */,
    double *x3 /* 0D_NOT_real in */,
    bool *allow_s_out_of_bounds /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */
);
GridFieldPt1Struct grid_field_interpolate(
    EleStruct &ele,
    CoordStruct &orbit,
    GridFieldStruct &grid,
    bool err_flag,
    double x1,
    std::optional<double> x2 = std::nullopt,
    std::optional<double> x3 = std::nullopt,
    std::optional<bool> allow_s_out_of_bounds = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_hard_multipole_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void hard_multipole_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_has_attribute(
    void *ele /* 0D_NOT_type inout */,
    const char *attrib /* 0D_NOT_character in */,
    bool &has_it /* 0D_NOT_logical in */
);
void has_attribute(EleStruct &ele, std::string attrib, bool has_it);
extern "C" bool
fortran_has_curvature(void *phot_ele /* 0D_NOT_type in */, bool &curved /* 0D_NOT_logical out */);
bool has_curvature(PhotonElementStruct &phot_ele);
extern "C" bool fortran_has_orientation_attributes(
    void *ele /* 0D_NOT_type in */,
    bool &has_attribs /* 0D_NOT_logical out */
);
bool has_orientation_attributes(EleStruct &ele);

// Skipped unusable routine hdf5_read_beam:
// - Untranslated type: pmd_header_struct (0D)

// Skipped unusable routine hdf5_read_grid_field:
// - Untranslated type: pmd_header_struct (0D)
extern "C" void fortran_hdf5_write_beam(
    const char *file_name /* 0D_NOT_character in */,
    Bmad::array_descriptor_t &bunches /* 1D_NOT_type inout */,
    bool &append /* 0D_NOT_logical in */,
    bool &error /* 0D_NOT_logical in */,
    void *lat /* 0D_NOT_type inout */,
    bool *alive_only /* 0D_NOT_logical in */
);
void hdf5_write_beam(
    std::string file_name,
    BunchStructArray1D bunches,
    bool append,
    bool error,
    optional_ref<LatStruct> lat = std::nullopt,
    std::optional<bool> alive_only = std::nullopt
);
extern "C" void fortran_hdf5_write_grid_field(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &g_field /* 1D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void hdf5_write_grid_field(
    std::string file_name,
    EleStruct &ele,
    GridFieldStructArray1D g_field,
    bool err_flag
);
extern "C" void fortran_hwang_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void hwang_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);

// Skipped unusable routine i_csr:
// - Untranslated type: csr_kick1_struct (0D)
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine ibs1:
// - Untranslated type: ibs_sim_param_struct (0D)
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine ibs_blowup1turn:
// - Untranslated type: ibs_sim_param_struct (0D)

// Skipped unusable routine ibs_delta_calc:
// - Untranslated type: ibs_sim_param_struct (0D)

// Skipped unusable routine ibs_equib_der:
// - Untranslated type: ibs_sim_param_struct (0D)

// Skipped unusable routine ibs_equib_rlx:
// - Untranslated type: ibs_sim_param_struct (0D)

// Skipped unusable routine ibs_lifetime:
// - Untranslated type: ibs_sim_param_struct (0D)
// - Untranslated type: ibs_maxratio_struct (0D)
// - Untranslated type: ibs_lifetime_struct (0D)
extern "C" bool fortran_ibs_matrix_c(
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real inout */,
    bool &tail_cut /* 0D_NOT_logical in */,
    double &tau /* 0D_NOT_real in */,
    double &energy /* 0D_NOT_real in */,
    double &n_part /* 0D_NOT_real in */,
    int &species /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &ibs_mat /* 2D_NOT_real inout */
);
void ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species,
    FixedArray2D<Real, 6, 6> ibs_mat
);

// Skipped unusable routine ibs_rates1turn:
// - Untranslated type: ibs_sim_param_struct (0D)
// - Untranslated type: ibs_struct (0D)
extern "C" bool fortran_igfcoulombfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &dx /* 0D_NOT_real in */,
    double &dy /* 0D_NOT_real in */,
    double &dz /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void igfcoulombfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
);
extern "C" bool fortran_igfexfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &dx /* 0D_NOT_real in */,
    double &dy /* 0D_NOT_real in */,
    double &dz /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void igfexfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
);
extern "C" bool fortran_igfeyfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &dx /* 0D_NOT_real in */,
    double &dy /* 0D_NOT_real in */,
    double &dz /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void igfeyfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
);
extern "C" bool fortran_igfezfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &dx /* 0D_NOT_real in */,
    double &dy /* 0D_NOT_real in */,
    double &dz /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void igfezfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res
);

// Skipped unusable routine image_charge_kick_calc:
// - Untranslated type: csr_kick1_struct (0D)
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine imageconvcorr3d:
// - Array bounds handling: "Enum 'RILO' found in bounds 'rilo' but not in provided map."
// - Array bounds handling: "Enum 'G1ILO' found in bounds 'g1ilo' but not in provided map."
// - Array bounds handling: "Enum 'CILO' found in bounds 'cilo' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_init_attribute_name1(
    bool &is_ok /* 0D_NOT_logical inout */,
    int &ix_key /* 0D_NOT_integer in */,
    int &ix_attrib /* 0D_NOT_integer in */,
    const char *name /* 0D_NOT_character in */,
    int *attrib_state /* 0D_NOT_integer in */,
    bool *override /* 0D_NOT_logical in */
);
void init_attribute_name1(
    bool &is_ok,
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state = std::nullopt,
    std::optional<bool> override = std::nullopt
);
extern "C" void fortran_init_attribute_name_array();
void init_attribute_name_array();
extern "C" void fortran_init_beam_distribution(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *beam_init /* 0D_NOT_type in */,
    void *beam /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *modes /* 0D_NOT_type in */,
    void *beam_init_set /* 0D_NOT_type out */,
    bool *print_p0c_shift_warning /* 0D_NOT_logical in */,
    bool *conserve_momentum /* 0D_NOT_logical in */
);
struct InitBeamDistribution {
  BeamStruct beam;
  bool err_flag;
  BeamInitStruct beam_init_set;
};
Bmad::InitBeamDistribution init_beam_distribution(
    EleStruct &ele,
    LatParamStruct &param,
    BeamInitStruct &beam_init,
    optional_ref<NormalModesStruct> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    std::optional<bool> conserve_momentum = std::nullopt
);
extern "C" void fortran_init_bmad();
void init_bmad();
extern "C" void fortran_init_bmad_parser_common(void *lat /* 0D_NOT_type inout */);
void init_bmad_parser_common(optional_ref<LatStruct> lat = std::nullopt);
extern "C" void fortran_init_bunch_distribution(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *beam_init /* 0D_NOT_type in */,
    int &ix_bunch /* 0D_NOT_integer in */,
    void *bunch /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *modes /* 0D_NOT_type in */,
    void *beam_init_used /* 0D_NOT_type out */,
    bool *print_p0c_shift_warning /* 0D_NOT_logical in */,
    bool *conserve_momentum /* 0D_NOT_logical in */
);
struct InitBunchDistribution {
  BunchStruct bunch;
  bool err_flag;
  BeamInitStruct beam_init_used;
};
Bmad::InitBunchDistribution init_bunch_distribution(
    EleStruct &ele,
    LatParamStruct &param,
    BeamInitStruct &beam_init,
    int ix_bunch,
    optional_ref<NormalModesStruct> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    std::optional<bool> conserve_momentum = std::nullopt
);
extern "C" void fortran_init_complex_taylor_series(
    void *complex_taylor /* 0D_NOT_type inout */,
    int &n_term /* 0D_NOT_integer in */,
    bool *save /* 0D_NOT_logical in */
);
void init_complex_taylor_series(
    ComplexTaylorStruct &complex_taylor,
    int n_term,
    std::optional<bool> save = std::nullopt
);
extern "C" void fortran_init_coord1(
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &vec /* 1D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    int *element_end /* 0D_NOT_integer in */,
    int *particle /* 0D_NOT_integer in */,
    int *direction /* 0D_NOT_integer in */,
    double *E_photon /* 0D_NOT_real in */,
    double *t_offset /* 0D_NOT_real in */,
    bool *shift_vec6 /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin /* 1D_NOT_real in */,
    double *s_pos /* 0D_NOT_real in */,
    bool *random_on /* 0D_NOT_logical in */
);
void init_coord(
    CoordStruct &orb,
    FixedArray1D<Real, 6> vec,
    optional_ref<EleStruct> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> random_on = std::nullopt
);
extern "C" void fortran_init_coord2(
    void *orb_out /* 0D_NOT_type out */,
    void *orb_in /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    int *element_end /* 0D_NOT_integer in */,
    int *particle /* 0D_NOT_integer in */,
    int *direction /* 0D_NOT_integer in */,
    double *E_photon /* 0D_NOT_real in */,
    double *t_offset /* 0D_NOT_real in */,
    bool *shift_vec6 /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin /* 1D_NOT_real in */,
    double *s_pos /* 0D_NOT_real in */,
    bool *random_on /* 0D_NOT_logical in */
);
CoordStruct init_coord(
    CoordStruct &orb_in,
    optional_ref<EleStruct> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> random_on = std::nullopt
);
extern "C" void fortran_init_coord3(
    void *orb /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    int *element_end /* 0D_NOT_integer in */,
    int *particle /* 0D_NOT_integer in */,
    int *direction /* 0D_NOT_integer in */,
    double *E_photon /* 0D_NOT_real in */,
    double *t_offset /* 0D_NOT_real in */,
    bool *shift_vec6 /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin /* 1D_NOT_real in */
);
void init_coord(
    CoordStruct &orb,
    optional_ref<EleStruct> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt
);
extern "C" void fortran_init_custom(void *lat /* 0D_NOT_type inout */);
void init_custom(LatStruct &lat);

// Skipped unusable routine init_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_init_ele(
    void *ele /* 0D_NOT_type out */,
    int *key /* 0D_NOT_integer in */,
    int *sub_key /* 0D_NOT_integer in */,
    int *ix_ele /* 0D_NOT_integer in */,
    void *branch /* 0D_NOT_type in */
);
EleStruct init_ele(
    std::optional<int> key = std::nullopt,
    std::optional<int> sub_key = std::nullopt,
    std::optional<int> ix_ele = std::nullopt,
    optional_ref<BranchStruct> branch = std::nullopt
);

// Skipped unusable routine init_fringe_info:
// - Untranslated type: fringe_field_info_struct (0D)
extern "C" void fortran_init_gg_taylor_series(
    void *gg_taylor /* 0D_NOT_type inout */,
    int &n_term /* 0D_NOT_integer in */,
    bool *save_old /* 0D_NOT_logical in */
);
void init_gg_taylor_series(
    GgTaylorStruct &gg_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt
);
extern "C" void fortran_init_lat(
    void *lat /* 0D_NOT_type out */,
    int *n /* 0D_NOT_integer in */,
    bool *init_beginning_ele /* 0D_NOT_logical in */
);
LatStruct init_lat(
    std::optional<int> n = std::nullopt,
    std::optional<bool> init_beginning_ele = std::nullopt
);
extern "C" void fortran_init_multipole_cache(void *ele /* 0D_NOT_type inout */);
void init_multipole_cache(EleStruct &ele);
extern "C" void fortran_init_photon_from_a_photon_init_ele(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type out */,
    bool *random_on /* 0D_NOT_logical in */
);
CoordStruct init_photon_from_a_photon_init_ele(
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> random_on = std::nullopt
);
extern "C" bool fortran_init_photon_integ_prob(
    double &gamma /* 0D_NOT_real in */,
    double &g /* 0D_NOT_real in */,
    double &E_min /* 0D_NOT_real in */,
    double &E_max /* 0D_NOT_real in */,
    double *vert_angle_min /* 0D_NOT_real in */,
    double *vert_angle_max /* 0D_NOT_real in */,
    bool *vert_angle_symmetric /* 0D_NOT_logical in */,
    double *energy_integ_prob /* 0D_NOT_real in */,
    double &E_photon /* 0D_NOT_real out */,
    double &integ_prob /* 0D_NOT_real out */
);
struct InitPhotonIntegProb {
  double E_photon;
  double integ_prob;
};
Bmad::InitPhotonIntegProb init_photon_integ_prob(
    double gamma,
    double g,
    double E_min,
    double E_max,
    std::optional<double> vert_angle_min = std::nullopt,
    std::optional<double> vert_angle_max = std::nullopt,
    std::optional<bool> vert_angle_symmetric = std::nullopt,
    std::optional<double> energy_integ_prob = std::nullopt
);
extern "C" void fortran_init_spin_distribution(
    void *beam_init /* 0D_NOT_type in */,
    void *bunch /* 0D_NOT_type out */,
    void *ele /* 0D_NOT_type inout */
);
BunchStruct init_spin_distribution(BeamInitStruct &beam_init, EleStruct &ele);
extern "C" void fortran_init_surface_segment(
    void *phot /* 0D_NOT_type in */,
    int &ix /* 0D_NOT_integer in */,
    int &iy /* 0D_NOT_integer in */
);
void init_surface_segment(PhotonElementStruct &phot, int ix, int iy);
extern "C" void fortran_init_taylor_series(
    void *bmad_taylor /* 0D_NOT_type inout */,
    int &n_term /* 0D_NOT_integer in */,
    bool *save_old /* 0D_NOT_logical in */
);
void init_taylor_series(
    TaylorStruct &bmad_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt
);
extern "C" void fortran_init_wake(
    void *wake /* 0D_PTR_type out */,
    int &n_sr_long /* 0D_NOT_integer in */,
    int &n_sr_trans /* 0D_NOT_integer in */,
    int &n_sr_z /* 0D_NOT_integer in */,
    int &n_lr_mode /* 0D_NOT_integer in */,
    bool *always_allocate /* 0D_NOT_logical in */
);
std::optional<WakeStruct> init_wake(
    int n_sr_long,
    int n_sr_trans,
    int n_sr_z,
    int n_lr_mode,
    std::optional<bool> always_allocate = std::nullopt
);
extern "C" void fortran_insert_element(
    void *lat /* 0D_NOT_type inout */,
    void *insert_ele /* 0D_NOT_type in */,
    int &ix_ele /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    void *orbit /* 1D_ALLOC_type inout */
);
void insert_element(
    LatStruct &lat,
    EleStruct &insert_ele,
    int ix_ele,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<CoordStructAlloc1D> orbit = std::nullopt
);

// Skipped unusable routine integrand:
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_integrand_base(
    double &t /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real inout */,
    double &func_retval__ /* 0D_NOT_real in */
);
void integrand_base(double t, FArray1D<Real> &args, double func_retval__);

// Skipped unusable routine integrand_base_cov:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine integrand_zap:
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_integrate_psi(
    double &bound /* 0D_NOT_real in */,
    double &p0 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    double &result /* 0D_NOT_real out */
);
double integrate_psi(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_integrated_mats(
    Bmad::array_descriptor_t &eles /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &coos /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &Lambda /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &Theta /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &Iota /* 2D_NOT_complex inout */,
    void *mode /* 0D_NOT_type inout */
);
void integrated_mats(
    EleStructArray1D eles,
    CoordStructArray1D coos,
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    NormalModesStruct &mode
);
extern "C" void fortran_integration_timer_ele(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type inout */,
    void *start /* 0D_NOT_type in */,
    void *orb_max /* 0D_NOT_type in */,
    double &tol /* 0D_NOT_real in */
);
void integration_timer(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &start,
    CoordStruct &orb_max,
    double tol
);
extern "C" void fortran_integration_timer_fibre(
    void *a_fibre /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &orbit /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &orbit_max /* 1D_NOT_real in */,
    double &tol_dp /* 0D_NOT_real in */
);
void integration_timer(
    Fibre &a_fibre,
    FixedArray1D<Real, 6> orbit,
    FixedArray1D<Real, 6> orbit_max,
    double tol_dp
);

// Skipped unusable routine interpolate_field:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine interpolate_field_batch:
// - Array bounds handling: "Enum 'N_PTS' found in bounds 'n_pts' but not in provided map."
// - Array bounds handling: "Enum 'N_PTS' found in bounds 'n_pts' but not in provided map."
// - Array bounds handling: "Enum 'N_PTS' found in bounds 'n_pts' but not in provided map."
// - Untranslated type: mesh3d_struct (0D)
// - Array bounds handling: "Enum 'N_PTS' found in bounds 'n_pts' but not in provided map."
// - Array bounds handling: "Enum 'N_PTS' found in bounds 'n_pts' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_ion_kick(
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &r_beam /* 1D_NOT_real in */,
    double &n_beam_part /* 0D_NOT_real in */,
    void *a_twiss /* 0D_NOT_type in */,
    void *b_twiss /* 0D_NOT_type in */,
    double &sig_ee /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &kick /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> ion_kick(
    CoordStruct &orbit,
    FixedArray1D<Real, 2> r_beam,
    double n_beam_part,
    TwissStruct &a_twiss,
    TwissStruct &b_twiss,
    double sig_ee
);
extern "C" bool fortran_is_attribute(
    int &ix_attrib /* 0D_NOT_integer in */,
    int &which /* 0D_NOT_integer in */,
    bool &is_attrib /* 0D_NOT_logical out */
);
bool is_attribute(int ix_attrib, int which);

// Skipped unusable routine jac:
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_key_name_to_key_index(
    const char *key_str /* 0D_NOT_character in */,
    bool *abbrev_allowed /* 0D_NOT_logical in */,
    int &key_index /* 0D_NOT_integer out */
);
int key_name_to_key_index(std::string key_str, std::optional<bool> abbrev_allowed = std::nullopt);
extern "C" void fortran_kick_vector_calc(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &s_body /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &dr_ds /* 1D_NOT_real out */,
    bool &err /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
struct KickVectorCalc {
  FixedArray1D<Real, 11> dr_ds;
  bool err;
};
Bmad::KickVectorCalc kick_vector_calc(
    EleStruct &ele,
    LatParamStruct &param,
    double s_body,
    CoordStruct &orbit,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void
fortran_kill_complex_taylor(Bmad::array_descriptor_t &complex_taylor /* 1D_NOT_type inout */);
void kill_complex_taylor(ComplexTaylorStructArray1D complex_taylor);
extern "C" void fortran_kill_ptc_layouts(void *lat /* 0D_NOT_type in */);
void kill_ptc_layouts(LatStruct &lat);
extern "C" void fortran_kill_taylor(Bmad::array_descriptor_t &bmad_taylor /* 1D_NOT_type inout */);
void kill_taylor(TaylorStructArray1D bmad_taylor);
extern "C" bool fortran_kind_name(
    int *this_kind /* 0D_PTR_integer in */,
    const char *kind_str /* 0D_NOT_character out */
);
std::string kind_name(int this_kind);
extern "C" bool fortran_knot_interpolate(
    Bmad::array_descriptor_t &x_knot /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y_knot /* 1D_NOT_real in */,
    double &x_pt /* 0D_NOT_real in */,
    int &interpolation /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &y_pt /* 0D_NOT_real out */
);
struct KnotInterpolate {
  bool err_flag;
  double y_pt;
};
Bmad::KnotInterpolate
knot_interpolate(FArray1D<Real> &x_knot, FArray1D<Real> &y_knot, double x_pt, int interpolation);
extern "C" bool fortran_knots_to_string(
    Bmad::array_descriptor_t &x_knot /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &y_knot /* 1D_NOT_real inout */,
    const char *str /* 0D_ALLOC_character in */
);
void knots_to_string(FArray1D<Real> &x_knot, FArray1D<Real> &y_knot, std::string str);

// Skipped unusable routine kubo_integrand:
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_lafun(
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    double &z /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void lafun(double x, double y, double z, double res);
extern "C" void fortran_lat_compute_ref_energy_and_time(
    void *lat /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool lat_compute_ref_energy_and_time(LatStruct &lat);
extern "C" void fortran_lat_ele_locator(
    const char *loc_str /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    void *eles /* 1D_ALLOC_type inout */,
    int &n_loc /* 0D_NOT_integer inout */,
    bool &err /* 0D_NOT_logical out */,
    bool *above_ubound_is_err /* 0D_NOT_logical in */,
    int *ix_dflt_branch /* 0D_NOT_integer in */,
    bool *order_by_index /* 0D_NOT_logical in */,
    bool *append_eles /* 0D_NOT_logical in */
);
bool lat_ele_locator(
    std::string loc_str,
    LatStruct &lat,
    ElePointerStructAlloc1D eles,
    int &n_loc,
    std::optional<bool> above_ubound_is_err = std::nullopt,
    std::optional<int> ix_dflt_branch = std::nullopt,
    std::optional<bool> order_by_index = std::nullopt,
    std::optional<bool> append_eles = std::nullopt
);
extern "C" void
fortran_lat_equal_lat(void *lat_out /* 0D_NOT_type inout */, void *lat_in /* 0D_NOT_type in */);
void lat_equal_lat(LatStruct &lat_out, LatStruct &lat_in);
extern "C" void fortran_lat_geometry(void *lat /* 0D_NOT_type inout */);
void lat_geometry(LatStruct &lat);
extern "C" void fortran_lat_make_mat6(
    void *lat /* 0D_NOT_type inout */,
    int *ix_ele /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &ref_orb /* 1D_NOT_type in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool lat_make_mat6(
    LatStruct &lat,
    std::optional<int> ix_ele = std::nullopt,
    std::optional<CoordStructArray1D> ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
);

// Skipped unusable routine lat_make_mat6_hook_def:
// - Routine in configuration skip list
extern "C" void
fortran_lat_sanity_check(void *lat /* 0D_NOT_type in */, bool &err_flag /* 0D_NOT_logical out */);
bool lat_sanity_check(LatStruct &lat);
extern "C" void fortran_lat_to_ptc_layout(void *lat /* 0D_NOT_type in */);
void lat_to_ptc_layout(LatStruct &lat);
extern "C" void fortran_lat_vec_equal_lat_vec(
    Bmad::array_descriptor_t &lat1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &lat2 /* 1D_NOT_type in */
);
void lat_vec_equal_lat_vec(LatStructArray1D lat1, LatStructArray1D lat2);
extern "C" void fortran_lattice_bookkeeper(
    void *lat /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool lattice_bookkeeper(LatStruct &lat);
extern "C" void fortran_lcavity_rf_step_setup(void *ele /* 0D_NOT_type inout */);
void lcavity_rf_step_setup(EleStruct &ele);
extern "C" void fortran_linear_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void linear_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_linear_coef(
    Bmad::array_descriptor_t &stack /* 1D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &coef /* 0D_NOT_real out */
);
struct LinearCoef {
  bool err_flag;
  double coef;
};
Bmad::LinearCoef linear_coef(ExpressionAtomStructArray1D stack);
extern "C" void fortran_linear_to_spin_taylor(
    Bmad::array_descriptor_t &q_map /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &spin_taylor /* 1D_NOT_type out */
);
TaylorStructArray1D linear_to_spin_taylor(FixedArray2D<Real, 4, 7> q_map);
extern "C" void fortran_load_parse_line(
    const char *action /* 0D_NOT_character in */,
    int &ix_start /* 0D_NOT_integer in */,
    bool &end_of_file /* 0D_NOT_logical out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct LoadParseLine {
  bool end_of_file;
  bool err_flag;
};
Bmad::LoadParseLine load_parse_line(std::string action, int ix_start);
extern "C" bool fortran_lord_edge_aligned(
    void *slave /* 0D_NOT_type in */,
    int &slave_edge /* 0D_NOT_integer in */,
    void *lord /* 0D_NOT_type in */,
    bool &is_aligned /* 0D_NOT_logical out */
);
bool lord_edge_aligned(EleStruct &slave, int slave_edge, EleStruct &lord);
extern "C" bool fortran_low_energy_z_correction(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    double &ds /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */,
    double &dz /* 0D_NOT_real out */
);
double low_energy_z_correction(
    CoordStruct &orbit,
    EleStruct &ele,
    double ds,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);

// Skipped unusable routine lsc_kick_params_calc:
// - Untranslated type: csr_struct (0D)
extern "C" void fortran_mad_add_offsets_and_multipoles(
    void *ele /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_add_offsets_and_multipoles(EleStruct &ele);
extern "C" void fortran_mad_concat_map2(
    void *map1 /* 0D_NOT_type in */,
    void *map2 /* 0D_NOT_type in */,
    void *map3 /* 0D_NOT_type out */
);
MadMapStruct mad_concat_map2(MadMapStruct &map1, MadMapStruct &map2);
extern "C" void fortran_mad_drift(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_drift(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_elsep(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_elsep(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_map_to_taylor(
    void *map /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &taylor /* 1D_NOT_type inout */
);
void mad_map_to_taylor(MadMapStruct &map, MadEnergyStruct &energy, TaylorStructArray1D taylor);
extern "C" void fortran_mad_quadrupole(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_quadrupole(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_rfcavity(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_rfcavity(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_sbend(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_sbend(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_sbend_body(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_sbend_body(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_sbend_fringe(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    bool &into /* 0D_NOT_logical in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_sbend_fringe(EleStruct &ele, MadEnergyStruct &energy, bool into);
extern "C" void fortran_mad_sextupole(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_sextupole(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_solenoid(
    void *ele /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct mad_solenoid(EleStruct &ele, MadEnergyStruct &energy);
extern "C" void fortran_mad_tmfoc(
    double &el /* 0D_NOT_real in */,
    double &sk1 /* 0D_NOT_real in */,
    double &c /* 0D_NOT_real out */,
    double &s /* 0D_NOT_real out */,
    double &d /* 0D_NOT_real out */,
    double &f /* 0D_NOT_real out */
);
struct MadTmfoc {
  double c;
  double s;
  double d;
  double f;
};
Bmad::MadTmfoc mad_tmfoc(double el, double sk1);
extern "C" void fortran_mad_tmsymm(Bmad::array_descriptor_t &te /* 3D_NOT_real inout */);
void mad_tmsymm(FixedArray3D<Real, 6, 6, 6> te);
extern "C" void
fortran_mad_tmtilt(void *map /* 0D_NOT_type inout */, double &tilt /* 0D_NOT_real in */);
void mad_tmtilt(MadMapStruct &map, double tilt);
extern "C" void fortran_mad_track1(
    void *c0 /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type in */,
    void *c1 /* 0D_NOT_type out */
);
CoordStruct mad_track1(CoordStruct &c0, MadMapStruct &map);
extern "C" void fortran_make_g2_mats(
    void *twiss /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &g2_mat /* 2D_NOT_real inout */,
    Bmad::array_descriptor_t &g2_inv_mat /* 2D_NOT_real inout */
);
void make_g2_mats(
    TwissStruct &twiss,
    FixedArray2D<Real, 2, 2> g2_mat,
    FixedArray2D<Real, 2, 2> g2_inv_mat
);
extern "C" void fortran_make_g_mats(
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &g_mat /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &g_inv_mat /* 2D_NOT_real out */
);
struct MakeGMats {
  FixedArray2D<Real, 4, 4> g_mat;
  FixedArray2D<Real, 4, 4> g_inv_mat;
};
Bmad::MakeGMats make_g_mats(EleStruct &ele);
extern "C" void fortran_make_hvbp(
    Bmad::array_descriptor_t &N /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &B /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &V /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &H /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &Vbar /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &Hbar /* 2D_NOT_real out */
);
struct MakeHvbp {
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> V;
  FixedArray2D<Real, 6, 6> H;
  std::optional<FixedArray2D<Real, 6, 6>> Vbar;
  std::optional<FixedArray2D<Real, 6, 6>> Hbar;
};
Bmad::MakeHvbp make_hvbp(FixedArray2D<Real, 6, 6> N);
extern "C" void fortran_make_hybrid_lat(
    void *lat_in /* 0D_NOT_type in */,
    void *lat_out /* 0D_NOT_type out */,
    bool *use_taylor /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &orb0_arr /* 1D_NOT_type in */
);
LatStruct make_hybrid_lat(
    LatStruct &lat_in,
    std::optional<bool> use_taylor = std::nullopt,
    std::optional<CoordArrayStructArray1D> orb0_arr = std::nullopt
);
extern "C" void fortran_make_mad_map(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *energy /* 0D_NOT_type out */,
    void *map /* 0D_NOT_type out */
);
struct MakeMadMap {
  MadEnergyStruct energy;
  MadMapStruct map;
};
Bmad::MakeMadMap make_mad_map(EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_make_mat6(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct MakeMat6 {
  CoordStruct end_orb;
  bool err_flag;
};
Bmad::MakeMat6 make_mat6(
    EleStruct &ele,
    LatParamStruct &param,
    optional_ref<CoordStruct> start_orb = std::nullopt
);
extern "C" void fortran_make_mat6_bmad(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    bool &err /* 0D_NOT_logical out */
);
struct MakeMat6Bmad {
  CoordStruct end_orb;
  bool err;
};
Bmad::MakeMat6Bmad make_mat6_bmad(EleStruct &ele, LatParamStruct &param, CoordStruct &start_orb);
extern "C" void fortran_make_mat6_bmad_photon(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    bool &err /* 0D_NOT_logical out */
);
struct MakeMat6BmadPhoton {
  CoordStruct end_orb;
  bool err;
};
Bmad::MakeMat6BmadPhoton
make_mat6_bmad_photon(EleStruct &ele, LatParamStruct &param, CoordStruct &start_orb);

// Skipped unusable routine make_mat6_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_make_mat6_high_energy_space_charge(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */
);
void make_mat6_high_energy_space_charge(EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_make_mat6_mad(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *c0 /* 0D_NOT_type in */,
    void *c1 /* 0D_NOT_type out */
);
CoordStruct make_mat6_mad(EleStruct &ele, LatParamStruct &param, CoordStruct &c0);
extern "C" void fortran_make_mat6_symp_lie_ptc(
    void *ele /* 0D_NOT_type inout */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */
);
CoordStruct make_mat6_symp_lie_ptc(EleStruct &ele, CoordStruct &start_orb);
extern "C" void fortran_make_mat6_taylor(
    void *ele /* 0D_NOT_type inout */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    bool *err_flag /* 0D_NOT_logical in */
);
CoordStruct make_mat6_taylor(
    EleStruct &ele,
    CoordStruct &start_orb,
    std::optional<bool> err_flag = std::nullopt
);
extern "C" void fortran_make_mat6_tracking(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *spin_only /* 0D_NOT_logical in */
);
struct MakeMat6Tracking {
  CoordStruct end_orb;
  bool err_flag;
};
Bmad::MakeMat6Tracking make_mat6_tracking(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &start_orb,
    std::optional<bool> spin_only = std::nullopt
);
extern "C" void fortran_make_n(
    Bmad::array_descriptor_t &t6 /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &N /* 2D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &tunes_out /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &U /* 2D_NOT_real out */
);
struct MakeN {
  FixedArray2D<Real, 6, 6> N;
  bool err_flag;
  FixedArray1D<Real, 3> tunes_out;
  std::optional<FixedArray2D<Real, 6, 6>> U;
};
Bmad::MakeN
make_n(FixedArray2D<Real, 6, 6> t6, std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt);
extern "C" void fortran_make_pbrh(
    Bmad::array_descriptor_t &M /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &P /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &Bp /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &R /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &H /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real in */
);
struct MakePbrh {
  FixedArray2D<Complex, 6, 6> P;
  FixedArray2D<Complex, 6, 6> Bp;
  FixedArray2D<Complex, 6, 6> R;
  FixedArray2D<Complex, 6, 6> H;
};
Bmad::MakePbrh make_pbrh(FixedArray2D<Real, 6, 6> M, FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_make_smat_from_abc(
    Bmad::array_descriptor_t &t6 /* 2D_NOT_real in */,
    void *mode /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &Nout /* 2D_NOT_real out */
);
struct MakeSmatFromAbc {
  FixedArray2D<Real, 6, 6> sigma_mat;
  bool err_flag;
  std::optional<FixedArray2D<Real, 6, 6>> Nout;
};
Bmad::MakeSmatFromAbc make_smat_from_abc(FixedArray2D<Real, 6, 6> t6, NormalModesStruct &mode);

// Skipped unusable routine make_sr_mats:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine make_srdt_cache:
// - Untranslated type: sliced_eles_struct (1D)
// - Variable inout sized array: 3D_ALLOC_complex
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_make_unit_mad_map(void *map /* 0D_NOT_type inout */);
void make_unit_mad_map(MadMapStruct &map);
extern "C" void fortran_make_v(
    Bmad::array_descriptor_t &M /* 2D_NOT_real inout */,
    Bmad::array_descriptor_t &V /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real inout */
);
void make_v(
    FixedArray2D<Real, 6, 6> M,
    FixedArray2D<Complex, 6, 6> V,
    FixedArray1D<Real, 3> abz_tunes
);
extern "C" void fortran_make_v_mats(
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &v_mat /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &v_inv_mat /* 2D_NOT_real out */
);
struct MakeVMats {
  std::optional<FixedArray2D<Real, 4, 4>> v_mat;
  std::optional<FixedArray2D<Real, 4, 4>> v_inv_mat;
};
Bmad::MakeVMats make_v_mats(EleStruct &ele);

// Skipped unusable routine make_ykick_mat:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_makeup_control_slave(
    void *lat /* 0D_NOT_type inout */,
    void *slave /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void makeup_control_slave(LatStruct &lat, EleStruct &slave, bool err_flag);
extern "C" void fortran_makeup_group_lord(
    void *lat /* 0D_NOT_type inout */,
    void *lord /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void makeup_group_lord(LatStruct &lat, EleStruct &lord, bool err_flag);
extern "C" void fortran_makeup_multipass_slave(
    void *lat /* 0D_NOT_type inout */,
    void *slave /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void makeup_multipass_slave(LatStruct &lat, EleStruct &slave, bool err_flag);
extern "C" void fortran_makeup_super_slave(
    void *lat /* 0D_NOT_type inout */,
    void *slave /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void makeup_super_slave(LatStruct &lat, EleStruct &slave, bool err_flag);
extern "C" void fortran_makeup_super_slave1(
    void *slave /* 0D_NOT_type inout */,
    void *lord /* 0D_NOT_type in */,
    double &offset /* 0D_NOT_real in */,
    void *param /* 0D_NOT_type in */,
    bool &include_upstream_end /* 0D_NOT_logical in */,
    bool &include_downstream_end /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool makeup_super_slave1(
    EleStruct &slave,
    EleStruct &lord,
    double offset,
    LatParamStruct &param,
    bool include_upstream_end,
    bool include_downstream_end
);
extern "C" bool
fortran_map1_inverse(void *map1 /* 0D_NOT_type in */, void *inv_map1 /* 0D_NOT_type out */);
SpinOrbitMap1Struct map1_inverse(SpinOrbitMap1Struct &map1);
extern "C" void fortran_map1_make_unit(void *map1 /* 0D_NOT_type out */);
SpinOrbitMap1Struct map1_make_unit();
extern "C" bool fortran_map1_times_map1(
    void *map2 /* 0D_NOT_type in */,
    void *map1 /* 0D_NOT_type in */,
    void *map_out /* 0D_NOT_type inout */
);
void map1_times_map1(
    SpinOrbitMap1Struct &map2,
    SpinOrbitMap1Struct &map1,
    SpinOrbitMap1Struct &map_out
);

// Skipped unusable routine map_coef:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_map_to_angle_coords(
    Bmad::array_descriptor_t &t_canon /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &t_angle /* 1D_NOT_type out */
);
TaylorStructArray1D map_to_angle_coords(TaylorStructArray1D t_canon);
extern "C" void fortran_mark_patch_regions(void *branch /* 0D_NOT_type inout */);
void mark_patch_regions(BranchStruct &branch);
extern "C" bool fortran_master_parameter_value(
    int &master_parameter /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    double &value /* 0D_NOT_real out */
);
double master_parameter_value(int master_parameter, EleStruct &ele);
extern "C" void fortran_mat4_multipole(
    double &knl /* 0D_NOT_real in */,
    double &tilt /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &kick_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 4, 4> mat4_multipole(double knl, double tilt, int n, CoordStruct &orbit);
extern "C" void
fortran_mat6_add_offsets(void *ele /* 0D_NOT_type inout */, void *param /* 0D_NOT_type in */);
void mat6_add_offsets(EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_mat6_add_pitch(
    double &x_pitch_tot /* 0D_NOT_real in */,
    double &y_pitch_tot /* 0D_NOT_real in */,
    int &orientation /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */
);
void mat6_add_pitch(
    double x_pitch_tot,
    double y_pitch_tot,
    int orientation,
    FixedArray2D<Real, 6, 6> mat6
);

// Skipped unusable routine mat6_from_s_to_s:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_mat6_to_complex_taylor(
    Bmad::array_descriptor_t &vec0 /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_complex in */,
    Bmad::array_descriptor_t &complex_taylor /* 1D_NOT_type out */
);
ComplexTaylorStructArray1D
mat6_to_complex_taylor(FixedArray1D<Complex, 6> vec0, FixedArray2D<Complex, 6, 6> mat6);
extern "C" void fortran_mat_symp_decouple(
    Bmad::array_descriptor_t &t0 /* 2D_NOT_real in */,
    int &stat /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t &U /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &V /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &Ubar /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &Vbar /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &G /* 2D_NOT_real out */,
    void *twiss1 /* 0D_NOT_type out */,
    void *twiss2 /* 0D_NOT_type out */,
    double &gamma /* 0D_NOT_real out */,
    bool &type_out /* 0D_NOT_logical in */
);
struct MatSympDecouple {
  int stat;
  FixedArray2D<Real, 4, 4> U;
  FixedArray2D<Real, 4, 4> V;
  FixedArray2D<Real, 4, 4> Ubar;
  FixedArray2D<Real, 4, 4> Vbar;
  FixedArray2D<Real, 4, 4> G;
  TwissStruct twiss1;
  TwissStruct twiss2;
  double gamma;
};
Bmad::MatSympDecouple mat_symp_decouple(FixedArray2D<Real, 4, 4> t0, bool type_out);
extern "C" void fortran_match_ele_to_mat6(
    void *ele /* 0D_NOT_type in */,
    void *start_orb /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &vec0 /* 1D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *include_delta_time /* 0D_NOT_logical in */,
    bool *set_trombone /* 0D_NOT_logical in */
);
struct MatchEleToMat6 {
  FixedArray2D<Real, 6, 6> mat6;
  FixedArray1D<Real, 6> vec0;
  bool err_flag;
};
Bmad::MatchEleToMat6 match_ele_to_mat6(
    EleStruct &ele,
    CoordStruct &start_orb,
    std::optional<bool> include_delta_time = std::nullopt,
    std::optional<bool> set_trombone = std::nullopt
);

// Skipped unusable routine mccfft1d:
// - Routine module (fft_interface_mod) in configuration skip list
extern "C" bool fortran_mexp(
    double &x /* 0D_NOT_real in */,
    int &m /* 0D_NOT_integer in */,
    double &this_exp /* 0D_NOT_real out */
);
double mexp(double x, int m);
extern "C" void fortran_mfft1(
    Bmad::array_descriptor_t &a /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &b /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &n /* 1D_NOT_integer in */,
    int &ndim /* 0D_NOT_integer in */,
    int &isn /* 0D_NOT_integer in */,
    int &ierr /* 0D_NOT_integer out */
);
int mfft1(FArray1D<Real> &a, FArray1D<Real> &b, FArray1D<Int> &n, int ndim, int isn);
extern "C" void fortran_misalign_ptc_fibre(
    void *ele /* 0D_NOT_type in */,
    bool &use_offsets /* 0D_NOT_logical in */,
    void *ptc_fibre /* 0D_PTR_type out */,
    bool &for_layout /* 0D_NOT_logical in */
);
std::optional<Fibre> misalign_ptc_fibre(EleStruct &ele, bool use_offsets, bool for_layout);
extern "C" bool fortran_momentum_compaction(
    void *branch /* 0D_NOT_type in */,
    double &mom_comp /* 0D_NOT_real out */
);
double momentum_compaction(BranchStruct &branch);

// Skipped unusable routine mpxx1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine mpxx_integrand:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine mpzt1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine multi_coulomb_log:
// - Untranslated type: ibs_sim_param_struct (0D)
extern "C" void fortran_multi_turn_tracking_analysis(
    Bmad::array_descriptor_t &track /* 1D_NOT_type in */,
    int &i_dim /* 0D_NOT_integer in */,
    void *track0 /* 0D_NOT_type out */,
    void *ele /* 0D_NOT_type out */,
    bool &stable /* 0D_NOT_logical out */,
    double &growth_rate /* 0D_NOT_real out */,
    double &chi /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct MultiTurnTrackingAnalysis {
  CoordStruct track0;
  EleStruct ele;
  bool stable;
  double growth_rate;
  double chi;
  bool err_flag;
};
Bmad::MultiTurnTrackingAnalysis multi_turn_tracking_analysis(CoordStructArray1D track, int i_dim);

// Skipped unusable routine multi_turn_tracking_to_mat:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_multilayer_type_to_multilayer_params(
    void *ele /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool multilayer_type_to_multilayer_params(EleStruct &ele);

// Skipped unusable routine multipass_all_info:
// - Untranslated type: multipass_all_info_struct (0D)
extern "C" void fortran_multipass_chain(
    void *ele /* 0D_NOT_type in */,
    int &ix_pass /* 0D_NOT_integer out */,
    int &n_links /* 0D_NOT_integer out */,
    void *chain_ele /* 1D_ALLOC_type out */,
    bool *use_super_lord /* 0D_NOT_logical in */
);
struct MultipassChain {
  int ix_pass;
  int n_links;
  ElePointerStructAlloc1D chain_ele;
};
Bmad::MultipassChain
multipass_chain(EleStruct &ele, std::optional<bool> use_super_lord = std::nullopt);

// Skipped unusable routine multipass_region_info:
// - Untranslated type: multipass_region_lat_struct (0D)
// - Untranslated type: multipass_all_info_struct (0D)
extern "C" void fortran_multipole1_ab_to_kt(
    double &an /* 0D_NOT_real in */,
    double &bn /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    double &knl /* 0D_NOT_real out */,
    double &tn /* 0D_NOT_real out */
);
struct Multipole1AbToKt {
  double knl;
  double tn;
};
Bmad::Multipole1AbToKt multipole1_ab_to_kt(double an, double bn, int n);
extern "C" void fortran_multipole1_kt_to_ab(
    double &knl /* 0D_NOT_real in */,
    double &knsl /* 0D_NOT_real in */,
    double &tn /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    double &an /* 0D_NOT_real out */,
    double &bn /* 0D_NOT_real out */
);
struct Multipole1KtToAb {
  double an;
  double bn;
};
Bmad::Multipole1KtToAb multipole1_kt_to_ab(double knl, double knsl, double tn, int n);
extern "C" void fortran_multipole_ab_to_kt(
    Bmad::array_descriptor_t &an /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &bn /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &knl /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &tn /* 1D_NOT_real inout */
);
void multipole_ab_to_kt(
    FArray1D<Real> &an,
    FArray1D<Real> &bn,
    FArray1D<Real> &knl,
    FArray1D<Real> &tn
);
extern "C" void fortran_multipole_ele_to_ab(
    void *ele /* 0D_NOT_type in */,
    bool &use_ele_tilt /* 0D_NOT_logical in */,
    int &ix_pole_max /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t &a /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &b /* 1D_NOT_real out */,
    int *pole_type /* 0D_NOT_integer in */,
    int *include_kicks /* 0D_NOT_integer in */,
    double &b1 /* 0D_NOT_real out */,
    bool *original /* 0D_NOT_logical in */
);
struct MultipoleEleToAb {
  int ix_pole_max;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> a;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> b;
  double b1;
};
Bmad::MultipoleEleToAb multipole_ele_to_ab(
    EleStruct &ele,
    bool use_ele_tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt,
    std::optional<bool> original = std::nullopt
);
extern "C" void fortran_multipole_ele_to_kt(
    void *ele /* 0D_NOT_type in */,
    bool &use_ele_tilt /* 0D_NOT_logical in */,
    int &ix_pole_max /* 0D_NOT_integer out */,
    Bmad::array_descriptor_t &knl /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &tilt /* 1D_NOT_real inout */,
    int *pole_type /* 0D_NOT_integer in */,
    int *include_kicks /* 0D_NOT_integer in */
);
int multipole_ele_to_kt(
    EleStruct &ele,
    bool use_ele_tilt,
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt
);
extern "C" void fortran_multipole_init(
    void *ele /* 0D_NOT_type out */,
    int &who /* 0D_NOT_integer in */,
    bool *zero /* 0D_NOT_logical in */
);
EleStruct multipole_init(int who, std::optional<bool> zero = std::nullopt);
extern "C" void fortran_multipole_kick(
    double &knl /* 0D_NOT_real in */,
    double &tilt /* 0D_NOT_real in */,
    int &n /* 0D_NOT_integer in */,
    int &ref_species /* 0D_NOT_integer in */,
    int &ele_orientation /* 0D_NOT_integer in */,
    void *coord /* 0D_NOT_type inout */,
    int *pole_type /* 0D_NOT_integer in */,
    bool *ref_orb_offset /* 0D_NOT_logical in */
);
void multipole_kick(
    double knl,
    double tilt,
    int n,
    int ref_species,
    int ele_orientation,
    CoordStruct &coord,
    std::optional<int> pole_type = std::nullopt,
    std::optional<bool> ref_orb_offset = std::nullopt
);
extern "C" void fortran_multipole_kick_mat(
    Bmad::array_descriptor_t &knl /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &tilt /* 1D_NOT_real in */,
    int &ref_species /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    double &factor /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */
);
FixedArray2D<Real, 6, 6> multipole_kick_mat(
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    int ref_species,
    EleStruct &ele,
    CoordStruct &orbit,
    double factor
);
extern "C" void fortran_multipole_kicks(
    Bmad::array_descriptor_t &knl /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &tilt /* 1D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    int *pole_type /* 0D_NOT_integer in */,
    bool *ref_orb_offset /* 0D_NOT_logical in */
);
void multipole_kicks(
    FArray1D<Real> &knl,
    FArray1D<Real> &tilt,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<int> pole_type = std::nullopt,
    std::optional<bool> ref_orb_offset = std::nullopt
);
extern "C" void fortran_multipole_kt_to_ab(
    Bmad::array_descriptor_t &knl /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &knsl /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &tn /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &an /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &bn /* 1D_NOT_real inout */
);
void multipole_kt_to_ab(
    FArray1D<Real> &knl,
    FArray1D<Real> &knsl,
    FArray1D<Real> &tn,
    FArray1D<Real> &an,
    FArray1D<Real> &bn
);
extern "C" void fortran_multipole_spin_tracking(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void multipole_spin_tracking(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" bool fortran_mytan(
    double &y /* 0D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    double &arg /* 0D_NOT_real in */
);
void mytan(double y, double x, double arg);
extern "C" bool fortran_n_attrib_string_max_len(int &max_len /* 0D_NOT_integer out */);
int n_attrib_string_max_len();
extern "C" void fortran_new_control(
    void *lat /* 0D_NOT_type in */,
    int &ix_ele /* 0D_NOT_integer out */,
    const char *ele_name /* 0D_NOT_character in */
);
int new_control(LatStruct &lat, std::optional<std::string> ele_name = std::nullopt);
extern "C" bool
fortran_nint_chk(double &re_val /* 0D_NOT_real in */, int &int_val /* 0D_NOT_integer out */);
int nint_chk(double re_val);
extern "C" void fortran_normal_form_complex_taylors(
    Bmad::array_descriptor_t &one_turn_taylor /* 1D_NOT_type inout */,
    bool &rf_on /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &F /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &L /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &A /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &A_inverse /* 1D_NOT_type inout */,
    int *order /* 0D_NOT_integer in */
);
void normal_form_complex_taylors(
    TaylorStructArray1D one_turn_taylor,
    bool rf_on,
    std::optional<ComplexTaylorStructArray1D> F = std::nullopt,
    std::optional<ComplexTaylorStructArray1D> L = std::nullopt,
    std::optional<TaylorStructArray1D> A = std::nullopt,
    std::optional<TaylorStructArray1D> A_inverse = std::nullopt,
    std::optional<int> order = std::nullopt
);

// Skipped unusable routine normal_form_rd_terms:
// - Untranslated type: probe_8 (0D)
extern "C" void fortran_normal_form_taylors(
    Bmad::array_descriptor_t &one_turn_taylor /* 1D_NOT_type in */,
    bool &rf_on /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &dhdj /* 1D_NOT_type out */,
    Bmad::array_descriptor_t &A /* 1D_NOT_type out */,
    Bmad::array_descriptor_t &A_inverse /* 1D_NOT_type out */
);
struct NormalFormTaylors {
  TaylorStructArray1D dhdj;
  TaylorStructArray1D A;
  TaylorStructArray1D A_inverse;
};
Bmad::NormalFormTaylors normal_form_taylors(TaylorStructArray1D one_turn_taylor, bool rf_on);
extern "C" void fortran_normal_mode3_calc(
    Bmad::array_descriptor_t &t6 /* 2D_NOT_real inout */,
    Bmad::array_descriptor_t &tune /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &B /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &HV /* 2D_NOT_real out */,
    bool *above_transition /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real in */
);
struct NormalMode3Calc {
  FixedArray1D<Real, 3> tune;
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> HV;
};
Bmad::NormalMode3Calc normal_mode3_calc(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<bool> above_transition = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt
);
extern "C" void fortran_normal_mode_dispersion(
    void *ele /* 0D_NOT_type inout */,
    bool *reverse /* 0D_NOT_logical in */
);
void normal_mode_dispersion(EleStruct &ele, std::optional<bool> reverse = std::nullopt);
extern "C" void fortran_normalize_evecs(
    Bmad::array_descriptor_t &evec /* 2D_NOT_complex inout */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool normalize_evecs(FixedArray2D<Complex, 6, 6> evec);
extern "C" bool
fortran_num_field_eles(void *ele /* 0D_NOT_type in */, int &n_field_ele /* 0D_NOT_integer out */);
int num_field_eles(EleStruct &ele);
extern "C" bool fortran_num_lords(
    void *slave /* 0D_NOT_type in */,
    int &lord_type /* 0D_NOT_integer in */,
    int &num /* 0D_NOT_integer out */
);
int num_lords(EleStruct &slave, int lord_type);
extern "C" void fortran_odeint_bmad(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &s1_body /* 0D_NOT_real in */,
    double &s2_body /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
struct OdeintBmad {
  bool err_flag;
  TrackStruct track;
};
Bmad::OdeintBmad odeint_bmad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    double s1_body,
    double s2_body,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_odeint_bmad_time(
    void *orb /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &t_dir /* 0D_NOT_integer in */,
    double &rf_time /* 0D_NOT_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type inout */,
    double *t_end /* 0D_NOT_real in */,
    double &dt_step /* 0D_NOT_real out */,
    void *extra_field /* 0D_NOT_type in */
);
struct OdeintBmadTime {
  bool err_flag;
  double dt_step;
};
Bmad::OdeintBmadTime odeint_bmad_time(
    CoordStruct &orb,
    EleStruct &ele,
    LatParamStruct &param,
    int t_dir,
    double &rf_time,
    optional_ref<TrackStruct> track = std::nullopt,
    std::optional<double> t_end = std::nullopt,
    optional_ref<EmFieldStruct> extra_field = std::nullopt
);
extern "C" void fortran_offset_particle(
    void *ele /* 0D_NOT_type in */,
    bool &set /* 0D_NOT_logical in */,
    void *orbit /* 0D_NOT_type inout */,
    bool *set_tilt /* 0D_NOT_logical in */,
    bool *set_hvkicks /* 0D_NOT_logical in */,
    int *drift_to_edge /* 0D_NOT_integer in */,
    double *s_pos /* 0D_NOT_real in */,
    double &s_out /* 0D_NOT_real out */,
    bool *set_spin /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin_qrot /* 1D_NOT_real out */,
    double *time /* 0D_NOT_real inout */
);
struct OffsetParticle {
  double s_out;
  FixedArray1D<Real, 4> spin_qrot;
};
Bmad::OffsetParticle offset_particle(
    EleStruct &ele,
    bool set,
    CoordStruct &orbit,
    std::optional<bool> set_tilt = std::nullopt,
    std::optional<bool> set_hvkicks = std::nullopt,
    std::optional<int> drift_to_edge = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> set_spin = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    optional_ref<double> time = std::nullopt
);
extern "C" void fortran_offset_photon(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    bool &set /* 0D_NOT_logical in */,
    bool *offset_position_only /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &rot_mat /* 2D_NOT_real in */
);
void offset_photon(
    EleStruct &ele,
    CoordStruct &orbit,
    bool set,
    std::optional<bool> offset_position_only = std::nullopt,
    std::optional<FixedArray2D<Real, 3, 3>> rot_mat = std::nullopt
);
extern "C" void fortran_one_turn_mat_at_ele(
    void *ele /* 0D_NOT_type in */,
    double &phi_a /* 0D_NOT_real in */,
    double &phi_b /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat4 /* 2D_NOT_real out */
);
FixedArray2D<Real, 4, 4> one_turn_mat_at_ele(EleStruct &ele, double phi_a, double phi_b);
extern "C" bool fortran_open_binary_file(
    const char *file_name /* 0D_NOT_character in */,
    const char *action /* 0D_NOT_character in */,
    int &iu /* 0D_NOT_integer out */,
    const char *r_name /* 0D_NOT_character in */,
    int &iver /* 0D_NOT_integer out */,
    bool &is_ok /* 0D_NOT_logical out */
);
struct OpenBinaryFile {
  int iu;
  int iver;
  bool is_ok;
};
Bmad::OpenBinaryFile
open_binary_file(std::string file_name, std::string action, std::string r_name);
extern "C" void fortran_orbit_amplitude_calc(
    void *ele /* 0D_NOT_type in */,
    void *orb /* 0D_NOT_type in */,
    double &amp_a /* 0D_NOT_real out */,
    double &amp_b /* 0D_NOT_real out */,
    double &amp_na /* 0D_NOT_real out */,
    double &amp_nb /* 0D_NOT_real out */
);
struct OrbitAmplitudeCalc {
  double amp_a;
  double amp_b;
  double amp_na;
  double amp_nb;
};
Bmad::OrbitAmplitudeCalc orbit_amplitude_calc(EleStruct &ele, CoordStruct &orb);
extern "C" void fortran_orbit_reference_energy_correction(
    void *orbit /* 0D_NOT_type inout */,
    double &p0c_new /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void orbit_reference_energy_correction(
    CoordStruct &orbit,
    double p0c_new,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_orbit_to_floor_phase_space(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &floor_phase_space /* 1D_NOT_real out */
);
FixedArray1D<Real, 6> orbit_to_floor_phase_space(CoordStruct &orbit, EleStruct &ele);
extern "C" bool fortran_orbit_to_local_curvilinear(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    int *z_direction /* 0D_NOT_integer in */,
    int *relative_to /* 0D_NOT_integer in */,
    void *local_position /* 0D_NOT_type out */
);
FloorPositionStruct orbit_to_local_curvilinear(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<int> z_direction = std::nullopt,
    std::optional<int> relative_to = std::nullopt
);
extern "C" bool fortran_orbit_too_large(
    void *orbit /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type out */,
    bool *check_momentum /* 0D_NOT_logical in */,
    bool &is_too_large /* 0D_NOT_logical out */
);
struct OrbitTooLarge {
  LatParamStruct param;
  bool is_too_large;
};
Bmad::OrbitTooLarge
orbit_too_large(CoordStruct &orbit, std::optional<bool> check_momentum = std::nullopt);
extern "C" void fortran_order_evecs_by_n_similarity(
    Bmad::array_descriptor_t &evec /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &eval /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &mat_tunes /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &Nmat /* 2D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct OrderEvecsByNSimilarity {
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
};
Bmad::OrderEvecsByNSimilarity order_evecs_by_n_similarity(
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray2D<Real, 6, 6> Nmat
);
extern "C" void fortran_order_evecs_by_plane_dominance(
    Bmad::array_descriptor_t &evec /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &eval /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &mat_tunes /* 1D_NOT_real inout */
);
void order_evecs_by_plane_dominance(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    std::optional<FixedArray1D<Real, 3>> mat_tunes = std::nullopt
);
extern "C" void fortran_order_evecs_by_tune(
    Bmad::array_descriptor_t &evec /* 2D_NOT_complex inout */,
    Bmad::array_descriptor_t &eval /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &mat_tunes /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool order_evecs_by_tune(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray1D<Real, 3> abz_tunes
);
extern "C" void fortran_order_particles_in_z(void *bunch /* 0D_NOT_type inout */);
void order_particles_in_z(BunchStruct &bunch);
extern "C" void fortran_order_super_lord_slaves(
    void *lat /* 0D_NOT_type inout */,
    int &ix_lord /* 0D_NOT_integer in */
);
void order_super_lord_slaves(LatStruct &lat, int ix_lord);
extern "C" void fortran_osc_alloc_freespace_array(
    Bmad::array_descriptor_t &nlo /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &nhi /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &npad /* 1D_NOT_integer in */
);
void osc_alloc_freespace_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
);
extern "C" void fortran_osc_alloc_image_array(
    Bmad::array_descriptor_t &nlo /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &nhi /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &npad /* 1D_NOT_integer in */
);
void osc_alloc_image_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
);
extern "C" void fortran_osc_alloc_rectpipe_arrays(
    Bmad::array_descriptor_t &nlo /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &nhi /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &npad /* 1D_NOT_integer in */
);
void osc_alloc_rectpipe_arrays(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad
);

// Skipped unusable routine osc_cathodeimages_solver:
// - Variable inout sized array: 3D_NOT_real
// - Variable inout sized array: 3D_NOT_real
// - Variable inout sized array: 4D_NOT_real
// - Variable inout sized array: 4D_NOT_real

// Skipped unusable routine osc_freespace_solver:
// - Variable in sized array: 3D_NOT_real
// - Variable inout sized array: 3D_NOT_real
// - Variable inout sized array: 4D_NOT_real
// - Variable inout sized array: 4D_NOT_real

// Skipped unusable routine osc_freespace_solver2:
// - Variable in sized array: 3D_NOT_real
// - Variable inout sized array: 4D_NOT_real
// - Variable inout sized array: 3D_NOT_real

// Skipped unusable routine osc_get_cgrn_freespace:
// - Variable inout sized array: 3D_NOT_complex

// Skipped unusable routine osc_getgrnfree:
// - Array bounds handling: "Enum 'ILO_GRN' found in bounds 'ilo_grn' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageconvcorr:
// - Array bounds handling: "Enum 'ILO_GRN' found in bounds 'ilo_grn' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageshift:
// - Array bounds handling: "Enum 'ILO_GRN' found in bounds 'ilo_grn' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_osc_getgrnpipe(
    double &gam /* 0D_NOT_real in */,
    double &a /* 0D_NOT_real in */,
    double &b /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &delta /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &umin /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &npad /* 1D_NOT_integer inout */
);
void osc_getgrnpipe(
    double gam,
    double a,
    double b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad
);
extern "C" void fortran_osc_read_rectpipe_grn();
void osc_read_rectpipe_grn();

// Skipped unusable routine osc_rectpipe_solver:
// - Variable in sized array: 3D_NOT_real
// - Variable inout sized array: 3D_NOT_real
// - Variable inout sized array: 4D_NOT_real
// - Variable inout sized array: 4D_NOT_real
extern "C" void fortran_osc_write_rectpipe_grn(
    double &apipe /* 0D_NOT_real in */,
    double &bpipe /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &delta /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &umin /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &umax /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &nlo /* 1D_NOT_integer inout */,
    Bmad::array_descriptor_t &nhi /* 1D_NOT_integer inout */,
    double &gamma /* 0D_NOT_real in */
);
void osc_write_rectpipe_grn(
    double apipe,
    double bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double gamma
);
extern "C" void fortran_parse_cartesian_map(
    void *ct_map /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parse_cartesian_map(
    CartesianMapStruct &ct_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
);
extern "C" void fortran_parse_cylindrical_map(
    void *cl_map /* 0D_PTR_type inout */,
    void *ele /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parse_cylindrical_map(
    CylindricalMapStruct &cl_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
);
extern "C" void fortran_parse_gen_grad_map(
    void *gg_map /* 0D_PTR_type inout */,
    void *ele /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parse_gen_grad_map(
    GenGradMapStruct &gg_map,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
);
extern "C" void fortran_parse_grid_field(
    void *g_field /* 0D_PTR_type inout */,
    void *ele /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parse_grid_field(
    GridFieldStruct &g_field,
    EleStruct &ele,
    LatStruct &lat,
    std::string delim,
    bool delim_found,
    bool err_flag
);
extern "C" bool fortran_parse_integer_list(
    const char *err_str /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &int_array /* 1D_NOT_integer inout */,
    bool &exact_size /* 0D_NOT_logical in */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    const char *open_delim /* 0D_NOT_character in */,
    const char *separator /* 0D_NOT_character in */,
    const char *close_delim /* 0D_NOT_character in */,
    int *default_value /* 0D_NOT_integer in */,
    bool &is_ok /* 0D_NOT_logical in */
);
void parse_integer_list(
    std::string err_str,
    LatStruct &lat,
    FArray1D<Int> &int_array,
    bool exact_size,
    std::string delim,
    bool delim_found,
    bool is_ok,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<int> default_value = std::nullopt
);
extern "C" bool fortran_parse_integer_list2(
    const char *err_str /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    void *int_array /* 1D_ALLOC_integer inout */,
    int &num_found /* 0D_NOT_integer out */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    int *num_expected /* 0D_NOT_integer in */,
    const char *open_delim /* 0D_NOT_character in */,
    const char *separator /* 0D_NOT_character in */,
    const char *close_delim /* 0D_NOT_character in */,
    int *default_value /* 0D_NOT_integer in */,
    bool &is_ok /* 0D_NOT_logical out */
);
struct ParseIntegerList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
Bmad::ParseIntegerList2 parse_integer_list2(
    std::string err_str,
    LatStruct &lat,
    IntAlloc1D &int_array,
    std::optional<int> num_expected = std::nullopt,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<int> default_value = std::nullopt
);

// Skipped unusable routine parse_line_or_list:
// - Untranslated type: seq_struct (1D)
extern "C" bool fortran_parse_real_list(
    void *lat /* 0D_NOT_type in */,
    const char *err_str /* 0D_NOT_character in */,
    Bmad::array_descriptor_t &real_array /* 1D_NOT_real inout */,
    bool &exact_size /* 0D_NOT_logical in */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    const char *open_delim /* 0D_NOT_character in */,
    const char *separator /* 0D_NOT_character in */,
    const char *close_delim /* 0D_NOT_character in */,
    double *default_value /* 0D_NOT_real in */,
    int &num_found /* 0D_NOT_integer out */,
    bool &is_ok /* 0D_NOT_logical in */
);
struct ParseRealList {
  std::string delim;
  bool delim_found;
  int num_found;
};
Bmad::ParseRealList parse_real_list(
    LatStruct &lat,
    std::string err_str,
    FArray1D<Real> &real_array,
    bool exact_size,
    bool is_ok,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<double> default_value = std::nullopt
);
extern "C" bool fortran_parse_real_list2(
    void *lat /* 0D_NOT_type in */,
    const char *err_str /* 0D_NOT_character in */,
    void *real_array /* 1D_ALLOC_real inout */,
    int &num_found /* 0D_NOT_integer out */,
    const char *delim /* 0D_NOT_character out */,
    bool &delim_found /* 0D_NOT_logical out */,
    int *num_expected /* 0D_NOT_integer in */,
    const char *open_brace /* 0D_NOT_character in */,
    const char *separator /* 0D_NOT_character in */,
    const char *close_brace /* 0D_NOT_character in */,
    double *default_value /* 0D_NOT_real in */,
    bool *single_value /* 0D_NOT_logical in */,
    bool &is_ok /* 0D_NOT_logical out */
);
struct ParseRealList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
Bmad::ParseRealList2 parse_real_list2(
    LatStruct &lat,
    std::string err_str,
    RealAlloc1D &real_array,
    std::optional<int> num_expected = std::nullopt,
    std::optional<std::string> open_brace = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_brace = std::nullopt,
    std::optional<double> default_value = std::nullopt,
    std::optional<bool> single_value = std::nullopt
);

// Skipped unusable routine parse_real_matrix:
// - Variable in sized array: 2D_ALLOC_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine parse_superimpose_command:
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser2_add_superimpose:
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser_add_branch:
// - Untranslated type: seq_struct (1D)
// - Variable-sized inout character array: 1D_ALLOC_character
// - Untranslated type: parser_lat_struct (0D)
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_parser_add_constant(
    const char *word /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */,
    bool &redef_is_error /* 0D_NOT_logical in */
);
void parser_add_constant(std::string word, LatStruct &lat, bool redef_is_error);

// Skipped unusable routine parser_add_lords:
// - Untranslated type: parser_lat_struct (0D)

// Skipped unusable routine parser_add_superimpose:
// - Untranslated type: parser_ele_struct (0D)
// - Untranslated type: parser_lat_struct (0D)
extern "C" void fortran_parser_call_check(
    const char *word /* 0D_NOT_character in */,
    int &ix_word /* 0D_NOT_integer in */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &call_found /* 0D_NOT_logical in */,
    bool *err_flag /* 0D_NOT_logical in */
);
void parser_call_check(
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool call_found,
    std::optional<bool> err_flag = std::nullopt
);

// Skipped unusable routine parser_debug_print_info:
// - Untranslated type: seq_struct (1D)

// Skipped unusable routine parser_error:
// - Untranslated type: seq_struct (0D)
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser_expand_line:
// - Untranslated type: seq_struct (1D)
// - Variable-sized in character array: 1D_ALLOC_character
// - Untranslated type: base_line_ele_struct (1D)
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_parser_fast_complex_read(
    Bmad::array_descriptor_t &cmplx_vec /* 1D_NOT_complex inout */,
    void *ele /* 0D_NOT_type in */,
    const char *delim /* 0D_NOT_character out */,
    const char *err_str /* 0D_NOT_character in */,
    bool &is_ok /* 0D_NOT_logical out */
);
struct ParserFastComplexRead {
  std::string delim;
  bool is_ok;
};
Bmad::ParserFastComplexRead
parser_fast_complex_read(FArray1D<Complex> &cmplx_vec, EleStruct &ele, std::string err_str);
extern "C" bool fortran_parser_fast_integer_read(
    Bmad::array_descriptor_t &int_vec /* 1D_NOT_integer inout */,
    void *ele /* 0D_NOT_type inout */,
    const char *delim_wanted /* 0D_NOT_character in */,
    const char *err_str /* 0D_NOT_character in */,
    bool &is_ok /* 0D_NOT_logical in */
);
void parser_fast_integer_read(
    FArray1D<Int> &int_vec,
    EleStruct &ele,
    std::string delim_wanted,
    std::string err_str,
    bool is_ok
);
extern "C" bool fortran_parser_fast_real_read(
    Bmad::array_descriptor_t &real_vec /* 1D_NOT_real inout */,
    void *ele /* 0D_NOT_type in */,
    const char *end_delims /* 0D_NOT_character in */,
    const char *delim /* 0D_NOT_character out */,
    const char *err_str /* 0D_NOT_character in */,
    bool *exact_size /* 0D_NOT_logical in */,
    int &n_real /* 0D_NOT_integer out */,
    bool &is_ok /* 0D_NOT_logical out */
);
struct ParserFastRealRead {
  std::string delim;
  int n_real;
  bool is_ok;
};
Bmad::ParserFastRealRead parser_fast_real_read(
    FArray1D<Real> &real_vec,
    EleStruct &ele,
    std::string end_delims,
    std::string err_str,
    std::optional<bool> exact_size = std::nullopt
);
extern "C" void fortran_parser_file_stack(
    const char *how /* 0D_NOT_character in */,
    const char *file_name_in /* 0D_NOT_character in */,
    bool *finished /* 0D_NOT_logical in */,
    bool *err /* 0D_NOT_logical in */,
    bool *open_file /* 0D_NOT_logical in */,
    bool *abort_on_open_error /* 0D_NOT_logical in */
);
void parser_file_stack(
    std::string how,
    std::optional<std::string> file_name_in = std::nullopt,
    std::optional<bool> finished = std::nullopt,
    std::optional<bool> err = std::nullopt,
    std::optional<bool> open_file = std::nullopt,
    std::optional<bool> abort_on_open_error = std::nullopt
);
extern "C" void fortran_parser_get_integer(
    int &int_val /* 0D_NOT_integer in */,
    const char *word /* 0D_NOT_character in */,
    int &ix_word /* 0D_NOT_integer in */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err /* 0D_NOT_logical in */,
    const char *str1 /* 0D_NOT_character in */,
    const char *str2 /* 0D_NOT_character in */
);
void parser_get_integer(
    int int_val,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool err,
    std::optional<std::string> str1 = std::nullopt,
    std::optional<std::string> str2 = std::nullopt
);
extern "C" void fortran_parser_get_logical(
    const char *attrib_name /* 0D_NOT_character in */,
    bool &this_logic /* 0D_NOT_logical in */,
    const char *ele_name /* 0D_NOT_character in */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err /* 0D_NOT_logical in */
);
void parser_get_logical(
    std::string attrib_name,
    bool this_logic,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool err
);
extern "C" void fortran_parser_identify_fork_to_element(void *lat /* 0D_NOT_type inout */);
void parser_identify_fork_to_element(LatStruct &lat);
extern "C" void fortran_parser_init_custom_elements(void *lat /* 0D_NOT_type inout */);
void parser_init_custom_elements(LatStruct &lat);
extern "C" void fortran_parser_print_line(
    void *lat /* 0D_NOT_type inout */,
    bool &end_of_file /* 0D_NOT_logical in */
);
void parser_print_line(LatStruct &lat, bool end_of_file);
extern "C" void fortran_parser_read_lr_wake(
    void *ele /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parser_read_lr_wake(EleStruct &ele, std::string delim, bool delim_found, bool err_flag);
extern "C" void fortran_parser_read_old_format_lr_wake(
    void *ele /* 0D_NOT_type inout */,
    const char *lr_file_name /* 0D_NOT_character in */
);
void parser_read_old_format_lr_wake(EleStruct &ele, std::string lr_file_name);
extern "C" void fortran_parser_read_old_format_sr_wake(
    void *ele /* 0D_NOT_type inout */,
    const char *sr_file_name /* 0D_NOT_character in */
);
void parser_read_old_format_sr_wake(EleStruct &ele, std::string sr_file_name);
extern "C" void fortran_parser_read_sr_wake(
    void *ele /* 0D_NOT_type inout */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void parser_read_sr_wake(EleStruct &ele, std::string delim, bool delim_found, bool err_flag);

// Skipped unusable routine parser_set_attribute:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_parser_transfer_control_struct(
    void *con_in /* 0D_NOT_type in */,
    void *con_out /* 0D_NOT_type out */,
    void *lord /* 0D_NOT_type in */,
    int &ix_var /* 0D_NOT_integer in */
);
ControlStruct parser_transfer_control_struct(ControlStruct &con_in, EleStruct &lord, int ix_var);
extern "C" bool fortran_particle_in_global_frame(
    void *orb /* 0D_NOT_type in */,
    void *branch /* 0D_NOT_type in */,
    bool *in_time_coordinates /* 0D_NOT_logical in */,
    bool *in_body_frame /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &w_mat_out /* 2D_NOT_real inout */,
    void *particle /* 0D_NOT_type out */
);
CoordStruct particle_in_global_frame(
    CoordStruct &orb,
    BranchStruct &branch,
    std::optional<bool> in_time_coordinates = std::nullopt,
    std::optional<bool> in_body_frame = std::nullopt,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat_out = std::nullopt
);
extern "C" bool fortran_particle_is_moving_backwards(
    void *orbit /* 0D_NOT_type in */,
    bool &is_moving_backwards /* 0D_NOT_logical out */
);
bool particle_is_moving_backwards(CoordStruct &orbit);
extern "C" bool fortran_particle_is_moving_forward(
    void *orbit /* 0D_NOT_type in */,
    int *dir /* 0D_NOT_integer in */,
    bool &is_moving_forward /* 0D_NOT_logical out */
);
bool particle_is_moving_forward(CoordStruct &orbit, std::optional<int> dir = std::nullopt);
extern "C" bool fortran_particle_rf_time(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    bool *reference_active_edge /* 0D_NOT_logical in */,
    double *s_rel /* 0D_NOT_real in */,
    bool *time_coords /* 0D_NOT_logical in */,
    double *rf_freq /* 0D_NOT_real in */,
    bool *abs_time /* 0D_NOT_logical in */,
    long double &time /* 0D_NOT_real16 out */
);
long double particle_rf_time(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<bool> reference_active_edge = std::nullopt,
    std::optional<double> s_rel = std::nullopt,
    std::optional<bool> time_coords = std::nullopt,
    std::optional<double> rf_freq = std::nullopt,
    std::optional<bool> abs_time = std::nullopt
);
extern "C" bool fortran_patch_flips_propagation_direction(
    double &x_pitch /* 0D_NOT_real in */,
    double &y_pitch /* 0D_NOT_real in */,
    bool &is_flip /* 0D_NOT_logical out */
);
bool patch_flips_propagation_direction(double x_pitch, double y_pitch);
extern "C" bool fortran_patch_length(
    void *patch /* 0D_NOT_type in */,
    int *ref_coords /* 0D_NOT_integer in */,
    double &length /* 0D_NOT_real out */
);
double patch_length(EleStruct &patch, std::optional<int> ref_coords = std::nullopt);
extern "C" void fortran_photon_absorption_and_phase_shift(
    const char *material /* 0D_NOT_character in */,
    double &Energy /* 0D_NOT_real in */,
    double &absorption /* 0D_NOT_real out */,
    double &phase_shift /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct PhotonAbsorptionAndPhaseShift {
  double absorption;
  double phase_shift;
  bool err_flag;
};
Bmad::PhotonAbsorptionAndPhaseShift
photon_absorption_and_phase_shift(std::string material, double Energy);
extern "C" void fortran_photon_add_to_detector_statistics(
    void *orbit0 /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type inout */,
    int &ix_pt /* 0D_NOT_integer out */,
    int &iy_pt /* 0D_NOT_integer out */,
    void *pixel_pt /* 0D_NOT_type in */
);
struct PhotonAddToDetectorStatistics {
  int ix_pt;
  int iy_pt;
};
Bmad::PhotonAddToDetectorStatistics photon_add_to_detector_statistics(
    CoordStruct &orbit0,
    CoordStruct &orbit,
    EleStruct &ele,
    optional_ref<PixelPtStruct> pixel_pt = std::nullopt
);

// Skipped unusable routine photon_diffuse_scattering:
// - Untranslated type: diffuse_param_struct (0D)

// Skipped unusable routine photon_read_spline:
// - Untranslated type: photon_init_splines_struct (0D)
extern "C" void fortran_photon_reflection(
    double &graze_angle_in /* 0D_NOT_real in */,
    double &energy /* 0D_NOT_real in */,
    void *surface /* 0D_NOT_type in */,
    double &graze_angle_out /* 0D_NOT_real out */,
    double &phi_out /* 0D_NOT_real out */
);
struct PhotonReflection {
  double graze_angle_out;
  double phi_out;
};
Bmad::PhotonReflection
photon_reflection(double graze_angle_in, double energy, PhotonReflectSurfaceStruct &surface);
extern "C" void fortran_photon_reflection_std_surface_init(void *surface /* 0D_NOT_type out */);
PhotonReflectSurfaceStruct photon_reflection_std_surface_init();
extern "C" void fortran_photon_reflectivity(
    double &angle /* 0D_NOT_real in */,
    double &energy /* 0D_NOT_real in */,
    void *surface /* 0D_NOT_type in */,
    double &p_reflect /* 0D_NOT_real out */,
    double &rel_p_specular /* 0D_NOT_real out */
);
struct PhotonReflectivity {
  double p_reflect;
  double rel_p_specular;
};
Bmad::PhotonReflectivity
photon_reflectivity(double angle, double energy, PhotonReflectSurfaceStruct &surface);
extern "C" void fortran_photon_target_corner_calc(
    void *aperture_ele /* 0D_NOT_type in */,
    double &x_lim /* 0D_NOT_real in */,
    double &y_lim /* 0D_NOT_real in */,
    double &z_lim /* 0D_NOT_real in */,
    void *source_ele /* 0D_NOT_type in */,
    void *corner /* 0D_NOT_type out */
);
TargetPointStruct photon_target_corner_calc(
    EleStruct &aperture_ele,
    double x_lim,
    double y_lim,
    double z_lim,
    EleStruct &source_ele
);
extern "C" void fortran_photon_target_setup(void *ele /* 0D_NOT_type inout */);
void photon_target_setup(EleStruct &ele);
extern "C" bool
fortran_photon_type(void *ele /* 0D_NOT_type in */, int &e_type /* 0D_NOT_integer out */);
int photon_type(EleStruct &ele);
extern "C" bool fortran_physical_ele_end(
    int &track_end /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type in */,
    int &ele_orientation /* 0D_NOT_integer in */,
    bool *return_stream_end /* 0D_NOT_logical in */,
    int &physical_end /* 0D_NOT_integer out */
);
int physical_ele_end(
    int track_end,
    CoordStruct &orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end = std::nullopt
);
extern "C" void fortran_point_photon_emission(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    int &direction /* 0D_NOT_integer in */,
    double &max_target_area /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &w_to_surface /* 2D_NOT_real in */
);
void point_photon_emission(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    int direction,
    double max_target_area,
    std::optional<FixedArray2D<Real, 3, 3>> w_to_surface = std::nullopt
);

// Skipped unusable routine pointer_to_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_pointer_to_branch_given_ele(
    void *ele /* 0D_NOT_type in */,
    void *branch_ptr /* 0D_PTR_type out */
);
std::optional<BranchStruct> pointer_to_branch(EleStruct &ele);
extern "C" bool fortran_pointer_to_branch_given_name(
    const char *branch_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    bool *parameter_is_branch0 /* 0D_NOT_logical in */,
    int *blank_branch /* 0D_NOT_integer in */,
    void *branch_ptr /* 0D_PTR_type out */
);
std::optional<BranchStruct> pointer_to_branch(
    std::string branch_name,
    LatStruct &lat,
    std::optional<bool> parameter_is_branch0 = std::nullopt,
    std::optional<int> blank_branch = std::nullopt
);
extern "C" bool fortran_pointer_to_ele1(
    void *lat /* 0D_NOT_type in */,
    int &ix_ele /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    void *ele_ptr /* 0D_PTR_type out */
);
std::optional<EleStruct>
pointer_to_ele(LatStruct &lat, int ix_ele, std::optional<int> ix_branch = std::nullopt);
extern "C" bool fortran_pointer_to_ele2(
    void *lat /* 0D_NOT_type in */,
    void *ele_loc /* 0D_NOT_type in */,
    void *ele_ptr /* 0D_PTR_type out */
);
std::optional<EleStruct> pointer_to_ele(LatStruct &lat, LatEleLocStruct &ele_loc);
extern "C" bool fortran_pointer_to_ele3(
    void *lat /* 0D_NOT_type in */,
    const char *ele_name /* 0D_NOT_character in */,
    void *ele_ptr /* 0D_PTR_type out */
);
std::optional<EleStruct> pointer_to_ele(LatStruct &lat, std::string ele_name);
extern "C" bool fortran_pointer_to_ele4(
    void *lat /* 0D_NOT_type in */,
    void *foreign_ele /* 0D_NOT_type in */,
    void *ele_ptr /* 0D_PTR_type out */
);
std::optional<EleStruct> pointer_to_ele(LatStruct &lat, EleStruct &foreign_ele);

// Skipped unusable routine pointer_to_ele_multipole:
// - Routine in configuration skip list
extern "C" bool fortran_pointer_to_element_at_s(
    void *branch /* 0D_NOT_type in */,
    double &s /* 0D_NOT_real in */,
    bool &choose_max /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &s_eff /* 0D_NOT_real out */,
    void *position /* 0D_NOT_type out */,
    bool *print_err /* 0D_NOT_logical in */,
    void *ele /* 0D_PTR_type out */
);
struct PointerToElementAtS {
  bool err_flag;
  double s_eff;
  CoordStruct position;
  std::optional<EleStruct> ele;
};
Bmad::PointerToElementAtS pointer_to_element_at_s(
    BranchStruct &branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool
fortran_pointer_to_fibre(void *ele /* 0D_NOT_type in */, void *assoc_fibre /* 0D_PTR_type out */);
std::optional<Fibre> pointer_to_fibre(EleStruct &ele);
extern "C" bool fortran_pointer_to_field_ele(
    void *ele /* 0D_NOT_type in */,
    int &ix_field_ele /* 0D_NOT_integer in */,
    double &dz_offset /* 0D_NOT_real out */,
    void *field_ele /* 0D_PTR_type out */
);
struct PointerToFieldEle {
  double dz_offset;
  std::optional<EleStruct> field_ele;
};
Bmad::PointerToFieldEle pointer_to_field_ele(EleStruct &ele, int ix_field_ele);
extern "C" bool fortran_pointer_to_girder(
    void *ele /* 0D_NOT_type in */,
    int &ix_slave_back /* 0D_NOT_integer out */,
    void *girder /* 0D_PTR_type out */
);
struct PointerToGirder {
  int ix_slave_back;
  std::optional<EleStruct> girder;
};
Bmad::PointerToGirder pointer_to_girder(EleStruct &ele);

// Skipped unusable routine pointer_to_indexed_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_pointer_to_lord(
    void *slave /* 0D_NOT_type in */,
    int &ix_lord /* 0D_NOT_integer in */,
    void *control /* 0D_PTR_type out */,
    int &ix_slave_back /* 0D_NOT_integer out */,
    int *lord_type /* 0D_NOT_integer in */,
    int &ix_control /* 0D_NOT_integer out */,
    int &ix_ic /* 0D_NOT_integer out */,
    void *lord_ptr /* 0D_PTR_type out */
);
struct PointerToLord {
  std::optional<ControlStruct> control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
  std::optional<EleStruct> lord_ptr;
};
Bmad::PointerToLord
pointer_to_lord(EleStruct &slave, int ix_lord, std::optional<int> lord_type = std::nullopt);
extern "C" bool fortran_pointer_to_multipass_lord(
    void *ele /* 0D_NOT_type in */,
    int &ix_pass /* 0D_NOT_integer out */,
    void *super_lord /* 0D_PTR_type out */,
    void *multi_lord /* 0D_PTR_type out */
);
struct PointerToMultipassLord {
  int ix_pass;
  std::optional<EleStruct> super_lord;
  std::optional<EleStruct> multi_lord;
};
Bmad::PointerToMultipassLord pointer_to_multipass_lord(EleStruct &ele);
extern "C" bool fortran_pointer_to_next_ele(
    void *this_ele /* 0D_NOT_type in */,
    int *offset /* 0D_NOT_integer in */,
    bool *skip_beginning /* 0D_NOT_logical in */,
    bool *follow_fork /* 0D_NOT_logical in */,
    void *next_ele /* 0D_PTR_type in */
);
void pointer_to_next_ele(
    EleStruct &this_ele,
    EleStruct &next_ele,
    std::optional<int> offset = std::nullopt,
    std::optional<bool> skip_beginning = std::nullopt,
    std::optional<bool> follow_fork = std::nullopt
);
extern "C" bool fortran_pointer_to_slave(
    void *lord /* 0D_NOT_type in */,
    int &ix_slave /* 0D_NOT_integer in */,
    void *control /* 0D_PTR_type out */,
    int *lord_type /* 0D_NOT_integer in */,
    int &ix_lord_back /* 0D_NOT_integer out */,
    int &ix_control /* 0D_NOT_integer out */,
    int &ix_ic /* 0D_NOT_integer out */,
    void *slave_ptr /* 0D_PTR_type out */
);
struct PointerToSlave {
  std::optional<ControlStruct> control;
  int ix_lord_back;
  int ix_control;
  int ix_ic;
  std::optional<EleStruct> slave_ptr;
};
Bmad::PointerToSlave
pointer_to_slave(EleStruct &lord, int ix_slave, std::optional<int> lord_type = std::nullopt);
extern "C" bool fortran_pointer_to_super_lord(
    void *slave /* 0D_NOT_type in */,
    void *control /* 0D_PTR_type out */,
    int &ix_slave_back /* 0D_NOT_integer out */,
    int &ix_control /* 0D_NOT_integer out */,
    int &ix_ic /* 0D_NOT_integer out */,
    int *lord_type /* 0D_NOT_integer in */,
    void *lord_ptr /* 0D_PTR_type out */
);
struct PointerToSuperLord {
  std::optional<ControlStruct> control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
  std::optional<EleStruct> lord_ptr;
};
Bmad::PointerToSuperLord
pointer_to_super_lord(EleStruct &slave, std::optional<int> lord_type = std::nullopt);
extern "C" bool fortran_pointer_to_surface_displacement_pt(
    void *ele /* 0D_NOT_type in */,
    bool &nearest /* 0D_NOT_logical in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    int &ix /* 0D_NOT_integer out */,
    int &iy /* 0D_NOT_integer out */,
    bool *extend_grid /* 0D_NOT_logical in */,
    double &xx /* 0D_NOT_real out */,
    double &yy /* 0D_NOT_real out */,
    void *pt /* 0D_PTR_type out */
);
struct PointerToSurfaceDisplacementPt {
  int ix;
  int iy;
  double xx;
  double yy;
  std::optional<SurfaceDisplacementPtStruct> pt;
};
Bmad::PointerToSurfaceDisplacementPt pointer_to_surface_displacement_pt(
    EleStruct &ele,
    bool nearest,
    double x,
    double y,
    std::optional<bool> extend_grid = std::nullopt
);
extern "C" bool fortran_pointer_to_surface_segmented_pt(
    void *ele /* 0D_NOT_type in */,
    bool &nearest /* 0D_NOT_logical in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    int &ix /* 0D_NOT_integer out */,
    int &iy /* 0D_NOT_integer out */,
    bool *extend_grid /* 0D_NOT_logical in */,
    double &xx /* 0D_NOT_real out */,
    double &yy /* 0D_NOT_real out */,
    void *pt /* 0D_PTR_type out */
);
struct PointerToSurfaceSegmentedPt {
  int ix;
  int iy;
  double xx;
  double yy;
  std::optional<SurfaceSegmentedPtStruct> pt;
};
Bmad::PointerToSurfaceSegmentedPt pointer_to_surface_segmented_pt(
    EleStruct &ele,
    bool nearest,
    double x,
    double y,
    std::optional<bool> extend_grid = std::nullopt
);
extern "C" bool fortran_pointer_to_wake_ele(
    void *ele /* 0D_NOT_type in */,
    double &delta_s /* 0D_NOT_real out */,
    void *wake_ele /* 0D_PTR_type out */
);
struct PointerToWakeEle {
  double delta_s;
  std::optional<EleStruct> wake_ele;
};
Bmad::PointerToWakeEle pointer_to_wake_ele(EleStruct &ele);
extern "C" bool fortran_pointer_to_wall3d(
    void *ele /* 0D_NOT_type in */,
    int *ix_wall /* 0D_NOT_integer in */,
    double &ds_offset /* 0D_NOT_real out */,
    bool &is_branch_wall /* 0D_NOT_logical out */,
    void *wall3d /* 0D_PTR_type out */
);
struct PointerToWall3d {
  double ds_offset;
  bool is_branch_wall;
  std::optional<Wall3dStruct> wall3d;
};
Bmad::PointerToWall3d pointer_to_wall3d(EleStruct &ele, std::optional<int> ix_wall = std::nullopt);

// Skipped unusable routine pointers_to_attribute:
// - Untranslated type: all_pointer_struct (1D)
extern "C" bool fortran_polar_to_spinor(
    void *polar /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &spinor /* 1D_NOT_complex out */
);
FixedArray1D<Complex, 2> polar_to_spinor(SpinPolarStruct &polar);
extern "C" bool fortran_polar_to_vec(
    void *polar /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &vec /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> polar_to_vec(SpinPolarStruct &polar);

// Skipped unusable routine print_mesh3d:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine prob_x_diffuse:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" void fortran_project_emit_to_xyz(
    void *ring /* 0D_NOT_type in */,
    int &ix /* 0D_NOT_integer in */,
    void *mode /* 0D_NOT_type in */,
    double &sigma_x /* 0D_NOT_real out */,
    double &sigma_y /* 0D_NOT_real out */,
    double &sigma_z /* 0D_NOT_real out */
);
struct ProjectEmitToXyz {
  double sigma_x;
  double sigma_y;
  double sigma_z;
};
Bmad::ProjectEmitToXyz project_emit_to_xyz(LatStruct &ring, int ix, NormalModesStruct &mode);

// Skipped unusable routine propagate_part_way:
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)

// Skipped unusable routine psi_prime:
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_psi_prime_sca(
    double &t /* 0D_NOT_real in */,
    double &p /* 0D_NOT_real in */,
    double &dpdt /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */
);
double psi_prime_sca(double t, double p, FixedArray1D<Real, 8> args);
extern "C" void fortran_ptc_bookkeeper(void *lat /* 0D_NOT_type inout */);
void ptc_bookkeeper(LatStruct &lat);
extern "C" void fortran_ptc_calculate_tracking_step_size(
    void *ptc_layout /* 0D_NOT_type inout */,
    double &kl_max /* 0D_NOT_real in */,
    double *ds_max /* 0D_NOT_real in */,
    void *even_steps /* 1D_ALLOC_logical in */,
    double *r_typical /* 0D_NOT_real in */,
    double *dx_tol_bend /* 0D_NOT_real in */,
    bool *use_2nd_order /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &crossover /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &crossover_wiggler /* 1D_NOT_integer in */
);
void ptc_calculate_tracking_step_size(
    Layout &ptc_layout,
    double kl_max,
    std::optional<double> ds_max = std::nullopt,
    optional_ref<BoolAlloc1D> even_steps = std::nullopt,
    std::optional<double> r_typical = std::nullopt,
    std::optional<double> dx_tol_bend = std::nullopt,
    std::optional<bool> use_2nd_order = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler = std::nullopt
);
extern "C" void fortran_ptc_check_for_lost_particle(
    int &state /* 0D_NOT_integer out */,
    void *ptc_fibre /* 0D_PTR_type out */,
    bool &do_reset /* 0D_NOT_logical in */
);
struct PtcCheckForLostParticle {
  int state;
  std::optional<Fibre> ptc_fibre;
};
Bmad::PtcCheckForLostParticle ptc_check_for_lost_particle(bool do_reset);
extern "C" void fortran_ptc_closed_orbit_calc(
    void *branch /* 0D_NOT_type in */,
    void *closed_orbit /* 1D_ALLOC_type out */,
    bool *radiation_damping_on /* 0D_NOT_logical in */
);
CoordStructAlloc1D ptc_closed_orbit_calc(
    BranchStruct &branch,
    std::optional<bool> radiation_damping_on = std::nullopt
);
extern "C" void fortran_ptc_emit_calc(
    void *ele /* 0D_NOT_type in */,
    void *norm_mode /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real inout */,
    void *closed_orb /* 0D_NOT_type out */
);
struct PtcEmitCalc {
  NormalModesStruct norm_mode;
  CoordStruct closed_orb;
};
Bmad::PtcEmitCalc ptc_emit_calc(EleStruct &ele, FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_kill_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_layouts_resplit(
    double &dKL_max /* 0D_NOT_real in */,
    double &l_max /* 0D_NOT_real in */,
    bool &l_max_drift_only /* 0D_NOT_logical in */,
    double &bend_dorb /* 0D_NOT_real in */,
    double &sex_dx /* 0D_NOT_real in */,
    bool *even /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &crossover /* 1D_NOT_integer in */,
    Bmad::array_descriptor_t &crossover_wiggler /* 1D_NOT_integer in */
);
void ptc_layouts_resplit(
    double dKL_max,
    double l_max,
    bool l_max_drift_only,
    double bend_dorb,
    double sex_dx,
    std::optional<bool> even = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler = std::nullopt
);

// Skipped unusable routine ptc_linear_isf_calc:
// - Untranslated type: linear_ele_isf_struct (1D)

// Skipped unusable routine ptc_map_to_normal_form:
// - Untranslated type: probe_8 (0D)

// Skipped unusable routine ptc_one_turn_map_at_ele:
// - Untranslated type: probe_8 (0D)
// - Untranslated type: internal_state (0D)
extern "C" void fortran_ptc_one_turn_mat_and_closed_orbit_calc(
    void *branch /* 0D_NOT_type inout */,
    double *pz /* 0D_NOT_real in */
);
void ptc_one_turn_mat_and_closed_orbit_calc(
    BranchStruct &branch,
    std::optional<double> pz = std::nullopt
);
extern "C" void fortran_ptc_ran_seed_put(int &iseed /* 0D_NOT_integer in */);
void ptc_ran_seed_put(int iseed);

// Skipped unusable routine ptc_read_flat_file:
// - Variable-sized in character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine ptc_read_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_set_rf_state_for_c_normal(bool &nocavity /* 0D_NOT_logical in */);
void ptc_set_rf_state_for_c_normal(bool nocavity);
extern "C" void fortran_ptc_set_taylor_order_if_needed();
void ptc_set_taylor_order_if_needed();

// Skipped unusable routine ptc_setup_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)

// Skipped unusable routine ptc_setup_tracking_with_damping_and_excitation:
// - Untranslated type: internal_state (0D)
extern "C" void fortran_ptc_spin_calc(
    void *ele /* 0D_NOT_type in */,
    void *norm_mode /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &sigma_mat /* 2D_NOT_real inout */,
    void *closed_orb /* 0D_NOT_type out */
);
struct PtcSpinCalc {
  NormalModesStruct norm_mode;
  CoordStruct closed_orb;
};
Bmad::PtcSpinCalc ptc_spin_calc(EleStruct &ele, FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_spin_matching_calc:
// - Untranslated type: spin_matching_struct (1D)

// Skipped unusable routine ptc_taylors_equal_bmad_taylors:
// - Untranslated type: taylor (1D)
extern "C" void fortran_ptc_track_all(
    void *branch /* 0D_NOT_type in */,
    void *orbit /* 1D_ALLOC_type inout */,
    int &track_state /* 0D_NOT_integer out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct PtcTrackAll {
  int track_state;
  bool err_flag;
};
Bmad::PtcTrackAll ptc_track_all(BranchStruct &branch, CoordStructAlloc1D orbit);

// Skipped unusable routine ptc_track_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_transfer_map_with_spin(
    void *branch /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &t_map /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &s_map /* 1D_NOT_type inout */,
    void *orb0 /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    int *ix1 /* 0D_NOT_integer in */,
    int *ix2 /* 0D_NOT_integer in */,
    bool *one_turn /* 0D_NOT_logical in */,
    bool *unit_start /* 0D_NOT_logical in */
);
bool ptc_transfer_map_with_spin(
    BranchStruct &branch,
    TaylorStructArray1D t_map,
    TaylorStructArray1D s_map,
    CoordStruct &orb0,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt
);

// Skipped unusable routine ptc_write_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)

// Skipped unusable routine ptwo:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" bool fortran_pwd_mat(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &t6 /* 2D_NOT_real in */,
    double &inductance /* 0D_NOT_real in */,
    double &sig_z /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &t6_pwd /* 2D_NOT_real out */
);
FixedArray2D<Real, 6, 6>
pwd_mat(LatStruct &lat, FixedArray2D<Real, 6, 6> t6, double inductance, double sig_z);

// Skipped unusable routine qromb_rad_int:
// - Array bounds handling: "Enum 'NUM_INT' found in bounds 'num_int' but not in provided map."
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)

// Skipped unusable routine quad_mat2_calc:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_rad1_damp_and_stoc_mats(
    void *ele /* 0D_NOT_type in */,
    bool &include_opening_angle /* 0D_NOT_logical in */,
    void *orb_in /* 0D_NOT_type in */,
    void *orb_out /* 0D_NOT_type in */,
    void *rad_map /* 0D_NOT_type out */,
    double &g2_tol /* 0D_NOT_real in */,
    double &g3_tol /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *ele0 /* 0D_NOT_type in */,
    void *rad_int1 /* 0D_NOT_type out */
);
struct Rad1DampAndStocMats {
  RadMapStruct rad_map;
  bool err_flag;
  RadInt1Struct rad_int1;
};
Bmad::Rad1DampAndStocMats rad1_damp_and_stoc_mats(
    EleStruct &ele,
    bool include_opening_angle,
    CoordStruct &orb_in,
    CoordStruct &orb_out,
    double g2_tol,
    double g3_tol,
    optional_ref<EleStruct> ele0 = std::nullopt
);
extern "C" void fortran_rad_damp_and_stoc_mats(
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool &include_opening_angle /* 0D_NOT_logical in */,
    void *rmap /* 0D_NOT_type out */,
    void *mode /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &xfer_nodamp_mat /* 2D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &closed_orbit /* 1D_NOT_type in */,
    void *rad_int_branch /* 0D_NOT_type out */
);
struct RadDampAndStocMats {
  RadMapStruct rmap;
  NormalModesStruct mode;
  FixedArray2D<Real, 6, 6> xfer_nodamp_mat;
  bool err_flag;
  RadIntBranchStruct rad_int_branch;
};
Bmad::RadDampAndStocMats rad_damp_and_stoc_mats(
    EleStruct &ele1,
    EleStruct &ele2,
    bool include_opening_angle,
    std::optional<CoordStructArray1D> closed_orbit = std::nullopt
);
extern "C" void fortran_rad_g_integrals(
    void *ele /* 0D_NOT_type in */,
    int &where /* 0D_NOT_integer in */,
    void *orb_in /* 0D_NOT_type in */,
    void *orb_out /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &int_g /* 1D_NOT_real out */,
    double &int_g2 /* 0D_NOT_real in */,
    double &int_g3 /* 0D_NOT_real out */,
    double &g_tol /* 0D_NOT_real in */,
    double &g2_tol /* 0D_NOT_real in */,
    double &g3_tol /* 0D_NOT_real in */
);
struct RadGIntegrals {
  FixedArray1D<Real, 2> int_g;
  double int_g3;
};
Bmad::RadGIntegrals rad_g_integrals(
    EleStruct &ele,
    int where,
    CoordStruct &orb_in,
    CoordStruct &orb_out,
    double int_g2,
    double g_tol,
    double g2_tol,
    double g3_tol
);
extern "C" void fortran_radiation_integrals(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &orbit /* 1D_NOT_type in */,
    void *mode /* 0D_NOT_type out */,
    int *ix_cache /* 0D_NOT_integer inout */,
    int *ix_branch /* 0D_NOT_integer in */,
    void *rad_int_by_ele /* 0D_NOT_type out */
);
struct RadiationIntegrals {
  NormalModesStruct mode;
  RadIntAllEleStruct rad_int_by_ele;
};
Bmad::RadiationIntegrals radiation_integrals(
    LatStruct &lat,
    CoordStructArray1D orbit,
    optional_ref<int> ix_cache = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
);

// Skipped unusable routine radiation_integrals_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_radiation_map_setup(
    void *ele /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *ref_orbit_in /* 0D_NOT_type inout */
);
bool radiation_map_setup(EleStruct &ele, optional_ref<CoordStruct> ref_orbit_in = std::nullopt);
extern "C" void fortran_ramper_slave_setup(
    void *lat /* 0D_NOT_type inout */,
    bool *force_setup /* 0D_NOT_logical in */
);
void ramper_slave_setup(LatStruct &lat, std::optional<bool> force_setup = std::nullopt);
extern "C" bool fortran_ramper_value(
    void *ramper /* 0D_NOT_type in */,
    void *r1 /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &value /* 0D_NOT_real out */
);
struct RamperValue {
  bool err_flag;
  double value;
};
Bmad::RamperValue ramper_value(EleStruct &ramper, ControlRamp1Struct &r1);
extern "C" void fortran_randomize_lr_wake_frequencies(
    void *ele /* 0D_NOT_type inout */,
    bool &set_done /* 0D_NOT_logical out */
);
bool randomize_lr_wake_frequencies(EleStruct &ele);

// Skipped unusable routine rb_field:
// - Array bounds handling: Calls in array bounds are not supported
extern "C" bool fortran_rchomp(
    double &rel /* 0D_NOT_real in */,
    int &plc /* 0D_NOT_integer in */,
    const char *out /* 0D_NOT_character in */
);
void rchomp(double rel, int plc, std::string out);

// Skipped unusable routine rclog_integrand:
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_re_allocate_eles(
    void *eles /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *save_old /* 0D_NOT_logical in */,
    bool *exact /* 0D_NOT_logical in */
);
void re_allocate_eles(
    ElePointerStructAlloc1D eles,
    int n,
    std::optional<bool> save_old = std::nullopt,
    std::optional<bool> exact = std::nullopt
);
extern "C" void fortran_re_allocate_wall3d_section_array(
    void *section /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void re_allocate(
    Wall3dSectionStructAlloc1D section,
    int n,
    std::optional<bool> exact = std::nullopt
);
extern "C" void fortran_re_allocate_wall3d_vertex_array(
    void *v /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void re_allocate(Wall3dVertexStructAlloc1D v, int n, std::optional<bool> exact = std::nullopt);
extern "C" void fortran_re_associate_node_array(
    void *tree /* 0D_NOT_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void re_associate_node_array(
    ExpressionTreeStruct &tree,
    int n,
    std::optional<bool> exact = std::nullopt
);
extern "C" bool fortran_re_str_qp(
    long double &rel /* 0D_NOT_real16 in */,
    const char *str_out /* 0D_NOT_character in */
);
void re_str(long double rel, std::string str_out);
extern "C" bool
fortran_re_str_rp(double &rel /* 0D_NOT_real in */, const char *str_out /* 0D_NOT_character in */);
void re_str(double rel, std::string str_out);
extern "C" void fortran_read_beam_ascii(
    const char *file_name /* 0D_NOT_character in */,
    void *beam /* 0D_NOT_type out */,
    void *beam_init /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ReadBeamAscii {
  BeamStruct beam;
  bool err_flag;
};
Bmad::ReadBeamAscii read_beam_ascii(std::string file_name, BeamInitStruct &beam_init);
extern "C" void fortran_read_beam_file(
    const char *file_name /* 0D_NOT_character in */,
    void *beam /* 0D_NOT_type out */,
    void *beam_init /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *ele /* 0D_NOT_type in */,
    bool *print_mom_shift_warning /* 0D_NOT_logical in */,
    bool *conserve_momentum /* 0D_NOT_logical in */
);
struct ReadBeamFile {
  BeamStruct beam;
  bool err_flag;
};
Bmad::ReadBeamFile read_beam_file(
    std::string file_name,
    BeamInitStruct &beam_init,
    optional_ref<EleStruct> ele = std::nullopt,
    std::optional<bool> print_mom_shift_warning = std::nullopt,
    std::optional<bool> conserve_momentum = std::nullopt
);
extern "C" void fortran_read_binary_cartesian_map(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *cart_map /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ReadBinaryCartesianMap {
  CartesianMapStruct cart_map;
  bool err_flag;
};
Bmad::ReadBinaryCartesianMap read_binary_cartesian_map(std::string file_name, EleStruct &ele);
extern "C" void fortran_read_binary_cylindrical_map(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *cl_map /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ReadBinaryCylindricalMap {
  CylindricalMapStruct cl_map;
  bool err_flag;
};
Bmad::ReadBinaryCylindricalMap read_binary_cylindrical_map(std::string file_name, EleStruct &ele);
extern "C" void fortran_read_binary_grid_field(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *g_field /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ReadBinaryGridField {
  GridFieldStruct g_field;
  bool err_flag;
};
Bmad::ReadBinaryGridField read_binary_grid_field(std::string file_name, EleStruct &ele);

// Skipped unusable routine read_digested_bmad_file:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_read_surface_reflection_file(
    const char *file_name /* 0D_NOT_character in */,
    void *surface /* 0D_NOT_type out */
);
PhotonReflectSurfaceStruct read_surface_reflection_file(std::string file_name);

// Skipped unusable routine real_8_to_taylor:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_reallocate_beam(
    void *beam /* 0D_NOT_type inout */,
    int &n_bunch /* 0D_NOT_integer in */,
    int *n_particle /* 0D_NOT_integer in */,
    bool *extend /* 0D_NOT_logical in */
);
void reallocate_beam(
    BeamStruct &beam,
    int n_bunch,
    std::optional<int> n_particle = std::nullopt,
    std::optional<bool> extend = std::nullopt
);
extern "C" void fortran_reallocate_bp_com_const();
void reallocate_bp_com_const();
extern "C" void fortran_reallocate_bunch(
    void *bunch /* 0D_NOT_type out */,
    int &n_particle /* 0D_NOT_integer in */,
    bool *save /* 0D_NOT_logical in */
);
BunchStruct reallocate_bunch(int n_particle, std::optional<bool> save = std::nullopt);
extern "C" void
fortran_reallocate_control(void *lat /* 0D_NOT_type inout */, int &n /* 0D_NOT_integer in */);
void reallocate_control(LatStruct &lat, int n);
extern "C" void fortran_reallocate_coord_array(
    void *coord_array /* 1D_ALLOC_type inout */,
    void *lat /* 0D_NOT_type in */
);
void reallocate_coord(CoordArrayStructAlloc1D coord_array, LatStruct &lat);
extern "C" void fortran_reallocate_coord_lat(
    void *coord /* 1D_ALLOC_type inout */,
    void *lat /* 0D_NOT_type in */,
    int *ix_branch /* 0D_NOT_integer in */
);
void reallocate_coord(
    CoordStructAlloc1D coord,
    LatStruct &lat,
    std::optional<int> ix_branch = std::nullopt
);
extern "C" void fortran_reallocate_coord_n(
    void *coord /* 1D_ALLOC_type inout */,
    int &n_coord /* 0D_NOT_integer in */
);
void reallocate_coord(CoordStructAlloc1D coord, int n_coord);
extern "C" void fortran_reallocate_expression_stack(
    void *stack /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void reallocate_expression_stack(
    ExpressionAtomStructAlloc1D stack,
    int n,
    std::optional<bool> exact = std::nullopt
);

// Skipped unusable routine reallocate_sequence:
// - Untranslated type: seq_struct (1D)

// Skipped unusable routine reals_8_equal_bmad_taylors:
// - Untranslated type: real_8 (1D)
extern "C" bool fortran_rel_tracking_charge_to_mass(
    void *orbit /* 0D_NOT_type in */,
    int &ref_species /* 0D_NOT_integer in */,
    double &rel_charge /* 0D_NOT_real out */
);
double rel_tracking_charge_to_mass(CoordStruct &orbit, int ref_species);
extern "C" bool fortran_relative_mode_flip(
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool &func_retval__ /* 0D_NOT_logical in */
);
void relative_mode_flip(EleStruct &ele1, EleStruct &ele2, bool func_retval__);
extern "C" void fortran_release_rad_int_cache(int &ix_cache /* 0D_NOT_integer inout */);
void release_rad_int_cache(int &ix_cache);
extern "C" void fortran_remove_constant_taylor(
    Bmad::array_descriptor_t &taylor_in /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &taylor_out /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &c0 /* 1D_NOT_real inout */,
    bool &remove_higher_order_terms /* 0D_NOT_logical in */
);
void remove_constant_taylor(
    TaylorStructArray1D taylor_in,
    TaylorStructArray1D taylor_out,
    FArray1D<Real> &c0,
    bool remove_higher_order_terms
);
extern "C" void fortran_remove_dead_from_bunch(
    void *bunch_in /* 0D_NOT_type in */,
    void *bunch_out /* 0D_NOT_type out */
);
BunchStruct remove_dead_from_bunch(BunchStruct &bunch_in);
extern "C" void fortran_remove_eles_from_lat(
    void *lat /* 0D_NOT_type inout */,
    bool *check_sanity /* 0D_NOT_logical in */
);
void remove_eles_from_lat(LatStruct &lat, std::optional<bool> check_sanity = std::nullopt);
extern "C" void fortran_remove_lord_slave_link(
    void *lord /* 0D_NOT_type inout */,
    void *slave /* 0D_NOT_type inout */
);
void remove_lord_slave_link(EleStruct &lord, EleStruct &slave);
extern "C" void fortran_reverse_lat(
    void *lat_in /* 0D_NOT_type in */,
    void *lat_rev /* 0D_NOT_type out */,
    bool *track_antiparticle /* 0D_NOT_logical in */
);
LatStruct reverse_lat(LatStruct &lat_in, std::optional<bool> track_antiparticle = std::nullopt);
extern "C" void fortran_rf_coupler_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    double &phase /* 0D_NOT_real in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void rf_coupler_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    double phase,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" bool fortran_rf_is_on(
    void *branch /* 0D_NOT_type in */,
    int *ix_ele1 /* 0D_NOT_integer in */,
    int *ix_ele2 /* 0D_NOT_integer in */,
    bool &is_on /* 0D_NOT_logical out */
);
bool rf_is_on(
    BranchStruct &branch,
    std::optional<int> ix_ele1 = std::nullopt,
    std::optional<int> ix_ele2 = std::nullopt
);
extern "C" bool fortran_rf_ref_time_offset(
    void *ele /* 0D_NOT_type in */,
    double *ds /* 0D_NOT_real in */,
    double &time /* 0D_NOT_real out */
);
double rf_ref_time_offset(EleStruct &ele, std::optional<double> ds = std::nullopt);
extern "C" bool fortran_rfun(
    double &u /* 0D_NOT_real in */,
    double &v /* 0D_NOT_real in */,
    double &w /* 0D_NOT_real in */,
    double &gam /* 0D_NOT_real in */,
    double &a /* 0D_NOT_real in */,
    double &b /* 0D_NOT_real in */,
    double &hz /* 0D_NOT_real in */,
    int &i /* 0D_NOT_integer in */,
    int &j /* 0D_NOT_integer in */,
    double &res /* 0D_NOT_real in */
);
void rfun(
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
);
extern "C" void fortran_rk_adaptive_time_step(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type inout */,
    void *orb /* 0D_NOT_type inout */,
    int &t_dir /* 0D_NOT_integer in */,
    double &rf_time /* 0D_NOT_real in */,
    double &dt_try /* 0D_NOT_real in */,
    double &dt_did /* 0D_NOT_real in */,
    double &dt_next /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical in */,
    void *extra_field /* 0D_NOT_type inout */
);
void rk_adaptive_time_step(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orb,
    int t_dir,
    double rf_time,
    double dt_try,
    double dt_did,
    double dt_next,
    bool err_flag,
    optional_ref<EmFieldStruct> extra_field = std::nullopt
);
extern "C" void fortran_rk_time_step1(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type inout */,
    double &rf_time /* 0D_NOT_real in */,
    void *orb /* 0D_NOT_type inout */,
    double &dt /* 0D_NOT_real in */,
    void *new_orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &r_err /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &dr_dt /* 1D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */,
    void *extra_field /* 0D_NOT_type inout */
);
FixedArray1D<Real, 10> rk_time_step1(
    EleStruct &ele,
    LatParamStruct &param,
    double rf_time,
    CoordStruct &orb,
    double dt,
    CoordStruct &new_orb,
    bool err_flag,
    std::optional<FixedArray1D<Real, 10>> dr_dt = std::nullopt,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EmFieldStruct> extra_field = std::nullopt
);
extern "C" bool fortran_rotate3(
    Bmad::array_descriptor_t &vec /* 1D_NOT_real inout */,
    double &angle /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &rvec /* 1D_NOT_real inout */
);
void rotate3(FixedArray1D<Real, 3> vec, double angle, FixedArray1D<Real, 3> rvec);
extern "C" void fortran_rotate_em_field(
    void *field /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &w_inv /* 2D_NOT_real in */,
    bool *calc_dfield /* 0D_NOT_logical in */,
    bool *calc_potential /* 0D_NOT_logical in */
);
void rotate_em_field(
    EmFieldStruct &field,
    FixedArray2D<Real, 3, 3> w_mat,
    FixedArray2D<Real, 3, 3> w_inv,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt
);
extern "C" void
fortran_rotate_field_zx(void *field /* 0D_NOT_type inout */, double &theta /* 0D_NOT_real in */);
void rotate_field_zx(EmFieldStruct &field, double theta);
extern "C" void fortran_rotate_for_curved_surface(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    bool &set /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &rot_mat /* 2D_NOT_real inout */
);
void rotate_for_curved_surface(
    EleStruct &ele,
    CoordStruct &orbit,
    bool set,
    FixedArray2D<Real, 3, 3> rot_mat
);
extern "C" void fortran_rotate_spin(
    Bmad::array_descriptor_t &rot_vec /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &spin /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &qrot /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> rotate_spin(FixedArray1D<Real, 3> rot_vec, FixedArray1D<Real, 3> spin);
extern "C" void fortran_rotate_spin_a_step(
    void *orbit /* 0D_NOT_type inout */,
    void *field /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    double &ds /* 0D_NOT_real in */
);
void rotate_spin_a_step(CoordStruct &orbit, EmFieldStruct &field, EleStruct &ele, double ds);
extern "C" void fortran_rotate_spin_given_field(
    void *orbit /* 0D_NOT_type inout */,
    int &sign_z_vel /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &BL /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &EL /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &qrot /* 1D_NOT_real inout */
);
void rotate_spin_given_field(
    CoordStruct &orbit,
    int sign_z_vel,
    std::optional<FixedArray1D<Real, 3>> BL = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> EL = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> qrot = std::nullopt
);
extern "C" bool fortran_s_body_calc(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    double &s_body /* 0D_NOT_real out */
);
double s_body_calc(CoordStruct &orbit, EleStruct &ele);
extern "C" void fortran_s_calc(void *lat /* 0D_NOT_type inout */);
void s_calc(LatStruct &lat);

// Skipped unusable routine s_ref_to_s_chord:
// - Untranslated type: csr_ele_info_struct (0D)

// Skipped unusable routine s_source_calc:
// - Untranslated type: csr_kick1_struct (0D)
// - Untranslated type: csr_struct (0D)
extern "C" void fortran_sad_mult_hard_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void sad_mult_hard_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_sad_soft_bend_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orb /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void sad_soft_bend_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_save_a_beam_step(
    void *ele /* 0D_NOT_type in */,
    void *beam /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &bunch_tracks /* 1D_NOT_type inout */,
    double *s_body /* 0D_NOT_real in */,
    bool *is_time_coords /* 0D_NOT_logical in */
);
void save_a_beam_step(
    EleStruct &ele,
    BeamStruct &beam,
    std::optional<BunchTrackStructArray1D> bunch_tracks = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt
);
extern "C" void fortran_save_a_bunch_step(
    void *ele /* 0D_NOT_type in */,
    void *bunch /* 0D_NOT_type in */,
    void *bunch_track /* 0D_NOT_type inout */,
    double *s_body /* 0D_NOT_real in */,
    bool *is_time_coords /* 0D_NOT_logical in */
);
void save_a_bunch_step(
    EleStruct &ele,
    BunchStruct &bunch,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt
);
extern "C" void fortran_save_a_step(
    void *track /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &local_ref_frame /* 0D_NOT_logical in */,
    void *orb /* 0D_NOT_type in */,
    double &s_rel /* 0D_NOT_real in */,
    bool *save_field /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real in */,
    bool *make_matrix /* 0D_NOT_logical in */,
    double *rf_time /* 0D_NOT_real in */,
    void *strong_beam /* 0D_NOT_type in */
);
void save_a_step(
    TrackStruct &track,
    EleStruct &ele,
    LatParamStruct &param,
    bool local_ref_frame,
    CoordStruct &orb,
    double s_rel,
    std::optional<bool> save_field = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<double> rf_time = std::nullopt,
    optional_ref<StrongBeamStruct> strong_beam = std::nullopt
);
extern "C" void fortran_sbend_body_with_k1_map(
    void *ele /* 0D_NOT_type in */,
    double &dg /* 0D_NOT_real in */,
    double &b1 /* 0D_NOT_real in */,
    void *param /* 0D_NOT_type in */,
    int &n_step /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void sbend_body_with_k1_map(
    EleStruct &ele,
    double dg,
    double b1,
    LatParamStruct &param,
    int n_step,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_sc_adaptive_step(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &include_image /* 0D_NOT_logical inout */,
    double &t_now /* 0D_NOT_real in */,
    double &dt_step /* 0D_NOT_real inout */,
    double &dt_next /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &sc_field /* 1D_NOT_type in */
);
double sc_adaptive_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool &include_image,
    double t_now,
    double &dt_step,
    EmFieldStructArray1D sc_field
);
extern "C" void fortran_sc_step(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &include_image /* 0D_NOT_logical inout */,
    double &t_end /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &sc_field /* 1D_NOT_type in */,
    int &n_emit /* 0D_NOT_integer out */
);
int sc_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool &include_image,
    double t_end,
    EmFieldStructArray1D sc_field
);
extern "C" void fortran_set_active_fixer(
    void *fixer /* 0D_NOT_type inout */,
    bool *turn_on /* 0D_NOT_logical in */,
    void *orbit /* 0D_NOT_type out */
);
CoordStruct set_active_fixer(EleStruct &fixer, std::optional<bool> turn_on = std::nullopt);

// Skipped unusable routine set_branch_and_ele_for_omp:
// - Untranslated type: lat_pointer_struct (1D)
extern "C" void fortran_set_custom_attribute_name(
    const char *custom_name /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */,
    int *custom_index /* 0D_NOT_integer in */
);
bool set_custom_attribute_name(
    std::string custom_name,
    std::optional<int> custom_index = std::nullopt
);
extern "C" void fortran_set_ele_attribute(
    void *ele /* 0D_NOT_type inout */,
    const char *set_string /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    bool *set_lords /* 0D_NOT_logical in */,
    int &err_id /* 0D_NOT_integer out */
);
struct SetEleAttribute {
  bool err_flag;
  int err_id;
};
Bmad::SetEleAttribute set_ele_attribute(
    EleStruct &ele,
    std::string set_string,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> set_lords = std::nullopt
);
extern "C" void fortran_set_ele_defaults(
    void *ele /* 0D_NOT_type inout */,
    bool *do_allocate /* 0D_NOT_logical in */
);
void set_ele_defaults(EleStruct &ele, std::optional<bool> do_allocate = std::nullopt);
extern "C" void
fortran_set_ele_name(void *ele /* 0D_NOT_type inout */, const char *name /* 0D_NOT_character in */);
void set_ele_name(EleStruct &ele, std::string name);
extern "C" void fortran_set_ele_real_attribute(
    void *ele /* 0D_NOT_type inout */,
    const char *attrib_name /* 0D_NOT_character in */,
    double &value /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *err_print_flag /* 0D_NOT_logical in */
);
bool set_ele_real_attribute(
    EleStruct &ele,
    std::string attrib_name,
    double value,
    std::optional<bool> err_print_flag = std::nullopt
);
extern "C" void fortran_set_ele_status_stale(
    void *ele /* 0D_NOT_type out */,
    int &status_group /* 0D_NOT_integer out */,
    bool &set_slaves /* 0D_NOT_logical out */
);
struct SetEleStatusStale {
  EleStruct ele;
  int status_group;
  bool set_slaves;
};
Bmad::SetEleStatusStale set_ele_status_stale();

// Skipped unusable routine set_flags_for_changed_all_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" void fortran_set_flags_for_changed_integer_attribute(
    void *ele /* 0D_NOT_type in */,
    int &attrib /* 0D_NOT_integer in */,
    bool *set_dependent /* 0D_NOT_logical in */
);
void set_flags_for_changed_attribute(
    EleStruct &ele,
    int attrib,
    std::optional<bool> set_dependent = std::nullopt
);
extern "C" void fortran_set_flags_for_changed_lat_attribute(
    void *lat /* 0D_NOT_type inout */,
    bool *set_dependent /* 0D_NOT_logical in */
);
void set_flags_for_changed_attribute(
    LatStruct &lat,
    std::optional<bool> set_dependent = std::nullopt
);
extern "C" void fortran_set_flags_for_changed_logical_attribute(
    void *ele /* 0D_NOT_type in */,
    bool &attrib /* 0D_NOT_logical in */,
    bool *set_dependent /* 0D_NOT_logical in */
);
void set_flags_for_changed_attribute(
    EleStruct &ele,
    bool attrib,
    std::optional<bool> set_dependent = std::nullopt
);
extern "C" void fortran_set_flags_for_changed_real_attribute(
    void *ele /* 0D_NOT_type in */,
    double *attrib /* 0D_NOT_real in */,
    bool *set_dependent /* 0D_NOT_logical in */
);
void set_flags_for_changed_attribute(
    EleStruct &ele,
    std::optional<double> attrib = std::nullopt,
    std::optional<bool> set_dependent = std::nullopt
);
extern "C" void fortran_set_fringe_on_off(
    double &fringe_at /* 0D_NOT_real inout */,
    int &ele_end /* 0D_NOT_integer in */,
    int &on_or_off /* 0D_NOT_integer in */
);
void set_fringe_on_off(double &fringe_at, int ele_end, int on_or_off);
extern "C" void fortran_set_lords_status_stale(
    void *ele /* 0D_NOT_type in */,
    int &stat_group /* 0D_NOT_integer in */,
    bool *control_bookkeeping /* 0D_NOT_logical in */,
    int *flag /* 0D_NOT_integer in */
);
void set_lords_status_stale(
    EleStruct &ele,
    int stat_group,
    std::optional<bool> control_bookkeeping = std::nullopt,
    std::optional<int> flag = std::nullopt
);
extern "C" void fortran_set_on_off(
    int &key /* 0D_NOT_integer in */,
    void *lat /* 0D_NOT_type inout */,
    int &switch_ /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &orb /* 1D_NOT_type in */,
    bool *use_ref_orb /* 0D_NOT_logical in */,
    int *ix_branch /* 0D_NOT_integer in */,
    void *saved_values /* 1D_ALLOC_real inout */,
    const char *attribute /* 0D_NOT_character in */,
    int *set_val /* 0D_NOT_integer in */
);
void set_on_off(
    int key,
    LatStruct &lat,
    int switch_,
    std::optional<CoordStructArray1D> orb = std::nullopt,
    std::optional<bool> use_ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<RealAlloc1D> saved_values = std::nullopt,
    std::optional<std::string> attribute = std::nullopt,
    std::optional<int> set_val = std::nullopt
);
extern "C" void fortran_set_orbit_to_zero(
    Bmad::array_descriptor_t &orbit /* 1D_NOT_type inout */,
    int &n1 /* 0D_NOT_integer in */,
    int &n2 /* 0D_NOT_integer in */,
    int *ix_noset /* 0D_NOT_integer in */
);
void set_orbit_to_zero(
    CoordStructArray1D orbit,
    int n1,
    int n2,
    std::optional<int> ix_noset = std::nullopt
);
extern "C" void fortran_set_ptc(
    double *e_tot /* 0D_NOT_real in */,
    int *particle /* 0D_NOT_integer in */,
    int *taylor_order /* 0D_NOT_integer in */,
    int *integ_order /* 0D_NOT_integer in */,
    int *n_step /* 0D_NOT_integer in */,
    bool *no_cavity /* 0D_NOT_logical in */,
    bool *force_init /* 0D_NOT_logical in */
);
void set_ptc(
    std::optional<double> e_tot = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> taylor_order = std::nullopt,
    std::optional<int> integ_order = std::nullopt,
    std::optional<int> n_step = std::nullopt,
    std::optional<bool> no_cavity = std::nullopt,
    std::optional<bool> force_init = std::nullopt
);
extern "C" void fortran_set_ptc_base_state(
    const char *component /* 0D_NOT_character in */,
    bool &set_val /* 0D_NOT_logical in */,
    bool &old_val /* 0D_NOT_logical out */
);
bool set_ptc_base_state(std::string component, bool set_val);
extern "C" void fortran_set_ptc_com_pointers();
void set_ptc_com_pointers();
extern "C" void fortran_set_ptc_quiet(
    int &channel /* 0D_NOT_integer in */,
    bool &set /* 0D_NOT_logical in */,
    int &old_val /* 0D_NOT_integer inout */
);
void set_ptc_quiet(int channel, bool set, int &old_val);
extern "C" void fortran_set_ptc_verbose(bool &on /* 0D_NOT_logical in */);
void set_ptc_verbose(bool on);
extern "C" void fortran_set_pwd_ele(
    void *lat /* 0D_NOT_type in */,
    void *mode0 /* 0D_NOT_type in */,
    double &inductance /* 0D_NOT_real in */
);
void set_pwd_ele(LatStruct &lat, NormalModesStruct &mode0, double inductance);
extern "C" void fortran_set_status_flags(
    void *bookkeeping_state /* 0D_NOT_type out */,
    int &stat /* 0D_NOT_integer in */
);
BookkeepingStateStruct set_status_flags(int stat);
extern "C" bool fortran_set_tune(
    double &phi_a_set /* 0D_NOT_real in */,
    double &phi_b_set /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &dk1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &eles /* 1D_NOT_type in */,
    void *branch /* 0D_NOT_type inout */,
    void *orb /* 1D_ALLOC_type inout */,
    bool *print_err /* 0D_NOT_logical in */,
    bool &ok /* 0D_NOT_logical out */
);
bool set_tune(
    double phi_a_set,
    double phi_b_set,
    FArray1D<Real> &dk1,
    ElePointerStructArray1D eles,
    BranchStruct &branch,
    CoordStructAlloc1D orb,
    std::optional<bool> print_err = std::nullopt
);

// Skipped unusable routine set_tune_via_group_knobs:
// - Routine in configuration skip list
extern "C" void fortran_set_twiss(
    void *branch /* 0D_NOT_type in */,
    void *twiss_ele /* 0D_NOT_type in */,
    int &ix_ele /* 0D_NOT_integer in */,
    bool &match_deta_ds /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */
);
void set_twiss(
    BranchStruct &branch,
    EleStruct &twiss_ele,
    int ix_ele,
    bool match_deta_ds,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_set_z_tune(
    void *branch /* 0D_NOT_type inout */,
    double &z_tune /* 0D_NOT_real in */,
    bool &ok /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
bool set_z_tune(BranchStruct &branch, double z_tune, std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_settable_dep_var_bookkeeping(void *ele /* 0D_NOT_type inout */);
void settable_dep_var_bookkeeping(EleStruct &ele);
extern "C" void fortran_setup_high_energy_space_charge_calc(
    bool &calc_on /* 0D_NOT_logical in */,
    void *branch /* 0D_NOT_type in */,
    double &n_part /* 0D_NOT_real in */,
    void *mode /* 0D_NOT_type in */,
    void *beam_init /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &closed_orb /* 1D_NOT_type in */
);
void setup_high_energy_space_charge_calc(
    bool calc_on,
    BranchStruct &branch,
    double n_part,
    NormalModesStruct &mode,
    optional_ref<BeamInitStruct> beam_init = std::nullopt,
    std::optional<CoordStructArray1D> closed_orb = std::nullopt
);

// Skipped unusable routine sfft:
// - Routine in configuration skip list
extern "C" void fortran_sigma_mat_ptc_to_bmad(
    Bmad::array_descriptor_t &sigma_mat_ptc /* 2D_NOT_real in */,
    double &beta0 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &sigma_mat_bmad /* 2D_NOT_real out */
);
FixedArray2D<Real, 6, 6>
sigma_mat_ptc_to_bmad(FixedArray2D<Real, 6, 6> sigma_mat_ptc, double beta0);
extern "C" bool fortran_significant_difference(
    double &value1 /* 0D_NOT_real in */,
    double &value2 /* 0D_NOT_real in */,
    double *abs_tol /* 0D_NOT_real in */,
    double *rel_tol /* 0D_NOT_real in */,
    bool &is_different /* 0D_NOT_logical out */
);
bool significant_difference(
    double value1,
    double value2,
    std::optional<double> abs_tol = std::nullopt,
    std::optional<double> rel_tol = std::nullopt
);
extern "C" bool
fortran_skip_ele_blender(void *ele /* 0D_NOT_type inout */, bool &skip /* 0D_NOT_logical in */);
void skip_ele_blender(EleStruct &ele, bool skip);
extern "C" void fortran_slice_lattice(
    void *lat /* 0D_NOT_type inout */,
    const char *ele_list /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    bool *do_bookkeeping /* 0D_NOT_logical in */
);
bool slice_lattice(
    LatStruct &lat,
    std::string ele_list,
    std::optional<bool> do_bookkeeping = std::nullopt
);
extern "C" void fortran_soft_quadrupole_edge_kick(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    int &particle_at /* 0D_NOT_integer in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void soft_quadrupole_edge_kick(
    EleStruct &ele,
    LatParamStruct &param,
    int particle_at,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_sol_quad_mat6_calc(
    double &ks_in /* 0D_NOT_real in */,
    double &k1_in /* 0D_NOT_real in */,
    double &tilt /* 0D_NOT_real in */,
    double &length /* 0D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void sol_quad_mat6_calc(
    double ks_in,
    double k1_in,
    double tilt,
    double length,
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);

// Skipped unusable routine solenoid_track_and_mat:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_solve_psi_adaptive(
    double &t0 /* 0D_NOT_real in */,
    double &t1 /* 0D_NOT_real in */,
    double &p0 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    double &p1 /* 0D_NOT_real out */
);
double solve_psi_adaptive(double t0, double t1, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_solve_psi_fixed_steps(
    double &t0 /* 0D_NOT_real in */,
    double &t1 /* 0D_NOT_real in */,
    double &p0 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &args /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &t /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &p /* 1D_NOT_real inout */
);
void solve_psi_fixed_steps(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args,
    FArray1D<Real> &t,
    FArray1D<Real> &p
);
extern "C" void fortran_sort_complex_taylor_terms(
    void *complex_taylor_in /* 0D_NOT_type in */,
    void *complex_taylor_sorted /* 0D_NOT_type out */
);
ComplexTaylorStruct sort_complex_taylor_terms(ComplexTaylorStruct &complex_taylor_in);

// Skipped unusable routine sort_universal_terms:
// - Untranslated type: universal_taylor (0D)
// - Untranslated type: universal_taylor (0D)

// Skipped unusable routine space_charge_3d:
// - Untranslated type: mesh3d_struct (0D)
// - Variable in sized array: 4D_ALLOC_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine space_charge_cathodeimages:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine space_charge_freespace:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine space_charge_rectpipe:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine spin_concat_linear_maps:
// - Routine in configuration skip list

// Skipped unusable routine spin_depolarization_rate:
// - Untranslated type: spin_matching_struct (1D)
extern "C" bool fortran_spin_dn_dpz_from_mat8(
    Bmad::array_descriptor_t &mat_1turn /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &dn_dpz_partial /* 2D_NOT_real in */,
    bool &error /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &dn_dpz /* 1D_NOT_real out */
);
struct SpinDnDpzFromMat8 {
  bool error;
  FixedArray1D<Real, 3> dn_dpz;
};
Bmad::SpinDnDpzFromMat8 spin_dn_dpz_from_mat8(
    FixedArray2D<Real, 8, 8> mat_1turn,
    std::optional<FixedArray2D<Real, 3, 3>> dn_dpz_partial = std::nullopt
);
extern "C" bool fortran_spin_dn_dpz_from_qmap(
    Bmad::array_descriptor_t &orb_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &q_map /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &dn_dpz_partial /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &dn_dpz_partial2 /* 2D_NOT_real in */,
    bool &error /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &n0 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &dn_dpz /* 1D_NOT_real out */
);
struct SpinDnDpzFromQmap {
  bool error;
  FixedArray1D<Real, 3> dn_dpz;
};
Bmad::SpinDnDpzFromQmap spin_dn_dpz_from_qmap(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> q_map,
    FixedArray2D<Real, 3, 3> dn_dpz_partial,
    FixedArray2D<Real, 3, 3> dn_dpz_partial2,
    std::optional<FixedArray1D<Real, 3>> n0 = std::nullopt
);
extern "C" void
fortran_spin_map1_normalize(Bmad::array_descriptor_t &spin1 /* 2D_NOT_real inout */);
void spin_map1_normalize(FixedArray2D<Real, 4, 7> spin1);
extern "C" void fortran_spin_mat8_resonance_strengths(
    Bmad::array_descriptor_t &orb_evec /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &mat8 /* 2D_NOT_real in */,
    double &xi_sum /* 0D_NOT_real out */,
    double &xi_diff /* 0D_NOT_real out */
);
struct SpinMat8ResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
Bmad::SpinMat8ResonanceStrengths
spin_mat8_resonance_strengths(FixedArray1D<Complex, 6> orb_evec, FixedArray2D<Real, 6, 6> mat8);
extern "C" void fortran_spin_mat_to_eigen(
    Bmad::array_descriptor_t &orb_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &spin_map /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &orb_eval /* 1D_NOT_complex out */,
    Bmad::array_descriptor_t &orb_evec /* 2D_NOT_complex out */,
    Bmad::array_descriptor_t &n0 /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &spin_evec /* 2D_NOT_complex out */,
    bool &error /* 0D_NOT_logical out */
);
struct SpinMatToEigen {
  FixedArray1D<Complex, 6> orb_eval;
  FixedArray2D<Complex, 6, 6> orb_evec;
  FixedArray1D<Real, 3> n0;
  FixedArray2D<Complex, 6, 3> spin_evec;
  bool error;
};
Bmad::SpinMatToEigen
spin_mat_to_eigen(FixedArray2D<Real, 6, 6> orb_mat, FixedArray2D<Real, 4, 7> spin_map);
extern "C" bool fortran_spin_omega(
    void *field /* 0D_NOT_type inout */,
    void *orbit /* 0D_NOT_type inout */,
    int &sign_z_vel /* 0D_NOT_integer in */,
    bool *phase_space_coords /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &omega /* 1D_NOT_real inout */
);
void spin_omega(
    EmFieldStruct &field,
    CoordStruct &orbit,
    int sign_z_vel,
    FixedArray1D<Real, 3> omega,
    std::optional<bool> phase_space_coords = std::nullopt
);
extern "C" void fortran_spin_quat_resonance_strengths(
    Bmad::array_descriptor_t &orb_evec /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &spin_q /* 2D_NOT_real in */,
    double &xi_sum /* 0D_NOT_real out */,
    double &xi_diff /* 0D_NOT_real out */
);
struct SpinQuatResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
Bmad::SpinQuatResonanceStrengths
spin_quat_resonance_strengths(FixedArray1D<Complex, 6> orb_evec, FixedArray2D<Real, 4, 7> spin_q);
extern "C" bool fortran_spin_taylor_to_linear(
    Bmad::array_descriptor_t &spin_taylor /* 1D_NOT_type in */,
    bool &normalize /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &dref_orb /* 1D_NOT_real in */,
    bool &is_on /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin_map1 /* 2D_NOT_real out */
);
FixedArray2D<Real, 4, 7> spin_taylor_to_linear(
    TaylorStructArray1D spin_taylor,
    bool normalize,
    FixedArray1D<Real, 6> dref_orb,
    bool is_on
);
extern "C" bool fortran_spinor_to_polar(
    Bmad::array_descriptor_t &spinor /* 1D_NOT_complex in */,
    void *polar /* 0D_NOT_type out */
);
SpinPolarStruct spinor_to_polar(FixedArray1D<Complex, 2> spinor);
extern "C" bool fortran_spinor_to_vec(
    Bmad::array_descriptor_t &spinor /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &vec /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> spinor_to_vec(FixedArray1D<Complex, 2> spinor);
extern "C" void fortran_spline_fit_orbit(
    void *start_orb /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &spline_x /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &spline_y /* 1D_NOT_real in */
);
void spline_fit_orbit(
    CoordStruct &start_orb,
    CoordStruct &end_orb,
    FixedArray1D<Real, 4> spline_x,
    FixedArray1D<Real, 4> spline_y
);

// Skipped unusable routine split_expression_string:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_split_lat(
    void *lat /* 0D_NOT_type inout */,
    double &s_split /* 0D_NOT_real in */,
    int &ix_branch /* 0D_NOT_integer in */,
    int &ix_split /* 0D_NOT_integer out */,
    bool &split_done /* 0D_NOT_logical out */,
    bool *add_suffix /* 0D_NOT_logical in */,
    bool *check_sanity /* 0D_NOT_logical in */,
    bool *save_null_drift /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *choose_max /* 0D_NOT_logical in */,
    int *ix_insert /* 0D_NOT_integer in */
);
struct SplitLat {
  int ix_split;
  bool split_done;
  bool err_flag;
};
Bmad::SplitLat split_lat(
    LatStruct &lat,
    double s_split,
    int ix_branch,
    std::optional<bool> add_suffix = std::nullopt,
    std::optional<bool> check_sanity = std::nullopt,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> choose_max = std::nullopt,
    std::optional<int> ix_insert = std::nullopt
);
extern "C" void fortran_sprint_spin_taylor_map(
    void *ele /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &start_orbit /* 1D_NOT_real in */
);
void sprint_spin_taylor_map(
    EleStruct &ele,
    std::optional<FixedArray1D<Real, 6>> start_orbit = std::nullopt
);
extern "C" void fortran_sr_longitudinal_wake_particle(
    void *ele /* 0D_NOT_type inout */,
    void *orbit /* 0D_NOT_type inout */
);
void sr_longitudinal_wake_particle(EleStruct &ele, CoordStruct &orbit);
extern "C" void fortran_sr_transverse_wake_particle(
    void *ele /* 0D_NOT_type inout */,
    void *orbit /* 0D_NOT_type inout */
);
void sr_transverse_wake_particle(EleStruct &ele, CoordStruct &orbit);
extern "C" void fortran_sr_z_long_wake(
    void *ele /* 0D_NOT_type in */,
    void *bunch /* 0D_NOT_type inout */,
    double &z_ave /* 0D_NOT_real in */
);
void sr_z_long_wake(EleStruct &ele, BunchStruct &bunch, double z_ave);
extern "C" void fortran_srdt_calc(
    void *lat /* 0D_NOT_type in */,
    void *srdt_sums /* 0D_NOT_type out */,
    int &order /* 0D_NOT_integer in */,
    int *n_slices_gen_opt /* 0D_NOT_integer in */,
    int *n_slices_sxt_opt /* 0D_NOT_integer in */,
    void *per_ele_out /* 1D_ALLOC_type inout */
);
SummationRdtStruct srdt_calc(
    LatStruct &lat,
    int order,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    std::optional<SummationRdtStructAlloc1D> per_ele_out = std::nullopt
);

// Skipped unusable routine srdt_calc_with_cache:
// - Variable inout sized array: 3D_ALLOC_complex
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_srdt_lsq_solution(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &var_indexes /* 1D_NOT_integer in */,
    void *ls_soln /* 1D_ALLOC_real out */,
    int *n_slices_gen_opt /* 0D_NOT_integer in */,
    int *n_slices_sxt_opt /* 0D_NOT_integer in */,
    double *chrom_set_x_opt /* 0D_NOT_real in */,
    double *chrom_set_y_opt /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &weight_in /* 1D_NOT_real in */
);
RealAlloc1D srdt_lsq_solution(
    LatStruct &lat,
    FArray1D<Int> &var_indexes,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    std::optional<double> chrom_set_x_opt = std::nullopt,
    std::optional<double> chrom_set_y_opt = std::nullopt,
    std::optional<FixedArray1D<Real, 10>> weight_in = std::nullopt
);
extern "C" void fortran_start_branch_at(
    void *lat /* 0D_NOT_type inout */,
    const char *ele_start /* 0D_NOT_character in */,
    bool &move_end_marker /* 0D_NOT_logical in */,
    bool &error /* 0D_NOT_logical out */
);
bool start_branch_at(LatStruct &lat, std::string ele_start, bool move_end_marker);
extern "C" bool fortran_stream_ele_end(
    int &physical_end /* 0D_NOT_integer in */,
    int &ele_orientation /* 0D_NOT_integer in */,
    int &stream_end /* 0D_NOT_integer out */
);
int stream_ele_end(int physical_end, int ele_orientation);
extern "C" void fortran_string_attrib(
    const char *attrib_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    const char *attrib_value /* 0D_NOT_character out */
);
std::string string_attrib(std::string attrib_name, EleStruct &ele);
extern "C" void fortran_strong_beam_sigma_calc(
    void *ele /* 0D_NOT_type in */,
    double &s_pos /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &sigma /* 1D_NOT_real out */,
    double &bbi_const /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &dsigma_ds /* 1D_NOT_real out */
);
struct StrongBeamSigmaCalc {
  FixedArray1D<Real, 2> sigma;
  double bbi_const;
  FixedArray1D<Real, 2> dsigma_ds;
};
Bmad::StrongBeamSigmaCalc strong_beam_sigma_calc(EleStruct &ele, double s_pos);
extern "C" bool fortran_strong_beam_strength(
    void *ele /* 0D_NOT_type in */,
    double &strength /* 0D_NOT_real out */
);
double strong_beam_strength(EleStruct &ele);
extern "C" void fortran_surface_grid_displacement(
    void *ele /* 0D_NOT_type in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &z /* 0D_NOT_real out */,
    Bmad::array_descriptor_t &dz_dxy /* 1D_NOT_real out */,
    bool *extend_grid /* 0D_NOT_logical in */
);
struct SurfaceGridDisplacement {
  bool err_flag;
  double z;
  FixedArray1D<Real, 2> dz_dxy;
};
Bmad::SurfaceGridDisplacement surface_grid_displacement(
    EleStruct &ele,
    double x,
    double y,
    std::optional<bool> extend_grid = std::nullopt
);

// Skipped unusable routine switch_attrib_value_name:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_symp_lie_bmad(
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    void *track /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */,
    bool *offset_ele /* 0D_NOT_logical in */
);
TrackStruct symp_lie_bmad(
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> offset_ele = std::nullopt
);
extern "C" void fortran_t6_to_b123(
    Bmad::array_descriptor_t &t6 /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &abz_tunes /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &B1 /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &B2 /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &B3 /* 2D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct T6ToB123 {
  FixedArray2D<Real, 6, 6> B1;
  FixedArray2D<Real, 6, 6> B2;
  FixedArray2D<Real, 6, 6> B3;
  bool err_flag;
};
Bmad::T6ToB123 t6_to_b123(FixedArray2D<Real, 6, 6> t6, FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_taper_mag_strengths(
    void *lat /* 0D_NOT_type inout */,
    void *ref_lat /* 0D_NOT_type in */,
    const char *except /* 0D_NOT_character in */,
    bool *err_flag /* 0D_NOT_logical in */
);
void taper_mag_strengths(
    LatStruct &lat,
    optional_ref<LatStruct> ref_lat = std::nullopt,
    std::optional<std::string> except = std::nullopt,
    std::optional<bool> err_flag = std::nullopt
);
extern "C" void fortran_target_min_max_calc(
    Bmad::array_descriptor_t &r_corner1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &r_corner2 /* 1D_NOT_real in */,
    double &y_min /* 0D_NOT_real inout */,
    double &y_max /* 0D_NOT_real inout */,
    double &phi_min /* 0D_NOT_real inout */,
    double &phi_max /* 0D_NOT_real inout */,
    bool *initial /* 0D_NOT_logical in */
);
void target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double &y_min,
    double &y_max,
    double &phi_min,
    double &phi_max,
    std::optional<bool> initial = std::nullopt
);
extern "C" void fortran_target_rot_mats(
    Bmad::array_descriptor_t &r_center /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &w_to_target /* 2D_NOT_real out */,
    Bmad::array_descriptor_t &w_to_ele /* 2D_NOT_real out */
);
struct TargetRotMats {
  FixedArray2D<Real, 3, 3> w_to_target;
  FixedArray2D<Real, 3, 3> w_to_ele;
};
Bmad::TargetRotMats target_rot_mats(FixedArray1D<Real, 3> r_center);
extern "C" void fortran_taylor_equal_taylor(
    void *taylor1 /* 0D_NOT_type inout */,
    void *taylor2 /* 0D_NOT_type in */
);
void taylor_equal_taylor(TaylorStruct &taylor1, TaylorStruct &taylor2);
extern "C" void fortran_taylor_inverse(
    Bmad::array_descriptor_t &taylor_in /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &taylor_inv /* 1D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */
);
bool taylor_inverse(TaylorStructArray1D taylor_in, TaylorStructArray1D taylor_inv);

// Skipped unusable routine taylor_minus_taylor:
// - Routine in configuration skip list

// Skipped unusable routine taylor_plus_taylor:
// - Routine in configuration skip list
extern "C" void fortran_taylor_propagate1(
    Bmad::array_descriptor_t &orb_taylor /* 1D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *ref_in /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &spin_taylor /* 1D_NOT_type inout */
);
bool taylor_propagate1(
    TaylorStructArray1D orb_taylor,
    EleStruct &ele,
    LatParamStruct &param,
    optional_ref<CoordStruct> ref_in = std::nullopt,
    std::optional<TaylorStructArray1D> spin_taylor = std::nullopt
);

// Skipped unusable routine taylor_to_genfield:
// - Untranslated type: genfield (0D)
extern "C" void fortran_taylor_to_mad_map(
    Bmad::array_descriptor_t &taylor /* 1D_NOT_type in */,
    void *energy /* 0D_NOT_type in */,
    void *map /* 0D_NOT_type out */
);
MadMapStruct taylor_to_mad_map(TaylorStructArray1D taylor, MadEnergyStruct &energy);

// Skipped unusable routine taylor_to_real_8:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_taylors_equal_taylors(
    Bmad::array_descriptor_t &taylor1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &taylor2 /* 1D_NOT_type in */
);
void taylors_equal_taylors(TaylorStructArray1D taylor1, TaylorStructArray1D taylor2);
extern "C" void fortran_tilt_coords(
    double &tilt_val /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &coord /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void tilt_coords(
    double tilt_val,
    FArray1D<Real> &coord,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_tilt_coords_photon(
    double &tilt_val /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &coord /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real inout */
);
void tilt_coords_photon(
    double tilt_val,
    FArray1D<Real> &coord,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat = std::nullopt
);
extern "C" void fortran_tilt_mat6(
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    double &tilt /* 0D_NOT_real in */
);
void tilt_mat6(FixedArray2D<Real, 6, 6> mat6, double tilt);

// Skipped unusable routine time_runge_kutta_periodic_kick_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_to_eta_reading(
    Bmad::array_descriptor_t &eta_actual /* 1D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    int &axis /* 0D_NOT_integer in */,
    bool &add_noise /* 0D_NOT_logical in */,
    double &reading /* 0D_NOT_real out */,
    bool &err /* 0D_NOT_logical out */
);
struct ToEtaReading {
  double reading;
  bool err;
};
Bmad::ToEtaReading
to_eta_reading(FArray1D<Real> &eta_actual, EleStruct &ele, int axis, bool add_noise);
extern "C" void fortran_to_fieldmap_coords(
    void *ele /* 0D_NOT_type in */,
    void *local_orb /* 0D_NOT_type in */,
    double &s_body /* 0D_NOT_real in */,
    int &ele_anchor_pt /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &r0 /* 1D_NOT_real in */,
    bool &curved_ref_frame /* 0D_NOT_logical in */,
    double &x /* 0D_NOT_real out */,
    double &y /* 0D_NOT_real out */,
    double &z /* 0D_NOT_real out */,
    double &cos_ang /* 0D_NOT_real out */,
    double &sin_ang /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct ToFieldmapCoords {
  double x;
  double y;
  double z;
  double cos_ang;
  double sin_ang;
  bool err_flag;
};
Bmad::ToFieldmapCoords to_fieldmap_coords(
    EleStruct &ele,
    CoordStruct &local_orb,
    double s_body,
    int ele_anchor_pt,
    FixedArray1D<Real, 3> r0,
    bool curved_ref_frame
);
extern "C" void fortran_to_orbit_reading(
    void *orb /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    int &axis /* 0D_NOT_integer in */,
    bool &add_noise /* 0D_NOT_logical in */,
    double &reading /* 0D_NOT_real out */,
    bool &err /* 0D_NOT_logical out */
);
struct ToOrbitReading {
  double reading;
  bool err;
};
Bmad::ToOrbitReading to_orbit_reading(CoordStruct &orb, EleStruct &ele, int axis, bool add_noise);
extern "C" void fortran_to_phase_and_coupling_reading(
    void *ele /* 0D_NOT_type in */,
    bool &add_noise /* 0D_NOT_logical in */,
    void *reading /* 0D_NOT_type out */,
    bool &err /* 0D_NOT_logical out */
);
struct ToPhaseAndCouplingReading {
  BpmPhaseCouplingStruct reading;
  bool err;
};
Bmad::ToPhaseAndCouplingReading to_phase_and_coupling_reading(EleStruct &ele, bool add_noise);
extern "C" bool fortran_to_photon_angle_coords(
    void *orb_in /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *orb_out /* 0D_NOT_type out */
);
CoordStruct to_photon_angle_coords(CoordStruct &orb_in, EleStruct &ele);
extern "C" void fortran_to_surface_coords(
    void *lab_orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *surface_orbit /* 0D_NOT_type out */
);
CoordStruct to_surface_coords(CoordStruct &lab_orbit, EleStruct &ele);
extern "C" void fortran_touschek_lifetime(
    void *mode /* 0D_NOT_type in */,
    double &Tl /* 0D_NOT_real out */,
    void *lat /* 0D_NOT_type in */
);
double touschek_lifetime(NormalModesStruct &mode, LatStruct &lat);

// Skipped unusable routine touschek_lifetime_ele_by_ele:
// - Untranslated type: momentum_aperture_struct (1D)

// Skipped unusable routine touschek_lifetime_with_aperture:
// - Untranslated type: momentum_aperture_struct (1D)
extern "C" void fortran_touschek_rate1(
    void *mode /* 0D_NOT_type in */,
    double &rate /* 0D_NOT_real out */,
    void *lat /* 0D_NOT_type in */,
    int *ix /* 0D_NOT_integer in */,
    double *s /* 0D_NOT_real in */
);
double touschek_rate1(
    NormalModesStruct &mode,
    LatStruct &lat,
    std::optional<int> ix = std::nullopt,
    std::optional<double> s = std::nullopt
);
extern "C" void fortran_touschek_rate1_zap(
    void *mode /* 0D_NOT_type inout */,
    double &rate /* 0D_NOT_real in */,
    void *lat /* 0D_NOT_type inout */,
    int *ix /* 0D_NOT_integer in */,
    double *s /* 0D_NOT_real in */
);
void touschek_rate1_zap(
    NormalModesStruct &mode,
    double rate,
    LatStruct &lat,
    std::optional<int> ix = std::nullopt,
    std::optional<double> s = std::nullopt
);
extern "C" void fortran_track1(
    void *start_orb /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    void *track /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *ignore_radiation /* 0D_NOT_logical in */,
    bool *make_map1 /* 0D_NOT_logical in */,
    bool *init_to_edge /* 0D_NOT_logical in */
);
struct Track1 {
  CoordStruct end_orb;
  bool err_flag;
};
Bmad::Track1 track1(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    optional_ref<TrackStruct> track = std::nullopt,
    std::optional<bool> ignore_radiation = std::nullopt,
    std::optional<bool> make_map1 = std::nullopt,
    std::optional<bool> init_to_edge = std::nullopt
);
extern "C" void fortran_track1_beam(
    void *beam /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    int *direction /* 0D_NOT_integer in */
);
bool track1_beam(
    BeamStruct &beam,
    EleStruct &ele,
    std::optional<CoordStructArray1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt
);
extern "C" void fortran_track1_bmad(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
struct Track1Bmad {
  bool err_flag;
  TrackStruct track;
};
Bmad::Track1Bmad track1_bmad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track1_bmad_photon(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool track1_bmad_photon(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_track1_bunch(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    int *direction /* 0D_NOT_integer in */,
    void *bunch_track /* 0D_NOT_type inout */
);
bool track1_bunch(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<CoordStructArray1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);
extern "C" void fortran_track1_bunch_csr(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    double *s_start /* 0D_NOT_real in */,
    double *s_end /* 0D_NOT_real in */,
    void *bunch_track /* 0D_NOT_type inout */
);
bool track1_bunch_csr(
    BunchStruct &bunch,
    EleStruct &ele,
    CoordStructArray1D centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);
extern "C" void fortran_track1_bunch_csr3d(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    double *s_start /* 0D_NOT_real in */,
    double *s_end /* 0D_NOT_real in */,
    void *bunch_track /* 0D_NOT_type inout */
);
bool track1_bunch_csr3d(
    BunchStruct &bunch,
    EleStruct &ele,
    CoordStructArray1D centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);
extern "C" void fortran_track1_bunch_hom(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    int *direction /* 0D_NOT_integer in */,
    void *bunch_track /* 0D_NOT_type inout */
);
void track1_bunch_hom(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);

// Skipped unusable routine track1_bunch_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_bunch_space_charge(
    void *bunch /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    bool *track_to_same_s /* 0D_NOT_logical in */,
    void *bunch_track /* 0D_NOT_type inout */
);
bool track1_bunch_space_charge(
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<bool> track_to_same_s = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);
extern "C" void fortran_track1_crystal(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_crystal(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);

// Skipped unusable routine track1_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_diffraction_plate_or_mask(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_diffraction_plate_or_mask(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_high_energy_space_charge(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_high_energy_space_charge(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_lens(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_lens(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_linear(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type inout */
);
void track1_linear(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param);
extern "C" void
fortran_track1_lr_wake(void *bunch /* 0D_NOT_type inout */, void *ele /* 0D_NOT_type inout */);
void track1_lr_wake(BunchStruct &bunch, EleStruct &ele);
extern "C" void fortran_track1_mad(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */
);
void track1_mad(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_track1_mirror(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_mirror(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_mosaic_crystal(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_mosaic_crystal(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_multilayer_mirror(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_multilayer_mirror(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);

// Skipped unusable routine track1_postprocess_def:
// - Routine in configuration skip list

// Skipped unusable routine track1_preprocess_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_radiation(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    int &edge /* 0D_NOT_integer in */
);
void track1_radiation(CoordStruct &orbit, EleStruct &ele, int edge);
extern "C" void fortran_track1_radiation_center(
    void *orbit /* 0D_NOT_type inout */,
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool *rad_damp /* 0D_NOT_logical in */,
    bool *rad_fluct /* 0D_NOT_logical in */
);
void track1_radiation_center(
    CoordStruct &orbit,
    EleStruct &ele1,
    EleStruct &ele2,
    std::optional<bool> rad_damp = std::nullopt,
    std::optional<bool> rad_fluct = std::nullopt
);
extern "C" void fortran_track1_runge_kutta(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
struct Track1RungeKutta {
  bool err_flag;
  TrackStruct track;
};
Bmad::Track1RungeKutta track1_runge_kutta(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track1_sample(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */
);
void track1_sample(EleStruct &ele, LatParamStruct &param, CoordStruct &orbit);
extern "C" void fortran_track1_spin(
    void *start_orb /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type inout */,
    bool *make_quaternion /* 0D_NOT_logical in */
);
void track1_spin(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &end_orb,
    std::optional<bool> make_quaternion = std::nullopt
);

// Skipped unusable routine track1_spin_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_spin_integration(
    void *start_orb /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type inout */
);
void track1_spin_integration(
    CoordStruct &start_orb,
    EleStruct &ele,
    LatParamStruct &param,
    CoordStruct &end_orb
);
extern "C" void fortran_track1_spin_taylor(
    void *start_orb /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */
);
CoordStruct track1_spin_taylor(CoordStruct &start_orb, EleStruct &ele, LatParamStruct &param);
extern "C" void
fortran_track1_sr_wake(void *bunch /* 0D_NOT_type inout */, void *ele /* 0D_NOT_type in */);
void track1_sr_wake(BunchStruct &bunch, EleStruct &ele);
extern "C" void fortran_track1_symp_lie_ptc(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *track /* 0D_NOT_type out */
);
TrackStruct track1_symp_lie_ptc(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_track1_taylor(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &taylor /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track1_taylor(
    CoordStruct &orbit,
    EleStruct &ele,
    std::optional<TaylorStructArray1D> taylor = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track1_time_runge_kutta(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type out */,
    double *t_end /* 0D_NOT_real in */,
    double *dt_step /* 0D_NOT_real inout */
);
struct Track1TimeRungeKutta {
  bool err_flag;
  TrackStruct track;
};
Bmad::Track1TimeRungeKutta track1_time_runge_kutta(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<double> t_end = std::nullopt,
    optional_ref<double> dt_step = std::nullopt
);

// Skipped unusable routine track1_wake_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track_a_beambeam(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    void *track /* 0D_NOT_type out */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
struct TrackABeambeam {
  TrackStruct track;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
Bmad::TrackABeambeam track_a_beambeam(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_bend(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_bend(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_bend_photon(
    void *orb /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    double &length /* 0D_NOT_real in */
);
void track_a_bend_photon(CoordStruct &orb, EleStruct &ele, double length);
extern "C" void
fortran_track_a_capillary(void *orb /* 0D_NOT_type inout */, void *ele /* 0D_NOT_type in */);
void track_a_capillary(CoordStruct &orb, EleStruct &ele);
extern "C" void fortran_track_a_converter(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_converter(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_crab_cavity(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_crab_cavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_drift(
    void *orb /* 0D_NOT_type inout */,
    double &length /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */,
    int *ele_orientation /* 0D_NOT_integer in */,
    bool *include_ref_motion /* 0D_NOT_logical in */,
    double *time /* 0D_NOT_real inout */
);
void track_a_drift(
    CoordStruct &orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<int> ele_orientation = std::nullopt,
    std::optional<bool> include_ref_motion = std::nullopt,
    optional_ref<double> time = std::nullopt
);
extern "C" void fortran_track_a_drift_photon(
    void *orb /* 0D_NOT_type inout */,
    double &length /* 0D_NOT_real in */,
    bool &phase_relative_to_ref /* 0D_NOT_logical in */
);
void track_a_drift_photon(CoordStruct &orb, double length, bool phase_relative_to_ref);
extern "C" void fortran_track_a_foil(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_foil(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_gkicker(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_gkicker(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_lcavity(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_lcavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_lcavity_old(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_lcavity_old(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_mask(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_mask(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_match(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool *err_flag /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_match(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_patch(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    bool *drift_to_exit /* 0D_NOT_logical in */,
    double &s_ent /* 0D_NOT_real out */,
    double &ds_ref /* 0D_NOT_real out */,
    bool *track_spin /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
struct TrackAPatch {
  double s_ent;
  double ds_ref;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
Bmad::TrackAPatch track_a_patch(
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> track_spin = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_patch_photon(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    bool *drift_to_exit /* 0D_NOT_logical in */,
    bool *use_z_pos /* 0D_NOT_logical in */
);
void track_a_patch_photon(
    EleStruct &ele,
    CoordStruct &orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> use_z_pos = std::nullopt
);
extern "C" void fortran_track_a_pickup(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool *err_flag /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_pickup(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_quadrupole(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_quadrupole(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_rfcavity(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_rfcavity(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_sad_mult(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_sad_mult(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_sol_quad(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_sol_quad(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_thick_multipole(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real inout */,
    bool *make_matrix /* 0D_NOT_logical in */
);
void track_a_thick_multipole(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_wiggler(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat6 /* 2D_NOT_real out */,
    bool *make_matrix /* 0D_NOT_logical in */
);
std::optional<FixedArray2D<Real, 6, 6>> track_a_wiggler(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<bool> make_matrix = std::nullopt
);
extern "C" void fortran_track_a_zero_length_element(
    void *orbit /* 0D_NOT_type inout */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *track /* 0D_NOT_type out */
);
struct TrackAZeroLengthElement {
  bool err_flag;
  TrackStruct track;
};
Bmad::TrackAZeroLengthElement
track_a_zero_length_element(CoordStruct &orbit, EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_track_all(
    void *lat /* 0D_NOT_type in */,
    void *orbit /* 1D_ALLOC_type inout */,
    int *ix_branch /* 0D_NOT_integer in */,
    int &track_state /* 0D_NOT_integer out */,
    bool &err_flag /* 0D_NOT_logical out */,
    void *orbit0 /* 1D_ALLOC_type out */,
    bool *init_lost /* 0D_NOT_logical in */
);
struct TrackAll {
  int track_state;
  bool err_flag;
  CoordStructAlloc1D orbit0;
};
Bmad::TrackAll track_all(
    LatStruct &lat,
    CoordStructAlloc1D orbit,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> init_lost = std::nullopt
);
extern "C" void fortran_track_beam(
    void *lat /* 0D_NOT_type in */,
    void *beam /* 0D_NOT_type inout */,
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    int *direction /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &bunch_tracks /* 1D_NOT_type inout */
);
bool track_beam(
    LatStruct &lat,
    BeamStruct &beam,
    optional_ref<EleStruct> ele1 = std::nullopt,
    optional_ref<EleStruct> ele2 = std::nullopt,
    std::optional<CoordStructArray1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<BunchTrackStructArray1D> bunch_tracks = std::nullopt
);
extern "C" void fortran_track_bunch(
    void *lat /* 0D_NOT_type in */,
    void *bunch /* 0D_NOT_type inout */,
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &centroid /* 1D_NOT_type in */,
    int *direction /* 0D_NOT_integer in */,
    void *bunch_track /* 0D_NOT_type inout */
);
bool track_bunch(
    LatStruct &lat,
    BunchStruct &bunch,
    optional_ref<EleStruct> ele1 = std::nullopt,
    optional_ref<EleStruct> ele2 = std::nullopt,
    std::optional<CoordStructArray1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackStruct> bunch_track = std::nullopt
);
extern "C" void fortran_track_bunch_time(
    void *bunch /* 0D_NOT_type inout */,
    void *branch /* 0D_NOT_type in */,
    double &t_end /* 0D_NOT_real in */,
    double &s_end /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &dt_step /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &extra_field /* 1D_NOT_type in */
);
void track_bunch_time(
    BunchStruct &bunch,
    BranchStruct &branch,
    double t_end,
    double s_end,
    std::optional<FArray1D<Real>> dt_step = std::nullopt,
    std::optional<EmFieldStructArray1D> extra_field = std::nullopt
);
extern "C" void fortran_track_bunch_to_s(
    void *bunch /* 0D_NOT_type inout */,
    double &s /* 0D_NOT_real in */,
    void *branch /* 0D_NOT_type in */
);
void track_bunch_to_s(BunchStruct &bunch, double s, BranchStruct &branch);
extern "C" void fortran_track_bunch_to_t(
    void *bunch /* 0D_NOT_type inout */,
    double &t_target /* 0D_NOT_real in */,
    void *branch /* 0D_NOT_type in */
);
void track_bunch_to_t(BunchStruct &bunch, double t_target, BranchStruct &branch);
extern "C" void fortran_track_complex_taylor(
    Bmad::array_descriptor_t &start_orb /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &complex_taylor /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &end_orb /* 1D_NOT_complex inout */
);
void track_complex_taylor(
    FArray1D<Complex> &start_orb,
    ComplexTaylorStructArray1D complex_taylor,
    FArray1D<Complex> &end_orb
);
extern "C" void fortran_track_from_s_to_s(
    void *lat /* 0D_NOT_type in */,
    double &s_start /* 0D_NOT_real in */,
    double &s_end /* 0D_NOT_real in */,
    void *orbit_start /* 0D_NOT_type in */,
    void *orbit_end /* 0D_NOT_type out */,
    void *all_orb /* 1D_ALLOC_type out */,
    int *ix_branch /* 0D_NOT_integer in */,
    int &track_state /* 0D_NOT_integer out */,
    int *ix_ele_end /* 0D_NOT_integer in */
);
struct TrackFromSToS {
  CoordStruct orbit_end;
  CoordStructAlloc1D all_orb;
  int track_state;
};
Bmad::TrackFromSToS track_from_s_to_s(
    LatStruct &lat,
    double s_start,
    double s_end,
    CoordStruct &orbit_start,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ix_ele_end = std::nullopt
);
extern "C" void fortran_track_many(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &orbit /* 1D_NOT_type inout */,
    int &ix_start /* 0D_NOT_integer in */,
    int &ix_end /* 0D_NOT_integer in */,
    int &direction /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    int &track_state /* 0D_NOT_integer out */
);
int track_many(
    LatStruct &lat,
    CoordStructArray1D orbit,
    int ix_start,
    int ix_end,
    int direction,
    std::optional<int> ix_branch = std::nullopt
);

// Skipped unusable routine track_many_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track_to_surface(
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type inout */,
    void *param /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &w_surface /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3>
track_to_surface(EleStruct &ele, CoordStruct &orbit, LatParamStruct &param);
extern "C" void fortran_track_until_dead(
    void *start_orb /* 0D_NOT_type in */,
    void *lat /* 0D_NOT_type in */,
    void *end_orb /* 0D_NOT_type out */,
    void *track /* 0D_NOT_type out */
);
struct TrackUntilDead {
  CoordStruct end_orb;
  TrackStruct track;
};
Bmad::TrackUntilDead track_until_dead(CoordStruct &start_orb, LatStruct &lat);
extern "C" void fortran_tracking_rad_map_setup(
    void *ele /* 0D_NOT_type in */,
    double &tollerance /* 0D_NOT_real in */,
    int &ref_edge /* 0D_NOT_integer in */,
    void *rad_map /* 0D_NOT_type out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct TrackingRadMapSetup {
  RadMapStruct rad_map;
  bool err_flag;
};
Bmad::TrackingRadMapSetup tracking_rad_map_setup(EleStruct &ele, double tollerance, int ref_edge);
extern "C" void
fortran_transfer_ac_kick(void *ac_in /* 0D_PTR_type in */, void *ac_out /* 0D_PTR_type out */);
std::optional<AcKickerStruct> transfer_ac_kick(AcKickerStruct &ac_in);
extern "C" void
fortran_transfer_branch(void *branch1 /* 0D_NOT_type in */, void *branch2 /* 0D_NOT_type out */);
BranchStruct transfer_branch(BranchStruct &branch1);
extern "C" void fortran_transfer_branch_parameters(
    void *branch_in /* 0D_NOT_type in */,
    void *branch_out /* 0D_NOT_type out */
);
BranchStruct transfer_branch_parameters(BranchStruct &branch_in);
extern "C" void fortran_transfer_branches(
    Bmad::array_descriptor_t &branch1 /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &branch2 /* 1D_NOT_type inout */
);
void transfer_branches(BranchStructArray1D branch1, BranchStructArray1D branch2);
extern "C" void fortran_transfer_ele(
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type out */,
    bool *nullify_pointers /* 0D_NOT_logical in */
);
EleStruct transfer_ele(EleStruct &ele1, std::optional<bool> nullify_pointers = std::nullopt);
extern "C" void fortran_transfer_ele_taylor(
    void *ele_in /* 0D_NOT_type in */,
    void *ele_out /* 0D_NOT_type out */,
    int *taylor_order /* 0D_NOT_integer in */
);
EleStruct transfer_ele_taylor(EleStruct &ele_in, std::optional<int> taylor_order = std::nullopt);
extern "C" void fortran_transfer_eles(
    Bmad::array_descriptor_t &ele1 /* 1D_NOT_type in */,
    Bmad::array_descriptor_t &ele2 /* 1D_NOT_type inout */
);
void transfer_eles(EleStructArray1D ele1, EleStructArray1D ele2);
extern "C" void fortran_transfer_fieldmap(
    void *ele_in /* 0D_NOT_type in */,
    void *ele_out /* 0D_NOT_type out */,
    int &who /* 0D_NOT_integer in */
);
EleStruct transfer_fieldmap(EleStruct &ele_in, int who);
extern "C" bool fortran_transfer_fixer_params(
    void *fixer /* 0D_NOT_type in */,
    bool &to_stored /* 0D_NOT_logical in */,
    void *orbit /* 0D_NOT_type inout */,
    const char *who /* 0D_NOT_character in */,
    bool &is_ok /* 0D_NOT_logical out */
);
bool transfer_fixer_params(
    EleStruct &fixer,
    bool to_stored,
    optional_ref<CoordStruct> orbit = std::nullopt,
    std::optional<std::string> who = std::nullopt
);
extern "C" void
fortran_transfer_lat(void *lat1 /* 0D_NOT_type in */, void *lat2 /* 0D_NOT_type out */);
LatStruct transfer_lat(LatStruct &lat1);
extern "C" void fortran_transfer_lat_parameters(
    void *lat_in /* 0D_NOT_type in */,
    void *lat_out /* 0D_NOT_type out */
);
LatStruct transfer_lat_parameters(LatStruct &lat_in);
extern "C" void fortran_transfer_map_calc(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &orb_map /* 1D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    int *ix1 /* 0D_NOT_integer in */,
    int *ix2 /* 0D_NOT_integer in */,
    void *ref_orb /* 0D_NOT_type in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *one_turn /* 0D_NOT_logical in */,
    bool *unit_start /* 0D_NOT_logical in */,
    bool *concat_if_possible /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin_map /* 1D_NOT_type inout */
);
bool transfer_map_calc(
    LatStruct &lat,
    TaylorStructArray1D orb_map,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    optional_ref<CoordStruct> ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt,
    std::optional<bool> concat_if_possible = std::nullopt,
    std::optional<TaylorStructArray1D> spin_map = std::nullopt
);
extern "C" void fortran_transfer_map_from_s_to_s(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &t_map /* 1D_NOT_type inout */,
    double *s1 /* 0D_NOT_real in */,
    double *s2 /* 0D_NOT_real in */,
    void *ref_orb_in /* 0D_NOT_type in */,
    void *ref_orb_out /* 0D_NOT_type out */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *one_turn /* 0D_NOT_logical in */,
    bool *unit_start /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *concat_if_possible /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &spin_map /* 1D_NOT_type inout */
);
struct TransferMapFromSToS {
  CoordStruct ref_orb_out;
  bool err_flag;
};
Bmad::TransferMapFromSToS transfer_map_from_s_to_s(
    LatStruct &lat,
    TaylorStructArray1D t_map,
    std::optional<double> s1 = std::nullopt,
    std::optional<double> s2 = std::nullopt,
    optional_ref<CoordStruct> ref_orb_in = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt,
    std::optional<bool> concat_if_possible = std::nullopt,
    std::optional<TaylorStructArray1D> spin_map = std::nullopt
);
extern "C" void fortran_transfer_mat2_from_twiss(
    void *twiss1 /* 0D_NOT_type in */,
    void *twiss2 /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 2, 2> transfer_mat2_from_twiss(TwissStruct &twiss1, TwissStruct &twiss2);
extern "C" void fortran_transfer_mat_from_twiss(
    void *ele1 /* 0D_NOT_type in */,
    void *ele2 /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &orb1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &orb2 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &m /* 2D_NOT_real out */
);
FixedArray2D<Real, 6, 6> transfer_mat_from_twiss(
    EleStruct &ele1,
    EleStruct &ele2,
    FixedArray1D<Real, 6> orb1,
    FixedArray1D<Real, 6> orb2
);
extern "C" void fortran_transfer_matrix_calc(
    void *lat /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &xfer_mat /* 2D_NOT_real inout */,
    Bmad::array_descriptor_t &xfer_vec /* 1D_NOT_real inout */,
    int *ix1 /* 0D_NOT_integer in */,
    int *ix2 /* 0D_NOT_integer in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *one_turn /* 0D_NOT_logical in */
);
void transfer_matrix_calc(
    LatStruct &lat,
    FixedArray2D<Real, 6, 6> xfer_mat,
    std::optional<FixedArray1D<Real, 6>> xfer_vec = std::nullopt,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt
);
extern "C" void fortran_transfer_twiss(
    void *ele_in /* 0D_NOT_type in */,
    void *ele_out /* 0D_NOT_type out */,
    bool *reverse /* 0D_NOT_logical in */
);
EleStruct transfer_twiss(EleStruct &ele_in, std::optional<bool> reverse = std::nullopt);
extern "C" void
fortran_transfer_wake(void *wake_in /* 0D_PTR_type in */, void *wake_out /* 0D_PTR_type out */);
std::optional<WakeStruct> transfer_wake(WakeStruct &wake_in);

// Skipped unusable routine transfer_wall3d:
// - Routine in configuration skip list

// Skipped unusable routine transport_with_sr:
// - Variable in sized array: 2D_NOT_real
// - Variable in sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine transport_with_sr_and_ibs:
// - Variable in sized array: 2D_NOT_real
// - Variable in sized array: 2D_NOT_real
// - Variable in sized array: 2D_NOT_real
extern "C" void fortran_truncate_complex_taylor_to_order(
    Bmad::array_descriptor_t &complex_taylor_in /* 1D_NOT_type in */,
    int &order /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &complex_taylor_out /* 1D_NOT_type inout */
);
void truncate_complex_taylor_to_order(
    ComplexTaylorStructArray1D complex_taylor_in,
    int order,
    ComplexTaylorStructArray1D complex_taylor_out
);
extern "C" void fortran_twiss1_propagate(
    void *twiss1 /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &mat2 /* 2D_NOT_real in */,
    int &ele_key /* 0D_NOT_integer in */,
    double &length /* 0D_NOT_real in */,
    void *twiss2 /* 0D_NOT_type out */,
    bool &err /* 0D_NOT_logical out */
);
struct Twiss1Propagate {
  TwissStruct twiss2;
  bool err;
};
Bmad::Twiss1Propagate
twiss1_propagate(TwissStruct &twiss1, FixedArray2D<Real, 2, 2> mat2, int ele_key, double length);
extern "C" void fortran_twiss3_at_start(
    void *lat /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */,
    int *ix_branch /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &tune3 /* 1D_NOT_real out */
);
FixedArray1D<Real, 3>
twiss3_at_start(LatStruct &lat, bool err_flag, std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss3_from_twiss2(void *ele /* 0D_NOT_type inout */);
void twiss3_from_twiss2(EleStruct &ele);
extern "C" void fortran_twiss3_propagate1(
    void *ele1 /* 0D_NOT_type inout */,
    void *ele2 /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */
);
void twiss3_propagate1(EleStruct &ele1, EleStruct &ele2, bool err_flag);
extern "C" void fortran_twiss3_propagate_all(
    void *lat /* 0D_NOT_type in */,
    int *ix_branch /* 0D_NOT_integer in */
);
void twiss3_propagate_all(LatStruct &lat, std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss_and_track_all(
    void *lat /* 0D_NOT_type inout */,
    void *orb_array /* 1D_ALLOC_type inout */,
    int &status /* 0D_NOT_integer out */,
    bool *print_err /* 0D_NOT_logical in */,
    bool *calc_chrom /* 0D_NOT_logical in */
);
int twiss_and_track(
    LatStruct &lat,
    CoordArrayStructAlloc1D orb_array,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> calc_chrom = std::nullopt
);
extern "C" void fortran_twiss_and_track_at_s(
    void *lat /* 0D_NOT_type in */,
    double &s /* 0D_NOT_real in */,
    void *ele_at_s /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &orb /* 1D_NOT_type in */,
    void *orb_at_s /* 0D_NOT_type inout */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical out */,
    bool *use_last /* 0D_NOT_logical in */,
    bool *compute_floor_coords /* 0D_NOT_logical in */
);
bool twiss_and_track_at_s(
    LatStruct &lat,
    double s,
    optional_ref<EleStruct> ele_at_s = std::nullopt,
    std::optional<CoordStructArray1D> orb = std::nullopt,
    optional_ref<CoordStruct> orb_at_s = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> use_last = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt
);
extern "C" void fortran_twiss_and_track_branch(
    void *lat /* 0D_NOT_type inout */,
    void *orb /* 1D_ALLOC_type inout */,
    int &status /* 0D_NOT_integer out */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */,
    bool *calc_chrom /* 0D_NOT_logical in */,
    void *orb_start /* 0D_NOT_type in */
);
int twiss_and_track(
    LatStruct &lat,
    CoordStructAlloc1D orb,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> calc_chrom = std::nullopt,
    optional_ref<CoordStruct> orb_start = std::nullopt
);
extern "C" void fortran_twiss_and_track_from_s_to_s(
    void *branch /* 0D_NOT_type in */,
    void *orbit_start /* 0D_NOT_type in */,
    double &s_end /* 0D_NOT_real in */,
    void *orbit_end /* 0D_NOT_type out */,
    void *ele_start /* 0D_NOT_type in */,
    void *ele_end /* 0D_NOT_type out */,
    bool &err /* 0D_NOT_logical out */,
    bool *compute_floor_coords /* 0D_NOT_logical in */,
    bool *compute_twiss /* 0D_NOT_logical in */
);
struct TwissAndTrackFromSToS {
  CoordStruct orbit_end;
  EleStruct ele_end;
  bool err;
};
Bmad::TwissAndTrackFromSToS twiss_and_track_from_s_to_s(
    BranchStruct &branch,
    CoordStruct &orbit_start,
    double s_end,
    optional_ref<EleStruct> ele_start = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt,
    std::optional<bool> compute_twiss = std::nullopt
);
extern "C" void fortran_twiss_and_track_intra_ele(
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &l_start /* 0D_NOT_real in */,
    double &l_end /* 0D_NOT_real in */,
    bool &track_upstream_end /* 0D_NOT_logical in */,
    bool &track_downstream_end /* 0D_NOT_logical in */,
    void *orbit_start /* 0D_NOT_type in */,
    void *orbit_end /* 0D_NOT_type out */,
    void *ele_start /* 0D_NOT_type in */,
    void *ele_end /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */,
    bool *compute_floor_coords /* 0D_NOT_logical in */,
    bool *compute_twiss /* 0D_NOT_logical in */,
    bool *reuse_ele_end /* 0D_NOT_logical in */
);
struct TwissAndTrackIntraEle {
  CoordStruct orbit_end;
  bool err;
};
Bmad::TwissAndTrackIntraEle twiss_and_track_intra_ele(
    EleStruct &ele,
    LatParamStruct &param,
    double l_start,
    double l_end,
    bool track_upstream_end,
    bool track_downstream_end,
    optional_ref<CoordStruct> orbit_start = std::nullopt,
    optional_ref<EleStruct> ele_start = std::nullopt,
    optional_ref<EleStruct> ele_end = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt,
    std::optional<bool> compute_twiss = std::nullopt,
    std::optional<bool> reuse_ele_end = std::nullopt
);
extern "C" void fortran_twiss_at_element(
    void *ele /* 0D_NOT_type in */,
    void *start /* 0D_NOT_type out */,
    void *end /* 0D_NOT_type out */,
    void *average /* 0D_NOT_type out */
);
struct TwissAtElement {
  EleStruct start;
  EleStruct end;
  EleStruct average;
};
Bmad::TwissAtElement twiss_at_element(EleStruct &ele);
extern "C" void fortran_twiss_at_start(
    void *lat /* 0D_NOT_type inout */,
    int &status /* 0D_NOT_integer out */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *type_out /* 0D_NOT_logical in */
);
int twiss_at_start(
    LatStruct &lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> type_out = std::nullopt
);

// Skipped unusable routine twiss_from_mat2:
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine twiss_from_mat6:
// - Variable in sized array: 2D_NOT_real
extern "C" void fortran_twiss_from_tracking(
    void *lat /* 0D_NOT_type inout */,
    void *ref_orb0 /* 0D_NOT_type in */,
    double &symp_err /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &d_orb /* 1D_NOT_real in */
);
struct TwissFromTracking {
  double symp_err;
  bool err_flag;
};
Bmad::TwissFromTracking twiss_from_tracking(
    LatStruct &lat,
    CoordStruct &ref_orb0,
    std::optional<FArray1D<Real>> d_orb = std::nullopt
);
extern "C" void fortran_twiss_propagate1(
    void *ele1 /* 0D_NOT_type inout */,
    void *ele2 /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *forward /* 0D_NOT_logical in */
);
bool twiss_propagate1(EleStruct &ele1, EleStruct &ele2, std::optional<bool> forward = std::nullopt);
extern "C" void fortran_twiss_propagate_all(
    void *lat /* 0D_NOT_type inout */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical out */,
    int *ie_start /* 0D_NOT_integer in */,
    int *ie_end /* 0D_NOT_integer in */
);
bool twiss_propagate_all(
    LatStruct &lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ie_start = std::nullopt,
    std::optional<int> ie_end = std::nullopt
);
extern "C" void fortran_twiss_to_1_turn_mat(
    void *twiss /* 0D_NOT_type in */,
    double &phi /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &mat2 /* 2D_NOT_real out */
);
FixedArray2D<Real, 2, 2> twiss_to_1_turn_mat(TwissStruct &twiss, double phi);

// Skipped unusable routine type_complex_taylors:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_type_coord(void *coord /* 0D_NOT_type in */);
void type_coord(CoordStruct &coord);

// Skipped unusable routine type_ele:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_end_stuff:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void
fortran_type_expression_tree(void *tree /* 0D_NOT_type in */, int *indent /* 0D_NOT_integer in */);
void type_expression_tree(ExpressionTreeStruct &tree, std::optional<int> indent = std::nullopt);

// Skipped unusable routine type_map:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_map1:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_ptc_fibre:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_ptc_internal_state:
// - Untranslated type: internal_state (0D)
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_type_ptc_layout(void *lay /* 0D_NOT_type inout */);
void type_ptc_layout(Layout &lay);

// Skipped unusable routine type_real_8_taylors:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_taylors:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_twiss:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine universal_equal_universal:
// - Untranslated type: universal_taylor (0D)
// - Untranslated type: universal_taylor (0D)

// Skipped unusable routine universal_to_bmad_taylor:
// - Untranslated type: universal_taylor (0D)

// Skipped unusable routine unlink_fieldmap:
// - Routine in configuration skip list

// Skipped unusable routine unlink_wall3d:
// - Routine in configuration skip list
extern "C" void fortran_update_ele_from_fibre(void *ele /* 0D_NOT_type inout */);
void update_ele_from_fibre(EleStruct &ele);
extern "C" void fortran_update_fibre_from_ele(
    void *ele /* 0D_NOT_type in */,
    bool &survey_needed /* 0D_NOT_logical out */
);
bool update_fibre_from_ele(EleStruct &ele);
extern "C" void
fortran_update_floor_angles(void *floor /* 0D_NOT_type inout */, void *floor0 /* 0D_NOT_type in */);
void update_floor_angles(
    FloorPositionStruct &floor,
    optional_ref<FloorPositionStruct> floor0 = std::nullopt
);
extern "C" bool fortran_valid_field_calc(
    void *ele /* 0D_NOT_type in */,
    int &field_calc /* 0D_NOT_integer in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool valid_field_calc(EleStruct &ele, int field_calc);
extern "C" bool fortran_valid_fringe_type(
    void *ele /* 0D_NOT_type in */,
    int &fringe_type /* 0D_NOT_integer in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool valid_fringe_type(EleStruct &ele, int fringe_type);
extern "C" bool fortran_valid_mat6_calc_method(
    void *ele /* 0D_NOT_type in */,
    int &species /* 0D_NOT_integer in */,
    int &mat6_calc_method /* 0D_NOT_integer in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool valid_mat6_calc_method(EleStruct &ele, int species, int mat6_calc_method);
extern "C" bool fortran_valid_spin_tracking_method(
    void *ele /* 0D_NOT_type in */,
    int &spin_tracking_method /* 0D_NOT_integer in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool valid_spin_tracking_method(EleStruct &ele, int spin_tracking_method);
extern "C" bool fortran_valid_tracking_method(
    void *ele /* 0D_NOT_type in */,
    int &species /* 0D_NOT_integer in */,
    int &tracking_method /* 0D_NOT_integer in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool valid_tracking_method(EleStruct &ele, int species, int tracking_method);
extern "C" bool fortran_value_of_attribute(
    void *ele /* 0D_NOT_type in */,
    const char *attrib_name /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    double *err_value /* 0D_NOT_real in */,
    double &value /* 0D_NOT_real out */
);
struct ValueOfAttribute {
  bool err_flag;
  double value;
};
Bmad::ValueOfAttribute value_of_attribute(
    EleStruct &ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<double> err_value = std::nullopt
);
extern "C" void fortran_value_to_line(
    const char *line /* 0D_NOT_character in */,
    double &value /* 0D_NOT_real in */,
    const char *str /* 0D_NOT_character in */,
    const char *typ /* 0D_NOT_character in */,
    bool *ignore_if_zero /* 0D_NOT_logical in */,
    bool *use_comma /* 0D_NOT_logical in */
);
void value_to_line(
    std::string line,
    double value,
    std::string str,
    std::string typ,
    std::optional<bool> ignore_if_zero = std::nullopt,
    std::optional<bool> use_comma = std::nullopt
);
extern "C" bool fortran_vec_to_polar(
    Bmad::array_descriptor_t &vec /* 1D_NOT_real in */,
    double *phase /* 0D_NOT_real in */,
    void *polar /* 0D_NOT_type out */
);
SpinPolarStruct vec_to_polar(FixedArray1D<Real, 3> vec, std::optional<double> phase = std::nullopt);
extern "C" bool fortran_vec_to_spinor(
    Bmad::array_descriptor_t &vec /* 1D_NOT_real in */,
    double *phase /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &spinor /* 1D_NOT_complex out */
);
FixedArray1D<Complex, 2>
vec_to_spinor(FixedArray1D<Real, 3> vec, std::optional<double> phase = std::nullopt);
extern "C" bool fortran_verify_valid_name(
    const char *name /* 0D_NOT_character in */,
    int &ix_name /* 0D_NOT_integer in */,
    bool *pure_name /* 0D_NOT_logical in */,
    bool *include_wild /* 0D_NOT_logical in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool verify_valid_name(
    std::string name,
    int ix_name,
    std::optional<bool> pure_name = std::nullopt,
    std::optional<bool> include_wild = std::nullopt
);
extern "C" bool fortran_w_mat_for_bend_angle(
    double &angle /* 0D_NOT_real in */,
    double &ref_tilt /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &r_vec /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3> w_mat_for_bend_angle(
    double angle,
    double ref_tilt,
    std::optional<FixedArray1D<Real, 3>> r_vec = std::nullopt
);
extern "C" bool fortran_w_mat_for_tilt(
    double &tilt /* 0D_NOT_real in */,
    bool *return_inverse /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3>
w_mat_for_tilt(double tilt, std::optional<bool> return_inverse = std::nullopt);
extern "C" bool fortran_w_mat_for_x_pitch(
    double &x_pitch /* 0D_NOT_real in */,
    bool *return_inverse /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3>
w_mat_for_x_pitch(double x_pitch, std::optional<bool> return_inverse = std::nullopt);
extern "C" bool fortran_w_mat_for_y_pitch(
    double &y_pitch /* 0D_NOT_real in */,
    bool *return_inverse /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3>
w_mat_for_y_pitch(double y_pitch, std::optional<bool> return_inverse = std::nullopt);
extern "C" bool fortran_wall3d_d_radius(
    Bmad::array_descriptor_t &position /* 1D_NOT_real in */,
    void *ele /* 0D_NOT_type in */,
    int *ix_wall /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &perp /* 1D_NOT_real out */,
    int &ix_section /* 0D_NOT_integer out */,
    bool &no_wall_here /* 0D_NOT_logical out */,
    Bmad::array_descriptor_t &origin /* 1D_NOT_real out */,
    double &radius_wall /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */,
    double &d_radius /* 0D_NOT_real out */
);
struct Wall3dDRadius {
  FixedArray1D<Real, 3> perp;
  int ix_section;
  bool no_wall_here;
  FixedArray1D<Real, 3> origin;
  double radius_wall;
  bool err_flag;
  double d_radius;
};
Bmad::Wall3dDRadius wall3d_d_radius(
    FArray1D<Real> &position,
    EleStruct &ele,
    std::optional<int> ix_wall = std::nullopt
);
extern "C" void fortran_wall3d_initializer(
    void *wall3d /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */
);
bool wall3d_initializer(Wall3dStruct &wall3d);
extern "C" void fortran_wall3d_section_initializer(
    void *section /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */
);
bool wall3d_section_initializer(Wall3dSectionStruct &section);
extern "C" bool fortran_wall3d_to_position(
    void *orbit /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &position /* 1D_NOT_real out */
);
FixedArray1D<Real, 6> wall3d_to_position(CoordStruct &orbit, EleStruct &ele);

// Skipped unusable routine wall_hit_handler_custom_def:
// - Routine in configuration skip list

// Skipped unusable routine wiggler_field:
// - Module name unset
extern "C" void fortran_word_to_value(
    const char *word /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */,
    double &value /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical in */,
    void *ele /* 0D_NOT_type inout */
);
void word_to_value(
    std::string word,
    LatStruct &lat,
    double value,
    bool err_flag,
    optional_ref<EleStruct> ele = std::nullopt
);

// Skipped unusable routine write_2d:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_write_ascii_beam_file(
    const char *file_name /* 0D_NOT_character in */,
    void *beam /* 0D_NOT_type in */,
    bool *new_file /* 0D_NOT_logical in */,
    bool *alive_only /* 0D_NOT_logical in */
);
void write_ascii_beam_file(
    std::string file_name,
    BeamStruct &beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<bool> alive_only = std::nullopt
);
extern "C" void fortran_write_astra_bend(
    int &iu /* 0D_NOT_integer in */,
    double &strength /* 0D_NOT_real in */,
    int &id /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &d1 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &d2 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &d3 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &d4 /* 1D_NOT_real inout */
);
void write_astra_bend(
    int iu,
    double strength,
    int id,
    FixedArray1D<Real, 2> d1,
    FixedArray1D<Real, 2> d2,
    FixedArray1D<Real, 2> d3,
    FixedArray1D<Real, 2> d4
);

// Skipped unusable routine write_astra_ele:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_write_astra_field_grid_file(
    int &astra_file_unit /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    double *dz /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteAstraFieldGridFile {
  double maxfield;
  bool err;
};
Bmad::WriteAstraFieldGridFile write_astra_field_grid_file(
    int astra_file_unit,
    EleStruct &ele,
    std::optional<double> dz = std::nullopt
);
extern "C" void fortran_write_astra_field_grid_file_3d(
    const char *base_filename /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    double *dz /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteAstraFieldGridFile3d {
  double maxfield;
  bool err;
};
Bmad::WriteAstraFieldGridFile3d write_astra_field_grid_file_3d(
    std::string base_filename,
    EleStruct &ele,
    std::optional<double> dz = std::nullopt
);

// Skipped unusable routine write_astra_lattice_file:
// - Untranslated type: astra_lattice_param_struct (0D)
extern "C" void fortran_write_beam_file(
    const char *file_name /* 0D_NOT_character in */,
    void *beam /* 0D_NOT_type in */,
    bool *new_file /* 0D_NOT_logical in */,
    int *file_format /* 0D_NOT_integer in */,
    void *lat /* 0D_NOT_type in */,
    bool *alive_only /* 0D_NOT_logical in */
);
void write_beam_file(
    std::string file_name,
    BeamStruct &beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<int> file_format = std::nullopt,
    optional_ref<LatStruct> lat = std::nullopt,
    std::optional<bool> alive_only = std::nullopt
);
extern "C" void fortran_write_beam_floor_positions(
    const char *file_name /* 0D_NOT_character in */,
    void *beam /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    bool *new_file /* 0D_NOT_logical in */
);
void write_beam_floor_positions(
    std::string file_name,
    BeamStruct &beam,
    EleStruct &ele,
    std::optional<bool> new_file = std::nullopt
);
extern "C" void fortran_write_binary_cartesian_map(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *cart_map /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool write_binary_cartesian_map(
    std::string file_name,
    EleStruct &ele,
    CartesianMapStruct &cart_map
);
extern "C" void fortran_write_binary_cylindrical_map(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *cl_map /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool write_binary_cylindrical_map(
    std::string file_name,
    EleStruct &ele,
    CylindricalMapStruct &cl_map
);
extern "C" void fortran_write_binary_grid_field(
    const char *file_name /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    void *g_field /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool write_binary_grid_field(std::string file_name, EleStruct &ele, GridFieldStruct &g_field);
extern "C" void fortran_write_blender_ele(
    int &iu /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type inout */,
    bool *old_format /* 0D_NOT_logical in */
);
void write_blender_ele(int iu, EleStruct &ele, std::optional<bool> old_format = std::nullopt);
extern "C" void fortran_write_blender_lat_layout(
    const char *file_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */
);
void write_blender_lat_layout(std::string file_name, LatStruct &lat);
extern "C" void fortran_write_bmad_lattice_file(
    const char *bmad_file /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */,
    int *output_form /* 0D_NOT_integer in */,
    void *orbit0 /* 0D_NOT_type in */
);
bool write_bmad_lattice_file(
    std::string bmad_file,
    LatStruct &lat,
    std::optional<int> output_form = std::nullopt,
    optional_ref<CoordStruct> orbit0 = std::nullopt
);

// Skipped unusable routine write_digested_bmad_file:
// - Variable-sized in character array: 1D_NOT_character
// - Untranslated type: extra_parsing_info_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine write_gpt_ele:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_write_gpt_field_grid_file_1d(
    int &gpt_file_unit /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    double &ref_time /* 0D_NOT_real out */,
    double *dz /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteGptFieldGridFile1d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile1d write_gpt_field_grid_file_1d(
    int gpt_file_unit,
    EleStruct &ele,
    std::optional<double> dz = std::nullopt
);
extern "C" void fortran_write_gpt_field_grid_file_2d(
    int &gpt_file_unit /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    double &ref_time /* 0D_NOT_real out */,
    double *dr /* 0D_NOT_real in */,
    double *dz /* 0D_NOT_real in */,
    double *r_max /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteGptFieldGridFile2d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile2d write_gpt_field_grid_file_2d(
    int gpt_file_unit,
    EleStruct &ele,
    std::optional<double> dr = std::nullopt,
    std::optional<double> dz = std::nullopt,
    std::optional<double> r_max = std::nullopt
);
extern "C" void fortran_write_gpt_field_grid_file_3d(
    const char *base_filename /* 0D_NOT_character in */,
    void *ele /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    double &ref_time /* 0D_NOT_real out */,
    double *dz /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteGptFieldGridFile3d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile3d write_gpt_field_grid_file_3d(
    std::string base_filename,
    EleStruct &ele,
    std::optional<double> dz = std::nullopt
);

// Skipped unusable routine write_gpt_lattice_file:
// - Untranslated type: gpt_lat_param_struct (0D)
extern "C" void fortran_write_lat_line(
    const char *line /* 0D_NOT_character inout */,
    int &iu /* 0D_NOT_integer in */,
    bool &end_is_neigh /* 0D_NOT_logical in */,
    bool *do_split /* 0D_NOT_logical in */,
    bool *ampersand_at_ends /* 0D_NOT_logical in */
);
void write_lat_line(
    std::string &line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split = std::nullopt,
    std::optional<bool> ampersand_at_ends = std::nullopt
);
extern "C" void fortran_write_lattice_in_elegant_format(
    const char *out_file_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    void *ref_orbit /* 1D_ALLOC_type in */,
    bool *use_matrix_model /* 0D_NOT_logical in */,
    bool *include_apertures /* 0D_NOT_logical in */,
    double *dr12_drift_max /* 0D_NOT_real in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical out */
);
bool write_lattice_in_elegant_format(
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
);
extern "C" void fortran_write_lattice_in_foreign_format(
    const char *out_type /* 0D_NOT_character in */,
    const char *out_file_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    void *ref_orbit /* 1D_ALLOC_type in */,
    bool *use_matrix_model /* 0D_NOT_logical in */,
    bool *include_apertures /* 0D_NOT_logical in */,
    double *dr12_drift_max /* 0D_NOT_real in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical out */
);
bool write_lattice_in_foreign_format(
    std::string out_type,
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
);
extern "C" void fortran_write_lattice_in_mad_format(
    const char *out_type /* 0D_NOT_character in */,
    const char *out_file_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    void *ref_orbit /* 1D_ALLOC_type in */,
    bool *use_matrix_model /* 0D_NOT_logical in */,
    bool *include_apertures /* 0D_NOT_logical in */,
    double *dr12_drift_max /* 0D_NOT_real in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical out */
);
bool write_lattice_in_mad_format(
    std::string out_type,
    std::string out_file_name,
    LatStruct &lat,
    std::optional<CoordStructAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
);
extern "C" void fortran_write_lattice_in_pals(
    const char *pals_file /* 0D_NOT_character out */,
    void *lat /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct WriteLatticeInPals {
  std::string pals_file;
  bool err_flag;
};
Bmad::WriteLatticeInPals write_lattice_in_pals(LatStruct &lat);
extern "C" void fortran_write_lattice_in_sad_format(
    const char *out_file_name /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type inout */,
    bool *include_apertures /* 0D_NOT_logical in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *err /* 0D_NOT_logical in */
);
void write_lattice_in_sad_format(
    std::string out_file_name,
    LatStruct &lat,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> err = std::nullopt
);
extern "C" void fortran_write_lattice_in_scibmad(
    const char *scibmad_file /* 0D_NOT_character out */,
    void *lat /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct WriteLatticeInScibmad {
  std::string scibmad_file;
  bool err_flag;
};
Bmad::WriteLatticeInScibmad write_lattice_in_scibmad(LatStruct &lat);
extern "C" void fortran_write_line_element(
    const char *line /* 0D_NOT_character in */,
    int &iu /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type inout */,
    void *lat /* 0D_NOT_type inout */
);
void write_line_element(std::string line, int iu, EleStruct &ele, LatStruct &lat);
extern "C" void fortran_write_opal_field_grid_file(
    int &opal_file_unit /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    void *param /* 0D_NOT_type in */,
    double &maxfield /* 0D_NOT_real out */,
    bool &err /* 0D_NOT_logical out */
);
struct WriteOpalFieldGridFile {
  double maxfield;
  bool err;
};
Bmad::WriteOpalFieldGridFile
write_opal_field_grid_file(int opal_file_unit, EleStruct &ele, LatParamStruct &param);
extern "C" void fortran_write_opal_lattice_file(
    int &opal_file_unit /* 0D_NOT_integer in */,
    void *lat /* 0D_NOT_type in */,
    bool &err /* 0D_NOT_logical out */
);
bool write_opal_lattice_file(int opal_file_unit, LatStruct &lat);
extern "C" void fortran_write_time_particle_distribution(
    int &time_file_unit /* 0D_NOT_integer in */,
    void *bunch /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    const char *style /* 0D_NOT_character in */,
    void *branch /* 0D_NOT_type in */,
    const char *format /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical out */
);
bool write_time_particle_distribution(
    int time_file_unit,
    BunchStruct &bunch,
    EleStruct &ele,
    std::optional<std::string> style = std::nullopt,
    optional_ref<BranchStruct> branch = std::nullopt,
    std::optional<std::string> format = std::nullopt
);
extern "C" bool fortran_xlafun(
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    double &z /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void xlafun(double x, double y, double z, double res);
extern "C" bool fortran_xraylib_nist_compound(
    const char *name /* 0D_NOT_character in */,
    int &indx /* 0D_NOT_integer out */
);
int xraylib_nist_compound(std::string name);

// Skipped unusable routine xsif_parser:
// - Routine in configuration skip list

// Skipped unusable routine xyz_to_action:
// - Untranslated type: * (0D)
extern "C" bool fortran_ylafun(
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    double &z /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void ylafun(double x, double y, double z, double res);
extern "C" bool fortran_z_at_surface(
    void *ele /* 0D_NOT_type in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *extend_grid /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &dz_dxy /* 1D_NOT_real out */,
    double &z /* 0D_NOT_real out */
);
struct ZAtSurface {
  bool err_flag;
  FixedArray1D<Real, 2> dz_dxy;
  double z;
};
Bmad::ZAtSurface
z_at_surface(EleStruct &ele, double x, double y, std::optional<bool> extend_grid = std::nullopt);
extern "C" void fortran_zero_ele_kicks(void *ele /* 0D_NOT_type inout */);
void zero_ele_kicks(EleStruct &ele);
extern "C" void fortran_zero_ele_offsets(void *ele /* 0D_NOT_type inout */);
void zero_ele_offsets(EleStruct &ele);
extern "C" void fortran_zero_lr_wakes_in_lat(void *lat /* 0D_NOT_type inout */);
void zero_lr_wakes_in_lat(LatStruct &lat);
extern "C" bool fortran_zlafun(
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */,
    double &z /* 0D_NOT_real in */,
    double &res /* 0D_NOT_real in */
);
void zlafun(double x, double y, double z, double res);

// Skipped unusable routine zot_integrand:
// - Untranslated type: c_ptr (0D)
} // namespace Bmad
