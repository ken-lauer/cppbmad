#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace Bmad {
extern "C" void fortran_ab_multipole_kick(
    double& a /* 0D_NOT_real in */,
    double& b /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer in */,
    int& ref_species /* 0D_NOT_integer in */,
    int& ele_orientation /* 0D_NOT_integer in */,
    void* coord /* 0D_NOT_type in */,
    double& kx /* 0D_NOT_real out */,
    double& ky /* 0D_NOT_real out */,
    double* dk /* 2D_NOT_real out */,
    int* pole_type /* 0D_NOT_integer in */,
    double* scale /* 0D_NOT_real in */);
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
    CoordProxy& coord,
    std::optional<int> pole_type = std::nullopt,
    std::optional<double> scale = std::nullopt);
extern "C" void fortran_ab_multipole_kicks(
    void* an /* 1D_ALLOC_real in */,
    void* bn /* 1D_ALLOC_real in */,
    int& ix_pole_max /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    int* pole_type /* 0D_NOT_integer in */,
    double* scale /* 0D_NOT_real in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void ab_multipole_kicks(
    RealAlloc1D& an,
    RealAlloc1D& bn,
    int ix_pole_max,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<int> pole_type = std::nullopt,
    std::optional<double> scale = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_absolute_photon_position(
    void* e_orb /* 0D_NOT_type in */,
    void* photon_orb /* 0D_NOT_type inout */);
void absolute_photon_position(CoordProxy& e_orb, CoordProxy& photon_orb);
extern "C" bool fortran_absolute_time_tracking(
    void* ele /* 0D_NOT_type in */,
    bool& is_abs_time /* 0D_NOT_logical inout */);
void absolute_time_tracking(EleProxy& ele, bool& is_abs_time);
extern "C" bool fortran_ac_kicker_amp(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    double* true_time /* 0D_NOT_real in */,
    double& ac_amp /* 0D_NOT_real inout */);
void ac_kicker_amp(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<double> true_time,
    double& ac_amp);
extern "C" void fortran_action_to_xyz(
    void* ring /* 0D_NOT_type in */,
    int& ix /* 0D_NOT_integer in */,
    double* J /* 1D_NOT_real in */,
    double* X /* 1D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ActionToXyz {
  FixedArray1D<Real, 6> X;
  bool err_flag;
};
Bmad::ActionToXyz action_to_xyz(
    LatProxy& ring,
    int ix,
    FixedArray1D<Real, 6> J);
extern "C" void fortran_add_lattice_control_structs(
    void* ele /* 0D_NOT_type in */,
    int* n_add_slave /* 0D_NOT_integer in */,
    int* n_add_lord /* 0D_NOT_integer in */,
    int* n_add_slave_field /* 0D_NOT_integer in */,
    int* n_add_lord_field /* 0D_NOT_integer in */,
    bool* add_at_end /* 0D_NOT_logical in */);
void add_lattice_control_structs(
    EleProxy& ele,
    std::optional<int> n_add_slave = std::nullopt,
    std::optional<int> n_add_lord = std::nullopt,
    std::optional<int> n_add_slave_field = std::nullopt,
    std::optional<int> n_add_lord_field = std::nullopt,
    std::optional<bool> add_at_end = std::nullopt);

// Skipped unusable routine add_ptc_layout_to_list:
// - Untranslated type: ptc_branch1_struct (0D)
// - Untranslated type: layout (0D)
extern "C" void fortran_add_superimpose(
    void* lat /* 0D_NOT_type inout */,
    void* super_ele_in /* 0D_NOT_type in */,
    int& ix_branch /* 0D_NOT_integer in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* super_ele_out /* 0D_PTR_type out */,
    bool* save_null_drift /* 0D_NOT_logical in */,
    bool* create_jumbo_slave /* 0D_NOT_logical in */,
    int* ix_insert /* 0D_NOT_integer in */,
    bool* mangle_slave_names /* 0D_NOT_logical in */,
    bool* wrap /* 0D_NOT_logical in */);
struct AddSuperimpose {
  bool err_flag;
  EleProxy super_ele_out;
};
Bmad::AddSuperimpose add_superimpose(
    LatProxy& lat,
    EleProxy& super_ele_in,
    int ix_branch,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> create_jumbo_slave = std::nullopt,
    std::optional<int> ix_insert = std::nullopt,
    std::optional<bool> mangle_slave_names = std::nullopt,
    std::optional<bool> wrap = std::nullopt);
extern "C" void fortran_add_this_multipass(
    void* lat /* 0D_NOT_type inout */,
    void* m_slaves /* 1D_ALLOC_type inout */,
    void* lord_in /* 0D_NOT_type inout */);
void add_this_multipass(
    LatProxy& lat,
    LatEleLocProxyAlloc1D& m_slaves,
    optional_ref<EleProxy> lord_in = std::nullopt);

// Skipped unusable routine add_this_name_to_list:
// - Variable-sized inout character array: names(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_add_this_taylor_term(
    void* ele /* 0D_NOT_type inout */,
    int& i_out /* 0D_NOT_integer inout */,
    double& coef /* 0D_NOT_real inout */,
    int* expn /* 1D_NOT_integer inout */);
void add_this_taylor_term(
    EleProxy& ele,
    int& i_out,
    double& coef,
    FixedArray1D<Int, 6> expn);
extern "C" void fortran_adjust_super_slave_names(
    void* lat /* 0D_NOT_type inout */,
    int& ix1_lord /* 0D_NOT_integer inout */,
    int& ix2_lord /* 0D_NOT_integer inout */,
    bool* first_time /* 0D_NOT_logical inout */);
void adjust_super_slave_names(
    LatProxy& lat,
    int& ix1_lord,
    int& ix2_lord,
    optional_ref<bool> first_time = std::nullopt);
extern "C" void fortran_allocate_branch_array(
    void* lat /* 0D_NOT_type inout */,
    int& upper_bound /* 0D_NOT_integer in */);
void allocate_branch_array(LatProxy& lat, int upper_bound);

// Skipped unusable routine allocate_element_array:
// - Routine in configuration skip list
extern "C" void fortran_allocate_lat_ele_array(
    void* lat /* 0D_NOT_type inout */,
    int* upper_bound /* 0D_NOT_integer in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* do_ramper_slave_setup /* 0D_NOT_logical in */);
void allocate_lat_ele_array(
    LatProxy& lat,
    std::optional<int> upper_bound = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> do_ramper_slave_setup = std::nullopt);

// Skipped unusable routine allocate_plat:
// - Untranslated type: parser_lat_struct (0D)

// Skipped unusable routine aml_parser:
// - Routine in configuration skip list
extern "C" bool fortran_angle_between_polars(
    void* polar1 /* 0D_NOT_type in */,
    void* polar2 /* 0D_NOT_type in */,
    double& angle /* 0D_NOT_real inout */);
void angle_between_polars(
    SpinPolarProxy& polar1,
    SpinPolarProxy& polar2,
    double& angle);
extern "C" void fortran_angle_to_canonical_coords(
    void* orbit /* 0D_NOT_type inout */,
    const char* coord_type /* 0D_NOT_character in */);
void angle_to_canonical_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type = std::nullopt);
extern "C" void fortran_aperture_bookkeeper(void* ele /* 0D_NOT_type inout */);
void aperture_bookkeeper(EleProxy& ele);
extern "C" void fortran_apply_all_rampers(
    void* lat /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool apply_all_rampers(LatProxy& lat);

// Skipped unusable routine apply_element_edge_kick:
// - Untranslated type: fringe_field_info_struct (0D)

// Skipped unusable routine apply_element_edge_kick_hook_def:
// - Untranslated type: fringe_field_info_struct (0D)
extern "C" void fortran_apply_energy_kick(
    double& dE /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type inout */,
    double* ddE_dr /* 1D_NOT_real in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void apply_energy_kick(
    double dE,
    CoordProxy& orbit,
    FixedArray1D<Real, 2> ddE_dr,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_apply_patch_to_ptc_fibre(
    void* ele /* 0D_NOT_type in */);
void apply_patch_to_ptc_fibre(EleProxy& ele);
extern "C" void fortran_apply_rampers_to_slave(
    void* slave /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool apply_rampers_to_slave(EleProxy& slave);
extern "C" bool fortran_array_re_str(
    void* arr /* 1D_ALLOC_real inout */,
    const char* parens_in /* 0D_NOT_character inout */,
    const char* str_out /* 0D_NOT_character inout */);
void array_re_str(
    RealAlloc1D& arr,
    optional_ref<std::string> parens_in,
    std::string& str_out);
extern "C" bool fortran_astra_max_field_reference(
    void* pt0 /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */,
    double& field_value /* 0D_NOT_real inout */);
void astra_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double& field_value);
extern "C" bool fortran_at_this_ele_end(
    int& now_at /* 0D_NOT_integer in */,
    int& where_at /* 0D_NOT_integer in */,
    bool& is_at_this_end /* 0D_NOT_logical inout */);
void at_this_ele_end(int now_at, int where_at, bool& is_at_this_end);
extern "C" void fortran_attribute_bookkeeper(
    void* ele /* 0D_NOT_type inout */,
    bool* force_bookkeeping /* 0D_NOT_logical in */);
void attribute_bookkeeper(
    EleProxy& ele,
    std::optional<bool> force_bookkeeping = std::nullopt);
extern "C" bool fortran_attribute_free1(
    int& ix_ele /* 0D_NOT_integer inout */,
    const char* attrib_name /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    bool* err_print_flag /* 0D_NOT_logical inout */,
    bool* except_overlay /* 0D_NOT_logical inout */,
    bool* dependent_attribs_free /* 0D_NOT_logical inout */,
    int* why_not_free /* 0D_NOT_integer inout */,
    bool& free /* 0D_NOT_logical inout */);
void attribute_free1(
    int& ix_ele,
    std::string& attrib_name,
    LatProxy& lat,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool& free);
extern "C" bool fortran_attribute_free2(
    void* ele /* 0D_NOT_type inout */,
    const char* attrib_name /* 0D_NOT_character inout */,
    bool* err_print_flag /* 0D_NOT_logical inout */,
    bool* except_overlay /* 0D_NOT_logical inout */,
    bool* dependent_attribs_free /* 0D_NOT_logical inout */,
    int* why_not_free /* 0D_NOT_integer inout */,
    bool& free /* 0D_NOT_logical inout */);
void attribute_free2(
    EleProxy& ele,
    std::string& attrib_name,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool& free);
extern "C" bool fortran_attribute_free3(
    int& ix_ele /* 0D_NOT_integer inout */,
    int& ix_branch /* 0D_NOT_integer inout */,
    const char* attrib_name /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    bool* err_print_flag /* 0D_NOT_logical inout */,
    bool* except_overlay /* 0D_NOT_logical inout */,
    bool* dependent_attribs_free /* 0D_NOT_logical inout */,
    int* why_not_free /* 0D_NOT_integer inout */,
    bool& free /* 0D_NOT_logical inout */);
void attribute_free3(
    int& ix_ele,
    int& ix_branch,
    std::string& attrib_name,
    LatProxy& lat,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool& free);
extern "C" bool fortran_attribute_index1(
    void* ele /* 0D_NOT_type inout */,
    const char* name /* 0D_NOT_character inout */,
    const char* full_name /* 0D_NOT_character inout */,
    bool* can_abbreviate /* 0D_NOT_logical inout */,
    bool* print_error /* 0D_NOT_logical inout */,
    int& attrib_index /* 0D_NOT_integer inout */);
void attribute_index1(
    EleProxy& ele,
    std::string& name,
    optional_ref<std::string> full_name,
    optional_ref<bool> can_abbreviate,
    optional_ref<bool> print_error,
    int& attrib_index);
extern "C" bool fortran_attribute_index2(
    int& key /* 0D_NOT_integer inout */,
    const char* name /* 0D_NOT_character inout */,
    const char* full_name /* 0D_NOT_character inout */,
    bool* can_abbreviate /* 0D_NOT_logical inout */,
    bool* print_error /* 0D_NOT_logical inout */,
    int& attrib_index /* 0D_NOT_integer inout */);
void attribute_index2(
    int& key,
    std::string& name,
    optional_ref<std::string> full_name,
    optional_ref<bool> can_abbreviate,
    optional_ref<bool> print_error,
    int& attrib_index);

// Skipped unusable routine attribute_info:
// - Untranslated type: ele_attribute_struct (0D)
extern "C" bool fortran_attribute_name1(
    int& key /* 0D_NOT_integer inout */,
    int& ix_att /* 0D_NOT_integer inout */,
    bool* show_private /* 0D_NOT_logical inout */,
    const char* attrib_name /* 0D_NOT_character inout */);
void attribute_name1(
    int& key,
    int& ix_att,
    optional_ref<bool> show_private,
    std::string& attrib_name);
extern "C" bool fortran_attribute_name2(
    void* ele /* 0D_NOT_type inout */,
    int& ix_att /* 0D_NOT_integer inout */,
    bool* show_private /* 0D_NOT_logical inout */,
    const char* attrib_name /* 0D_NOT_character inout */);
void attribute_name2(
    EleProxy& ele,
    int& ix_att,
    optional_ref<bool> show_private,
    std::string& attrib_name);

// Skipped unusable routine attribute_set_bookkeeping:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_attribute_type(
    const char* attrib_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    int& attrib_type /* 0D_NOT_integer out */);
int attribute_type(
    std::string attrib_name,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" bool fortran_attribute_units(
    const char* attrib_name /* 0D_NOT_character in */,
    const char* unrecognized_units /* 0D_NOT_character in */,
    const char* attrib_units /* 0D_NOT_character out */);
std::string attribute_units(
    std::string attrib_name,
    std::optional<std::string> unrecognized_units = std::nullopt);
extern "C" void fortran_autoscale_phase_and_amp(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* scale_phase /* 0D_NOT_logical in */,
    bool* scale_amp /* 0D_NOT_logical in */,
    bool* call_bookkeeper /* 0D_NOT_logical in */);
bool autoscale_phase_and_amp(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> scale_phase = std::nullopt,
    std::optional<bool> scale_amp = std::nullopt,
    std::optional<bool> call_bookkeeper = std::nullopt);
extern "C" bool fortran_average_twiss(
    double& frac1 /* 0D_NOT_real in */,
    void* twiss1 /* 0D_NOT_type in */,
    void* twiss2 /* 0D_NOT_type inout */,
    void* ave_twiss /* 0D_NOT_type inout */);
void average_twiss(
    double frac1,
    TwissProxy& twiss1,
    TwissProxy& twiss2,
    TwissProxy& ave_twiss);

// Skipped unusable routine bane1:
// - Untranslated type: ibs_struct (0D)
extern "C" void fortran_bbi_kick(
    double& x /* 0D_NOT_real in */,
    double& y /* 0D_NOT_real in */,
    double* sigma /* 1D_NOT_real in */,
    double* nk /* 1D_NOT_real out */,
    double* dnk /* 2D_NOT_real out */);
struct BbiKick {
  FixedArray1D<Real, 2> nk;
  FixedArray2D<Real, 2, 2> dnk;
};
Bmad::BbiKick bbi_kick(double x, double y, FixedArray1D<Real, 2> sigma);
extern "C" void fortran_bbi_slice_calc(
    void* ele /* 0D_NOT_type in */,
    int& n_slice /* 0D_NOT_integer in */,
    void* z_slice /* 1D_ALLOC_real out */);
RealAlloc1D bbi_slice_calc(EleProxy& ele, int n_slice);
extern "C" void fortran_beam_envelope_ibs(
    double* sigma_mat /* 2D_NOT_real in */,
    double* ibs_mat /* 2D_NOT_real out */,
    bool& tail_cut /* 0D_NOT_logical in */,
    double& tau /* 0D_NOT_real in */,
    double& energy /* 0D_NOT_real in */,
    double& n_part /* 0D_NOT_real in */,
    int& species /* 0D_NOT_integer in */);
FixedArray2D<Real, 6, 6> beam_envelope_ibs(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species);
extern "C" void fortran_beam_equal_beam(
    void* beam1 /* 0D_NOT_type in */,
    void* beam2 /* 0D_NOT_type in */);
void beam_equal_beam(BeamProxy& beam1, BeamProxy& beam2);
extern "C" void fortran_beam_tilts(
    double* S /* 2D_NOT_real in */,
    double& angle_xy /* 0D_NOT_real out */,
    double& angle_xz /* 0D_NOT_real out */,
    double& angle_yz /* 0D_NOT_real out */,
    double& angle_xpz /* 0D_NOT_real out */,
    double& angle_ypz /* 0D_NOT_real out */);
struct BeamTilts {
  double angle_xy;
  double angle_xz;
  double angle_yz;
  double angle_xpz;
  double angle_ypz;
};
Bmad::BeamTilts beam_tilts(FixedArray2D<Real, 6, 6> S);

// Skipped unusable routine beambeam_fibre_setup:
// - Untranslated type: fibre (0D)
extern "C" void fortran_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orb /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */,
    bool* track_spin /* 0D_NOT_logical in */);
void bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> track_spin = std::nullopt);
extern "C" void fortran_bend_exact_multipole_field(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    bool& local_ref_frame /* 0D_NOT_logical in */,
    void* field /* 0D_NOT_type out */,
    bool* calc_dfield /* 0D_NOT_logical in */,
    bool* calc_potential /* 0D_NOT_logical in */);
EmFieldProxy bend_exact_multipole_field(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" bool fortran_bend_length_has_been_set(
    void* ele /* 0D_NOT_type in */,
    bool& is_set /* 0D_NOT_logical inout */);
void bend_length_has_been_set(EleProxy& ele, bool& is_set);
extern "C" bool fortran_bend_photon_e_rel_init(
    double* r_in /* 0D_NOT_real in */,
    double& E_rel /* 0D_NOT_real out */);
double bend_photon_e_rel_init(std::optional<double> r_in = std::nullopt);
extern "C" bool fortran_bend_photon_energy_integ_prob(
    double& E_photon /* 0D_NOT_real in */,
    double& g_bend /* 0D_NOT_real in */,
    double& gamma /* 0D_NOT_real in */,
    double& integ_prob /* 0D_NOT_real out */);
double bend_photon_energy_integ_prob(
    double E_photon,
    double g_bend,
    double gamma);
extern "C" bool fortran_bend_photon_energy_normalized_probability(
    double& E_rel /* 0D_NOT_real in */,
    double& prob /* 0D_NOT_real out */);
double bend_photon_energy_normalized_probability(double E_rel);
extern "C" void fortran_bend_photon_init(
    double& g_bend_x /* 0D_NOT_real in */,
    double& g_bend_y /* 0D_NOT_real in */,
    double& gamma /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type out */,
    double* E_min /* 0D_NOT_real in */,
    double* E_max /* 0D_NOT_real in */,
    double* E_integ_prob /* 0D_NOT_real in */,
    double* vert_angle_min /* 0D_NOT_real in */,
    double* vert_angle_max /* 0D_NOT_real in */,
    bool* vert_angle_symmetric /* 0D_NOT_logical in */,
    double* emit_probability /* 0D_NOT_real in */);
CoordProxy bend_photon_init(
    double g_bend_x,
    double g_bend_y,
    double gamma,
    std::optional<double> E_min = std::nullopt,
    std::optional<double> E_max = std::nullopt,
    std::optional<double> E_integ_prob = std::nullopt,
    std::optional<double> vert_angle_min = std::nullopt,
    std::optional<double> vert_angle_max = std::nullopt,
    std::optional<bool> vert_angle_symmetric = std::nullopt,
    std::optional<double> emit_probability = std::nullopt);
extern "C" void fortran_bend_photon_polarization_init(
    double& g_bend_x /* 0D_NOT_real in */,
    double& g_bend_y /* 0D_NOT_real in */,
    double& E_rel /* 0D_NOT_real in */,
    double& gamma_phi /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type out */);
CoordProxy bend_photon_polarization_init(
    double g_bend_x,
    double g_bend_y,
    double E_rel,
    double gamma_phi);
extern "C" bool fortran_bend_photon_vert_angle_init(
    double& E_rel /* 0D_NOT_real in */,
    double& gamma /* 0D_NOT_real in */,
    double* r_in /* 0D_NOT_real in */,
    bool* invert /* 0D_NOT_logical in */,
    double& phi /* 0D_NOT_real out */);
double bend_photon_vert_angle_init(
    double E_rel,
    double gamma,
    std::optional<double> r_in = std::nullopt,
    std::optional<bool> invert = std::nullopt);
extern "C" bool fortran_bend_shift(
    void* position1 /* 0D_NOT_type in */,
    double& g /* 0D_NOT_real in */,
    double& delta_s /* 0D_NOT_real in */,
    double* w_mat /* 2D_NOT_real out */,
    double* ref_tilt /* 0D_NOT_real in */,
    void* position2 /* 0D_NOT_type inout */);
FixedArray2D<Real, 3, 3> bend_shift(
    FloorPositionProxy& position1,
    double g,
    double delta_s,
    std::optional<double> ref_tilt,
    FloorPositionProxy& position2);
extern "C" bool fortran_bend_vert_angle_integ_prob(
    double& vert_angle /* 0D_NOT_real in */,
    double& E_rel /* 0D_NOT_real in */,
    double& gamma /* 0D_NOT_real in */,
    double& integ_prob /* 0D_NOT_real out */);
double bend_vert_angle_integ_prob(
    double vert_angle,
    double E_rel,
    double gamma);

// Skipped unusable routine bjmt1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine bjmt_integrand:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine bl_via_mat:
// - Untranslated type: ibs_sim_param_struct (0D)
extern "C" void fortran_bl_via_vlassov(
    double& current /* 0D_NOT_real in */,
    double& alpha /* 0D_NOT_real in */,
    double& Energy /* 0D_NOT_real in */,
    double& sigma_p /* 0D_NOT_real in */,
    double& Vrf /* 0D_NOT_real in */,
    double& omega /* 0D_NOT_real in */,
    double& U0 /* 0D_NOT_real in */,
    double& circ /* 0D_NOT_real in */,
    double& R /* 0D_NOT_real in */,
    double& L /* 0D_NOT_real in */,
    double& sigma_z /* 0D_NOT_real out */);
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
    double L);

// Skipped unusable routine bmad_and_xsif_parser:
// - Routine in configuration skip list
extern "C" void fortran_bmad_parser(
    const char* lat_file /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type out */,
    bool* make_mats6 /* 0D_NOT_logical in */,
    bool& digested_read_ok /* 0D_NOT_logical out */,
    const char* use_line /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* parse_lat /* 0D_NOT_type out */);
struct BmadParser {
  LatProxy lat;
  bool digested_read_ok;
  bool err_flag;
  LatProxy parse_lat;
};
Bmad::BmadParser bmad_parser(
    std::string lat_file,
    std::optional<bool> make_mats6 = std::nullopt,
    std::optional<std::string> use_line = std::nullopt);
extern "C" void fortran_bmad_parser2(
    const char* lat_file /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type inout */,
    void* orbit /* 1D_ALLOC_type in */,
    bool* make_mats6 /* 0D_NOT_logical in */,
    bool* err_flag /* 0D_NOT_logical inout */,
    void* parse_lat /* 0D_NOT_type in */);
void bmad_parser2(
    std::string lat_file,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> orbit = std::nullopt,
    std::optional<bool> make_mats6 = std::nullopt,
    optional_ref<bool> err_flag = std::nullopt,
    optional_ref<LatProxy> parse_lat = std::nullopt);

// Skipped unusable routine bmad_parser_string_attribute_set:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_bmad_patch_parameters_to_ptc(
    double* ang /* 1D_NOT_real inout */,
    double* exi /* 2D_NOT_real inout */);
void bmad_patch_parameters_to_ptc(
    FixedArray1D<Real, 3> ang,
    FixedArray2D<Real, 3, 3> exi);

// Skipped unusable routine bmad_taylor_equal_damap:
// - Untranslated type: damap (0D)

// Skipped unusable routine bmad_taylors_equal_ptc_taylors:
// - Untranslated type: taylor (1D)

// Skipped unusable routine bmad_taylors_equal_reals_8:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_bp_set_ran_status();
void bp_set_ran_status();
extern "C" void fortran_branch_equal_branch(
    void* branch1 /* 0D_NOT_type out */,
    void* branch2 /* 0D_NOT_type in */);
BranchProxy branch_equal_branch(BranchProxy& branch2);
extern "C" bool fortran_branch_name(
    void* branch /* 0D_NOT_type in */,
    const char* name /* 0D_NOT_character inout */);
void branch_name(BranchProxy& branch, std::string& name);
extern "C" void fortran_branch_to_ptc_m_u(void* branch /* 0D_NOT_type in */);
void branch_to_ptc_m_u(BranchProxy& branch);
extern "C" void fortran_bunch_equal_bunch(
    void* bunch1 /* 0D_NOT_type in */,
    void* bunch2 /* 0D_NOT_type in */);
void bunch_equal_bunch(BunchProxy& bunch1, BunchProxy& bunch2);

// Skipped unusable routine c_multi:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_c_to_cbar(
    void* ele /* 0D_NOT_type in */,
    double* cbar_mat /* 2D_NOT_real out */);
FixedArray2D<Real, 2, 2> c_to_cbar(EleProxy& ele);
extern "C" void fortran_calc_bunch_params(
    void* bunch /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type in */,
    bool& error /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */,
    double* n_mat /* 2D_NOT_real in */,
    bool* is_time_coords /* 0D_NOT_logical in */,
    void* ele /* 0D_NOT_type in */);
void calc_bunch_params(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    bool error,
    std::optional<bool> print_err = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> n_mat = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" void fortran_calc_bunch_params_slice(
    void* bunch /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type inout */,
    int& plane /* 0D_NOT_integer in */,
    double& slice_center /* 0D_NOT_real in */,
    double& slice_spread /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */,
    bool* is_time_coords /* 0D_NOT_logical in */,
    void* ele /* 0D_NOT_type in */);
void calc_bunch_params_slice(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    int plane,
    double slice_center,
    double slice_spread,
    bool err,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" void fortran_calc_bunch_params_z_slice(
    void* bunch /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type inout */,
    double* slice_bounds /* 1D_NOT_real in */,
    bool& err /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */,
    bool* is_time_coords /* 0D_NOT_logical in */,
    void* ele /* 0D_NOT_type in */);
void calc_bunch_params_z_slice(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    FixedArray1D<Real, 2> slice_bounds,
    bool err,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" void fortran_calc_bunch_sigma_matrix_etc(
    void* particle /* 1D_ALLOC_type in */,
    void* charge /* 1D_ALLOC_real in */,
    void* bunch_params /* 0D_NOT_type out */,
    bool* is_time_coords /* 0D_NOT_logical inout */,
    void* ele /* 0D_NOT_type inout */);
BunchParamsProxy calc_bunch_sigma_matrix_etc(
    CoordProxyAlloc1D& particle,
    RealAlloc1D& charge,
    optional_ref<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);

// Skipped unusable routine calc_density_derivative_complex:
// - Variable in sized array: density(:,:,:) 3D_NOT_real
// - Variable inout sized array: density_prime(:,:,:) 3D_NOT_complex
extern "C" void fortran_calc_emittances_and_twiss_from_sigma_matrix(
    double* sigma_mat /* 2D_NOT_real in */,
    void* bunch_params /* 0D_NOT_type out */,
    bool& error /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */,
    double* n_mat /* 2D_NOT_real out */);
struct CalcEmittancesAndTwissFromSigmaMatrix {
  BunchParamsProxy bunch_params;
  bool error;
  std::optional<FixedArray2D<Real, 6, 6>> n_mat;
};
Bmad::CalcEmittancesAndTwissFromSigmaMatrix
calc_emittances_and_twiss_from_sigma_matrix(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine calc_next_fringe_edge:
// - Untranslated type: fringe_field_info_struct (0D)
extern "C" void fortran_calc_spin_params(
    void* bunch /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type out */);
BunchParamsProxy calc_spin_params(BunchProxy& bunch);
extern "C" void fortran_calc_super_slave_key(
    void* lord1 /* 0D_NOT_type in */,
    void* lord2 /* 0D_NOT_type in */,
    void* slave /* 0D_NOT_type out */,
    bool* create_jumbo_slave /* 0D_NOT_logical in */);
EleProxy calc_super_slave_key(
    EleProxy& lord1,
    EleProxy& lord2,
    std::optional<bool> create_jumbo_slave = std::nullopt);
extern "C" void fortran_calc_wall_radius(
    void* v /* 1D_ALLOC_type in */,
    double& cos_ang /* 0D_NOT_real in */,
    double& sin_ang /* 0D_NOT_real in */,
    double& r_wall /* 0D_NOT_real out */,
    double& dr_dtheta /* 0D_NOT_real out */,
    int& ix_vertex /* 0D_NOT_integer out */);
struct CalcWallRadius {
  double r_wall;
  double dr_dtheta;
  int ix_vertex;
};
Bmad::CalcWallRadius calc_wall_radius(
    Wall3dVertexProxyAlloc1D& v,
    double cos_ang,
    double sin_ang);

// Skipped unusable routine calc_wiggler_g_params:
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)
extern "C" void fortran_calc_z_tune(void* branch /* 0D_NOT_type inout */);
void calc_z_tune(BranchProxy& branch);
extern "C" void fortran_canonical_to_angle_coords(
    void* orbit /* 0D_NOT_type inout */,
    const char* coord_type /* 0D_NOT_character in */);
void canonical_to_angle_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type = std::nullopt);

// Skipped unusable routine capillary_photon_hit_spot_calc:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_propagate_photon_a_step:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_reflect_photon:
// - Untranslated type: photon_track_struct (0D)

// Skipped unusable routine capillary_track_photon_to_wall:
// - Untranslated type: photon_track_struct (0D)
extern "C" void fortran_cbar_to_c(
    double* cbar_mat /* 2D_NOT_real in */,
    void* a /* 0D_NOT_type in */,
    void* b /* 0D_NOT_type in */,
    double* c_mat /* 2D_NOT_real out */);
FixedArray2D<Real, 2, 2> cbar_to_c(
    FixedArray2D<Real, 2, 2> cbar_mat,
    TwissProxy& a,
    TwissProxy& b);

// Skipped unusable routine ccfft3d:
// - Routine module (fft_interface_mod) in configuration skip list

// Skipped unusable routine ccfftam:
// - Routine in configuration skip list
extern "C" void fortran_check_aperture_limit(
    void* orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* param /* 0D_NOT_type inout */,
    void* old_orb /* 0D_NOT_type in */,
    bool* check_momentum /* 0D_NOT_logical in */);
void check_aperture_limit(
    CoordProxy& orb,
    EleProxy& ele,
    int particle_at,
    LatParamProxy& param,
    optional_ref<CoordProxy> old_orb = std::nullopt,
    std::optional<bool> check_momentum = std::nullopt);

// Skipped unusable routine check_aperture_limit_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_check_controller_controls(
    int& ele_key /* 0D_NOT_integer in */,
    void* contrl /* 1D_ALLOC_type in */,
    const char* name /* 0D_NOT_character in */,
    bool& err /* 0D_NOT_logical out */);
bool check_controller_controls(
    int ele_key,
    ControlProxyAlloc1D& contrl,
    std::string name);
extern "C" void fortran_check_for_superimpose_problem(
    void* branch /* 0D_NOT_type inout */,
    void* super_ele /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */,
    void* ref_ele /* 0D_NOT_type inout */,
    bool& wrap /* 0D_NOT_logical inout */);
void check_for_superimpose_problem(
    BranchProxy& branch,
    EleProxy& super_ele,
    bool& err_flag,
    optional_ref<EleProxy> ref_ele,
    bool& wrap);
extern "C" void fortran_check_if_s_in_bounds(
    void* branch /* 0D_NOT_type in */,
    double& s /* 0D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& translated_s /* 0D_NOT_real out */,
    bool* print_err /* 0D_NOT_logical in */);
struct CheckIfSInBounds {
  bool err_flag;
  double translated_s;
};
Bmad::CheckIfSInBounds check_if_s_in_bounds(
    BranchProxy& branch,
    double s,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_choose_quads_for_set_tune(
    void* branch /* 0D_NOT_type in */,
    void* dk1 /* 1D_ALLOC_real out */,
    void* eles /* 1D_ALLOC_type out */,
    const char* mask /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ChooseQuadsForSetTune {
  RealAlloc1D dk1;
  ElePointerProxyAlloc1D eles;
  bool err_flag;
};
Bmad::ChooseQuadsForSetTune choose_quads_for_set_tune(
    BranchProxy& branch,
    std::optional<std::string> mask = std::nullopt);
extern "C" void fortran_chrom_calc(
    void* lat /* 0D_NOT_type in */,
    double& delta_e /* 0D_NOT_real inout */,
    double& chrom_a /* 0D_NOT_real out */,
    double& chrom_b /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    double* pz /* 0D_NOT_real in */,
    void* low_E_lat /* 0D_NOT_type out */,
    void* high_E_lat /* 0D_NOT_type out */,
    void* low_E_orb /* 1D_ALLOC_type out */,
    void* high_E_orb /* 1D_ALLOC_type out */,
    int* ix_branch /* 0D_NOT_integer in */,
    void* orb0 /* 0D_NOT_type in */);
struct ChromCalc {
  double chrom_a;
  double chrom_b;
  bool err_flag;
  LatProxy low_E_lat;
  LatProxy high_E_lat;
  CoordProxyAlloc1D low_E_orb;
  CoordProxyAlloc1D high_E_orb;
};
Bmad::ChromCalc chrom_calc(
    LatProxy& lat,
    double& delta_e,
    std::optional<double> pz = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<CoordProxy> orb0 = std::nullopt);
extern "C" void fortran_chrom_tune(
    void* lat /* 0D_NOT_type inout */,
    double& delta_e /* 0D_NOT_real inout */,
    double& target_x /* 0D_NOT_real in */,
    double& target_y /* 0D_NOT_real in */,
    double& err_tol /* 0D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool chrom_tune(
    LatProxy& lat,
    double& delta_e,
    double target_x,
    double target_y,
    double err_tol);

// Skipped unusable routine cimp1:
// - Untranslated type: ibs_struct (0D)
extern "C" bool fortran_classical_radius(
    int& species /* 0D_NOT_integer in */,
    double& radius /* 0D_NOT_real inout */);
void classical_radius(int species, double& radius);
extern "C" void fortran_clear_lat_1turn_mats(void* lat /* 0D_NOT_type out */);
LatProxy clear_lat_1turn_mats();
extern "C" void fortran_clear_taylor_maps_from_elements(
    void* lat /* 0D_NOT_type inout */);
void clear_taylor_maps_from_elements(LatProxy& lat);
extern "C" void fortran_closed_orbit_calc(
    void* lat /* 0D_NOT_type inout */,
    void* closed_orb /* 1D_ALLOC_type inout */,
    int* i_dim /* 0D_NOT_integer in */,
    int* direction /* 0D_NOT_integer in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */);
bool closed_orbit_calc(
    LatProxy& lat,
    CoordProxyAlloc1D& closed_orb,
    std::optional<int> i_dim = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_closed_orbit_from_tracking(
    void* lat /* 0D_NOT_type in */,
    void* closed_orb /* 1D_ALLOC_type out */,
    int& i_dim /* 0D_NOT_integer in */,
    void* eps_rel /* 1D_ALLOC_real in */,
    void* eps_abs /* 1D_ALLOC_real in */,
    void* init_guess /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ClosedOrbitFromTracking {
  CoordProxyAlloc1D closed_orb;
  bool err_flag;
};
Bmad::ClosedOrbitFromTracking closed_orbit_from_tracking(
    LatProxy& lat,
    int i_dim,
    optional_ref<RealAlloc1D> eps_rel = std::nullopt,
    optional_ref<RealAlloc1D> eps_abs = std::nullopt,
    optional_ref<CoordProxy> init_guess = std::nullopt);
extern "C" bool fortran_cmplx_re_str(
    std::complex<double>& cmp /* 0D_NOT_complex inout */,
    const char* str_out /* 0D_NOT_character inout */);
void cmplx_re_str(std::complex<double>& cmp, std::string& str_out);
extern "C" void fortran_combine_consecutive_elements(
    void* lat /* 0D_NOT_type inout */,
    bool& error /* 0D_NOT_logical out */);
bool combine_consecutive_elements(LatProxy& lat);
extern "C" void fortran_complex_taylor_clean(
    void* complex_taylor /* 0D_NOT_type inout */);
void complex_taylor_clean(ComplexTaylorProxy& complex_taylor);

// Skipped unusable routine complex_taylor_equal_c_taylor:
// - Untranslated type: c_taylor (0D)
extern "C" void fortran_complex_taylor_equal_complex_taylor(
    void* complex_taylor1 /* 0D_NOT_type out */,
    void* complex_taylor2 /* 0D_NOT_type in */);
ComplexTaylorProxy complex_taylor_equal_complex_taylor(
    ComplexTaylorProxy& complex_taylor2);
extern "C" bool fortran_complex_taylor_exponent_index(
    int* expn /* 1D_NOT_integer in */,
    int& index /* 0D_NOT_integer out */);
int complex_taylor_exponent_index(FixedArray1D<Int, 6> expn);
extern "C" void fortran_complex_taylor_make_unit(
    void* complex_taylor /* 1D_ALLOC_type out */);
ComplexTaylorProxyAlloc1D complex_taylor_make_unit();
extern "C" void fortran_complex_taylor_to_mat6(
    void* a_complex_taylor /* 1D_NOT_type in */,
    void* r_in /* 1D_ALLOC_complex in */,
    std::complex<double>* vec0 /* 1D_NOT_complex out */,
    std::complex<double>* mat6 /* 2D_NOT_complex out */,
    void* r_out /* 1D_ALLOC_complex out */);
struct ComplexTaylorToMat6 {
  FixedArray1D<Complex, 6> vec0;
  FixedArray2D<Complex, 6, 6> mat6;
  ComplexAlloc1D r_out;
};
Bmad::ComplexTaylorToMat6 complex_taylor_to_mat6(
    FixedArray1D<ComplexTaylorProxy, 6> a_complex_taylor,
    ComplexAlloc1D& r_in);

// Skipped unusable routine complex_taylors_equal_c_taylors:
// - Untranslated type: c_taylor (1D)
extern "C" void fortran_complex_taylors_equal_complex_taylors(
    void* complex_taylor1 /* 1D_ALLOC_type out */,
    void* complex_taylor2 /* 1D_ALLOC_type in */);
ComplexTaylorProxyAlloc1D complex_taylors_equal_complex_taylors(
    ComplexTaylorProxyAlloc1D& complex_taylor2);
extern "C" void fortran_compute_slave_coupler(
    void* slave /* 0D_NOT_type inout */);
void compute_slave_coupler(EleProxy& slave);

// Skipped unusable routine compute_super_lord_s:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_concat_ele_taylor(
    void* orb_taylor /* 1D_ALLOC_type in */,
    void* ele /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */,
    void* spin_taylor /* 1D_ALLOC_type in */);
void concat_ele_taylor(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    bool err_flag,
    optional_ref<TaylorProxyAlloc1D> spin_taylor = std::nullopt);

// Skipped unusable routine concat_real_8:
// - Untranslated type: real_8 (1D)
// - Untranslated type: real_8 (1D)
// - Untranslated type: real_8 (1D)
extern "C" void fortran_concat_taylor(
    void* taylor1 /* 1D_ALLOC_type in */,
    void* taylor2 /* 1D_ALLOC_type in */,
    void* taylor3 /* 1D_ALLOC_type in */);
void concat_taylor(
    TaylorProxyAlloc1D& taylor1,
    TaylorProxyAlloc1D& taylor2,
    TaylorProxyAlloc1D& taylor3);
extern "C" void fortran_concat_transfer_mat(
    double* mat_1 /* 2D_NOT_real in */,
    double* vec_1 /* 1D_NOT_real inout */,
    double* mat_0 /* 2D_NOT_real in */,
    double* vec_0 /* 1D_NOT_real inout */,
    double* mat_out /* 2D_NOT_real out */,
    double* vec_out /* 1D_NOT_real inout */);
FixedArray2D<Real, 6, 6> concat_transfer_mat(
    FixedArray2D<Real, 6, 6> mat_1,
    FixedArray1D<Real, 6> vec_1,
    FixedArray2D<Real, 6, 6> mat_0,
    FixedArray1D<Real, 6> vec_0,
    FixedArray1D<Real, 6> vec_out);
extern "C" void fortran_control_bookkeeper(
    void* lat /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    bool* err_flag /* 0D_NOT_logical in */);
void control_bookkeeper(
    LatProxy& lat,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<bool> err_flag = std::nullopt);

// Skipped unusable routine conv3d:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_convert_bend_exact_multipole(
    double& g /* 0D_NOT_real in */,
    int& out_type /* 0D_NOT_integer in */,
    double* an /* 1D_NOT_real inout */,
    double* bn /* 1D_NOT_real inout */);
void convert_bend_exact_multipole(
    double g,
    int out_type,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> an,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> bn);
extern "C" void fortran_convert_coords(
    const char* in_type_str /* 0D_NOT_character in */,
    void* coord_in /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    const char* out_type_str /* 0D_NOT_character out */,
    void* coord_out /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ConvertCoords {
  std::string out_type_str;
  CoordProxy coord_out;
  bool err_flag;
};
Bmad::ConvertCoords convert_coords(
    std::string in_type_str,
    CoordProxy& coord_in,
    EleProxy& ele);
extern "C" void fortran_convert_field_ele_to_lab(
    void* ele /* 0D_NOT_type in */,
    double& s_here /* 0D_NOT_real in */,
    bool& forward_transform /* 0D_NOT_logical in */,
    void* field /* 0D_NOT_type out */,
    bool* calc_dfield /* 0D_NOT_logical in */,
    bool* calc_potential /* 0D_NOT_logical in */);
EmFieldProxy convert_field_ele_to_lab(
    EleProxy& ele,
    double s_here,
    bool forward_transform,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" void fortran_convert_local_cartesian_to_local_curvilinear(
    double& x /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& g /* 0D_NOT_real inout */,
    double& xout /* 0D_NOT_real inout */,
    double& sout /* 0D_NOT_real inout */);
void convert_local_cartesian_to_local_curvilinear(
    double& x,
    double& z,
    double& g,
    double& xout,
    double& sout);
extern "C" void fortran_convert_local_curvilinear_to_local_cartesian(
    double& x /* 0D_NOT_real inout */,
    double& s /* 0D_NOT_real inout */,
    double& g /* 0D_NOT_real inout */,
    double& xout /* 0D_NOT_real inout */,
    double& zout /* 0D_NOT_real inout */);
void convert_local_curvilinear_to_local_cartesian(
    double& x,
    double& s,
    double& g,
    double& xout,
    double& zout);
extern "C" void fortran_convert_particle_coordinates_s_to_t(
    void* particle /* 0D_NOT_type inout */,
    double& s_body /* 0D_NOT_real in */,
    int& orientation /* 0D_NOT_integer in */);
void convert_particle_coordinates_s_to_t(
    CoordProxy& particle,
    double s_body,
    int orientation);
extern "C" void fortran_convert_particle_coordinates_t_to_s(
    void* particle /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    double& s_body /* 0D_NOT_real out */,
    bool* use_downstream_p0c /* 0D_NOT_logical in */);
double convert_particle_coordinates_t_to_s(
    CoordProxy& particle,
    EleProxy& ele,
    std::optional<bool> use_downstream_p0c = std::nullopt);
extern "C" void fortran_convert_pc_to(
    double& pc /* 0D_NOT_real in */,
    int& particle /* 0D_NOT_integer in */,
    double& E_tot /* 0D_NOT_real out */,
    double& gamma /* 0D_NOT_real out */,
    double& kinetic /* 0D_NOT_real out */,
    double& beta /* 0D_NOT_real out */,
    double& brho /* 0D_NOT_real out */,
    double& beta1 /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
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
    double& E_tot /* 0D_NOT_real in */,
    int& particle /* 0D_NOT_integer in */,
    double& gamma /* 0D_NOT_real out */,
    double& kinetic /* 0D_NOT_real out */,
    double& beta /* 0D_NOT_real out */,
    double& pc /* 0D_NOT_real out */,
    double& brho /* 0D_NOT_real out */,
    double& beta1 /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */);
struct ConvertTotalEnergyTo {
  double gamma;
  double kinetic;
  double beta;
  double pc;
  double brho;
  double beta1;
  bool err_flag;
};
Bmad::ConvertTotalEnergyTo convert_total_energy_to(
    double E_tot,
    int particle,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_converter_distribution_parser(
    void* ele /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character out */,
    bool& delim_found /* 0D_NOT_logical out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ConverterDistributionParser {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
Bmad::ConverterDistributionParser converter_distribution_parser(EleProxy& ele);
extern "C" void fortran_coord_equal_coord(
    void* coord1 /* 0D_NOT_type out */,
    void* coord2 /* 0D_NOT_type in */);
CoordProxy coord_equal_coord(CoordProxy& coord2);
extern "C" bool fortran_coord_state_name(
    int& coord_state /* 0D_NOT_integer in */,
    bool* one_word /* 0D_NOT_logical inout */,
    const char* state_str /* 0D_NOT_character out */);
std::string coord_state_name(
    int coord_state,
    optional_ref<bool> one_word = std::nullopt);
extern "C" bool fortran_coords_body_to_local(
    void* body_position /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double* w_mat /* 2D_NOT_real in */,
    bool* calculate_angles /* 0D_NOT_logical in */,
    void* local_position /* 0D_NOT_type inout */);
void coords_body_to_local(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_body_to_rel_exit(
    void* body_position /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double* w_mat /* 2D_NOT_real in */,
    bool* calculate_angles /* 0D_NOT_logical in */,
    void* rel_exit /* 0D_NOT_type inout */);
void coords_body_to_rel_exit(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& rel_exit);
extern "C" bool fortran_coords_curvilinear_to_floor(
    double* xys /* 1D_NOT_real in */,
    void* branch /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* global /* 0D_NOT_type inout */);
bool coords_curvilinear_to_floor(
    FixedArray1D<Real, 3> xys,
    BranchProxy& branch,
    FloorPositionProxy& global);
extern "C" bool fortran_coords_floor_to_curvilinear(
    void* floor_coords /* 0D_NOT_type in */,
    void* ele0 /* 0D_NOT_type in */,
    void* ele1 /* 0D_PTR_type out */,
    int& status /* 0D_NOT_integer out */,
    double* w_mat /* 2D_NOT_real out */,
    void* local_coords /* 0D_NOT_type inout */);
struct CoordsFloorToCurvilinear {
  EleProxy ele1;
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
};
Bmad::CoordsFloorToCurvilinear coords_floor_to_curvilinear(
    FloorPositionProxy& floor_coords,
    EleProxy& ele0,
    FloorPositionProxy& local_coords);
extern "C" bool fortran_coords_floor_to_local_curvilinear(
    void* global_position /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    int& status /* 0D_NOT_integer out */,
    double* w_mat /* 2D_NOT_real out */,
    int* relative_to /* 0D_NOT_integer in */,
    void* local_position /* 0D_NOT_type inout */);
struct CoordsFloorToLocalCurvilinear {
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
};
Bmad::CoordsFloorToLocalCurvilinear coords_floor_to_local_curvilinear(
    FloorPositionProxy& global_position,
    EleProxy& ele,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_floor_to_relative(
    void* floor0 /* 0D_NOT_type in */,
    void* global_position /* 0D_NOT_type in */,
    bool* calculate_angles /* 0D_NOT_logical in */,
    bool* is_delta_position /* 0D_NOT_logical in */,
    void* local_position /* 0D_NOT_type inout */);
void coords_floor_to_relative(
    FloorPositionProxy& floor0,
    FloorPositionProxy& global_position,
    std::optional<bool> calculate_angles,
    std::optional<bool> is_delta_position,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_local_curvilinear_to_body(
    void* local_position /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double* w_mat /* 2D_NOT_real in */,
    bool* calculate_angles /* 0D_NOT_logical in */,
    void* body_position /* 0D_NOT_type inout */);
void coords_local_curvilinear_to_body(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& body_position);
extern "C" bool fortran_coords_local_curvilinear_to_floor(
    void* local_position /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    bool* in_body_frame /* 0D_NOT_logical in */,
    double* w_mat /* 2D_NOT_real out */,
    bool* calculate_angles /* 0D_NOT_logical in */,
    int* relative_to /* 0D_NOT_integer in */,
    void* global_position /* 0D_NOT_type inout */);
FixedArray2D<Real, 3, 3> coords_local_curvilinear_to_floor(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<bool> in_body_frame,
    std::optional<bool> calculate_angles,
    std::optional<int> relative_to,
    FloorPositionProxy& global_position);
extern "C" bool fortran_coords_relative_to_floor(
    void* floor0 /* 0D_NOT_type in */,
    double* dr /* 1D_NOT_real in */,
    double* theta /* 0D_NOT_real inout */,
    double* phi /* 0D_NOT_real inout */,
    double* psi /* 0D_NOT_real inout */,
    void* floor1 /* 0D_NOT_type inout */);
void coords_relative_to_floor(
    FloorPositionProxy& floor0,
    FixedArray1D<Real, 3> dr,
    optional_ref<double> theta,
    optional_ref<double> phi,
    optional_ref<double> psi,
    FloorPositionProxy& floor1);

// Skipped unusable routine cos_phi:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" bool fortran_coulombfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void coulombfun(double& u, double& v, double& w, double& gam, double& res);
extern "C" void fortran_create_concatenated_wall3d(
    void* lat /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical inout */);
void create_concatenated_wall3d(LatProxy& lat, bool& err);
extern "C" void fortran_create_element_slice(
    void* sliced_ele /* 0D_NOT_type out */,
    void* ele_in /* 0D_NOT_type in */,
    double& l_slice /* 0D_NOT_real in */,
    double& offset /* 0D_NOT_real in */,
    void* param /* 0D_NOT_type in */,
    bool& include_upstream_end /* 0D_NOT_logical in */,
    bool& include_downstream_end /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* old_slice /* 0D_NOT_type in */,
    void* orb_in /* 0D_NOT_type in */);
struct CreateElementSlice {
  EleProxy sliced_ele;
  bool err_flag;
};
Bmad::CreateElementSlice create_element_slice(
    EleProxy& ele_in,
    double l_slice,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end,
    optional_ref<EleProxy> old_slice = std::nullopt,
    optional_ref<CoordProxy> orb_in = std::nullopt);

// Skipped unusable routine create_feedback:
// - Variable-sized in character array: input(:) 1D_ALLOC_character
// - Variable-sized in character array: output(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_create_field_overlap(
    void* lat /* 0D_NOT_type inout */,
    const char* lord_name /* 0D_NOT_character in */,
    const char* slave_name /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool create_field_overlap(
    LatProxy& lat,
    std::string lord_name,
    std::string slave_name);
extern "C" void fortran_create_girder(
    void* lat /* 0D_NOT_type inout */,
    int& ix_girder /* 0D_NOT_integer in */,
    void* contrl /* 1D_ALLOC_type in */,
    void* girder_info /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical inout */);
void create_girder(
    LatProxy& lat,
    int ix_girder,
    ControlProxyAlloc1D& contrl,
    EleProxy& girder_info,
    bool& err_flag);
extern "C" void fortran_create_group(
    void* lord /* 0D_NOT_type inout */,
    void* contrl /* 1D_ALLOC_type in */,
    bool& err /* 0D_NOT_logical in */);
void create_group(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);
extern "C" void fortran_create_lat_ele_nametable(
    void* lat /* 0D_NOT_type in */,
    void* nametable /* 0D_NOT_type in */);
void create_lat_ele_nametable(LatProxy& lat, NametableProxy& nametable);
extern "C" void fortran_create_overlay(
    void* lord /* 0D_NOT_type inout */,
    void* contrl /* 1D_ALLOC_type in */,
    bool& err /* 0D_NOT_logical in */);
void create_overlay(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);
extern "C" void fortran_create_planar_wiggler_model(
    void* wiggler_in /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */);
struct CreatePlanarWigglerModel {
  LatProxy lat;
  bool err_flag;
};
Bmad::CreatePlanarWigglerModel create_planar_wiggler_model(
    EleProxy& wiggler_in,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_create_ramper(
    void* lord /* 0D_NOT_type inout */,
    void* contrl /* 1D_ALLOC_type in */,
    bool& err /* 0D_NOT_logical in */);
void create_ramper(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);
extern "C" void fortran_create_sol_quad_model(
    void* sol_quad /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */);
void create_sol_quad_model(EleProxy& sol_quad, LatProxy& lat);
extern "C" void fortran_create_unique_ele_names(
    void* lat /* 0D_NOT_type inout */,
    int& key /* 0D_NOT_integer in */,
    const char* suffix /* 0D_NOT_character in */);
void create_unique_ele_names(LatProxy& lat, int key, std::string suffix);
extern "C" void fortran_create_wiggler_cartesian_map(
    void* ele /* 0D_NOT_type in */,
    void* cart_map /* 0D_NOT_type out */);
CartesianMapProxy create_wiggler_cartesian_map(EleProxy& ele);
extern "C" void fortran_crystal_attribute_bookkeeper(
    void* ele /* 0D_NOT_type in */);
void crystal_attribute_bookkeeper(EleProxy& ele);

// Skipped unusable routine crystal_diffraction_field_calc:
// - Untranslated type: crystal_param_struct (0D)
extern "C" void fortran_crystal_h_misalign(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    double* h_vec /* 1D_NOT_real inout */);
void crystal_h_misalign(
    EleProxy& ele,
    CoordProxy& orbit,
    FixedArray1D<Real, 3> h_vec);
extern "C" void fortran_crystal_type_to_crystal_params(
    void* ele /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool crystal_type_to_crystal_params(EleProxy& ele);

// Skipped unusable routine csr3d_steady_state_solver:
// - Variable in sized array: density(:,:,:) 3D_NOT_real
// - Variable inout sized array: wake(:,:,:,:) 4D_NOT_real

// Skipped unusable routine csr_and_sc_apply_kicks:
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine csr_bin_kicks:
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine csr_bin_particles:
// - Untranslated type: csr_struct (0D)
extern "C" bool fortran_custom_attribute_ubound_index(
    int& ele_class /* 0D_NOT_integer in */,
    int& ix_ubound /* 0D_NOT_integer out */);
int custom_attribute_ubound_index(int ele_class);

// Skipped unusable routine custom_ele_attrib_name_list:
// - Variable-sized out character array: name_list(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine damap_equal_bmad_taylor:
// - Untranslated type: damap (0D)
extern "C" bool fortran_damping_matrix_d(
    double& gamma /* 0D_NOT_real inout */,
    double& g_tot /* 0D_NOT_real inout */,
    double& B0 /* 0D_NOT_real inout */,
    double& B1 /* 0D_NOT_real inout */,
    double& delta /* 0D_NOT_real inout */,
    int& species /* 0D_NOT_integer inout */,
    double* mat /* 2D_NOT_real inout */);
void damping_matrix_d(
    double& gamma,
    double& g_tot,
    double& B0,
    double& B1,
    double& delta,
    int& species,
    FixedArray2D<Real, 6, 6> mat);

// Skipped unusable routine deallocate_ele_array_pointers:
// - Routine in configuration skip list
extern "C" void fortran_deallocate_ele_pointers(
    void* ele /* 0D_NOT_type inout */,
    bool* nullify_only /* 0D_NOT_logical in */,
    bool* nullify_branch /* 0D_NOT_logical in */,
    bool* dealloc_poles /* 0D_NOT_logical in */);
void deallocate_ele_pointers(
    EleProxy& ele,
    std::optional<bool> nullify_only = std::nullopt,
    std::optional<bool> nullify_branch = std::nullopt,
    std::optional<bool> dealloc_poles = std::nullopt);
extern "C" void fortran_deallocate_expression_tree(
    void* tree /* 0D_NOT_type inout */);
void deallocate_expression_tree(ExpressionTreeProxy& tree);
extern "C" void fortran_deallocate_lat_pointers(
    void* lat /* 0D_NOT_type inout */);
void deallocate_lat_pointers(LatProxy& lat);
extern "C" bool fortran_default_tracking_species(
    void* param /* 0D_NOT_type in */,
    int& species /* 0D_NOT_integer inout */);
void default_tracking_species(LatParamProxy& param, int& species);

// Skipped unusable routine deposit_particles:
// - Untranslated type: mesh3d_struct (0D)
extern "C" bool fortran_detector_pixel_pt(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    int* ix_pix /* 1D_NOT_integer out */);
FixedArray1D<Int, 2> detector_pixel_pt(CoordProxy& orbit, EleProxy& ele);
extern "C" bool fortran_diffraction_plate_or_mask_hit_spot(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    int& ix_section /* 0D_NOT_integer inout */);
void diffraction_plate_or_mask_hit_spot(
    EleProxy& ele,
    CoordProxy& orbit,
    int& ix_section);
extern "C" bool fortran_diffusion_matrix_b(
    double& gamma /* 0D_NOT_real inout */,
    double& g_tot /* 0D_NOT_real inout */,
    int& species /* 0D_NOT_integer inout */,
    double* mat /* 2D_NOT_real inout */);
void diffusion_matrix_b(
    double& gamma,
    double& g_tot,
    int& species,
    FixedArray2D<Real, 6, 6> mat);
extern "C" bool fortran_distance_to_aperture(
    void* orbit /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    bool& no_aperture_here /* 0D_NOT_logical out */,
    double& dist /* 0D_NOT_real inout */);
bool distance_to_aperture(
    CoordProxy& orbit,
    int particle_at,
    EleProxy& ele,
    double& dist);

// Skipped unusable routine distance_to_aperture_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_do_mode_flip(
    void* ele /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool do_mode_flip(EleProxy& ele);
extern "C" bool fortran_dpc_given_de(
    double& pc_old /* 0D_NOT_real inout */,
    double& mass /* 0D_NOT_real inout */,
    double& dE /* 0D_NOT_real inout */,
    double& dpc /* 0D_NOT_real inout */);
void dpc_given_de(double& pc_old, double& mass, double& dE, double& dpc);
extern "C" void fortran_drift_and_pipe_track_methods_adjustment(
    void* lat /* 0D_NOT_type inout */);
void drift_and_pipe_track_methods_adjustment(LatProxy& lat);
extern "C" void fortran_drift_multipass_name_correction(
    void* lat /* 0D_NOT_type inout */);
void drift_multipass_name_correction(LatProxy& lat);
extern "C" void fortran_drift_orbit_time(
    void* orbit /* 0D_NOT_type inout */,
    double& beta0 /* 0D_NOT_real in */,
    double* delta_s /* 0D_NOT_real in */,
    double* delta_t /* 0D_NOT_real in */);
void drift_orbit_time(
    CoordProxy& orbit,
    double beta0,
    std::optional<double> delta_s = std::nullopt,
    std::optional<double> delta_t = std::nullopt);
extern "C" void fortran_drift_particle_to_s(
    void* p /* 0D_NOT_type inout */,
    double& s /* 0D_NOT_real in */,
    void* branch /* 0D_NOT_type in */);
void drift_particle_to_s(CoordProxy& p, double s, BranchProxy& branch);
extern "C" void fortran_drift_particle_to_t(
    void* p /* 0D_NOT_type inout */,
    double& t /* 0D_NOT_real in */,
    void* branch /* 0D_NOT_type in */);
void drift_particle_to_t(CoordProxy& p, double t, BranchProxy& branch);
extern "C" bool fortran_dspline_len(
    double& s_chord0 /* 0D_NOT_real in */,
    double& s_chord1 /* 0D_NOT_real in */,
    void* spline /* 0D_NOT_type in */,
    double* dtheta_ref /* 0D_NOT_real in */,
    double& dlen /* 0D_NOT_real out */);
double dspline_len(
    double s_chord0,
    double s_chord1,
    SplineProxy& spline,
    std::optional<double> dtheta_ref = std::nullopt);
extern "C" void fortran_dynamic_aperture_point(
    void* branch /* 0D_NOT_type in */,
    void* ele0 /* 0D_NOT_type in */,
    void* orb0 /* 0D_NOT_type in */,
    double& theta_xy /* 0D_NOT_real in */,
    void* ap_param /* 0D_NOT_type in */,
    void* ap_point /* 0D_NOT_type out */,
    bool* check_xy_init /* 0D_NOT_logical in */);
AperturePointProxy dynamic_aperture_point(
    BranchProxy& branch,
    EleProxy& ele0,
    CoordProxy& orb0,
    double theta_xy,
    ApertureParamProxy& ap_param,
    std::optional<bool> check_xy_init = std::nullopt);
extern "C" void fortran_dynamic_aperture_scan(
    void* aperture_scan /* 1D_ALLOC_type out */,
    void* aperture_param /* 0D_NOT_type in */,
    void* pz_start /* 1D_ALLOC_real in */,
    void* lat /* 0D_NOT_type in */,
    bool* print_timing /* 0D_NOT_logical in */);
ApertureScanProxyAlloc1D dynamic_aperture_scan(
    ApertureParamProxy& aperture_param,
    RealAlloc1D& pz_start,
    LatProxy& lat,
    std::optional<bool> print_timing = std::nullopt);
extern "C" bool fortran_e_accel_field(
    void* ele /* 0D_NOT_type in */,
    int& voltage_or_gradient /* 0D_NOT_integer in */,
    bool* bmad_standard_tracking /* 0D_NOT_logical in */,
    double& field /* 0D_NOT_real inout */);
void e_accel_field(
    EleProxy& ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking,
    double& field);
extern "C" bool fortran_e_crit_photon(
    double& gamma /* 0D_NOT_real in */,
    double& g_bend /* 0D_NOT_real in */,
    double& E_crit /* 0D_NOT_real out */);
double e_crit_photon(double gamma, double g_bend);
extern "C" void fortran_eigen_decomp_6mat(
    double* mat /* 2D_NOT_real in */,
    std::complex<double>* eval /* 1D_NOT_complex out */,
    std::complex<double>* evec /* 2D_NOT_complex out */,
    bool& err_flag /* 0D_NOT_logical out */,
    double* tunes /* 1D_NOT_real out */);
struct EigenDecomp6mat {
  FixedArray1D<Complex, 6> eval;
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
  FixedArray1D<Real, 3> tunes;
};
Bmad::EigenDecomp6mat eigen_decomp_6mat(FixedArray2D<Real, 6, 6> mat);

// Skipped unusable routine eigensys:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_ele_compute_ref_energy_and_time(
    void* ele0 /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void ele_compute_ref_energy_and_time(
    EleProxy& ele0,
    EleProxy& ele,
    LatParamProxy& param,
    bool err_flag);
extern "C" void fortran_ele_equal_ele(
    void* ele_out /* 0D_NOT_type out */,
    void* ele_in /* 0D_NOT_type in */);
EleProxy ele_equal_ele(EleProxy& ele_in);
extern "C" void fortran_ele_equals_ele(
    void* ele_out /* 0D_NOT_type out */,
    void* ele_in /* 0D_NOT_type in */,
    bool& update_nametable /* 0D_NOT_logical in */);
EleProxy ele_equals_ele(EleProxy& ele_in, bool update_nametable);
extern "C" void fortran_ele_finalizer(void* ele /* 0D_NOT_type inout */);
void ele_finalizer(EleProxy& ele);
extern "C" bool fortran_ele_full_name(
    void* ele /* 0D_NOT_type in */,
    const char* template_ /* 0D_NOT_character in */,
    const char* str /* 0D_ALLOC_character inout */);
void ele_full_name(
    EleProxy& ele,
    std::optional<std::string> template_,
    std::string& str);
extern "C" void fortran_ele_geometry(
    void* floor_start /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* floor_end /* 0D_NOT_type out */,
    double* len_scale /* 0D_NOT_real in */,
    bool* ignore_patch_err /* 0D_NOT_logical in */);
FloorPositionProxy ele_geometry(
    FloorPositionProxy& floor_start,
    EleProxy& ele,
    std::optional<double> len_scale = std::nullopt,
    std::optional<bool> ignore_patch_err = std::nullopt);

// Skipped unusable routine ele_geometry_hook_def:
// - Routine in configuration skip list
extern "C" bool fortran_ele_geometry_with_misalignments(
    void* ele /* 0D_NOT_type in */,
    double* len_scale /* 0D_NOT_real in */,
    void* floor /* 0D_NOT_type inout */);
void ele_geometry_with_misalignments(
    EleProxy& ele,
    std::optional<double> len_scale,
    FloorPositionProxy& floor);
extern "C" bool fortran_ele_has_constant_ds_dt_ref(
    void* ele /* 0D_NOT_type in */,
    bool& is_const /* 0D_NOT_logical inout */);
void ele_has_constant_ds_dt_ref(EleProxy& ele, bool& is_const);
extern "C" bool fortran_ele_has_nonzero_kick(
    void* ele /* 0D_NOT_type out */,
    bool& has_kick /* 0D_NOT_logical inout */);
EleProxy ele_has_nonzero_kick(bool& has_kick);
extern "C" bool fortran_ele_has_nonzero_offset(
    void* ele /* 0D_NOT_type inout */,
    bool& has_offset /* 0D_NOT_logical inout */);
void ele_has_nonzero_offset(EleProxy& ele, bool& has_offset);
extern "C" bool fortran_ele_is_monitor(
    void* ele /* 0D_NOT_type in */,
    bool* print_warning /* 0D_NOT_logical in */,
    bool& is_monitor /* 0D_NOT_logical out */);
bool ele_is_monitor(
    EleProxy& ele,
    std::optional<bool> print_warning = std::nullopt);
extern "C" bool fortran_ele_loc(
    void* ele /* 0D_NOT_type in */,
    void* loc /* 0D_NOT_type inout */);
void ele_loc(EleProxy& ele, LatEleLocProxy& loc);
extern "C" bool fortran_ele_loc_name(
    void* ele /* 0D_NOT_type in */,
    bool* show_branch0 /* 0D_NOT_logical in */,
    const char* parens /* 0D_NOT_character in */,
    const char* str /* 0D_NOT_character inout */);
void ele_loc_name(
    EleProxy& ele,
    std::optional<bool> show_branch0,
    std::optional<std::string> parens,
    std::string& str);
extern "C" void fortran_ele_misalignment_l_s_calc(
    void* ele /* 0D_NOT_type in */,
    double* L_mis /* 1D_NOT_real out */,
    double* S_mis /* 2D_NOT_real out */);
struct EleMisalignmentLSCalc {
  FixedArray1D<Real, 3> L_mis;
  FixedArray2D<Real, 3, 3> S_mis;
};
Bmad::EleMisalignmentLSCalc ele_misalignment_l_s_calc(EleProxy& ele);
extern "C" bool fortran_ele_nametable_index(
    void* ele /* 0D_NOT_type in */,
    int& ix_nt /* 0D_NOT_integer inout */);
void ele_nametable_index(EleProxy& ele, int& ix_nt);
extern "C" void fortran_ele_order_calc(
    void* lat /* 0D_NOT_type in */,
    void* order /* 0D_NOT_type out */);
LatEleOrderProxy ele_order_calc(LatProxy& lat);
extern "C" void fortran_ele_reference_energy_correction(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    int& particle_at /* 0D_NOT_integer in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void ele_reference_energy_correction(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_ele_rf_step_index(
    double& E_ref /* 0D_NOT_real in */,
    double& s_rel /* 0D_NOT_real in */,
    void* ele /* 0D_NOT_type in */,
    int& ix_step /* 0D_NOT_integer inout */);
void ele_rf_step_index(double E_ref, double s_rel, EleProxy& ele, int& ix_step);

// Skipped unusable routine ele_to_fibre:
// - Untranslated type: fibre (0D)

// Skipped unusable routine ele_to_fibre_hook_def:
// - Untranslated type: fibre (0D)
extern "C" void fortran_ele_to_ptc_magnetic_bn_an(
    void* ele /* 0D_NOT_type in */,
    void* bn /* 1D_ALLOC_real out */,
    void* an /* 1D_ALLOC_real out */,
    int& n_max /* 0D_NOT_integer out */);
struct EleToPtcMagneticBnAn {
  RealAlloc1D bn;
  RealAlloc1D an;
  int n_max;
};
Bmad::EleToPtcMagneticBnAn ele_to_ptc_magnetic_bn_an(EleProxy& ele);
extern "C" void fortran_ele_to_spin_taylor(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* orb0 /* 0D_NOT_type in */);
void ele_to_spin_taylor(EleProxy& ele, LatParamProxy& param, CoordProxy& orb0);
extern "C" void fortran_ele_to_taylor(
    void* ele /* 0D_NOT_type in */,
    void* orb0 /* 0D_NOT_type in */,
    bool* taylor_map_includes_offsets /* 0D_NOT_logical in */,
    bool* include_damping /* 0D_NOT_logical in */,
    void* orbital_taylor /* 1D_NOT_type out */,
    void* spin_taylor /* 1D_NOT_type out */);
struct EleToTaylor {
  TaylorProxyArray1D orbital_taylor;
  TaylorProxyArray1D spin_taylor;
};
Bmad::EleToTaylor ele_to_taylor(
    EleProxy& ele,
    optional_ref<CoordProxy> orb0 = std::nullopt,
    std::optional<bool> taylor_map_includes_offsets = std::nullopt,
    std::optional<bool> include_damping = std::nullopt);
extern "C" bool fortran_ele_unique_name(
    void* ele /* 0D_NOT_type in */,
    void* order /* 0D_NOT_type in */,
    const char* unique_name /* 0D_NOT_character inout */);
void ele_unique_name(
    EleProxy& ele,
    LatEleOrderProxy& order,
    std::string& unique_name);
extern "C" bool fortran_ele_value_has_changed(
    void* ele /* 0D_NOT_type inout */,
    void* list /* 1D_ALLOC_integer in */,
    void* abs_tol /* 1D_ALLOC_real in */,
    bool& set_old /* 0D_NOT_logical in */,
    bool& has_changed /* 0D_NOT_logical inout */);
void ele_value_has_changed(
    EleProxy& ele,
    IntAlloc1D& list,
    RealAlloc1D& abs_tol,
    bool set_old,
    bool& has_changed);
extern "C" void fortran_ele_vec_equal_ele_vec(
    void* ele1 /* 1D_ALLOC_type out */,
    void* ele2 /* 1D_ALLOC_type in */);
EleProxyAlloc1D ele_vec_equal_ele_vec(EleProxyAlloc1D& ele2);
extern "C" void fortran_elec_multipole_field(
    double& a /* 0D_NOT_real in */,
    double& b /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer in */,
    void* coord /* 0D_NOT_type in */,
    double& Ex /* 0D_NOT_real out */,
    double& Ey /* 0D_NOT_real out */,
    double* dE /* 2D_NOT_real out */,
    bool& compute_dE /* 0D_NOT_logical out */);
struct ElecMultipoleField {
  double Ex;
  double Ey;
  std::optional<FixedArray2D<Real, 2, 2>> dE;
  bool compute_dE;
};
Bmad::ElecMultipoleField elec_multipole_field(
    double a,
    double b,
    int n,
    CoordProxy& coord);
extern "C" void fortran_element_slice_iterator(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& i_slice /* 0D_NOT_integer in */,
    int& n_slice_tot /* 0D_NOT_integer in */,
    void* sliced_ele /* 0D_NOT_type inout */,
    double* s_start /* 0D_NOT_real in */,
    double* s_end /* 0D_NOT_real in */);
void element_slice_iterator(
    EleProxy& ele,
    LatParamProxy& param,
    int i_slice,
    int n_slice_tot,
    EleProxy& sliced_ele,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt);
extern "C" void fortran_ellipinc_test();
void ellipinc_test();
extern "C" void fortran_em_field_calc(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& s_pos /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type in */,
    bool& local_ref_frame /* 0D_NOT_logical in */,
    void* field /* 0D_NOT_type out */,
    bool* calc_dfield /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* calc_potential /* 0D_NOT_logical in */,
    bool* use_overlap /* 0D_NOT_logical in */,
    bool* grid_allow_s_out_of_bounds /* 0D_NOT_logical in */,
    double* rf_time /* 0D_NOT_real in */,
    void* used_eles /* 1D_ALLOC_type in */,
    bool* print_err /* 0D_NOT_logical in */,
    void* original_ele /* 0D_NOT_type in */);
struct EmFieldCalc {
  EmFieldProxy field;
  bool err_flag;
};
Bmad::EmFieldCalc em_field_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_pos,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt,
    std::optional<bool> use_overlap = std::nullopt,
    std::optional<bool> grid_allow_s_out_of_bounds = std::nullopt,
    std::optional<double> rf_time = std::nullopt,
    optional_ref<ElePointerProxyAlloc1D> used_eles = std::nullopt,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EleProxy> original_ele = std::nullopt);

// Skipped unusable routine em_field_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_em_field_derivatives(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    double& s_pos /* 0D_NOT_real inout */,
    void* orbit /* 0D_NOT_type inout */,
    bool& local_ref_frame /* 0D_NOT_logical inout */,
    void* dfield /* 0D_NOT_type out */,
    bool* grid_allow_s_out_of_bounds /* 0D_NOT_logical inout */,
    double* rf_time /* 0D_NOT_real inout */);
EmFieldProxy em_field_derivatives(
    EleProxy& ele,
    LatParamProxy& param,
    double& s_pos,
    CoordProxy& orbit,
    bool& local_ref_frame,
    optional_ref<bool> grid_allow_s_out_of_bounds = std::nullopt,
    optional_ref<double> rf_time = std::nullopt);
extern "C" void fortran_em_field_kick_vector_time(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& rf_time /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type in */,
    double* dvec_dt /* 1D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */,
    void* extra_field /* 0D_NOT_type in */);
FixedArray1D<Real, 10> em_field_kick_vector_time(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orbit,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" bool fortran_em_field_plus_em_field(
    void* field1 /* 0D_NOT_type in */,
    void* field2 /* 0D_NOT_type in */,
    void* field_tot /* 0D_NOT_type out */);
EmFieldProxy em_field_plus_em_field(EmFieldProxy& field1, EmFieldProxy& field2);
extern "C" void fortran_em_taylor_equal_em_taylor(
    void* em_taylor1 /* 0D_NOT_type out */,
    void* em_taylor2 /* 0D_NOT_type in */);
EmTaylorProxy em_taylor_equal_em_taylor(EmTaylorProxy& em_taylor2);
extern "C" void fortran_em_taylors_equal_em_taylors(
    void* em_taylor1 /* 1D_ALLOC_type out */,
    void* em_taylor2 /* 1D_ALLOC_type in */);
EmTaylorProxyAlloc1D em_taylors_equal_em_taylors(
    EmTaylorProxyAlloc1D& em_taylor2);
extern "C" void fortran_emit_6d(
    void* ele_ref /* 0D_NOT_type in */,
    bool& include_opening_angle /* 0D_NOT_logical in */,
    void* mode /* 0D_NOT_type out */,
    double* sigma_mat /* 2D_NOT_real out */,
    void* closed_orbit /* 1D_ALLOC_type in */,
    void* rad_int_by_ele /* 0D_NOT_type out */);
struct Emit6d {
  NormalModesProxy mode;
  FixedArray2D<Real, 6, 6> sigma_mat;
  RadIntAllEleProxy rad_int_by_ele;
};
Bmad::Emit6d emit_6d(
    EleProxy& ele_ref,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit = std::nullopt);
extern "C" bool fortran_entering_element(
    void* orbit /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    bool& is_entering /* 0D_NOT_logical inout */);
void entering_element(CoordProxy& orbit, int particle_at, bool& is_entering);
extern "C" void fortran_envelope_radints(
    std::complex<double>* Lambda /* 2D_NOT_complex inout */,
    std::complex<double>* Theta /* 2D_NOT_complex inout */,
    std::complex<double>* Iota /* 2D_NOT_complex inout */,
    double* alpha /* 1D_NOT_real inout */,
    double* emit /* 1D_NOT_real inout */);
void envelope_radints(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    FixedArray1D<Real, 3> alpha,
    FixedArray1D<Real, 3> emit);
extern "C" void fortran_envelope_radints_ibs(
    std::complex<double>* Lambda /* 2D_NOT_complex in */,
    std::complex<double>* Theta /* 2D_NOT_complex in */,
    std::complex<double>* Iota /* 2D_NOT_complex in */,
    void* eles /* 1D_ALLOC_type in */,
    double* alpha /* 1D_NOT_real out */,
    double* emit /* 1D_NOT_real out */,
    void* mode /* 0D_NOT_type in */,
    bool& tail_cut /* 0D_NOT_logical in */,
    double& npart /* 0D_NOT_real in */,
    int& species /* 0D_NOT_integer in */);
struct EnvelopeRadintsIbs {
  FixedArray1D<Real, 3> alpha;
  FixedArray1D<Real, 3> emit;
};
Bmad::EnvelopeRadintsIbs envelope_radints_ibs(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    EleProxyAlloc1D& eles,
    NormalModesProxy& mode,
    bool tail_cut,
    double npart,
    int species);
extern "C" bool fortran_eq_ac_kicker(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ac_kicker(AcKickerProxy& f1, AcKickerProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_ac_kicker_freq(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ac_kicker_freq(
    AcKickerFreqProxy& f1,
    AcKickerFreqProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_ac_kicker_time(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ac_kicker_time(
    AcKickerTimeProxy& f1,
    AcKickerTimeProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_anormal_mode(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_anormal_mode(AnormalModeProxy& f1, AnormalModeProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_aperture_param(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_aperture_param(
    ApertureParamProxy& f1,
    ApertureParamProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_aperture_point(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_aperture_point(
    AperturePointProxy& f1,
    AperturePointProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_aperture_scan(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_aperture_scan(
    ApertureScanProxy& f1,
    ApertureScanProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_beam(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_beam(BeamProxy& f1, BeamProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_beam_init(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_beam_init(BeamInitProxy& f1, BeamInitProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_bmad_common(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_bmad_common(BmadCommonProxy& f1, BmadCommonProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_bookkeeping_state(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_bookkeeping_state(
    BookkeepingStateProxy& f1,
    BookkeepingStateProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_bpm_phase_coupling(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_bpm_phase_coupling(
    BpmPhaseCouplingProxy& f1,
    BpmPhaseCouplingProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_branch(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_branch(BranchProxy& f1, BranchProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_bunch(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_bunch(BunchProxy& f1, BunchProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_bunch_params(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_bunch_params(BunchParamsProxy& f1, BunchParamsProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_cartesian_map(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cartesian_map(
    CartesianMapProxy& f1,
    CartesianMapProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_cartesian_map_term(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cartesian_map_term(
    CartesianMapTermProxy& f1,
    CartesianMapTermProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_cartesian_map_term1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cartesian_map_term1(
    CartesianMapTerm1Proxy& f1,
    CartesianMapTerm1Proxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_complex_taylor(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_complex_taylor(
    ComplexTaylorProxy& f1,
    ComplexTaylorProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_complex_taylor_term(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_complex_taylor_term(
    ComplexTaylorTermProxy& f1,
    ComplexTaylorTermProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_control(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_control(ControlProxy& f1, ControlProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_control_ramp1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_control_ramp1(
    ControlRamp1Proxy& f1,
    ControlRamp1Proxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_control_var1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_control_var1(ControlVar1Proxy& f1, ControlVar1Proxy& f2, bool& is_eq);
extern "C" bool fortran_eq_controller(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_controller(ControllerProxy& f1, ControllerProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_coord(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_coord(CoordProxy& f1, CoordProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_coord_array(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_coord_array(CoordArrayProxy& f1, CoordArrayProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_cylindrical_map(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cylindrical_map(
    CylindricalMapProxy& f1,
    CylindricalMapProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_cylindrical_map_term(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cylindrical_map_term(
    CylindricalMapTermProxy& f1,
    CylindricalMapTermProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_cylindrical_map_term1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_cylindrical_map_term1(
    CylindricalMapTerm1Proxy& f1,
    CylindricalMapTerm1Proxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_ele(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ele(EleProxy& f1, EleProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_ellipse_beam_init(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ellipse_beam_init(
    EllipseBeamInitProxy& f1,
    EllipseBeamInitProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_em_field(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_em_field(EmFieldProxy& f1, EmFieldProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_em_taylor(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_em_taylor(EmTaylorProxy& f1, EmTaylorProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_em_taylor_term(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_em_taylor_term(
    EmTaylorTermProxy& f1,
    EmTaylorTermProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_expression_atom(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_expression_atom(
    ExpressionAtomProxy& f1,
    ExpressionAtomProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_floor_position(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_floor_position(
    FloorPositionProxy& f1,
    FloorPositionProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_gen_grad1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_gen_grad1(GenGrad1Proxy& f1, GenGrad1Proxy& f2, bool& is_eq);
extern "C" bool fortran_eq_gen_grad_map(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_gen_grad_map(GenGradMapProxy& f1, GenGradMapProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_grid_beam_init(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_grid_beam_init(
    GridBeamInitProxy& f1,
    GridBeamInitProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_grid_field(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_grid_field(GridFieldProxy& f1, GridFieldProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_grid_field_pt(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_grid_field_pt(GridFieldPtProxy& f1, GridFieldPtProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_grid_field_pt1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_grid_field_pt1(
    GridFieldPt1Proxy& f1,
    GridFieldPt1Proxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_high_energy_space_charge(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_high_energy_space_charge(
    HighEnergySpaceChargeProxy& f1,
    HighEnergySpaceChargeProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_interval1_coef(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_interval1_coef(
    Interval1CoefProxy& f1,
    Interval1CoefProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_kv_beam_init(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_kv_beam_init(KvBeamInitProxy& f1, KvBeamInitProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_lat(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_lat(LatProxy& f1, LatProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_lat_ele_loc(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_lat_ele_loc(LatEleLocProxy& f1, LatEleLocProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_lat_param(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_lat_param(LatParamProxy& f1, LatParamProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_linac_normal_mode(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_linac_normal_mode(
    LinacNormalModeProxy& f1,
    LinacNormalModeProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_mode3(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_mode3(Mode3Proxy& f1, Mode3Proxy& f2, bool& is_eq);
extern "C" bool fortran_eq_mode_info(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_mode_info(ModeInfoProxy& f1, ModeInfoProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_normal_modes(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_normal_modes(NormalModesProxy& f1, NormalModesProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_photon_element(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_photon_element(
    PhotonElementProxy& f1,
    PhotonElementProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_photon_material(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_photon_material(
    PhotonMaterialProxy& f1,
    PhotonMaterialProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_photon_reflect_surface(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_photon_reflect_surface(
    PhotonReflectSurfaceProxy& f1,
    PhotonReflectSurfaceProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_photon_reflect_table(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_photon_reflect_table(
    PhotonReflectTableProxy& f1,
    PhotonReflectTableProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_photon_target(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_photon_target(
    PhotonTargetProxy& f1,
    PhotonTargetProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_pixel_detec(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_pixel_detec(PixelDetecProxy& f1, PixelDetecProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_pixel_pt(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_pixel_pt(PixelPtProxy& f1, PixelPtProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_pre_tracker(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_pre_tracker(PreTrackerProxy& f1, PreTrackerProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_rad_int1(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_rad_int1(RadInt1Proxy& f1, RadInt1Proxy& f2, bool& is_eq);
extern "C" bool fortran_eq_rad_int_all_ele(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_rad_int_all_ele(
    RadIntAllEleProxy& f1,
    RadIntAllEleProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_rad_int_branch(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_rad_int_branch(
    RadIntBranchProxy& f1,
    RadIntBranchProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_rad_map(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_rad_map(RadMapProxy& f1, RadMapProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_rad_map_ele(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_rad_map_ele(RadMapEleProxy& f1, RadMapEleProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_ramper_lord(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_ramper_lord(RamperLordProxy& f1, RamperLordProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_space_charge_common(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_space_charge_common(
    SpaceChargeCommonProxy& f1,
    SpaceChargeCommonProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_spin_polar(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_spin_polar(SpinPolarProxy& f1, SpinPolarProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_spline(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_spline(SplineProxy& f1, SplineProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_strong_beam(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_strong_beam(StrongBeamProxy& f1, StrongBeamProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_surface_curvature(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_curvature(
    SurfaceCurvatureProxy& f1,
    SurfaceCurvatureProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_displacement(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_displacement(
    SurfaceDisplacementProxy& f1,
    SurfaceDisplacementProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_displacement_pt(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_displacement_pt(
    SurfaceDisplacementPtProxy& f1,
    SurfaceDisplacementPtProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_h_misalign(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_h_misalign(
    SurfaceHMisalignProxy& f1,
    SurfaceHMisalignProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_h_misalign_pt(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_h_misalign_pt(
    SurfaceHMisalignPtProxy& f1,
    SurfaceHMisalignPtProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_segmented(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_segmented(
    SurfaceSegmentedProxy& f1,
    SurfaceSegmentedProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_surface_segmented_pt(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_surface_segmented_pt(
    SurfaceSegmentedPtProxy& f1,
    SurfaceSegmentedPtProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_target_point(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_target_point(TargetPointProxy& f1, TargetPointProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_taylor(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_taylor(TaylorProxy& f1, TaylorProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_taylor_term(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_taylor_term(TaylorTermProxy& f1, TaylorTermProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_track(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_track(TrackProxy& f1, TrackProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_track_point(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_track_point(TrackPointProxy& f1, TrackPointProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_twiss(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_twiss(TwissProxy& f1, TwissProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake(WakeProxy& f1, WakeProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake_lr(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake_lr(WakeLrProxy& f1, WakeLrProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake_lr_mode(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake_lr_mode(WakeLrModeProxy& f1, WakeLrModeProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake_sr(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake_sr(WakeSrProxy& f1, WakeSrProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake_sr_mode(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake_sr_mode(WakeSrModeProxy& f1, WakeSrModeProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wake_sr_z_long(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wake_sr_z_long(WakeSrZLongProxy& f1, WakeSrZLongProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wall3d(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wall3d(Wall3dProxy& f1, Wall3dProxy& f2, bool& is_eq);
extern "C" bool fortran_eq_wall3d_section(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wall3d_section(
    Wall3dSectionProxy& f1,
    Wall3dSectionProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_wall3d_vertex(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_wall3d_vertex(
    Wall3dVertexProxy& f1,
    Wall3dVertexProxy& f2,
    bool& is_eq);
extern "C" bool fortran_eq_xy_disp(
    void* f1 /* 0D_NOT_type in */,
    void* f2 /* 0D_NOT_type in */,
    bool& is_eq /* 0D_NOT_logical inout */);
void eq_xy_disp(XyDispProxy& f1, XyDispProxy& f2, bool& is_eq);
extern "C" bool fortran_equal_sign_here(
    void* ele /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& is_here /* 0D_NOT_logical inout */);
void equal_sign_here(EleProxy& ele, std::string& delim, bool& is_here);
extern "C" bool fortran_equivalent_taylor_attributes(
    void* ele_taylor /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    bool& equiv /* 0D_NOT_logical inout */);
void equivalent_taylor_attributes(
    EleProxy& ele_taylor,
    EleProxy& ele2,
    bool& equiv);
extern "C" void fortran_etdiv(
    double& A /* 0D_NOT_real inout */,
    double& B /* 0D_NOT_real inout */,
    double& C /* 0D_NOT_real inout */,
    double& D /* 0D_NOT_real inout */,
    double& E /* 0D_NOT_real inout */,
    double& F /* 0D_NOT_real inout */);
void etdiv(double& A, double& B, double& C, double& D, double& E, double& F);

// Skipped unusable routine ety:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine ety2:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine etyt:
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_evaluate_array_index(
    bool& err_flag /* 0D_NOT_logical out */,
    const char* delim_list1 /* 0D_NOT_character in */,
    const char* word2 /* 0D_NOT_character out */,
    const char* delim_list2 /* 0D_NOT_character in */,
    const char* delim2 /* 0D_NOT_character out */,
    int& this_index /* 0D_NOT_integer out */);
struct EvaluateArrayIndex {
  bool err_flag;
  std::string word2;
  std::string delim2;
  int this_index;
};
Bmad::EvaluateArrayIndex evaluate_array_index(
    std::string delim_list1,
    std::string delim_list2);
extern "C" bool fortran_evaluate_logical(
    const char* word /* 0D_NOT_character in */,
    int& iostat /* 0D_NOT_integer out */,
    bool& this_logic /* 0D_NOT_logical out */);
struct EvaluateLogical {
  int iostat;
  bool this_logic;
};
Bmad::EvaluateLogical evaluate_logical(std::string word);
extern "C" void fortran_exact_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orb /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void exact_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_exp_bessi0(
    double& t /* 0D_NOT_real in */,
    double& B1 /* 0D_NOT_real in */,
    double& B2 /* 0D_NOT_real in */,
    double& func_retval__ /* 0D_NOT_real out */);
double exp_bessi0(double t, double B1, double B2);
extern "C" bool fortran_expect_one_of(
    const char* delim_list /* 0D_NOT_character in */,
    bool& check_input_delim /* 0D_NOT_logical in */,
    const char* ele_name /* 0D_NOT_character in */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& is_ok /* 0D_NOT_logical inout */);
void expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string& delim,
    bool& delim_found,
    bool& is_ok);
extern "C" bool fortran_expect_this(
    const char* expecting /* 0D_NOT_character in */,
    bool& check_delim /* 0D_NOT_logical in */,
    bool& call_check /* 0D_NOT_logical in */,
    const char* err_str /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    const char* delim /* 0D_NOT_character out */,
    bool& delim_found /* 0D_NOT_logical out */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ExpectThis {
  std::string delim;
  bool delim_found;
  bool is_ok;
};
Bmad::ExpectThis expect_this(
    std::string expecting,
    bool check_delim,
    bool call_check,
    std::string err_str,
    EleProxy& ele);
extern "C" bool fortran_expression_stack_to_string(
    void* stack /* 1D_ALLOC_type in */,
    bool* polish /* 0D_NOT_logical in */,
    const char* str /* 0D_ALLOC_character out */);
std::string expression_stack_to_string(
    ExpressionAtomProxyAlloc1D& stack,
    std::optional<bool> polish = std::nullopt);
extern "C" bool fortran_expression_stack_value(
    void* stack /* 1D_ALLOC_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    const char* err_str /* 0D_NOT_character out */,
    void* var /* 1D_ALLOC_type in */,
    bool* use_old /* 0D_NOT_logical in */,
    double& value /* 0D_NOT_real out */);
struct ExpressionStackValue {
  bool err_flag;
  std::string err_str;
  double value;
};
Bmad::ExpressionStackValue expression_stack_value(
    ExpressionAtomProxyAlloc1D& stack,
    optional_ref<ControlVar1ProxyAlloc1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt);
extern "C" void fortran_expression_string_to_stack(
    const char* string /* 0D_NOT_character in */,
    void* stack /* 1D_ALLOC_type out */,
    int& n_stack /* 0D_NOT_integer out */,
    bool& err_flag /* 0D_NOT_logical out */,
    const char* err_str /* 0D_NOT_character out */);
struct ExpressionStringToStack {
  ExpressionAtomProxyAlloc1D stack;
  int n_stack;
  bool err_flag;
  std::string err_str;
};
Bmad::ExpressionStringToStack expression_string_to_stack(std::string string);
extern "C" void fortran_expression_string_to_tree(
    const char* string /* 0D_NOT_character in */,
    void* root_tree /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    const char* err_str /* 0D_NOT_character out */);
struct ExpressionStringToTree {
  bool err_flag;
  std::string err_str;
};
Bmad::ExpressionStringToTree expression_string_to_tree(
    std::string string,
    ExpressionTreeProxy& root_tree);
extern "C" bool fortran_expression_tree_to_string(
    void* tree /* 0D_NOT_type in */,
    bool* include_root /* 0D_NOT_logical in */,
    int* n_node /* 0D_NOT_integer in */,
    void* parent /* 0D_NOT_type in */,
    const char* str_out /* 0D_ALLOC_character out */);
std::string expression_tree_to_string(
    ExpressionTreeProxy& tree,
    std::optional<bool> include_root = std::nullopt,
    std::optional<int> n_node = std::nullopt,
    optional_ref<ExpressionTreeProxy> parent = std::nullopt);
extern "C" bool fortran_expression_value(
    const char* expression /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */,
    const char* err_str /* 0D_NOT_character out */,
    void* var /* 1D_ALLOC_type in */,
    bool* use_old /* 0D_NOT_logical in */,
    double& value /* 0D_NOT_real out */);
struct ExpressionValue {
  bool err_flag;
  std::string err_str;
  double value;
};
Bmad::ExpressionValue expression_value(
    std::string expression,
    optional_ref<ControlVar1ProxyAlloc1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt);
extern "C" void fortran_fft1(
    void* a /* 1D_ALLOC_real inout */,
    void* b /* 1D_ALLOC_real inout */,
    int& n /* 0D_NOT_integer in */,
    int& isn /* 0D_NOT_integer in */,
    int& ierr /* 0D_NOT_integer out */);
int fft1(RealAlloc1D& a, RealAlloc1D& b, int n, int isn);

// Skipped unusable routine fftconvcorr3d:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine fibre_to_ele:
// - Untranslated type: fibre (0D)
extern "C" bool fortran_field_attribute_free(
    void* ele /* 0D_NOT_type in */,
    const char* attrib_name /* 0D_NOT_character in */,
    bool& free /* 0D_NOT_logical out */);
bool field_attribute_free(EleProxy& ele, std::string attrib_name);

// Skipped unusable routine field_interpolate_3d:
// - Variable in sized array: field_mesh(0:,0:,0:) 3D_NOT_real
extern "C" void fortran_finalize_reflectivity_table(
    void* table /* 0D_NOT_type inout */,
    bool& in_degrees /* 0D_NOT_logical in */);
void finalize_reflectivity_table(
    PhotonReflectTableProxy& table,
    bool in_degrees);
extern "C" void fortran_find_element_ends(
    void* ele /* 0D_NOT_type in */,
    void* ele1 /* 0D_PTR_type out */,
    void* ele2 /* 0D_PTR_type out */,
    int* ix_multipass /* 0D_NOT_integer in */);
struct FindElementEnds {
  EleProxy ele1;
  EleProxy ele2;
};
Bmad::FindElementEnds find_element_ends(
    EleProxy& ele,
    std::optional<int> ix_multipass = std::nullopt);
extern "C" void fortran_find_fwhm(
    double& bound /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    double& fwhm /* 0D_NOT_real out */);
double find_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_find_matching_fieldmap(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    int& fm_type /* 0D_NOT_integer in */,
    void* match_ele /* 0D_PTR_type out */,
    int& ix_field /* 0D_NOT_integer out */,
    bool* ignore_slaves /* 0D_NOT_logical in */);
struct FindMatchingFieldmap {
  EleProxy match_ele;
  int ix_field;
};
Bmad::FindMatchingFieldmap find_matching_fieldmap(
    std::string file_name,
    EleProxy& ele,
    int fm_type,
    std::optional<bool> ignore_slaves = std::nullopt);
extern "C" void fortran_find_normalization(
    double& bound /* 0D_NOT_real in */,
    double& p0 /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    double& pnrml /* 0D_NOT_real out */);
double find_normalization(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_floor_angles_to_w_mat(
    double& theta /* 0D_NOT_real in */,
    double& phi /* 0D_NOT_real in */,
    double& psi /* 0D_NOT_real in */,
    double* w_mat /* 2D_NOT_real out */,
    double* w_mat_inv /* 2D_NOT_real out */);
struct FloorAnglesToWMat {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat_inv;
};
Bmad::FloorAnglesToWMat floor_angles_to_w_mat(
    double theta,
    double phi,
    double psi);
extern "C" void fortran_floor_w_mat_to_angles(
    double* w_mat /* 2D_NOT_real in */,
    double& theta /* 0D_NOT_real out */,
    double& phi /* 0D_NOT_real out */,
    double& psi /* 0D_NOT_real out */,
    void* floor0 /* 0D_NOT_type in */);
struct FloorWMatToAngles {
  double theta;
  double phi;
  double psi;
};
Bmad::FloorWMatToAngles floor_w_mat_to_angles(
    FixedArray2D<Real, 3, 3> w_mat,
    optional_ref<FloorPositionProxy> floor0 = std::nullopt);
extern "C" void fortran_form_complex_taylor(
    void* re_taylor /* 0D_NOT_type in */,
    void* im_taylor /* 0D_NOT_type in */,
    void* complex_taylor /* 0D_NOT_type out */);
ComplexTaylorProxy form_complex_taylor(
    TaylorProxy& re_taylor,
    TaylorProxy& im_taylor);
extern "C" void fortran_form_digested_bmad_file_name(
    const char* lat_file /* 0D_NOT_character in */,
    const char* digested_file /* 0D_NOT_character out */,
    const char* full_lat_file /* 0D_NOT_character out */,
    const char* use_line /* 0D_NOT_character in */);
struct FormDigestedBmadFileName {
  std::string digested_file;
  std::string full_lat_file;
};
Bmad::FormDigestedBmadFileName form_digested_bmad_file_name(
    std::string lat_file,
    std::optional<std::string> use_line = std::nullopt);
extern "C" bool fortran_fringe_here(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    bool& is_here /* 0D_NOT_logical inout */);
void fringe_here(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    bool& is_here);
extern "C" bool fortran_g_bend_from_em_field(
    double* b /* 1D_NOT_real in */,
    double* e /* 1D_NOT_real in */,
    void* orbit /* 0D_NOT_type in */,
    double* g_bend /* 1D_NOT_real out */);
FixedArray1D<Real, 3> g_bend_from_em_field(
    FixedArray1D<Real, 3> b,
    FixedArray1D<Real, 3> e,
    CoordProxy& orbit);
extern "C" void fortran_g_bending_strength_from_em_field(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& s_rel /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type in */,
    bool& local_ref_frame /* 0D_NOT_logical in */,
    double* g /* 1D_NOT_real out */,
    double* dg /* 2D_NOT_real out */);
struct GBendingStrengthFromEmField {
  FixedArray1D<Real, 3> g;
  std::optional<FixedArray2D<Real, 3, 3>> dg;
};
Bmad::GBendingStrengthFromEmField g_bending_strength_from_em_field(
    EleProxy& ele,
    LatParamProxy& param,
    double s_rel,
    CoordProxy& orbit,
    bool local_ref_frame);
extern "C" void fortran_g_integrals_calc(void* lat /* 0D_NOT_type inout */);
void g_integrals_calc(LatProxy& lat);
extern "C" bool fortran_gamma_ref(
    void* ele /* 0D_NOT_type in */,
    double& gamma /* 0D_NOT_real inout */);
void gamma_ref(EleProxy& ele, double& gamma);
extern "C" void fortran_gen_grad1_to_em_taylor(
    void* ele /* 0D_NOT_type in */,
    void* gen_grad /* 0D_NOT_type in */,
    int& iz /* 0D_NOT_integer in */,
    void* em_taylor /* 1D_NOT_type out */);
EmTaylorProxyArray1D gen_grad1_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    int iz);
extern "C" void fortran_gen_grad_at_s_to_em_taylor(
    void* ele /* 0D_NOT_type in */,
    void* gen_grad /* 0D_NOT_type in */,
    double& s_pos /* 0D_NOT_real in */,
    void* em_taylor /* 1D_NOT_type out */);
EmTaylorProxyArray1D gen_grad_at_s_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    double s_pos);
extern "C" bool fortran_gen_grad_field(
    void* deriv /* 1D_ALLOC_real inout */,
    void* gg /* 0D_NOT_type inout */,
    double& rho /* 0D_NOT_real inout */,
    double& theta /* 0D_NOT_real inout */,
    double* field /* 1D_NOT_real inout */);
void gen_grad_field(
    RealAlloc1D& deriv,
    GenGrad1Proxy& gg,
    double& rho,
    double& theta,
    FixedArray1D<Real, 3> field);

// Skipped unusable routine get_astra_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_get_bl_from_fwhm(
    double& bound /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    double& sigma /* 0D_NOT_real out */);
double get_bl_from_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_get_called_file(
    const char* delim /* 0D_NOT_character inout */,
    const char* call_file /* 0D_NOT_character inout */,
    bool& err /* 0D_NOT_logical inout */);
void get_called_file(std::string& delim, std::string& call_file, bool& err);

// Skipped unusable routine get_cgrn_csr3d:
// - Variable inout sized array: cgrn(0:,0:,0:) 3D_NOT_complex
extern "C" void fortran_get_emit_from_sigma_mat(
    double* sigma_mat /* 2D_NOT_real in */,
    double* normal /* 1D_NOT_real out */,
    double* Nmat /* 2D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct GetEmitFromSigmaMat {
  FixedArray1D<Real, 3> normal;
  bool err_flag;
};
Bmad::GetEmitFromSigmaMat get_emit_from_sigma_mat(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<FixedArray2D<Real, 6, 6>> Nmat = std::nullopt);

// Skipped unusable routine get_gpt_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)

// Skipped unusable routine get_list_of_names:
// - Variable-sized inout character array: name_list(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_get_next_word(
    const char* word /* 0D_NOT_character in */,
    int& ix_word /* 0D_NOT_integer in */,
    const char* delim_list /* 0D_NOT_character in */,
    const char* delim /* 0D_NOT_character in */,
    bool& delim_found /* 0D_NOT_logical in */,
    bool* upper_case_word /* 0D_NOT_logical in */,
    bool* call_check /* 0D_NOT_logical in */,
    bool* err_flag /* 0D_NOT_logical in */);
void get_next_word(
    std::string word,
    int ix_word,
    std::string delim_list,
    std::string delim,
    bool delim_found,
    std::optional<bool> upper_case_word = std::nullopt,
    std::optional<bool> call_check = std::nullopt,
    std::optional<bool> err_flag = std::nullopt);

// Skipped unusable routine get_opal_fieldgrid_name_and_scaling:
// - Untranslated type: str_index_struct (0D)

// Skipped unusable routine get_overlay_group_names:
// - Untranslated type: parser_ele_struct (0D)
// - Variable-sized inout character array: names_out(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine get_sequence_args:
// - Variable-sized inout character array: arg_list(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_get_slave_list(
    void* lord /* 0D_NOT_type in */,
    void* slaves /* 1D_ALLOC_type out */,
    int& n_slave /* 0D_NOT_integer out */);
struct GetSlaveList {
  ElePointerProxyAlloc1D slaves;
  int n_slave;
};
Bmad::GetSlaveList get_slave_list(EleProxy& lord);

// Skipped unusable routine get_switch:
// - Variable-sized inout character array: name_list(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine getrhotilde:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_gpt_field_grid_scaling(
    void* ele /* 0D_NOT_type inout */,
    int& dimensions /* 0D_NOT_integer inout */,
    double& field_scale /* 0D_NOT_real inout */,
    double& ref_time /* 0D_NOT_real inout */);
void gpt_field_grid_scaling(
    EleProxy& ele,
    int& dimensions,
    double& field_scale,
    double& ref_time);
extern "C" bool fortran_gpt_max_field_reference(
    void* pt0 /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */,
    double& field_value /* 0D_NOT_real inout */);
void gpt_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double& field_value);
extern "C" void fortran_gpt_to_particle_bunch(
    const char* gpt_file /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* bunch /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct GptToParticleBunch {
  BunchProxy bunch;
  bool err_flag;
};
Bmad::GptToParticleBunch gpt_to_particle_bunch(
    std::string gpt_file,
    EleProxy& ele);
extern "C" bool fortran_gradient_shift_sr_wake(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& grad_shift /* 0D_NOT_real inout */);
void gradient_shift_sr_wake(
    EleProxy& ele,
    LatParamProxy& param,
    double& grad_shift);
extern "C" void fortran_grid_field_interpolate(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    void* grid /* 0D_NOT_type in */,
    void* g_field /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical in */,
    double& x1 /* 0D_NOT_real in */,
    double* x2 /* 0D_NOT_real in */,
    double* x3 /* 0D_NOT_real in */,
    bool* allow_s_out_of_bounds /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */);
GridFieldPt1Proxy grid_field_interpolate(
    EleProxy& ele,
    CoordProxy& orbit,
    GridFieldProxy& grid,
    bool err_flag,
    double x1,
    std::optional<double> x2 = std::nullopt,
    std::optional<double> x3 = std::nullopt,
    std::optional<bool> allow_s_out_of_bounds = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_hard_multipole_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void hard_multipole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_has_attribute(
    void* ele /* 0D_NOT_type inout */,
    const char* attrib /* 0D_NOT_character inout */,
    bool& has_it /* 0D_NOT_logical inout */);
void has_attribute(EleProxy& ele, std::string& attrib, bool& has_it);
extern "C" bool fortran_has_curvature(
    void* phot_ele /* 0D_NOT_type in */,
    bool& curved /* 0D_NOT_logical out */);
bool has_curvature(PhotonElementProxy& phot_ele);
extern "C" bool fortran_has_orientation_attributes(
    void* ele /* 0D_NOT_type in */,
    bool& has_attribs /* 0D_NOT_logical out */);
bool has_orientation_attributes(EleProxy& ele);

// Skipped unusable routine hdf5_read_beam:
// - Untranslated type: pmd_header_struct (0D)

// Skipped unusable routine hdf5_read_grid_field:
// - Untranslated type: pmd_header_struct (0D)
extern "C" void fortran_hdf5_write_beam(
    const char* file_name /* 0D_NOT_character inout */,
    void* bunches /* 1D_ALLOC_type inout */,
    bool& append /* 0D_NOT_logical inout */,
    bool& error /* 0D_NOT_logical inout */,
    void* lat /* 0D_NOT_type inout */,
    bool* alive_only /* 0D_NOT_logical inout */);
void hdf5_write_beam(
    std::string& file_name,
    BunchProxyAlloc1D& bunches,
    bool& append,
    bool& error,
    optional_ref<LatProxy> lat = std::nullopt,
    optional_ref<bool> alive_only = std::nullopt);
extern "C" void fortran_hdf5_write_grid_field(
    const char* file_name /* 0D_NOT_character inout */,
    void* ele /* 0D_NOT_type inout */,
    void* g_field /* 1D_ALLOC_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void hdf5_write_grid_field(
    std::string& file_name,
    EleProxy& ele,
    GridFieldProxyAlloc1D& g_field,
    bool& err_flag);
extern "C" void fortran_hwang_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orb /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void hwang_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);

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
    double* sigma_mat /* 2D_NOT_real inout */,
    bool& tail_cut /* 0D_NOT_logical inout */,
    double& tau /* 0D_NOT_real inout */,
    double& energy /* 0D_NOT_real inout */,
    double& n_part /* 0D_NOT_real inout */,
    int& species /* 0D_NOT_integer inout */,
    double* ibs_mat /* 2D_NOT_real inout */);
void ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool& tail_cut,
    double& tau,
    double& energy,
    double& n_part,
    int& species,
    FixedArray2D<Real, 6, 6> ibs_mat);

// Skipped unusable routine ibs_rates1turn:
// - Untranslated type: ibs_sim_param_struct (0D)
// - Untranslated type: ibs_struct (0D)
extern "C" bool fortran_igfcoulombfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& dx /* 0D_NOT_real inout */,
    double& dy /* 0D_NOT_real inout */,
    double& dz /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void igfcoulombfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res);
extern "C" bool fortran_igfexfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& dx /* 0D_NOT_real inout */,
    double& dy /* 0D_NOT_real inout */,
    double& dz /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void igfexfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res);
extern "C" bool fortran_igfeyfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& dx /* 0D_NOT_real inout */,
    double& dy /* 0D_NOT_real inout */,
    double& dz /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void igfeyfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res);
extern "C" bool fortran_igfezfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& dx /* 0D_NOT_real inout */,
    double& dy /* 0D_NOT_real inout */,
    double& dz /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void igfezfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& dx,
    double& dy,
    double& dz,
    double& res);

// Skipped unusable routine image_charge_kick_calc:
// - Untranslated type: csr_kick1_struct (0D)
// - Untranslated type: csr_struct (0D)

// Skipped unusable routine imageconvcorr3d:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_init_attribute_name1(
    int& ix_key /* 0D_NOT_integer in */,
    int& ix_attrib /* 0D_NOT_integer in */,
    const char* name /* 0D_NOT_character in */,
    int* attrib_state /* 0D_NOT_integer in */,
    bool* override /* 0D_NOT_logical in */);
void init_attribute_name1(
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state = std::nullopt,
    std::optional<bool> override = std::nullopt);
extern "C" void fortran_init_attribute_name_array();
void init_attribute_name_array();
extern "C" void fortran_init_beam_distribution(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* beam_init /* 0D_NOT_type in */,
    void* beam /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* modes /* 0D_NOT_type in */,
    void* beam_init_set /* 0D_NOT_type out */,
    bool* print_p0c_shift_warning /* 0D_NOT_logical in */,
    bool* conserve_momentum /* 0D_NOT_logical inout */);
struct InitBeamDistribution {
  BeamProxy beam;
  bool err_flag;
  BeamInitProxy beam_init_set;
};
Bmad::InitBeamDistribution init_beam_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_init_bmad();
void init_bmad();
extern "C" void fortran_init_bmad_parser_common(
    void* lat /* 0D_NOT_type inout */);
void init_bmad_parser_common(optional_ref<LatProxy> lat = std::nullopt);
extern "C" void fortran_init_bunch_distribution(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* beam_init /* 0D_NOT_type in */,
    int& ix_bunch /* 0D_NOT_integer in */,
    void* bunch /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* modes /* 0D_NOT_type in */,
    void* beam_init_used /* 0D_NOT_type out */,
    bool* print_p0c_shift_warning /* 0D_NOT_logical in */,
    bool* conserve_momentum /* 0D_NOT_logical inout */);
struct InitBunchDistribution {
  BunchProxy bunch;
  bool err_flag;
  BeamInitProxy beam_init_used;
};
Bmad::InitBunchDistribution init_bunch_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    int ix_bunch,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_init_complex_taylor_series(
    void* complex_taylor /* 0D_NOT_type inout */,
    int& n_term /* 0D_NOT_integer in */,
    bool* save /* 0D_NOT_logical in */);
void init_complex_taylor_series(
    ComplexTaylorProxy& complex_taylor,
    int n_term,
    std::optional<bool> save = std::nullopt);
extern "C" void fortran_init_coord1(
    void* orb /* 0D_NOT_type inout */,
    double* vec /* 1D_NOT_real in */,
    void* ele /* 0D_NOT_type in */,
    int* element_end /* 0D_NOT_integer in */,
    int* particle /* 0D_NOT_integer in */,
    int* direction /* 0D_NOT_integer in */,
    double* E_photon /* 0D_NOT_real in */,
    double* t_offset /* 0D_NOT_real in */,
    bool* shift_vec6 /* 0D_NOT_logical in */,
    double* spin /* 1D_NOT_real in */,
    double* s_pos /* 0D_NOT_real in */,
    bool* random_on /* 0D_NOT_logical in */);
void init_coord(
    CoordProxy& orb,
    FixedArray1D<Real, 6> vec,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> random_on = std::nullopt);
extern "C" void fortran_init_coord2(
    void* orb_out /* 0D_NOT_type out */,
    void* orb_in /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    int* element_end /* 0D_NOT_integer in */,
    int* particle /* 0D_NOT_integer in */,
    int* direction /* 0D_NOT_integer in */,
    double* E_photon /* 0D_NOT_real in */,
    double* t_offset /* 0D_NOT_real in */,
    bool* shift_vec6 /* 0D_NOT_logical in */,
    double* spin /* 1D_NOT_real in */,
    double* s_pos /* 0D_NOT_real in */,
    bool* random_on /* 0D_NOT_logical in */);
CoordProxy init_coord(
    CoordProxy& orb_in,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> random_on = std::nullopt);
extern "C" void fortran_init_coord3(
    void* orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    int* element_end /* 0D_NOT_integer in */,
    int* particle /* 0D_NOT_integer in */,
    int* direction /* 0D_NOT_integer in */,
    double* E_photon /* 0D_NOT_real in */,
    double* t_offset /* 0D_NOT_real in */,
    bool* shift_vec6 /* 0D_NOT_logical in */,
    double* spin /* 1D_NOT_real in */);
void init_coord(
    CoordProxy& orb,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt);
extern "C" void fortran_init_custom(void* lat /* 0D_NOT_type inout */);
void init_custom(LatProxy& lat);

// Skipped unusable routine init_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_init_ele(
    void* ele /* 0D_NOT_type out */,
    int* key /* 0D_NOT_integer in */,
    int* sub_key /* 0D_NOT_integer in */,
    int* ix_ele /* 0D_NOT_integer in */,
    void* branch /* 0D_NOT_type in */);
EleProxy init_ele(
    std::optional<int> key = std::nullopt,
    std::optional<int> sub_key = std::nullopt,
    std::optional<int> ix_ele = std::nullopt,
    optional_ref<BranchProxy> branch = std::nullopt);
extern "C" void fortran_init_em_taylor_series(
    void* em_taylor /* 0D_NOT_type inout */,
    int& n_term /* 0D_NOT_integer in */,
    bool* save_old /* 0D_NOT_logical in */);
void init_em_taylor_series(
    EmTaylorProxy& em_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt);

// Skipped unusable routine init_fringe_info:
// - Untranslated type: fringe_field_info_struct (0D)
extern "C" void fortran_init_lat(
    void* lat /* 0D_NOT_type out */,
    int* n /* 0D_NOT_integer in */,
    bool* init_beginning_ele /* 0D_NOT_logical in */);
LatProxy init_lat(
    std::optional<int> n = std::nullopt,
    std::optional<bool> init_beginning_ele = std::nullopt);
extern "C" void fortran_init_multipole_cache(void* ele /* 0D_NOT_type inout */);
void init_multipole_cache(EleProxy& ele);
extern "C" void fortran_init_photon_from_a_photon_init_ele(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type out */,
    bool* random_on /* 0D_NOT_logical in */);
CoordProxy init_photon_from_a_photon_init_ele(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> random_on = std::nullopt);
extern "C" bool fortran_init_photon_integ_prob(
    double& gamma /* 0D_NOT_real in */,
    double& g /* 0D_NOT_real in */,
    double& E_min /* 0D_NOT_real in */,
    double& E_max /* 0D_NOT_real in */,
    double* vert_angle_min /* 0D_NOT_real in */,
    double* vert_angle_max /* 0D_NOT_real in */,
    bool* vert_angle_symmetric /* 0D_NOT_logical in */,
    double* energy_integ_prob /* 0D_NOT_real in */,
    double& E_photon /* 0D_NOT_real out */,
    double& integ_prob /* 0D_NOT_real out */);
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
    std::optional<double> energy_integ_prob = std::nullopt);
extern "C" void fortran_init_spin_distribution(
    void* beam_init /* 0D_NOT_type in */,
    void* bunch /* 0D_NOT_type out */,
    void* ele /* 0D_NOT_type inout */);
BunchProxy init_spin_distribution(BeamInitProxy& beam_init, EleProxy& ele);
extern "C" void fortran_init_surface_segment(
    void* phot /* 0D_NOT_type in */,
    int& ix /* 0D_NOT_integer inout */,
    int& iy /* 0D_NOT_integer inout */);
void init_surface_segment(PhotonElementProxy& phot, int& ix, int& iy);
extern "C" void fortran_init_taylor_series(
    void* bmad_taylor /* 0D_NOT_type inout */,
    int& n_term /* 0D_NOT_integer in */,
    bool* save_old /* 0D_NOT_logical in */);
void init_taylor_series(
    TaylorProxy& bmad_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt);
extern "C" void fortran_init_wake(
    void* wake /* 0D_PTR_type out */,
    int& n_sr_long /* 0D_NOT_integer in */,
    int& n_sr_trans /* 0D_NOT_integer in */,
    int& n_sr_z /* 0D_NOT_integer in */,
    int& n_lr_mode /* 0D_NOT_integer in */,
    bool* always_allocate /* 0D_NOT_logical in */);
WakeProxy init_wake(
    int n_sr_long,
    int n_sr_trans,
    int n_sr_z,
    int n_lr_mode,
    std::optional<bool> always_allocate = std::nullopt);
extern "C" void fortran_insert_element(
    void* lat /* 0D_NOT_type inout */,
    void* insert_ele /* 0D_NOT_type in */,
    int& ix_ele /* 0D_NOT_integer in */,
    int* ix_branch /* 0D_NOT_integer in */,
    void* orbit /* 1D_ALLOC_type inout */);
void insert_element(
    LatProxy& lat,
    EleProxy& insert_ele,
    int ix_ele,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<CoordProxyAlloc1D> orbit = std::nullopt);

// Skipped unusable routine integrand:
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_integrand_base(
    double& t /* 0D_NOT_real in */,
    void* args /* 1D_ALLOC_real inout */,
    double& func_retval__ /* 0D_NOT_real inout */);
void integrand_base(double t, RealAlloc1D& args, double& func_retval__);

// Skipped unusable routine integrand_base_cov:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine integrand_zap:
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_integrate_psi(
    double& bound /* 0D_NOT_real in */,
    double& p0 /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    double& result /* 0D_NOT_real out */);
double integrate_psi(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_integrated_mats(
    void* eles /* 1D_ALLOC_type inout */,
    void* coos /* 1D_ALLOC_type inout */,
    std::complex<double>* Lambda /* 2D_NOT_complex inout */,
    std::complex<double>* Theta /* 2D_NOT_complex inout */,
    std::complex<double>* Iota /* 2D_NOT_complex inout */,
    void* mode /* 0D_NOT_type inout */);
void integrated_mats(
    EleProxyAlloc1D& eles,
    CoordProxyAlloc1D& coos,
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    NormalModesProxy& mode);
extern "C" void fortran_integration_timer_ele(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    void* start /* 0D_NOT_type in */,
    void* orb_max /* 0D_NOT_type in */,
    double& tol /* 0D_NOT_real inout */);
void integration_timer_ele(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start,
    CoordProxy& orb_max,
    double& tol);

// Skipped unusable routine integration_timer_fibre:
// - Untranslated type: fibre (0D)

// Skipped unusable routine interpolate_field:
// - Untranslated type: mesh3d_struct (0D)
extern "C" void fortran_ion_kick(
    void* orbit /* 0D_NOT_type in */,
    double* r_beam /* 1D_NOT_real in */,
    double& n_beam_part /* 0D_NOT_real in */,
    void* a_twiss /* 0D_NOT_type in */,
    void* b_twiss /* 0D_NOT_type in */,
    double& sig_ee /* 0D_NOT_real in */,
    double* kick /* 1D_NOT_real out */);
FixedArray1D<Real, 3> ion_kick(
    CoordProxy& orbit,
    FixedArray1D<Real, 2> r_beam,
    double n_beam_part,
    TwissProxy& a_twiss,
    TwissProxy& b_twiss,
    double sig_ee);
extern "C" bool fortran_is_attribute(
    int& ix_attrib /* 0D_NOT_integer in */,
    int& which /* 0D_NOT_integer in */,
    bool& is_attrib /* 0D_NOT_logical out */);
bool is_attribute(int ix_attrib, int which);

// Skipped unusable routine jac:
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_key_name_to_key_index(
    const char* key_str /* 0D_NOT_character in */,
    bool* abbrev_allowed /* 0D_NOT_logical in */,
    int& key_index /* 0D_NOT_integer inout */);
void key_name_to_key_index(
    std::string key_str,
    std::optional<bool> abbrev_allowed,
    int& key_index);
extern "C" void fortran_kick_vector_calc(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& s_body /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type in */,
    double* dr_ds /* 1D_NOT_real out */,
    bool& err /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical inout */);
struct KickVectorCalc {
  FixedArray1D<Real, 11> dr_ds;
  bool err;
};
Bmad::KickVectorCalc kick_vector_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_body,
    CoordProxy& orbit,
    optional_ref<bool> print_err = std::nullopt);
extern "C" void fortran_kill_complex_taylor(
    void* complex_taylor /* 1D_ALLOC_type inout */);
void kill_complex_taylor(ComplexTaylorProxyAlloc1D& complex_taylor);
extern "C" void fortran_kill_ptc_layouts(void* lat /* 0D_NOT_type in */);
void kill_ptc_layouts(LatProxy& lat);
extern "C" void fortran_kill_taylor(
    void* bmad_taylor /* 1D_ALLOC_type inout */);
void kill_taylor(TaylorProxyAlloc1D& bmad_taylor);
extern "C" bool fortran_kind_name(
    int* this_kind /* 0D_PTR_integer in */,
    const char* kind_str /* 0D_NOT_character out */);
std::string kind_name(int this_kind);
extern "C" bool fortran_knot_interpolate(
    void* x_knot /* 1D_ALLOC_real in */,
    void* y_knot /* 1D_ALLOC_real in */,
    double& x_pt /* 0D_NOT_real in */,
    int& interpolation /* 0D_NOT_integer in */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& y_pt /* 0D_NOT_real inout */);
bool knot_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x_pt,
    int interpolation,
    double& y_pt);
extern "C" bool fortran_knots_to_string(
    void* x_knot /* 1D_ALLOC_real inout */,
    void* y_knot /* 1D_ALLOC_real inout */,
    const char* str /* 0D_ALLOC_character inout */);
void knots_to_string(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    std::string& str);

// Skipped unusable routine kubo_integrand:
// - Untranslated type: c_ptr (0D)
extern "C" bool fortran_lafun(
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void lafun(double& x, double& y, double& z, double& res);
extern "C" void fortran_lat_compute_ref_energy_and_time(
    void* lat /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool lat_compute_ref_energy_and_time(LatProxy& lat);
extern "C" void fortran_lat_ele_locator(
    const char* loc_str /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    void* eles /* 1D_ALLOC_type inout */,
    int& n_loc /* 0D_NOT_integer inout */,
    bool& err /* 0D_NOT_logical out */,
    bool* above_ubound_is_err /* 0D_NOT_logical in */,
    int* ix_dflt_branch /* 0D_NOT_integer in */,
    bool* order_by_index /* 0D_NOT_logical in */,
    bool* append_eles /* 0D_NOT_logical in */);
bool lat_ele_locator(
    std::string loc_str,
    LatProxy& lat,
    ElePointerProxyAlloc1D& eles,
    int& n_loc,
    std::optional<bool> above_ubound_is_err = std::nullopt,
    std::optional<int> ix_dflt_branch = std::nullopt,
    std::optional<bool> order_by_index = std::nullopt,
    std::optional<bool> append_eles = std::nullopt);
extern "C" void fortran_lat_equal_lat(
    void* lat_out /* 0D_NOT_type out */,
    void* lat_in /* 0D_NOT_type in */);
LatProxy lat_equal_lat(LatProxy& lat_in);
extern "C" void fortran_lat_geometry(void* lat /* 0D_NOT_type inout */);
void lat_geometry(LatProxy& lat);
extern "C" void fortran_lat_make_mat6(
    void* lat /* 0D_NOT_type inout */,
    int* ix_ele /* 0D_NOT_integer in */,
    void* ref_orb /* 1D_ALLOC_type in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool lat_make_mat6(
    LatProxy& lat,
    std::optional<int> ix_ele = std::nullopt,
    optional_ref<CoordProxyAlloc1D> ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine lat_make_mat6_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_lat_sanity_check(
    void* lat /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool lat_sanity_check(LatProxy& lat);
extern "C" void fortran_lat_to_ptc_layout(void* lat /* 0D_NOT_type in */);
void lat_to_ptc_layout(LatProxy& lat);
extern "C" void fortran_lat_vec_equal_lat_vec(
    void* lat1 /* 1D_ALLOC_type out */,
    void* lat2 /* 1D_ALLOC_type in */);
LatProxyAlloc1D lat_vec_equal_lat_vec(LatProxyAlloc1D& lat2);
extern "C" void fortran_lattice_bookkeeper(
    void* lat /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool lattice_bookkeeper(LatProxy& lat);
extern "C" void fortran_lcavity_rf_step_setup(
    void* ele /* 0D_NOT_type inout */);
void lcavity_rf_step_setup(EleProxy& ele);
extern "C" void fortran_linear_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orb /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void linear_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_linear_coef(
    void* stack /* 1D_ALLOC_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& coef /* 0D_NOT_real out */);
struct LinearCoef {
  bool err_flag;
  double coef;
};
Bmad::LinearCoef linear_coef(ExpressionAtomProxyAlloc1D& stack);
extern "C" void fortran_linear_to_spin_taylor(
    double* q_map /* 2D_NOT_real in */,
    void* spin_taylor /* 1D_NOT_type out */);
TaylorProxyArray1D linear_to_spin_taylor(FixedArray2D<Real, 4, 7> q_map);
extern "C" void fortran_load_parse_line(
    const char* action /* 0D_NOT_character in */,
    int& ix_start /* 0D_NOT_integer in */,
    bool& end_of_file /* 0D_NOT_logical out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct LoadParseLine {
  bool end_of_file;
  bool err_flag;
};
Bmad::LoadParseLine load_parse_line(std::string action, int ix_start);
extern "C" bool fortran_lord_edge_aligned(
    void* slave /* 0D_NOT_type in */,
    int& slave_edge /* 0D_NOT_integer in */,
    void* lord /* 0D_NOT_type in */,
    bool& is_aligned /* 0D_NOT_logical inout */);
void lord_edge_aligned(
    EleProxy& slave,
    int slave_edge,
    EleProxy& lord,
    bool& is_aligned);
extern "C" bool fortran_low_energy_z_correction(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double& ds /* 0D_NOT_real in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */,
    double& dz /* 0D_NOT_real inout */);
void low_energy_z_correction(
    CoordProxy& orbit,
    EleProxy& ele,
    double ds,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    double& dz);

// Skipped unusable routine lsc_kick_params_calc:
// - Untranslated type: csr_struct (0D)
extern "C" void fortran_mad_add_offsets_and_multipoles(
    void* ele /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_add_offsets_and_multipoles(EleProxy& ele);
extern "C" void fortran_mad_concat_map2(
    void* map1 /* 0D_NOT_type in */,
    void* map2 /* 0D_NOT_type in */,
    void* map3 /* 0D_NOT_type out */);
MadMapProxy mad_concat_map2(MadMapProxy& map1, MadMapProxy& map2);
extern "C" void fortran_mad_drift(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_drift(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_elsep(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_elsep(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_map_to_taylor(
    void* map /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* taylor /* 1D_ALLOC_type out */);
TaylorProxyAlloc1D mad_map_to_taylor(MadMapProxy& map, MadEnergyProxy& energy);
extern "C" void fortran_mad_quadrupole(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_quadrupole(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_rfcavity(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_rfcavity(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_sbend(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_sbend(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_sbend_body(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_sbend_body(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_sbend_fringe(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    bool& into /* 0D_NOT_logical in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_sbend_fringe(EleProxy& ele, MadEnergyProxy& energy, bool into);
extern "C" void fortran_mad_sextupole(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_sextupole(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_solenoid(
    void* ele /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy mad_solenoid(EleProxy& ele, MadEnergyProxy& energy);
extern "C" void fortran_mad_tmfoc(
    double& el /* 0D_NOT_real in */,
    double& sk1 /* 0D_NOT_real in */,
    double& c /* 0D_NOT_real out */,
    double& s /* 0D_NOT_real out */,
    double& d /* 0D_NOT_real out */,
    double& f /* 0D_NOT_real out */);
struct MadTmfoc {
  double c;
  double s;
  double d;
  double f;
};
Bmad::MadTmfoc mad_tmfoc(double el, double sk1);
extern "C" void fortran_mad_tmsymm(double* te /* 3D_NOT_real inout */);
void mad_tmsymm(FixedArray3D<Real, 6, 6, 6> te);
extern "C" void fortran_mad_tmtilt(
    void* map /* 0D_NOT_type inout */,
    double& tilt /* 0D_NOT_real in */);
void mad_tmtilt(MadMapProxy& map, double tilt);
extern "C" void fortran_mad_track1(
    void* c0 /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type in */,
    void* c1 /* 0D_NOT_type out */);
CoordProxy mad_track1(CoordProxy& c0, MadMapProxy& map);
extern "C" void fortran_make_g2_mats(
    void* twiss /* 0D_NOT_type in */,
    double* g2_mat /* 2D_NOT_real inout */,
    double* g2_inv_mat /* 2D_NOT_real inout */);
void make_g2_mats(
    TwissProxy& twiss,
    FixedArray2D<Real, 2, 2> g2_mat,
    FixedArray2D<Real, 2, 2> g2_inv_mat);
extern "C" void fortran_make_g_mats(
    void* ele /* 0D_NOT_type in */,
    double* g_mat /* 2D_NOT_real out */,
    double* g_inv_mat /* 2D_NOT_real out */);
struct MakeGMats {
  FixedArray2D<Real, 4, 4> g_mat;
  FixedArray2D<Real, 4, 4> g_inv_mat;
};
Bmad::MakeGMats make_g_mats(EleProxy& ele);
extern "C" void fortran_make_hvbp(
    double* N /* 2D_NOT_real in */,
    double* B /* 2D_NOT_real out */,
    double* V /* 2D_NOT_real out */,
    double* H /* 2D_NOT_real out */,
    double* Vbar /* 2D_NOT_real out */,
    double* Hbar /* 2D_NOT_real out */);
struct MakeHvbp {
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> V;
  FixedArray2D<Real, 6, 6> H;
  std::optional<FixedArray2D<Real, 6, 6>> Vbar;
  std::optional<FixedArray2D<Real, 6, 6>> Hbar;
};
Bmad::MakeHvbp make_hvbp(FixedArray2D<Real, 6, 6> N);
extern "C" void fortran_make_hybrid_lat(
    void* lat_in /* 0D_NOT_type in */,
    void* lat_out /* 0D_NOT_type out */,
    bool* use_taylor /* 0D_NOT_logical in */,
    void* orb0_arr /* 1D_ALLOC_type in */);
LatProxy make_hybrid_lat(
    LatProxy& lat_in,
    std::optional<bool> use_taylor = std::nullopt,
    optional_ref<CoordArrayProxyAlloc1D> orb0_arr = std::nullopt);
extern "C" void fortran_make_mad_map(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* energy /* 0D_NOT_type out */,
    void* map /* 0D_NOT_type out */);
struct MakeMadMap {
  MadEnergyProxy energy;
  MadMapProxy map;
};
Bmad::MakeMadMap make_mad_map(EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_make_mat6(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct MakeMat6 {
  CoordProxy end_orb;
  bool err_flag;
};
Bmad::MakeMat6 make_mat6(
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> start_orb = std::nullopt);
extern "C" void fortran_make_mat6_bmad(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    bool& err /* 0D_NOT_logical out */);
struct MakeMat6Bmad {
  CoordProxy end_orb;
  bool err;
};
Bmad::MakeMat6Bmad make_mat6_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb);
extern "C" void fortran_make_mat6_bmad_photon(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    bool& err /* 0D_NOT_logical out */);
struct MakeMat6BmadPhoton {
  CoordProxy end_orb;
  bool err;
};
Bmad::MakeMat6BmadPhoton make_mat6_bmad_photon(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb);

// Skipped unusable routine make_mat6_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_make_mat6_high_energy_space_charge(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */);
void make_mat6_high_energy_space_charge(EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_make_mat6_mad(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* c0 /* 0D_NOT_type in */,
    void* c1 /* 0D_NOT_type out */);
CoordProxy make_mat6_mad(EleProxy& ele, LatParamProxy& param, CoordProxy& c0);
extern "C" void fortran_make_mat6_symp_lie_ptc(
    void* ele /* 0D_NOT_type inout */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */);
CoordProxy make_mat6_symp_lie_ptc(EleProxy& ele, CoordProxy& start_orb);
extern "C" void fortran_make_mat6_taylor(
    void* ele /* 0D_NOT_type inout */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    bool* err_flag /* 0D_NOT_logical inout */);
CoordProxy make_mat6_taylor(
    EleProxy& ele,
    CoordProxy& start_orb,
    optional_ref<bool> err_flag = std::nullopt);
extern "C" void fortran_make_mat6_tracking(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* spin_only /* 0D_NOT_logical in */);
struct MakeMat6Tracking {
  CoordProxy end_orb;
  bool err_flag;
};
Bmad::MakeMat6Tracking make_mat6_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb,
    std::optional<bool> spin_only = std::nullopt);
extern "C" void fortran_make_n(
    double* t6 /* 2D_NOT_real in */,
    double* N /* 2D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    double* abz_tunes /* 1D_NOT_real in */,
    double* tunes_out /* 1D_NOT_real out */,
    double* U /* 2D_NOT_real out */);
struct MakeN {
  FixedArray2D<Real, 6, 6> N;
  bool err_flag;
  FixedArray1D<Real, 3> tunes_out;
  std::optional<FixedArray2D<Real, 6, 6>> U;
};
Bmad::MakeN make_n(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt);
extern "C" void fortran_make_pbrh(
    double* M /* 2D_NOT_real in */,
    std::complex<double>* P /* 2D_NOT_complex out */,
    std::complex<double>* Bp /* 2D_NOT_complex out */,
    std::complex<double>* R /* 2D_NOT_complex out */,
    std::complex<double>* H /* 2D_NOT_complex out */,
    double* abz_tunes /* 1D_NOT_real in */);
struct MakePbrh {
  FixedArray2D<Complex, 6, 6> P;
  FixedArray2D<Complex, 6, 6> Bp;
  FixedArray2D<Complex, 6, 6> R;
  FixedArray2D<Complex, 6, 6> H;
};
Bmad::MakePbrh make_pbrh(
    FixedArray2D<Real, 6, 6> M,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_make_smat_from_abc(
    double* t6 /* 2D_NOT_real in */,
    void* mode /* 0D_NOT_type in */,
    double* sigma_mat /* 2D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    double* Nout /* 2D_NOT_real out */);
struct MakeSmatFromAbc {
  FixedArray2D<Real, 6, 6> sigma_mat;
  bool err_flag;
  std::optional<FixedArray2D<Real, 6, 6>> Nout;
};
Bmad::MakeSmatFromAbc make_smat_from_abc(
    FixedArray2D<Real, 6, 6> t6,
    NormalModesProxy& mode);

// Skipped unusable routine make_sr_mats:
// - Variable out sized array: M(:,:) 2D_NOT_real
// - Variable out sized array: Bone(:,:) 2D_NOT_real
// - Variable out sized array: Done(:,:) 2D_NOT_real

// Skipped unusable routine make_srdt_cache:
// - Untranslated type: sliced_eles_struct (1D)
// - Variable inout sized array: cache(:,:,:) 3D_ALLOC_complex
extern "C" void fortran_make_unit_mad_map(void* map /* 0D_NOT_type inout */);
void make_unit_mad_map(MadMapProxy& map);
extern "C" void fortran_make_v(
    double* M /* 2D_NOT_real inout */,
    std::complex<double>* V /* 2D_NOT_complex inout */,
    double* abz_tunes /* 1D_NOT_real inout */);
void make_v(
    FixedArray2D<Real, 6, 6> M,
    FixedArray2D<Complex, 6, 6> V,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_make_v_mats(
    void* ele /* 0D_NOT_type in */,
    double* v_mat /* 2D_NOT_real out */,
    double* v_inv_mat /* 2D_NOT_real out */);
struct MakeVMats {
  std::optional<FixedArray2D<Real, 4, 4>> v_mat;
  std::optional<FixedArray2D<Real, 4, 4>> v_inv_mat;
};
Bmad::MakeVMats make_v_mats(EleProxy& ele);

// Skipped unusable routine make_ykick_mat:
// - Variable inout sized array: Yone(:,:) 2D_NOT_real
extern "C" void fortran_makeup_control_slave(
    void* lat /* 0D_NOT_type inout */,
    void* slave /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void makeup_control_slave(LatProxy& lat, EleProxy& slave, bool& err_flag);
extern "C" void fortran_makeup_group_lord(
    void* lat /* 0D_NOT_type inout */,
    void* lord /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void makeup_group_lord(LatProxy& lat, EleProxy& lord, bool& err_flag);
extern "C" void fortran_makeup_multipass_slave(
    void* lat /* 0D_NOT_type inout */,
    void* slave /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void makeup_multipass_slave(LatProxy& lat, EleProxy& slave, bool& err_flag);
extern "C" void fortran_makeup_super_slave(
    void* lat /* 0D_NOT_type inout */,
    void* slave /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void makeup_super_slave(LatProxy& lat, EleProxy& slave, bool& err_flag);
extern "C" void fortran_makeup_super_slave1(
    void* slave /* 0D_NOT_type inout */,
    void* lord /* 0D_NOT_type in */,
    double& offset /* 0D_NOT_real in */,
    void* param /* 0D_NOT_type in */,
    bool& include_upstream_end /* 0D_NOT_logical in */,
    bool& include_downstream_end /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool makeup_super_slave1(
    EleProxy& slave,
    EleProxy& lord,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end);
extern "C" bool fortran_map1_inverse(
    void* map1 /* 0D_NOT_type in */,
    void* inv_map1 /* 0D_NOT_type inout */);
void map1_inverse(SpinOrbitMap1Proxy& map1, SpinOrbitMap1Proxy& inv_map1);
extern "C" void fortran_map1_make_unit(void* map1 /* 0D_NOT_type out */);
SpinOrbitMap1Proxy map1_make_unit();
extern "C" bool fortran_map1_times_map1(
    void* map2 /* 0D_NOT_type in */,
    void* map1 /* 0D_NOT_type in */,
    void* map_out /* 0D_NOT_type out */);
SpinOrbitMap1Proxy map1_times_map1(
    SpinOrbitMap1Proxy& map2,
    SpinOrbitMap1Proxy& map1);

// Skipped unusable routine map_coef:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_map_to_angle_coords(
    void* t_canon /* 1D_NOT_type in */,
    void* t_angle /* 1D_NOT_type out */);
TaylorProxyArray1D map_to_angle_coords(FixedArray1D<TaylorProxy, 6> t_canon);
extern "C" void fortran_mark_patch_regions(
    void* branch /* 0D_NOT_type inout */);
void mark_patch_regions(BranchProxy& branch);
extern "C" bool fortran_master_parameter_value(
    int& master_parameter /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    double& value /* 0D_NOT_real inout */);
void master_parameter_value(int master_parameter, EleProxy& ele, double& value);
extern "C" void fortran_mat4_multipole(
    double& knl /* 0D_NOT_real in */,
    double& tilt /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer inout */,
    void* orbit /* 0D_NOT_type in */,
    double* kick_mat /* 2D_NOT_real out */);
FixedArray2D<Real, 4, 4> mat4_multipole(
    double knl,
    double tilt,
    int& n,
    CoordProxy& orbit);
extern "C" void fortran_mat6_add_offsets(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */);
void mat6_add_offsets(EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_mat6_add_pitch(
    double& x_pitch_tot /* 0D_NOT_real in */,
    double& y_pitch_tot /* 0D_NOT_real in */,
    int& orientation /* 0D_NOT_integer in */,
    double* mat6 /* 2D_NOT_real inout */);
void mat6_add_pitch(
    double x_pitch_tot,
    double y_pitch_tot,
    int orientation,
    FixedArray2D<Real, 6, 6> mat6);

// Skipped unusable routine mat6_from_s_to_s:
// - Variable inout sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_mat6_to_complex_taylor(
    std::complex<double>* vec0 /* 1D_NOT_complex in */,
    std::complex<double>* mat6 /* 2D_NOT_complex in */,
    void* complex_taylor /* 1D_NOT_type out */);
ComplexTaylorProxyArray1D mat6_to_complex_taylor(
    FixedArray1D<Complex, 6> vec0,
    FixedArray2D<Complex, 6, 6> mat6);
extern "C" void fortran_mat_symp_decouple(
    double* t0 /* 2D_NOT_real in */,
    int& stat /* 0D_NOT_integer out */,
    double* U /* 2D_NOT_real inout */,
    double* V /* 2D_NOT_real inout */,
    double* Ubar /* 2D_NOT_real inout */,
    double* Vbar /* 2D_NOT_real inout */,
    double* G /* 2D_NOT_real inout */,
    void* twiss1 /* 0D_NOT_type out */,
    void* twiss2 /* 0D_NOT_type out */,
    double& gamma /* 0D_NOT_real out */,
    bool& type_out /* 0D_NOT_logical in */);
struct MatSympDecouple {
  int stat;
  TwissProxy twiss1;
  TwissProxy twiss2;
  double gamma;
};
Bmad::MatSympDecouple mat_symp_decouple(
    FixedArray2D<Real, 4, 4> t0,
    FixedArray2D<Real, 4, 4> U,
    FixedArray2D<Real, 4, 4> V,
    FixedArray2D<Real, 4, 4> Ubar,
    FixedArray2D<Real, 4, 4> Vbar,
    FixedArray2D<Real, 4, 4> G,
    bool type_out);
extern "C" void fortran_match_ele_to_mat6(
    void* ele /* 0D_NOT_type in */,
    void* start_orb /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    double* vec0 /* 1D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* include_delta_time /* 0D_NOT_logical in */,
    bool* set_trombone /* 0D_NOT_logical in */);
struct MatchEleToMat6 {
  FixedArray2D<Real, 6, 6> mat6;
  FixedArray1D<Real, 6> vec0;
  bool err_flag;
};
Bmad::MatchEleToMat6 match_ele_to_mat6(
    EleProxy& ele,
    CoordProxy& start_orb,
    std::optional<bool> include_delta_time = std::nullopt,
    std::optional<bool> set_trombone = std::nullopt);

// Skipped unusable routine mccfft1d:
// - Routine module (fft_interface_mod) in configuration skip list
extern "C" bool fortran_mexp(
    double& x /* 0D_NOT_real in */,
    int& m /* 0D_NOT_integer in */,
    double& this_exp /* 0D_NOT_real inout */);
void mexp(double x, int m, double& this_exp);
extern "C" void fortran_mfft1(
    void* a /* 1D_ALLOC_real inout */,
    void* b /* 1D_ALLOC_real inout */,
    void* n /* 1D_ALLOC_integer in */,
    int& ndim /* 0D_NOT_integer in */,
    int& isn /* 0D_NOT_integer in */,
    int& ierr /* 0D_NOT_integer out */);
int mfft1(RealAlloc1D& a, RealAlloc1D& b, IntAlloc1D& n, int ndim, int isn);

// Skipped unusable routine misalign_ptc_fibre:
// - Untranslated type: fibre (0D)
extern "C" bool fortran_momentum_compaction(
    void* branch /* 0D_NOT_type in */,
    double& mom_comp /* 0D_NOT_real inout */);
void momentum_compaction(BranchProxy& branch, double& mom_comp);

// Skipped unusable routine mpxx1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine mpxx_integrand:
// - Untranslated type: c_ptr (0D)

// Skipped unusable routine mpzt1:
// - Untranslated type: ibs_struct (0D)

// Skipped unusable routine multi_coulomb_log:
// - Untranslated type: ibs_sim_param_struct (0D)
extern "C" void fortran_multi_turn_tracking_analysis(
    void* track /* 1D_ALLOC_type in */,
    int& i_dim /* 0D_NOT_integer in */,
    void* track0 /* 0D_NOT_type out */,
    void* ele /* 0D_NOT_type out */,
    bool& stable /* 0D_NOT_logical out */,
    double& growth_rate /* 0D_NOT_real out */,
    double& chi /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct MultiTurnTrackingAnalysis {
  CoordProxy track0;
  EleProxy ele;
  bool stable;
  double growth_rate;
  double chi;
  bool err_flag;
};
Bmad::MultiTurnTrackingAnalysis multi_turn_tracking_analysis(
    CoordProxyAlloc1D& track,
    int i_dim);

// Skipped unusable routine multi_turn_tracking_to_mat:
// - Variable out sized array: map1(:,:) 2D_NOT_real
extern "C" void fortran_multilayer_type_to_multilayer_params(
    void* ele /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool multilayer_type_to_multilayer_params(EleProxy& ele);

// Skipped unusable routine multipass_all_info:
// - Untranslated type: multipass_all_info_struct (0D)
extern "C" void fortran_multipass_chain(
    void* ele /* 0D_NOT_type in */,
    int& ix_pass /* 0D_NOT_integer in */,
    int& n_links /* 0D_NOT_integer in */,
    void* chain_ele /* 1D_ALLOC_type in */,
    bool* use_super_lord /* 0D_NOT_logical in */);
void multipass_chain(
    EleProxy& ele,
    int ix_pass,
    int n_links,
    optional_ref<ElePointerProxyAlloc1D> chain_ele = std::nullopt,
    std::optional<bool> use_super_lord = std::nullopt);

// Skipped unusable routine multipass_region_info:
// - Untranslated type: multipass_region_lat_struct (0D)
// - Untranslated type: multipass_all_info_struct (0D)
extern "C" void fortran_multipole1_ab_to_kt(
    double& an /* 0D_NOT_real in */,
    double& bn /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer in */,
    double& knl /* 0D_NOT_real out */,
    double& tn /* 0D_NOT_real out */);
struct Multipole1AbToKt {
  double knl;
  double tn;
};
Bmad::Multipole1AbToKt multipole1_ab_to_kt(double an, double bn, int n);
extern "C" void fortran_multipole1_kt_to_ab(
    double& knl /* 0D_NOT_real in */,
    double& knsl /* 0D_NOT_real in */,
    double& tn /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer in */,
    double& an /* 0D_NOT_real out */,
    double& bn /* 0D_NOT_real out */);
struct Multipole1KtToAb {
  double an;
  double bn;
};
Bmad::Multipole1KtToAb multipole1_kt_to_ab(
    double knl,
    double knsl,
    double tn,
    int n);
extern "C" void fortran_multipole_ab_to_kt(
    void* an /* 1D_ALLOC_real in */,
    void* bn /* 1D_ALLOC_real in */,
    void* knl /* 1D_ALLOC_real out */,
    void* tn /* 1D_ALLOC_real out */);
struct MultipoleAbToKt {
  RealAlloc1D knl;
  RealAlloc1D tn;
};
Bmad::MultipoleAbToKt multipole_ab_to_kt(RealAlloc1D& an, RealAlloc1D& bn);
extern "C" void fortran_multipole_ele_to_ab(
    void* ele /* 0D_NOT_type in */,
    bool& use_ele_tilt /* 0D_NOT_logical in */,
    int& ix_pole_max /* 0D_NOT_integer out */,
    double* a /* 1D_NOT_real out */,
    double* b /* 1D_NOT_real out */,
    int* pole_type /* 0D_NOT_integer in */,
    int* include_kicks /* 0D_NOT_integer in */,
    double& b1 /* 0D_NOT_real out */,
    bool* original /* 0D_NOT_logical in */);
struct MultipoleEleToAb {
  int ix_pole_max;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> a;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> b;
  double b1;
};
Bmad::MultipoleEleToAb multipole_ele_to_ab(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt,
    std::optional<bool> original = std::nullopt);
extern "C" void fortran_multipole_ele_to_kt(
    void* ele /* 0D_NOT_type in */,
    bool& use_ele_tilt /* 0D_NOT_logical in */,
    int& ix_pole_max /* 0D_NOT_integer out */,
    void* knl /* 1D_ALLOC_real out */,
    void* tilt /* 1D_ALLOC_real out */,
    int* pole_type /* 0D_NOT_integer in */,
    int* include_kicks /* 0D_NOT_integer in */);
struct MultipoleEleToKt {
  int ix_pole_max;
  RealAlloc1D knl;
  RealAlloc1D tilt;
};
Bmad::MultipoleEleToKt multipole_ele_to_kt(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt);
extern "C" void fortran_multipole_init(
    void* ele /* 0D_NOT_type out */,
    int& who /* 0D_NOT_integer in */,
    bool* zero /* 0D_NOT_logical in */);
EleProxy multipole_init(int who, std::optional<bool> zero = std::nullopt);
extern "C" void fortran_multipole_kick(
    double& knl /* 0D_NOT_real in */,
    double& tilt /* 0D_NOT_real in */,
    int& n /* 0D_NOT_integer in */,
    int& ref_species /* 0D_NOT_integer in */,
    int& ele_orientation /* 0D_NOT_integer in */,
    void* coord /* 0D_NOT_type inout */,
    int* pole_type /* 0D_NOT_integer in */,
    bool* ref_orb_offset /* 0D_NOT_logical in */);
void multipole_kick(
    double knl,
    double tilt,
    int n,
    int ref_species,
    int ele_orientation,
    CoordProxy& coord,
    std::optional<int> pole_type = std::nullopt,
    std::optional<bool> ref_orb_offset = std::nullopt);
extern "C" void fortran_multipole_kick_mat(
    void* knl /* 1D_ALLOC_real in */,
    void* tilt /* 1D_ALLOC_real in */,
    int& ref_species /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    double& factor /* 0D_NOT_real in */,
    double* mat6 /* 2D_NOT_real out */);
FixedArray2D<Real, 6, 6> multipole_kick_mat(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    int ref_species,
    EleProxy& ele,
    CoordProxy& orbit,
    double factor);
extern "C" void fortran_multipole_kicks(
    void* knl /* 1D_ALLOC_real in */,
    void* tilt /* 1D_ALLOC_real in */,
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    int* pole_type /* 0D_NOT_integer in */,
    bool* ref_orb_offset /* 0D_NOT_logical in */);
void multipole_kicks(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<int> pole_type = std::nullopt,
    std::optional<bool> ref_orb_offset = std::nullopt);
extern "C" void fortran_multipole_kt_to_ab(
    void* knl /* 1D_ALLOC_real in */,
    void* knsl /* 1D_ALLOC_real in */,
    void* tn /* 1D_ALLOC_real in */,
    void* an /* 1D_ALLOC_real out */,
    void* bn /* 1D_ALLOC_real out */);
struct MultipoleKtToAb {
  RealAlloc1D an;
  RealAlloc1D bn;
};
Bmad::MultipoleKtToAb multipole_kt_to_ab(
    RealAlloc1D& knl,
    RealAlloc1D& knsl,
    RealAlloc1D& tn);
extern "C" void fortran_multipole_spin_tracking(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void multipole_spin_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" bool fortran_mytan(
    double& y /* 0D_NOT_real inout */,
    double& x /* 0D_NOT_real inout */,
    double& arg /* 0D_NOT_real inout */);
void mytan(double& y, double& x, double& arg);
extern "C" bool fortran_n_attrib_string_max_len(
    int& max_len /* 0D_NOT_integer out */);
int n_attrib_string_max_len();
extern "C" void fortran_new_control(
    void* lat /* 0D_NOT_type in */,
    int& ix_ele /* 0D_NOT_integer in */,
    const char* ele_name /* 0D_NOT_character in */);
void new_control(
    LatProxy& lat,
    int ix_ele,
    std::optional<std::string> ele_name = std::nullopt);
extern "C" bool fortran_nint_chk(
    double& re_val /* 0D_NOT_real in */,
    int& int_val /* 0D_NOT_integer out */);
int nint_chk(double re_val);
extern "C" void fortran_normal_form_complex_taylors(
    void* one_turn_taylor /* 1D_NOT_type inout */,
    bool& rf_on /* 0D_NOT_logical inout */,
    void* F /* 1D_NOT_type inout */,
    void* L /* 1D_NOT_type inout */,
    void* A /* 1D_NOT_type inout */,
    void* A_inverse /* 1D_NOT_type inout */,
    int* order /* 0D_NOT_integer inout */);
void normal_form_complex_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool& rf_on,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> F = std::nullopt,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> L = std::nullopt,
    std::optional<FixedArray1D<TaylorProxy, 6>> A = std::nullopt,
    std::optional<FixedArray1D<TaylorProxy, 6>> A_inverse = std::nullopt,
    optional_ref<int> order = std::nullopt);

// Skipped unusable routine normal_form_rd_terms:
// - Untranslated type: probe_8 (0D)
extern "C" void fortran_normal_form_taylors(
    void* one_turn_taylor /* 1D_NOT_type in */,
    bool& rf_on /* 0D_NOT_logical in */,
    void* dhdj /* 1D_NOT_type out */,
    void* A /* 1D_NOT_type out */,
    void* A_inverse /* 1D_NOT_type out */);
struct NormalFormTaylors {
  TaylorProxyArray1D dhdj;
  TaylorProxyArray1D A;
  TaylorProxyArray1D A_inverse;
};
Bmad::NormalFormTaylors normal_form_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool rf_on);
extern "C" void fortran_normal_mode3_calc(
    double* t6 /* 2D_NOT_real inout */,
    double* tune /* 1D_NOT_real out */,
    double* B /* 2D_NOT_real out */,
    double* HV /* 2D_NOT_real out */,
    bool* above_transition /* 0D_NOT_logical in */,
    double* abz_tunes /* 1D_NOT_real in */);
struct NormalMode3Calc {
  FixedArray1D<Real, 3> tune;
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> HV;
};
Bmad::NormalMode3Calc normal_mode3_calc(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<bool> above_transition = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt);
extern "C" void fortran_normal_mode_dispersion(
    void* ele /* 0D_NOT_type inout */,
    bool* reverse /* 0D_NOT_logical in */);
void normal_mode_dispersion(
    EleProxy& ele,
    std::optional<bool> reverse = std::nullopt);
extern "C" void fortran_normalize_evecs(
    std::complex<double>* evec /* 2D_NOT_complex inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool normalize_evecs(FixedArray2D<Complex, 6, 6> evec);
extern "C" bool fortran_num_field_eles(
    void* ele /* 0D_NOT_type in */,
    int& n_field_ele /* 0D_NOT_integer inout */);
void num_field_eles(EleProxy& ele, int& n_field_ele);
extern "C" bool fortran_num_lords(
    void* slave /* 0D_NOT_type in */,
    int& lord_type /* 0D_NOT_integer in */,
    int& num /* 0D_NOT_integer inout */);
void num_lords(EleProxy& slave, int lord_type, int& num);
extern "C" void fortran_odeint_bmad(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& s1_body /* 0D_NOT_real in */,
    double& s2_body /* 0D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type out */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
struct OdeintBmad {
  bool err_flag;
  TrackProxy track;
};
Bmad::OdeintBmad odeint_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    double s1_body,
    double s2_body,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_odeint_bmad_time(
    void* orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& t_dir /* 0D_NOT_integer in */,
    double& rf_time /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type inout */,
    double* t_end /* 0D_NOT_real in */,
    double& dt_step /* 0D_NOT_real out */,
    void* extra_field /* 0D_NOT_type in */);
struct OdeintBmadTime {
  bool err_flag;
  double dt_step;
};
Bmad::OdeintBmadTime odeint_bmad_time(
    CoordProxy& orb,
    EleProxy& ele,
    LatParamProxy& param,
    int t_dir,
    double& rf_time,
    optional_ref<TrackProxy> track = std::nullopt,
    std::optional<double> t_end = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" void fortran_offset_particle(
    void* ele /* 0D_NOT_type in */,
    bool& set /* 0D_NOT_logical in */,
    void* orbit /* 0D_NOT_type inout */,
    bool* set_tilt /* 0D_NOT_logical in */,
    bool* set_hvkicks /* 0D_NOT_logical in */,
    int* drift_to_edge /* 0D_NOT_integer in */,
    double* s_pos /* 0D_NOT_real in */,
    double& s_out /* 0D_NOT_real out */,
    bool* set_spin /* 0D_NOT_logical in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */,
    double* spin_qrot /* 1D_NOT_real out */,
    double* time /* 0D_NOT_real inout */);
struct OffsetParticle {
  double s_out;
  FixedArray1D<Real, 4> spin_qrot;
};
Bmad::OffsetParticle offset_particle(
    EleProxy& ele,
    bool set,
    CoordProxy& orbit,
    std::optional<bool> set_tilt = std::nullopt,
    std::optional<bool> set_hvkicks = std::nullopt,
    std::optional<int> drift_to_edge = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> set_spin = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    optional_ref<double> time = std::nullopt);
extern "C" void fortran_offset_photon(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    bool& set /* 0D_NOT_logical in */,
    bool* offset_position_only /* 0D_NOT_logical in */,
    double* rot_mat /* 2D_NOT_real in */);
void offset_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    std::optional<bool> offset_position_only = std::nullopt,
    std::optional<FixedArray2D<Real, 3, 3>> rot_mat = std::nullopt);
extern "C" void fortran_one_turn_mat_at_ele(
    void* ele /* 0D_NOT_type in */,
    double& phi_a /* 0D_NOT_real in */,
    double& phi_b /* 0D_NOT_real in */,
    double* mat4 /* 2D_NOT_real out */);
FixedArray2D<Real, 4, 4> one_turn_mat_at_ele(
    EleProxy& ele,
    double phi_a,
    double phi_b);
extern "C" bool fortran_open_binary_file(
    const char* file_name /* 0D_NOT_character in */,
    const char* action /* 0D_NOT_character in */,
    int& iu /* 0D_NOT_integer out */,
    const char* r_name /* 0D_NOT_character in */,
    int& iver /* 0D_NOT_integer out */,
    bool& is_ok /* 0D_NOT_logical out */);
struct OpenBinaryFile {
  int iu;
  int iver;
  bool is_ok;
};
Bmad::OpenBinaryFile open_binary_file(
    std::string file_name,
    std::string action,
    std::string r_name);
extern "C" void fortran_orbit_amplitude_calc(
    void* ele /* 0D_NOT_type in */,
    void* orb /* 0D_NOT_type in */,
    double& amp_a /* 0D_NOT_real out */,
    double& amp_b /* 0D_NOT_real out */,
    double& amp_na /* 0D_NOT_real out */,
    double& amp_nb /* 0D_NOT_real out */);
struct OrbitAmplitudeCalc {
  double amp_a;
  double amp_b;
  double amp_na;
  double amp_nb;
};
Bmad::OrbitAmplitudeCalc orbit_amplitude_calc(EleProxy& ele, CoordProxy& orb);
extern "C" void fortran_orbit_reference_energy_correction(
    void* orbit /* 0D_NOT_type inout */,
    double& p0c_new /* 0D_NOT_real in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void orbit_reference_energy_correction(
    CoordProxy& orbit,
    double p0c_new,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_orbit_to_floor_phase_space(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double* floor_phase_space /* 1D_NOT_real inout */);
void orbit_to_floor_phase_space(
    CoordProxy& orbit,
    EleProxy& ele,
    FixedArray1D<Real, 6> floor_phase_space);
extern "C" bool fortran_orbit_to_local_curvilinear(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    int* z_direction /* 0D_NOT_integer in */,
    int* relative_to /* 0D_NOT_integer in */,
    void* local_position /* 0D_NOT_type inout */);
void orbit_to_local_curvilinear(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<int> z_direction,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position);
extern "C" bool fortran_orbit_too_large(
    void* orbit /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type out */,
    bool* check_momentum /* 0D_NOT_logical in */,
    bool& is_too_large /* 0D_NOT_logical inout */);
LatParamProxy orbit_too_large(
    CoordProxy& orbit,
    std::optional<bool> check_momentum,
    bool& is_too_large);
extern "C" void fortran_order_evecs_by_n_similarity(
    std::complex<double>* evec /* 2D_NOT_complex out */,
    std::complex<double>* eval /* 1D_NOT_complex inout */,
    double* mat_tunes /* 1D_NOT_real inout */,
    double* Nmat /* 2D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct OrderEvecsByNSimilarity {
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
};
Bmad::OrderEvecsByNSimilarity order_evecs_by_n_similarity(
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray2D<Real, 6, 6> Nmat);
extern "C" void fortran_order_evecs_by_plane_dominance(
    std::complex<double>* evec /* 2D_NOT_complex inout */,
    std::complex<double>* eval /* 1D_NOT_complex inout */,
    double* mat_tunes /* 1D_NOT_real inout */);
void order_evecs_by_plane_dominance(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    std::optional<FixedArray1D<Real, 3>> mat_tunes = std::nullopt);
extern "C" void fortran_order_evecs_by_tune(
    std::complex<double>* evec /* 2D_NOT_complex inout */,
    std::complex<double>* eval /* 1D_NOT_complex inout */,
    double* mat_tunes /* 1D_NOT_real in */,
    double* abz_tunes /* 1D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool order_evecs_by_tune(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_order_particles_in_z(
    void* bunch /* 0D_NOT_type inout */);
void order_particles_in_z(BunchProxy& bunch);
extern "C" void fortran_order_super_lord_slaves(
    void* lat /* 0D_NOT_type in */,
    int& ix_lord /* 0D_NOT_integer in */);
void order_super_lord_slaves(LatProxy& lat, int ix_lord);
extern "C" void fortran_osc_alloc_freespace_array(
    int* nlo /* 1D_NOT_integer in */,
    int* nhi /* 1D_NOT_integer in */,
    int* npad /* 1D_NOT_integer in */);
void osc_alloc_freespace_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_alloc_image_array(
    int* nlo /* 1D_NOT_integer in */,
    int* nhi /* 1D_NOT_integer in */,
    int* npad /* 1D_NOT_integer in */);
void osc_alloc_image_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_alloc_rectpipe_arrays(
    int* nlo /* 1D_NOT_integer in */,
    int* nhi /* 1D_NOT_integer in */,
    int* npad /* 1D_NOT_integer in */);
void osc_alloc_rectpipe_arrays(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);

// Skipped unusable routine osc_cathodeimages_solver:
// - Variable inout sized array: rho(:,:,:) 3D_NOT_real
// - Variable out sized array: phi(:,:,:) 3D_NOT_real
// - Variable out sized array: efield(:,:,:,:) 4D_NOT_real
// - Variable out sized array: bfield(:,:,:,:) 4D_NOT_real

// Skipped unusable routine osc_freespace_solver:
// - Variable in sized array: rho(:,:,:) 3D_NOT_real
// - Variable out sized array: phi(:,:,:) 3D_NOT_real
// - Variable out sized array: efield(:,:,:,:) 4D_NOT_real
// - Variable out sized array: bfield(:,:,:,:) 4D_NOT_real

// Skipped unusable routine osc_freespace_solver2:
// - Variable in sized array: rho(:,:,:) 3D_NOT_real
// - Variable inout sized array: efield(:,:,:,:) 4D_NOT_real
// - Variable inout sized array: phi(:,:,:) 3D_NOT_real

// Skipped unusable routine osc_get_cgrn_freespace:
// - Variable inout sized array: cgrn(:,:,:) 3D_NOT_complex

// Skipped unusable routine osc_getgrnfree:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageconvcorr:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageshift:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_osc_getgrnpipe(
    double& gam /* 0D_NOT_real inout */,
    double& a /* 0D_NOT_real inout */,
    double& b /* 0D_NOT_real inout */,
    double* delta /* 1D_NOT_real inout */,
    double* umin /* 1D_NOT_real inout */,
    int* npad /* 1D_NOT_integer inout */);
void osc_getgrnpipe(
    double& gam,
    double& a,
    double& b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_read_rectpipe_grn();
void osc_read_rectpipe_grn();

// Skipped unusable routine osc_rectpipe_solver:
// - Variable in sized array: rho(:,:,:) 3D_NOT_real
// - Variable out sized array: phi(:,:,:) 3D_NOT_real
// - Variable out sized array: efield(:,:,:,:) 4D_NOT_real
// - Variable out sized array: bfield(:,:,:,:) 4D_NOT_real
extern "C" void fortran_osc_write_rectpipe_grn(
    double& apipe /* 0D_NOT_real inout */,
    double& bpipe /* 0D_NOT_real inout */,
    double* delta /* 1D_NOT_real inout */,
    double* umin /* 1D_NOT_real inout */,
    double* umax /* 1D_NOT_real inout */,
    int* nlo /* 1D_NOT_integer inout */,
    int* nhi /* 1D_NOT_integer inout */,
    double& gamma /* 0D_NOT_real inout */);
void osc_write_rectpipe_grn(
    double& apipe,
    double& bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double& gamma);
extern "C" void fortran_parse_cartesian_map(
    void* ct_map /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parse_cartesian_map(
    CartesianMapProxy& ct_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);
extern "C" void fortran_parse_cylindrical_map(
    void* cl_map /* 0D_PTR_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parse_cylindrical_map(
    CylindricalMapProxy& cl_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);
extern "C" void fortran_parse_gen_grad_map(
    void* gg_map /* 0D_PTR_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parse_gen_grad_map(
    GenGradMapProxy& gg_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);
extern "C" void fortran_parse_grid_field(
    void* g_field /* 0D_PTR_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parse_grid_field(
    GridFieldProxy& g_field,
    EleProxy& ele,
    LatProxy& lat,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);
extern "C" bool fortran_parse_integer_list(
    const char* err_str /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    void* int_array /* 1D_ALLOC_integer inout */,
    bool& exact_size /* 0D_NOT_logical inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    const char* open_delim /* 0D_NOT_character inout */,
    const char* separator /* 0D_NOT_character inout */,
    const char* close_delim /* 0D_NOT_character inout */,
    int* default_value /* 0D_NOT_integer inout */,
    bool& is_ok /* 0D_NOT_logical inout */);
void parse_integer_list(
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
    bool& is_ok);
extern "C" bool fortran_parse_integer_list2(
    const char* err_str /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    void* int_array /* 1D_ALLOC_integer inout */,
    int& num_found /* 0D_NOT_integer out */,
    const char* delim /* 0D_NOT_character out */,
    bool& delim_found /* 0D_NOT_logical out */,
    int* num_expected /* 0D_NOT_integer inout */,
    const char* open_delim /* 0D_NOT_character inout */,
    const char* separator /* 0D_NOT_character inout */,
    const char* close_delim /* 0D_NOT_character inout */,
    int* default_value /* 0D_NOT_integer inout */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ParseIntegerList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
Bmad::ParseIntegerList2 parse_integer_list2(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    optional_ref<int> num_expected = std::nullopt,
    optional_ref<std::string> open_delim = std::nullopt,
    optional_ref<std::string> separator = std::nullopt,
    optional_ref<std::string> close_delim = std::nullopt,
    optional_ref<int> default_value = std::nullopt);

// Skipped unusable routine parse_line_or_list:
// - Untranslated type: seq_struct (1D)
extern "C" bool fortran_parse_real_list(
    void* lat /* 0D_NOT_type in */,
    const char* err_str /* 0D_NOT_character in */,
    void* real_array /* 1D_ALLOC_real out */,
    bool& exact_size /* 0D_NOT_logical in */,
    const char* delim /* 0D_NOT_character out */,
    bool& delim_found /* 0D_NOT_logical out */,
    const char* open_delim /* 0D_NOT_character in */,
    const char* separator /* 0D_NOT_character in */,
    const char* close_delim /* 0D_NOT_character in */,
    double* default_value /* 0D_NOT_real in */,
    int& num_found /* 0D_NOT_integer out */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ParseRealList {
  RealAlloc1D real_array;
  std::string delim;
  bool delim_found;
  int num_found;
  bool is_ok;
};
Bmad::ParseRealList parse_real_list(
    LatProxy& lat,
    std::string err_str,
    bool exact_size,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<double> default_value = std::nullopt);
extern "C" bool fortran_parse_real_list2(
    void* lat /* 0D_NOT_type in */,
    const char* err_str /* 0D_NOT_character in */,
    void* real_array /* 1D_ALLOC_real inout */,
    int& num_found /* 0D_NOT_integer out */,
    const char* delim /* 0D_NOT_character out */,
    bool& delim_found /* 0D_NOT_logical out */,
    int* num_expected /* 0D_NOT_integer inout */,
    const char* open_brace /* 0D_NOT_character inout */,
    const char* separator /* 0D_NOT_character inout */,
    const char* close_brace /* 0D_NOT_character inout */,
    double* default_value /* 0D_NOT_real inout */,
    bool* single_value /* 0D_NOT_logical inout */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ParseRealList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
Bmad::ParseRealList2 parse_real_list2(
    LatProxy& lat,
    std::string err_str,
    RealAlloc1D& real_array,
    optional_ref<int> num_expected = std::nullopt,
    optional_ref<std::string> open_brace = std::nullopt,
    optional_ref<std::string> separator = std::nullopt,
    optional_ref<std::string> close_brace = std::nullopt,
    optional_ref<double> default_value = std::nullopt,
    optional_ref<bool> single_value = std::nullopt);

// Skipped unusable routine parse_real_matrix:
// - Variable in sized array: table(:,:) 2D_ALLOC_real

// Skipped unusable routine parse_superimpose_command:
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser2_add_superimpose:
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser_add_branch:
// - Untranslated type: seq_struct (1D)
// - Variable-sized inout character array: seq_name(:) 1D_ALLOC_character
// - Untranslated type: parser_lat_struct (0D)
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_parser_add_constant(
    const char* word /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    bool& redef_is_error /* 0D_NOT_logical inout */);
void parser_add_constant(
    std::string& word,
    LatProxy& lat,
    bool& redef_is_error);

// Skipped unusable routine parser_add_lords:
// - Untranslated type: parser_lat_struct (0D)

// Skipped unusable routine parser_add_superimpose:
// - Untranslated type: parser_ele_struct (0D)
// - Untranslated type: parser_lat_struct (0D)
extern "C" void fortran_parser_call_check(
    const char* word /* 0D_NOT_character inout */,
    int& ix_word /* 0D_NOT_integer inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& call_found /* 0D_NOT_logical inout */,
    bool* err_flag /* 0D_NOT_logical inout */);
void parser_call_check(
    std::string& word,
    int& ix_word,
    std::string& delim,
    bool& delim_found,
    bool& call_found,
    optional_ref<bool> err_flag = std::nullopt);

// Skipped unusable routine parser_debug_print_info:
// - Untranslated type: seq_struct (1D)

// Skipped unusable routine parser_error:
// - Untranslated type: seq_struct (0D)
// - Untranslated type: parser_ele_struct (0D)

// Skipped unusable routine parser_expand_line:
// - Untranslated type: seq_struct (1D)
// - Variable-sized in character array: seq_name(:) 1D_ALLOC_character
// - Untranslated type: base_line_ele_struct (1D)
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_parser_fast_complex_read(
    void* cmplx_vec /* 1D_ALLOC_complex out */,
    void* ele /* 0D_NOT_type in */,
    const char* delim /* 0D_NOT_character out */,
    const char* err_str /* 0D_NOT_character in */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ParserFastComplexRead {
  ComplexAlloc1D cmplx_vec;
  std::string delim;
  bool is_ok;
};
Bmad::ParserFastComplexRead parser_fast_complex_read(
    EleProxy& ele,
    std::string err_str);
extern "C" bool fortran_parser_fast_integer_read(
    void* int_vec /* 1D_ALLOC_integer inout */,
    void* ele /* 0D_NOT_type inout */,
    const char* delim_wanted /* 0D_NOT_character inout */,
    const char* err_str /* 0D_NOT_character inout */,
    bool& is_ok /* 0D_NOT_logical inout */);
void parser_fast_integer_read(
    IntAlloc1D& int_vec,
    EleProxy& ele,
    std::string& delim_wanted,
    std::string& err_str,
    bool& is_ok);
extern "C" bool fortran_parser_fast_real_read(
    void* real_vec /* 1D_ALLOC_real out */,
    void* ele /* 0D_NOT_type in */,
    const char* end_delims /* 0D_NOT_character in */,
    const char* delim /* 0D_NOT_character out */,
    const char* err_str /* 0D_NOT_character in */,
    bool* exact_size /* 0D_NOT_logical in */,
    int& n_real /* 0D_NOT_integer out */,
    bool& is_ok /* 0D_NOT_logical out */);
struct ParserFastRealRead {
  RealAlloc1D real_vec;
  std::string delim;
  int n_real;
  bool is_ok;
};
Bmad::ParserFastRealRead parser_fast_real_read(
    EleProxy& ele,
    std::string end_delims,
    std::string err_str,
    std::optional<bool> exact_size = std::nullopt);
extern "C" void fortran_parser_file_stack(
    const char* how /* 0D_NOT_character inout */,
    const char* file_name_in /* 0D_NOT_character inout */,
    bool* finished /* 0D_NOT_logical inout */,
    bool* err /* 0D_NOT_logical inout */,
    bool* open_file /* 0D_NOT_logical inout */,
    bool* abort_on_open_error /* 0D_NOT_logical inout */);
void parser_file_stack(
    std::string& how,
    optional_ref<std::string> file_name_in = std::nullopt,
    optional_ref<bool> finished = std::nullopt,
    optional_ref<bool> err = std::nullopt,
    optional_ref<bool> open_file = std::nullopt,
    optional_ref<bool> abort_on_open_error = std::nullopt);
extern "C" void fortran_parser_get_integer(
    int& int_val /* 0D_NOT_integer inout */,
    const char* word /* 0D_NOT_character inout */,
    int& ix_word /* 0D_NOT_integer inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err /* 0D_NOT_logical inout */,
    const char* str1 /* 0D_NOT_character inout */,
    const char* str2 /* 0D_NOT_character inout */);
void parser_get_integer(
    int& int_val,
    std::string& word,
    int& ix_word,
    std::string& delim,
    bool& delim_found,
    bool& err,
    optional_ref<std::string> str1 = std::nullopt,
    optional_ref<std::string> str2 = std::nullopt);
extern "C" void fortran_parser_get_logical(
    const char* attrib_name /* 0D_NOT_character inout */,
    bool& this_logic /* 0D_NOT_logical inout */,
    const char* ele_name /* 0D_NOT_character inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err /* 0D_NOT_logical inout */);
void parser_get_logical(
    std::string& attrib_name,
    bool& this_logic,
    std::string& ele_name,
    std::string& delim,
    bool& delim_found,
    bool& err);
extern "C" void fortran_parser_identify_fork_to_element(
    void* lat /* 0D_NOT_type inout */);
void parser_identify_fork_to_element(LatProxy& lat);
extern "C" void fortran_parser_init_custom_elements(
    void* lat /* 0D_NOT_type inout */);
void parser_init_custom_elements(LatProxy& lat);
extern "C" void fortran_parser_print_line(
    void* lat /* 0D_NOT_type inout */,
    bool& end_of_file /* 0D_NOT_logical inout */);
void parser_print_line(LatProxy& lat, bool& end_of_file);
extern "C" void fortran_parser_read_lr_wake(
    void* ele /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parser_read_lr_wake(
    EleProxy& ele,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);
extern "C" void fortran_parser_read_old_format_lr_wake(
    void* ele /* 0D_NOT_type inout */,
    const char* lr_file_name /* 0D_NOT_character in */);
void parser_read_old_format_lr_wake(EleProxy& ele, std::string lr_file_name);
extern "C" void fortran_parser_read_old_format_sr_wake(
    void* ele /* 0D_NOT_type inout */,
    const char* sr_file_name /* 0D_NOT_character in */);
void parser_read_old_format_sr_wake(EleProxy& ele, std::string sr_file_name);
extern "C" void fortran_parser_read_sr_wake(
    void* ele /* 0D_NOT_type inout */,
    const char* delim /* 0D_NOT_character inout */,
    bool& delim_found /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void parser_read_sr_wake(
    EleProxy& ele,
    std::string& delim,
    bool& delim_found,
    bool& err_flag);

// Skipped unusable routine parser_set_attribute:
// - Untranslated type: parser_ele_struct (0D)
extern "C" void fortran_parser_transfer_control_struct(
    void* con_in /* 0D_NOT_type in */,
    void* con_out /* 0D_NOT_type out */,
    void* lord /* 0D_NOT_type in */,
    int& ix_var /* 0D_NOT_integer in */);
ControlProxy parser_transfer_control_struct(
    ControlProxy& con_in,
    EleProxy& lord,
    int ix_var);
extern "C" bool fortran_particle_in_global_frame(
    void* orb /* 0D_NOT_type in */,
    void* branch /* 0D_NOT_type in */,
    bool* in_time_coordinates /* 0D_NOT_logical in */,
    bool* in_body_frame /* 0D_NOT_logical in */,
    double* w_mat_out /* 2D_NOT_real inout */,
    void* particle /* 0D_NOT_type in */);
void particle_in_global_frame(
    CoordProxy& orb,
    BranchProxy& branch,
    std::optional<bool> in_time_coordinates,
    std::optional<bool> in_body_frame,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat_out,
    CoordProxy& particle);
extern "C" bool fortran_particle_is_moving_backwards(
    void* orbit /* 0D_NOT_type in */,
    bool& is_moving_backwards /* 0D_NOT_logical inout */);
void particle_is_moving_backwards(CoordProxy& orbit, bool& is_moving_backwards);
extern "C" bool fortran_particle_is_moving_forward(
    void* orbit /* 0D_NOT_type in */,
    int* dir /* 0D_NOT_integer in */,
    bool& is_moving_forward /* 0D_NOT_logical inout */);
void particle_is_moving_forward(
    CoordProxy& orbit,
    std::optional<int> dir,
    bool& is_moving_forward);
extern "C" bool fortran_particle_rf_time(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    bool* reference_active_edge /* 0D_NOT_logical in */,
    double* s_rel /* 0D_NOT_real in */,
    bool* time_coords /* 0D_NOT_logical in */,
    double* rf_freq /* 0D_NOT_real in */,
    bool* abs_time /* 0D_NOT_logical in */,
    long double& time /* 0D_NOT_real16 inout */);
void particle_rf_time(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<bool> reference_active_edge,
    std::optional<double> s_rel,
    std::optional<bool> time_coords,
    std::optional<double> rf_freq,
    std::optional<bool> abs_time,
    long double& time);
extern "C" bool fortran_patch_flips_propagation_direction(
    double& x_pitch /* 0D_NOT_real in */,
    double& y_pitch /* 0D_NOT_real in */,
    bool& is_flip /* 0D_NOT_logical inout */);
void patch_flips_propagation_direction(
    double x_pitch,
    double y_pitch,
    bool& is_flip);
extern "C" bool fortran_patch_length(
    void* patch /* 0D_NOT_type in */,
    int* ref_coords /* 0D_NOT_integer in */,
    double& length /* 0D_NOT_real inout */);
void patch_length(
    EleProxy& patch,
    std::optional<int> ref_coords,
    double& length);
extern "C" void fortran_photon_absorption_and_phase_shift(
    const char* material /* 0D_NOT_character in */,
    double& Energy /* 0D_NOT_real in */,
    double& absorption /* 0D_NOT_real out */,
    double& phase_shift /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct PhotonAbsorptionAndPhaseShift {
  double absorption;
  double phase_shift;
  bool err_flag;
};
Bmad::PhotonAbsorptionAndPhaseShift photon_absorption_and_phase_shift(
    std::string material,
    double Energy);
extern "C" void fortran_photon_add_to_detector_statistics(
    void* orbit0 /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type inout */,
    int* ix_pt /* 0D_NOT_integer inout */,
    int* iy_pt /* 0D_NOT_integer inout */,
    void* pixel_pt /* 0D_NOT_type in */);
void photon_add_to_detector_statistics(
    CoordProxy& orbit0,
    CoordProxy& orbit,
    EleProxy& ele,
    optional_ref<int> ix_pt = std::nullopt,
    optional_ref<int> iy_pt = std::nullopt,
    optional_ref<PixelPtProxy> pixel_pt = std::nullopt);

// Skipped unusable routine photon_diffuse_scattering:
// - Untranslated type: diffuse_param_struct (0D)

// Skipped unusable routine photon_read_spline:
// - Untranslated type: photon_init_splines_struct (0D)
extern "C" void fortran_photon_reflection(
    double& graze_angle_in /* 0D_NOT_real in */,
    double& energy /* 0D_NOT_real in */,
    void* surface /* 0D_NOT_type in */,
    double& graze_angle_out /* 0D_NOT_real out */,
    double& phi_out /* 0D_NOT_real out */);
struct PhotonReflection {
  double graze_angle_out;
  double phi_out;
};
Bmad::PhotonReflection photon_reflection(
    double graze_angle_in,
    double energy,
    PhotonReflectSurfaceProxy& surface);
extern "C" void fortran_photon_reflection_std_surface_init(
    void* surface /* 0D_NOT_type out */);
PhotonReflectSurfaceProxy photon_reflection_std_surface_init();
extern "C" void fortran_photon_reflectivity(
    double& angle /* 0D_NOT_real in */,
    double& energy /* 0D_NOT_real in */,
    void* surface /* 0D_NOT_type in */,
    double& p_reflect /* 0D_NOT_real out */,
    double& rel_p_specular /* 0D_NOT_real out */);
struct PhotonReflectivity {
  double p_reflect;
  double rel_p_specular;
};
Bmad::PhotonReflectivity photon_reflectivity(
    double angle,
    double energy,
    PhotonReflectSurfaceProxy& surface);
extern "C" void fortran_photon_target_corner_calc(
    void* aperture_ele /* 0D_NOT_type in */,
    double& x_lim /* 0D_NOT_real inout */,
    double& y_lim /* 0D_NOT_real inout */,
    double& z_lim /* 0D_NOT_real inout */,
    void* source_ele /* 0D_NOT_type in */,
    void* corner /* 0D_NOT_type out */);
TargetPointProxy photon_target_corner_calc(
    EleProxy& aperture_ele,
    double& x_lim,
    double& y_lim,
    double& z_lim,
    EleProxy& source_ele);
extern "C" void fortran_photon_target_setup(void* ele /* 0D_NOT_type inout */);
void photon_target_setup(EleProxy& ele);
extern "C" bool fortran_photon_type(
    void* ele /* 0D_NOT_type in */,
    int& e_type /* 0D_NOT_integer out */);
int photon_type(EleProxy& ele);
extern "C" bool fortran_physical_ele_end(
    int& track_end /* 0D_NOT_integer in */,
    void* orbit /* 0D_NOT_type in */,
    int& ele_orientation /* 0D_NOT_integer in */,
    bool* return_stream_end /* 0D_NOT_logical in */,
    int& physical_end /* 0D_NOT_integer inout */);
void physical_ele_end(
    int track_end,
    CoordProxy& orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end,
    int& physical_end);
extern "C" void fortran_point_photon_emission(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    int& direction /* 0D_NOT_integer in */,
    double& max_target_area /* 0D_NOT_real in */,
    double* w_to_surface /* 2D_NOT_real in */);
void point_photon_emission(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    int direction,
    double max_target_area,
    std::optional<FixedArray2D<Real, 3, 3>> w_to_surface = std::nullopt);

// Skipped unusable routine pointer_to_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_pointer_to_branch_given_ele(
    void* ele /* 0D_NOT_type in */,
    void* branch_ptr /* 0D_PTR_type out */);
BranchProxy pointer_to_branch(EleProxy& ele);
extern "C" bool fortran_pointer_to_branch_given_name(
    const char* branch_name /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    bool* parameter_is_branch0 /* 0D_NOT_logical in */,
    int* blank_branch /* 0D_NOT_integer in */,
    void* branch_ptr /* 0D_PTR_type out */);
BranchProxy pointer_to_branch(
    std::string branch_name,
    LatProxy& lat,
    std::optional<bool> parameter_is_branch0 = std::nullopt,
    std::optional<int> blank_branch = std::nullopt);
extern "C" bool fortran_pointer_to_ele1(
    void* lat /* 0D_NOT_type inout */,
    int& ix_ele /* 0D_NOT_integer inout */,
    int* ix_branch /* 0D_NOT_integer inout */,
    void* ele_ptr /* 0D_PTR_type inout */);
void pointer_to_ele1(
    LatProxy& lat,
    int& ix_ele,
    optional_ref<int> ix_branch,
    EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele2(
    void* lat /* 0D_NOT_type inout */,
    void* ele_loc /* 0D_NOT_type inout */,
    void* ele_ptr /* 0D_PTR_type inout */);
void pointer_to_ele2(LatProxy& lat, LatEleLocProxy& ele_loc, EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele3(
    void* lat /* 0D_NOT_type inout */,
    const char* ele_name /* 0D_NOT_character inout */,
    void* ele_ptr /* 0D_PTR_type inout */);
void pointer_to_ele3(LatProxy& lat, std::string& ele_name, EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele4(
    void* lat /* 0D_NOT_type inout */,
    void* foreign_ele /* 0D_NOT_type inout */,
    void* ele_ptr /* 0D_PTR_type inout */);
void pointer_to_ele4(LatProxy& lat, EleProxy& foreign_ele, EleProxy& ele_ptr);

// Skipped unusable routine pointer_to_ele_multipole:
// - Routine in configuration skip list
extern "C" bool fortran_pointer_to_element_at_s(
    void* branch /* 0D_NOT_type in */,
    double& s /* 0D_NOT_real in */,
    bool& choose_max /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& s_eff /* 0D_NOT_real out */,
    void* position /* 0D_NOT_type out */,
    bool* print_err /* 0D_NOT_logical in */,
    void* ele /* 0D_PTR_type out */);
struct PointerToElementAtS {
  bool err_flag;
  double s_eff;
  CoordProxy position;
  EleProxy ele;
};
Bmad::PointerToElementAtS pointer_to_element_at_s(
    BranchProxy& branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine pointer_to_fibre:
// - Untranslated type: fibre (0D)
extern "C" bool fortran_pointer_to_field_ele(
    void* ele /* 0D_NOT_type in */,
    int& ix_field_ele /* 0D_NOT_integer in */,
    double& dz_offset /* 0D_NOT_real out */,
    void* field_ele /* 0D_PTR_type inout */);
double pointer_to_field_ele(
    EleProxy& ele,
    int ix_field_ele,
    EleProxy& field_ele);
extern "C" bool fortran_pointer_to_girder(
    void* ele /* 0D_NOT_type in */,
    int& ix_slave_back /* 0D_NOT_integer out */,
    void* girder /* 0D_PTR_type inout */);
int pointer_to_girder(EleProxy& ele, EleProxy& girder);

// Skipped unusable routine pointer_to_indexed_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_pointer_to_lord(
    void* slave /* 0D_NOT_type in */,
    int& ix_lord /* 0D_NOT_integer in */,
    void* control /* 0D_PTR_type out */,
    int& ix_slave_back /* 0D_NOT_integer out */,
    int* lord_type /* 0D_NOT_integer in */,
    int& ix_control /* 0D_NOT_integer out */,
    int& ix_ic /* 0D_NOT_integer out */,
    void* lord_ptr /* 0D_PTR_type inout */);
struct PointerToLord {
  ControlProxy control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
};
Bmad::PointerToLord pointer_to_lord(
    EleProxy& slave,
    int ix_lord,
    std::optional<int> lord_type,
    EleProxy& lord_ptr);
extern "C" bool fortran_pointer_to_multipass_lord(
    void* ele /* 0D_NOT_type in */,
    int& ix_pass /* 0D_NOT_integer out */,
    void* super_lord /* 0D_PTR_type out */,
    void* multi_lord /* 0D_PTR_type inout */);
struct PointerToMultipassLord {
  int ix_pass;
  EleProxy super_lord;
};
Bmad::PointerToMultipassLord pointer_to_multipass_lord(
    EleProxy& ele,
    EleProxy& multi_lord);
extern "C" bool fortran_pointer_to_next_ele(
    void* this_ele /* 0D_NOT_type in */,
    int* offset /* 0D_NOT_integer in */,
    bool* skip_beginning /* 0D_NOT_logical in */,
    bool* follow_fork /* 0D_NOT_logical in */,
    void* next_ele /* 0D_PTR_type inout */);
void pointer_to_next_ele(
    EleProxy& this_ele,
    std::optional<int> offset,
    std::optional<bool> skip_beginning,
    std::optional<bool> follow_fork,
    EleProxy& next_ele);
extern "C" bool fortran_pointer_to_slave(
    void* lord /* 0D_NOT_type in */,
    int& ix_slave /* 0D_NOT_integer in */,
    void* control /* 0D_PTR_type out */,
    int* lord_type /* 0D_NOT_integer in */,
    int& ix_lord_back /* 0D_NOT_integer out */,
    int& ix_control /* 0D_NOT_integer out */,
    int& ix_ic /* 0D_NOT_integer out */,
    void* slave_ptr /* 0D_PTR_type out */);
struct PointerToSlave {
  ControlProxy control;
  int ix_lord_back;
  int ix_control;
  int ix_ic;
  EleProxy slave_ptr;
};
Bmad::PointerToSlave pointer_to_slave(
    EleProxy& lord,
    int ix_slave,
    std::optional<int> lord_type = std::nullopt);
extern "C" bool fortran_pointer_to_super_lord(
    void* slave /* 0D_NOT_type in */,
    void* control /* 0D_PTR_type out */,
    int& ix_slave_back /* 0D_NOT_integer out */,
    int& ix_control /* 0D_NOT_integer out */,
    int& ix_ic /* 0D_NOT_integer out */,
    int* lord_type /* 0D_NOT_integer in */,
    void* lord_ptr /* 0D_PTR_type inout */);
struct PointerToSuperLord {
  ControlProxy control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
};
Bmad::PointerToSuperLord pointer_to_super_lord(
    EleProxy& slave,
    std::optional<int> lord_type,
    EleProxy& lord_ptr);
extern "C" bool fortran_pointer_to_surface_displacement_pt(
    void* ele /* 0D_NOT_type in */,
    bool& nearest /* 0D_NOT_logical in */,
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    int* ix /* 0D_NOT_integer inout */,
    int* iy /* 0D_NOT_integer inout */,
    bool* extend_grid /* 0D_NOT_logical in */,
    double* xx /* 0D_NOT_real inout */,
    double* yy /* 0D_NOT_real inout */,
    void* pt /* 0D_PTR_type out */);
SurfaceDisplacementPtProxy pointer_to_surface_displacement_pt(
    EleProxy& ele,
    bool nearest,
    double& x,
    double& y,
    optional_ref<int> ix = std::nullopt,
    optional_ref<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    optional_ref<double> xx = std::nullopt,
    optional_ref<double> yy = std::nullopt);
extern "C" bool fortran_pointer_to_surface_segmented_pt(
    void* ele /* 0D_NOT_type in */,
    bool& nearest /* 0D_NOT_logical in */,
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    int* ix /* 0D_NOT_integer inout */,
    int* iy /* 0D_NOT_integer inout */,
    bool* extend_grid /* 0D_NOT_logical in */,
    double* xx /* 0D_NOT_real inout */,
    double* yy /* 0D_NOT_real inout */,
    void* pt /* 0D_PTR_type out */);
SurfaceSegmentedPtProxy pointer_to_surface_segmented_pt(
    EleProxy& ele,
    bool nearest,
    double& x,
    double& y,
    optional_ref<int> ix = std::nullopt,
    optional_ref<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    optional_ref<double> xx = std::nullopt,
    optional_ref<double> yy = std::nullopt);
extern "C" bool fortran_pointer_to_wake_ele(
    void* ele /* 0D_NOT_type in */,
    double& delta_s /* 0D_NOT_real out */,
    void* wake_ele /* 0D_PTR_type inout */);
double pointer_to_wake_ele(EleProxy& ele, EleProxy& wake_ele);
extern "C" bool fortran_pointer_to_wall3d(
    void* ele /* 0D_NOT_type in */,
    int* ix_wall /* 0D_NOT_integer in */,
    double& ds_offset /* 0D_NOT_real out */,
    bool& is_branch_wall /* 0D_NOT_logical out */,
    void* wall3d /* 0D_PTR_type out */);
struct PointerToWall3d {
  double ds_offset;
  bool is_branch_wall;
  Wall3dProxy wall3d;
};
Bmad::PointerToWall3d pointer_to_wall3d(
    EleProxy& ele,
    std::optional<int> ix_wall = std::nullopt);

// Skipped unusable routine pointers_to_attribute:
// - Untranslated type: all_pointer_struct (1D)
extern "C" bool fortran_polar_to_spinor(
    void* polar /* 0D_NOT_type in */,
    std::complex<double>* spinor /* 1D_NOT_complex inout */);
void polar_to_spinor(SpinPolarProxy& polar, FixedArray1D<Complex, 2> spinor);
extern "C" bool fortran_polar_to_vec(
    void* polar /* 0D_NOT_type in */,
    double* vec /* 1D_NOT_real inout */);
void polar_to_vec(SpinPolarProxy& polar, FixedArray1D<Real, 3> vec);

// Skipped unusable routine print_mesh3d:
// - Untranslated type: mesh3d_struct (0D)

// Skipped unusable routine prob_x_diffuse:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" void fortran_project_emit_to_xyz(
    void* ring /* 0D_NOT_type in */,
    int& ix /* 0D_NOT_integer in */,
    void* mode /* 0D_NOT_type in */,
    double& sigma_x /* 0D_NOT_real out */,
    double& sigma_y /* 0D_NOT_real out */,
    double& sigma_z /* 0D_NOT_real out */);
struct ProjectEmitToXyz {
  double sigma_x;
  double sigma_y;
  double sigma_z;
};
Bmad::ProjectEmitToXyz project_emit_to_xyz(
    LatProxy& ring,
    int ix,
    NormalModesProxy& mode);

// Skipped unusable routine propagate_part_way:
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)

// Skipped unusable routine psi_prime:
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_psi_prime_sca(
    double& t /* 0D_NOT_real in */,
    double& p /* 0D_NOT_real in */,
    double& dpdt /* 0D_NOT_real out */,
    double* args /* 1D_NOT_real in */);
double psi_prime_sca(double t, double p, FixedArray1D<Real, 8> args);
extern "C" void fortran_ptc_bookkeeper(void* lat /* 0D_NOT_type inout */);
void ptc_bookkeeper(LatProxy& lat);

// Skipped unusable routine ptc_calculate_tracking_step_size:
// - Untranslated type: layout (0D)

// Skipped unusable routine ptc_check_for_lost_particle:
// - Untranslated type: fibre (0D)
extern "C" void fortran_ptc_closed_orbit_calc(
    void* branch /* 0D_NOT_type in */,
    void* closed_orbit /* 1D_ALLOC_type out */,
    bool* radiation_damping_on /* 0D_NOT_logical in */);
CoordProxyAlloc1D ptc_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<bool> radiation_damping_on = std::nullopt);
extern "C" void fortran_ptc_emit_calc(
    void* ele /* 0D_NOT_type in */,
    void* norm_mode /* 0D_NOT_type out */,
    double* sigma_mat /* 2D_NOT_real inout */,
    void* closed_orb /* 0D_NOT_type out */);
struct PtcEmitCalc {
  NormalModesProxy norm_mode;
  CoordProxy closed_orb;
};
Bmad::PtcEmitCalc ptc_emit_calc(
    EleProxy& ele,
    FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_kill_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_layouts_resplit(
    double& dKL_max /* 0D_NOT_real in */,
    double& l_max /* 0D_NOT_real in */,
    bool& l_max_drift_only /* 0D_NOT_logical in */,
    double& bend_dorb /* 0D_NOT_real in */,
    double& sex_dx /* 0D_NOT_real in */,
    bool* even /* 0D_NOT_logical in */,
    int* crossover /* 1D_NOT_integer in */,
    int* crossover_wiggler /* 1D_NOT_integer inout */);
void ptc_layouts_resplit(
    double dKL_max,
    double l_max,
    bool l_max_drift_only,
    double bend_dorb,
    double sex_dx,
    std::optional<bool> even = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover = std::nullopt,
    std::optional<FixedArray1D<Int, 2>> crossover_wiggler = std::nullopt);

// Skipped unusable routine ptc_linear_isf_calc:
// - Untranslated type: linear_ele_isf_struct (1D)

// Skipped unusable routine ptc_map_to_normal_form:
// - Untranslated type: probe_8 (0D)
// - Untranslated type: c_normal_form (0D)
// - Untranslated type: c_taylor (1D)
// - Untranslated type: c_taylor (0D)

// Skipped unusable routine ptc_one_turn_map_at_ele:
// - Untranslated type: probe_8 (0D)
// - Untranslated type: internal_state (0D)
extern "C" void fortran_ptc_one_turn_mat_and_closed_orbit_calc(
    void* branch /* 0D_NOT_type inout */,
    double* pz /* 0D_NOT_real in */);
void ptc_one_turn_mat_and_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<double> pz = std::nullopt);
extern "C" void fortran_ptc_ran_seed_put(int& iseed /* 0D_NOT_integer in */);
void ptc_ran_seed_put(int iseed);

// Skipped unusable routine ptc_read_flat_file:
// - Variable-sized in character array: flat_file(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine ptc_read_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_set_rf_state_for_c_normal(
    bool& nocavity /* 0D_NOT_logical in */);
void ptc_set_rf_state_for_c_normal(bool nocavity);
extern "C" void fortran_ptc_set_taylor_order_if_needed();
void ptc_set_taylor_order_if_needed();

// Skipped unusable routine ptc_setup_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)

// Skipped unusable routine ptc_setup_tracking_with_damping_and_excitation:
// - Untranslated type: internal_state (0D)
extern "C" void fortran_ptc_spin_calc(
    void* ele /* 0D_NOT_type in */,
    void* norm_mode /* 0D_NOT_type out */,
    double* sigma_mat /* 2D_NOT_real inout */,
    void* closed_orb /* 0D_NOT_type out */);
struct PtcSpinCalc {
  NormalModesProxy norm_mode;
  CoordProxy closed_orb;
};
Bmad::PtcSpinCalc ptc_spin_calc(
    EleProxy& ele,
    FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_spin_matching_calc:
// - Untranslated type: spin_matching_struct (1D)

// Skipped unusable routine ptc_taylors_equal_bmad_taylors:
// - Untranslated type: taylor (1D)
extern "C" void fortran_ptc_track_all(
    void* branch /* 0D_NOT_type in */,
    void* orbit /* 1D_ALLOC_type inout */,
    int& track_state /* 0D_NOT_integer out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct PtcTrackAll {
  int track_state;
  bool err_flag;
};
Bmad::PtcTrackAll ptc_track_all(BranchProxy& branch, CoordProxyAlloc1D& orbit);

// Skipped unusable routine ptc_track_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)
extern "C" void fortran_ptc_transfer_map_with_spin(
    void* branch /* 0D_NOT_type in */,
    void* t_map /* 1D_NOT_type inout */,
    void* s_map /* 1D_NOT_type inout */,
    void* orb0 /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    int* ix1 /* 0D_NOT_integer in */,
    int* ix2 /* 0D_NOT_integer in */,
    bool* one_turn /* 0D_NOT_logical in */,
    bool* unit_start /* 0D_NOT_logical in */);
bool ptc_transfer_map_with_spin(
    BranchProxy& branch,
    FixedArray1D<TaylorProxy, 6> t_map,
    FixedArray1D<TaylorProxy, 4> s_map,
    CoordProxy& orb0,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt);

// Skipped unusable routine ptc_write_map_with_radiation:
// - Untranslated type: ptc_rad_map_struct (0D)

// Skipped unusable routine ptwo:
// - Untranslated type: diffuse_param_struct (0D)
extern "C" bool fortran_pwd_mat(
    void* lat /* 0D_NOT_type in */,
    double* t6 /* 2D_NOT_real in */,
    double& inductance /* 0D_NOT_real in */,
    double& sig_z /* 0D_NOT_real in */,
    double* t6_pwd /* 2D_NOT_real out */);
FixedArray2D<Real, 6, 6> pwd_mat(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> t6,
    double inductance,
    double sig_z);

// Skipped unusable routine qromb_rad_int:
// - Untranslated type: rad_int_track_point_struct (0D)
// - Untranslated type: rad_int_info_struct (0D)

// Skipped unusable routine quad_mat2_calc:
// - Variable out sized array: mat2(:,:) 2D_NOT_real
extern "C" void fortran_rad1_damp_and_stoc_mats(
    void* ele /* 0D_NOT_type in */,
    bool& include_opening_angle /* 0D_NOT_logical in */,
    void* orb_in /* 0D_NOT_type in */,
    void* orb_out /* 0D_NOT_type in */,
    void* rad_map /* 0D_NOT_type out */,
    double& g2_tol /* 0D_NOT_real in */,
    double& g3_tol /* 0D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* ele0 /* 0D_NOT_type in */,
    void* rad_int1 /* 0D_NOT_type out */);
struct Rad1DampAndStocMats {
  RadMapProxy rad_map;
  bool err_flag;
  RadInt1Proxy rad_int1;
};
Bmad::Rad1DampAndStocMats rad1_damp_and_stoc_mats(
    EleProxy& ele,
    bool include_opening_angle,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double g2_tol,
    double g3_tol,
    optional_ref<EleProxy> ele0 = std::nullopt);
extern "C" void fortran_rad_damp_and_stoc_mats(
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    bool& include_opening_angle /* 0D_NOT_logical in */,
    void* rmap /* 0D_NOT_type out */,
    void* mode /* 0D_NOT_type out */,
    double* xfer_nodamp_mat /* 2D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* closed_orbit /* 1D_ALLOC_type in */,
    void* rad_int_branch /* 0D_NOT_type out */);
struct RadDampAndStocMats {
  RadMapProxy rmap;
  NormalModesProxy mode;
  FixedArray2D<Real, 6, 6> xfer_nodamp_mat;
  bool err_flag;
  RadIntBranchProxy rad_int_branch;
};
Bmad::RadDampAndStocMats rad_damp_and_stoc_mats(
    EleProxy& ele1,
    EleProxy& ele2,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit = std::nullopt);
extern "C" void fortran_rad_g_integrals(
    void* ele /* 0D_NOT_type in */,
    int& where /* 0D_NOT_integer in */,
    void* orb_in /* 0D_NOT_type in */,
    void* orb_out /* 0D_NOT_type in */,
    double* int_g /* 1D_NOT_real out */,
    double& int_g2 /* 0D_NOT_real inout */,
    double& int_g3 /* 0D_NOT_real inout */,
    double& g_tol /* 0D_NOT_real in */,
    double& g2_tol /* 0D_NOT_real in */,
    double& g3_tol /* 0D_NOT_real in */);
FixedArray1D<Real, 2> rad_g_integrals(
    EleProxy& ele,
    int where,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double& int_g2,
    double& int_g3,
    double g_tol,
    double g2_tol,
    double g3_tol);
extern "C" void fortran_radiation_integrals(
    void* lat /* 0D_NOT_type in */,
    void* orbit /* 1D_ALLOC_type in */,
    void* mode /* 0D_NOT_type out */,
    int* ix_cache /* 0D_NOT_integer inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    void* rad_int_by_ele /* 0D_NOT_type out */);
struct RadiationIntegrals {
  NormalModesProxy mode;
  RadIntAllEleProxy rad_int_by_ele;
};
Bmad::RadiationIntegrals radiation_integrals(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    optional_ref<int> ix_cache = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine radiation_integrals_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_radiation_map_setup(
    void* ele /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* ref_orbit_in /* 0D_NOT_type inout */);
bool radiation_map_setup(
    EleProxy& ele,
    optional_ref<CoordProxy> ref_orbit_in = std::nullopt);
extern "C" void fortran_ramper_slave_setup(
    void* lat /* 0D_NOT_type inout */,
    bool* force_setup /* 0D_NOT_logical in */);
void ramper_slave_setup(
    LatProxy& lat,
    std::optional<bool> force_setup = std::nullopt);
extern "C" bool fortran_ramper_value(
    void* ramper /* 0D_NOT_type in */,
    void* r1 /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& value /* 0D_NOT_real inout */);
bool ramper_value(EleProxy& ramper, ControlRamp1Proxy& r1, double& value);
extern "C" void fortran_randomize_lr_wake_frequencies(
    void* ele /* 0D_NOT_type inout */,
    bool& set_done /* 0D_NOT_logical out */);
bool randomize_lr_wake_frequencies(EleProxy& ele);
extern "C" bool fortran_rchomp(
    double& rel /* 0D_NOT_real inout */,
    int& plc /* 0D_NOT_integer inout */,
    const char* out /* 0D_NOT_character inout */);
void rchomp(double& rel, int& plc, std::string& out);

// Skipped unusable routine rclog_integrand:
// - Untranslated type: c_ptr (0D)
extern "C" void fortran_re_allocate_eles(
    void* eles /* 1D_ALLOC_type inout */,
    int& n /* 0D_NOT_integer in */,
    bool* save_old /* 0D_NOT_logical in */,
    bool* exact /* 0D_NOT_logical in */);
void re_allocate_eles(
    ElePointerProxyAlloc1D& eles,
    int n,
    std::optional<bool> save_old = std::nullopt,
    std::optional<bool> exact = std::nullopt);
extern "C" void fortran_re_allocate_wall3d_section_array(
    void* section /* 1D_ALLOC_type inout */,
    int& n /* 0D_NOT_integer in */,
    bool* exact /* 0D_NOT_logical in */);
void re_allocate_wall3d_section_array(
    Wall3dSectionProxyAlloc1D& section,
    int n,
    std::optional<bool> exact = std::nullopt);
extern "C" void fortran_re_allocate_wall3d_vertex_array(
    void* v /* 1D_ALLOC_type inout */,
    int& n /* 0D_NOT_integer in */,
    bool* exact /* 0D_NOT_logical in */);
void re_allocate_wall3d_vertex_array(
    Wall3dVertexProxyAlloc1D& v,
    int n,
    std::optional<bool> exact = std::nullopt);
extern "C" void fortran_re_associate_node_array(
    void* tree /* 0D_NOT_type inout */,
    int& n /* 0D_NOT_integer in */,
    bool* exact /* 0D_NOT_logical in */);
void re_associate_node_array(
    ExpressionTreeProxy& tree,
    int n,
    std::optional<bool> exact = std::nullopt);
extern "C" bool fortran_re_str_qp(
    long double& rel /* 0D_NOT_real16 inout */,
    const char* str_out /* 0D_NOT_character inout */);
void re_str_qp(long double& rel, std::string& str_out);
extern "C" bool fortran_re_str_rp(
    double& rel /* 0D_NOT_real inout */,
    const char* str_out /* 0D_NOT_character inout */);
void re_str_rp(double& rel, std::string& str_out);
extern "C" void fortran_read_beam_ascii(
    const char* file_name /* 0D_NOT_character in */,
    void* beam /* 0D_NOT_type out */,
    void* beam_init /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct ReadBeamAscii {
  BeamProxy beam;
  bool err_flag;
};
Bmad::ReadBeamAscii read_beam_ascii(
    std::string file_name,
    BeamInitProxy& beam_init);
extern "C" void fortran_read_beam_file(
    const char* file_name /* 0D_NOT_character in */,
    void* beam /* 0D_NOT_type out */,
    void* beam_init /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* ele /* 0D_NOT_type in */,
    bool* print_mom_shift_warning /* 0D_NOT_logical in */,
    bool* conserve_momentum /* 0D_NOT_logical inout */);
struct ReadBeamFile {
  BeamProxy beam;
  bool err_flag;
};
Bmad::ReadBeamFile read_beam_file(
    std::string file_name,
    BeamInitProxy& beam_init,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<bool> print_mom_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_read_binary_cartesian_map(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* cart_map /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void read_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag);
extern "C" void fortran_read_binary_cylindrical_map(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* cl_map /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void read_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag);
extern "C" void fortran_read_binary_grid_field(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* g_field /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void read_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag);

// Skipped unusable routine read_digested_bmad_file:
// - Variable-sized out character array: lat_files(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_read_surface_reflection_file(
    const char* file_name /* 0D_NOT_character in */,
    void* surface /* 0D_NOT_type out */);
PhotonReflectSurfaceProxy read_surface_reflection_file(std::string file_name);

// Skipped unusable routine real_8_to_taylor:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_reallocate_beam(
    void* beam /* 0D_NOT_type inout */,
    int& n_bunch /* 0D_NOT_integer in */,
    int* n_particle /* 0D_NOT_integer in */,
    bool* extend /* 0D_NOT_logical inout */);
void reallocate_beam(
    BeamProxy& beam,
    int n_bunch,
    std::optional<int> n_particle = std::nullopt,
    optional_ref<bool> extend = std::nullopt);
extern "C" void fortran_reallocate_bp_com_const();
void reallocate_bp_com_const();
extern "C" void fortran_reallocate_bunch(
    void* bunch /* 0D_NOT_type out */,
    int& n_particle /* 0D_NOT_integer in */,
    bool* save /* 0D_NOT_logical in */);
BunchProxy reallocate_bunch(
    int n_particle,
    std::optional<bool> save = std::nullopt);
extern "C" void fortran_reallocate_control(
    void* lat /* 0D_NOT_type inout */,
    int& n /* 0D_NOT_integer in */);
void reallocate_control(LatProxy& lat, int n);
extern "C" void fortran_reallocate_coord_array(
    void* coord_array /* 1D_ALLOC_type inout */,
    void* lat /* 0D_NOT_type in */);
void reallocate_coord(CoordArrayProxyAlloc1D& coord_array, LatProxy& lat);
extern "C" void fortran_reallocate_coord_lat(
    void* coord /* 1D_ALLOC_type inout */,
    void* lat /* 0D_NOT_type in */,
    int* ix_branch /* 0D_NOT_integer in */);
void reallocate_coord(
    CoordProxyAlloc1D& coord,
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_reallocate_coord_n(
    void* coord /* 1D_ALLOC_type inout */,
    int& n_coord /* 0D_NOT_integer in */);
void reallocate_coord(CoordProxyAlloc1D& coord, int n_coord);
extern "C" void fortran_reallocate_expression_stack(
    void* stack /* 1D_ALLOC_type inout */,
    int& n /* 0D_NOT_integer in */,
    bool* exact /* 0D_NOT_logical in */);
void reallocate_expression_stack(
    ExpressionAtomProxyAlloc1D& stack,
    int n,
    std::optional<bool> exact = std::nullopt);

// Skipped unusable routine reallocate_sequence:
// - Untranslated type: seq_struct (1D)

// Skipped unusable routine reals_8_equal_bmad_taylors:
// - Untranslated type: real_8 (1D)
extern "C" bool fortran_rel_tracking_charge_to_mass(
    void* orbit /* 0D_NOT_type in */,
    int& ref_species /* 0D_NOT_integer in */,
    double& rel_charge /* 0D_NOT_real inout */);
void rel_tracking_charge_to_mass(
    CoordProxy& orbit,
    int ref_species,
    double& rel_charge);
extern "C" bool fortran_relative_mode_flip(
    void* ele1 /* 0D_NOT_type inout */,
    void* ele2 /* 0D_NOT_type inout */,
    bool& func_retval__ /* 0D_NOT_logical inout */);
void relative_mode_flip(EleProxy& ele1, EleProxy& ele2, bool& func_retval__);
extern "C" void fortran_release_rad_int_cache(
    int& ix_cache /* 0D_NOT_integer inout */);
void release_rad_int_cache(int& ix_cache);
extern "C" void fortran_remove_constant_taylor(
    void* taylor_in /* 1D_ALLOC_type in */,
    void* taylor_out /* 1D_ALLOC_type out */,
    void* c0 /* 1D_ALLOC_real out */,
    bool& remove_higher_order_terms /* 0D_NOT_logical in */);
struct RemoveConstantTaylor {
  TaylorProxyAlloc1D taylor_out;
  RealAlloc1D c0;
};
Bmad::RemoveConstantTaylor remove_constant_taylor(
    TaylorProxyAlloc1D& taylor_in,
    bool remove_higher_order_terms);
extern "C" void fortran_remove_dead_from_bunch(
    void* bunch_in /* 0D_NOT_type in */,
    void* bunch_out /* 0D_NOT_type out */);
BunchProxy remove_dead_from_bunch(BunchProxy& bunch_in);
extern "C" void fortran_remove_eles_from_lat(
    void* lat /* 0D_NOT_type inout */,
    bool* check_sanity /* 0D_NOT_logical in */);
void remove_eles_from_lat(
    LatProxy& lat,
    std::optional<bool> check_sanity = std::nullopt);
extern "C" void fortran_remove_lord_slave_link(
    void* lord /* 0D_NOT_type inout */,
    void* slave /* 0D_NOT_type inout */);
void remove_lord_slave_link(EleProxy& lord, EleProxy& slave);
extern "C" void fortran_reverse_lat(
    void* lat_in /* 0D_NOT_type in */,
    void* lat_rev /* 0D_NOT_type out */,
    bool* track_antiparticle /* 0D_NOT_logical in */);
LatProxy reverse_lat(
    LatProxy& lat_in,
    std::optional<bool> track_antiparticle = std::nullopt);
extern "C" void fortran_rf_coupler_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    double& phase /* 0D_NOT_real in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void rf_coupler_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    double phase,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_rf_is_on(
    void* branch /* 0D_NOT_type in */,
    int* ix_ele1 /* 0D_NOT_integer in */,
    int* ix_ele2 /* 0D_NOT_integer in */,
    bool& is_on /* 0D_NOT_logical inout */);
void rf_is_on(
    BranchProxy& branch,
    std::optional<int> ix_ele1,
    std::optional<int> ix_ele2,
    bool& is_on);
extern "C" bool fortran_rf_ref_time_offset(
    void* ele /* 0D_NOT_type in */,
    double* ds /* 0D_NOT_real in */,
    double& time /* 0D_NOT_real inout */);
void rf_ref_time_offset(EleProxy& ele, std::optional<double> ds, double& time);
extern "C" bool fortran_rfun(
    double& u /* 0D_NOT_real inout */,
    double& v /* 0D_NOT_real inout */,
    double& w /* 0D_NOT_real inout */,
    double& gam /* 0D_NOT_real inout */,
    double& a /* 0D_NOT_real inout */,
    double& b /* 0D_NOT_real inout */,
    double& hz /* 0D_NOT_real inout */,
    int& i /* 0D_NOT_integer inout */,
    int& j /* 0D_NOT_integer inout */,
    double& res /* 0D_NOT_real inout */);
void rfun(
    double& u,
    double& v,
    double& w,
    double& gam,
    double& a,
    double& b,
    double& hz,
    int& i,
    int& j,
    double& res);
extern "C" void fortran_rk_adaptive_time_step(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    void* orb /* 0D_NOT_type inout */,
    int& t_dir /* 0D_NOT_integer inout */,
    double& rf_time /* 0D_NOT_real inout */,
    double& dt_try /* 0D_NOT_real inout */,
    double& dt_did /* 0D_NOT_real inout */,
    double& dt_next /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical inout */,
    void* extra_field /* 0D_NOT_type inout */);
void rk_adaptive_time_step(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orb,
    int& t_dir,
    double& rf_time,
    double& dt_try,
    double& dt_did,
    double& dt_next,
    bool& err_flag,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" void fortran_rk_time_step1(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    double& rf_time /* 0D_NOT_real in */,
    void* orb /* 0D_NOT_type inout */,
    double& dt /* 0D_NOT_real in */,
    void* new_orb /* 0D_NOT_type inout */,
    double* r_err /* 1D_NOT_real out */,
    double* dr_dt /* 1D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical inout */,
    bool* print_err /* 0D_NOT_logical inout */,
    void* extra_field /* 0D_NOT_type inout */);
FixedArray1D<Real, 10> rk_time_step1(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orb,
    double dt,
    CoordProxy& new_orb,
    std::optional<FixedArray1D<Real, 10>> dr_dt,
    bool& err_flag,
    optional_ref<bool> print_err = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" bool fortran_rotate3(
    double* vec /* 1D_NOT_real inout */,
    double& angle /* 0D_NOT_real inout */,
    double* rvec /* 1D_NOT_real inout */);
void rotate3(
    FixedArray1D<Real, 3> vec,
    double& angle,
    FixedArray1D<Real, 3> rvec);
extern "C" void fortran_rotate_em_field(
    void* field /* 0D_NOT_type inout */,
    double* w_mat /* 2D_NOT_real in */,
    double* w_inv /* 2D_NOT_real in */,
    bool* calc_dfield /* 0D_NOT_logical in */,
    bool* calc_potential /* 0D_NOT_logical in */);
void rotate_em_field(
    EmFieldProxy& field,
    FixedArray2D<Real, 3, 3> w_mat,
    FixedArray2D<Real, 3, 3> w_inv,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" void fortran_rotate_field_zx(
    void* field /* 0D_NOT_type inout */,
    double& theta /* 0D_NOT_real inout */);
void rotate_field_zx(EmFieldProxy& field, double& theta);
extern "C" void fortran_rotate_for_curved_surface(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    bool& set /* 0D_NOT_logical in */,
    double* rot_mat /* 2D_NOT_real inout */);
void rotate_for_curved_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    FixedArray2D<Real, 3, 3> rot_mat);
extern "C" void fortran_rotate_spin(
    double* rot_vec /* 1D_NOT_real in */,
    double* spin /* 1D_NOT_real inout */,
    double* qrot /* 1D_NOT_real out */);
FixedArray1D<Real, 4> rotate_spin(
    FixedArray1D<Real, 3> rot_vec,
    FixedArray1D<Real, 3> spin);
extern "C" void fortran_rotate_spin_a_step(
    void* orbit /* 0D_NOT_type inout */,
    void* field /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double& ds /* 0D_NOT_real in */);
void rotate_spin_a_step(
    CoordProxy& orbit,
    EmFieldProxy& field,
    EleProxy& ele,
    double ds);
extern "C" void fortran_rotate_spin_given_field(
    void* orbit /* 0D_NOT_type inout */,
    int& sign_z_vel /* 0D_NOT_integer in */,
    double* BL /* 1D_NOT_real in */,
    double* EL /* 1D_NOT_real in */,
    double* qrot /* 1D_NOT_real inout */);
void rotate_spin_given_field(
    CoordProxy& orbit,
    int sign_z_vel,
    std::optional<FixedArray1D<Real, 3>> BL = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> EL = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> qrot = std::nullopt);
extern "C" bool fortran_s_body_calc(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double& s_body /* 0D_NOT_real inout */);
void s_body_calc(CoordProxy& orbit, EleProxy& ele, double& s_body);
extern "C" void fortran_s_calc(void* lat /* 0D_NOT_type inout */);
void s_calc(LatProxy& lat);

// Skipped unusable routine s_ref_to_s_chord:
// - Untranslated type: csr_ele_info_struct (0D)

// Skipped unusable routine s_source_calc:
// - Untranslated type: csr_kick1_struct (0D)
// - Untranslated type: csr_struct (0D)
extern "C" void fortran_sad_mult_hard_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void sad_mult_hard_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_sad_soft_bend_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orb /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void sad_soft_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_save_a_beam_step(
    void* ele /* 0D_NOT_type in */,
    void* beam /* 0D_NOT_type in */,
    void* bunch_tracks /* 1D_ALLOC_type in */,
    double* s_body /* 0D_NOT_real in */,
    bool* is_time_coords /* 0D_NOT_logical in */);
void save_a_beam_step(
    EleProxy& ele,
    BeamProxy& beam,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt);
extern "C" void fortran_save_a_bunch_step(
    void* ele /* 0D_NOT_type in */,
    void* bunch /* 0D_NOT_type in */,
    void* bunch_track /* 0D_NOT_type in */,
    double* s_body /* 0D_NOT_real in */,
    bool* is_time_coords /* 0D_NOT_logical in */);
void save_a_bunch_step(
    EleProxy& ele,
    BunchProxy& bunch,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt);
extern "C" void fortran_save_a_step(
    void* track /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& local_ref_frame /* 0D_NOT_logical in */,
    void* orb /* 0D_NOT_type in */,
    double& s_rel /* 0D_NOT_real in */,
    bool* save_field /* 0D_NOT_logical in */,
    double* mat6 /* 2D_NOT_real in */,
    bool* make_matrix /* 0D_NOT_logical in */,
    double* rf_time /* 0D_NOT_real in */,
    void* strong_beam /* 0D_NOT_type in */);
void save_a_step(
    TrackProxy& track,
    EleProxy& ele,
    LatParamProxy& param,
    bool local_ref_frame,
    CoordProxy& orb,
    double s_rel,
    std::optional<bool> save_field = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<double> rf_time = std::nullopt,
    optional_ref<StrongBeamProxy> strong_beam = std::nullopt);
extern "C" void fortran_sbend_body_with_k1_map(
    void* ele /* 0D_NOT_type in */,
    double& dg /* 0D_NOT_real in */,
    double& b1 /* 0D_NOT_real in */,
    void* param /* 0D_NOT_type in */,
    int& n_step /* 0D_NOT_integer in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void sbend_body_with_k1_map(
    EleProxy& ele,
    double dg,
    double b1,
    LatParamProxy& param,
    int n_step,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_sc_adaptive_step(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    bool& include_image /* 0D_NOT_logical inout */,
    double& t_now /* 0D_NOT_real in */,
    double& dt_step /* 0D_NOT_real inout */,
    double& dt_next /* 0D_NOT_real out */,
    void* sc_field /* 1D_ALLOC_type in */);
double sc_adaptive_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool& include_image,
    double t_now,
    double& dt_step,
    EmFieldProxyAlloc1D& sc_field);
extern "C" void fortran_sc_step(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    bool& include_image /* 0D_NOT_logical inout */,
    double& t_end /* 0D_NOT_real in */,
    void* sc_field /* 1D_ALLOC_type in */,
    int& n_emit /* 0D_NOT_integer out */);
int sc_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool& include_image,
    double t_end,
    EmFieldProxyAlloc1D& sc_field);
extern "C" void fortran_set_active_fixer(
    void* fixer /* 0D_NOT_type inout */,
    bool* turn_on /* 0D_NOT_logical in */,
    void* orbit /* 0D_NOT_type out */);
CoordProxy set_active_fixer(
    EleProxy& fixer,
    std::optional<bool> turn_on = std::nullopt);

// Skipped unusable routine set_branch_and_ele_for_omp:
// - Untranslated type: lat_pointer_struct (1D)
extern "C" void fortran_set_custom_attribute_name(
    const char* custom_name /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */,
    int* custom_index /* 0D_NOT_integer in */);
bool set_custom_attribute_name(
    std::string custom_name,
    std::optional<int> custom_index = std::nullopt);
extern "C" void fortran_set_ele_attribute(
    void* ele /* 0D_NOT_type inout */,
    const char* set_string /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* err_print_flag /* 0D_NOT_logical in */,
    bool* set_lords /* 0D_NOT_logical in */,
    int& err_id /* 0D_NOT_integer out */);
struct SetEleAttribute {
  bool err_flag;
  int err_id;
};
Bmad::SetEleAttribute set_ele_attribute(
    EleProxy& ele,
    std::string set_string,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> set_lords = std::nullopt);
extern "C" void fortran_set_ele_defaults(
    void* ele /* 0D_NOT_type inout */,
    bool* do_allocate /* 0D_NOT_logical in */);
void set_ele_defaults(
    EleProxy& ele,
    std::optional<bool> do_allocate = std::nullopt);
extern "C" void fortran_set_ele_name(
    void* ele /* 0D_NOT_type inout */,
    const char* name /* 0D_NOT_character in */);
void set_ele_name(EleProxy& ele, std::string name);
extern "C" void fortran_set_ele_real_attribute(
    void* ele /* 0D_NOT_type inout */,
    const char* attrib_name /* 0D_NOT_character in */,
    double& value /* 0D_NOT_real in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* err_print_flag /* 0D_NOT_logical in */);
bool set_ele_real_attribute(
    EleProxy& ele,
    std::string attrib_name,
    double value,
    std::optional<bool> err_print_flag = std::nullopt);
extern "C" void fortran_set_ele_status_stale(
    void* ele /* 0D_NOT_type out */,
    int& status_group /* 0D_NOT_integer out */,
    bool& set_slaves /* 0D_NOT_logical out */);
struct SetEleStatusStale {
  EleProxy ele;
  int status_group;
  bool set_slaves;
};
Bmad::SetEleStatusStale set_ele_status_stale();
extern "C" bool fortran_set_emit_from_beam_init(
    void* beam_init_in /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    int& species /* 0D_NOT_integer in */,
    void* modes /* 0D_NOT_type in */,
    bool* err_flag /* 0D_NOT_logical in */,
    void* beam_init_set /* 0D_NOT_type inout */);
void set_emit_from_beam_init(
    BeamInitProxy& beam_init_in,
    EleProxy& ele,
    int species,
    optional_ref<NormalModesProxy> modes,
    std::optional<bool> err_flag,
    BeamInitProxy& beam_init_set);

// Skipped unusable routine set_flags_for_changed_all_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" void fortran_set_flags_for_changed_integer_attribute(
    void* ele /* 0D_NOT_type in */,
    int& attrib /* 0D_NOT_integer inout */,
    bool* set_dependent /* 0D_NOT_logical in */);
void set_flags_for_changed_attribute(
    EleProxy& ele,
    int& attrib,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_lat_attribute(
    void* lat /* 0D_NOT_type inout */,
    bool* set_dependent /* 0D_NOT_logical in */);
void set_flags_for_changed_attribute(
    LatProxy& lat,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_logical_attribute(
    void* ele /* 0D_NOT_type in */,
    bool& attrib /* 0D_NOT_logical inout */,
    bool* set_dependent /* 0D_NOT_logical in */);
void set_flags_for_changed_attribute(
    EleProxy& ele,
    bool& attrib,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_real_attribute(
    void* ele /* 0D_NOT_type in */,
    double* attrib /* 0D_NOT_real inout */,
    bool* set_dependent /* 0D_NOT_logical in */);
void set_flags_for_changed_attribute(
    EleProxy& ele,
    optional_ref<double> attrib = std::nullopt,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_fringe_on_off(
    double& fringe_at /* 0D_NOT_real inout */,
    int& ele_end /* 0D_NOT_integer in */,
    int& on_or_off /* 0D_NOT_integer in */);
void set_fringe_on_off(double& fringe_at, int ele_end, int on_or_off);
extern "C" void fortran_set_lords_status_stale(
    void* ele /* 0D_NOT_type in */,
    int& stat_group /* 0D_NOT_integer in */,
    bool* control_bookkeeping /* 0D_NOT_logical in */,
    int* flag /* 0D_NOT_integer in */);
void set_lords_status_stale(
    EleProxy& ele,
    int stat_group,
    std::optional<bool> control_bookkeeping = std::nullopt,
    std::optional<int> flag = std::nullopt);
extern "C" void fortran_set_on_off(
    int& key /* 0D_NOT_integer in */,
    void* lat /* 0D_NOT_type inout */,
    int& switch_ /* 0D_NOT_integer in */,
    void* orb /* 1D_ALLOC_type in */,
    bool* use_ref_orb /* 0D_NOT_logical in */,
    int* ix_branch /* 0D_NOT_integer in */,
    void* saved_values /* 1D_ALLOC_real inout */,
    const char* attribute /* 0D_NOT_character in */,
    int* set_val /* 0D_NOT_integer in */);
void set_on_off(
    int key,
    LatProxy& lat,
    int switch_,
    optional_ref<CoordProxyAlloc1D> orb = std::nullopt,
    std::optional<bool> use_ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<RealAlloc1D> saved_values = std::nullopt,
    std::optional<std::string> attribute = std::nullopt,
    std::optional<int> set_val = std::nullopt);
extern "C" void fortran_set_orbit_to_zero(
    void* orbit /* 1D_ALLOC_type out */,
    int& n1 /* 0D_NOT_integer in */,
    int& n2 /* 0D_NOT_integer in */,
    int* ix_noset /* 0D_NOT_integer in */);
CoordProxyAlloc1D set_orbit_to_zero(
    int n1,
    int n2,
    std::optional<int> ix_noset = std::nullopt);
extern "C" void fortran_set_ptc(
    double* e_tot /* 0D_NOT_real in */,
    int* particle /* 0D_NOT_integer in */,
    int* taylor_order /* 0D_NOT_integer in */,
    int* integ_order /* 0D_NOT_integer in */,
    int* n_step /* 0D_NOT_integer in */,
    bool* no_cavity /* 0D_NOT_logical in */,
    bool* force_init /* 0D_NOT_logical in */);
void set_ptc(
    std::optional<double> e_tot = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> taylor_order = std::nullopt,
    std::optional<int> integ_order = std::nullopt,
    std::optional<int> n_step = std::nullopt,
    std::optional<bool> no_cavity = std::nullopt,
    std::optional<bool> force_init = std::nullopt);
extern "C" void fortran_set_ptc_base_state(
    const char* component /* 0D_NOT_character in */,
    bool& set_val /* 0D_NOT_logical in */,
    bool& old_val /* 0D_NOT_logical out */);
bool set_ptc_base_state(std::string component, bool set_val);
extern "C" void fortran_set_ptc_com_pointers();
void set_ptc_com_pointers();
extern "C" void fortran_set_ptc_quiet(
    int& channel /* 0D_NOT_integer in */,
    bool& set /* 0D_NOT_logical in */,
    int& old_val /* 0D_NOT_integer inout */);
void set_ptc_quiet(int channel, bool set, int& old_val);
extern "C" void fortran_set_ptc_verbose(bool& on /* 0D_NOT_logical inout */);
void set_ptc_verbose(bool& on);
extern "C" void fortran_set_pwd_ele(
    void* lat /* 0D_NOT_type in */,
    void* mode0 /* 0D_NOT_type in */,
    double& inductance /* 0D_NOT_real in */);
void set_pwd_ele(LatProxy& lat, NormalModesProxy& mode0, double inductance);
extern "C" void fortran_set_status_flags(
    void* bookkeeping_state /* 0D_NOT_type out */,
    int& stat /* 0D_NOT_integer in */);
BookkeepingStateProxy set_status_flags(int stat);
extern "C" bool fortran_set_tune(
    double& phi_a_set /* 0D_NOT_real in */,
    double& phi_b_set /* 0D_NOT_real in */,
    void* dk1 /* 1D_ALLOC_real in */,
    void* eles /* 1D_ALLOC_type in */,
    void* branch /* 0D_NOT_type inout */,
    void* orb /* 1D_ALLOC_type inout */,
    bool* print_err /* 0D_NOT_logical in */,
    bool& ok /* 0D_NOT_logical inout */);
void set_tune(
    double phi_a_set,
    double phi_b_set,
    RealAlloc1D& dk1,
    ElePointerProxyAlloc1D& eles,
    BranchProxy& branch,
    CoordProxyAlloc1D& orb,
    std::optional<bool> print_err,
    bool& ok);

// Skipped unusable routine set_tune_via_group_knobs:
// - Routine in configuration skip list
extern "C" void fortran_set_twiss(
    void* branch /* 0D_NOT_type in */,
    void* twiss_ele /* 0D_NOT_type in */,
    int& ix_ele /* 0D_NOT_integer in */,
    bool& match_deta_ds /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */);
void set_twiss(
    BranchProxy& branch,
    EleProxy& twiss_ele,
    int ix_ele,
    bool match_deta_ds,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_set_z_tune(
    void* branch /* 0D_NOT_type inout */,
    double& z_tune /* 0D_NOT_real in */,
    bool& ok /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */);
bool set_z_tune(
    BranchProxy& branch,
    double z_tune,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_settable_dep_var_bookkeeping(
    void* ele /* 0D_NOT_type inout */);
void settable_dep_var_bookkeeping(EleProxy& ele);
extern "C" void fortran_setup_high_energy_space_charge_calc(
    bool& calc_on /* 0D_NOT_logical in */,
    void* branch /* 0D_NOT_type in */,
    double& n_part /* 0D_NOT_real in */,
    void* mode /* 0D_NOT_type in */,
    void* closed_orb /* 1D_ALLOC_type in */);
void setup_high_energy_space_charge_calc(
    bool calc_on,
    BranchProxy& branch,
    double n_part,
    NormalModesProxy& mode,
    optional_ref<CoordProxyAlloc1D> closed_orb = std::nullopt);

// Skipped unusable routine sfft:
// - Routine in configuration skip list
extern "C" void fortran_sigma_mat_ptc_to_bmad(
    double* sigma_mat_ptc /* 2D_NOT_real in */,
    double& beta0 /* 0D_NOT_real in */,
    double* sigma_mat_bmad /* 2D_NOT_real out */);
FixedArray2D<Real, 6, 6> sigma_mat_ptc_to_bmad(
    FixedArray2D<Real, 6, 6> sigma_mat_ptc,
    double beta0);
extern "C" bool fortran_significant_difference(
    double& value1 /* 0D_NOT_real in */,
    double& value2 /* 0D_NOT_real in */,
    double* abs_tol /* 0D_NOT_real in */,
    double* rel_tol /* 0D_NOT_real in */,
    bool& is_different /* 0D_NOT_logical inout */);
void significant_difference(
    double value1,
    double value2,
    std::optional<double> abs_tol,
    std::optional<double> rel_tol,
    bool& is_different);
extern "C" bool fortran_skip_ele_blender(
    void* ele /* 0D_NOT_type inout */,
    bool& skip /* 0D_NOT_logical inout */);
void skip_ele_blender(EleProxy& ele, bool& skip);
extern "C" void fortran_slice_lattice(
    void* lat /* 0D_NOT_type inout */,
    const char* ele_list /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    bool* do_bookkeeping /* 0D_NOT_logical in */);
bool slice_lattice(
    LatProxy& lat,
    std::string ele_list,
    std::optional<bool> do_bookkeeping = std::nullopt);
extern "C" void fortran_soft_quadrupole_edge_kick(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    int& particle_at /* 0D_NOT_integer in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void soft_quadrupole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_sol_quad_mat6_calc(
    double& ks_in /* 0D_NOT_real inout */,
    double& k1_in /* 0D_NOT_real inout */,
    double& tilt /* 0D_NOT_real in */,
    double& length /* 0D_NOT_real in */,
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void sol_quad_mat6_calc(
    double& ks_in,
    double& k1_in,
    double tilt,
    double length,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);

// Skipped unusable routine solenoid_track_and_mat:
// - Variable inout sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_solve_psi_adaptive(
    double& t0 /* 0D_NOT_real in */,
    double& t1 /* 0D_NOT_real in */,
    double& p0 /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    double& p1 /* 0D_NOT_real out */);
double solve_psi_adaptive(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args);
extern "C" void fortran_solve_psi_fixed_steps(
    double& t0 /* 0D_NOT_real in */,
    double& t1 /* 0D_NOT_real in */,
    double& p0 /* 0D_NOT_real in */,
    double* args /* 1D_NOT_real in */,
    void* t /* 1D_ALLOC_real out */,
    void* p /* 1D_ALLOC_real out */);
struct SolvePsiFixedSteps {
  RealAlloc1D t;
  RealAlloc1D p;
};
Bmad::SolvePsiFixedSteps solve_psi_fixed_steps(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args);
extern "C" void fortran_sort_complex_taylor_terms(
    void* complex_taylor_in /* 0D_NOT_type in */,
    void* complex_taylor_sorted /* 0D_NOT_type out */);
ComplexTaylorProxy sort_complex_taylor_terms(
    ComplexTaylorProxy& complex_taylor_in);

// Skipped unusable routine sort_universal_terms:
// - Untranslated type: universal_taylor (0D)
// - Untranslated type: universal_taylor (0D)

// Skipped unusable routine space_charge_3d:
// - Untranslated type: mesh3d_struct (0D)
// - Variable in sized array: image_efield(:,:,:,:) 4D_ALLOC_real

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
    double* mat_1turn /* 2D_NOT_real in */,
    double* dn_dpz_partial /* 2D_NOT_real in */,
    bool& error /* 0D_NOT_logical out */,
    double* dn_dpz /* 1D_NOT_real inout */);
bool spin_dn_dpz_from_mat8(
    FixedArray2D<Real, 8, 8> mat_1turn,
    std::optional<FixedArray2D<Real, 3, 3>> dn_dpz_partial,
    FixedArray1D<Real, 3> dn_dpz);
extern "C" bool fortran_spin_dn_dpz_from_qmap(
    double* orb_mat /* 2D_NOT_real in */,
    double* q_map /* 2D_NOT_real in */,
    double* dn_dpz_partial /* 2D_NOT_real in */,
    double* dn_dpz_partial2 /* 2D_NOT_real in */,
    bool& error /* 0D_NOT_logical out */,
    double* n0 /* 1D_NOT_real in */,
    double* dn_dpz /* 1D_NOT_real inout */);
bool spin_dn_dpz_from_qmap(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> q_map,
    FixedArray2D<Real, 3, 3> dn_dpz_partial,
    FixedArray2D<Real, 3, 3> dn_dpz_partial2,
    std::optional<FixedArray1D<Real, 3>> n0,
    FixedArray1D<Real, 3> dn_dpz);
extern "C" void fortran_spin_map1_normalize(
    double* spin1 /* 2D_NOT_real inout */);
void spin_map1_normalize(FixedArray2D<Real, 4, 7> spin1);
extern "C" void fortran_spin_mat8_resonance_strengths(
    std::complex<double>* orb_evec /* 1D_NOT_complex in */,
    double* mat8 /* 2D_NOT_real in */,
    double& xi_sum /* 0D_NOT_real out */,
    double& xi_diff /* 0D_NOT_real out */);
struct SpinMat8ResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
Bmad::SpinMat8ResonanceStrengths spin_mat8_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 6, 6> mat8);
extern "C" void fortran_spin_mat_to_eigen(
    double* orb_mat /* 2D_NOT_real in */,
    double* spin_map /* 2D_NOT_real in */,
    std::complex<double>* orb_eval /* 1D_NOT_complex out */,
    std::complex<double>* orb_evec /* 2D_NOT_complex out */,
    double* n0 /* 1D_NOT_real out */,
    std::complex<double>* spin_evec /* 2D_NOT_complex out */,
    bool& error /* 0D_NOT_logical out */);
struct SpinMatToEigen {
  FixedArray1D<Complex, 6> orb_eval;
  FixedArray2D<Complex, 6, 6> orb_evec;
  FixedArray1D<Real, 3> n0;
  FixedArray2D<Complex, 6, 3> spin_evec;
  bool error;
};
Bmad::SpinMatToEigen spin_mat_to_eigen(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> spin_map);
extern "C" bool fortran_spin_omega(
    void* field /* 0D_NOT_type inout */,
    void* orbit /* 0D_NOT_type inout */,
    int& sign_z_vel /* 0D_NOT_integer inout */,
    bool* phase_space_coords /* 0D_NOT_logical inout */,
    double* omega /* 1D_NOT_real inout */);
void spin_omega(
    EmFieldProxy& field,
    CoordProxy& orbit,
    int& sign_z_vel,
    optional_ref<bool> phase_space_coords,
    FixedArray1D<Real, 3> omega);
extern "C" void fortran_spin_quat_resonance_strengths(
    std::complex<double>* orb_evec /* 1D_NOT_complex in */,
    double* spin_q /* 2D_NOT_real in */,
    double& xi_sum /* 0D_NOT_real out */,
    double& xi_diff /* 0D_NOT_real out */);
struct SpinQuatResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
Bmad::SpinQuatResonanceStrengths spin_quat_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 4, 7> spin_q);
extern "C" bool fortran_spin_taylor_to_linear(
    void* spin_taylor /* 1D_NOT_type in */,
    bool& normalize /* 0D_NOT_logical in */,
    double* dref_orb /* 1D_NOT_real in */,
    bool& is_on /* 0D_NOT_logical in */,
    double* spin_map1 /* 2D_NOT_real inout */);
void spin_taylor_to_linear(
    FixedArray1D<TaylorProxy, 4> spin_taylor,
    bool normalize,
    FixedArray1D<Real, 6> dref_orb,
    bool is_on,
    FixedArray2D<Real, 4, 7> spin_map1);
extern "C" bool fortran_spinor_to_polar(
    std::complex<double>* spinor /* 1D_NOT_complex in */,
    void* polar /* 0D_NOT_type inout */);
void spinor_to_polar(FixedArray1D<Complex, 2> spinor, SpinPolarProxy& polar);
extern "C" bool fortran_spinor_to_vec(
    std::complex<double>* spinor /* 1D_NOT_complex in */,
    double* vec /* 1D_NOT_real inout */);
void spinor_to_vec(FixedArray1D<Complex, 2> spinor, FixedArray1D<Real, 3> vec);
extern "C" void fortran_spline_fit_orbit(
    void* start_orb /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type in */,
    double* spline_x /* 1D_NOT_real in */,
    double* spline_y /* 1D_NOT_real in */);
void spline_fit_orbit(
    CoordProxy& start_orb,
    CoordProxy& end_orb,
    FixedArray1D<Real, 4> spline_x,
    FixedArray1D<Real, 4> spline_y);

// Skipped unusable routine split_expression_string:
// - Variable-sized out character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_split_lat(
    void* lat /* 0D_NOT_type inout */,
    double& s_split /* 0D_NOT_real in */,
    int& ix_branch /* 0D_NOT_integer in */,
    int& ix_split /* 0D_NOT_integer out */,
    bool& split_done /* 0D_NOT_logical out */,
    bool* add_suffix /* 0D_NOT_logical in */,
    bool* check_sanity /* 0D_NOT_logical in */,
    bool* save_null_drift /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* choose_max /* 0D_NOT_logical in */,
    int* ix_insert /* 0D_NOT_integer in */);
struct SplitLat {
  int ix_split;
  bool split_done;
  bool err_flag;
};
Bmad::SplitLat split_lat(
    LatProxy& lat,
    double s_split,
    int ix_branch,
    std::optional<bool> add_suffix = std::nullopt,
    std::optional<bool> check_sanity = std::nullopt,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> choose_max = std::nullopt,
    std::optional<int> ix_insert = std::nullopt);
extern "C" void fortran_sprint_spin_taylor_map(
    void* ele /* 0D_NOT_type inout */,
    double* start_orbit /* 1D_NOT_real in */);
void sprint_spin_taylor_map(
    EleProxy& ele,
    std::optional<FixedArray1D<Real, 6>> start_orbit = std::nullopt);
extern "C" void fortran_sr_longitudinal_wake_particle(
    void* ele /* 0D_NOT_type inout */,
    void* orbit /* 0D_NOT_type inout */);
void sr_longitudinal_wake_particle(EleProxy& ele, CoordProxy& orbit);
extern "C" void fortran_sr_transverse_wake_particle(
    void* ele /* 0D_NOT_type inout */,
    void* orbit /* 0D_NOT_type inout */);
void sr_transverse_wake_particle(EleProxy& ele, CoordProxy& orbit);
extern "C" void fortran_sr_z_long_wake(
    void* ele /* 0D_NOT_type in */,
    void* bunch /* 0D_NOT_type inout */,
    double& z_ave /* 0D_NOT_real in */);
void sr_z_long_wake(EleProxy& ele, BunchProxy& bunch, double z_ave);
extern "C" void fortran_srdt_calc(
    void* lat /* 0D_NOT_type in */,
    void* srdt_sums /* 0D_NOT_type out */,
    int& order /* 0D_NOT_integer in */,
    int* n_slices_gen_opt /* 0D_NOT_integer in */,
    int* n_slices_sxt_opt /* 0D_NOT_integer in */,
    void* per_ele_out /* 1D_ALLOC_type inout */);
SummationRdtProxy srdt_calc(
    LatProxy& lat,
    int order,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    optional_ref<SummationRdtProxyAlloc1D> per_ele_out = std::nullopt);

// Skipped unusable routine srdt_calc_with_cache:
// - Variable inout sized array: cache(:,:,:) 3D_ALLOC_complex
extern "C" void fortran_srdt_lsq_solution(
    void* lat /* 0D_NOT_type in */,
    void* var_indexes /* 1D_ALLOC_integer in */,
    void* ls_soln /* 1D_ALLOC_real out */,
    int* n_slices_gen_opt /* 0D_NOT_integer in */,
    int* n_slices_sxt_opt /* 0D_NOT_integer in */,
    double* chrom_set_x_opt /* 0D_NOT_real in */,
    double* chrom_set_y_opt /* 0D_NOT_real in */,
    double* weight_in /* 1D_NOT_real in */);
RealAlloc1D srdt_lsq_solution(
    LatProxy& lat,
    IntAlloc1D& var_indexes,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    std::optional<double> chrom_set_x_opt = std::nullopt,
    std::optional<double> chrom_set_y_opt = std::nullopt,
    std::optional<FixedArray1D<Real, 10>> weight_in = std::nullopt);
extern "C" void fortran_start_branch_at(
    void* lat /* 0D_NOT_type inout */,
    const char* ele_start /* 0D_NOT_character in */,
    bool& move_end_marker /* 0D_NOT_logical in */,
    bool& error /* 0D_NOT_logical out */);
bool start_branch_at(
    LatProxy& lat,
    std::string ele_start,
    bool move_end_marker);
extern "C" bool fortran_stream_ele_end(
    int& physical_end /* 0D_NOT_integer in */,
    int& ele_orientation /* 0D_NOT_integer in */,
    int& stream_end /* 0D_NOT_integer inout */);
void stream_ele_end(int physical_end, int ele_orientation, int& stream_end);
extern "C" void fortran_string_attrib(
    const char* attrib_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    const char* attrib_value /* 0D_NOT_character out */);
std::string string_attrib(std::string attrib_name, EleProxy& ele);
extern "C" void fortran_strong_beam_sigma_calc(
    void* ele /* 0D_NOT_type in */,
    double& s_pos /* 0D_NOT_real in */,
    double* sigma /* 1D_NOT_real out */,
    double& bbi_const /* 0D_NOT_real out */,
    double* dsigma_ds /* 1D_NOT_real out */);
struct StrongBeamSigmaCalc {
  FixedArray1D<Real, 2> sigma;
  double bbi_const;
  FixedArray1D<Real, 2> dsigma_ds;
};
Bmad::StrongBeamSigmaCalc strong_beam_sigma_calc(EleProxy& ele, double s_pos);
extern "C" bool fortran_strong_beam_strength(
    void* ele /* 0D_NOT_type in */,
    double& strength /* 0D_NOT_real inout */);
void strong_beam_strength(EleProxy& ele, double& strength);
extern "C" void fortran_surface_grid_displacement(
    void* ele /* 0D_NOT_type in */,
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical in */,
    double& z /* 0D_NOT_real in */,
    double* dz_dxy /* 1D_NOT_real in */,
    bool* extend_grid /* 0D_NOT_logical in */);
void surface_grid_displacement(
    EleProxy& ele,
    double& x,
    double& y,
    bool err_flag,
    double z,
    std::optional<FixedArray1D<Real, 2>> dz_dxy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt);

// Skipped unusable routine switch_attrib_value_name:
// - Variable-sized out character array: name_list(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_symp_lie_bmad(
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    void* track /* 0D_NOT_type out */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */,
    bool* offset_ele /* 0D_NOT_logical in */);
TrackProxy symp_lie_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> offset_ele = std::nullopt);
extern "C" void fortran_t6_to_b123(
    double* t6 /* 2D_NOT_real in */,
    double* abz_tunes /* 1D_NOT_real in */,
    double* B1 /* 2D_NOT_real out */,
    double* B2 /* 2D_NOT_real out */,
    double* B3 /* 2D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct T6ToB123 {
  FixedArray2D<Real, 6, 6> B1;
  FixedArray2D<Real, 6, 6> B2;
  FixedArray2D<Real, 6, 6> B3;
  bool err_flag;
};
Bmad::T6ToB123 t6_to_b123(
    FixedArray2D<Real, 6, 6> t6,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_taper_mag_strengths(
    void* lat /* 0D_NOT_type inout */,
    void* ref_lat /* 0D_NOT_type in */,
    const char* except /* 0D_NOT_character in */,
    bool* err_flag /* 0D_NOT_logical inout */);
void taper_mag_strengths(
    LatProxy& lat,
    optional_ref<LatProxy> ref_lat = std::nullopt,
    std::optional<std::string> except = std::nullopt,
    optional_ref<bool> err_flag = std::nullopt);
extern "C" void fortran_target_min_max_calc(
    double* r_corner1 /* 1D_NOT_real in */,
    double* r_corner2 /* 1D_NOT_real in */,
    double& y_min /* 0D_NOT_real inout */,
    double& y_max /* 0D_NOT_real inout */,
    double& phi_min /* 0D_NOT_real inout */,
    double& phi_max /* 0D_NOT_real inout */,
    bool* initial /* 0D_NOT_logical in */);
void target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double& y_min,
    double& y_max,
    double& phi_min,
    double& phi_max,
    std::optional<bool> initial = std::nullopt);
extern "C" void fortran_target_rot_mats(
    double* r_center /* 1D_NOT_real in */,
    double* w_to_target /* 2D_NOT_real out */,
    double* w_to_ele /* 2D_NOT_real out */);
struct TargetRotMats {
  FixedArray2D<Real, 3, 3> w_to_target;
  FixedArray2D<Real, 3, 3> w_to_ele;
};
Bmad::TargetRotMats target_rot_mats(FixedArray1D<Real, 3> r_center);
extern "C" void fortran_taylor_equal_taylor(
    void* taylor1 /* 0D_NOT_type out */,
    void* taylor2 /* 0D_NOT_type in */);
TaylorProxy taylor_equal_taylor(TaylorProxy& taylor2);
extern "C" void fortran_taylor_inverse(
    void* taylor_in /* 1D_ALLOC_type in */,
    void* taylor_inv /* 1D_ALLOC_type out */,
    bool& err /* 0D_NOT_logical out */);
struct TaylorInverse {
  TaylorProxyAlloc1D taylor_inv;
  bool err;
};
Bmad::TaylorInverse taylor_inverse(TaylorProxyAlloc1D& taylor_in);

// Skipped unusable routine taylor_minus_taylor:
// - Routine in configuration skip list

// Skipped unusable routine taylor_plus_taylor:
// - Routine in configuration skip list
extern "C" void fortran_taylor_propagate1(
    void* orb_taylor /* 1D_ALLOC_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* ref_in /* 0D_NOT_type in */,
    void* spin_taylor /* 1D_ALLOC_type inout */);
bool taylor_propagate1(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> ref_in = std::nullopt,
    optional_ref<TaylorProxyAlloc1D> spin_taylor = std::nullopt);

// Skipped unusable routine taylor_to_genfield:
// - Untranslated type: genfield (0D)
extern "C" void fortran_taylor_to_mad_map(
    void* taylor /* 1D_ALLOC_type in */,
    void* energy /* 0D_NOT_type in */,
    void* map /* 0D_NOT_type out */);
MadMapProxy taylor_to_mad_map(
    TaylorProxyAlloc1D& taylor,
    MadEnergyProxy& energy);

// Skipped unusable routine taylor_to_real_8:
// - Untranslated type: real_8 (1D)
extern "C" void fortran_taylors_equal_taylors(
    void* taylor1 /* 1D_ALLOC_type out */,
    void* taylor2 /* 1D_ALLOC_type in */);
TaylorProxyAlloc1D taylors_equal_taylors(TaylorProxyAlloc1D& taylor2);
extern "C" void fortran_tilt_coords(
    double& tilt_val /* 0D_NOT_real in */,
    void* coord /* 1D_ALLOC_real inout */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void tilt_coords(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_tilt_coords_photon(
    double& tilt_val /* 0D_NOT_real in */,
    void* coord /* 1D_ALLOC_real inout */,
    double* w_mat /* 2D_NOT_real inout */);
void tilt_coords_photon(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat = std::nullopt);
extern "C" void fortran_tilt_mat6(
    double* mat6 /* 2D_NOT_real inout */,
    double& tilt /* 0D_NOT_real in */);
void tilt_mat6(FixedArray2D<Real, 6, 6> mat6, double tilt);

// Skipped unusable routine time_runge_kutta_periodic_kick_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_to_eta_reading(
    void* eta_actual /* 1D_ALLOC_real in */,
    void* ele /* 0D_NOT_type in */,
    int& axis /* 0D_NOT_integer in */,
    bool& add_noise /* 0D_NOT_logical in */,
    double& reading /* 0D_NOT_real out */,
    bool& err /* 0D_NOT_logical out */);
struct ToEtaReading {
  double reading;
  bool err;
};
Bmad::ToEtaReading to_eta_reading(
    RealAlloc1D& eta_actual,
    EleProxy& ele,
    int axis,
    bool add_noise);
extern "C" void fortran_to_fieldmap_coords(
    void* ele /* 0D_NOT_type in */,
    void* local_orb /* 0D_NOT_type in */,
    double& s_body /* 0D_NOT_real in */,
    int& ele_anchor_pt /* 0D_NOT_integer in */,
    double* r0 /* 1D_NOT_real in */,
    bool& curved_ref_frame /* 0D_NOT_logical in */,
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& cos_ang /* 0D_NOT_real inout */,
    double& sin_ang /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical in */);
void to_fieldmap_coords(
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
    bool err_flag);
extern "C" void fortran_to_orbit_reading(
    void* orb /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    int& axis /* 0D_NOT_integer in */,
    bool& add_noise /* 0D_NOT_logical in */,
    double& reading /* 0D_NOT_real out */,
    bool& err /* 0D_NOT_logical out */);
struct ToOrbitReading {
  double reading;
  bool err;
};
Bmad::ToOrbitReading to_orbit_reading(
    CoordProxy& orb,
    EleProxy& ele,
    int axis,
    bool add_noise);
extern "C" void fortran_to_phase_and_coupling_reading(
    void* ele /* 0D_NOT_type in */,
    bool& add_noise /* 0D_NOT_logical in */,
    void* reading /* 0D_NOT_type out */,
    bool& err /* 0D_NOT_logical out */);
struct ToPhaseAndCouplingReading {
  BpmPhaseCouplingProxy reading;
  bool err;
};
Bmad::ToPhaseAndCouplingReading to_phase_and_coupling_reading(
    EleProxy& ele,
    bool add_noise);
extern "C" bool fortran_to_photon_angle_coords(
    void* orb_in /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* orb_out /* 0D_NOT_type out */);
CoordProxy to_photon_angle_coords(CoordProxy& orb_in, EleProxy& ele);
extern "C" void fortran_to_surface_coords(
    void* lab_orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* surface_orbit /* 0D_NOT_type out */);
CoordProxy to_surface_coords(CoordProxy& lab_orbit, EleProxy& ele);
extern "C" void fortran_touschek_lifetime(
    void* mode /* 0D_NOT_type in */,
    double& Tl /* 0D_NOT_real out */,
    void* lat /* 0D_NOT_type in */);
double touschek_lifetime(NormalModesProxy& mode, LatProxy& lat);

// Skipped unusable routine touschek_lifetime_ele_by_ele:
// - Untranslated type: momentum_aperture_struct (1D)

// Skipped unusable routine touschek_lifetime_with_aperture:
// - Untranslated type: momentum_aperture_struct (1D)
extern "C" void fortran_touschek_rate1(
    void* mode /* 0D_NOT_type in */,
    double& rate /* 0D_NOT_real out */,
    void* lat /* 0D_NOT_type in */,
    int* ix /* 0D_NOT_integer in */,
    double* s /* 0D_NOT_real in */);
double touschek_rate1(
    NormalModesProxy& mode,
    LatProxy& lat,
    std::optional<int> ix = std::nullopt,
    std::optional<double> s = std::nullopt);
extern "C" void fortran_touschek_rate1_zap(
    void* mode /* 0D_NOT_type inout */,
    double& rate /* 0D_NOT_real inout */,
    void* lat /* 0D_NOT_type inout */,
    int* ix /* 0D_NOT_integer inout */,
    double* s /* 0D_NOT_real inout */);
void touschek_rate1_zap(
    NormalModesProxy& mode,
    double& rate,
    LatProxy& lat,
    optional_ref<int> ix = std::nullopt,
    optional_ref<double> s = std::nullopt);
extern "C" void fortran_track1(
    void* start_orb /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    void* track /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* ignore_radiation /* 0D_NOT_logical in */,
    bool* make_map1 /* 0D_NOT_logical in */,
    bool* init_to_edge /* 0D_NOT_logical in */);
struct Track1 {
  CoordProxy end_orb;
  bool err_flag;
};
Bmad::Track1 track1(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<TrackProxy> track = std::nullopt,
    std::optional<bool> ignore_radiation = std::nullopt,
    std::optional<bool> make_map1 = std::nullopt,
    std::optional<bool> init_to_edge = std::nullopt);
extern "C" void fortran_track1_beam(
    void* beam /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    void* centroid /* 1D_ALLOC_type in */,
    int* direction /* 0D_NOT_integer in */);
bool track1_beam(
    BeamProxy& beam,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt);
extern "C" void fortran_track1_bmad(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type out */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
struct Track1Bmad {
  bool err_flag;
  TrackProxy track;
};
Bmad::Track1Bmad track1_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_bmad_photon(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool track1_bmad_photon(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_bunch(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    void* centroid /* 1D_ALLOC_type in */,
    int* direction /* 0D_NOT_integer in */,
    void* bunch_track /* 0D_NOT_type inout */);
bool track1_bunch(
    BunchProxy& bunch,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_csr(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* centroid /* 1D_ALLOC_type in */,
    bool& err /* 0D_NOT_logical out */,
    double* s_start /* 0D_NOT_real in */,
    double* s_end /* 0D_NOT_real in */,
    void* bunch_track /* 0D_NOT_type inout */);
bool track1_bunch_csr(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_csr3d(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* centroid /* 1D_ALLOC_type in */,
    bool& err /* 0D_NOT_logical out */,
    double* s_start /* 0D_NOT_real in */,
    double* s_end /* 0D_NOT_real in */,
    void* bunch_track /* 0D_NOT_type inout */);
bool track1_bunch_csr3d(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_hom(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    int* direction /* 0D_NOT_integer in */,
    void* bunch_track /* 0D_NOT_type inout */);
void track1_bunch_hom(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);

// Skipped unusable routine track1_bunch_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_bunch_space_charge(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    bool* track_to_same_s /* 0D_NOT_logical in */,
    void* bunch_track /* 0D_NOT_type inout */);
bool track1_bunch_space_charge(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<bool> track_to_same_s = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_crystal(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_crystal(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);

// Skipped unusable routine track1_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_diffraction_plate_or_mask(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_diffraction_plate_or_mask(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_high_energy_space_charge(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_high_energy_space_charge(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_lens(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_lens(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_linear(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type inout */);
void track1_linear(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_lr_wake(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */);
void track1_lr_wake(BunchProxy& bunch, EleProxy& ele);
extern "C" void fortran_track1_mad(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */);
void track1_mad(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_mirror(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_mirror(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_mosaic_crystal(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_mosaic_crystal(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_multilayer_mirror(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_multilayer_mirror(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);

// Skipped unusable routine track1_postprocess_def:
// - Routine in configuration skip list

// Skipped unusable routine track1_preprocess_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_radiation(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    int& edge /* 0D_NOT_integer in */);
void track1_radiation(CoordProxy& orbit, EleProxy& ele, int edge);
extern "C" void fortran_track1_radiation_center(
    void* orbit /* 0D_NOT_type inout */,
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    bool* rad_damp /* 0D_NOT_logical in */,
    bool* rad_fluct /* 0D_NOT_logical in */);
void track1_radiation_center(
    CoordProxy& orbit,
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> rad_damp = std::nullopt,
    std::optional<bool> rad_fluct = std::nullopt);
extern "C" void fortran_track1_runge_kutta(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type out */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
struct Track1RungeKutta {
  bool err_flag;
  TrackProxy track;
};
Bmad::Track1RungeKutta track1_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_sample(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */);
void track1_sample(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_spin(
    void* start_orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type out */,
    void* param /* 0D_NOT_type inout */,
    void* end_orb /* 0D_NOT_type out */,
    bool* make_quaternion /* 0D_NOT_logical inout */);
struct Track1Spin {
  EleProxy ele;
  CoordProxy end_orb;
};
Bmad::Track1Spin track1_spin(
    CoordProxy& start_orb,
    LatParamProxy& param,
    optional_ref<bool> make_quaternion = std::nullopt);

// Skipped unusable routine track1_spin_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_track1_spin_integration(
    void* start_orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    void* end_orb /* 0D_NOT_type out */);
CoordProxy track1_spin_integration(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_spin_taylor(
    void* start_orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type inout */,
    void* end_orb /* 0D_NOT_type out */);
CoordProxy track1_spin_taylor(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_sr_wake(
    void* bunch /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */);
void track1_sr_wake(BunchProxy& bunch, EleProxy& ele);
extern "C" void fortran_track1_symp_lie_ptc(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* track /* 0D_NOT_type out */);
TrackProxy track1_symp_lie_ptc(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_taylor(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* taylor /* 1D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track1_taylor(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<FixedArray1D<TaylorProxy, 6>> taylor = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_time_runge_kutta(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type out */,
    double* t_end /* 0D_NOT_real in */,
    double* dt_step /* 0D_NOT_real inout */);
struct Track1TimeRungeKutta {
  bool err_flag;
  TrackProxy track;
};
Bmad::Track1TimeRungeKutta track1_time_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<double> t_end = std::nullopt,
    optional_ref<double> dt_step = std::nullopt);

// Skipped unusable routine track1_wake_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track_a_beambeam(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    void* track /* 0D_NOT_type out */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
struct TrackABeambeam {
  TrackProxy track;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
Bmad::TrackABeambeam track_a_beambeam(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_bend(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_bend(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_bend_photon(
    void* orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    double& length /* 0D_NOT_real in */);
void track_a_bend_photon(CoordProxy& orb, EleProxy& ele, double length);
extern "C" void fortran_track_a_capillary(
    void* orb /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */);
void track_a_capillary(CoordProxy& orb, EleProxy& ele);
extern "C" void fortran_track_a_converter(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_converter(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_crab_cavity(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_crab_cavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_drift(
    void* orb /* 0D_NOT_type inout */,
    double& length /* 0D_NOT_real in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */,
    int* ele_orientation /* 0D_NOT_integer in */,
    bool* include_ref_motion /* 0D_NOT_logical in */,
    double* time /* 0D_NOT_real inout */);
void track_a_drift(
    CoordProxy& orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<int> ele_orientation = std::nullopt,
    std::optional<bool> include_ref_motion = std::nullopt,
    optional_ref<double> time = std::nullopt);
extern "C" void fortran_track_a_drift_photon(
    void* orb /* 0D_NOT_type inout */,
    double& length /* 0D_NOT_real in */,
    bool& phase_relative_to_ref /* 0D_NOT_logical in */);
void track_a_drift_photon(
    CoordProxy& orb,
    double length,
    bool phase_relative_to_ref);
extern "C" void fortran_track_a_foil(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_foil(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_gkicker(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_gkicker(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_lcavity(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_lcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_lcavity_old(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_lcavity_old(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_mask(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_mask(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_match(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool* err_flag /* 0D_NOT_logical inout */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_match(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_patch(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    bool* drift_to_exit /* 0D_NOT_logical in */,
    double& s_ent /* 0D_NOT_real out */,
    double& ds_ref /* 0D_NOT_real out */,
    bool* track_spin /* 0D_NOT_logical in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
struct TrackAPatch {
  double s_ent;
  double ds_ref;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
Bmad::TrackAPatch track_a_patch(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> track_spin = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_patch_photon(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    bool* drift_to_exit /* 0D_NOT_logical in */,
    bool* use_z_pos /* 0D_NOT_logical in */);
void track_a_patch_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> use_z_pos = std::nullopt);
extern "C" void fortran_track_a_pickup(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool* err_flag /* 0D_NOT_logical inout */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_pickup(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_quadrupole(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_quadrupole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_rfcavity(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_rfcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_sad_mult(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_sad_mult(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_sol_quad(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_sol_quad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_thick_multipole(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real inout */,
    bool* make_matrix /* 0D_NOT_logical in */);
void track_a_thick_multipole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_wiggler(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double* mat6 /* 2D_NOT_real out */,
    bool* make_matrix /* 0D_NOT_logical in */);
FixedArray2D<Real, 6, 6> track_a_wiggler(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_zero_length_element(
    void* orbit /* 0D_NOT_type inout */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* track /* 0D_NOT_type out */);
struct TrackAZeroLengthElement {
  bool err_flag;
  TrackProxy track;
};
Bmad::TrackAZeroLengthElement track_a_zero_length_element(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track_all(
    void* lat /* 0D_NOT_type in */,
    void* orbit /* 1D_ALLOC_type inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    int& track_state /* 0D_NOT_integer out */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* orbit0 /* 1D_ALLOC_type out */,
    bool* init_lost /* 0D_NOT_logical in */);
struct TrackAll {
  int track_state;
  bool err_flag;
  CoordProxyAlloc1D orbit0;
};
Bmad::TrackAll track_all(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> init_lost = std::nullopt);
extern "C" void fortran_track_beam(
    void* lat /* 0D_NOT_type in */,
    void* beam /* 0D_NOT_type inout */,
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    void* centroid /* 1D_ALLOC_type in */,
    int* direction /* 0D_NOT_integer in */,
    void* bunch_tracks /* 1D_ALLOC_type inout */);
bool track_beam(
    LatProxy& lat,
    BeamProxy& beam,
    optional_ref<EleProxy> ele1 = std::nullopt,
    optional_ref<EleProxy> ele2 = std::nullopt,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks = std::nullopt);
extern "C" void fortran_track_bunch(
    void* lat /* 0D_NOT_type in */,
    void* bunch /* 0D_NOT_type inout */,
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    void* centroid /* 1D_ALLOC_type in */,
    int* direction /* 0D_NOT_integer in */,
    void* bunch_track /* 0D_NOT_type inout */);
bool track_bunch(
    LatProxy& lat,
    BunchProxy& bunch,
    optional_ref<EleProxy> ele1 = std::nullopt,
    optional_ref<EleProxy> ele2 = std::nullopt,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track_bunch_time(
    void* bunch /* 0D_NOT_type inout */,
    void* branch /* 0D_NOT_type in */,
    double& t_end /* 0D_NOT_real in */,
    double& s_end /* 0D_NOT_real in */,
    void* dt_step /* 1D_ALLOC_real inout */,
    void* extra_field /* 1D_ALLOC_type in */);
void track_bunch_time(
    BunchProxy& bunch,
    BranchProxy& branch,
    double t_end,
    double s_end,
    optional_ref<RealAlloc1D> dt_step = std::nullopt,
    optional_ref<EmFieldProxyAlloc1D> extra_field = std::nullopt);
extern "C" void fortran_track_bunch_to_s(
    void* bunch /* 0D_NOT_type inout */,
    double& s /* 0D_NOT_real in */,
    void* branch /* 0D_NOT_type in */);
void track_bunch_to_s(BunchProxy& bunch, double s, BranchProxy& branch);
extern "C" void fortran_track_bunch_to_t(
    void* bunch /* 0D_NOT_type inout */,
    double& t_target /* 0D_NOT_real in */,
    void* branch /* 0D_NOT_type in */);
void track_bunch_to_t(BunchProxy& bunch, double t_target, BranchProxy& branch);
extern "C" void fortran_track_complex_taylor(
    void* start_orb /* 1D_ALLOC_complex in */,
    void* complex_taylor /* 1D_ALLOC_type in */,
    void* end_orb /* 1D_ALLOC_complex out */);
ComplexAlloc1D track_complex_taylor(
    ComplexAlloc1D& start_orb,
    ComplexTaylorProxyAlloc1D& complex_taylor);
extern "C" void fortran_track_from_s_to_s(
    void* lat /* 0D_NOT_type in */,
    double& s_start /* 0D_NOT_real in */,
    double& s_end /* 0D_NOT_real in */,
    void* orbit_start /* 0D_NOT_type in */,
    void* orbit_end /* 0D_NOT_type out */,
    void* all_orb /* 1D_ALLOC_type out */,
    int* ix_branch /* 0D_NOT_integer in */,
    int& track_state /* 0D_NOT_integer out */,
    int* ix_ele_end /* 0D_NOT_integer in */);
struct TrackFromSToS {
  CoordProxy orbit_end;
  CoordProxyAlloc1D all_orb;
  int track_state;
};
Bmad::TrackFromSToS track_from_s_to_s(
    LatProxy& lat,
    double s_start,
    double s_end,
    CoordProxy& orbit_start,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ix_ele_end = std::nullopt);
extern "C" void fortran_track_many(
    void* lat /* 0D_NOT_type in */,
    void* orbit /* 1D_ALLOC_type inout */,
    int& ix_start /* 0D_NOT_integer in */,
    int& ix_end /* 0D_NOT_integer in */,
    int& direction /* 0D_NOT_integer in */,
    int* ix_branch /* 0D_NOT_integer in */,
    int& track_state /* 0D_NOT_integer out */);
int track_many(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    int ix_start,
    int ix_end,
    int direction,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine track_many_hook_def:
// - Routine in configuration skip list
extern "C" void fortran_track_to_surface(
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type inout */,
    void* param /* 0D_NOT_type in */,
    double* w_surface /* 2D_NOT_real out */);
FixedArray2D<Real, 3, 3> track_to_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    LatParamProxy& param);
extern "C" void fortran_track_until_dead(
    void* start_orb /* 0D_NOT_type in */,
    void* lat /* 0D_NOT_type in */,
    void* end_orb /* 0D_NOT_type out */,
    void* track /* 0D_NOT_type out */);
struct TrackUntilDead {
  CoordProxy end_orb;
  TrackProxy track;
};
Bmad::TrackUntilDead track_until_dead(CoordProxy& start_orb, LatProxy& lat);
extern "C" void fortran_tracking_rad_map_setup(
    void* ele /* 0D_NOT_type in */,
    double& tollerance /* 0D_NOT_real in */,
    int& ref_edge /* 0D_NOT_integer in */,
    void* rad_map /* 0D_NOT_type out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct TrackingRadMapSetup {
  RadMapProxy rad_map;
  bool err_flag;
};
Bmad::TrackingRadMapSetup tracking_rad_map_setup(
    EleProxy& ele,
    double tollerance,
    int ref_edge);
extern "C" void fortran_transfer_ac_kick(
    void* ac_in /* 0D_PTR_type in */,
    void* ac_out /* 0D_PTR_type out */);
AcKickerProxy transfer_ac_kick(AcKickerProxy& ac_in);
extern "C" void fortran_transfer_branch(
    void* branch1 /* 0D_NOT_type in */,
    void* branch2 /* 0D_NOT_type out */);
BranchProxy transfer_branch(BranchProxy& branch1);
extern "C" void fortran_transfer_branch_parameters(
    void* branch_in /* 0D_NOT_type in */,
    void* branch_out /* 0D_NOT_type out */);
BranchProxy transfer_branch_parameters(BranchProxy& branch_in);
extern "C" void fortran_transfer_branches(
    void* branch1 /* 1D_ALLOC_type in */,
    void* branch2 /* 1D_ALLOC_type out */);
BranchProxyAlloc1D transfer_branches(BranchProxyAlloc1D& branch1);
extern "C" void fortran_transfer_ele(
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type out */,
    bool* nullify_pointers /* 0D_NOT_logical in */);
EleProxy transfer_ele(
    EleProxy& ele1,
    std::optional<bool> nullify_pointers = std::nullopt);
extern "C" void fortran_transfer_ele_taylor(
    void* ele_in /* 0D_NOT_type in */,
    void* ele_out /* 0D_NOT_type out */,
    int* taylor_order /* 0D_NOT_integer in */);
EleProxy transfer_ele_taylor(
    EleProxy& ele_in,
    std::optional<int> taylor_order = std::nullopt);
extern "C" void fortran_transfer_eles(
    void* ele1 /* 1D_ALLOC_type in */,
    void* ele2 /* 1D_ALLOC_type out */);
EleProxyAlloc1D transfer_eles(EleProxyAlloc1D& ele1);
extern "C" void fortran_transfer_fieldmap(
    void* ele_in /* 0D_NOT_type in */,
    void* ele_out /* 0D_NOT_type out */,
    int& who /* 0D_NOT_integer in */);
EleProxy transfer_fieldmap(EleProxy& ele_in, int who);
extern "C" bool fortran_transfer_fixer_params(
    void* fixer /* 0D_NOT_type in */,
    bool& to_stored /* 0D_NOT_logical in */,
    void* orbit /* 0D_NOT_type in */,
    const char* who /* 0D_NOT_character in */,
    bool& is_ok /* 0D_NOT_logical out */);
bool transfer_fixer_params(
    EleProxy& fixer,
    bool to_stored,
    optional_ref<CoordProxy> orbit = std::nullopt,
    std::optional<std::string> who = std::nullopt);
extern "C" void fortran_transfer_lat(
    void* lat1 /* 0D_NOT_type in */,
    void* lat2 /* 0D_NOT_type out */);
LatProxy transfer_lat(LatProxy& lat1);
extern "C" void fortran_transfer_lat_parameters(
    void* lat_in /* 0D_NOT_type in */,
    void* lat_out /* 0D_NOT_type out */);
LatProxy transfer_lat_parameters(LatProxy& lat_in);
extern "C" void fortran_transfer_map_calc(
    void* lat /* 0D_NOT_type in */,
    void* orb_map /* 1D_ALLOC_type inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    int* ix1 /* 0D_NOT_integer in */,
    int* ix2 /* 0D_NOT_integer in */,
    void* ref_orb /* 0D_NOT_type in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* one_turn /* 0D_NOT_logical in */,
    bool* unit_start /* 0D_NOT_logical in */,
    bool* concat_if_possible /* 0D_NOT_logical in */,
    void* spin_map /* 1D_ALLOC_type inout */);
bool transfer_map_calc(
    LatProxy& lat,
    TaylorProxyAlloc1D& orb_map,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    optional_ref<CoordProxy> ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt,
    std::optional<bool> concat_if_possible = std::nullopt,
    optional_ref<TaylorProxyAlloc1D> spin_map = std::nullopt);
extern "C" void fortran_transfer_map_from_s_to_s(
    void* lat /* 0D_NOT_type in */,
    void* t_map /* 1D_ALLOC_type inout */,
    double* s1 /* 0D_NOT_real in */,
    double* s2 /* 0D_NOT_real in */,
    void* ref_orb_in /* 0D_NOT_type in */,
    void* ref_orb_out /* 0D_NOT_type out */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* one_turn /* 0D_NOT_logical in */,
    bool* unit_start /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* concat_if_possible /* 0D_NOT_logical in */,
    void* spin_map /* 1D_ALLOC_type inout */);
struct TransferMapFromSToS {
  CoordProxy ref_orb_out;
  bool err_flag;
};
Bmad::TransferMapFromSToS transfer_map_from_s_to_s(
    LatProxy& lat,
    TaylorProxyAlloc1D& t_map,
    std::optional<double> s1 = std::nullopt,
    std::optional<double> s2 = std::nullopt,
    optional_ref<CoordProxy> ref_orb_in = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt,
    std::optional<bool> unit_start = std::nullopt,
    std::optional<bool> concat_if_possible = std::nullopt,
    optional_ref<TaylorProxyAlloc1D> spin_map = std::nullopt);
extern "C" void fortran_transfer_mat2_from_twiss(
    void* twiss1 /* 0D_NOT_type in */,
    void* twiss2 /* 0D_NOT_type in */,
    double* mat /* 2D_NOT_real out */);
FixedArray2D<Real, 2, 2> transfer_mat2_from_twiss(
    TwissProxy& twiss1,
    TwissProxy& twiss2);
extern "C" void fortran_transfer_mat_from_twiss(
    void* ele1 /* 0D_NOT_type in */,
    void* ele2 /* 0D_NOT_type in */,
    double* orb1 /* 1D_NOT_real in */,
    double* orb2 /* 1D_NOT_real in */,
    double* m /* 2D_NOT_real out */);
FixedArray2D<Real, 6, 6> transfer_mat_from_twiss(
    EleProxy& ele1,
    EleProxy& ele2,
    FixedArray1D<Real, 6> orb1,
    FixedArray1D<Real, 6> orb2);
extern "C" void fortran_transfer_matrix_calc(
    void* lat /* 0D_NOT_type in */,
    double* xfer_mat /* 2D_NOT_real inout */,
    double* xfer_vec /* 1D_NOT_real inout */,
    int* ix1 /* 0D_NOT_integer in */,
    int* ix2 /* 0D_NOT_integer in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* one_turn /* 0D_NOT_logical in */);
void transfer_matrix_calc(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> xfer_mat,
    std::optional<FixedArray1D<Real, 6>> xfer_vec = std::nullopt,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt);
extern "C" void fortran_transfer_twiss(
    void* ele_in /* 0D_NOT_type in */,
    void* ele_out /* 0D_NOT_type out */,
    bool* reverse /* 0D_NOT_logical in */);
EleProxy transfer_twiss(
    EleProxy& ele_in,
    std::optional<bool> reverse = std::nullopt);
extern "C" void fortran_transfer_wake(
    void* wake_in /* 0D_PTR_type in */,
    void* wake_out /* 0D_PTR_type out */);
WakeProxy transfer_wake(WakeProxy& wake_in);

// Skipped unusable routine transfer_wall3d:
// - Routine in configuration skip list

// Skipped unusable routine transport_with_sr:
// - Variable in sized array: M(:,:) 2D_NOT_real
// - Variable in sized array: Bone(:,:) 2D_NOT_real
// - Variable inout sized array: Yone(:,:) 2D_NOT_real

// Skipped unusable routine transport_with_sr_and_ibs:
// - Variable in sized array: M(:,:) 2D_NOT_real
// - Variable in sized array: Bone(:,:) 2D_NOT_real
// - Variable in sized array: Yone(:,:) 2D_NOT_real
extern "C" void fortran_truncate_complex_taylor_to_order(
    void* complex_taylor_in /* 1D_ALLOC_type in */,
    int& order /* 0D_NOT_integer in */,
    void* complex_taylor_out /* 1D_ALLOC_type out */);
ComplexTaylorProxyAlloc1D truncate_complex_taylor_to_order(
    ComplexTaylorProxyAlloc1D& complex_taylor_in,
    int order);
extern "C" void fortran_twiss1_propagate(
    void* twiss1 /* 0D_NOT_type in */,
    double* mat2 /* 2D_NOT_real in */,
    int& ele_key /* 0D_NOT_integer in */,
    double& length /* 0D_NOT_real in */,
    void* twiss2 /* 0D_NOT_type out */,
    bool& err /* 0D_NOT_logical out */);
struct Twiss1Propagate {
  TwissProxy twiss2;
  bool err;
};
Bmad::Twiss1Propagate twiss1_propagate(
    TwissProxy& twiss1,
    FixedArray2D<Real, 2, 2> mat2,
    int ele_key,
    double length);
extern "C" void fortran_twiss3_at_start(
    void* lat /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    double* tune3 /* 1D_NOT_real out */);
FixedArray1D<Real, 3> twiss3_at_start(
    LatProxy& lat,
    bool& err_flag,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss3_from_twiss2(void* ele /* 0D_NOT_type inout */);
void twiss3_from_twiss2(EleProxy& ele);
extern "C" void fortran_twiss3_propagate1(
    void* ele1 /* 0D_NOT_type inout */,
    void* ele2 /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void twiss3_propagate1(EleProxy& ele1, EleProxy& ele2, bool& err_flag);
extern "C" void fortran_twiss3_propagate_all(
    void* lat /* 0D_NOT_type in */,
    int* ix_branch /* 0D_NOT_integer in */);
void twiss3_propagate_all(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss_and_track_at_s(
    void* lat /* 0D_NOT_type in */,
    double& s /* 0D_NOT_real in */,
    void* ele_at_s /* 0D_NOT_type inout */,
    void* orb /* 1D_ALLOC_type in */,
    void* orb_at_s /* 0D_NOT_type inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err /* 0D_NOT_logical out */,
    bool* use_last /* 0D_NOT_logical in */,
    bool* compute_floor_coords /* 0D_NOT_logical in */);
bool twiss_and_track_at_s(
    LatProxy& lat,
    double s,
    optional_ref<EleProxy> ele_at_s = std::nullopt,
    optional_ref<CoordProxyAlloc1D> orb = std::nullopt,
    optional_ref<CoordProxy> orb_at_s = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> use_last = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt);
extern "C" void fortran_twiss_and_track_from_s_to_s(
    void* branch /* 0D_NOT_type in */,
    void* orbit_start /* 0D_NOT_type in */,
    double& s_end /* 0D_NOT_real in */,
    void* orbit_end /* 0D_NOT_type out */,
    void* ele_start /* 0D_NOT_type in */,
    void* ele_end /* 0D_NOT_type out */,
    bool& err /* 0D_NOT_logical out */,
    bool* compute_floor_coords /* 0D_NOT_logical in */,
    bool* compute_twiss /* 0D_NOT_logical in */);
struct TwissAndTrackFromSToS {
  CoordProxy orbit_end;
  EleProxy ele_end;
  bool err;
};
Bmad::TwissAndTrackFromSToS twiss_and_track_from_s_to_s(
    BranchProxy& branch,
    CoordProxy& orbit_start,
    double s_end,
    optional_ref<EleProxy> ele_start = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt,
    std::optional<bool> compute_twiss = std::nullopt);
extern "C" void fortran_twiss_and_track_intra_ele(
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& l_start /* 0D_NOT_real in */,
    double& l_end /* 0D_NOT_real in */,
    bool& track_upstream_end /* 0D_NOT_logical in */,
    bool& track_downstream_end /* 0D_NOT_logical in */,
    void* orbit_start /* 0D_NOT_type in */,
    void* orbit_end /* 0D_NOT_type out */,
    void* ele_start /* 0D_NOT_type in */,
    void* ele_end /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical out */,
    bool* compute_floor_coords /* 0D_NOT_logical in */,
    bool* compute_twiss /* 0D_NOT_logical in */,
    bool* reuse_ele_end /* 0D_NOT_logical in */);
struct TwissAndTrackIntraEle {
  CoordProxy orbit_end;
  bool err;
};
Bmad::TwissAndTrackIntraEle twiss_and_track_intra_ele(
    EleProxy& ele,
    LatParamProxy& param,
    double l_start,
    double l_end,
    bool track_upstream_end,
    bool track_downstream_end,
    optional_ref<CoordProxy> orbit_start = std::nullopt,
    optional_ref<EleProxy> ele_start = std::nullopt,
    optional_ref<EleProxy> ele_end = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt,
    std::optional<bool> compute_twiss = std::nullopt,
    std::optional<bool> reuse_ele_end = std::nullopt);
extern "C" void fortran_twiss_at_element(
    void* ele /* 0D_NOT_type in */,
    void* start /* 0D_NOT_type out */,
    void* end /* 0D_NOT_type out */,
    void* average /* 0D_NOT_type out */);
struct TwissAtElement {
  EleProxy start;
  EleProxy end;
  EleProxy average;
};
Bmad::TwissAtElement twiss_at_element(EleProxy& ele);
extern "C" void fortran_twiss_at_start(
    void* lat /* 0D_NOT_type inout */,
    int& status /* 0D_NOT_integer out */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* type_out /* 0D_NOT_logical in */);
int twiss_at_start(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> type_out = std::nullopt);

// Skipped unusable routine twiss_from_mat2:
// - Variable inout sized array: mat_in(:,:) 2D_NOT_real

// Skipped unusable routine twiss_from_mat6:
// - Variable in sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_twiss_from_tracking(
    void* lat /* 0D_NOT_type inout */,
    void* ref_orb0 /* 0D_NOT_type in */,
    double& symp_err /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    void* d_orb /* 1D_ALLOC_real in */);
struct TwissFromTracking {
  double symp_err;
  bool err_flag;
};
Bmad::TwissFromTracking twiss_from_tracking(
    LatProxy& lat,
    CoordProxy& ref_orb0,
    optional_ref<RealAlloc1D> d_orb = std::nullopt);
extern "C" void fortran_twiss_propagate1(
    void* ele1 /* 0D_NOT_type inout */,
    void* ele2 /* 0D_NOT_type inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* forward /* 0D_NOT_logical in */);
bool twiss_propagate1(
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> forward = std::nullopt);
extern "C" void fortran_twiss_propagate_all(
    void* lat /* 0D_NOT_type inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err_flag /* 0D_NOT_logical out */,
    int* ie_start /* 0D_NOT_integer in */,
    int* ie_end /* 0D_NOT_integer in */);
bool twiss_propagate_all(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ie_start = std::nullopt,
    std::optional<int> ie_end = std::nullopt);
extern "C" void fortran_twiss_to_1_turn_mat(
    void* twiss /* 0D_NOT_type in */,
    double& phi /* 0D_NOT_real in */,
    double* mat2 /* 2D_NOT_real out */);
FixedArray2D<Real, 2, 2> twiss_to_1_turn_mat(TwissProxy& twiss, double phi);

// Skipped unusable routine type_complex_taylors:
// - Variable-sized out character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_type_coord(void* coord /* 0D_NOT_type in */);
void type_coord(CoordProxy& coord);

// Skipped unusable routine type_ele:
// - Variable-sized in character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_end_stuff:
// - Variable-sized inout character array: li(:) 1D_ALLOC_character
// - Variable-sized inout character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_type_expression_tree(
    void* tree /* 0D_NOT_type in */,
    int* indent /* 0D_NOT_integer in */);
void type_expression_tree(
    ExpressionTreeProxy& tree,
    std::optional<int> indent = std::nullopt);

// Skipped unusable routine type_map:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_map1:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_ptc_fibre:
// - Untranslated type: fibre (0D)
// - Variable-sized out character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_ptc_internal_state:
// - Untranslated type: internal_state (0D)
// - Variable-sized out character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_ptc_layout:
// - Untranslated type: layout (0D)

// Skipped unusable routine type_real_8_taylors:
// - Untranslated type: real_8 (1D)

// Skipped unusable routine type_taylors:
// - Variable-sized inout character array: lines(:) 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_twiss:
// - Variable-sized out character array: lines(:) 1D_ALLOC_character
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
extern "C" void fortran_update_ele_from_fibre(
    void* ele /* 0D_NOT_type inout */);
void update_ele_from_fibre(EleProxy& ele);
extern "C" void fortran_update_fibre_from_ele(
    void* ele /* 0D_NOT_type in */,
    bool& survey_needed /* 0D_NOT_logical out */);
bool update_fibre_from_ele(EleProxy& ele);
extern "C" void fortran_update_floor_angles(
    void* floor /* 0D_NOT_type inout */,
    void* floor0 /* 0D_NOT_type in */);
void update_floor_angles(
    FloorPositionProxy& floor,
    optional_ref<FloorPositionProxy> floor0 = std::nullopt);
extern "C" bool fortran_valid_field_calc(
    void* ele /* 0D_NOT_type in */,
    int& field_calc /* 0D_NOT_integer in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void valid_field_calc(EleProxy& ele, int field_calc, bool& is_valid);
extern "C" bool fortran_valid_fringe_type(
    void* ele /* 0D_NOT_type in */,
    int& fringe_type /* 0D_NOT_integer in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void valid_fringe_type(EleProxy& ele, int fringe_type, bool& is_valid);
extern "C" bool fortran_valid_mat6_calc_method(
    void* ele /* 0D_NOT_type in */,
    int& species /* 0D_NOT_integer in */,
    int& mat6_calc_method /* 0D_NOT_integer in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void valid_mat6_calc_method(
    EleProxy& ele,
    int species,
    int mat6_calc_method,
    bool& is_valid);
extern "C" bool fortran_valid_spin_tracking_method(
    void* ele /* 0D_NOT_type in */,
    int& spin_tracking_method /* 0D_NOT_integer in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void valid_spin_tracking_method(
    EleProxy& ele,
    int spin_tracking_method,
    bool& is_valid);
extern "C" bool fortran_valid_tracking_method(
    void* ele /* 0D_NOT_type in */,
    int& species /* 0D_NOT_integer in */,
    int& tracking_method /* 0D_NOT_integer in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void valid_tracking_method(
    EleProxy& ele,
    int species,
    int tracking_method,
    bool& is_valid);
extern "C" bool fortran_value_of_attribute(
    void* ele /* 0D_NOT_type in */,
    const char* attrib_name /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* err_print_flag /* 0D_NOT_logical in */,
    double* err_value /* 0D_NOT_real in */,
    double& value /* 0D_NOT_real inout */);
bool value_of_attribute(
    EleProxy& ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<double> err_value,
    double& value);
extern "C" void fortran_value_to_line(
    const char* line /* 0D_NOT_character inout */,
    double& value /* 0D_NOT_real inout */,
    const char* str /* 0D_NOT_character inout */,
    const char* typ /* 0D_NOT_character inout */,
    bool* ignore_if_zero /* 0D_NOT_logical inout */,
    bool* use_comma /* 0D_NOT_logical inout */);
void value_to_line(
    std::string& line,
    double& value,
    std::string& str,
    std::string& typ,
    optional_ref<bool> ignore_if_zero = std::nullopt,
    optional_ref<bool> use_comma = std::nullopt);
extern "C" bool fortran_vec_to_polar(
    double* vec /* 1D_NOT_real in */,
    double* phase /* 0D_NOT_real in */,
    void* polar /* 0D_NOT_type inout */);
void vec_to_polar(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    SpinPolarProxy& polar);
extern "C" bool fortran_vec_to_spinor(
    double* vec /* 1D_NOT_real in */,
    double* phase /* 0D_NOT_real in */,
    std::complex<double>* spinor /* 1D_NOT_complex inout */);
void vec_to_spinor(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    FixedArray1D<Complex, 2> spinor);
extern "C" bool fortran_verify_valid_name(
    const char* name /* 0D_NOT_character in */,
    int& ix_name /* 0D_NOT_integer in */,
    bool* pure_name /* 0D_NOT_logical in */,
    bool* include_wild /* 0D_NOT_logical in */,
    bool& is_valid /* 0D_NOT_logical out */);
bool verify_valid_name(
    std::string name,
    int ix_name,
    std::optional<bool> pure_name = std::nullopt,
    std::optional<bool> include_wild = std::nullopt);
extern "C" bool fortran_w_mat_for_bend_angle(
    double& angle /* 0D_NOT_real in */,
    double& ref_tilt /* 0D_NOT_real in */,
    double* r_vec /* 1D_NOT_real inout */,
    double* w_mat /* 2D_NOT_real inout */);
void w_mat_for_bend_angle(
    double angle,
    double ref_tilt,
    std::optional<FixedArray1D<Real, 3>> r_vec,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_tilt(
    double& tilt /* 0D_NOT_real in */,
    bool* return_inverse /* 0D_NOT_logical in */,
    double* w_mat /* 2D_NOT_real inout */);
void w_mat_for_tilt(
    double tilt,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_x_pitch(
    double& x_pitch /* 0D_NOT_real in */,
    bool* return_inverse /* 0D_NOT_logical in */,
    double* w_mat /* 2D_NOT_real inout */);
void w_mat_for_x_pitch(
    double x_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_y_pitch(
    double& y_pitch /* 0D_NOT_real in */,
    bool* return_inverse /* 0D_NOT_logical in */,
    double* w_mat /* 2D_NOT_real inout */);
void w_mat_for_y_pitch(
    double y_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_wall3d_d_radius(
    void* position /* 1D_ALLOC_real in */,
    void* ele /* 0D_NOT_type in */,
    int* ix_wall /* 0D_NOT_integer in */,
    double* perp /* 1D_NOT_real out */,
    int& ix_section /* 0D_NOT_integer out */,
    bool& no_wall_here /* 0D_NOT_logical out */,
    double* origin /* 1D_NOT_real out */,
    double& radius_wall /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */,
    double& d_radius /* 0D_NOT_real out */);
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
    RealAlloc1D& position,
    EleProxy& ele,
    std::optional<int> ix_wall = std::nullopt);
extern "C" void fortran_wall3d_initializer(
    void* wall3d /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical out */);
bool wall3d_initializer(Wall3dProxy& wall3d);
extern "C" void fortran_wall3d_section_initializer(
    void* section /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical out */);
bool wall3d_section_initializer(Wall3dSectionProxy& section);
extern "C" bool fortran_wall3d_to_position(
    void* orbit /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double* position /* 1D_NOT_real out */);
FixedArray1D<Real, 6> wall3d_to_position(CoordProxy& orbit, EleProxy& ele);

// Skipped unusable routine wall_hit_handler_custom_def:
// - Routine in configuration skip list
extern "C" void fortran_word_to_value(
    const char* word /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    double& value /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical inout */,
    void* ele /* 0D_NOT_type inout */);
void word_to_value(
    std::string& word,
    LatProxy& lat,
    double& value,
    bool& err_flag,
    optional_ref<EleProxy> ele = std::nullopt);

// Skipped unusable routine write_2d:
// - Variable inout sized array: grid(:,:) 2D_NOT_real
extern "C" void fortran_write_ascii_beam_file(
    const char* file_name /* 0D_NOT_character in */,
    void* beam /* 0D_NOT_type in */,
    bool* new_file /* 0D_NOT_logical in */,
    bool* alive_only /* 0D_NOT_logical in */);
void write_ascii_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<bool> alive_only = std::nullopt);
extern "C" void fortran_write_astra_bend(
    int& iu /* 0D_NOT_integer inout */,
    double& strength /* 0D_NOT_real inout */,
    int& id /* 0D_NOT_integer inout */,
    double* d1 /* 1D_NOT_real inout */,
    double* d2 /* 1D_NOT_real inout */,
    double* d3 /* 1D_NOT_real inout */,
    double* d4 /* 1D_NOT_real inout */);
void write_astra_bend(
    int& iu,
    double& strength,
    int& id,
    FixedArray1D<Real, 2> d1,
    FixedArray1D<Real, 2> d2,
    FixedArray1D<Real, 2> d3,
    FixedArray1D<Real, 2> d4);

// Skipped unusable routine write_astra_ele:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_write_astra_field_grid_file(
    int& astra_file_unit /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    double* dz /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */);
struct WriteAstraFieldGridFile {
  double maxfield;
  bool err;
};
Bmad::WriteAstraFieldGridFile write_astra_field_grid_file(
    int astra_file_unit,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);
extern "C" void fortran_write_astra_field_grid_file_3d(
    const char* base_filename /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    double* dz /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */);
struct WriteAstraFieldGridFile3d {
  double maxfield;
  bool err;
};
Bmad::WriteAstraFieldGridFile3d write_astra_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);

// Skipped unusable routine write_astra_lattice_file:
// - Untranslated type: astra_lattice_param_struct (0D)
extern "C" void fortran_write_beam_file(
    const char* file_name /* 0D_NOT_character in */,
    void* beam /* 0D_NOT_type in */,
    bool* new_file /* 0D_NOT_logical in */,
    int* file_format /* 0D_NOT_integer in */,
    void* lat /* 0D_NOT_type in */,
    bool* alive_only /* 0D_NOT_logical in */);
void write_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<int> file_format = std::nullopt,
    optional_ref<LatProxy> lat = std::nullopt,
    std::optional<bool> alive_only = std::nullopt);
extern "C" void fortran_write_beam_floor_positions(
    const char* file_name /* 0D_NOT_character in */,
    void* beam /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    bool* new_file /* 0D_NOT_logical in */);
void write_beam_floor_positions(
    std::string file_name,
    BeamProxy& beam,
    EleProxy& ele,
    std::optional<bool> new_file = std::nullopt);
extern "C" void fortran_write_binary_cartesian_map(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* cart_map /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void write_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag);
extern "C" void fortran_write_binary_cylindrical_map(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* cl_map /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void write_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag);
extern "C" void fortran_write_binary_grid_field(
    const char* file_name /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    void* g_field /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical in */);
void write_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag);
extern "C" void fortran_write_blender_ele(
    int& iu /* 0D_NOT_integer inout */,
    void* ele /* 0D_NOT_type inout */,
    bool* old_format /* 0D_NOT_logical inout */);
void write_blender_ele(
    int& iu,
    EleProxy& ele,
    optional_ref<bool> old_format = std::nullopt);
extern "C" void fortran_write_blender_lat_layout(
    const char* file_name /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */);
void write_blender_lat_layout(std::string& file_name, LatProxy& lat);
extern "C" void fortran_write_bmad_lattice_file(
    const char* bmad_file /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */,
    int* output_form /* 0D_NOT_integer in */,
    void* orbit0 /* 0D_NOT_type in */);
bool write_bmad_lattice_file(
    std::string bmad_file,
    LatProxy& lat,
    std::optional<int> output_form = std::nullopt,
    optional_ref<CoordProxy> orbit0 = std::nullopt);

// Skipped unusable routine write_digested_bmad_file:
// - Variable-sized in character array: file_names(:) 1D_ALLOC_character
// - Untranslated type: extra_parsing_info_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine write_gpt_ele:
// - Untranslated type: str_index_struct (0D)
extern "C" void fortran_write_gpt_field_grid_file_1d(
    int& gpt_file_unit /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    double& ref_time /* 0D_NOT_real out */,
    double* dz /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */);
struct WriteGptFieldGridFile1d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile1d write_gpt_field_grid_file_1d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);
extern "C" void fortran_write_gpt_field_grid_file_2d(
    int& gpt_file_unit /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    double& ref_time /* 0D_NOT_real out */,
    double* dr /* 0D_NOT_real in */,
    double* dz /* 0D_NOT_real in */,
    double* r_max /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */);
struct WriteGptFieldGridFile2d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile2d write_gpt_field_grid_file_2d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dr = std::nullopt,
    std::optional<double> dz = std::nullopt,
    std::optional<double> r_max = std::nullopt);
extern "C" void fortran_write_gpt_field_grid_file_3d(
    const char* base_filename /* 0D_NOT_character in */,
    void* ele /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    double& ref_time /* 0D_NOT_real out */,
    double* dz /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */);
struct WriteGptFieldGridFile3d {
  double maxfield;
  double ref_time;
  bool err;
};
Bmad::WriteGptFieldGridFile3d write_gpt_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);

// Skipped unusable routine write_gpt_lattice_file:
// - Untranslated type: gpt_lat_param_struct (0D)
extern "C" void fortran_write_lat_line(
    const char* line /* 0D_NOT_character inout */,
    int& iu /* 0D_NOT_integer in */,
    bool& end_is_neigh /* 0D_NOT_logical in */,
    bool* do_split /* 0D_NOT_logical in */,
    bool* scibmad /* 0D_NOT_logical in */);
void write_lat_line(
    std::string& line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split = std::nullopt,
    std::optional<bool> scibmad = std::nullopt);
extern "C" void fortran_write_lattice_in_elegant_format(
    const char* out_file_name /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    void* ref_orbit /* 1D_ALLOC_type in */,
    bool* use_matrix_model /* 0D_NOT_logical in */,
    bool* include_apertures /* 0D_NOT_logical in */,
    double* dr12_drift_max /* 0D_NOT_real in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err /* 0D_NOT_logical out */);
bool write_lattice_in_elegant_format(
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_write_lattice_in_foreign_format(
    const char* out_type /* 0D_NOT_character in */,
    const char* out_file_name /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    void* ref_orbit /* 1D_ALLOC_type in */,
    bool* use_matrix_model /* 0D_NOT_logical in */,
    bool* include_apertures /* 0D_NOT_logical in */,
    double* dr12_drift_max /* 0D_NOT_real in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err /* 0D_NOT_logical out */);
bool write_lattice_in_foreign_format(
    std::string out_type,
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_write_lattice_in_mad_format(
    const char* out_type /* 0D_NOT_character in */,
    const char* out_file_name /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    void* ref_orbit /* 1D_ALLOC_type in */,
    bool* use_matrix_model /* 0D_NOT_logical in */,
    bool* include_apertures /* 0D_NOT_logical in */,
    double* dr12_drift_max /* 0D_NOT_real in */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool& err /* 0D_NOT_logical out */);
bool write_lattice_in_mad_format(
    std::string out_type,
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_write_lattice_in_sad_format(
    const char* out_file_name /* 0D_NOT_character inout */,
    void* lat /* 0D_NOT_type inout */,
    bool* include_apertures /* 0D_NOT_logical inout */,
    int* ix_branch /* 0D_NOT_integer inout */,
    bool* err /* 0D_NOT_logical inout */);
void write_lattice_in_sad_format(
    std::string& out_file_name,
    LatProxy& lat,
    optional_ref<bool> include_apertures = std::nullopt,
    optional_ref<int> ix_branch = std::nullopt,
    optional_ref<bool> err = std::nullopt);
extern "C" void fortran_write_lattice_in_scibmad(
    const char* scibmad_file /* 0D_NOT_character out */,
    void* lat /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
struct WriteLatticeInScibmad {
  std::string scibmad_file;
  bool err_flag;
};
Bmad::WriteLatticeInScibmad write_lattice_in_scibmad(LatProxy& lat);
extern "C" void fortran_write_line_element(
    const char* line /* 0D_NOT_character inout */,
    int& iu /* 0D_NOT_integer inout */,
    void* ele /* 0D_NOT_type inout */,
    void* lat /* 0D_NOT_type inout */);
void write_line_element(
    std::string& line,
    int& iu,
    EleProxy& ele,
    LatProxy& lat);
extern "C" void fortran_write_opal_field_grid_file(
    int& opal_file_unit /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* param /* 0D_NOT_type in */,
    double& maxfield /* 0D_NOT_real out */,
    bool& err /* 0D_NOT_logical out */);
struct WriteOpalFieldGridFile {
  double maxfield;
  bool err;
};
Bmad::WriteOpalFieldGridFile write_opal_field_grid_file(
    int opal_file_unit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_write_opal_lattice_file(
    int& opal_file_unit /* 0D_NOT_integer in */,
    void* lat /* 0D_NOT_type in */,
    bool& err /* 0D_NOT_logical out */);
bool write_opal_lattice_file(int opal_file_unit, LatProxy& lat);
extern "C" void fortran_write_time_particle_distribution(
    int& time_file_unit /* 0D_NOT_integer in */,
    void* bunch /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    const char* style /* 0D_NOT_character in */,
    void* branch /* 0D_NOT_type in */,
    const char* format /* 0D_NOT_character in */,
    bool& err /* 0D_NOT_logical out */);
bool write_time_particle_distribution(
    int time_file_unit,
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<std::string> style = std::nullopt,
    optional_ref<BranchProxy> branch = std::nullopt,
    std::optional<std::string> format = std::nullopt);
extern "C" bool fortran_xlafun(
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void xlafun(double& x, double& y, double& z, double& res);
extern "C" bool fortran_xraylib_nist_compound(
    const char* name /* 0D_NOT_character in */,
    int& indx /* 0D_NOT_integer out */);
int xraylib_nist_compound(std::string name);

// Skipped unusable routine xsif_parser:
// - Routine in configuration skip list

// Skipped unusable routine xyz_to_action:
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_ylafun(
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void ylafun(double& x, double& y, double& z, double& res);
extern "C" bool fortran_z_at_surface(
    void* ele /* 0D_NOT_type in */,
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    bool& err_flag /* 0D_NOT_logical out */,
    bool* extend_grid /* 0D_NOT_logical in */,
    double* dz_dxy /* 1D_NOT_real out */,
    double& z /* 0D_NOT_real out */);
struct ZAtSurface {
  bool err_flag;
  FixedArray1D<Real, 2> dz_dxy;
  double z;
};
Bmad::ZAtSurface z_at_surface(
    EleProxy& ele,
    double& x,
    double& y,
    std::optional<bool> extend_grid = std::nullopt);
extern "C" void fortran_zero_ele_kicks(void* ele /* 0D_NOT_type out */);
EleProxy zero_ele_kicks();
extern "C" void fortran_zero_ele_offsets(void* ele /* 0D_NOT_type out */);
EleProxy zero_ele_offsets();
extern "C" void fortran_zero_lr_wakes_in_lat(void* lat /* 0D_NOT_type inout */);
void zero_lr_wakes_in_lat(LatProxy& lat);
extern "C" bool fortran_zlafun(
    double& x /* 0D_NOT_real inout */,
    double& y /* 0D_NOT_real inout */,
    double& z /* 0D_NOT_real inout */,
    double& res /* 0D_NOT_real inout */);
void zlafun(double& x, double& y, double& z, double& res);

// Skipped unusable routine zot_integrand:
// - Untranslated type: c_ptr (0D)
} // namespace Bmad
