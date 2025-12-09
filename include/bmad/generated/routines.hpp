#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.h"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace Bmad {
extern "C" void fortran_ab_multipole_kick(
    c_Real& a /* 0D_NOT_real */,
    c_Real& b /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Int& ref_species /* 0D_NOT_integer */,
    c_Int& ele_orientation /* 0D_NOT_integer */,
    void* coord /* 0D_NOT_type */,
    c_Real& kx /* 0D_NOT_real */,
    c_Real& ky /* 0D_NOT_real */,
    c_RealArr dk /* 2D_NOT_real */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Real* scale /* 0D_NOT_real */);
struct AbMultipoleKick {
  double kx;
  double ky;
  std::optional<FixedArray2D<Real, 2, 2>> dk;
};
AbMultipoleKick ab_multipole_kick(
    double a,
    double b,
    int n,
    int ref_species,
    int ele_orientation,
    CoordProxy& coord,
    std::optional<int> pole_type = std::nullopt,
    std::optional<double> scale = std::nullopt);
extern "C" void fortran_ab_multipole_kicks(
    void* an /* 1D_ALLOC_real */,
    void* bn /* 1D_ALLOC_real */,
    c_Int& ix_pole_max /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Real* scale /* 0D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
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
    void* e_orb /* 0D_NOT_type */,
    void* photon_orb /* 0D_NOT_type */);
void absolute_photon_position(CoordProxy& e_orb, CoordProxy& photon_orb);
extern "C" bool fortran_absolute_time_tracking(
    void* ele /* 0D_NOT_type */,
    c_Bool& is_abs_time /* 0D_NOT_logical */);
void absolute_time_tracking(EleProxy& ele, bool is_abs_time);
extern "C" bool fortran_ac_kicker_amp(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Real* true_time /* 0D_NOT_real */,
    c_Real& ac_amp /* 0D_NOT_real */);
void ac_kicker_amp(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<double> true_time,
    double ac_amp);

// Skipped unusable routine ac_kicker_freq_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ac_kicker_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ac_kicker_time_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_action_to_xyz(
    void* ring /* 0D_NOT_type */,
    c_Int& ix /* 0D_NOT_integer */,
    c_RealArr J /* 1D_NOT_real */,
    c_RealArr X /* 1D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ActionToXyz {
  FixedArray1D<Real, 6> X;
  bool err_flag;
};
ActionToXyz action_to_xyz(LatProxy& ring, int ix, FixedArray1D<Real, 6> J);
extern "C" void fortran_add_lattice_control_structs(
    void* ele /* 0D_NOT_type */,
    c_Int* n_add_slave /* 0D_NOT_integer */,
    c_Int* n_add_lord /* 0D_NOT_integer */,
    c_Int* n_add_slave_field /* 0D_NOT_integer */,
    c_Int* n_add_lord_field /* 0D_NOT_integer */,
    c_Bool* add_at_end /* 0D_NOT_logical */);
void add_lattice_control_structs(
    EleProxy& ele,
    std::optional<int> n_add_slave = std::nullopt,
    std::optional<int> n_add_lord = std::nullopt,
    std::optional<int> n_add_slave_field = std::nullopt,
    std::optional<int> n_add_lord_field = std::nullopt,
    std::optional<bool> add_at_end = std::nullopt);

// Skipped unusable routine add_ptc_layout_to_list:
// Untranslated type: PtcBranch1Proxy (0D_NOT_type)
// Untranslated type: LayoutProxy (0D_NOT_type)
extern "C" void fortran_add_superimpose(
    void* lat /* 0D_NOT_type */,
    void* super_ele_in /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* super_ele_out /* 0D_PTR_type */,
    c_Bool* save_null_drift /* 0D_NOT_logical */,
    c_Bool* create_jumbo_slave /* 0D_NOT_logical */,
    c_Int* ix_insert /* 0D_NOT_integer */,
    c_Bool* mangle_slave_names /* 0D_NOT_logical */,
    c_Bool* wrap /* 0D_NOT_logical */);
struct AddSuperimpose {
  bool err_flag;
  EleProxy super_ele_out;
};
AddSuperimpose add_superimpose(
    LatProxy& lat,
    EleProxy& super_ele_in,
    int ix_branch,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> create_jumbo_slave = std::nullopt,
    std::optional<int> ix_insert = std::nullopt,
    std::optional<bool> mangle_slave_names = std::nullopt,
    std::optional<bool> wrap = std::nullopt);
extern "C" void fortran_add_this_multipass(
    void* lat /* 0D_NOT_type */,
    void* m_slaves /* 1D_ALLOC_type */,
    void* lord_in /* 0D_NOT_type */);
void add_this_multipass(
    LatProxy& lat,
    LatEleLocProxyAlloc1D& m_slaves,
    optional_ref<EleProxy> lord_in = std::nullopt);

// Skipped unusable routine add_this_name_to_list:
// Variable-sized inout character array: names(:) 1D_ALLOC_character
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_add_this_taylor_term(
    void* ele /* 0D_NOT_type */,
    c_Int& i_out /* 0D_NOT_integer */,
    c_Real& coef /* 0D_NOT_real */,
    c_IntArr expn /* 1D_NOT_integer */);
void add_this_taylor_term(
    EleProxy& ele,
    int i_out,
    double coef,
    FixedArray1D<Int, 6> expn);
extern "C" void fortran_adjust_super_slave_names(
    void* lat /* 0D_NOT_type */,
    c_Int& ix1_lord /* 0D_NOT_integer */,
    c_Int& ix2_lord /* 0D_NOT_integer */,
    c_Bool* first_time /* 0D_NOT_logical */);
void adjust_super_slave_names(
    LatProxy& lat,
    int ix1_lord,
    int ix2_lord,
    optional_ref<bool> first_time = std::nullopt);
extern "C" void fortran_allocate_branch_array(
    void* lat /* 0D_NOT_type */,
    c_Int& upper_bound /* 0D_NOT_integer */);
void allocate_branch_array(LatProxy& lat, int upper_bound);

// Skipped unusable routine allocate_element_array:
// Routine in configuration skip list
extern "C" void fortran_allocate_lat_ele_array(
    void* lat /* 0D_NOT_type */,
    c_Int* upper_bound /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* do_ramper_slave_setup /* 0D_NOT_logical */);
void allocate_lat_ele_array(
    LatProxy& lat,
    std::optional<int> upper_bound = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> do_ramper_slave_setup = std::nullopt);

// Skipped unusable routine allocate_plat:
// Untranslated type: ParserLatProxy (0D_NOT_type)

// Skipped unusable routine aml_parser:
// Routine in configuration skip list
extern "C" bool fortran_angle_between_polars(
    void* polar1 /* 0D_NOT_type */,
    void* polar2 /* 0D_NOT_type */,
    c_Real& angle /* 0D_NOT_real */);
void angle_between_polars(
    SpinPolarProxy& polar1,
    SpinPolarProxy& polar2,
    double angle);
extern "C" void fortran_angle_to_canonical_coords(
    void* orbit /* 0D_NOT_type */,
    c_Char coord_type /* 0D_NOT_character */);
void angle_to_canonical_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type = std::nullopt);

// Skipped unusable routine anormal_mode_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_aperture_bookkeeper(void* ele /* 0D_NOT_type */);
void aperture_bookkeeper(EleProxy& ele);

// Skipped unusable routine aperture_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine aperture_point_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine aperture_scan_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_apply_all_rampers(
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool apply_all_rampers(LatProxy& lat);

// Skipped unusable routine apply_element_edge_kick:
// Untranslated type: FringeFieldInfoProxy (0D_NOT_type)

// Skipped unusable routine apply_element_edge_kick_hook_def:
// Untranslated type: FringeFieldInfoProxy (0D_NOT_type)
extern "C" void fortran_apply_energy_kick(
    c_Real& dE /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr ddE_dr /* 1D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void apply_energy_kick(
    double dE,
    CoordProxy& orbit,
    FixedArray1D<Real, 2> ddE_dr,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_apply_patch_to_ptc_fibre(void* ele /* 0D_NOT_type */);
void apply_patch_to_ptc_fibre(EleProxy& ele);
extern "C" void fortran_apply_rampers_to_slave(
    void* slave /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool apply_rampers_to_slave(EleProxy& slave);
extern "C" bool fortran_array_re_str(
    void* arr /* 1D_ALLOC_real */,
    c_Char parens_in /* 0D_NOT_character */,
    c_Char str_out /* 0D_NOT_character */);
void array_re_str(
    RealAlloc1D& arr,
    optional_ref<std::string> parens_in,
    std::string str_out);

// Skipped unusable routine astra_lattice_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_astra_max_field_reference(
    void* pt0 /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& field_value /* 0D_NOT_real */);
void astra_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double field_value);
extern "C" bool fortran_at_this_ele_end(
    c_Int& now_at /* 0D_NOT_integer */,
    c_Int& where_at /* 0D_NOT_integer */,
    c_Bool& is_at_this_end /* 0D_NOT_logical */);
void at_this_ele_end(int now_at, int where_at, bool is_at_this_end);
extern "C" void fortran_attribute_bookkeeper(
    void* ele /* 0D_NOT_type */,
    c_Bool* force_bookkeeping /* 0D_NOT_logical */);
void attribute_bookkeeper(
    EleProxy& ele,
    std::optional<bool> force_bookkeeping = std::nullopt);
extern "C" bool fortran_attribute_free1(
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Char attrib_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Bool* except_overlay /* 0D_NOT_logical */,
    c_Bool* dependent_attribs_free /* 0D_NOT_logical */,
    c_Int* why_not_free /* 0D_NOT_integer */,
    c_Bool& free /* 0D_NOT_logical */);
void attribute_free1(
    int ix_ele,
    std::string attrib_name,
    LatProxy& lat,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool free);
extern "C" bool fortran_attribute_free2(
    void* ele /* 0D_NOT_type */,
    c_Char attrib_name /* 0D_NOT_character */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Bool* except_overlay /* 0D_NOT_logical */,
    c_Bool* dependent_attribs_free /* 0D_NOT_logical */,
    c_Int* why_not_free /* 0D_NOT_integer */,
    c_Bool& free /* 0D_NOT_logical */);
void attribute_free2(
    EleProxy& ele,
    std::string attrib_name,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool free);
extern "C" bool fortran_attribute_free3(
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Char attrib_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Bool* except_overlay /* 0D_NOT_logical */,
    c_Bool* dependent_attribs_free /* 0D_NOT_logical */,
    c_Int* why_not_free /* 0D_NOT_integer */,
    c_Bool& free /* 0D_NOT_logical */);
void attribute_free3(
    int ix_ele,
    int ix_branch,
    std::string attrib_name,
    LatProxy& lat,
    optional_ref<bool> err_print_flag,
    optional_ref<bool> except_overlay,
    optional_ref<bool> dependent_attribs_free,
    optional_ref<int> why_not_free,
    bool free);
extern "C" bool fortran_attribute_index1(
    void* ele /* 0D_NOT_type */,
    c_Char name /* 0D_NOT_character */,
    c_Char full_name /* 0D_NOT_character */,
    c_Bool* can_abbreviate /* 0D_NOT_logical */,
    c_Bool* print_error /* 0D_NOT_logical */,
    c_Int& attrib_index /* 0D_NOT_integer */);
void attribute_index1(
    EleProxy& ele,
    std::string name,
    optional_ref<std::string> full_name,
    optional_ref<bool> can_abbreviate,
    optional_ref<bool> print_error,
    int attrib_index);
extern "C" bool fortran_attribute_index2(
    c_Int& key /* 0D_NOT_integer */,
    c_Char name /* 0D_NOT_character */,
    c_Char full_name /* 0D_NOT_character */,
    c_Bool* can_abbreviate /* 0D_NOT_logical */,
    c_Bool* print_error /* 0D_NOT_logical */,
    c_Int& attrib_index /* 0D_NOT_integer */);
void attribute_index2(
    int key,
    std::string name,
    optional_ref<std::string> full_name,
    optional_ref<bool> can_abbreviate,
    optional_ref<bool> print_error,
    int attrib_index);

// Skipped unusable routine attribute_info:
// Untranslated type: EleAttributeProxy (0D_NOT_type)
extern "C" bool fortran_attribute_name1(
    c_Int& key /* 0D_NOT_integer */,
    c_Int& ix_att /* 0D_NOT_integer */,
    c_Bool* show_private /* 0D_NOT_logical */,
    c_Char attrib_name /* 0D_NOT_character */);
void attribute_name1(
    int key,
    int ix_att,
    optional_ref<bool> show_private,
    std::string attrib_name);
extern "C" bool fortran_attribute_name2(
    void* ele /* 0D_NOT_type */,
    c_Int& ix_att /* 0D_NOT_integer */,
    c_Bool* show_private /* 0D_NOT_logical */,
    c_Char attrib_name /* 0D_NOT_character */);
void attribute_name2(
    EleProxy& ele,
    int ix_att,
    optional_ref<bool> show_private,
    std::string attrib_name);

// Skipped unusable routine attribute_set_bookkeeping:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" bool fortran_attribute_type(
    c_Char attrib_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Int& attrib_type /* 0D_NOT_integer */);
int attribute_type(
    std::string attrib_name,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" bool fortran_attribute_units(
    c_Char attrib_name /* 0D_NOT_character */,
    c_Char unrecognized_units /* 0D_NOT_character */,
    c_Char attrib_units /* 0D_NOT_character */);
std::string attribute_units(
    std::string attrib_name,
    std::optional<std::string> unrecognized_units = std::nullopt);
extern "C" void fortran_attributes_need_bookkeeping(
    void* ele /* 0D_NOT_type */,
    void* dval /* 1D_ALLOC_real */);
RealAlloc1D attributes_need_bookkeeping(EleProxy& ele);
extern "C" void fortran_autoscale_phase_and_amp(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* scale_phase /* 0D_NOT_logical */,
    c_Bool* scale_amp /* 0D_NOT_logical */,
    c_Bool* call_bookkeeper /* 0D_NOT_logical */);
bool autoscale_phase_and_amp(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> scale_phase = std::nullopt,
    std::optional<bool> scale_amp = std::nullopt,
    std::optional<bool> call_bookkeeper = std::nullopt);
extern "C" bool fortran_average_twiss(
    c_Real& frac1 /* 0D_NOT_real */,
    void* twiss1 /* 0D_NOT_type */,
    void* twiss2 /* 0D_NOT_type */,
    void* ave_twiss /* 0D_NOT_type */);
void average_twiss(
    double frac1,
    TwissProxy& twiss1,
    TwissProxy& twiss2,
    TwissProxy& ave_twiss);

// Skipped unusable routine bane1:
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine base_line_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_bbi_kick(
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_RealArr sigma /* 1D_NOT_real */,
    c_RealArr nk /* 1D_NOT_real */,
    c_RealArr dnk /* 2D_NOT_real */);
struct BbiKick {
  FixedArray1D<Real, 2> nk;
  FixedArray2D<Real, 2, 2> dnk;
};
BbiKick bbi_kick(double x, double y, FixedArray1D<Real, 2> sigma);
extern "C" void fortran_bbi_slice_calc(
    void* ele /* 0D_NOT_type */,
    c_Int& n_slice /* 0D_NOT_integer */,
    void* z_slice /* 1D_ALLOC_real */);
RealAlloc1D bbi_slice_calc(EleProxy& ele, int n_slice);
extern "C" void fortran_beam_envelope_ibs(
    c_RealArr sigma_mat /* 2D_NOT_real */,
    c_RealArr ibs_mat /* 2D_NOT_real */,
    c_Bool& tail_cut /* 0D_NOT_logical */,
    c_Real& tau /* 0D_NOT_real */,
    c_Real& energy /* 0D_NOT_real */,
    c_Real& n_part /* 0D_NOT_real */,
    c_Int& species /* 0D_NOT_integer */);
FixedArray2D<Real, 6, 6> beam_envelope_ibs(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species);
extern "C" void fortran_beam_equal_beam(
    void* beam1 /* 0D_NOT_type */,
    void* beam2 /* 0D_NOT_type */);
void beam_equal_beam(BeamProxy& beam1, BeamProxy& beam2);

// Skipped unusable routine beam_init_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine beam_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_beam_tilts(
    c_RealArr S /* 2D_NOT_real */,
    c_Real& angle_xy /* 0D_NOT_real */,
    c_Real& angle_xz /* 0D_NOT_real */,
    c_Real& angle_yz /* 0D_NOT_real */,
    c_Real& angle_xpz /* 0D_NOT_real */,
    c_Real& angle_ypz /* 0D_NOT_real */);
struct BeamTilts {
  double angle_xy;
  double angle_xz;
  double angle_yz;
  double angle_xpz;
  double angle_ypz;
};
BeamTilts beam_tilts(FixedArray2D<Real, 6, 6> S);

// Skipped unusable routine beambeam_fibre_setup:
// Untranslated type: FibreProxy (0D_NOT_type)
extern "C" void fortran_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_Bool* track_spin /* 0D_NOT_logical */);
void bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> track_spin = std::nullopt);
extern "C" void fortran_bend_exact_multipole_field(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& local_ref_frame /* 0D_NOT_logical */,
    void* field /* 0D_NOT_type */,
    c_Bool* calc_dfield /* 0D_NOT_logical */,
    c_Bool* calc_potential /* 0D_NOT_logical */);
EmFieldProxy bend_exact_multipole_field(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" bool fortran_bend_length_has_been_set(
    void* ele /* 0D_NOT_type */,
    c_Bool& is_set /* 0D_NOT_logical */);
void bend_length_has_been_set(EleProxy& ele, bool is_set);
extern "C" bool fortran_bend_photon_e_rel_init(
    c_Real* r_in /* 0D_NOT_real */,
    c_Real& E_rel /* 0D_NOT_real */);
double bend_photon_e_rel_init(std::optional<double> r_in = std::nullopt);
extern "C" bool fortran_bend_photon_energy_integ_prob(
    c_Real& E_photon /* 0D_NOT_real */,
    c_Real& g_bend /* 0D_NOT_real */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& integ_prob /* 0D_NOT_real */);
double bend_photon_energy_integ_prob(
    double E_photon,
    double g_bend,
    double gamma);
extern "C" bool fortran_bend_photon_energy_normalized_probability(
    c_Real& E_rel /* 0D_NOT_real */,
    c_Real& prob /* 0D_NOT_real */);
double bend_photon_energy_normalized_probability(double E_rel);
extern "C" void fortran_bend_photon_init(
    c_Real& g_bend_x /* 0D_NOT_real */,
    c_Real& g_bend_y /* 0D_NOT_real */,
    c_Real& gamma /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_Real* E_min /* 0D_NOT_real */,
    c_Real* E_max /* 0D_NOT_real */,
    c_Real* E_integ_prob /* 0D_NOT_real */,
    c_Real* vert_angle_min /* 0D_NOT_real */,
    c_Real* vert_angle_max /* 0D_NOT_real */,
    c_Bool* vert_angle_symmetric /* 0D_NOT_logical */,
    c_Real* emit_probability /* 0D_NOT_real */);
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
    c_Real& g_bend_x /* 0D_NOT_real */,
    c_Real& g_bend_y /* 0D_NOT_real */,
    c_Real& E_rel /* 0D_NOT_real */,
    c_Real& gamma_phi /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */);
CoordProxy bend_photon_polarization_init(
    double g_bend_x,
    double g_bend_y,
    double E_rel,
    double gamma_phi);
extern "C" bool fortran_bend_photon_vert_angle_init(
    c_Real& E_rel /* 0D_NOT_real */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Real* r_in /* 0D_NOT_real */,
    c_Bool* invert /* 0D_NOT_logical */,
    c_Real& phi /* 0D_NOT_real */);
double bend_photon_vert_angle_init(
    double E_rel,
    double gamma,
    std::optional<double> r_in = std::nullopt,
    std::optional<bool> invert = std::nullopt);
extern "C" bool fortran_bend_shift(
    void* position1 /* 0D_NOT_type */,
    c_Real& g /* 0D_NOT_real */,
    c_Real& delta_s /* 0D_NOT_real */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Real* ref_tilt /* 0D_NOT_real */,
    void* position2 /* 0D_NOT_type */);
FixedArray2D<Real, 3, 3> bend_shift(
    FloorPositionProxy& position1,
    double g,
    double delta_s,
    std::optional<double> ref_tilt,
    FloorPositionProxy& position2);
extern "C" bool fortran_bend_vert_angle_integ_prob(
    c_Real& vert_angle /* 0D_NOT_real */,
    c_Real& E_rel /* 0D_NOT_real */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& integ_prob /* 0D_NOT_real */);
double bend_vert_angle_integ_prob(
    double vert_angle,
    double E_rel,
    double gamma);

// Skipped unusable routine bjmt1:
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine bjmt_integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)

// Skipped unusable routine bl_via_mat:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)
extern "C" void fortran_bl_via_vlassov(
    c_Real& current /* 0D_NOT_real */,
    c_Real& alpha /* 0D_NOT_real */,
    c_Real& Energy /* 0D_NOT_real */,
    c_Real& sigma_p /* 0D_NOT_real */,
    c_Real& Vrf /* 0D_NOT_real */,
    c_Real& omega /* 0D_NOT_real */,
    c_Real& U0 /* 0D_NOT_real */,
    c_Real& circ /* 0D_NOT_real */,
    c_Real& R /* 0D_NOT_real */,
    c_Real& L /* 0D_NOT_real */,
    c_Real& sigma_z /* 0D_NOT_real */);
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
// Routine in configuration skip list

// Skipped unusable routine bmad_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bmad_normal_form_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_bmad_parser(
    c_Char lat_file /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool* make_mats6 /* 0D_NOT_logical */,
    c_Bool& digested_read_ok /* 0D_NOT_logical */,
    c_Char use_line /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* parse_lat /* 0D_NOT_type */);
struct BmadParser {
  LatProxy lat;
  bool digested_read_ok;
  bool err_flag;
  LatProxy parse_lat;
};
BmadParser bmad_parser(
    std::string lat_file,
    std::optional<bool> make_mats6 = std::nullopt,
    std::optional<std::string> use_line = std::nullopt);
extern "C" void fortran_bmad_parser2(
    c_Char lat_file /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */,
    c_Bool* make_mats6 /* 0D_NOT_logical */,
    c_Bool* err_flag /* 0D_NOT_logical */,
    void* parse_lat /* 0D_NOT_type */);
void bmad_parser2(
    std::string lat_file,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> orbit = std::nullopt,
    std::optional<bool> make_mats6 = std::nullopt,
    optional_ref<bool> err_flag = std::nullopt,
    optional_ref<LatProxy> parse_lat = std::nullopt);

// Skipped unusable routine bmad_parser_string_attribute_set:
// Untranslated type: ParserEleProxy (0D_NOT_type)
extern "C" void fortran_bmad_patch_parameters_to_ptc(
    c_RealArr ang /* 1D_NOT_real */,
    c_RealArr exi /* 2D_NOT_real */);
void bmad_patch_parameters_to_ptc(
    FixedArray1D<Real, 3> ang,
    FixedArray2D<Real, 3, 3> exi);

// Skipped unusable routine bmad_private_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bmad_taylor_equal_damap:
// Untranslated type: DamapProxy (0D_NOT_type)

// Skipped unusable routine bmad_taylors_equal_ptc_taylors:
// Untranslated type: TaylorProxy (1D_ALLOC_type)

// Skipped unusable routine bmad_taylors_equal_reals_8:
// Untranslated type: Real8Proxy (1D_ALLOC_type)

// Skipped unusable routine bookkeeping_state_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bp_common2_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bp_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bp_const_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_bp_set_ran_status();
void bp_set_ran_status();

// Skipped unusable routine bpm_phase_coupling_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_branch_equal_branch(
    void* branch1 /* 0D_NOT_type */,
    void* branch2 /* 0D_NOT_type */);
BranchProxy branch_equal_branch(BranchProxy& branch2);
extern "C" bool fortran_branch_name(
    void* branch /* 0D_NOT_type */,
    c_Char name /* 0D_NOT_character */);
void branch_name(BranchProxy& branch, std::string name);

// Skipped unusable routine branch_pointer_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine branch_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_branch_to_ptc_m_u(void* branch /* 0D_NOT_type */);
void branch_to_ptc_m_u(BranchProxy& branch);
extern "C" void fortran_bunch_equal_bunch(
    void* bunch1 /* 0D_NOT_type */,
    void* bunch2 /* 0D_NOT_type */);
void bunch_equal_bunch(BunchProxy& bunch1, BunchProxy& bunch2);

// Skipped unusable routine bunch_params_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bunch_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine bunch_track_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine c_multi:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_c_to_cbar(
    void* ele /* 0D_NOT_type */,
    c_RealArr cbar_mat /* 2D_NOT_real */);
FixedArray2D<Real, 2, 2> c_to_cbar(EleProxy& ele);
extern "C" void fortran_calc_bunch_params(
    void* bunch /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */,
    c_Bool& error /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_RealArr n_mat /* 2D_NOT_real */,
    c_Bool* is_time_coords /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */);
void calc_bunch_params(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    bool error,
    std::optional<bool> print_err = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> n_mat = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" void fortran_calc_bunch_params_slice(
    void* bunch /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */,
    c_Int& plane /* 0D_NOT_integer */,
    c_Real& slice_center /* 0D_NOT_real */,
    c_Real& slice_spread /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Bool* is_time_coords /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */);
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
    void* bunch /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */,
    c_RealArr slice_bounds /* 1D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Bool* is_time_coords /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */);
void calc_bunch_params_z_slice(
    BunchProxy& bunch,
    BunchParamsProxy& bunch_params,
    FixedArray1D<Real, 2> slice_bounds,
    bool err,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);
extern "C" void fortran_calc_bunch_sigma_matrix_etc(
    void* particle /* 1D_ALLOC_type */,
    void* charge /* 1D_ALLOC_real */,
    void* bunch_params /* 0D_NOT_type */,
    c_Bool* is_time_coords /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */);
BunchParamsProxy calc_bunch_sigma_matrix_etc(
    CoordProxyAlloc1D& particle,
    RealAlloc1D& charge,
    optional_ref<bool> is_time_coords = std::nullopt,
    optional_ref<EleProxy> ele = std::nullopt);

// Skipped unusable routine calc_density_derivative_complex:
// Variable inout sized array: density(:,:,:) 3D_NOT_real
// Variable inout sized array: density_prime(:,:,:) 3D_NOT_complex
extern "C" void fortran_calc_emittances_and_twiss_from_sigma_matrix(
    c_RealArr sigma_mat /* 2D_NOT_real */,
    void* bunch_params /* 0D_NOT_type */,
    c_Bool& error /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_RealArr n_mat /* 2D_NOT_real */);
struct CalcEmittancesAndTwissFromSigmaMatrix {
  BunchParamsProxy bunch_params;
  bool error;
  std::optional<FixedArray2D<Real, 6, 6>> n_mat;
};
CalcEmittancesAndTwissFromSigmaMatrix
calc_emittances_and_twiss_from_sigma_matrix(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine calc_next_fringe_edge:
// Untranslated type: FringeFieldInfoProxy (0D_NOT_type)
extern "C" void fortran_calc_spin_params(
    void* bunch /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */);
BunchParamsProxy calc_spin_params(BunchProxy& bunch);
extern "C" void fortran_calc_super_slave_key(
    void* lord1 /* 0D_NOT_type */,
    void* lord2 /* 0D_NOT_type */,
    void* slave /* 0D_NOT_type */,
    c_Bool* create_jumbo_slave /* 0D_NOT_logical */);
EleProxy calc_super_slave_key(
    EleProxy& lord1,
    EleProxy& lord2,
    std::optional<bool> create_jumbo_slave = std::nullopt);
extern "C" void fortran_calc_wall_radius(
    void* v /* 1D_ALLOC_type */,
    c_Real& cos_ang /* 0D_NOT_real */,
    c_Real& sin_ang /* 0D_NOT_real */,
    c_Real& r_wall /* 0D_NOT_real */,
    c_Real& dr_dtheta /* 0D_NOT_real */,
    c_Int& ix_vertex /* 0D_NOT_integer */);
struct CalcWallRadius {
  double r_wall;
  double dr_dtheta;
  int ix_vertex;
};
CalcWallRadius calc_wall_radius(
    Wall3dVertexProxyAlloc1D& v,
    double cos_ang,
    double sin_ang);

// Skipped unusable routine calc_wiggler_g_params:
// Untranslated type: RadIntTrackPointProxy (0D_NOT_type)
// Untranslated type: RadIntInfoProxy (0D_NOT_type)
extern "C" void fortran_calc_z_tune(void* branch /* 0D_NOT_type */);
void calc_z_tune(BranchProxy& branch);
extern "C" void fortran_canonical_to_angle_coords(
    void* orbit /* 0D_NOT_type */,
    c_Char coord_type /* 0D_NOT_character */);
void canonical_to_angle_coords(
    CoordProxy& orbit,
    std::optional<std::string> coord_type = std::nullopt);

// Skipped unusable routine capillary_photon_hit_spot_calc:
// Untranslated type: PhotonTrackProxy (0D_NOT_type)

// Skipped unusable routine capillary_propagate_photon_a_step:
// Untranslated type: PhotonTrackProxy (0D_NOT_type)

// Skipped unusable routine capillary_reflect_photon:
// Untranslated type: PhotonTrackProxy (0D_NOT_type)

// Skipped unusable routine capillary_track_photon_to_wall:
// Untranslated type: PhotonTrackProxy (0D_NOT_type)

// Skipped unusable routine cartesian_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine cartesian_map_term1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine cartesian_map_term_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_cbar_to_c(
    c_RealArr cbar_mat /* 2D_NOT_real */,
    void* a /* 0D_NOT_type */,
    void* b /* 0D_NOT_type */,
    c_RealArr c_mat /* 2D_NOT_real */);
FixedArray2D<Real, 2, 2> cbar_to_c(
    FixedArray2D<Real, 2, 2> cbar_mat,
    TwissProxy& a,
    TwissProxy& b);

// Skipped unusable routine ccfft3d:
// Routine module (fft_interface_mod) in configuration skip list

// Skipped unusable routine ccfftam:
// Routine in configuration skip list

// Skipped unusable routine cheb_diffuse_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_check_aperture_limit(
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* param /* 0D_NOT_type */,
    void* old_orb /* 0D_NOT_type */,
    c_Bool* check_momentum /* 0D_NOT_logical */);
void check_aperture_limit(
    CoordProxy& orb,
    EleProxy& ele,
    int particle_at,
    LatParamProxy& param,
    optional_ref<CoordProxy> old_orb = std::nullopt,
    std::optional<bool> check_momentum = std::nullopt);

// Skipped unusable routine check_aperture_limit_custom_def:
// Routine in configuration skip list
extern "C" void fortran_check_controller_controls(
    c_Int& ele_key /* 0D_NOT_integer */,
    void* contrl /* 1D_ALLOC_type */,
    c_Char name /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
bool check_controller_controls(
    int ele_key,
    ControlProxyAlloc1D& contrl,
    std::string name);
extern "C" void fortran_check_for_superimpose_problem(
    void* branch /* 0D_NOT_type */,
    void* super_ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ref_ele /* 0D_NOT_type */,
    c_Bool& wrap /* 0D_NOT_logical */);
void check_for_superimpose_problem(
    BranchProxy& branch,
    EleProxy& super_ele,
    bool err_flag,
    optional_ref<EleProxy> ref_ele,
    bool wrap);
extern "C" void fortran_check_if_s_in_bounds(
    void* branch /* 0D_NOT_type */,
    c_Real& s /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& translated_s /* 0D_NOT_real */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct CheckIfSInBounds {
  bool err_flag;
  double translated_s;
};
CheckIfSInBounds check_if_s_in_bounds(
    BranchProxy& branch,
    double s,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine choose_quads_for_set_tune:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" void fortran_chrom_calc(
    void* lat /* 0D_NOT_type */,
    c_Real& delta_e /* 0D_NOT_real */,
    c_Real& chrom_a /* 0D_NOT_real */,
    c_Real& chrom_b /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real* pz /* 0D_NOT_real */,
    void* low_E_lat /* 0D_NOT_type */,
    void* high_E_lat /* 0D_NOT_type */,
    void* low_E_orb /* 1D_ALLOC_type */,
    void* high_E_orb /* 1D_ALLOC_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    void* orb0 /* 0D_NOT_type */);
struct ChromCalc {
  double chrom_a;
  double chrom_b;
  bool err_flag;
  LatProxy low_E_lat;
  LatProxy high_E_lat;
  CoordProxyAlloc1D low_E_orb;
  CoordProxyAlloc1D high_E_orb;
};
ChromCalc chrom_calc(
    LatProxy& lat,
    double delta_e,
    std::optional<double> pz = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<CoordProxy> orb0 = std::nullopt);
extern "C" void fortran_chrom_tune(
    void* lat /* 0D_NOT_type */,
    c_Real& delta_e /* 0D_NOT_real */,
    c_Real& target_x /* 0D_NOT_real */,
    c_Real& target_y /* 0D_NOT_real */,
    c_Real& err_tol /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool chrom_tune(
    LatProxy& lat,
    double delta_e,
    double target_x,
    double target_y,
    double err_tol);

// Skipped unusable routine cimp1:
// Untranslated type: IbsProxy (0D_NOT_type)
extern "C" bool fortran_classical_radius(
    c_Int& species /* 0D_NOT_integer */,
    c_Real& radius /* 0D_NOT_real */);
void classical_radius(int species, double radius);
extern "C" void fortran_clear_lat_1turn_mats(void* lat /* 0D_NOT_type */);
LatProxy clear_lat_1turn_mats();
extern "C" void fortran_clear_taylor_maps_from_elements(
    void* lat /* 0D_NOT_type */);
void clear_taylor_maps_from_elements(LatProxy& lat);
extern "C" void fortran_closed_orbit_calc(
    void* lat /* 0D_NOT_type */,
    void* closed_orb /* 1D_ALLOC_type */,
    c_Int* i_dim /* 0D_NOT_integer */,
    c_Int* direction /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
bool closed_orbit_calc(
    LatProxy& lat,
    CoordProxyAlloc1D& closed_orb,
    std::optional<int> i_dim = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_closed_orbit_from_tracking(
    void* lat /* 0D_NOT_type */,
    void* closed_orb /* 1D_ALLOC_type */,
    c_Int& i_dim /* 0D_NOT_integer */,
    void* eps_rel /* 1D_ALLOC_real */,
    void* eps_abs /* 1D_ALLOC_real */,
    void* init_guess /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ClosedOrbitFromTracking {
  CoordProxyAlloc1D closed_orb;
  bool err_flag;
};
ClosedOrbitFromTracking closed_orbit_from_tracking(
    LatProxy& lat,
    int i_dim,
    optional_ref<RealAlloc1D> eps_rel = std::nullopt,
    optional_ref<RealAlloc1D> eps_abs = std::nullopt,
    optional_ref<CoordProxy> init_guess = std::nullopt);
extern "C" bool fortran_cmplx_re_str(
    c_Complex& cmp /* 0D_NOT_complex */,
    c_Char str_out /* 0D_NOT_character */);
void cmplx_re_str(std::complex<double> cmp, std::string str_out);
extern "C" void fortran_combine_consecutive_elements(
    void* lat /* 0D_NOT_type */,
    c_Bool& error /* 0D_NOT_logical */);
bool combine_consecutive_elements(LatProxy& lat);
extern "C" void fortran_complex_taylor_clean(
    void* complex_taylor /* 0D_NOT_type */);
void complex_taylor_clean(ComplexTaylorProxy& complex_taylor);

// Skipped unusable routine complex_taylor_equal_c_taylor:
// Untranslated type: CTaylorProxy (0D_NOT_type)
extern "C" void fortran_complex_taylor_equal_complex_taylor(
    void* complex_taylor1 /* 0D_NOT_type */,
    void* complex_taylor2 /* 0D_NOT_type */);
ComplexTaylorProxy complex_taylor_equal_complex_taylor(
    ComplexTaylorProxy& complex_taylor2);
extern "C" bool fortran_complex_taylor_exponent_index(
    c_IntArr expn /* 1D_NOT_integer */,
    c_Int& index /* 0D_NOT_integer */);
int complex_taylor_exponent_index(FixedArray1D<Int, 6> expn);
extern "C" void fortran_complex_taylor_make_unit(
    void* complex_taylor /* 1D_ALLOC_type */);
ComplexTaylorProxyAlloc1D complex_taylor_make_unit();

// Skipped unusable routine complex_taylor_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine complex_taylor_term_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_complex_taylor_to_mat6(
    void* a_complex_taylor /* 1D_NOT_type */,
    void* r_in /* 1D_ALLOC_complex */,
    c_ComplexArr vec0 /* 1D_NOT_complex */,
    c_ComplexArr mat6 /* 2D_NOT_complex */,
    void* r_out /* 1D_ALLOC_complex */);
struct ComplexTaylorToMat6 {
  FixedArray1D<Complex, 6> vec0;
  FixedArray2D<Complex, 6, 6> mat6;
  ComplexAlloc1D r_out;
};
ComplexTaylorToMat6 complex_taylor_to_mat6(
    FixedArray1D<ComplexTaylorProxy, 6> a_complex_taylor,
    ComplexAlloc1D& r_in);

// Skipped unusable routine complex_taylors_equal_c_taylors:
// Untranslated type: CTaylorProxy (1D_ALLOC_type)
extern "C" void fortran_complex_taylors_equal_complex_taylors(
    void* complex_taylor1 /* 1D_ALLOC_type */,
    void* complex_taylor2 /* 1D_ALLOC_type */);
ComplexTaylorProxyAlloc1D complex_taylors_equal_complex_taylors(
    ComplexTaylorProxyAlloc1D& complex_taylor2);

// Skipped unusable routine complex_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine compounddatanist_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_compute_slave_coupler(void* slave /* 0D_NOT_type */);
void compute_slave_coupler(EleProxy& slave);

// Skipped unusable routine compute_super_lord_s:
// Untranslated type: ParserEleProxy (0D_NOT_type)
extern "C" void fortran_concat_ele_taylor(
    void* orb_taylor /* 1D_ALLOC_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* spin_taylor /* 1D_ALLOC_type */);
void concat_ele_taylor(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    bool err_flag,
    optional_ref<TaylorProxyAlloc1D> spin_taylor = std::nullopt);

// Skipped unusable routine concat_real_8:
// Untranslated type: Real8Proxy (1D_ALLOC_type)
// Untranslated type: Real8Proxy (1D_ALLOC_type)
// Untranslated type: Real8Proxy (1D_ALLOC_type)
extern "C" void fortran_concat_taylor(
    void* taylor1 /* 1D_ALLOC_type */,
    void* taylor2 /* 1D_ALLOC_type */,
    void* taylor3 /* 1D_ALLOC_type */);
void concat_taylor(
    TaylorProxyAlloc1D& taylor1,
    TaylorProxyAlloc1D& taylor2,
    TaylorProxyAlloc1D& taylor3);
extern "C" void fortran_concat_transfer_mat(
    c_RealArr mat_1 /* 2D_NOT_real */,
    c_RealArr vec_1 /* 1D_NOT_real */,
    c_RealArr mat_0 /* 2D_NOT_real */,
    c_RealArr vec_0 /* 1D_NOT_real */,
    c_RealArr mat_out /* 2D_NOT_real */,
    c_RealArr vec_out /* 1D_NOT_real */);
FixedArray2D<Real, 6, 6> concat_transfer_mat(
    FixedArray2D<Real, 6, 6> mat_1,
    FixedArray1D<Real, 6> vec_1,
    FixedArray2D<Real, 6, 6> mat_0,
    FixedArray1D<Real, 6> vec_0,
    FixedArray1D<Real, 6> vec_out);
extern "C" void fortran_control_bookkeeper(
    void* lat /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool* err_flag /* 0D_NOT_logical */);
void control_bookkeeper(
    LatProxy& lat,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<bool> err_flag = std::nullopt);

// Skipped unusable routine control_ramp1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine control_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine control_var1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine controller_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine conv3d:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_convert_bend_exact_multipole(
    c_Real& g /* 0D_NOT_real */,
    c_Int& out_type /* 0D_NOT_integer */,
    c_RealArr an /* 1D_NOT_real */,
    c_RealArr bn /* 1D_NOT_real */);
void convert_bend_exact_multipole(
    double g,
    int out_type,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> an,
    FixedArray1D<Real, Bmad::N_POLE_MAXX> bn);
extern "C" void fortran_convert_coords(
    c_Char in_type_str /* 0D_NOT_character */,
    void* coord_in /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Char out_type_str /* 0D_NOT_character */,
    void* coord_out /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ConvertCoords {
  std::string out_type_str;
  CoordProxy coord_out;
  bool err_flag;
};
ConvertCoords convert_coords(
    std::string in_type_str,
    CoordProxy& coord_in,
    EleProxy& ele);
extern "C" void fortran_convert_field_ele_to_lab(
    void* ele /* 0D_NOT_type */,
    c_Real& s_here /* 0D_NOT_real */,
    c_Bool& forward_transform /* 0D_NOT_logical */,
    void* field /* 0D_NOT_type */,
    c_Bool* calc_dfield /* 0D_NOT_logical */,
    c_Bool* calc_potential /* 0D_NOT_logical */);
EmFieldProxy convert_field_ele_to_lab(
    EleProxy& ele,
    double s_here,
    bool forward_transform,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" void fortran_convert_local_cartesian_to_local_curvilinear(
    c_Real& x /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& g /* 0D_NOT_real */,
    c_Real& xout /* 0D_NOT_real */,
    c_Real& sout /* 0D_NOT_real */);
void convert_local_cartesian_to_local_curvilinear(
    double x,
    double z,
    double g,
    double xout,
    double sout);
extern "C" void fortran_convert_local_curvilinear_to_local_cartesian(
    c_Real& x /* 0D_NOT_real */,
    c_Real& s /* 0D_NOT_real */,
    c_Real& g /* 0D_NOT_real */,
    c_Real& xout /* 0D_NOT_real */,
    c_Real& zout /* 0D_NOT_real */);
void convert_local_curvilinear_to_local_cartesian(
    double x,
    double s,
    double g,
    double xout,
    double zout);
extern "C" void fortran_convert_particle_coordinates_s_to_t(
    void* particle /* 0D_NOT_type */,
    c_Real& s_body /* 0D_NOT_real */,
    c_Int& orientation /* 0D_NOT_integer */);
void convert_particle_coordinates_s_to_t(
    CoordProxy& particle,
    double s_body,
    int orientation);
extern "C" void fortran_convert_particle_coordinates_t_to_s(
    void* particle /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& s_body /* 0D_NOT_real */,
    c_Bool* use_downstream_p0c /* 0D_NOT_logical */);
double convert_particle_coordinates_t_to_s(
    CoordProxy& particle,
    EleProxy& ele,
    std::optional<bool> use_downstream_p0c = std::nullopt);
extern "C" void fortran_convert_pc_to(
    c_Real& pc /* 0D_NOT_real */,
    c_Int& particle /* 0D_NOT_integer */,
    c_Real& E_tot /* 0D_NOT_real */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& kinetic /* 0D_NOT_real */,
    c_Real& beta /* 0D_NOT_real */,
    c_Real& brho /* 0D_NOT_real */,
    c_Real& beta1 /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ConvertPcTo {
  double E_tot;
  double gamma;
  double kinetic;
  double beta;
  double brho;
  double beta1;
  bool err_flag;
};
ConvertPcTo convert_pc_to(double pc, int particle);
extern "C" void fortran_convert_total_energy_to(
    c_Real& E_tot /* 0D_NOT_real */,
    c_Int& particle /* 0D_NOT_integer */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& kinetic /* 0D_NOT_real */,
    c_Real& beta /* 0D_NOT_real */,
    c_Real& pc /* 0D_NOT_real */,
    c_Real& brho /* 0D_NOT_real */,
    c_Real& beta1 /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct ConvertTotalEnergyTo {
  double gamma;
  double kinetic;
  double beta;
  double pc;
  double brho;
  double beta1;
  bool err_flag;
};
ConvertTotalEnergyTo convert_total_energy_to(
    double E_tot,
    int particle,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine converter_dir_1d_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_dir_2d_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_dir_coef_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_direction_out_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_converter_distribution_parser(
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ConverterDistributionParser {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
ConverterDistributionParser converter_distribution_parser(EleProxy& ele);

// Skipped unusable routine converter_distribution_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_prob_pc_r_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine converter_sub_distribution_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine coord_array_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_coord_equal_coord(
    void* coord1 /* 0D_NOT_type */,
    void* coord2 /* 0D_NOT_type */);
CoordProxy coord_equal_coord(CoordProxy& coord2);
extern "C" bool fortran_coord_state_name(
    c_Int& coord_state /* 0D_NOT_integer */,
    c_Bool* one_word /* 0D_NOT_logical */,
    c_Char state_str /* 0D_NOT_character */);
std::string coord_state_name(
    int coord_state,
    optional_ref<bool> one_word = std::nullopt);

// Skipped unusable routine coord_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_coords_body_to_local(
    void* body_position /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Bool* calculate_angles /* 0D_NOT_logical */,
    void* local_position /* 0D_NOT_type */);
void coords_body_to_local(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_body_to_rel_exit(
    void* body_position /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Bool* calculate_angles /* 0D_NOT_logical */,
    void* rel_exit /* 0D_NOT_type */);
void coords_body_to_rel_exit(
    FloorPositionProxy& body_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& rel_exit);
extern "C" bool fortran_coords_curvilinear_to_floor(
    c_RealArr xys /* 1D_NOT_real */,
    void* branch /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* global /* 0D_NOT_type */);
bool coords_curvilinear_to_floor(
    FixedArray1D<Real, 3> xys,
    BranchProxy& branch,
    FloorPositionProxy& global);
extern "C" bool fortran_coords_floor_to_curvilinear(
    void* floor_coords /* 0D_NOT_type */,
    void* ele0 /* 0D_NOT_type */,
    void* ele1 /* 0D_PTR_type */,
    c_Int& status /* 0D_NOT_integer */,
    c_RealArr w_mat /* 2D_NOT_real */,
    void* local_coords /* 0D_NOT_type */);
struct CoordsFloorToCurvilinear {
  EleProxy ele1;
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
};
CoordsFloorToCurvilinear coords_floor_to_curvilinear(
    FloorPositionProxy& floor_coords,
    EleProxy& ele0,
    FloorPositionProxy& local_coords);
extern "C" bool fortran_coords_floor_to_local_curvilinear(
    void* global_position /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int& status /* 0D_NOT_integer */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Int* relative_to /* 0D_NOT_integer */,
    void* local_position /* 0D_NOT_type */);
struct CoordsFloorToLocalCurvilinear {
  int status;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
};
CoordsFloorToLocalCurvilinear coords_floor_to_local_curvilinear(
    FloorPositionProxy& global_position,
    EleProxy& ele,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_floor_to_relative(
    void* floor0 /* 0D_NOT_type */,
    void* global_position /* 0D_NOT_type */,
    c_Bool* calculate_angles /* 0D_NOT_logical */,
    c_Bool* is_delta_position /* 0D_NOT_logical */,
    void* local_position /* 0D_NOT_type */);
void coords_floor_to_relative(
    FloorPositionProxy& floor0,
    FloorPositionProxy& global_position,
    std::optional<bool> calculate_angles,
    std::optional<bool> is_delta_position,
    FloorPositionProxy& local_position);
extern "C" bool fortran_coords_local_curvilinear_to_body(
    void* local_position /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Bool* calculate_angles /* 0D_NOT_logical */,
    void* body_position /* 0D_NOT_type */);
void coords_local_curvilinear_to_body(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat,
    std::optional<bool> calculate_angles,
    FloorPositionProxy& body_position);
extern "C" bool fortran_coords_local_curvilinear_to_floor(
    void* local_position /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool* in_body_frame /* 0D_NOT_logical */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Bool* calculate_angles /* 0D_NOT_logical */,
    c_Int* relative_to /* 0D_NOT_integer */,
    void* global_position /* 0D_NOT_type */);
FixedArray2D<Real, 3, 3> coords_local_curvilinear_to_floor(
    FloorPositionProxy& local_position,
    EleProxy& ele,
    std::optional<bool> in_body_frame,
    std::optional<bool> calculate_angles,
    std::optional<int> relative_to,
    FloorPositionProxy& global_position);
extern "C" bool fortran_coords_relative_to_floor(
    void* floor0 /* 0D_NOT_type */,
    c_RealArr dr /* 1D_NOT_real */,
    c_Real* theta /* 0D_NOT_real */,
    c_Real* phi /* 0D_NOT_real */,
    c_Real* psi /* 0D_NOT_real */,
    void* floor1 /* 0D_NOT_type */);
void coords_relative_to_floor(
    FloorPositionProxy& floor0,
    FixedArray1D<Real, 3> dr,
    optional_ref<double> theta,
    optional_ref<double> phi,
    optional_ref<double> psi,
    FloorPositionProxy& floor1);

// Skipped unusable routine cos_phi:
// Untranslated type: DiffuseParamProxy (0D_NOT_type)
extern "C" bool fortran_coulombfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void coulombfun(double u, double v, double w, double gam, double res);
extern "C" void fortran_create_concatenated_wall3d(
    void* lat /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
void create_concatenated_wall3d(LatProxy& lat, bool err);
extern "C" void fortran_create_element_slice(
    void* sliced_ele /* 0D_NOT_type */,
    void* ele_in /* 0D_NOT_type */,
    c_Real& l_slice /* 0D_NOT_real */,
    c_Real& offset /* 0D_NOT_real */,
    void* param /* 0D_NOT_type */,
    c_Bool& include_upstream_end /* 0D_NOT_logical */,
    c_Bool& include_downstream_end /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* old_slice /* 0D_NOT_type */,
    void* orb_in /* 0D_NOT_type */);
struct CreateElementSlice {
  EleProxy sliced_ele;
  bool err_flag;
};
CreateElementSlice create_element_slice(
    EleProxy& ele_in,
    double l_slice,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end,
    optional_ref<EleProxy> old_slice = std::nullopt,
    optional_ref<CoordProxy> orb_in = std::nullopt);

// Skipped unusable routine create_feedback:
// Variable-sized in character array: input(:) 1D_ALLOC_character
// Variable-sized in character array: output(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_create_field_overlap(
    void* lat /* 0D_NOT_type */,
    c_Char lord_name /* 0D_NOT_character */,
    c_Char slave_name /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool create_field_overlap(
    LatProxy& lat,
    std::string lord_name,
    std::string slave_name);
extern "C" void fortran_create_girder(
    void* lat /* 0D_NOT_type */,
    c_Int& ix_girder /* 0D_NOT_integer */,
    void* contrl /* 1D_ALLOC_type */,
    void* girder_info /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void create_girder(
    LatProxy& lat,
    int ix_girder,
    ControlProxyAlloc1D& contrl,
    EleProxy& girder_info,
    bool err_flag);
extern "C" void fortran_create_group(
    void* lord /* 0D_NOT_type */,
    void* contrl /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */);
void create_group(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);

// Skipped unusable routine create_lat_ele_nametable:
// Untranslated type: NametableProxy (0D_NOT_type)
extern "C" void fortran_create_overlay(
    void* lord /* 0D_NOT_type */,
    void* contrl /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */);
void create_overlay(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);
extern "C" void fortran_create_planar_wiggler_model(
    void* wiggler_in /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct CreatePlanarWigglerModel {
  LatProxy lat;
  bool err_flag;
};
CreatePlanarWigglerModel create_planar_wiggler_model(
    EleProxy& wiggler_in,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_create_ramper(
    void* lord /* 0D_NOT_type */,
    void* contrl /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */);
void create_ramper(EleProxy& lord, ControlProxyAlloc1D& contrl, bool err);
extern "C" void fortran_create_sol_quad_model(
    void* sol_quad /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */);
void create_sol_quad_model(EleProxy& sol_quad, LatProxy& lat);
extern "C" void fortran_create_unique_ele_names(
    void* lat /* 0D_NOT_type */,
    c_Int& key /* 0D_NOT_integer */,
    c_Char suffix /* 0D_NOT_character */);
void create_unique_ele_names(LatProxy& lat, int key, std::string suffix);
extern "C" void fortran_create_wiggler_cartesian_map(
    void* ele /* 0D_NOT_type */,
    void* cart_map /* 0D_NOT_type */);
CartesianMapProxy create_wiggler_cartesian_map(EleProxy& ele);
extern "C" void fortran_crystal_attribute_bookkeeper(
    void* ele /* 0D_NOT_type */);
void crystal_attribute_bookkeeper(EleProxy& ele);

// Skipped unusable routine crystal_diffraction_field_calc:
// Untranslated type: CrystalParamProxy (0D_NOT_type)
extern "C" void fortran_crystal_h_misalign(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr h_vec /* 1D_NOT_real */);
void crystal_h_misalign(
    EleProxy& ele,
    CoordProxy& orbit,
    FixedArray1D<Real, 3> h_vec);

// Skipped unusable routine crystal_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine crystal_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_crystal_type_to_crystal_params(
    void* ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool crystal_type_to_crystal_params(EleProxy& ele);

// Skipped unusable routine csr3d_steady_state_solver:
// Variable in sized array: density(:,:,:) 3D_NOT_real
// Variable inout sized array: wake(:,:,:,:) 4D_NOT_real

// Skipped unusable routine csr_and_sc_apply_kicks:
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine csr_bin_kicks:
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine csr_bin_particles:
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine csr_bunch_slice_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine csr_ele_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine csr_kick1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine csr_particle_position_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine csr_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_custom_attribute_ubound_index(
    c_Int& ele_class /* 0D_NOT_integer */,
    c_Int& ix_ubound /* 0D_NOT_integer */);
int custom_attribute_ubound_index(int ele_class);

// Skipped unusable routine custom_ele_attrib_name_list:
// Variable-sized out character array: name_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine cylindrical_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine cylindrical_map_term1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine cylindrical_map_term_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine damap_equal_bmad_taylor:
// Untranslated type: DamapProxy (0D_NOT_type)
extern "C" bool fortran_damping_matrix_d(
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& g_tot /* 0D_NOT_real */,
    c_Real& B0 /* 0D_NOT_real */,
    c_Real& B1 /* 0D_NOT_real */,
    c_Real& delta /* 0D_NOT_real */,
    c_Int& species /* 0D_NOT_integer */,
    c_RealArr mat /* 2D_NOT_real */);
void damping_matrix_d(
    double gamma,
    double g_tot,
    double B0,
    double B1,
    double delta,
    int species,
    FixedArray2D<Real, 6, 6> mat);

// Skipped unusable routine deallocate_ele_array_pointers:
// Routine in configuration skip list
extern "C" void fortran_deallocate_ele_pointers(
    void* ele /* 0D_NOT_type */,
    c_Bool* nullify_only /* 0D_NOT_logical */,
    c_Bool* nullify_branch /* 0D_NOT_logical */,
    c_Bool* dealloc_poles /* 0D_NOT_logical */);
void deallocate_ele_pointers(
    EleProxy& ele,
    std::optional<bool> nullify_only = std::nullopt,
    std::optional<bool> nullify_branch = std::nullopt,
    std::optional<bool> dealloc_poles = std::nullopt);

// Skipped unusable routine deallocate_expression_tree:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)
extern "C" void fortran_deallocate_lat_pointers(void* lat /* 0D_NOT_type */);
void deallocate_lat_pointers(LatProxy& lat);

// Skipped unusable routine deallocate_tree:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)
extern "C" bool fortran_default_tracking_species(
    void* param /* 0D_NOT_type */,
    c_Int& species /* 0D_NOT_integer */);
void default_tracking_species(LatParamProxy& param, int species);

// Skipped unusable routine deposit_particles:
// Untranslated type: Mesh3dProxy (0D_NOT_type)
extern "C" bool fortran_detector_pixel_pt(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_IntArr ix_pix /* 1D_NOT_integer */);
FixedArray1D<Int, 2> detector_pixel_pt(CoordProxy& orbit, EleProxy& ele);
extern "C" bool fortran_diffraction_plate_or_mask_hit_spot(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int& ix_section /* 0D_NOT_integer */);
void diffraction_plate_or_mask_hit_spot(
    EleProxy& ele,
    CoordProxy& orbit,
    int ix_section);

// Skipped unusable routine diffuse_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine diffuse_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_diffusion_matrix_b(
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& g_tot /* 0D_NOT_real */,
    c_Int& species /* 0D_NOT_integer */,
    c_RealArr mat /* 2D_NOT_real */);
void diffusion_matrix_b(
    double gamma,
    double g_tot,
    int species,
    FixedArray2D<Real, 6, 6> mat);
extern "C" bool fortran_distance_to_aperture(
    void* orbit /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Bool& no_aperture_here /* 0D_NOT_logical */,
    c_Real& dist /* 0D_NOT_real */);
bool distance_to_aperture(
    CoordProxy& orbit,
    int particle_at,
    EleProxy& ele,
    double dist);

// Skipped unusable routine distance_to_aperture_custom_def:
// Routine in configuration skip list
extern "C" void fortran_do_mode_flip(
    void* ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool do_mode_flip(EleProxy& ele);
extern "C" bool fortran_dpc_given_de(
    c_Real& pc_old /* 0D_NOT_real */,
    c_Real& mass /* 0D_NOT_real */,
    c_Real& dE /* 0D_NOT_real */,
    c_Real& dpc /* 0D_NOT_real */);
void dpc_given_de(double pc_old, double mass, double dE, double dpc);
extern "C" void fortran_drift_and_pipe_track_methods_adjustment(
    void* lat /* 0D_NOT_type */);
void drift_and_pipe_track_methods_adjustment(LatProxy& lat);
extern "C" void fortran_drift_multipass_name_correction(
    void* lat /* 0D_NOT_type */);
void drift_multipass_name_correction(LatProxy& lat);
extern "C" void fortran_drift_orbit_time(
    void* orbit /* 0D_NOT_type */,
    c_Real& beta0 /* 0D_NOT_real */,
    c_Real* delta_s /* 0D_NOT_real */,
    c_Real* delta_t /* 0D_NOT_real */);
void drift_orbit_time(
    CoordProxy& orbit,
    double beta0,
    std::optional<double> delta_s = std::nullopt,
    std::optional<double> delta_t = std::nullopt);
extern "C" void fortran_drift_particle_to_s(
    void* p /* 0D_NOT_type */,
    c_Real& s /* 0D_NOT_real */,
    void* branch /* 0D_NOT_type */);
void drift_particle_to_s(CoordProxy& p, double s, BranchProxy& branch);
extern "C" void fortran_drift_particle_to_t(
    void* p /* 0D_NOT_type */,
    c_Real& t /* 0D_NOT_real */,
    void* branch /* 0D_NOT_type */);
void drift_particle_to_t(CoordProxy& p, double t, BranchProxy& branch);
extern "C" bool fortran_dspline_len(
    c_Real& s_chord0 /* 0D_NOT_real */,
    c_Real& s_chord1 /* 0D_NOT_real */,
    void* spline /* 0D_NOT_type */,
    c_Real* dtheta_ref /* 0D_NOT_real */,
    c_Real& dlen /* 0D_NOT_real */);
double dspline_len(
    double s_chord0,
    double s_chord1,
    SplineProxy& spline,
    std::optional<double> dtheta_ref = std::nullopt);
extern "C" void fortran_dynamic_aperture_point(
    void* branch /* 0D_NOT_type */,
    void* ele0 /* 0D_NOT_type */,
    void* orb0 /* 0D_NOT_type */,
    c_Real& theta_xy /* 0D_NOT_real */,
    void* ap_param /* 0D_NOT_type */,
    void* ap_point /* 0D_NOT_type */,
    c_Bool* check_xy_init /* 0D_NOT_logical */);
AperturePointProxy dynamic_aperture_point(
    BranchProxy& branch,
    EleProxy& ele0,
    CoordProxy& orb0,
    double theta_xy,
    ApertureParamProxy& ap_param,
    std::optional<bool> check_xy_init = std::nullopt);
extern "C" void fortran_dynamic_aperture_scan(
    void* aperture_scan /* 1D_ALLOC_type */,
    void* aperture_param /* 0D_NOT_type */,
    void* pz_start /* 1D_ALLOC_real */,
    void* lat /* 0D_NOT_type */,
    c_Bool* print_timing /* 0D_NOT_logical */);
ApertureScanProxyAlloc1D dynamic_aperture_scan(
    ApertureParamProxy& aperture_param,
    RealAlloc1D& pz_start,
    LatProxy& lat,
    std::optional<bool> print_timing = std::nullopt);
extern "C" bool fortran_e_accel_field(
    void* ele /* 0D_NOT_type */,
    c_Int& voltage_or_gradient /* 0D_NOT_integer */,
    c_Bool* bmad_standard_tracking /* 0D_NOT_logical */,
    c_Real& field /* 0D_NOT_real */);
void e_accel_field(
    EleProxy& ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking,
    double field);
extern "C" bool fortran_e_crit_photon(
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& g_bend /* 0D_NOT_real */,
    c_Real& E_crit /* 0D_NOT_real */);
double e_crit_photon(double gamma, double g_bend);
extern "C" void fortran_eigen_decomp_6mat(
    c_RealArr mat /* 2D_NOT_real */,
    c_ComplexArr eval /* 1D_NOT_complex */,
    c_ComplexArr evec /* 2D_NOT_complex */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_RealArr tunes /* 1D_NOT_real */);
struct EigenDecomp6mat {
  FixedArray1D<Complex, 6> eval;
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
  FixedArray1D<Real, 3> tunes;
};
EigenDecomp6mat eigen_decomp_6mat(FixedArray2D<Real, 6, 6> mat);

// Skipped unusable routine eigensys:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine ele_attribute_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ele_compute_ref_energy_and_time(
    void* ele0 /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void ele_compute_ref_energy_and_time(
    EleProxy& ele0,
    EleProxy& ele,
    LatParamProxy& param,
    bool err_flag);
extern "C" void fortran_ele_equal_ele(
    void* ele_out /* 0D_NOT_type */,
    void* ele_in /* 0D_NOT_type */);
EleProxy ele_equal_ele(EleProxy& ele_in);
extern "C" void fortran_ele_equals_ele(
    void* ele_out /* 0D_NOT_type */,
    void* ele_in /* 0D_NOT_type */,
    c_Bool& update_nametable /* 0D_NOT_logical */);
EleProxy ele_equals_ele(EleProxy& ele_in, bool update_nametable);
extern "C" void fortran_ele_finalizer(void* ele /* 0D_NOT_type */);
void ele_finalizer(EleProxy& ele);
extern "C" bool fortran_ele_full_name(
    void* ele /* 0D_NOT_type */,
    c_Char template_ /* 0D_NOT_character */,
    c_Char str /* 0D_ALLOC_character */);
void ele_full_name(
    EleProxy& ele,
    std::optional<std::string> template_,
    std::string str);
extern "C" void fortran_ele_geometry(
    void* floor_start /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* floor_end /* 0D_NOT_type */,
    c_Real* len_scale /* 0D_NOT_real */,
    c_Bool* ignore_patch_err /* 0D_NOT_logical */);
FloorPositionProxy ele_geometry(
    FloorPositionProxy& floor_start,
    EleProxy& ele,
    std::optional<double> len_scale = std::nullopt,
    std::optional<bool> ignore_patch_err = std::nullopt);

// Skipped unusable routine ele_geometry_hook_def:
// Routine in configuration skip list
extern "C" bool fortran_ele_geometry_with_misalignments(
    void* ele /* 0D_NOT_type */,
    c_Real* len_scale /* 0D_NOT_real */,
    void* floor /* 0D_NOT_type */);
void ele_geometry_with_misalignments(
    EleProxy& ele,
    std::optional<double> len_scale,
    FloorPositionProxy& floor);
extern "C" bool fortran_ele_has_constant_ds_dt_ref(
    void* ele /* 0D_NOT_type */,
    c_Bool& is_const /* 0D_NOT_logical */);
void ele_has_constant_ds_dt_ref(EleProxy& ele, bool is_const);
extern "C" bool fortran_ele_has_nonzero_kick(
    void* ele /* 0D_NOT_type */,
    c_Bool& has_kick /* 0D_NOT_logical */);
EleProxy ele_has_nonzero_kick(bool has_kick);
extern "C" bool fortran_ele_has_nonzero_offset(
    void* ele /* 0D_NOT_type */,
    c_Bool& has_offset /* 0D_NOT_logical */);
void ele_has_nonzero_offset(EleProxy& ele, bool has_offset);
extern "C" bool fortran_ele_is_monitor(
    void* ele /* 0D_NOT_type */,
    c_Bool* print_warning /* 0D_NOT_logical */,
    c_Bool& is_monitor /* 0D_NOT_logical */);
bool ele_is_monitor(
    EleProxy& ele,
    std::optional<bool> print_warning = std::nullopt);
extern "C" bool fortran_ele_loc(
    void* ele /* 0D_NOT_type */,
    void* loc /* 0D_NOT_type */);
void ele_loc(EleProxy& ele, LatEleLocProxy& loc);
extern "C" bool fortran_ele_loc_name(
    void* ele /* 0D_NOT_type */,
    c_Bool* show_branch0 /* 0D_NOT_logical */,
    c_Char parens /* 0D_NOT_character */,
    c_Char str /* 0D_NOT_character */);
void ele_loc_name(
    EleProxy& ele,
    std::optional<bool> show_branch0,
    std::optional<std::string> parens,
    std::string str);
extern "C" void fortran_ele_misalignment_l_s_calc(
    void* ele /* 0D_NOT_type */,
    c_RealArr L_mis /* 1D_NOT_real */,
    c_RealArr S_mis /* 2D_NOT_real */);
struct EleMisalignmentLSCalc {
  FixedArray1D<Real, 3> L_mis;
  FixedArray2D<Real, 3, 3> S_mis;
};
EleMisalignmentLSCalc ele_misalignment_l_s_calc(EleProxy& ele);
extern "C" bool fortran_ele_nametable_index(
    void* ele /* 0D_NOT_type */,
    c_Int& ix_nt /* 0D_NOT_integer */);
void ele_nametable_index(EleProxy& ele, int ix_nt);
extern "C" void fortran_ele_order_calc(
    void* lat /* 0D_NOT_type */,
    void* order /* 0D_NOT_type */);
LatEleOrderProxy ele_order_calc(LatProxy& lat);

// Skipped unusable routine ele_pointer_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ele_reference_energy_correction(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void ele_reference_energy_correction(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_ele_rf_step_index(
    c_Real& E_ref /* 0D_NOT_real */,
    c_Real& s_rel /* 0D_NOT_real */,
    void* ele /* 0D_NOT_type */,
    c_Int& ix_step /* 0D_NOT_integer */);
void ele_rf_step_index(double E_ref, double s_rel, EleProxy& ele, int ix_step);

// Skipped unusable routine ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ele_to_fibre:
// Untranslated type: FibreProxy (0D_PTR_type)

// Skipped unusable routine ele_to_fibre_hook_def:
// Untranslated type: FibreProxy (0D_NOT_type)
extern "C" void fortran_ele_to_ptc_magnetic_bn_an(
    void* ele /* 0D_NOT_type */,
    void* bn /* 1D_ALLOC_real */,
    void* an /* 1D_ALLOC_real */,
    c_Int& n_max /* 0D_NOT_integer */);
struct EleToPtcMagneticBnAn {
  RealAlloc1D bn;
  RealAlloc1D an;
  int n_max;
};
EleToPtcMagneticBnAn ele_to_ptc_magnetic_bn_an(EleProxy& ele);
extern "C" void fortran_ele_to_spin_taylor(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orb0 /* 0D_NOT_type */);
void ele_to_spin_taylor(EleProxy& ele, LatParamProxy& param, CoordProxy& orb0);
extern "C" void fortran_ele_to_taylor(
    void* ele /* 0D_NOT_type */,
    void* orb0 /* 0D_NOT_type */,
    c_Bool* taylor_map_includes_offsets /* 0D_NOT_logical */,
    c_Bool* include_damping /* 0D_NOT_logical */,
    void* orbital_taylor /* 1D_NOT_type */,
    void* spin_taylor /* 1D_NOT_type */);
struct EleToTaylor {
  TaylorProxyArray1D orbital_taylor;
  TaylorProxyArray1D spin_taylor;
};
EleToTaylor ele_to_taylor(
    EleProxy& ele,
    optional_ref<CoordProxy> orb0 = std::nullopt,
    std::optional<bool> taylor_map_includes_offsets = std::nullopt,
    std::optional<bool> include_damping = std::nullopt);
extern "C" bool fortran_ele_unique_name(
    void* ele /* 0D_NOT_type */,
    void* order /* 0D_NOT_type */,
    c_Char unique_name /* 0D_NOT_character */);
void ele_unique_name(
    EleProxy& ele,
    LatEleOrderProxy& order,
    std::string unique_name);
extern "C" bool fortran_ele_value_has_changed(
    void* ele /* 0D_NOT_type */,
    void* list /* 1D_ALLOC_integer */,
    void* abs_tol /* 1D_ALLOC_real */,
    c_Bool& set_old /* 0D_NOT_logical */,
    c_Bool& has_changed /* 0D_NOT_logical */);
void ele_value_has_changed(
    EleProxy& ele,
    IntAlloc1D& list,
    RealAlloc1D& abs_tol,
    bool set_old,
    bool has_changed);
extern "C" void fortran_ele_vec_equal_ele_vec(
    void* ele1 /* 1D_ALLOC_type */,
    void* ele2 /* 1D_ALLOC_type */);
EleProxyAlloc1D ele_vec_equal_ele_vec(EleProxyAlloc1D& ele2);
extern "C" void fortran_elec_multipole_field(
    c_Real& a /* 0D_NOT_real */,
    c_Real& b /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    void* coord /* 0D_NOT_type */,
    c_Real& Ex /* 0D_NOT_real */,
    c_Real& Ey /* 0D_NOT_real */,
    c_RealArr dE /* 2D_NOT_real */,
    c_Bool& compute_dE /* 0D_NOT_logical */);
struct ElecMultipoleField {
  double Ex;
  double Ey;
  std::optional<FixedArray2D<Real, 2, 2>> dE;
  bool compute_dE;
};
ElecMultipoleField elec_multipole_field(
    double a,
    double b,
    int n,
    CoordProxy& coord);
extern "C" void fortran_element_slice_iterator(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& i_slice /* 0D_NOT_integer */,
    c_Int& n_slice_tot /* 0D_NOT_integer */,
    void* sliced_ele /* 0D_NOT_type */,
    c_Real* s_start /* 0D_NOT_real */,
    c_Real* s_end /* 0D_NOT_real */);
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

// Skipped unusable routine ellipse_beam_init_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine em_field_calc:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine em_field_custom_def:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" void fortran_em_field_derivatives(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& s_pos /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& local_ref_frame /* 0D_NOT_logical */,
    void* dfield /* 0D_NOT_type */,
    c_Bool* grid_allow_s_out_of_bounds /* 0D_NOT_logical */,
    c_Real* rf_time /* 0D_NOT_real */);
EmFieldProxy em_field_derivatives(
    EleProxy& ele,
    LatParamProxy& param,
    double s_pos,
    CoordProxy& orbit,
    bool local_ref_frame,
    optional_ref<bool> grid_allow_s_out_of_bounds = std::nullopt,
    optional_ref<double> rf_time = std::nullopt);
extern "C" void fortran_em_field_kick_vector_time(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& rf_time /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr dvec_dt /* 1D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    void* extra_field /* 0D_NOT_type */);
FixedArray1D<Real, 10> em_field_kick_vector_time(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orbit,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" bool fortran_em_field_plus_em_field(
    void* field1 /* 0D_NOT_type */,
    void* field2 /* 0D_NOT_type */,
    void* field_tot /* 0D_NOT_type */);
EmFieldProxy em_field_plus_em_field(EmFieldProxy& field1, EmFieldProxy& field2);

// Skipped unusable routine em_field_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_em_taylor_equal_em_taylor(
    void* em_taylor1 /* 0D_NOT_type */,
    void* em_taylor2 /* 0D_NOT_type */);
EmTaylorProxy em_taylor_equal_em_taylor(EmTaylorProxy& em_taylor2);

// Skipped unusable routine em_taylor_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine em_taylor_term_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_em_taylors_equal_em_taylors(
    void* em_taylor1 /* 1D_ALLOC_type */,
    void* em_taylor2 /* 1D_ALLOC_type */);
EmTaylorProxyAlloc1D em_taylors_equal_em_taylors(
    EmTaylorProxyAlloc1D& em_taylor2);
extern "C" void fortran_emit_6d(
    void* ele_ref /* 0D_NOT_type */,
    c_Bool& include_opening_angle /* 0D_NOT_logical */,
    void* mode /* 0D_NOT_type */,
    c_RealArr sigma_mat /* 2D_NOT_real */,
    void* closed_orbit /* 1D_ALLOC_type */,
    void* rad_int_by_ele /* 0D_NOT_type */);
struct Emit6d {
  NormalModesProxy mode;
  FixedArray2D<Real, 6, 6> sigma_mat;
  RadIntAllEleProxy rad_int_by_ele;
};
Emit6d emit_6d(
    EleProxy& ele_ref,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit = std::nullopt);
extern "C" bool fortran_entering_element(
    void* orbit /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    c_Bool& is_entering /* 0D_NOT_logical */);
void entering_element(CoordProxy& orbit, int particle_at, bool is_entering);
extern "C" void fortran_envelope_radints(
    c_ComplexArr Lambda /* 2D_NOT_complex */,
    c_ComplexArr Theta /* 2D_NOT_complex */,
    c_ComplexArr Iota /* 2D_NOT_complex */,
    c_RealArr alpha /* 1D_NOT_real */,
    c_RealArr emit /* 1D_NOT_real */);
void envelope_radints(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    FixedArray1D<Real, 3> alpha,
    FixedArray1D<Real, 3> emit);
extern "C" void fortran_envelope_radints_ibs(
    c_ComplexArr Lambda /* 2D_NOT_complex */,
    c_ComplexArr Theta /* 2D_NOT_complex */,
    c_ComplexArr Iota /* 2D_NOT_complex */,
    void* eles /* 1D_ALLOC_type */,
    c_RealArr alpha /* 1D_NOT_real */,
    c_RealArr emit /* 1D_NOT_real */,
    void* mode /* 0D_NOT_type */,
    c_Bool& tail_cut /* 0D_NOT_logical */,
    c_Real& npart /* 0D_NOT_real */,
    c_Int& species /* 0D_NOT_integer */);
struct EnvelopeRadintsIbs {
  FixedArray1D<Real, 3> alpha;
  FixedArray1D<Real, 3> emit;
};
EnvelopeRadintsIbs envelope_radints_ibs(
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    EleProxyAlloc1D& eles,
    NormalModesProxy& mode,
    bool tail_cut,
    double npart,
    int species);
extern "C" bool fortran_eq_ac_kicker(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ac_kicker(AcKickerProxy& f1, AcKickerProxy& f2, bool is_eq);
extern "C" bool fortran_eq_ac_kicker_freq(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ac_kicker_freq(
    AcKickerFreqProxy& f1,
    AcKickerFreqProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_ac_kicker_time(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ac_kicker_time(
    AcKickerTimeProxy& f1,
    AcKickerTimeProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_anormal_mode(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_anormal_mode(AnormalModeProxy& f1, AnormalModeProxy& f2, bool is_eq);
extern "C" bool fortran_eq_aperture_param(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_aperture_param(
    ApertureParamProxy& f1,
    ApertureParamProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_aperture_point(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_aperture_point(
    AperturePointProxy& f1,
    AperturePointProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_aperture_scan(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_aperture_scan(ApertureScanProxy& f1, ApertureScanProxy& f2, bool is_eq);
extern "C" bool fortran_eq_beam(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_beam(BeamProxy& f1, BeamProxy& f2, bool is_eq);
extern "C" bool fortran_eq_beam_init(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_beam_init(BeamInitProxy& f1, BeamInitProxy& f2, bool is_eq);
extern "C" bool fortran_eq_bmad_common(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_bmad_common(BmadCommonProxy& f1, BmadCommonProxy& f2, bool is_eq);
extern "C" bool fortran_eq_bookkeeping_state(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_bookkeeping_state(
    BookkeepingStateProxy& f1,
    BookkeepingStateProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_bpm_phase_coupling(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_bpm_phase_coupling(
    BpmPhaseCouplingProxy& f1,
    BpmPhaseCouplingProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_branch(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_branch(BranchProxy& f1, BranchProxy& f2, bool is_eq);
extern "C" bool fortran_eq_bunch(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_bunch(BunchProxy& f1, BunchProxy& f2, bool is_eq);
extern "C" bool fortran_eq_bunch_params(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_bunch_params(BunchParamsProxy& f1, BunchParamsProxy& f2, bool is_eq);
extern "C" bool fortran_eq_cartesian_map(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cartesian_map(CartesianMapProxy& f1, CartesianMapProxy& f2, bool is_eq);
extern "C" bool fortran_eq_cartesian_map_term(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cartesian_map_term(
    CartesianMapTermProxy& f1,
    CartesianMapTermProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_cartesian_map_term1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cartesian_map_term1(
    CartesianMapTerm1Proxy& f1,
    CartesianMapTerm1Proxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_complex_taylor(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_complex_taylor(
    ComplexTaylorProxy& f1,
    ComplexTaylorProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_complex_taylor_term(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_complex_taylor_term(
    ComplexTaylorTermProxy& f1,
    ComplexTaylorTermProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_control(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_control(ControlProxy& f1, ControlProxy& f2, bool is_eq);
extern "C" bool fortran_eq_control_ramp1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_control_ramp1(ControlRamp1Proxy& f1, ControlRamp1Proxy& f2, bool is_eq);
extern "C" bool fortran_eq_control_var1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_control_var1(ControlVar1Proxy& f1, ControlVar1Proxy& f2, bool is_eq);
extern "C" bool fortran_eq_controller(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_controller(ControllerProxy& f1, ControllerProxy& f2, bool is_eq);
extern "C" bool fortran_eq_coord(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_coord(CoordProxy& f1, CoordProxy& f2, bool is_eq);
extern "C" bool fortran_eq_coord_array(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_coord_array(CoordArrayProxy& f1, CoordArrayProxy& f2, bool is_eq);
extern "C" bool fortran_eq_cylindrical_map(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cylindrical_map(
    CylindricalMapProxy& f1,
    CylindricalMapProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_cylindrical_map_term(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cylindrical_map_term(
    CylindricalMapTermProxy& f1,
    CylindricalMapTermProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_cylindrical_map_term1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_cylindrical_map_term1(
    CylindricalMapTerm1Proxy& f1,
    CylindricalMapTerm1Proxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_ele(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ele(EleProxy& f1, EleProxy& f2, bool is_eq);
extern "C" bool fortran_eq_ellipse_beam_init(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ellipse_beam_init(
    EllipseBeamInitProxy& f1,
    EllipseBeamInitProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_em_field(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_em_field(EmFieldProxy& f1, EmFieldProxy& f2, bool is_eq);
extern "C" bool fortran_eq_em_taylor(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_em_taylor(EmTaylorProxy& f1, EmTaylorProxy& f2, bool is_eq);
extern "C" bool fortran_eq_em_taylor_term(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_em_taylor_term(
    EmTaylorTermProxy& f1,
    EmTaylorTermProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_expression_atom(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_expression_atom(
    ExpressionAtomProxy& f1,
    ExpressionAtomProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_floor_position(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_floor_position(
    FloorPositionProxy& f1,
    FloorPositionProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_gen_grad1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_gen_grad1(GenGrad1Proxy& f1, GenGrad1Proxy& f2, bool is_eq);
extern "C" bool fortran_eq_gen_grad_map(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_gen_grad_map(GenGradMapProxy& f1, GenGradMapProxy& f2, bool is_eq);
extern "C" bool fortran_eq_grid_beam_init(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_grid_beam_init(
    GridBeamInitProxy& f1,
    GridBeamInitProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_grid_field(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_grid_field(GridFieldProxy& f1, GridFieldProxy& f2, bool is_eq);
extern "C" bool fortran_eq_grid_field_pt(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_grid_field_pt(GridFieldPtProxy& f1, GridFieldPtProxy& f2, bool is_eq);
extern "C" bool fortran_eq_grid_field_pt1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_grid_field_pt1(
    GridFieldPt1Proxy& f1,
    GridFieldPt1Proxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_high_energy_space_charge(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_high_energy_space_charge(
    HighEnergySpaceChargeProxy& f1,
    HighEnergySpaceChargeProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_interval1_coef(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_interval1_coef(
    Interval1CoefProxy& f1,
    Interval1CoefProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_kv_beam_init(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_kv_beam_init(KvBeamInitProxy& f1, KvBeamInitProxy& f2, bool is_eq);
extern "C" bool fortran_eq_lat(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_lat(LatProxy& f1, LatProxy& f2, bool is_eq);
extern "C" bool fortran_eq_lat_ele_loc(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_lat_ele_loc(LatEleLocProxy& f1, LatEleLocProxy& f2, bool is_eq);
extern "C" bool fortran_eq_lat_param(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_lat_param(LatParamProxy& f1, LatParamProxy& f2, bool is_eq);
extern "C" bool fortran_eq_linac_normal_mode(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_linac_normal_mode(
    LinacNormalModeProxy& f1,
    LinacNormalModeProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_mode3(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_mode3(Mode3Proxy& f1, Mode3Proxy& f2, bool is_eq);
extern "C" bool fortran_eq_mode_info(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_mode_info(ModeInfoProxy& f1, ModeInfoProxy& f2, bool is_eq);
extern "C" bool fortran_eq_normal_modes(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_normal_modes(NormalModesProxy& f1, NormalModesProxy& f2, bool is_eq);
extern "C" bool fortran_eq_photon_element(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_photon_element(
    PhotonElementProxy& f1,
    PhotonElementProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_photon_material(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_photon_material(
    PhotonMaterialProxy& f1,
    PhotonMaterialProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_photon_reflect_surface(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_photon_reflect_surface(
    PhotonReflectSurfaceProxy& f1,
    PhotonReflectSurfaceProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_photon_reflect_table(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_photon_reflect_table(
    PhotonReflectTableProxy& f1,
    PhotonReflectTableProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_photon_target(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_photon_target(PhotonTargetProxy& f1, PhotonTargetProxy& f2, bool is_eq);
extern "C" bool fortran_eq_pixel_detec(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_pixel_detec(PixelDetecProxy& f1, PixelDetecProxy& f2, bool is_eq);
extern "C" bool fortran_eq_pixel_pt(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_pixel_pt(PixelPtProxy& f1, PixelPtProxy& f2, bool is_eq);
extern "C" bool fortran_eq_pre_tracker(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_pre_tracker(PreTrackerProxy& f1, PreTrackerProxy& f2, bool is_eq);
extern "C" bool fortran_eq_rad_int1(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_rad_int1(RadInt1Proxy& f1, RadInt1Proxy& f2, bool is_eq);
extern "C" bool fortran_eq_rad_int_all_ele(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_rad_int_all_ele(
    RadIntAllEleProxy& f1,
    RadIntAllEleProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_rad_int_branch(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_rad_int_branch(
    RadIntBranchProxy& f1,
    RadIntBranchProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_rad_map(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_rad_map(RadMapProxy& f1, RadMapProxy& f2, bool is_eq);
extern "C" bool fortran_eq_rad_map_ele(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_rad_map_ele(RadMapEleProxy& f1, RadMapEleProxy& f2, bool is_eq);
extern "C" bool fortran_eq_ramper_lord(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_ramper_lord(RamperLordProxy& f1, RamperLordProxy& f2, bool is_eq);
extern "C" bool fortran_eq_space_charge_common(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_space_charge_common(
    SpaceChargeCommonProxy& f1,
    SpaceChargeCommonProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_spin_polar(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_spin_polar(SpinPolarProxy& f1, SpinPolarProxy& f2, bool is_eq);
extern "C" bool fortran_eq_spline(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_spline(SplineProxy& f1, SplineProxy& f2, bool is_eq);
extern "C" bool fortran_eq_strong_beam(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_strong_beam(StrongBeamProxy& f1, StrongBeamProxy& f2, bool is_eq);
extern "C" bool fortran_eq_surface_curvature(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_curvature(
    SurfaceCurvatureProxy& f1,
    SurfaceCurvatureProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_displacement(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_displacement(
    SurfaceDisplacementProxy& f1,
    SurfaceDisplacementProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_displacement_pt(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_displacement_pt(
    SurfaceDisplacementPtProxy& f1,
    SurfaceDisplacementPtProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_h_misalign(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_h_misalign(
    SurfaceHMisalignProxy& f1,
    SurfaceHMisalignProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_h_misalign_pt(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_h_misalign_pt(
    SurfaceHMisalignPtProxy& f1,
    SurfaceHMisalignPtProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_segmented(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_segmented(
    SurfaceSegmentedProxy& f1,
    SurfaceSegmentedProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_surface_segmented_pt(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_surface_segmented_pt(
    SurfaceSegmentedPtProxy& f1,
    SurfaceSegmentedPtProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_target_point(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_target_point(TargetPointProxy& f1, TargetPointProxy& f2, bool is_eq);
extern "C" bool fortran_eq_taylor(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_taylor(TaylorProxy& f1, TaylorProxy& f2, bool is_eq);
extern "C" bool fortran_eq_taylor_term(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_taylor_term(TaylorTermProxy& f1, TaylorTermProxy& f2, bool is_eq);
extern "C" bool fortran_eq_track(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_track(TrackProxy& f1, TrackProxy& f2, bool is_eq);
extern "C" bool fortran_eq_track_point(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_track_point(TrackPointProxy& f1, TrackPointProxy& f2, bool is_eq);
extern "C" bool fortran_eq_twiss(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_twiss(TwissProxy& f1, TwissProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake(WakeProxy& f1, WakeProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake_lr(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake_lr(WakeLrProxy& f1, WakeLrProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake_lr_mode(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake_lr_mode(WakeLrModeProxy& f1, WakeLrModeProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake_sr(WakeSrProxy& f1, WakeSrProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr_mode(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake_sr_mode(WakeSrModeProxy& f1, WakeSrModeProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wake_sr_z_long(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wake_sr_z_long(WakeSrZLongProxy& f1, WakeSrZLongProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wall3d(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wall3d(Wall3dProxy& f1, Wall3dProxy& f2, bool is_eq);
extern "C" bool fortran_eq_wall3d_section(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wall3d_section(
    Wall3dSectionProxy& f1,
    Wall3dSectionProxy& f2,
    bool is_eq);
extern "C" bool fortran_eq_wall3d_vertex(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_wall3d_vertex(Wall3dVertexProxy& f1, Wall3dVertexProxy& f2, bool is_eq);
extern "C" bool fortran_eq_xy_disp(
    void* f1 /* 0D_NOT_type */,
    void* f2 /* 0D_NOT_type */,
    c_Bool& is_eq /* 0D_NOT_logical */);
void eq_xy_disp(XyDispProxy& f1, XyDispProxy& f2, bool is_eq);
extern "C" bool fortran_equal_sign_here(
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& is_here /* 0D_NOT_logical */);
void equal_sign_here(EleProxy& ele, std::string delim, bool is_here);
extern "C" bool fortran_equivalent_taylor_attributes(
    void* ele_taylor /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& equiv /* 0D_NOT_logical */);
void equivalent_taylor_attributes(
    EleProxy& ele_taylor,
    EleProxy& ele2,
    bool equiv);
extern "C" void fortran_etdiv(
    c_Real& A /* 0D_NOT_real */,
    c_Real& B /* 0D_NOT_real */,
    c_Real& C /* 0D_NOT_real */,
    c_Real& D /* 0D_NOT_real */,
    c_Real& E /* 0D_NOT_real */,
    c_Real& F /* 0D_NOT_real */);
void etdiv(double A, double B, double C, double D, double E, double F);

// Skipped unusable routine ety:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine ety2:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine etyt:
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_evaluate_array_index(
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Char delim_list1 /* 0D_NOT_character */,
    c_Char word2 /* 0D_NOT_character */,
    c_Char delim_list2 /* 0D_NOT_character */,
    c_Char delim2 /* 0D_NOT_character */,
    c_Int& this_index /* 0D_NOT_integer */);
struct EvaluateArrayIndex {
  bool err_flag;
  std::string word2;
  std::string delim2;
  int this_index;
};
EvaluateArrayIndex evaluate_array_index(
    std::string delim_list1,
    std::string delim_list2);
extern "C" bool fortran_evaluate_logical(
    c_Char word /* 0D_NOT_character */,
    c_Int& iostat /* 0D_NOT_integer */,
    c_Bool& this_logic /* 0D_NOT_logical */);
struct EvaluateLogical {
  int iostat;
  bool this_logic;
};
EvaluateLogical evaluate_logical(std::string word);
extern "C" void fortran_exact_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void exact_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_exp_bessi0(
    c_Real& t /* 0D_NOT_real */,
    c_Real& B1 /* 0D_NOT_real */,
    c_Real& B2 /* 0D_NOT_real */,
    c_Real& func_retval__ /* 0D_NOT_real */);
double exp_bessi0(double t, double B1, double B2);
extern "C" bool fortran_expect_one_of(
    c_Char delim_list /* 0D_NOT_character */,
    c_Bool& check_input_delim /* 0D_NOT_logical */,
    c_Char ele_name /* 0D_NOT_character */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& is_ok /* 0D_NOT_logical */);
void expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool is_ok);
extern "C" bool fortran_expect_this(
    c_Char expecting /* 0D_NOT_character */,
    c_Bool& check_delim /* 0D_NOT_logical */,
    c_Bool& call_check /* 0D_NOT_logical */,
    c_Char err_str /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ExpectThis {
  std::string delim;
  bool delim_found;
  bool is_ok;
};
ExpectThis expect_this(
    std::string expecting,
    bool check_delim,
    bool call_check,
    std::string err_str,
    EleProxy& ele);

// Skipped unusable routine expression_atom_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_expression_stack_to_string(
    void* stack /* 1D_ALLOC_type */,
    c_Bool* polish /* 0D_NOT_logical */,
    c_Char str /* 0D_ALLOC_character */);
std::string expression_stack_to_string(
    ExpressionAtomProxyAlloc1D& stack,
    std::optional<bool> polish = std::nullopt);
extern "C" bool fortran_expression_stack_value(
    void* stack /* 1D_ALLOC_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Char err_str /* 0D_NOT_character */,
    void* var /* 1D_ALLOC_type */,
    c_Bool* use_old /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
struct ExpressionStackValue {
  bool err_flag;
  std::string err_str;
  double value;
};
ExpressionStackValue expression_stack_value(
    ExpressionAtomProxyAlloc1D& stack,
    optional_ref<ControlVar1ProxyAlloc1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt);
extern "C" void fortran_expression_string_to_stack(
    c_Char string /* 0D_NOT_character */,
    void* stack /* 1D_ALLOC_type */,
    c_Int& n_stack /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Char err_str /* 0D_NOT_character */);
struct ExpressionStringToStack {
  ExpressionAtomProxyAlloc1D stack;
  int n_stack;
  bool err_flag;
  std::string err_str;
};
ExpressionStringToStack expression_string_to_stack(std::string string);

// Skipped unusable routine expression_string_to_tree:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)

// Skipped unusable routine expression_tree_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine expression_tree_to_string:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)
extern "C" bool fortran_expression_value(
    c_Char expression /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Char err_str /* 0D_NOT_character */,
    void* var /* 1D_ALLOC_type */,
    c_Bool* use_old /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
struct ExpressionValue {
  bool err_flag;
  std::string err_str;
  double value;
};
ExpressionValue expression_value(
    std::string expression,
    optional_ref<ControlVar1ProxyAlloc1D> var = std::nullopt,
    std::optional<bool> use_old = std::nullopt);

// Skipped unusable routine extra_parsing_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_fft1(
    void* a /* 1D_ALLOC_real */,
    void* b /* 1D_ALLOC_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Int& isn /* 0D_NOT_integer */,
    c_Int& ierr /* 0D_NOT_integer */);
void fft1(RealAlloc1D& a, RealAlloc1D& b, int n, int isn, int ierr);

// Skipped unusable routine fftconvcorr3d:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine fibre_to_ele:
// Untranslated type: FibreProxy (0D_NOT_type)
extern "C" bool fortran_field_attribute_free(
    void* ele /* 0D_NOT_type */,
    c_Char attrib_name /* 0D_NOT_character */,
    c_Bool& free /* 0D_NOT_logical */);
bool field_attribute_free(EleProxy& ele, std::string attrib_name);

// Skipped unusable routine field_interpolate_3d:
// Variable in sized array: field_mesh(0:,0:,0:) 3D_NOT_real
extern "C" void fortran_finalize_reflectivity_table(
    void* table /* 0D_NOT_type */,
    c_Bool& in_degrees /* 0D_NOT_logical */);
void finalize_reflectivity_table(
    PhotonReflectTableProxy& table,
    bool in_degrees);
extern "C" void fortran_find_element_ends(
    void* ele /* 0D_NOT_type */,
    void* ele1 /* 0D_PTR_type */,
    void* ele2 /* 0D_PTR_type */,
    c_Int* ix_multipass /* 0D_NOT_integer */);
struct FindElementEnds {
  EleProxy ele1;
  EleProxy ele2;
};
FindElementEnds find_element_ends(
    EleProxy& ele,
    std::optional<int> ix_multipass = std::nullopt);
extern "C" void fortran_find_fwhm(
    c_Real& bound /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    c_Real& fwhm /* 0D_NOT_real */);
double find_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_find_matching_fieldmap(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Int& fm_type /* 0D_NOT_integer */,
    void* match_ele /* 0D_PTR_type */,
    c_Int& ix_field /* 0D_NOT_integer */,
    c_Bool* ignore_slaves /* 0D_NOT_logical */);
struct FindMatchingFieldmap {
  EleProxy match_ele;
  int ix_field;
};
FindMatchingFieldmap find_matching_fieldmap(
    std::string file_name,
    EleProxy& ele,
    int fm_type,
    std::optional<bool> ignore_slaves = std::nullopt);
extern "C" void fortran_find_normalization(
    c_Real& bound /* 0D_NOT_real */,
    c_Real& p0 /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    c_Real& pnrml /* 0D_NOT_real */);
double find_normalization(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_floor_angles_to_w_mat(
    c_Real& theta /* 0D_NOT_real */,
    c_Real& phi /* 0D_NOT_real */,
    c_Real& psi /* 0D_NOT_real */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_RealArr w_mat_inv /* 2D_NOT_real */);
struct FloorAnglesToWMat {
  std::optional<FixedArray2D<Real, 3, 3>> w_mat;
  std::optional<FixedArray2D<Real, 3, 3>> w_mat_inv;
};
FloorAnglesToWMat floor_angles_to_w_mat(double theta, double phi, double psi);

// Skipped unusable routine floor_position_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_floor_w_mat_to_angles(
    c_RealArr w_mat /* 2D_NOT_real */,
    c_Real& theta /* 0D_NOT_real */,
    c_Real& phi /* 0D_NOT_real */,
    c_Real& psi /* 0D_NOT_real */,
    void* floor0 /* 0D_NOT_type */);
struct FloorWMatToAngles {
  double theta;
  double phi;
  double psi;
};
FloorWMatToAngles floor_w_mat_to_angles(
    FixedArray2D<Real, 3, 3> w_mat,
    optional_ref<FloorPositionProxy> floor0 = std::nullopt);

// Skipped unusable routine foil_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_form_complex_taylor(
    void* re_taylor /* 0D_NOT_type */,
    void* im_taylor /* 0D_NOT_type */,
    void* complex_taylor /* 0D_NOT_type */);
ComplexTaylorProxy form_complex_taylor(
    TaylorProxy& re_taylor,
    TaylorProxy& im_taylor);
extern "C" void fortran_form_digested_bmad_file_name(
    c_Char lat_file /* 0D_NOT_character */,
    c_Char digested_file /* 0D_NOT_character */,
    c_Char full_lat_file /* 0D_NOT_character */,
    c_Char use_line /* 0D_NOT_character */);
struct FormDigestedBmadFileName {
  std::string digested_file;
  std::string full_lat_file;
};
FormDigestedBmadFileName form_digested_bmad_file_name(
    std::string lat_file,
    std::optional<std::string> use_line = std::nullopt);

// Skipped unusable routine fringe_field_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_fringe_here(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    c_Bool& is_here /* 0D_NOT_logical */);
void fringe_here(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    bool is_here);
extern "C" bool fortran_g_bend_from_em_field(
    c_RealArr b /* 1D_NOT_real */,
    c_RealArr e /* 1D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr g_bend /* 1D_NOT_real */);
FixedArray1D<Real, 3> g_bend_from_em_field(
    FixedArray1D<Real, 3> b,
    FixedArray1D<Real, 3> e,
    CoordProxy& orbit);
extern "C" void fortran_g_bending_strength_from_em_field(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& s_rel /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& local_ref_frame /* 0D_NOT_logical */,
    c_RealArr g /* 1D_NOT_real */,
    c_RealArr dg /* 2D_NOT_real */);
struct GBendingStrengthFromEmField {
  FixedArray1D<Real, 3> g;
  std::optional<FixedArray2D<Real, 3, 3>> dg;
};
GBendingStrengthFromEmField g_bending_strength_from_em_field(
    EleProxy& ele,
    LatParamProxy& param,
    double s_rel,
    CoordProxy& orbit,
    bool local_ref_frame);
extern "C" void fortran_g_integrals_calc(void* lat /* 0D_NOT_type */);
void g_integrals_calc(LatProxy& lat);
extern "C" bool fortran_gamma_ref(
    void* ele /* 0D_NOT_type */,
    c_Real& gamma /* 0D_NOT_real */);
void gamma_ref(EleProxy& ele, double gamma);

// Skipped unusable routine gen_grad1_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_gen_grad1_to_em_taylor(
    void* ele /* 0D_NOT_type */,
    void* gen_grad /* 0D_NOT_type */,
    c_Int& iz /* 0D_NOT_integer */,
    void* em_taylor /* 1D_NOT_type */);
EmTaylorProxyArray1D gen_grad1_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    int iz);
extern "C" void fortran_gen_grad_at_s_to_em_taylor(
    void* ele /* 0D_NOT_type */,
    void* gen_grad /* 0D_NOT_type */,
    c_Real& s_pos /* 0D_NOT_real */,
    void* em_taylor /* 1D_NOT_type */);
EmTaylorProxyArray1D gen_grad_at_s_to_em_taylor(
    EleProxy& ele,
    GenGradMapProxy& gen_grad,
    double s_pos);
extern "C" bool fortran_gen_grad_field(
    void* deriv /* 1D_ALLOC_real */,
    void* gg /* 0D_NOT_type */,
    c_Real& rho /* 0D_NOT_real */,
    c_Real& theta /* 0D_NOT_real */,
    c_RealArr field /* 1D_NOT_real */);
void gen_grad_field(
    RealAlloc1D& deriv,
    GenGrad1Proxy& gg,
    double rho,
    double theta,
    FixedArray1D<Real, 3> field);

// Skipped unusable routine gen_grad_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine get_astra_fieldgrid_name_and_scaling:
// Untranslated type: StrIndexProxy (0D_NOT_type)
extern "C" void fortran_get_bl_from_fwhm(
    c_Real& bound /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    c_Real& sigma /* 0D_NOT_real */);
double get_bl_from_fwhm(double bound, FixedArray1D<Real, 8> args);
extern "C" void fortran_get_called_file(
    c_Char delim /* 0D_NOT_character */,
    c_Char call_file /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
void get_called_file(std::string delim, std::string call_file, bool err);

// Skipped unusable routine get_cgrn_csr3d:
// Variable inout sized array: cgrn(0:,0:,0:) 3D_NOT_complex
extern "C" void fortran_get_emit_from_sigma_mat(
    c_RealArr sigma_mat /* 2D_NOT_real */,
    c_RealArr normal /* 1D_NOT_real */,
    c_RealArr Nmat /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct GetEmitFromSigmaMat {
  FixedArray1D<Real, 3> normal;
  bool err_flag;
};
GetEmitFromSigmaMat get_emit_from_sigma_mat(
    FixedArray2D<Real, 6, 6> sigma_mat,
    std::optional<FixedArray2D<Real, 6, 6>> Nmat = std::nullopt);

// Skipped unusable routine get_gpt_fieldgrid_name_and_scaling:
// Untranslated type: StrIndexProxy (0D_NOT_type)

// Skipped unusable routine get_list_of_names:
// Variable-sized inout character array: name_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_get_next_word(
    c_Char word /* 0D_NOT_character */,
    c_Int& ix_word /* 0D_NOT_integer */,
    c_Char delim_list /* 0D_NOT_character */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool* upper_case_word /* 0D_NOT_logical */,
    c_Bool* call_check /* 0D_NOT_logical */,
    c_Bool* err_flag /* 0D_NOT_logical */);
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
// Untranslated type: StrIndexProxy (0D_NOT_type)

// Skipped unusable routine get_overlay_group_names:
// Untranslated type: ParserEleProxy (0D_NOT_type)
// Variable-sized inout character array: names_out(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine get_sequence_args:
// Variable-sized inout character array: arg_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine get_slave_list:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine get_switch:
// Variable-sized inout character array: name_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine getrhotilde:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_gpt_field_grid_scaling(
    void* ele /* 0D_NOT_type */,
    c_Int& dimensions /* 0D_NOT_integer */,
    c_Real& field_scale /* 0D_NOT_real */,
    c_Real& ref_time /* 0D_NOT_real */);
void gpt_field_grid_scaling(
    EleProxy& ele,
    int dimensions,
    double field_scale,
    double ref_time);

// Skipped unusable routine gpt_lat_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_gpt_max_field_reference(
    void* pt0 /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& field_value /* 0D_NOT_real */);
void gpt_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double field_value);
extern "C" void fortran_gpt_to_particle_bunch(
    c_Char gpt_file /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* bunch /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct GptToParticleBunch {
  BunchProxy bunch;
  bool err_flag;
};
GptToParticleBunch gpt_to_particle_bunch(std::string gpt_file, EleProxy& ele);
extern "C" bool fortran_gradient_shift_sr_wake(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& grad_shift /* 0D_NOT_real */);
void gradient_shift_sr_wake(
    EleProxy& ele,
    LatParamProxy& param,
    double grad_shift);

// Skipped unusable routine grid_beam_init_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_grid_field_interpolate(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    void* grid /* 0D_NOT_type */,
    void* g_field /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& x1 /* 0D_NOT_real */,
    c_Real* x2 /* 0D_NOT_real */,
    c_Real* x3 /* 0D_NOT_real */,
    c_Bool* allow_s_out_of_bounds /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
void grid_field_interpolate(
    EleProxy& ele,
    CoordProxy& orbit,
    GridFieldProxy& grid,
    GridFieldPt1Proxy& g_field,
    bool err_flag,
    double x1,
    std::optional<double> x2 = std::nullopt,
    std::optional<double> x3 = std::nullopt,
    std::optional<bool> allow_s_out_of_bounds = std::nullopt,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine grid_field_pt1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine grid_field_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine grid_field_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_hard_multipole_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void hard_multipole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_has_attribute(
    void* ele /* 0D_NOT_type */,
    c_Char attrib /* 0D_NOT_character */,
    c_Bool& has_it /* 0D_NOT_logical */);
void has_attribute(EleProxy& ele, std::string attrib, bool has_it);
extern "C" bool fortran_has_curvature(
    void* phot_ele /* 0D_NOT_type */,
    c_Bool& curved /* 0D_NOT_logical */);
bool has_curvature(PhotonElementProxy& phot_ele);
extern "C" bool fortran_has_orientation_attributes(
    void* ele /* 0D_NOT_type */,
    c_Bool& has_attribs /* 0D_NOT_logical */);
bool has_orientation_attributes(EleProxy& ele);

// Skipped unusable routine hdf5_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine hdf5_read_beam:
// Untranslated type: PmdHeaderProxy (0D_NOT_type)

// Skipped unusable routine hdf5_read_grid_field:
// Untranslated type: PmdHeaderProxy (0D_NOT_type)
extern "C" void fortran_hdf5_write_beam(
    c_Char file_name /* 0D_NOT_character */,
    void* bunches /* 1D_ALLOC_type */,
    c_Bool& append /* 0D_NOT_logical */,
    c_Bool& error /* 0D_NOT_logical */,
    void* lat /* 0D_NOT_type */,
    c_Bool* alive_only /* 0D_NOT_logical */);
void hdf5_write_beam(
    std::string file_name,
    BunchProxyAlloc1D& bunches,
    bool append,
    bool error,
    optional_ref<LatProxy> lat = std::nullopt,
    optional_ref<bool> alive_only = std::nullopt);
extern "C" void fortran_hdf5_write_grid_field(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* g_field /* 1D_ALLOC_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void hdf5_write_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxyAlloc1D& g_field,
    bool err_flag);

// Skipped unusable routine high_energy_space_charge_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_hwang_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void hwang_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);

// Skipped unusable routine i_csr:
// Untranslated type: CsrKick1Proxy (0D_NOT_type)
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine ibs1:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine ibs_blowup1turn:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)

// Skipped unusable routine ibs_delta_calc:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)

// Skipped unusable routine ibs_equib_der:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)

// Skipped unusable routine ibs_equib_rlx:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)

// Skipped unusable routine ibs_lifetime:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)
// Untranslated type: IbsMaxratioProxy (0D_NOT_type)
// Untranslated type: IbsLifetimeProxy (0D_NOT_type)

// Skipped unusable routine ibs_lifetime_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_ibs_matrix_c(
    c_RealArr sigma_mat /* 2D_NOT_real */,
    c_Bool& tail_cut /* 0D_NOT_logical */,
    c_Real& tau /* 0D_NOT_real */,
    c_Real& energy /* 0D_NOT_real */,
    c_Real& n_part /* 0D_NOT_real */,
    c_Int& species /* 0D_NOT_integer */,
    c_RealArr ibs_mat /* 2D_NOT_real */);
void ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species,
    FixedArray2D<Real, 6, 6> ibs_mat);

// Skipped unusable routine ibs_maxratio_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ibs_rates1turn:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine ibs_sim_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ibs_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_igfcoulombfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& dx /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */,
    c_Real& dz /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void igfcoulombfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res);
extern "C" bool fortran_igfexfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& dx /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */,
    c_Real& dz /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void igfexfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res);
extern "C" bool fortran_igfeyfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& dx /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */,
    c_Real& dz /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void igfeyfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res);
extern "C" bool fortran_igfezfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& dx /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */,
    c_Real& dz /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void igfezfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res);

// Skipped unusable routine image_charge_kick_calc:
// Untranslated type: CsrKick1Proxy (0D_NOT_type)
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine imageconvcorr3d:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_init_attribute_name1(
    c_Int& ix_key /* 0D_NOT_integer */,
    c_Int& ix_attrib /* 0D_NOT_integer */,
    c_Char name /* 0D_NOT_character */,
    c_Int* attrib_state /* 0D_NOT_integer */,
    c_Bool* override /* 0D_NOT_logical */);
void init_attribute_name1(
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state = std::nullopt,
    std::optional<bool> override = std::nullopt);
extern "C" void fortran_init_attribute_name_array();
void init_attribute_name_array();
extern "C" void fortran_init_beam_distribution(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* beam_init /* 0D_NOT_type */,
    void* beam /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* modes /* 0D_NOT_type */,
    void* beam_init_set /* 0D_NOT_type */,
    c_Bool* print_p0c_shift_warning /* 0D_NOT_logical */,
    c_Bool* conserve_momentum /* 0D_NOT_logical */);
struct InitBeamDistribution {
  BeamProxy beam;
  bool err_flag;
  BeamInitProxy beam_init_set;
};
InitBeamDistribution init_beam_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_init_bmad();
void init_bmad();
extern "C" void fortran_init_bmad_parser_common(void* lat /* 0D_NOT_type */);
void init_bmad_parser_common(optional_ref<LatProxy> lat = std::nullopt);
extern "C" void fortran_init_bunch_distribution(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* beam_init /* 0D_NOT_type */,
    c_Int& ix_bunch /* 0D_NOT_integer */,
    void* bunch /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* modes /* 0D_NOT_type */,
    void* beam_init_used /* 0D_NOT_type */,
    c_Bool* print_p0c_shift_warning /* 0D_NOT_logical */,
    c_Bool* conserve_momentum /* 0D_NOT_logical */);
struct InitBunchDistribution {
  BunchProxy bunch;
  bool err_flag;
  BeamInitProxy beam_init_used;
};
InitBunchDistribution init_bunch_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    int ix_bunch,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_init_complex_taylor_series(
    void* complex_taylor /* 0D_NOT_type */,
    c_Int& n_term /* 0D_NOT_integer */,
    c_Bool* save /* 0D_NOT_logical */);
void init_complex_taylor_series(
    ComplexTaylorProxy& complex_taylor,
    int n_term,
    std::optional<bool> save = std::nullopt);
extern "C" void fortran_init_coord1(
    void* orb /* 0D_NOT_type */,
    c_RealArr vec /* 1D_NOT_real */,
    void* ele /* 0D_NOT_type */,
    c_Int* element_end /* 0D_NOT_integer */,
    c_Int* particle /* 0D_NOT_integer */,
    c_Int* direction /* 0D_NOT_integer */,
    c_Real* E_photon /* 0D_NOT_real */,
    c_Real* t_offset /* 0D_NOT_real */,
    c_Bool* shift_vec6 /* 0D_NOT_logical */,
    c_RealArr spin /* 1D_NOT_real */,
    c_Real* s_pos /* 0D_NOT_real */,
    c_Bool* random_on /* 0D_NOT_logical */);
void init_coord1(
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
    void* orb_out /* 0D_NOT_type */,
    void* orb_in /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int* element_end /* 0D_NOT_integer */,
    c_Int* particle /* 0D_NOT_integer */,
    c_Int* direction /* 0D_NOT_integer */,
    c_Real* E_photon /* 0D_NOT_real */,
    c_Real* t_offset /* 0D_NOT_real */,
    c_Bool* shift_vec6 /* 0D_NOT_logical */,
    c_RealArr spin /* 1D_NOT_real */,
    c_Real* s_pos /* 0D_NOT_real */,
    c_Bool* random_on /* 0D_NOT_logical */);
CoordProxy init_coord2(
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
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int* element_end /* 0D_NOT_integer */,
    c_Int* particle /* 0D_NOT_integer */,
    c_Int* direction /* 0D_NOT_integer */,
    c_Real* E_photon /* 0D_NOT_real */,
    c_Real* t_offset /* 0D_NOT_real */,
    c_Bool* shift_vec6 /* 0D_NOT_logical */,
    c_RealArr spin /* 1D_NOT_real */);
void init_coord3(
    CoordProxy& orb,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<int> element_end = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> direction = std::nullopt,
    std::optional<double> E_photon = std::nullopt,
    std::optional<double> t_offset = std::nullopt,
    std::optional<bool> shift_vec6 = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> spin = std::nullopt);
extern "C" void fortran_init_custom(void* lat /* 0D_NOT_type */);
void init_custom(LatProxy& lat);

// Skipped unusable routine init_custom_def:
// Routine in configuration skip list
extern "C" void fortran_init_ele(
    void* ele /* 0D_NOT_type */,
    c_Int* key /* 0D_NOT_integer */,
    c_Int* sub_key /* 0D_NOT_integer */,
    c_Int* ix_ele /* 0D_NOT_integer */,
    void* branch /* 0D_NOT_type */);
EleProxy init_ele(
    std::optional<int> key = std::nullopt,
    std::optional<int> sub_key = std::nullopt,
    std::optional<int> ix_ele = std::nullopt,
    optional_ref<BranchProxy> branch = std::nullopt);
extern "C" void fortran_init_em_taylor_series(
    void* em_taylor /* 0D_NOT_type */,
    c_Int& n_term /* 0D_NOT_integer */,
    c_Bool* save_old /* 0D_NOT_logical */);
void init_em_taylor_series(
    EmTaylorProxy& em_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt);

// Skipped unusable routine init_fringe_info:
// Untranslated type: FringeFieldInfoProxy (0D_NOT_type)
extern "C" void fortran_init_lat(
    void* lat /* 0D_NOT_type */,
    c_Int* n /* 0D_NOT_integer */,
    c_Bool* init_beginning_ele /* 0D_NOT_logical */);
LatProxy init_lat(
    std::optional<int> n = std::nullopt,
    std::optional<bool> init_beginning_ele = std::nullopt);
extern "C" void fortran_init_multipole_cache(void* ele /* 0D_NOT_type */);
void init_multipole_cache(EleProxy& ele);
extern "C" void fortran_init_photon_from_a_photon_init_ele(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool* random_on /* 0D_NOT_logical */);
CoordProxy init_photon_from_a_photon_init_ele(
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> random_on = std::nullopt);
extern "C" bool fortran_init_photon_integ_prob(
    c_Real& gamma /* 0D_NOT_real */,
    c_Real& g /* 0D_NOT_real */,
    c_Real& E_min /* 0D_NOT_real */,
    c_Real& E_max /* 0D_NOT_real */,
    c_Real* vert_angle_min /* 0D_NOT_real */,
    c_Real* vert_angle_max /* 0D_NOT_real */,
    c_Bool* vert_angle_symmetric /* 0D_NOT_logical */,
    c_Real* energy_integ_prob /* 0D_NOT_real */,
    c_Real& E_photon /* 0D_NOT_real */,
    c_Real& integ_prob /* 0D_NOT_real */);
struct InitPhotonIntegProb {
  double E_photon;
  double integ_prob;
};
InitPhotonIntegProb init_photon_integ_prob(
    double gamma,
    double g,
    double E_min,
    double E_max,
    std::optional<double> vert_angle_min = std::nullopt,
    std::optional<double> vert_angle_max = std::nullopt,
    std::optional<bool> vert_angle_symmetric = std::nullopt,
    std::optional<double> energy_integ_prob = std::nullopt);
extern "C" void fortran_init_spin_distribution(
    void* beam_init /* 0D_NOT_type */,
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */);
BunchProxy init_spin_distribution(BeamInitProxy& beam_init, EleProxy& ele);
extern "C" void fortran_init_surface_segment(
    void* phot /* 0D_NOT_type */,
    c_Int& ix /* 0D_NOT_integer */,
    c_Int& iy /* 0D_NOT_integer */);
void init_surface_segment(PhotonElementProxy& phot, int ix, int iy);
extern "C" void fortran_init_taylor_series(
    void* bmad_taylor /* 0D_NOT_type */,
    c_Int& n_term /* 0D_NOT_integer */,
    c_Bool* save_old /* 0D_NOT_logical */);
void init_taylor_series(
    TaylorProxy& bmad_taylor,
    int n_term,
    std::optional<bool> save_old = std::nullopt);
extern "C" void fortran_init_wake(
    void* wake /* 0D_PTR_type */,
    c_Int& n_sr_long /* 0D_NOT_integer */,
    c_Int& n_sr_trans /* 0D_NOT_integer */,
    c_Int& n_sr_z /* 0D_NOT_integer */,
    c_Int& n_lr_mode /* 0D_NOT_integer */,
    c_Bool* always_allocate /* 0D_NOT_logical */);
WakeProxy init_wake(
    int n_sr_long,
    int n_sr_trans,
    int n_sr_z,
    int n_lr_mode,
    std::optional<bool> always_allocate = std::nullopt);
extern "C" void fortran_insert_element(
    void* lat /* 0D_NOT_type */,
    void* insert_ele /* 0D_NOT_type */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    void* orbit /* 1D_ALLOC_type */);
void insert_element(
    LatProxy& lat,
    EleProxy& insert_ele,
    int ix_ele,
    std::optional<int> ix_branch = std::nullopt,
    optional_ref<CoordProxyAlloc1D> orbit = std::nullopt);

// Skipped unusable routine integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)
extern "C" bool fortran_integrand_base(
    c_Real& t /* 0D_NOT_real */,
    void* args /* 1D_ALLOC_real */,
    c_Real& func_retval__ /* 0D_NOT_real */);
void integrand_base(double t, RealAlloc1D& args, double func_retval__);

// Skipped unusable routine integrand_base_cov:
// Untranslated type: CPtrProxy (0D_NOT_type)

// Skipped unusable routine integrand_zap:
// Untranslated type: CPtrProxy (0D_NOT_type)
extern "C" void fortran_integrate_psi(
    c_Real& bound /* 0D_NOT_real */,
    c_Real& p0 /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    c_Real& result /* 0D_NOT_real */);
double integrate_psi(double bound, double p0, FixedArray1D<Real, 8> args);
extern "C" void fortran_integrated_mats(
    void* eles /* 1D_ALLOC_type */,
    void* coos /* 1D_ALLOC_type */,
    c_ComplexArr Lambda /* 2D_NOT_complex */,
    c_ComplexArr Theta /* 2D_NOT_complex */,
    c_ComplexArr Iota /* 2D_NOT_complex */,
    void* mode /* 0D_NOT_type */);
void integrated_mats(
    EleProxyAlloc1D& eles,
    CoordProxyAlloc1D& coos,
    FixedArray2D<Complex, 6, 6> Lambda,
    FixedArray2D<Complex, 6, 6> Theta,
    FixedArray2D<Complex, 6, 6> Iota,
    NormalModesProxy& mode);
extern "C" void fortran_integration_timer_ele(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* start /* 0D_NOT_type */,
    void* orb_max /* 0D_NOT_type */,
    c_Real& tol /* 0D_NOT_real */);
void integration_timer_ele(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start,
    CoordProxy& orb_max,
    double tol);

// Skipped unusable routine integration_timer_fibre:
// Untranslated type: FibreProxy (0D_NOT_type)

// Skipped unusable routine interpolate_field:
// Untranslated type: Mesh3dProxy (0D_NOT_type)

// Skipped unusable routine interval1_coef_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ion_kick(
    void* orbit /* 0D_NOT_type */,
    c_RealArr r_beam /* 1D_NOT_real */,
    c_Real& n_beam_part /* 0D_NOT_real */,
    void* a_twiss /* 0D_NOT_type */,
    void* b_twiss /* 0D_NOT_type */,
    c_Real& sig_ee /* 0D_NOT_real */,
    c_RealArr kick /* 1D_NOT_real */);
FixedArray1D<Real, 3> ion_kick(
    CoordProxy& orbit,
    FixedArray1D<Real, 2> r_beam,
    double n_beam_part,
    TwissProxy& a_twiss,
    TwissProxy& b_twiss,
    double sig_ee);
extern "C" bool fortran_is_attribute(
    c_Int& ix_attrib /* 0D_NOT_integer */,
    c_Int& which /* 0D_NOT_integer */,
    c_Bool& is_attrib /* 0D_NOT_logical */);
bool is_attribute(int ix_attrib, int which);

// Skipped unusable routine jac:
// Untranslated type: CPtrProxy (0D_NOT_type)
// Untranslated type: CPtrProxy (0D_NOT_type)
// Untranslated type: CPtrProxy (0D_NOT_type)
// Untranslated type: CPtrProxy (0D_NOT_type)
extern "C" bool fortran_key_name_to_key_index(
    c_Char key_str /* 0D_NOT_character */,
    c_Bool* abbrev_allowed /* 0D_NOT_logical */,
    c_Int& key_index /* 0D_NOT_integer */);
void key_name_to_key_index(
    std::string key_str,
    std::optional<bool> abbrev_allowed,
    int key_index);
extern "C" void fortran_kick_vector_calc(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& s_body /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr dr_ds /* 1D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct KickVectorCalc {
  FixedArray1D<Real, 11> dr_ds;
  bool err;
};
KickVectorCalc kick_vector_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_body,
    CoordProxy& orbit,
    optional_ref<bool> print_err = std::nullopt);
extern "C" void fortran_kill_complex_taylor(
    void* complex_taylor /* 1D_ALLOC_type */);
void kill_complex_taylor(ComplexTaylorProxyAlloc1D& complex_taylor);
extern "C" void fortran_kill_ptc_layouts(void* lat /* 0D_NOT_type */);
void kill_ptc_layouts(LatProxy& lat);
extern "C" void fortran_kill_taylor(void* bmad_taylor /* 1D_ALLOC_type */);
void kill_taylor(TaylorProxyAlloc1D& bmad_taylor);
extern "C" bool fortran_kind_name(
    c_IntArr this_kind /* 0D_PTR_integer */,
    c_Char kind_str /* 0D_NOT_character */);
std::string kind_name(int this_kind);
extern "C" bool fortran_knot_interpolate(
    void* x_knot /* 1D_ALLOC_real */,
    void* y_knot /* 1D_ALLOC_real */,
    c_Real& x_pt /* 0D_NOT_real */,
    c_Int& interpolation /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& y_pt /* 0D_NOT_real */);
bool knot_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x_pt,
    int interpolation,
    double y_pt);
extern "C" bool fortran_knots_to_string(
    void* x_knot /* 1D_ALLOC_real */,
    void* y_knot /* 1D_ALLOC_real */,
    c_Char str /* 0D_ALLOC_character */);
void knots_to_string(RealAlloc1D& x_knot, RealAlloc1D& y_knot, std::string str);

// Skipped unusable routine kubo_integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)

// Skipped unusable routine kv_beam_init_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_lafun(
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void lafun(double x, double y, double z, double res);
extern "C" void fortran_lat_compute_ref_energy_and_time(
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool lat_compute_ref_energy_and_time(LatProxy& lat);

// Skipped unusable routine lat_ele_loc_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine lat_ele_locator:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine lat_ele_order1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine lat_ele_order_array_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine lat_ele_order_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_lat_equal_lat(
    void* lat_out /* 0D_NOT_type */,
    void* lat_in /* 0D_NOT_type */);
LatProxy lat_equal_lat(LatProxy& lat_in);
extern "C" void fortran_lat_geometry(void* lat /* 0D_NOT_type */);
void lat_geometry(LatProxy& lat);
extern "C" void fortran_lat_make_mat6(
    void* lat /* 0D_NOT_type */,
    c_Int* ix_ele /* 0D_NOT_integer */,
    void* ref_orb /* 1D_ALLOC_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool lat_make_mat6(
    LatProxy& lat,
    std::optional<int> ix_ele = std::nullopt,
    optional_ref<CoordProxyAlloc1D> ref_orb = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine lat_make_mat6_hook_def:
// Routine in configuration skip list

// Skipped unusable routine lat_param_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine lat_pointer_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_lat_sanity_check(
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool lat_sanity_check(LatProxy& lat);

// Skipped unusable routine lat_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_lat_to_ptc_layout(void* lat /* 0D_NOT_type */);
void lat_to_ptc_layout(LatProxy& lat);
extern "C" void fortran_lat_vec_equal_lat_vec(
    void* lat1 /* 1D_ALLOC_type */,
    void* lat2 /* 1D_ALLOC_type */);
LatProxyAlloc1D lat_vec_equal_lat_vec(LatProxyAlloc1D& lat2);
extern "C" void fortran_lattice_bookkeeper(
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool lattice_bookkeeper(LatProxy& lat);
extern "C" void fortran_lcavity_rf_step_setup(void* ele /* 0D_NOT_type */);
void lcavity_rf_step_setup(EleProxy& ele);

// Skipped unusable routine linac_normal_mode_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_linear_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void linear_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_linear_coef(
    void* stack /* 1D_ALLOC_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& coef /* 0D_NOT_real */);
struct LinearCoef {
  bool err_flag;
  double coef;
};
LinearCoef linear_coef(ExpressionAtomProxyAlloc1D& stack);

// Skipped unusable routine linear_ele_isf_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine linear_isf1_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_linear_to_spin_taylor(
    c_RealArr q_map /* 2D_NOT_real */,
    void* spin_taylor /* 1D_NOT_type */);
TaylorProxyArray1D linear_to_spin_taylor(FixedArray2D<Real, 4, 7> q_map);
extern "C" void fortran_load_parse_line(
    c_Char action /* 0D_NOT_character */,
    c_Int& ix_start /* 0D_NOT_integer */,
    c_Bool& end_of_file /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct LoadParseLine {
  bool end_of_file;
  bool err_flag;
};
LoadParseLine load_parse_line(std::string action, int ix_start);
extern "C" bool fortran_lord_edge_aligned(
    void* slave /* 0D_NOT_type */,
    c_Int& slave_edge /* 0D_NOT_integer */,
    void* lord /* 0D_NOT_type */,
    c_Bool& is_aligned /* 0D_NOT_logical */);
void lord_edge_aligned(
    EleProxy& slave,
    int slave_edge,
    EleProxy& lord,
    bool is_aligned);
extern "C" bool fortran_low_energy_z_correction(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& ds /* 0D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_Real& dz /* 0D_NOT_real */);
void low_energy_z_correction(
    CoordProxy& orbit,
    EleProxy& ele,
    double ds,
    std::optional<FixedArray2D<Real, 6, 6>> mat6,
    std::optional<bool> make_matrix,
    double dz);

// Skipped unusable routine lsc_kick_params_calc:
// Untranslated type: CsrProxy (0D_NOT_type)

// Skipped unusable routine mad_add_offsets_and_multipoles:
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_concat_map2:
// Untranslated type: MadMapProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_drift:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_elsep:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_energy_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine mad_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine mad_map_to_taylor:
// Untranslated type: MadMapProxy (0D_NOT_type)
// Untranslated type: MadEnergyProxy (0D_NOT_type)

// Skipped unusable routine mad_quadrupole:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_rfcavity:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_sbend:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_sbend_body:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_sbend_fringe:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_sextupole:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_solenoid:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)
extern "C" void fortran_mad_tmfoc(
    c_Real& el /* 0D_NOT_real */,
    c_Real& sk1 /* 0D_NOT_real */,
    c_Real& c /* 0D_NOT_real */,
    c_Real& s /* 0D_NOT_real */,
    c_Real& d /* 0D_NOT_real */,
    c_Real& f /* 0D_NOT_real */);
struct MadTmfoc {
  double c;
  double s;
  double d;
  double f;
};
MadTmfoc mad_tmfoc(double el, double sk1);
extern "C" void fortran_mad_tmsymm(c_RealArr te /* 3D_NOT_real */);
void mad_tmsymm(FixedArray3D<Real, 6, 6, 6> te);

// Skipped unusable routine mad_tmtilt:
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine mad_track1:
// Untranslated type: MadMapProxy (0D_NOT_type)
extern "C" void fortran_make_g2_mats(
    void* twiss /* 0D_NOT_type */,
    c_RealArr g2_mat /* 2D_NOT_real */,
    c_RealArr g2_inv_mat /* 2D_NOT_real */);
void make_g2_mats(
    TwissProxy& twiss,
    FixedArray2D<Real, 2, 2> g2_mat,
    FixedArray2D<Real, 2, 2> g2_inv_mat);
extern "C" void fortran_make_g_mats(
    void* ele /* 0D_NOT_type */,
    c_RealArr g_mat /* 2D_NOT_real */,
    c_RealArr g_inv_mat /* 2D_NOT_real */);
struct MakeGMats {
  FixedArray2D<Real, 4, 4> g_mat;
  FixedArray2D<Real, 4, 4> g_inv_mat;
};
MakeGMats make_g_mats(EleProxy& ele);
extern "C" void fortran_make_hvbp(
    c_RealArr N /* 2D_NOT_real */,
    c_RealArr B /* 2D_NOT_real */,
    c_RealArr V /* 2D_NOT_real */,
    c_RealArr H /* 2D_NOT_real */,
    c_RealArr Vbar /* 2D_NOT_real */,
    c_RealArr Hbar /* 2D_NOT_real */);
struct MakeHvbp {
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> V;
  FixedArray2D<Real, 6, 6> H;
  std::optional<FixedArray2D<Real, 6, 6>> Vbar;
  std::optional<FixedArray2D<Real, 6, 6>> Hbar;
};
MakeHvbp make_hvbp(FixedArray2D<Real, 6, 6> N);
extern "C" void fortran_make_hybrid_lat(
    void* lat_in /* 0D_NOT_type */,
    void* lat_out /* 0D_NOT_type */,
    c_Bool* use_taylor /* 0D_NOT_logical */,
    void* orb0_arr /* 1D_ALLOC_type */);
LatProxy make_hybrid_lat(
    LatProxy& lat_in,
    std::optional<bool> use_taylor = std::nullopt,
    optional_ref<CoordArrayProxyAlloc1D> orb0_arr = std::nullopt);

// Skipped unusable routine make_mad_map:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)
extern "C" void fortran_make_mat6(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct MakeMat6 {
  CoordProxy end_orb;
  bool err_flag;
};
MakeMat6 make_mat6(
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> start_orb = std::nullopt);
extern "C" void fortran_make_mat6_bmad(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
struct MakeMat6Bmad {
  CoordProxy end_orb;
  bool err;
};
MakeMat6Bmad make_mat6_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb);
extern "C" void fortran_make_mat6_bmad_photon(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
struct MakeMat6BmadPhoton {
  CoordProxy end_orb;
  bool err;
};
MakeMat6BmadPhoton make_mat6_bmad_photon(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb);

// Skipped unusable routine make_mat6_custom_def:
// Routine in configuration skip list
extern "C" void fortran_make_mat6_high_energy_space_charge(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */);
void make_mat6_high_energy_space_charge(EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_make_mat6_mad(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* c0 /* 0D_NOT_type */,
    void* c1 /* 0D_NOT_type */);
CoordProxy make_mat6_mad(EleProxy& ele, LatParamProxy& param, CoordProxy& c0);
extern "C" void fortran_make_mat6_symp_lie_ptc(
    void* ele /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */);
CoordProxy make_mat6_symp_lie_ptc(EleProxy& ele, CoordProxy& start_orb);
extern "C" void fortran_make_mat6_taylor(
    void* ele /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool* err_flag /* 0D_NOT_logical */);
CoordProxy make_mat6_taylor(
    EleProxy& ele,
    CoordProxy& start_orb,
    optional_ref<bool> err_flag = std::nullopt);
extern "C" void fortran_make_mat6_tracking(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* spin_only /* 0D_NOT_logical */);
struct MakeMat6Tracking {
  CoordProxy end_orb;
  bool err_flag;
};
MakeMat6Tracking make_mat6_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start_orb,
    std::optional<bool> spin_only = std::nullopt);
extern "C" void fortran_make_n(
    c_RealArr t6 /* 2D_NOT_real */,
    c_RealArr N /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_RealArr abz_tunes /* 1D_NOT_real */,
    c_RealArr tunes_out /* 1D_NOT_real */,
    c_RealArr U /* 2D_NOT_real */);
struct MakeN {
  FixedArray2D<Real, 6, 6> N;
  bool err_flag;
  FixedArray1D<Real, 3> tunes_out;
  std::optional<FixedArray2D<Real, 6, 6>> U;
};
MakeN make_n(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt);
extern "C" void fortran_make_pbrh(
    c_RealArr M /* 2D_NOT_real */,
    c_ComplexArr P /* 2D_NOT_complex */,
    c_ComplexArr Bp /* 2D_NOT_complex */,
    c_ComplexArr R /* 2D_NOT_complex */,
    c_ComplexArr H /* 2D_NOT_complex */,
    c_RealArr abz_tunes /* 1D_NOT_real */);
struct MakePbrh {
  FixedArray2D<Complex, 6, 6> P;
  FixedArray2D<Complex, 6, 6> Bp;
  FixedArray2D<Complex, 6, 6> R;
  FixedArray2D<Complex, 6, 6> H;
};
MakePbrh make_pbrh(FixedArray2D<Real, 6, 6> M, FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_make_smat_from_abc(
    c_RealArr t6 /* 2D_NOT_real */,
    void* mode /* 0D_NOT_type */,
    c_RealArr sigma_mat /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_RealArr Nout /* 2D_NOT_real */);
struct MakeSmatFromAbc {
  FixedArray2D<Real, 6, 6> sigma_mat;
  bool err_flag;
  std::optional<FixedArray2D<Real, 6, 6>> Nout;
};
MakeSmatFromAbc make_smat_from_abc(
    FixedArray2D<Real, 6, 6> t6,
    NormalModesProxy& mode);

// Skipped unusable routine make_sr_mats:
// Variable out sized array: M(:,:) 2D_NOT_real
// Variable out sized array: Bone(:,:) 2D_NOT_real
// Variable out sized array: Done(:,:) 2D_NOT_real

// Skipped unusable routine make_srdt_cache:
// Untranslated type: SlicedElesProxy (1D_ALLOC_type)
// Variable inout sized array: cache(:,:,:) 3D_ALLOC_complex

// Skipped unusable routine make_unit_mad_map:
// Untranslated type: MadMapProxy (0D_NOT_type)
extern "C" void fortran_make_v(
    c_RealArr M /* 2D_NOT_real */,
    c_ComplexArr V /* 2D_NOT_complex */,
    c_RealArr abz_tunes /* 1D_NOT_real */);
void make_v(
    FixedArray2D<Real, 6, 6> M,
    FixedArray2D<Complex, 6, 6> V,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_make_v_mats(
    void* ele /* 0D_NOT_type */,
    c_RealArr v_mat /* 2D_NOT_real */,
    c_RealArr v_inv_mat /* 2D_NOT_real */);
struct MakeVMats {
  std::optional<FixedArray2D<Real, 4, 4>> v_mat;
  std::optional<FixedArray2D<Real, 4, 4>> v_inv_mat;
};
MakeVMats make_v_mats(EleProxy& ele);

// Skipped unusable routine make_ykick_mat:
// Variable inout sized array: Yone(:,:) 2D_NOT_real
extern "C" void fortran_makeup_control_slave(
    void* lat /* 0D_NOT_type */,
    void* slave /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void makeup_control_slave(LatProxy& lat, EleProxy& slave, bool err_flag);
extern "C" void fortran_makeup_group_lord(
    void* lat /* 0D_NOT_type */,
    void* lord /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void makeup_group_lord(LatProxy& lat, EleProxy& lord, bool err_flag);
extern "C" void fortran_makeup_multipass_slave(
    void* lat /* 0D_NOT_type */,
    void* slave /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void makeup_multipass_slave(LatProxy& lat, EleProxy& slave, bool err_flag);
extern "C" void fortran_makeup_super_slave(
    void* lat /* 0D_NOT_type */,
    void* slave /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void makeup_super_slave(LatProxy& lat, EleProxy& slave, bool err_flag);
extern "C" void fortran_makeup_super_slave1(
    void* slave /* 0D_NOT_type */,
    void* lord /* 0D_NOT_type */,
    c_Real& offset /* 0D_NOT_real */,
    void* param /* 0D_NOT_type */,
    c_Bool& include_upstream_end /* 0D_NOT_logical */,
    c_Bool& include_downstream_end /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool makeup_super_slave1(
    EleProxy& slave,
    EleProxy& lord,
    double offset,
    LatParamProxy& param,
    bool include_upstream_end,
    bool include_downstream_end);
extern "C" bool fortran_map1_inverse(
    void* map1 /* 0D_NOT_type */,
    void* inv_map1 /* 0D_NOT_type */);
void map1_inverse(SpinOrbitMap1Proxy& map1, SpinOrbitMap1Proxy& inv_map1);
extern "C" void fortran_map1_make_unit(void* map1 /* 0D_NOT_type */);
SpinOrbitMap1Proxy map1_make_unit();
extern "C" bool fortran_map1_times_map1(
    void* map2 /* 0D_NOT_type */,
    void* map1 /* 0D_NOT_type */,
    void* map_out /* 0D_NOT_type */);
SpinOrbitMap1Proxy map1_times_map1(
    SpinOrbitMap1Proxy& map2,
    SpinOrbitMap1Proxy& map1);

// Skipped unusable routine map_coef:
// Untranslated type: Real8Proxy (1D_ALLOC_type)
extern "C" void fortran_map_to_angle_coords(
    void* t_canon /* 1D_NOT_type */,
    void* t_angle /* 1D_NOT_type */);
TaylorProxyArray1D map_to_angle_coords(FixedArray1D<TaylorProxy, 6> t_canon);
extern "C" void fortran_mark_patch_regions(void* branch /* 0D_NOT_type */);
void mark_patch_regions(BranchProxy& branch);
extern "C" bool fortran_master_parameter_value(
    c_Int& master_parameter /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Real& value /* 0D_NOT_real */);
void master_parameter_value(int master_parameter, EleProxy& ele, double value);
extern "C" void fortran_mat4_multipole(
    c_Real& knl /* 0D_NOT_real */,
    c_Real& tilt /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr kick_mat /* 2D_NOT_real */);
FixedArray2D<Real, 4, 4> mat4_multipole(
    double knl,
    double tilt,
    int n,
    CoordProxy& orbit);
extern "C" void fortran_mat6_add_offsets(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */);
void mat6_add_offsets(EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_mat6_add_pitch(
    c_Real& x_pitch_tot /* 0D_NOT_real */,
    c_Real& y_pitch_tot /* 0D_NOT_real */,
    c_Int& orientation /* 0D_NOT_integer */,
    c_RealArr mat6 /* 2D_NOT_real */);
void mat6_add_pitch(
    double x_pitch_tot,
    double y_pitch_tot,
    int orientation,
    FixedArray2D<Real, 6, 6> mat6);

// Skipped unusable routine mat6_from_s_to_s:
// Variable inout sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_mat6_to_complex_taylor(
    c_ComplexArr vec0 /* 1D_NOT_complex */,
    c_ComplexArr mat6 /* 2D_NOT_complex */,
    void* complex_taylor /* 1D_NOT_type */);
ComplexTaylorProxyArray1D mat6_to_complex_taylor(
    FixedArray1D<Complex, 6> vec0,
    FixedArray2D<Complex, 6, 6> mat6);
extern "C" void fortran_mat_symp_decouple(
    c_RealArr t0 /* 2D_NOT_real */,
    c_Int& stat /* 0D_NOT_integer */,
    c_RealArr U /* 2D_NOT_real */,
    c_RealArr V /* 2D_NOT_real */,
    c_RealArr Ubar /* 2D_NOT_real */,
    c_RealArr Vbar /* 2D_NOT_real */,
    c_RealArr G /* 2D_NOT_real */,
    void* twiss1 /* 0D_NOT_type */,
    void* twiss2 /* 0D_NOT_type */,
    c_Real& gamma /* 0D_NOT_real */,
    c_Bool& type_out /* 0D_NOT_logical */);
struct MatSympDecouple {
  int stat;
  TwissProxy twiss1;
  TwissProxy twiss2;
  double gamma;
};
MatSympDecouple mat_symp_decouple(
    FixedArray2D<Real, 4, 4> t0,
    FixedArray2D<Real, 4, 4> U,
    FixedArray2D<Real, 4, 4> V,
    FixedArray2D<Real, 4, 4> Ubar,
    FixedArray2D<Real, 4, 4> Vbar,
    FixedArray2D<Real, 4, 4> G,
    bool type_out);
extern "C" void fortran_match_ele_to_mat6(
    void* ele /* 0D_NOT_type */,
    void* start_orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_RealArr vec0 /* 1D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* include_delta_time /* 0D_NOT_logical */,
    c_Bool* set_trombone /* 0D_NOT_logical */);
struct MatchEleToMat6 {
  FixedArray2D<Real, 6, 6> mat6;
  FixedArray1D<Real, 6> vec0;
  bool err_flag;
};
MatchEleToMat6 match_ele_to_mat6(
    EleProxy& ele,
    CoordProxy& start_orb,
    std::optional<bool> include_delta_time = std::nullopt,
    std::optional<bool> set_trombone = std::nullopt);

// Skipped unusable routine material_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine mccfft1d:
// Routine module (fft_interface_mod) in configuration skip list

// Skipped unusable routine mesh3d_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_mexp(
    c_Real& x /* 0D_NOT_real */,
    c_Int& m /* 0D_NOT_integer */,
    c_Real& this_exp /* 0D_NOT_real */);
void mexp(double x, int m, double this_exp);
extern "C" void fortran_mfft1(
    void* a /* 1D_ALLOC_real */,
    void* b /* 1D_ALLOC_real */,
    void* n /* 1D_ALLOC_integer */,
    c_Int& ndim /* 0D_NOT_integer */,
    c_Int& isn /* 0D_NOT_integer */,
    c_Int& ierr /* 0D_NOT_integer */);
void mfft1(
    RealAlloc1D& a,
    RealAlloc1D& b,
    IntAlloc1D& n,
    int ndim,
    int isn,
    int ierr);

// Skipped unusable routine misalign_ptc_fibre:
// Untranslated type: FibreProxy (0D_PTR_type)

// Skipped unusable routine mode3_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine mode_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine momentum_aperture_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_momentum_compaction(
    void* branch /* 0D_NOT_type */,
    c_Real& mom_comp /* 0D_NOT_real */);
void momentum_compaction(BranchProxy& branch, double mom_comp);

// Skipped unusable routine mpxx1:
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine mpxx_integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)

// Skipped unusable routine mpzt1:
// Untranslated type: IbsProxy (0D_NOT_type)

// Skipped unusable routine multi_coulomb_log:
// Untranslated type: IbsSimParamProxy (0D_NOT_type)
extern "C" void fortran_multi_turn_tracking_analysis(
    void* track /* 1D_ALLOC_type */,
    c_Int& i_dim /* 0D_NOT_integer */,
    void* track0 /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& stable /* 0D_NOT_logical */,
    c_Real& growth_rate /* 0D_NOT_real */,
    c_Real& chi /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct MultiTurnTrackingAnalysis {
  CoordProxy track0;
  EleProxy ele;
  bool stable;
  double growth_rate;
  double chi;
  bool err_flag;
};
MultiTurnTrackingAnalysis multi_turn_tracking_analysis(
    CoordProxyAlloc1D& track,
    int i_dim);

// Skipped unusable routine multi_turn_tracking_to_mat:
// Variable out sized array: map1(:,:) 2D_NOT_real
extern "C" void fortran_multilayer_type_to_multilayer_params(
    void* ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool multilayer_type_to_multilayer_params(EleProxy& ele);

// Skipped unusable routine multipass_all_info:
// Untranslated type: MultipassAllInfoProxy (0D_NOT_type)

// Skipped unusable routine multipass_all_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_branch_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_chain:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine multipass_ele_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_lord_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_region_branch_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_region_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine multipass_region_info:
// Untranslated type: MultipassRegionLatProxy (0D_NOT_type)
// Untranslated type: MultipassAllInfoProxy (0D_NOT_type)

// Skipped unusable routine multipass_region_lat_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_multipole1_ab_to_kt(
    c_Real& an /* 0D_NOT_real */,
    c_Real& bn /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Real& knl /* 0D_NOT_real */,
    c_Real& tn /* 0D_NOT_real */);
struct Multipole1AbToKt {
  double knl;
  double tn;
};
Multipole1AbToKt multipole1_ab_to_kt(double an, double bn, int n);
extern "C" void fortran_multipole1_kt_to_ab(
    c_Real& knl /* 0D_NOT_real */,
    c_Real& knsl /* 0D_NOT_real */,
    c_Real& tn /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Real& an /* 0D_NOT_real */,
    c_Real& bn /* 0D_NOT_real */);
struct Multipole1KtToAb {
  double an;
  double bn;
};
Multipole1KtToAb multipole1_kt_to_ab(double knl, double knsl, double tn, int n);
extern "C" void fortran_multipole_ab_to_kt(
    void* an /* 1D_ALLOC_real */,
    void* bn /* 1D_ALLOC_real */,
    void* knl /* 1D_ALLOC_real */,
    void* tn /* 1D_ALLOC_real */);
struct MultipoleAbToKt {
  RealAlloc1D knl;
  RealAlloc1D tn;
};
MultipoleAbToKt multipole_ab_to_kt(RealAlloc1D& an, RealAlloc1D& bn);

// Skipped unusable routine multipole_cache_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_multipole_ele_to_ab(
    void* ele /* 0D_NOT_type */,
    c_Bool& use_ele_tilt /* 0D_NOT_logical */,
    c_Int& ix_pole_max /* 0D_NOT_integer */,
    c_RealArr a /* 1D_NOT_real */,
    c_RealArr b /* 1D_NOT_real */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Int* include_kicks /* 0D_NOT_integer */,
    c_Real& b1 /* 0D_NOT_real */,
    c_Bool* original /* 0D_NOT_logical */);
struct MultipoleEleToAb {
  int ix_pole_max;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> a;
  FixedArray1D<Real, Bmad::N_POLE_MAXX> b;
  double b1;
};
MultipoleEleToAb multipole_ele_to_ab(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt,
    std::optional<bool> original = std::nullopt);
extern "C" void fortran_multipole_ele_to_kt(
    void* ele /* 0D_NOT_type */,
    c_Bool& use_ele_tilt /* 0D_NOT_logical */,
    c_Int& ix_pole_max /* 0D_NOT_integer */,
    void* knl /* 1D_ALLOC_real */,
    void* tilt /* 1D_ALLOC_real */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Int* include_kicks /* 0D_NOT_integer */);
struct MultipoleEleToKt {
  int ix_pole_max;
  RealAlloc1D knl;
  RealAlloc1D tilt;
};
MultipoleEleToKt multipole_ele_to_kt(
    EleProxy& ele,
    bool use_ele_tilt,
    std::optional<int> pole_type = std::nullopt,
    std::optional<int> include_kicks = std::nullopt);
extern "C" void fortran_multipole_init(
    void* ele /* 0D_NOT_type */,
    c_Int& who /* 0D_NOT_integer */,
    c_Bool* zero /* 0D_NOT_logical */);
EleProxy multipole_init(int who, std::optional<bool> zero = std::nullopt);
extern "C" void fortran_multipole_kick(
    c_Real& knl /* 0D_NOT_real */,
    c_Real& tilt /* 0D_NOT_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Int& ref_species /* 0D_NOT_integer */,
    c_Int& ele_orientation /* 0D_NOT_integer */,
    void* coord /* 0D_NOT_type */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Bool* ref_orb_offset /* 0D_NOT_logical */);
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
    void* knl /* 1D_ALLOC_real */,
    void* tilt /* 1D_ALLOC_real */,
    c_Int& ref_species /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Real& factor /* 0D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */);
FixedArray2D<Real, 6, 6> multipole_kick_mat(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    int ref_species,
    EleProxy& ele,
    CoordProxy& orbit,
    double factor);
extern "C" void fortran_multipole_kicks(
    void* knl /* 1D_ALLOC_real */,
    void* tilt /* 1D_ALLOC_real */,
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int* pole_type /* 0D_NOT_integer */,
    c_Bool* ref_orb_offset /* 0D_NOT_logical */);
void multipole_kicks(
    RealAlloc1D& knl,
    RealAlloc1D& tilt,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<int> pole_type = std::nullopt,
    std::optional<bool> ref_orb_offset = std::nullopt);
extern "C" void fortran_multipole_kt_to_ab(
    void* knl /* 1D_ALLOC_real */,
    void* knsl /* 1D_ALLOC_real */,
    void* tn /* 1D_ALLOC_real */,
    void* an /* 1D_ALLOC_real */,
    void* bn /* 1D_ALLOC_real */);
struct MultipoleKtToAb {
  RealAlloc1D an;
  RealAlloc1D bn;
};
MultipoleKtToAb multipole_kt_to_ab(
    RealAlloc1D& knl,
    RealAlloc1D& knsl,
    RealAlloc1D& tn);
extern "C" void fortran_multipole_spin_tracking(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void multipole_spin_tracking(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" bool fortran_mytan(
    c_Real& y /* 0D_NOT_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& arg /* 0D_NOT_real */);
void mytan(double y, double x, double arg);
extern "C" bool fortran_n_attrib_string_max_len(
    c_Int& max_len /* 0D_NOT_integer */);
int n_attrib_string_max_len();
extern "C" void fortran_new_control(
    void* lat /* 0D_NOT_type */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Char ele_name /* 0D_NOT_character */);
void new_control(
    LatProxy& lat,
    int ix_ele,
    std::optional<std::string> ele_name = std::nullopt);
extern "C" bool fortran_nint_chk(
    c_Real& re_val /* 0D_NOT_real */,
    c_Int& int_val /* 0D_NOT_integer */);
int nint_chk(double re_val);
extern "C" void fortran_normal_form_complex_taylors(
    void* one_turn_taylor /* 1D_NOT_type */,
    c_Bool& rf_on /* 0D_NOT_logical */,
    void* F /* 1D_NOT_type */,
    void* L /* 1D_NOT_type */,
    void* A /* 1D_NOT_type */,
    void* A_inverse /* 1D_NOT_type */,
    c_Int* order /* 0D_NOT_integer */);
void normal_form_complex_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool rf_on,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> F = std::nullopt,
    std::optional<FixedArray1D<ComplexTaylorProxy, 6>> L = std::nullopt,
    std::optional<FixedArray1D<TaylorProxy, 6>> A = std::nullopt,
    std::optional<FixedArray1D<TaylorProxy, 6>> A_inverse = std::nullopt,
    optional_ref<int> order = std::nullopt);

// Skipped unusable routine normal_form_rd_terms:
// Untranslated type: Probe8Proxy (0D_NOT_type)
extern "C" void fortran_normal_form_taylors(
    void* one_turn_taylor /* 1D_NOT_type */,
    c_Bool& rf_on /* 0D_NOT_logical */,
    void* dhdj /* 1D_NOT_type */,
    void* A /* 1D_NOT_type */,
    void* A_inverse /* 1D_NOT_type */);
struct NormalFormTaylors {
  TaylorProxyArray1D dhdj;
  TaylorProxyArray1D A;
  TaylorProxyArray1D A_inverse;
};
NormalFormTaylors normal_form_taylors(
    FixedArray1D<TaylorProxy, 6> one_turn_taylor,
    bool rf_on);
extern "C" void fortran_normal_mode3_calc(
    c_RealArr t6 /* 2D_NOT_real */,
    c_RealArr tune /* 1D_NOT_real */,
    c_RealArr B /* 2D_NOT_real */,
    c_RealArr HV /* 2D_NOT_real */,
    c_Bool* above_transition /* 0D_NOT_logical */,
    c_RealArr abz_tunes /* 1D_NOT_real */);
struct NormalMode3Calc {
  FixedArray1D<Real, 3> tune;
  FixedArray2D<Real, 6, 6> B;
  FixedArray2D<Real, 6, 6> HV;
};
NormalMode3Calc normal_mode3_calc(
    FixedArray2D<Real, 6, 6> t6,
    std::optional<bool> above_transition = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> abz_tunes = std::nullopt);
extern "C" void fortran_normal_mode_dispersion(
    void* ele /* 0D_NOT_type */,
    c_Bool* reverse /* 0D_NOT_logical */);
void normal_mode_dispersion(
    EleProxy& ele,
    std::optional<bool> reverse = std::nullopt);

// Skipped unusable routine normal_modes_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_normalize_evecs(
    c_ComplexArr evec /* 2D_NOT_complex */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool normalize_evecs(FixedArray2D<Complex, 6, 6> evec);
extern "C" bool fortran_num_field_eles(
    void* ele /* 0D_NOT_type */,
    c_Int& n_field_ele /* 0D_NOT_integer */);
void num_field_eles(EleProxy& ele, int n_field_ele);
extern "C" bool fortran_num_lords(
    void* slave /* 0D_NOT_type */,
    c_Int& lord_type /* 0D_NOT_integer */,
    c_Int& num /* 0D_NOT_integer */);
void num_lords(EleProxy& slave, int lord_type, int num);
extern "C" void fortran_odeint_bmad(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& s1_body /* 0D_NOT_real */,
    c_Real& s2_body /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
struct OdeintBmad {
  bool err_flag;
  TrackProxy track;
};
OdeintBmad odeint_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    double s1_body,
    double s2_body,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_odeint_bmad_time(
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& t_dir /* 0D_NOT_integer */,
    c_Real& rf_time /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */,
    c_Real* t_end /* 0D_NOT_real */,
    c_Real& dt_step /* 0D_NOT_real */,
    void* extra_field /* 0D_NOT_type */);
struct OdeintBmadTime {
  bool err_flag;
  double dt_step;
};
OdeintBmadTime odeint_bmad_time(
    CoordProxy& orb,
    EleProxy& ele,
    LatParamProxy& param,
    int t_dir,
    double rf_time,
    optional_ref<TrackProxy> track = std::nullopt,
    std::optional<double> t_end = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" void fortran_offset_particle(
    void* ele /* 0D_NOT_type */,
    c_Bool& set /* 0D_NOT_logical */,
    void* orbit /* 0D_NOT_type */,
    c_Bool* set_tilt /* 0D_NOT_logical */,
    c_Bool* set_hvkicks /* 0D_NOT_logical */,
    c_Int* drift_to_edge /* 0D_NOT_integer */,
    c_Real* s_pos /* 0D_NOT_real */,
    c_Real& s_out /* 0D_NOT_real */,
    c_Bool* set_spin /* 0D_NOT_logical */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_RealArr spin_qrot /* 1D_NOT_real */,
    c_Real* time /* 0D_NOT_real */);
struct OffsetParticle {
  double s_out;
  FixedArray1D<Real, 4> spin_qrot;
};
OffsetParticle offset_particle(
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
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& set /* 0D_NOT_logical */,
    c_Bool* offset_position_only /* 0D_NOT_logical */,
    c_RealArr rot_mat /* 2D_NOT_real */);
void offset_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    std::optional<bool> offset_position_only = std::nullopt,
    std::optional<FixedArray2D<Real, 3, 3>> rot_mat = std::nullopt);
extern "C" void fortran_one_turn_mat_at_ele(
    void* ele /* 0D_NOT_type */,
    c_Real& phi_a /* 0D_NOT_real */,
    c_Real& phi_b /* 0D_NOT_real */,
    c_RealArr mat4 /* 2D_NOT_real */);
FixedArray2D<Real, 4, 4> one_turn_mat_at_ele(
    EleProxy& ele,
    double phi_a,
    double phi_b);
extern "C" bool fortran_open_binary_file(
    c_Char file_name /* 0D_NOT_character */,
    c_Char action /* 0D_NOT_character */,
    c_Int& iu /* 0D_NOT_integer */,
    c_Char r_name /* 0D_NOT_character */,
    c_Int& iver /* 0D_NOT_integer */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct OpenBinaryFile {
  int iu;
  int iver;
  bool is_ok;
};
OpenBinaryFile open_binary_file(
    std::string file_name,
    std::string action,
    std::string r_name);
extern "C" void fortran_orbit_amplitude_calc(
    void* ele /* 0D_NOT_type */,
    void* orb /* 0D_NOT_type */,
    c_Real& amp_a /* 0D_NOT_real */,
    c_Real& amp_b /* 0D_NOT_real */,
    c_Real& amp_na /* 0D_NOT_real */,
    c_Real& amp_nb /* 0D_NOT_real */);
struct OrbitAmplitudeCalc {
  double amp_a;
  double amp_b;
  double amp_na;
  double amp_nb;
};
OrbitAmplitudeCalc orbit_amplitude_calc(EleProxy& ele, CoordProxy& orb);
extern "C" void fortran_orbit_reference_energy_correction(
    void* orbit /* 0D_NOT_type */,
    c_Real& p0c_new /* 0D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void orbit_reference_energy_correction(
    CoordProxy& orbit,
    double p0c_new,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" bool fortran_orbit_to_floor_phase_space(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_RealArr floor_phase_space /* 1D_NOT_real */);
void orbit_to_floor_phase_space(
    CoordProxy& orbit,
    EleProxy& ele,
    FixedArray1D<Real, 6> floor_phase_space);
extern "C" bool fortran_orbit_to_local_curvilinear(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int* z_direction /* 0D_NOT_integer */,
    c_Int* relative_to /* 0D_NOT_integer */,
    void* local_position /* 0D_NOT_type */);
void orbit_to_local_curvilinear(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<int> z_direction,
    std::optional<int> relative_to,
    FloorPositionProxy& local_position);
extern "C" bool fortran_orbit_too_large(
    void* orbit /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool* check_momentum /* 0D_NOT_logical */,
    c_Bool& is_too_large /* 0D_NOT_logical */);
LatParamProxy orbit_too_large(
    CoordProxy& orbit,
    std::optional<bool> check_momentum,
    bool is_too_large);
extern "C" void fortran_order_evecs_by_n_similarity(
    c_ComplexArr evec /* 2D_NOT_complex */,
    c_ComplexArr eval /* 1D_NOT_complex */,
    c_RealArr mat_tunes /* 1D_NOT_real */,
    c_RealArr Nmat /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct OrderEvecsByNSimilarity {
  FixedArray2D<Complex, 6, 6> evec;
  bool err_flag;
};
OrderEvecsByNSimilarity order_evecs_by_n_similarity(
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray2D<Real, 6, 6> Nmat);
extern "C" void fortran_order_evecs_by_plane_dominance(
    c_ComplexArr evec /* 2D_NOT_complex */,
    c_ComplexArr eval /* 1D_NOT_complex */,
    c_RealArr mat_tunes /* 1D_NOT_real */);
void order_evecs_by_plane_dominance(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    std::optional<FixedArray1D<Real, 3>> mat_tunes = std::nullopt);
extern "C" void fortran_order_evecs_by_tune(
    c_ComplexArr evec /* 2D_NOT_complex */,
    c_ComplexArr eval /* 1D_NOT_complex */,
    c_RealArr mat_tunes /* 1D_NOT_real */,
    c_RealArr abz_tunes /* 1D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool order_evecs_by_tune(
    FixedArray2D<Complex, 6, 6> evec,
    FixedArray1D<Complex, 6> eval,
    FixedArray1D<Real, 3> mat_tunes,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_order_particles_in_z(void* bunch /* 0D_NOT_type */);
void order_particles_in_z(BunchProxy& bunch);
extern "C" void fortran_order_super_lord_slaves(
    void* lat /* 0D_NOT_type */,
    c_Int& ix_lord /* 0D_NOT_integer */);
void order_super_lord_slaves(LatProxy& lat, int ix_lord);
extern "C" void fortran_osc_alloc_freespace_array(
    c_IntArr nlo /* 1D_NOT_integer */,
    c_IntArr nhi /* 1D_NOT_integer */,
    c_IntArr npad /* 1D_NOT_integer */);
void osc_alloc_freespace_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_alloc_image_array(
    c_IntArr nlo /* 1D_NOT_integer */,
    c_IntArr nhi /* 1D_NOT_integer */,
    c_IntArr npad /* 1D_NOT_integer */);
void osc_alloc_image_array(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_alloc_rectpipe_arrays(
    c_IntArr nlo /* 1D_NOT_integer */,
    c_IntArr nhi /* 1D_NOT_integer */,
    c_IntArr npad /* 1D_NOT_integer */);
void osc_alloc_rectpipe_arrays(
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    FixedArray1D<Int, 3> npad);

// Skipped unusable routine osc_cathodeimages_solver:
// Variable inout sized array: rho(:,:,:) 3D_NOT_real
// Variable inout sized array: phi(:,:,:) 3D_NOT_real
// Variable inout sized array: efield(:,:,:,:) 4D_NOT_real
// Variable inout sized array: bfield(:,:,:,:) 4D_NOT_real

// Skipped unusable routine osc_freespace_solver:
// Variable inout sized array: rho(:,:,:) 3D_NOT_real
// Variable inout sized array: phi(:,:,:) 3D_NOT_real
// Variable inout sized array: efield(:,:,:,:) 4D_NOT_real
// Variable inout sized array: bfield(:,:,:,:) 4D_NOT_real

// Skipped unusable routine osc_freespace_solver2:
// Variable in sized array: rho(:,:,:) 3D_NOT_real
// Variable inout sized array: efield(:,:,:,:) 4D_NOT_real
// Variable inout sized array: phi(:,:,:) 3D_NOT_real

// Skipped unusable routine osc_get_cgrn_freespace:
// Variable inout sized array: cgrn(:,:,:) 3D_NOT_complex

// Skipped unusable routine osc_getgrnfree:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageconvcorr:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine osc_getgrnimageshift:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_osc_getgrnpipe(
    c_Real& gam /* 0D_NOT_real */,
    c_Real& a /* 0D_NOT_real */,
    c_Real& b /* 0D_NOT_real */,
    c_RealArr delta /* 1D_NOT_real */,
    c_RealArr umin /* 1D_NOT_real */,
    c_IntArr npad /* 1D_NOT_integer */);
void osc_getgrnpipe(
    double gam,
    double a,
    double b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad);
extern "C" void fortran_osc_read_rectpipe_grn();
void osc_read_rectpipe_grn();

// Skipped unusable routine osc_rectpipe_solver:
// Variable inout sized array: rho(:,:,:) 3D_NOT_real
// Variable inout sized array: phi(:,:,:) 3D_NOT_real
// Variable inout sized array: efield(:,:,:,:) 4D_NOT_real
// Variable inout sized array: bfield(:,:,:,:) 4D_NOT_real
extern "C" void fortran_osc_write_rectpipe_grn(
    c_Real& apipe /* 0D_NOT_real */,
    c_Real& bpipe /* 0D_NOT_real */,
    c_RealArr delta /* 1D_NOT_real */,
    c_RealArr umin /* 1D_NOT_real */,
    c_RealArr umax /* 1D_NOT_real */,
    c_IntArr nlo /* 1D_NOT_integer */,
    c_IntArr nhi /* 1D_NOT_integer */,
    c_Real& gamma /* 0D_NOT_real */);
void osc_write_rectpipe_grn(
    double apipe,
    double bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double gamma);
extern "C" void fortran_parse_cartesian_map(
    void* ct_map /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parse_cartesian_map(
    CartesianMapProxy& ct_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag);
extern "C" void fortran_parse_cylindrical_map(
    void* cl_map /* 0D_PTR_type */,
    void* ele /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parse_cylindrical_map(
    CylindricalMapProxy& cl_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag);
extern "C" void fortran_parse_gen_grad_map(
    void* gg_map /* 0D_PTR_type */,
    void* ele /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parse_gen_grad_map(
    GenGradMapProxy& gg_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag);
extern "C" void fortran_parse_grid_field(
    void* g_field /* 0D_PTR_type */,
    void* ele /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parse_grid_field(
    GridFieldProxy& g_field,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag);
extern "C" bool fortran_parse_integer_list(
    c_Char err_str /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* int_array /* 1D_ALLOC_integer */,
    c_Bool& exact_size /* 0D_NOT_logical */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Char open_delim /* 0D_NOT_character */,
    c_Char separator /* 0D_NOT_character */,
    c_Char close_delim /* 0D_NOT_character */,
    c_Int* default_value /* 0D_NOT_integer */,
    c_Bool& is_ok /* 0D_NOT_logical */);
void parse_integer_list(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    bool exact_size,
    std::string delim,
    bool delim_found,
    optional_ref<std::string> open_delim,
    optional_ref<std::string> separator,
    optional_ref<std::string> close_delim,
    optional_ref<int> default_value,
    bool is_ok);
extern "C" bool fortran_parse_integer_list2(
    c_Char err_str /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* int_array /* 1D_ALLOC_integer */,
    c_Int& num_found /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Int* num_expected /* 0D_NOT_integer */,
    c_Char open_delim /* 0D_NOT_character */,
    c_Char separator /* 0D_NOT_character */,
    c_Char close_delim /* 0D_NOT_character */,
    c_Int* default_value /* 0D_NOT_integer */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ParseIntegerList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
ParseIntegerList2 parse_integer_list2(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    optional_ref<int> num_expected = std::nullopt,
    optional_ref<std::string> open_delim = std::nullopt,
    optional_ref<std::string> separator = std::nullopt,
    optional_ref<std::string> close_delim = std::nullopt,
    optional_ref<int> default_value = std::nullopt);

// Skipped unusable routine parse_line_or_list:
// Untranslated type: SeqProxy (1D_ALLOC_type)
extern "C" bool fortran_parse_real_list(
    void* lat /* 0D_NOT_type */,
    c_Char err_str /* 0D_NOT_character */,
    void* real_array /* 1D_ALLOC_real */,
    c_Bool& exact_size /* 0D_NOT_logical */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Char open_delim /* 0D_NOT_character */,
    c_Char separator /* 0D_NOT_character */,
    c_Char close_delim /* 0D_NOT_character */,
    c_Real* default_value /* 0D_NOT_real */,
    c_Int& num_found /* 0D_NOT_integer */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ParseRealList {
  RealAlloc1D real_array;
  std::string delim;
  bool delim_found;
  int num_found;
  bool is_ok;
};
ParseRealList parse_real_list(
    LatProxy& lat,
    std::string err_str,
    bool exact_size,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<double> default_value = std::nullopt);
extern "C" bool fortran_parse_real_list2(
    void* lat /* 0D_NOT_type */,
    c_Char err_str /* 0D_NOT_character */,
    void* real_array /* 1D_ALLOC_real */,
    c_Int& num_found /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Int* num_expected /* 0D_NOT_integer */,
    c_Char open_brace /* 0D_NOT_character */,
    c_Char separator /* 0D_NOT_character */,
    c_Char close_brace /* 0D_NOT_character */,
    c_Real* default_value /* 0D_NOT_real */,
    c_Bool* single_value /* 0D_NOT_logical */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ParseRealList2 {
  int num_found;
  std::string delim;
  bool delim_found;
  bool is_ok;
};
ParseRealList2 parse_real_list2(
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
// Variable in sized array: table(:,:) 2D_ALLOC_real

// Skipped unusable routine parse_superimpose_command:
// Untranslated type: ParserEleProxy (0D_NOT_type)

// Skipped unusable routine parser2_add_superimpose:
// Untranslated type: ParserEleProxy (0D_NOT_type)

// Skipped unusable routine parser_add_branch:
// Untranslated type: SeqProxy (1D_ALLOC_type)
// Variable-sized inout character array: seq_name(:) 1D_ALLOC_character
// Untranslated type: ParserLatProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_parser_add_constant(
    c_Char word /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool& redef_is_error /* 0D_NOT_logical */);
void parser_add_constant(std::string word, LatProxy& lat, bool redef_is_error);

// Skipped unusable routine parser_add_lords:
// Untranslated type: ParserLatProxy (0D_NOT_type)

// Skipped unusable routine parser_add_superimpose:
// Untranslated type: ParserEleProxy (0D_NOT_type)
// Untranslated type: ParserLatProxy (0D_NOT_type)
extern "C" void fortran_parser_call_check(
    c_Char word /* 0D_NOT_character */,
    c_Int& ix_word /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& call_found /* 0D_NOT_logical */,
    c_Bool* err_flag /* 0D_NOT_logical */);
void parser_call_check(
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool call_found,
    optional_ref<bool> err_flag = std::nullopt);

// Skipped unusable routine parser_controller_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine parser_debug_print_info:
// Untranslated type: SeqProxy (1D_ALLOC_type)

// Skipped unusable routine parser_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine parser_error:
// Untranslated type: SeqProxy (0D_NOT_type)
// Untranslated type: ParserEleProxy (0D_NOT_type)

// Skipped unusable routine parser_expand_line:
// Untranslated type: SeqProxy (1D_ALLOC_type)
// Variable-sized in character array: seq_name(:) 1D_ALLOC_character
// Untranslated type: BaseLineEleProxy (1D_ALLOC_type)
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_parser_fast_complex_read(
    void* cmplx_vec /* 1D_ALLOC_complex */,
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Char err_str /* 0D_NOT_character */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ParserFastComplexRead {
  ComplexAlloc1D cmplx_vec;
  std::string delim;
  bool is_ok;
};
ParserFastComplexRead parser_fast_complex_read(
    EleProxy& ele,
    std::string err_str);
extern "C" bool fortran_parser_fast_integer_read(
    void* int_vec /* 1D_ALLOC_integer */,
    void* ele /* 0D_NOT_type */,
    c_Char delim_wanted /* 0D_NOT_character */,
    c_Char err_str /* 0D_NOT_character */,
    c_Bool& is_ok /* 0D_NOT_logical */);
void parser_fast_integer_read(
    IntAlloc1D& int_vec,
    EleProxy& ele,
    std::string delim_wanted,
    std::string err_str,
    bool is_ok);
extern "C" bool fortran_parser_fast_real_read(
    void* real_vec /* 1D_ALLOC_real */,
    void* ele /* 0D_NOT_type */,
    c_Char end_delims /* 0D_NOT_character */,
    c_Char delim /* 0D_NOT_character */,
    c_Char err_str /* 0D_NOT_character */,
    c_Bool* exact_size /* 0D_NOT_logical */,
    c_Int& n_real /* 0D_NOT_integer */,
    c_Bool& is_ok /* 0D_NOT_logical */);
struct ParserFastRealRead {
  RealAlloc1D real_vec;
  std::string delim;
  int n_real;
  bool is_ok;
};
ParserFastRealRead parser_fast_real_read(
    EleProxy& ele,
    std::string end_delims,
    std::string err_str,
    std::optional<bool> exact_size = std::nullopt);
extern "C" void fortran_parser_file_stack(
    c_Char how /* 0D_NOT_character */,
    c_Char file_name_in /* 0D_NOT_character */,
    c_Bool* finished /* 0D_NOT_logical */,
    c_Bool* err /* 0D_NOT_logical */,
    c_Bool* open_file /* 0D_NOT_logical */,
    c_Bool* abort_on_open_error /* 0D_NOT_logical */);
void parser_file_stack(
    std::string how,
    optional_ref<std::string> file_name_in = std::nullopt,
    optional_ref<bool> finished = std::nullopt,
    optional_ref<bool> err = std::nullopt,
    optional_ref<bool> open_file = std::nullopt,
    optional_ref<bool> abort_on_open_error = std::nullopt);
extern "C" void fortran_parser_get_integer(
    c_Int& int_val /* 0D_NOT_integer */,
    c_Char word /* 0D_NOT_character */,
    c_Int& ix_word /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Char str1 /* 0D_NOT_character */,
    c_Char str2 /* 0D_NOT_character */);
void parser_get_integer(
    int int_val,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool err,
    optional_ref<std::string> str1 = std::nullopt,
    optional_ref<std::string> str2 = std::nullopt);
extern "C" void fortran_parser_get_logical(
    c_Char attrib_name /* 0D_NOT_character */,
    c_Bool& this_logic /* 0D_NOT_logical */,
    c_Char ele_name /* 0D_NOT_character */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err /* 0D_NOT_logical */);
void parser_get_logical(
    std::string attrib_name,
    bool this_logic,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool err);
extern "C" void fortran_parser_identify_fork_to_element(
    void* lat /* 0D_NOT_type */);
void parser_identify_fork_to_element(LatProxy& lat);
extern "C" void fortran_parser_init_custom_elements(
    void* lat /* 0D_NOT_type */);
void parser_init_custom_elements(LatProxy& lat);

// Skipped unusable routine parser_lat_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_parser_print_line(
    void* lat /* 0D_NOT_type */,
    c_Bool& end_of_file /* 0D_NOT_logical */);
void parser_print_line(LatProxy& lat, bool end_of_file);
extern "C" void fortran_parser_read_lr_wake(
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parser_read_lr_wake(
    EleProxy& ele,
    std::string delim,
    bool delim_found,
    bool err_flag);
extern "C" void fortran_parser_read_old_format_lr_wake(
    void* ele /* 0D_NOT_type */,
    c_Char lr_file_name /* 0D_NOT_character */);
void parser_read_old_format_lr_wake(EleProxy& ele, std::string lr_file_name);
extern "C" void fortran_parser_read_old_format_sr_wake(
    void* ele /* 0D_NOT_type */,
    c_Char sr_file_name /* 0D_NOT_character */);
void parser_read_old_format_sr_wake(EleProxy& ele, std::string sr_file_name);
extern "C" void fortran_parser_read_sr_wake(
    void* ele /* 0D_NOT_type */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void parser_read_sr_wake(
    EleProxy& ele,
    std::string delim,
    bool delim_found,
    bool err_flag);

// Skipped unusable routine parser_set_attribute:
// Untranslated type: ParserEleProxy (0D_NOT_type)
extern "C" void fortran_parser_transfer_control_struct(
    void* con_in /* 0D_NOT_type */,
    void* con_out /* 0D_NOT_type */,
    void* lord /* 0D_NOT_type */,
    c_Int& ix_var /* 0D_NOT_integer */);
ControlProxy parser_transfer_control_struct(
    ControlProxy& con_in,
    EleProxy& lord,
    int ix_var);
extern "C" bool fortran_particle_in_global_frame(
    void* orb /* 0D_NOT_type */,
    void* branch /* 0D_NOT_type */,
    c_Bool* in_time_coordinates /* 0D_NOT_logical */,
    c_Bool* in_body_frame /* 0D_NOT_logical */,
    c_RealArr w_mat_out /* 2D_NOT_real */,
    void* particle /* 0D_NOT_type */);
void particle_in_global_frame(
    CoordProxy& orb,
    BranchProxy& branch,
    std::optional<bool> in_time_coordinates,
    std::optional<bool> in_body_frame,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat_out,
    CoordProxy& particle);
extern "C" bool fortran_particle_is_moving_backwards(
    void* orbit /* 0D_NOT_type */,
    c_Bool& is_moving_backwards /* 0D_NOT_logical */);
void particle_is_moving_backwards(CoordProxy& orbit, bool is_moving_backwards);
extern "C" bool fortran_particle_is_moving_forward(
    void* orbit /* 0D_NOT_type */,
    c_Int* dir /* 0D_NOT_integer */,
    c_Bool& is_moving_forward /* 0D_NOT_logical */);
void particle_is_moving_forward(
    CoordProxy& orbit,
    std::optional<int> dir,
    bool is_moving_forward);
extern "C" bool fortran_particle_rf_time(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool* reference_active_edge /* 0D_NOT_logical */,
    c_Real* s_rel /* 0D_NOT_real */,
    c_Bool* time_coords /* 0D_NOT_logical */,
    c_Real* rf_freq /* 0D_NOT_real */,
    c_Bool* abs_time /* 0D_NOT_logical */,
    const long double& time /* 0D_NOT_real16 */);
void particle_rf_time(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<bool> reference_active_edge,
    std::optional<double> s_rel,
    std::optional<bool> time_coords,
    std::optional<double> rf_freq,
    std::optional<bool> abs_time,
    long double time);
extern "C" bool fortran_patch_flips_propagation_direction(
    c_Real& x_pitch /* 0D_NOT_real */,
    c_Real& y_pitch /* 0D_NOT_real */,
    c_Bool& is_flip /* 0D_NOT_logical */);
void patch_flips_propagation_direction(
    double x_pitch,
    double y_pitch,
    bool is_flip);
extern "C" bool fortran_patch_length(
    void* patch /* 0D_NOT_type */,
    c_Int* ref_coords /* 0D_NOT_integer */,
    c_Real& length /* 0D_NOT_real */);
void patch_length(
    EleProxy& patch,
    std::optional<int> ref_coords,
    double length);

// Skipped unusable routine pauli_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_photon_absorption_and_phase_shift(
    c_Char material /* 0D_NOT_character */,
    c_Real& Energy /* 0D_NOT_real */,
    c_Real& absorption /* 0D_NOT_real */,
    c_Real& phase_shift /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct PhotonAbsorptionAndPhaseShift {
  double absorption;
  double phase_shift;
  bool err_flag;
};
PhotonAbsorptionAndPhaseShift photon_absorption_and_phase_shift(
    std::string material,
    double Energy);
extern "C" void fortran_photon_add_to_detector_statistics(
    void* orbit0 /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int* ix_pt /* 0D_NOT_integer */,
    c_Int* iy_pt /* 0D_NOT_integer */,
    void* pixel_pt /* 0D_NOT_type */);
void photon_add_to_detector_statistics(
    CoordProxy& orbit0,
    CoordProxy& orbit,
    EleProxy& ele,
    optional_ref<int> ix_pt = std::nullopt,
    optional_ref<int> iy_pt = std::nullopt,
    optional_ref<PixelPtProxy> pixel_pt = std::nullopt);

// Skipped unusable routine photon_coord_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_diffuse_scattering:
// Untranslated type: DiffuseParamProxy (0D_NOT_type)

// Skipped unusable routine photon_element_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_init_spline_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_init_spline_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_init_splines_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_init_x_angle_spline_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_init_y_angle_spline_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_material_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_read_spline:
// Untranslated type: PhotonInitSplinesProxy (0D_NOT_type)

// Skipped unusable routine photon_reflect_surface_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_reflect_table_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_photon_reflection(
    c_Real& graze_angle_in /* 0D_NOT_real */,
    c_Real& energy /* 0D_NOT_real */,
    void* surface /* 0D_NOT_type */,
    c_Real& graze_angle_out /* 0D_NOT_real */,
    c_Real& phi_out /* 0D_NOT_real */);
struct PhotonReflection {
  double graze_angle_out;
  double phi_out;
};
PhotonReflection photon_reflection(
    double graze_angle_in,
    double energy,
    PhotonReflectSurfaceProxy& surface);
extern "C" void fortran_photon_reflection_std_surface_init(
    void* surface /* 0D_NOT_type */);
PhotonReflectSurfaceProxy photon_reflection_std_surface_init();
extern "C" void fortran_photon_reflectivity(
    c_Real& angle /* 0D_NOT_real */,
    c_Real& energy /* 0D_NOT_real */,
    void* surface /* 0D_NOT_type */,
    c_Real& p_reflect /* 0D_NOT_real */,
    c_Real& rel_p_specular /* 0D_NOT_real */);
struct PhotonReflectivity {
  double p_reflect;
  double rel_p_specular;
};
PhotonReflectivity photon_reflectivity(
    double angle,
    double energy,
    PhotonReflectSurfaceProxy& surface);
extern "C" void fortran_photon_target_corner_calc(
    void* aperture_ele /* 0D_NOT_type */,
    c_Real& x_lim /* 0D_NOT_real */,
    c_Real& y_lim /* 0D_NOT_real */,
    c_Real& z_lim /* 0D_NOT_real */,
    void* source_ele /* 0D_NOT_type */,
    void* corner /* 0D_NOT_type */);
TargetPointProxy photon_target_corner_calc(
    EleProxy& aperture_ele,
    double x_lim,
    double y_lim,
    double z_lim,
    EleProxy& source_ele);
extern "C" void fortran_photon_target_setup(void* ele /* 0D_NOT_type */);
void photon_target_setup(EleProxy& ele);

// Skipped unusable routine photon_target_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine photon_track_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_photon_type(
    void* ele /* 0D_NOT_type */,
    c_Int& e_type /* 0D_NOT_integer */);
int photon_type(EleProxy& ele);
extern "C" bool fortran_physical_ele_end(
    c_Int& track_end /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_Int& ele_orientation /* 0D_NOT_integer */,
    c_Bool* return_stream_end /* 0D_NOT_logical */,
    c_Int& physical_end /* 0D_NOT_integer */);
void physical_ele_end(
    int track_end,
    CoordProxy& orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end,
    int physical_end);

// Skipped unusable routine pixel_detec_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine pixel_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine pmd_header_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine pmd_unit_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_point_photon_emission(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int& direction /* 0D_NOT_integer */,
    c_Real& max_target_area /* 0D_NOT_real */,
    c_RealArr w_to_surface /* 2D_NOT_real */);
void point_photon_emission(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    int direction,
    double max_target_area,
    std::optional<FixedArray2D<Real, 3, 3>> w_to_surface = std::nullopt);

// Skipped unusable routine pointer_to_attribute:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" bool fortran_pointer_to_branch_given_ele(
    void* ele /* 0D_NOT_type */,
    void* branch_ptr /* 0D_PTR_type */);
BranchProxy pointer_to_branch_given_ele(EleProxy& ele);
extern "C" bool fortran_pointer_to_branch_given_name(
    c_Char branch_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool* parameter_is_branch0 /* 0D_NOT_logical */,
    c_Int* blank_branch /* 0D_NOT_integer */,
    void* branch_ptr /* 0D_PTR_type */);
BranchProxy pointer_to_branch_given_name(
    std::string branch_name,
    LatProxy& lat,
    std::optional<bool> parameter_is_branch0 = std::nullopt,
    std::optional<int> blank_branch = std::nullopt);
extern "C" bool fortran_pointer_to_ele1(
    void* lat /* 0D_NOT_type */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    void* ele_ptr /* 0D_PTR_type */);
void pointer_to_ele1(
    LatProxy& lat,
    int ix_ele,
    optional_ref<int> ix_branch,
    EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele2(
    void* lat /* 0D_NOT_type */,
    void* ele_loc /* 0D_NOT_type */,
    void* ele_ptr /* 0D_PTR_type */);
void pointer_to_ele2(LatProxy& lat, LatEleLocProxy& ele_loc, EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele3(
    void* lat /* 0D_NOT_type */,
    c_Char ele_name /* 0D_NOT_character */,
    void* ele_ptr /* 0D_PTR_type */);
void pointer_to_ele3(LatProxy& lat, std::string ele_name, EleProxy& ele_ptr);
extern "C" bool fortran_pointer_to_ele4(
    void* lat /* 0D_NOT_type */,
    void* foreign_ele /* 0D_NOT_type */,
    void* ele_ptr /* 0D_PTR_type */);
void pointer_to_ele4(LatProxy& lat, EleProxy& foreign_ele, EleProxy& ele_ptr);

// Skipped unusable routine pointer_to_ele_multipole:
// Routine in configuration skip list
extern "C" bool fortran_pointer_to_element_at_s(
    void* branch /* 0D_NOT_type */,
    c_Real& s /* 0D_NOT_real */,
    c_Bool& choose_max /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& s_eff /* 0D_NOT_real */,
    void* position /* 0D_NOT_type */,
    c_Bool* print_err /* 0D_NOT_logical */,
    void* ele /* 0D_PTR_type */);
struct PointerToElementAtS {
  bool err_flag;
  double s_eff;
  CoordProxy position;
  EleProxy ele;
};
PointerToElementAtS pointer_to_element_at_s(
    BranchProxy& branch,
    double s,
    bool choose_max,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine pointer_to_fibre:
// Untranslated type: FibreProxy (0D_PTR_type)
extern "C" bool fortran_pointer_to_field_ele(
    void* ele /* 0D_NOT_type */,
    c_Int& ix_field_ele /* 0D_NOT_integer */,
    c_Real& dz_offset /* 0D_NOT_real */,
    void* field_ele /* 0D_PTR_type */);
double pointer_to_field_ele(
    EleProxy& ele,
    int ix_field_ele,
    EleProxy& field_ele);
extern "C" bool fortran_pointer_to_girder(
    void* ele /* 0D_NOT_type */,
    c_Int& ix_slave_back /* 0D_NOT_integer */,
    void* girder /* 0D_PTR_type */);
int pointer_to_girder(EleProxy& ele, EleProxy& girder);

// Skipped unusable routine pointer_to_indexed_attribute:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" bool fortran_pointer_to_lord(
    void* slave /* 0D_NOT_type */,
    c_Int& ix_lord /* 0D_NOT_integer */,
    void* control /* 0D_PTR_type */,
    c_Int& ix_slave_back /* 0D_NOT_integer */,
    c_Int* lord_type /* 0D_NOT_integer */,
    c_Int& ix_control /* 0D_NOT_integer */,
    c_Int& ix_ic /* 0D_NOT_integer */,
    void* lord_ptr /* 0D_PTR_type */);
struct PointerToLord {
  ControlProxy control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
};
PointerToLord pointer_to_lord(
    EleProxy& slave,
    int ix_lord,
    std::optional<int> lord_type,
    EleProxy& lord_ptr);
extern "C" bool fortran_pointer_to_multipass_lord(
    void* ele /* 0D_NOT_type */,
    c_Int& ix_pass /* 0D_NOT_integer */,
    void* super_lord /* 0D_PTR_type */,
    void* multi_lord /* 0D_PTR_type */);
struct PointerToMultipassLord {
  int ix_pass;
  EleProxy super_lord;
};
PointerToMultipassLord pointer_to_multipass_lord(
    EleProxy& ele,
    EleProxy& multi_lord);
extern "C" bool fortran_pointer_to_next_ele(
    void* this_ele /* 0D_NOT_type */,
    c_Int* offset /* 0D_NOT_integer */,
    c_Bool* skip_beginning /* 0D_NOT_logical */,
    c_Bool* follow_fork /* 0D_NOT_logical */,
    void* next_ele /* 0D_PTR_type */);
void pointer_to_next_ele(
    EleProxy& this_ele,
    std::optional<int> offset,
    std::optional<bool> skip_beginning,
    std::optional<bool> follow_fork,
    EleProxy& next_ele);
extern "C" bool fortran_pointer_to_slave(
    void* lord /* 0D_NOT_type */,
    c_Int& ix_slave /* 0D_NOT_integer */,
    void* control /* 0D_PTR_type */,
    c_Int* lord_type /* 0D_NOT_integer */,
    c_Int& ix_lord_back /* 0D_NOT_integer */,
    c_Int& ix_control /* 0D_NOT_integer */,
    c_Int& ix_ic /* 0D_NOT_integer */,
    void* slave_ptr /* 0D_PTR_type */);
struct PointerToSlave {
  ControlProxy control;
  int ix_lord_back;
  int ix_control;
  int ix_ic;
  EleProxy slave_ptr;
};
PointerToSlave pointer_to_slave(
    EleProxy& lord,
    int ix_slave,
    std::optional<int> lord_type = std::nullopt);
extern "C" bool fortran_pointer_to_super_lord(
    void* slave /* 0D_NOT_type */,
    void* control /* 0D_PTR_type */,
    c_Int& ix_slave_back /* 0D_NOT_integer */,
    c_Int& ix_control /* 0D_NOT_integer */,
    c_Int& ix_ic /* 0D_NOT_integer */,
    c_Int* lord_type /* 0D_NOT_integer */,
    void* lord_ptr /* 0D_PTR_type */);
struct PointerToSuperLord {
  ControlProxy control;
  int ix_slave_back;
  int ix_control;
  int ix_ic;
};
PointerToSuperLord pointer_to_super_lord(
    EleProxy& slave,
    std::optional<int> lord_type,
    EleProxy& lord_ptr);
extern "C" bool fortran_pointer_to_surface_displacement_pt(
    void* ele /* 0D_NOT_type */,
    c_Bool& nearest /* 0D_NOT_logical */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Int* ix /* 0D_NOT_integer */,
    c_Int* iy /* 0D_NOT_integer */,
    c_Bool* extend_grid /* 0D_NOT_logical */,
    c_Real* xx /* 0D_NOT_real */,
    c_Real* yy /* 0D_NOT_real */,
    void* pt /* 0D_PTR_type */);
SurfaceDisplacementPtProxy pointer_to_surface_displacement_pt(
    EleProxy& ele,
    bool nearest,
    double x,
    double y,
    optional_ref<int> ix = std::nullopt,
    optional_ref<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    optional_ref<double> xx = std::nullopt,
    optional_ref<double> yy = std::nullopt);
extern "C" bool fortran_pointer_to_surface_segmented_pt(
    void* ele /* 0D_NOT_type */,
    c_Bool& nearest /* 0D_NOT_logical */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Int* ix /* 0D_NOT_integer */,
    c_Int* iy /* 0D_NOT_integer */,
    c_Bool* extend_grid /* 0D_NOT_logical */,
    c_Real* xx /* 0D_NOT_real */,
    c_Real* yy /* 0D_NOT_real */,
    void* pt /* 0D_PTR_type */);
SurfaceSegmentedPtProxy pointer_to_surface_segmented_pt(
    EleProxy& ele,
    bool nearest,
    double x,
    double y,
    optional_ref<int> ix = std::nullopt,
    optional_ref<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    optional_ref<double> xx = std::nullopt,
    optional_ref<double> yy = std::nullopt);
extern "C" bool fortran_pointer_to_wake_ele(
    void* ele /* 0D_NOT_type */,
    c_Real& delta_s /* 0D_NOT_real */,
    void* wake_ele /* 0D_PTR_type */);
double pointer_to_wake_ele(EleProxy& ele, EleProxy& wake_ele);
extern "C" bool fortran_pointer_to_wall3d(
    void* ele /* 0D_NOT_type */,
    c_Int* ix_wall /* 0D_NOT_integer */,
    c_Real& ds_offset /* 0D_NOT_real */,
    c_Bool& is_branch_wall /* 0D_NOT_logical */,
    void* wall3d /* 0D_PTR_type */);
struct PointerToWall3d {
  double ds_offset;
  bool is_branch_wall;
  Wall3dProxy wall3d;
};
PointerToWall3d pointer_to_wall3d(
    EleProxy& ele,
    std::optional<int> ix_wall = std::nullopt);

// Skipped unusable routine pointers_to_attribute:
// Untranslated type: AllPointerProxy (1D_ALLOC_type)
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" bool fortran_polar_to_spinor(
    void* polar /* 0D_NOT_type */,
    c_ComplexArr spinor /* 1D_NOT_complex */);
void polar_to_spinor(SpinPolarProxy& polar, FixedArray1D<Complex, 2> spinor);
extern "C" bool fortran_polar_to_vec(
    void* polar /* 0D_NOT_type */,
    c_RealArr vec /* 1D_NOT_real */);
void polar_to_vec(SpinPolarProxy& polar, FixedArray1D<Real, 3> vec);

// Skipped unusable routine pre_tracker_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine print_mesh3d:
// Untranslated type: Mesh3dProxy (0D_NOT_type)

// Skipped unusable routine prob_x_diffuse:
// Untranslated type: DiffuseParamProxy (0D_NOT_type)
extern "C" void fortran_project_emit_to_xyz(
    void* ring /* 0D_NOT_type */,
    c_Int& ix /* 0D_NOT_integer */,
    void* mode /* 0D_NOT_type */,
    c_Real& sigma_x /* 0D_NOT_real */,
    c_Real& sigma_y /* 0D_NOT_real */,
    c_Real& sigma_z /* 0D_NOT_real */);
struct ProjectEmitToXyz {
  double sigma_x;
  double sigma_y;
  double sigma_z;
};
ProjectEmitToXyz project_emit_to_xyz(
    LatProxy& ring,
    int ix,
    NormalModesProxy& mode);

// Skipped unusable routine propagate_part_way:
// Untranslated type: RadIntTrackPointProxy (0D_NOT_type)
// Untranslated type: RadIntInfoProxy (0D_NOT_type)

// Skipped unusable routine psi_prime:
// Untranslated type: CPtrProxy (0D_NOT_type)
// Untranslated type: CPtrProxy (0D_NOT_type)
// Untranslated type: CPtrProxy (0D_NOT_type)
extern "C" void fortran_psi_prime_sca(
    c_Real& t /* 0D_NOT_real */,
    c_Real& p /* 0D_NOT_real */,
    c_Real& dpdt /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */);
double psi_prime_sca(double t, double p, FixedArray1D<Real, 8> args);
extern "C" void fortran_ptc_bookkeeper(void* lat /* 0D_NOT_type */);
void ptc_bookkeeper(LatProxy& lat);

// Skipped unusable routine ptc_branch1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ptc_calculate_tracking_step_size:
// Untranslated type: LayoutProxy (0D_NOT_type)

// Skipped unusable routine ptc_check_for_lost_particle:
// Untranslated type: FibreProxy (0D_PTR_type)
extern "C" void fortran_ptc_closed_orbit_calc(
    void* branch /* 0D_NOT_type */,
    void* closed_orbit /* 1D_ALLOC_type */,
    c_Bool* radiation_damping_on /* 0D_NOT_logical */);
CoordProxyAlloc1D ptc_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<bool> radiation_damping_on = std::nullopt);

// Skipped unusable routine ptc_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ptc_emit_calc(
    void* ele /* 0D_NOT_type */,
    void* norm_mode /* 0D_NOT_type */,
    c_RealArr sigma_mat /* 2D_NOT_real */,
    void* closed_orb /* 0D_NOT_type */);
struct PtcEmitCalc {
  NormalModesProxy norm_mode;
  CoordProxy closed_orb;
};
PtcEmitCalc ptc_emit_calc(EleProxy& ele, FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_kill_map_with_radiation:
// Untranslated type: PtcRadMapProxy (0D_NOT_type)

// Skipped unusable routine ptc_layout_pointer_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ptc_layouts_resplit(
    c_Real& dKL_max /* 0D_NOT_real */,
    c_Real& l_max /* 0D_NOT_real */,
    c_Bool& l_max_drift_only /* 0D_NOT_logical */,
    c_Real& bend_dorb /* 0D_NOT_real */,
    c_Real& sex_dx /* 0D_NOT_real */,
    c_Bool* even /* 0D_NOT_logical */,
    c_IntArr crossover /* 1D_NOT_integer */,
    c_IntArr crossover_wiggler /* 1D_NOT_integer */);
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
// Untranslated type: LinearEleIsfProxy (1D_ALLOC_type)

// Skipped unusable routine ptc_map_to_normal_form:
// Untranslated type: Probe8Proxy (0D_NOT_type)
// Untranslated type: CNormalFormProxy (0D_NOT_type)
// Untranslated type: CTaylorProxy (1D_NOT_type)
// Untranslated type: CTaylorProxy (0D_NOT_type)

// Skipped unusable routine ptc_normal_form_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ptc_one_turn_map_at_ele:
// Untranslated type: Probe8Proxy (0D_NOT_type)
// Untranslated type: InternalStateProxy (0D_NOT_type)
extern "C" void fortran_ptc_one_turn_mat_and_closed_orbit_calc(
    void* branch /* 0D_NOT_type */,
    c_Real* pz /* 0D_NOT_real */);
void ptc_one_turn_mat_and_closed_orbit_calc(
    BranchProxy& branch,
    std::optional<double> pz = std::nullopt);

// Skipped unusable routine ptc_private_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine ptc_rad_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ptc_ran_seed_put(c_Int& iseed /* 0D_NOT_integer */);
void ptc_ran_seed_put(int iseed);

// Skipped unusable routine ptc_read_flat_file:
// Variable-sized in character array: flat_file(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine ptc_read_map_with_radiation:
// Untranslated type: PtcRadMapProxy (0D_NOT_type)
extern "C" void fortran_ptc_set_rf_state_for_c_normal(
    c_Bool& nocavity /* 0D_NOT_logical */);
void ptc_set_rf_state_for_c_normal(bool nocavity);
extern "C" void fortran_ptc_set_taylor_order_if_needed();
void ptc_set_taylor_order_if_needed();

// Skipped unusable routine ptc_setup_map_with_radiation:
// Untranslated type: PtcRadMapProxy (0D_NOT_type)

// Skipped unusable routine ptc_setup_tracking_with_damping_and_excitation:
// Untranslated type: InternalStateProxy (0D_NOT_type)
extern "C" void fortran_ptc_spin_calc(
    void* ele /* 0D_NOT_type */,
    void* norm_mode /* 0D_NOT_type */,
    c_RealArr sigma_mat /* 2D_NOT_real */,
    void* closed_orb /* 0D_NOT_type */);
struct PtcSpinCalc {
  NormalModesProxy norm_mode;
  CoordProxy closed_orb;
};
PtcSpinCalc ptc_spin_calc(EleProxy& ele, FixedArray2D<Real, 6, 6> sigma_mat);

// Skipped unusable routine ptc_spin_matching_calc:
// Untranslated type: SpinMatchingProxy (1D_ALLOC_type)

// Skipped unusable routine ptc_taylors_equal_bmad_taylors:
// Untranslated type: TaylorProxy (1D_ALLOC_type)
extern "C" void fortran_ptc_track_all(
    void* branch /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */,
    c_Int& track_state /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct PtcTrackAll {
  int track_state;
  bool err_flag;
};
PtcTrackAll ptc_track_all(BranchProxy& branch, CoordProxyAlloc1D& orbit);

// Skipped unusable routine ptc_track_map_with_radiation:
// Untranslated type: PtcRadMapProxy (0D_NOT_type)
extern "C" void fortran_ptc_transfer_map_with_spin(
    void* branch /* 0D_NOT_type */,
    void* t_map /* 1D_NOT_type */,
    void* s_map /* 1D_NOT_type */,
    void* orb0 /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Int* ix1 /* 0D_NOT_integer */,
    c_Int* ix2 /* 0D_NOT_integer */,
    c_Bool* one_turn /* 0D_NOT_logical */,
    c_Bool* unit_start /* 0D_NOT_logical */);
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
// Untranslated type: PtcRadMapProxy (0D_NOT_type)

// Skipped unusable routine ptwo:
// Untranslated type: DiffuseParamProxy (0D_NOT_type)
extern "C" bool fortran_pwd_mat(
    void* lat /* 0D_NOT_type */,
    c_RealArr t6 /* 2D_NOT_real */,
    c_Real& inductance /* 0D_NOT_real */,
    c_Real& sig_z /* 0D_NOT_real */,
    c_RealArr t6_pwd /* 2D_NOT_real */);
FixedArray2D<Real, 6, 6> pwd_mat(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> t6,
    double inductance,
    double sig_z);

// Skipped unusable routine qromb_rad_int:
// Untranslated type: RadIntTrackPointProxy (0D_NOT_type)
// Untranslated type: RadIntInfoProxy (0D_NOT_type)

// Skipped unusable routine quad_mat2_calc:
// Variable out sized array: mat2(:,:) 2D_NOT_real
extern "C" void fortran_rad1_damp_and_stoc_mats(
    void* ele /* 0D_NOT_type */,
    c_Bool& include_opening_angle /* 0D_NOT_logical */,
    void* orb_in /* 0D_NOT_type */,
    void* orb_out /* 0D_NOT_type */,
    void* rad_map /* 0D_NOT_type */,
    c_Real& g2_tol /* 0D_NOT_real */,
    c_Real& g3_tol /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ele0 /* 0D_NOT_type */,
    void* rad_int1 /* 0D_NOT_type */);
struct Rad1DampAndStocMats {
  RadMapProxy rad_map;
  bool err_flag;
  RadInt1Proxy rad_int1;
};
Rad1DampAndStocMats rad1_damp_and_stoc_mats(
    EleProxy& ele,
    bool include_opening_angle,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double g2_tol,
    double g3_tol,
    optional_ref<EleProxy> ele0 = std::nullopt);
extern "C" void fortran_rad_damp_and_stoc_mats(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& include_opening_angle /* 0D_NOT_logical */,
    void* rmap /* 0D_NOT_type */,
    void* mode /* 0D_NOT_type */,
    c_RealArr xfer_nodamp_mat /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* closed_orbit /* 1D_ALLOC_type */,
    void* rad_int_branch /* 0D_NOT_type */);
struct RadDampAndStocMats {
  RadMapProxy rmap;
  NormalModesProxy mode;
  FixedArray2D<Real, 6, 6> xfer_nodamp_mat;
  bool err_flag;
  RadIntBranchProxy rad_int_branch;
};
RadDampAndStocMats rad_damp_and_stoc_mats(
    EleProxy& ele1,
    EleProxy& ele2,
    bool include_opening_angle,
    optional_ref<CoordProxyAlloc1D> closed_orbit = std::nullopt);
extern "C" void fortran_rad_g_integrals(
    void* ele /* 0D_NOT_type */,
    c_Int& where /* 0D_NOT_integer */,
    void* orb_in /* 0D_NOT_type */,
    void* orb_out /* 0D_NOT_type */,
    c_RealArr int_g /* 1D_NOT_real */,
    c_Real& int_g2 /* 0D_NOT_real */,
    c_Real& int_g3 /* 0D_NOT_real */,
    c_Real& g_tol /* 0D_NOT_real */,
    c_Real& g2_tol /* 0D_NOT_real */,
    c_Real& g3_tol /* 0D_NOT_real */);
FixedArray1D<Real, 2> rad_g_integrals(
    EleProxy& ele,
    int where,
    CoordProxy& orb_in,
    CoordProxy& orb_out,
    double int_g2,
    double int_g3,
    double g_tol,
    double g2_tol,
    double g3_tol);

// Skipped unusable routine rad_int1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_all_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_branch_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_cache1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_cache_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_info_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_int_track_point_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_map_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine rad_map_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_radiation_integrals(
    void* lat /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */,
    void* mode /* 0D_NOT_type */,
    c_Int* ix_cache /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    void* rad_int_by_ele /* 0D_NOT_type */);
struct RadiationIntegrals {
  NormalModesProxy mode;
  RadIntAllEleProxy rad_int_by_ele;
};
RadiationIntegrals radiation_integrals(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    optional_ref<int> ix_cache = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine radiation_integrals_custom_def:
// Routine in configuration skip list
extern "C" void fortran_radiation_map_setup(
    void* ele /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ref_orbit_in /* 0D_NOT_type */);
bool radiation_map_setup(
    EleProxy& ele,
    optional_ref<CoordProxy> ref_orbit_in = std::nullopt);

// Skipped unusable routine ramper_lord_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_ramper_slave_setup(
    void* lat /* 0D_NOT_type */,
    c_Bool* force_setup /* 0D_NOT_logical */);
void ramper_slave_setup(
    LatProxy& lat,
    std::optional<bool> force_setup = std::nullopt);
extern "C" bool fortran_ramper_value(
    void* ramper /* 0D_NOT_type */,
    void* r1 /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
bool ramper_value(EleProxy& ramper, ControlRamp1Proxy& r1, double value);
extern "C" void fortran_randomize_lr_wake_frequencies(
    void* ele /* 0D_NOT_type */,
    c_Bool& set_done /* 0D_NOT_logical */);
bool randomize_lr_wake_frequencies(EleProxy& ele);
extern "C" bool fortran_rchomp(
    c_Real& rel /* 0D_NOT_real */,
    c_Int& plc /* 0D_NOT_integer */,
    c_Char out /* 0D_NOT_character */);
void rchomp(double rel, int plc, std::string out);

// Skipped unusable routine rclog_integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)

// Skipped unusable routine re_allocate_eles:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" void fortran_re_allocate_wall3d_section_array(
    void* section /* 1D_ALLOC_type */,
    c_Int& n /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */);
void re_allocate_wall3d_section_array(
    Wall3dSectionProxyAlloc1D& section,
    int n,
    std::optional<bool> exact = std::nullopt);
extern "C" void fortran_re_allocate_wall3d_vertex_array(
    void* v /* 1D_ALLOC_type */,
    c_Int& n /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */);
void re_allocate_wall3d_vertex_array(
    Wall3dVertexProxyAlloc1D& v,
    int n,
    std::optional<bool> exact = std::nullopt);

// Skipped unusable routine re_associate_node_array:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)
extern "C" bool fortran_re_str_qp(
    const long double& rel /* 0D_NOT_real16 */,
    c_Char str_out /* 0D_NOT_character */);
void re_str_qp(long double rel, std::string str_out);
extern "C" bool fortran_re_str_rp(
    c_Real& rel /* 0D_NOT_real */,
    c_Char str_out /* 0D_NOT_character */);
void re_str_rp(double rel, std::string str_out);
extern "C" void fortran_read_beam_ascii(
    c_Char file_name /* 0D_NOT_character */,
    void* beam /* 0D_NOT_type */,
    void* beam_init /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct ReadBeamAscii {
  BeamProxy beam;
  bool err_flag;
};
ReadBeamAscii read_beam_ascii(std::string file_name, BeamInitProxy& beam_init);
extern "C" void fortran_read_beam_file(
    c_Char file_name /* 0D_NOT_character */,
    void* beam /* 0D_NOT_type */,
    void* beam_init /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */,
    c_Bool* print_mom_shift_warning /* 0D_NOT_logical */,
    c_Bool* conserve_momentum /* 0D_NOT_logical */);
struct ReadBeamFile {
  BeamProxy beam;
  bool err_flag;
};
ReadBeamFile read_beam_file(
    std::string file_name,
    BeamInitProxy& beam_init,
    optional_ref<EleProxy> ele = std::nullopt,
    std::optional<bool> print_mom_shift_warning = std::nullopt,
    optional_ref<bool> conserve_momentum = std::nullopt);
extern "C" void fortran_read_binary_cartesian_map(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* cart_map /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void read_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag);
extern "C" void fortran_read_binary_cylindrical_map(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* cl_map /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void read_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag);
extern "C" void fortran_read_binary_grid_field(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* g_field /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void read_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag);

// Skipped unusable routine read_digested_bmad_file:
// Variable-sized out character array: lat_files(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_read_surface_reflection_file(
    c_Char file_name /* 0D_NOT_character */,
    void* surface /* 0D_NOT_type */);
PhotonReflectSurfaceProxy read_surface_reflection_file(std::string file_name);

// Skipped unusable routine real_8_to_taylor:
// Untranslated type: Real8Proxy (1D_ALLOC_type)
extern "C" void fortran_reallocate_beam(
    void* beam /* 0D_NOT_type */,
    c_Int& n_bunch /* 0D_NOT_integer */,
    c_Int* n_particle /* 0D_NOT_integer */,
    c_Bool* extend /* 0D_NOT_logical */);
void reallocate_beam(
    BeamProxy& beam,
    int n_bunch,
    std::optional<int> n_particle = std::nullopt,
    optional_ref<bool> extend = std::nullopt);
extern "C" void fortran_reallocate_bp_com_const();
void reallocate_bp_com_const();
extern "C" void fortran_reallocate_bunch(
    void* bunch /* 0D_NOT_type */,
    c_Int& n_particle /* 0D_NOT_integer */,
    c_Bool* save /* 0D_NOT_logical */);
BunchProxy reallocate_bunch(
    int n_particle,
    std::optional<bool> save = std::nullopt);
extern "C" void fortran_reallocate_control(
    void* lat /* 0D_NOT_type */,
    c_Int& n /* 0D_NOT_integer */);
void reallocate_control(LatProxy& lat, int n);
extern "C" void fortran_reallocate_coord_array(
    void* coord_array /* 1D_ALLOC_type */,
    void* lat /* 0D_NOT_type */);
void reallocate_coord_array(CoordArrayProxyAlloc1D& coord_array, LatProxy& lat);
extern "C" void fortran_reallocate_coord_lat(
    void* coord /* 1D_ALLOC_type */,
    void* lat /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */);
void reallocate_coord_lat(
    CoordProxyAlloc1D& coord,
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_reallocate_coord_n(
    void* coord /* 1D_ALLOC_type */,
    c_Int& n_coord /* 0D_NOT_integer */);
void reallocate_coord_n(CoordProxyAlloc1D& coord, int n_coord);
extern "C" void fortran_reallocate_expression_stack(
    void* stack /* 1D_ALLOC_type */,
    c_Int& n /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */);
void reallocate_expression_stack(
    ExpressionAtomProxyAlloc1D& stack,
    int n,
    std::optional<bool> exact = std::nullopt);

// Skipped unusable routine reallocate_sequence:
// Untranslated type: SeqProxy (1D_ALLOC_type)

// Skipped unusable routine reals_8_equal_bmad_taylors:
// Untranslated type: Real8Proxy (1D_ALLOC_type)
extern "C" bool fortran_rel_tracking_charge_to_mass(
    void* orbit /* 0D_NOT_type */,
    c_Int& ref_species /* 0D_NOT_integer */,
    c_Real& rel_charge /* 0D_NOT_real */);
void rel_tracking_charge_to_mass(
    CoordProxy& orbit,
    int ref_species,
    double rel_charge);
extern "C" bool fortran_relative_mode_flip(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& func_retval__ /* 0D_NOT_logical */);
void relative_mode_flip(EleProxy& ele1, EleProxy& ele2, bool func_retval__);
extern "C" void fortran_release_rad_int_cache(
    c_Int& ix_cache /* 0D_NOT_integer */);
void release_rad_int_cache(int ix_cache);
extern "C" void fortran_remove_constant_taylor(
    void* taylor_in /* 1D_ALLOC_type */,
    void* taylor_out /* 1D_ALLOC_type */,
    void* c0 /* 1D_ALLOC_real */,
    c_Bool& remove_higher_order_terms /* 0D_NOT_logical */);
struct RemoveConstantTaylor {
  TaylorProxyAlloc1D taylor_out;
  RealAlloc1D c0;
};
RemoveConstantTaylor remove_constant_taylor(
    TaylorProxyAlloc1D& taylor_in,
    bool remove_higher_order_terms);
extern "C" void fortran_remove_dead_from_bunch(
    void* bunch_in /* 0D_NOT_type */,
    void* bunch_out /* 0D_NOT_type */);
BunchProxy remove_dead_from_bunch(BunchProxy& bunch_in);
extern "C" void fortran_remove_eles_from_lat(
    void* lat /* 0D_NOT_type */,
    c_Bool* check_sanity /* 0D_NOT_logical */);
void remove_eles_from_lat(
    LatProxy& lat,
    std::optional<bool> check_sanity = std::nullopt);
extern "C" void fortran_remove_lord_slave_link(
    void* lord /* 0D_NOT_type */,
    void* slave /* 0D_NOT_type */);
void remove_lord_slave_link(EleProxy& lord, EleProxy& slave);

// Skipped unusable routine resonance_h_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_reverse_lat(
    void* lat_in /* 0D_NOT_type */,
    void* lat_rev /* 0D_NOT_type */,
    c_Bool* track_antiparticle /* 0D_NOT_logical */);
LatProxy reverse_lat(
    LatProxy& lat_in,
    std::optional<bool> track_antiparticle = std::nullopt);
extern "C" void fortran_rf_coupler_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    c_Real& phase /* 0D_NOT_real */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void rf_coupler_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    double phase,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);

// Skipped unusable routine rf_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_rf_is_on(
    void* branch /* 0D_NOT_type */,
    c_Int* ix_ele1 /* 0D_NOT_integer */,
    c_Int* ix_ele2 /* 0D_NOT_integer */,
    c_Bool& is_on /* 0D_NOT_logical */);
void rf_is_on(
    BranchProxy& branch,
    std::optional<int> ix_ele1,
    std::optional<int> ix_ele2,
    bool is_on);
extern "C" bool fortran_rf_ref_time_offset(
    void* ele /* 0D_NOT_type */,
    c_Real* ds /* 0D_NOT_real */,
    c_Real& time /* 0D_NOT_real */);
void rf_ref_time_offset(EleProxy& ele, std::optional<double> ds, double time);

// Skipped unusable routine rf_stair_step_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_rfun(
    c_Real& u /* 0D_NOT_real */,
    c_Real& v /* 0D_NOT_real */,
    c_Real& w /* 0D_NOT_real */,
    c_Real& gam /* 0D_NOT_real */,
    c_Real& a /* 0D_NOT_real */,
    c_Real& b /* 0D_NOT_real */,
    c_Real& hz /* 0D_NOT_real */,
    c_Int& i /* 0D_NOT_integer */,
    c_Int& j /* 0D_NOT_integer */,
    c_Real& res /* 0D_NOT_real */);
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
    double res);
extern "C" void fortran_rk_adaptive_time_step(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orb /* 0D_NOT_type */,
    c_Int& t_dir /* 0D_NOT_integer */,
    c_Real& rf_time /* 0D_NOT_real */,
    c_Real& dt_try /* 0D_NOT_real */,
    c_Real& dt_did /* 0D_NOT_real */,
    c_Real& dt_next /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* extra_field /* 0D_NOT_type */);
void rk_adaptive_time_step(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orb,
    int t_dir,
    double rf_time,
    double dt_try,
    double dt_did,
    double dt_next,
    bool err_flag,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" void fortran_rk_time_step1(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& rf_time /* 0D_NOT_real */,
    void* orb /* 0D_NOT_type */,
    c_Real& dt /* 0D_NOT_real */,
    void* new_orb /* 0D_NOT_type */,
    c_RealArr r_err /* 1D_NOT_real */,
    c_RealArr dr_dt /* 1D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    void* extra_field /* 0D_NOT_type */);
void rk_time_step1(
    EleProxy& ele,
    LatParamProxy& param,
    double rf_time,
    CoordProxy& orb,
    double dt,
    CoordProxy& new_orb,
    FixedArray1D<Real, 10> r_err,
    std::optional<FixedArray1D<Real, 10>> dr_dt,
    bool err_flag,
    optional_ref<bool> print_err = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt);
extern "C" bool fortran_rotate3(
    c_RealArr vec /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */,
    c_RealArr rvec /* 1D_NOT_real */);
void rotate3(
    FixedArray1D<Real, 3> vec,
    double angle,
    FixedArray1D<Real, 3> rvec);
extern "C" void fortran_rotate_em_field(
    void* field /* 0D_NOT_type */,
    c_RealArr w_mat /* 2D_NOT_real */,
    c_RealArr w_inv /* 2D_NOT_real */,
    c_Bool* calc_dfield /* 0D_NOT_logical */,
    c_Bool* calc_potential /* 0D_NOT_logical */);
void rotate_em_field(
    EmFieldProxy& field,
    FixedArray2D<Real, 3, 3> w_mat,
    FixedArray2D<Real, 3, 3> w_inv,
    std::optional<bool> calc_dfield = std::nullopt,
    std::optional<bool> calc_potential = std::nullopt);
extern "C" void fortran_rotate_field_zx(
    void* field /* 0D_NOT_type */,
    c_Real& theta /* 0D_NOT_real */);
void rotate_field_zx(EmFieldProxy& field, double theta);
extern "C" void fortran_rotate_for_curved_surface(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& set /* 0D_NOT_logical */,
    c_RealArr rot_mat /* 2D_NOT_real */);
void rotate_for_curved_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    bool set,
    FixedArray2D<Real, 3, 3> rot_mat);
extern "C" void fortran_rotate_spin(
    c_RealArr rot_vec /* 1D_NOT_real */,
    c_RealArr spin /* 1D_NOT_real */,
    c_RealArr qrot /* 1D_NOT_real */);
FixedArray1D<Real, 4> rotate_spin(
    FixedArray1D<Real, 3> rot_vec,
    FixedArray1D<Real, 3> spin);
extern "C" void fortran_rotate_spin_a_step(
    void* orbit /* 0D_NOT_type */,
    void* field /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& ds /* 0D_NOT_real */);
void rotate_spin_a_step(
    CoordProxy& orbit,
    EmFieldProxy& field,
    EleProxy& ele,
    double ds);
extern "C" void fortran_rotate_spin_given_field(
    void* orbit /* 0D_NOT_type */,
    c_Int& sign_z_vel /* 0D_NOT_integer */,
    c_RealArr BL /* 1D_NOT_real */,
    c_RealArr EL /* 1D_NOT_real */,
    c_RealArr qrot /* 1D_NOT_real */);
void rotate_spin_given_field(
    CoordProxy& orbit,
    int sign_z_vel,
    std::optional<FixedArray1D<Real, 3>> BL = std::nullopt,
    std::optional<FixedArray1D<Real, 3>> EL = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> qrot = std::nullopt);

// Skipped unusable routine runge_kutta_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_s_body_calc(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& s_body /* 0D_NOT_real */);
void s_body_calc(CoordProxy& orbit, EleProxy& ele, double s_body);
extern "C" void fortran_s_calc(void* lat /* 0D_NOT_type */);
void s_calc(LatProxy& lat);

// Skipped unusable routine s_ref_to_s_chord:
// Untranslated type: CsrEleInfoProxy (0D_NOT_type)

// Skipped unusable routine s_source_calc:
// Untranslated type: CsrKick1Proxy (0D_NOT_type)
// Untranslated type: CsrProxy (0D_NOT_type)
extern "C" void fortran_sad_mult_hard_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void sad_mult_hard_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_sad_soft_bend_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orb /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void sad_soft_bend_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orb,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_save_a_beam_step(
    void* ele /* 0D_NOT_type */,
    void* beam /* 0D_NOT_type */,
    void* bunch_tracks /* 1D_ALLOC_type */,
    c_Real* s_body /* 0D_NOT_real */,
    c_Bool* is_time_coords /* 0D_NOT_logical */);
void save_a_beam_step(
    EleProxy& ele,
    BeamProxy& beam,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt);
extern "C" void fortran_save_a_bunch_step(
    void* ele /* 0D_NOT_type */,
    void* bunch /* 0D_NOT_type */,
    void* bunch_track /* 0D_NOT_type */,
    c_Real* s_body /* 0D_NOT_real */,
    c_Bool* is_time_coords /* 0D_NOT_logical */);
void save_a_bunch_step(
    EleProxy& ele,
    BunchProxy& bunch,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt,
    std::optional<double> s_body = std::nullopt,
    std::optional<bool> is_time_coords = std::nullopt);
extern "C" void fortran_save_a_step(
    void* track /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& local_ref_frame /* 0D_NOT_logical */,
    void* orb /* 0D_NOT_type */,
    c_Real& s_rel /* 0D_NOT_real */,
    c_Bool* save_field /* 0D_NOT_logical */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_Real* rf_time /* 0D_NOT_real */,
    void* strong_beam /* 0D_NOT_type */);
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
    void* ele /* 0D_NOT_type */,
    c_Real& dg /* 0D_NOT_real */,
    c_Real& b1 /* 0D_NOT_real */,
    void* param /* 0D_NOT_type */,
    c_Int& n_step /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
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
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& include_image /* 0D_NOT_logical */,
    c_Real& t_now /* 0D_NOT_real */,
    c_Real& dt_step /* 0D_NOT_real */,
    c_Real& dt_next /* 0D_NOT_real */,
    void* sc_field /* 1D_ALLOC_type */);
double sc_adaptive_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool include_image,
    double t_now,
    double dt_step,
    EmFieldProxyAlloc1D& sc_field);
extern "C" void fortran_sc_step(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& include_image /* 0D_NOT_logical */,
    c_Real& t_end /* 0D_NOT_real */,
    void* sc_field /* 1D_ALLOC_type */,
    c_Int& n_emit /* 0D_NOT_integer */);
int sc_step(
    BunchProxy& bunch,
    EleProxy& ele,
    bool include_image,
    double t_end,
    EmFieldProxyAlloc1D& sc_field);

// Skipped unusable routine seq_ele_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine seq_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_set_active_fixer(
    void* fixer /* 0D_NOT_type */,
    c_Bool* is_on /* 0D_NOT_logical */,
    void* orbit /* 0D_NOT_type */);
CoordProxy set_active_fixer(
    EleProxy& fixer,
    std::optional<bool> is_on = std::nullopt);

// Skipped unusable routine set_branch_and_ele_for_omp:
// Untranslated type: LatPointerProxy (1D_ALLOC_type)
extern "C" void fortran_set_custom_attribute_name(
    c_Char custom_name /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Int* custom_index /* 0D_NOT_integer */);
bool set_custom_attribute_name(
    std::string custom_name,
    std::optional<int> custom_index = std::nullopt);
extern "C" void fortran_set_ele_attribute(
    void* ele /* 0D_NOT_type */,
    c_Char set_string /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Bool* set_lords /* 0D_NOT_logical */,
    c_Int& err_id /* 0D_NOT_integer */);
struct SetEleAttribute {
  bool err_flag;
  int err_id;
};
SetEleAttribute set_ele_attribute(
    EleProxy& ele,
    std::string set_string,
    std::optional<bool> err_print_flag = std::nullopt,
    std::optional<bool> set_lords = std::nullopt);
extern "C" void fortran_set_ele_defaults(
    void* ele /* 0D_NOT_type */,
    c_Bool* do_allocate /* 0D_NOT_logical */);
void set_ele_defaults(
    EleProxy& ele,
    std::optional<bool> do_allocate = std::nullopt);
extern "C" void fortran_set_ele_name(
    void* ele /* 0D_NOT_type */,
    c_Char name /* 0D_NOT_character */);
void set_ele_name(EleProxy& ele, std::string name);
extern "C" void fortran_set_ele_real_attribute(
    void* ele /* 0D_NOT_type */,
    c_Char attrib_name /* 0D_NOT_character */,
    c_Real& value /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* err_print_flag /* 0D_NOT_logical */);
bool set_ele_real_attribute(
    EleProxy& ele,
    std::string attrib_name,
    double value,
    std::optional<bool> err_print_flag = std::nullopt);
extern "C" void fortran_set_ele_status_stale(
    void* ele /* 0D_NOT_type */,
    c_Int& status_group /* 0D_NOT_integer */,
    c_Bool& set_slaves /* 0D_NOT_logical */);
struct SetEleStatusStale {
  EleProxy ele;
  int status_group;
  bool set_slaves;
};
SetEleStatusStale set_ele_status_stale();
extern "C" bool fortran_set_emit_from_beam_init(
    void* beam_init_in /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int& species /* 0D_NOT_integer */,
    void* modes /* 0D_NOT_type */,
    c_Bool* err_flag /* 0D_NOT_logical */,
    void* beam_init_set /* 0D_NOT_type */);
void set_emit_from_beam_init(
    BeamInitProxy& beam_init_in,
    EleProxy& ele,
    int species,
    optional_ref<NormalModesProxy> modes,
    std::optional<bool> err_flag,
    BeamInitProxy& beam_init_set);

// Skipped unusable routine set_flags_for_changed_all_attribute:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" void fortran_set_flags_for_changed_integer_attribute(
    void* ele /* 0D_NOT_type */,
    c_Int& attrib /* 0D_NOT_integer */,
    c_Bool* set_dependent /* 0D_NOT_logical */);
void set_flags_for_changed_integer_attribute(
    EleProxy& ele,
    int attrib,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_lat_attribute(
    void* lat /* 0D_NOT_type */,
    c_Bool* set_dependent /* 0D_NOT_logical */);
void set_flags_for_changed_lat_attribute(
    LatProxy& lat,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_logical_attribute(
    void* ele /* 0D_NOT_type */,
    c_Bool& attrib /* 0D_NOT_logical */,
    c_Bool* set_dependent /* 0D_NOT_logical */);
void set_flags_for_changed_logical_attribute(
    EleProxy& ele,
    bool attrib,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_flags_for_changed_real_attribute(
    void* ele /* 0D_NOT_type */,
    c_Real* attrib /* 0D_NOT_real */,
    c_Bool* set_dependent /* 0D_NOT_logical */);
void set_flags_for_changed_real_attribute(
    EleProxy& ele,
    optional_ref<double> attrib = std::nullopt,
    std::optional<bool> set_dependent = std::nullopt);
extern "C" void fortran_set_fringe_on_off(
    c_Real& fringe_at /* 0D_NOT_real */,
    c_Int& ele_end /* 0D_NOT_integer */,
    c_Int& on_or_off /* 0D_NOT_integer */);
void set_fringe_on_off(double fringe_at, int ele_end, int on_or_off);
extern "C" void fortran_set_lords_status_stale(
    void* ele /* 0D_NOT_type */,
    c_Int& stat_group /* 0D_NOT_integer */,
    c_Bool* control_bookkeeping /* 0D_NOT_logical */,
    c_Int* flag /* 0D_NOT_integer */);
void set_lords_status_stale(
    EleProxy& ele,
    int stat_group,
    std::optional<bool> control_bookkeeping = std::nullopt,
    std::optional<int> flag = std::nullopt);
extern "C" void fortran_set_on_off(
    c_Int& key /* 0D_NOT_integer */,
    void* lat /* 0D_NOT_type */,
    c_Int& switch_ /* 0D_NOT_integer */,
    void* orb /* 1D_ALLOC_type */,
    c_Bool* use_ref_orb /* 0D_NOT_logical */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    void* saved_values /* 1D_ALLOC_real */,
    c_Char attribute /* 0D_NOT_character */,
    c_Int* set_val /* 0D_NOT_integer */);
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
    void* orbit /* 1D_ALLOC_type */,
    c_Int& n1 /* 0D_NOT_integer */,
    c_Int& n2 /* 0D_NOT_integer */,
    c_Int* ix_noset /* 0D_NOT_integer */);
CoordProxyAlloc1D set_orbit_to_zero(
    int n1,
    int n2,
    std::optional<int> ix_noset = std::nullopt);
extern "C" void fortran_set_ptc(
    c_Real* e_tot /* 0D_NOT_real */,
    c_Int* particle /* 0D_NOT_integer */,
    c_Int* taylor_order /* 0D_NOT_integer */,
    c_Int* integ_order /* 0D_NOT_integer */,
    c_Int* n_step /* 0D_NOT_integer */,
    c_Bool* no_cavity /* 0D_NOT_logical */,
    c_Bool* force_init /* 0D_NOT_logical */);
void set_ptc(
    std::optional<double> e_tot = std::nullopt,
    std::optional<int> particle = std::nullopt,
    std::optional<int> taylor_order = std::nullopt,
    std::optional<int> integ_order = std::nullopt,
    std::optional<int> n_step = std::nullopt,
    std::optional<bool> no_cavity = std::nullopt,
    std::optional<bool> force_init = std::nullopt);
extern "C" void fortran_set_ptc_base_state(
    c_Char component /* 0D_NOT_character */,
    c_Bool& set_val /* 0D_NOT_logical */,
    c_Bool& old_val /* 0D_NOT_logical */);
bool set_ptc_base_state(std::string component, bool set_val);
extern "C" void fortran_set_ptc_com_pointers();
void set_ptc_com_pointers();
extern "C" void fortran_set_ptc_quiet(
    c_Int& channel /* 0D_NOT_integer */,
    c_Bool& set /* 0D_NOT_logical */,
    c_Int& old_val /* 0D_NOT_integer */);
void set_ptc_quiet(int channel, bool set, int old_val);
extern "C" void fortran_set_ptc_verbose(c_Bool& on /* 0D_NOT_logical */);
void set_ptc_verbose(bool on);
extern "C" void fortran_set_pwd_ele(
    void* lat /* 0D_NOT_type */,
    void* mode0 /* 0D_NOT_type */,
    c_Real& inductance /* 0D_NOT_real */);
void set_pwd_ele(LatProxy& lat, NormalModesProxy& mode0, double inductance);
extern "C" void fortran_set_status_flags(
    void* bookkeeping_state /* 0D_NOT_type */,
    c_Int& stat /* 0D_NOT_integer */);
BookkeepingStateProxy set_status_flags(int stat);

// Skipped unusable routine set_tune:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine set_tune_via_group_knobs:
// Routine in configuration skip list
extern "C" void fortran_set_twiss(
    void* branch /* 0D_NOT_type */,
    void* twiss_ele /* 0D_NOT_type */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Bool& match_deta_ds /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
void set_twiss(
    BranchProxy& branch,
    EleProxy& twiss_ele,
    int ix_ele,
    bool match_deta_ds,
    bool err_flag,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_set_z_tune(
    void* branch /* 0D_NOT_type */,
    c_Real& z_tune /* 0D_NOT_real */,
    c_Bool& ok /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
bool set_z_tune(
    BranchProxy& branch,
    double z_tune,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_settable_dep_var_bookkeeping(
    void* ele /* 0D_NOT_type */);
void settable_dep_var_bookkeeping(EleProxy& ele);
extern "C" void fortran_setup_high_energy_space_charge_calc(
    c_Bool& calc_on /* 0D_NOT_logical */,
    void* branch /* 0D_NOT_type */,
    c_Real& n_part /* 0D_NOT_real */,
    void* mode /* 0D_NOT_type */,
    void* closed_orb /* 1D_ALLOC_type */);
void setup_high_energy_space_charge_calc(
    bool calc_on,
    BranchProxy& branch,
    double n_part,
    NormalModesProxy& mode,
    optional_ref<CoordProxyAlloc1D> closed_orb = std::nullopt);

// Skipped unusable routine sfft:
// Routine in configuration skip list
extern "C" void fortran_sigma_mat_ptc_to_bmad(
    c_RealArr sigma_mat_ptc /* 2D_NOT_real */,
    c_Real& beta0 /* 0D_NOT_real */,
    c_RealArr sigma_mat_bmad /* 2D_NOT_real */);
FixedArray2D<Real, 6, 6> sigma_mat_ptc_to_bmad(
    FixedArray2D<Real, 6, 6> sigma_mat_ptc,
    double beta0);
extern "C" bool fortran_significant_difference(
    c_Real& value1 /* 0D_NOT_real */,
    c_Real& value2 /* 0D_NOT_real */,
    c_Real* abs_tol /* 0D_NOT_real */,
    c_Real* rel_tol /* 0D_NOT_real */,
    c_Bool& is_different /* 0D_NOT_logical */);
void significant_difference(
    double value1,
    double value2,
    std::optional<double> abs_tol,
    std::optional<double> rel_tol,
    bool is_different);
extern "C" bool fortran_skip_ele_blender(
    void* ele /* 0D_NOT_type */,
    c_Bool& skip /* 0D_NOT_logical */);
void skip_ele_blender(EleProxy& ele, bool skip);
extern "C" void fortran_slice_lattice(
    void* lat /* 0D_NOT_type */,
    c_Char ele_list /* 0D_NOT_character */,
    c_Bool& error /* 0D_NOT_logical */,
    c_Bool* do_bookkeeping /* 0D_NOT_logical */);
bool slice_lattice(
    LatProxy& lat,
    std::string ele_list,
    std::optional<bool> do_bookkeeping = std::nullopt);

// Skipped unusable routine sliced_eles_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_soft_quadrupole_edge_kick(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Int& particle_at /* 0D_NOT_integer */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void soft_quadrupole_edge_kick(
    EleProxy& ele,
    LatParamProxy& param,
    int particle_at,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_sol_quad_mat6_calc(
    c_Real& ks_in /* 0D_NOT_real */,
    c_Real& k1_in /* 0D_NOT_real */,
    c_Real& tilt /* 0D_NOT_real */,
    c_Real& length /* 0D_NOT_real */,
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void sol_quad_mat6_calc(
    double ks_in,
    double k1_in,
    double tilt,
    double length,
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);

// Skipped unusable routine solenoid_track_and_mat:
// Variable inout sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_solve_psi_adaptive(
    c_Real& t0 /* 0D_NOT_real */,
    c_Real& t1 /* 0D_NOT_real */,
    c_Real& p0 /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    c_Real& p1 /* 0D_NOT_real */);
double solve_psi_adaptive(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args);
extern "C" void fortran_solve_psi_fixed_steps(
    c_Real& t0 /* 0D_NOT_real */,
    c_Real& t1 /* 0D_NOT_real */,
    c_Real& p0 /* 0D_NOT_real */,
    c_RealArr args /* 1D_NOT_real */,
    void* t /* 1D_ALLOC_real */,
    void* p /* 1D_ALLOC_real */);
struct SolvePsiFixedSteps {
  RealAlloc1D t;
  RealAlloc1D p;
};
SolvePsiFixedSteps solve_psi_fixed_steps(
    double t0,
    double t1,
    double p0,
    FixedArray1D<Real, 8> args);
extern "C" void fortran_sort_complex_taylor_terms(
    void* complex_taylor_in /* 0D_NOT_type */,
    void* complex_taylor_sorted /* 0D_NOT_type */);
ComplexTaylorProxy sort_complex_taylor_terms(
    ComplexTaylorProxy& complex_taylor_in);

// Skipped unusable routine sort_universal_terms:
// Untranslated type: UniversalTaylorProxy (0D_NOT_type)
// Untranslated type: UniversalTaylorProxy (0D_NOT_type)

// Skipped unusable routine space_charge_3d:
// Untranslated type: Mesh3dProxy (0D_NOT_type)
// Variable in sized array: image_efield(:,:,:,:) 4D_ALLOC_real

// Skipped unusable routine space_charge_cathodeimages:
// Untranslated type: Mesh3dProxy (0D_NOT_type)

// Skipped unusable routine space_charge_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine space_charge_freespace:
// Untranslated type: Mesh3dProxy (0D_NOT_type)

// Skipped unusable routine space_charge_rectpipe:
// Untranslated type: Mesh3dProxy (0D_NOT_type)

// Skipped unusable routine spin_axis_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine spin_concat_linear_maps:
// Routine in configuration skip list

// Skipped unusable routine spin_depolarization_rate:
// Untranslated type: SpinMatchingProxy (1D_ALLOC_type)
extern "C" bool fortran_spin_dn_dpz_from_mat8(
    c_RealArr mat_1turn /* 2D_NOT_real */,
    c_RealArr dn_dpz_partial /* 2D_NOT_real */,
    c_Bool& error /* 0D_NOT_logical */,
    c_RealArr dn_dpz /* 1D_NOT_real */);
bool spin_dn_dpz_from_mat8(
    FixedArray2D<Real, 8, 8> mat_1turn,
    std::optional<FixedArray2D<Real, 3, 3>> dn_dpz_partial,
    FixedArray1D<Real, 3> dn_dpz);
extern "C" bool fortran_spin_dn_dpz_from_qmap(
    c_RealArr orb_mat /* 2D_NOT_real */,
    c_RealArr q_map /* 2D_NOT_real */,
    c_RealArr dn_dpz_partial /* 2D_NOT_real */,
    c_RealArr dn_dpz_partial2 /* 2D_NOT_real */,
    c_Bool& error /* 0D_NOT_logical */,
    c_RealArr n0 /* 1D_NOT_real */,
    c_RealArr dn_dpz /* 1D_NOT_real */);
bool spin_dn_dpz_from_qmap(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> q_map,
    FixedArray2D<Real, 3, 3> dn_dpz_partial,
    FixedArray2D<Real, 3, 3> dn_dpz_partial2,
    std::optional<FixedArray1D<Real, 3>> n0,
    FixedArray1D<Real, 3> dn_dpz);

// Skipped unusable routine spin_eigen_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_spin_map1_normalize(c_RealArr spin1 /* 2D_NOT_real */);
void spin_map1_normalize(FixedArray2D<Real, 4, 7> spin1);
extern "C" void fortran_spin_mat8_resonance_strengths(
    c_ComplexArr orb_evec /* 1D_NOT_complex */,
    c_RealArr mat8 /* 2D_NOT_real */,
    c_Real& xi_sum /* 0D_NOT_real */,
    c_Real& xi_diff /* 0D_NOT_real */);
struct SpinMat8ResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
SpinMat8ResonanceStrengths spin_mat8_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 6, 6> mat8);
extern "C" void fortran_spin_mat_to_eigen(
    c_RealArr orb_mat /* 2D_NOT_real */,
    c_RealArr spin_map /* 2D_NOT_real */,
    c_ComplexArr orb_eval /* 1D_NOT_complex */,
    c_ComplexArr orb_evec /* 2D_NOT_complex */,
    c_RealArr n0 /* 1D_NOT_real */,
    c_ComplexArr spin_evec /* 2D_NOT_complex */,
    c_Bool& error /* 0D_NOT_logical */);
struct SpinMatToEigen {
  FixedArray1D<Complex, 6> orb_eval;
  FixedArray2D<Complex, 6, 6> orb_evec;
  FixedArray1D<Real, 3> n0;
  FixedArray2D<Complex, 6, 3> spin_evec;
  bool error;
};
SpinMatToEigen spin_mat_to_eigen(
    FixedArray2D<Real, 6, 6> orb_mat,
    FixedArray2D<Real, 4, 7> spin_map);

// Skipped unusable routine spin_matching_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_spin_omega(
    void* field /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Int& sign_z_vel /* 0D_NOT_integer */,
    c_Bool* phase_space_coords /* 0D_NOT_logical */,
    c_RealArr omega /* 1D_NOT_real */);
void spin_omega(
    EmFieldProxy& field,
    CoordProxy& orbit,
    int sign_z_vel,
    optional_ref<bool> phase_space_coords,
    FixedArray1D<Real, 3> omega);

// Skipped unusable routine spin_orbit_map1_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine spin_polar_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_spin_quat_resonance_strengths(
    c_ComplexArr orb_evec /* 1D_NOT_complex */,
    c_RealArr spin_q /* 2D_NOT_real */,
    c_Real& xi_sum /* 0D_NOT_real */,
    c_Real& xi_diff /* 0D_NOT_real */);
struct SpinQuatResonanceStrengths {
  double xi_sum;
  double xi_diff;
};
SpinQuatResonanceStrengths spin_quat_resonance_strengths(
    FixedArray1D<Complex, 6> orb_evec,
    FixedArray2D<Real, 4, 7> spin_q);
extern "C" bool fortran_spin_taylor_to_linear(
    void* spin_taylor /* 1D_NOT_type */,
    c_Bool& normalize /* 0D_NOT_logical */,
    c_RealArr dref_orb /* 1D_NOT_real */,
    c_Bool& is_on /* 0D_NOT_logical */,
    c_RealArr spin_map1 /* 2D_NOT_real */);
void spin_taylor_to_linear(
    FixedArray1D<TaylorProxy, 4> spin_taylor,
    bool normalize,
    FixedArray1D<Real, 6> dref_orb,
    bool is_on,
    FixedArray2D<Real, 4, 7> spin_map1);
extern "C" bool fortran_spinor_to_polar(
    c_ComplexArr spinor /* 1D_NOT_complex */,
    void* polar /* 0D_NOT_type */);
void spinor_to_polar(FixedArray1D<Complex, 2> spinor, SpinPolarProxy& polar);
extern "C" bool fortran_spinor_to_vec(
    c_ComplexArr spinor /* 1D_NOT_complex */,
    c_RealArr vec /* 1D_NOT_real */);
void spinor_to_vec(FixedArray1D<Complex, 2> spinor, FixedArray1D<Real, 3> vec);
extern "C" void fortran_spline_fit_orbit(
    void* start_orb /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_RealArr spline_x /* 1D_NOT_real */,
    c_RealArr spline_y /* 1D_NOT_real */);
void spline_fit_orbit(
    CoordProxy& start_orb,
    CoordProxy& end_orb,
    FixedArray1D<Real, 4> spline_x,
    FixedArray1D<Real, 4> spline_y);

// Skipped unusable routine split_expression_string:
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_split_lat(
    void* lat /* 0D_NOT_type */,
    c_Real& s_split /* 0D_NOT_real */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Int& ix_split /* 0D_NOT_integer */,
    c_Bool& split_done /* 0D_NOT_logical */,
    c_Bool* add_suffix /* 0D_NOT_logical */,
    c_Bool* check_sanity /* 0D_NOT_logical */,
    c_Bool* save_null_drift /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* choose_max /* 0D_NOT_logical */,
    c_Int* ix_insert /* 0D_NOT_integer */);
struct SplitLat {
  int ix_split;
  bool split_done;
  bool err_flag;
};
SplitLat split_lat(
    LatProxy& lat,
    double s_split,
    int ix_branch,
    std::optional<bool> add_suffix = std::nullopt,
    std::optional<bool> check_sanity = std::nullopt,
    std::optional<bool> save_null_drift = std::nullopt,
    std::optional<bool> choose_max = std::nullopt,
    std::optional<int> ix_insert = std::nullopt);
extern "C" void fortran_sprint_spin_taylor_map(
    void* ele /* 0D_NOT_type */,
    c_RealArr start_orbit /* 1D_NOT_real */);
void sprint_spin_taylor_map(
    EleProxy& ele,
    std::optional<FixedArray1D<Real, 6>> start_orbit = std::nullopt);
extern "C" void fortran_sr_longitudinal_wake_particle(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void sr_longitudinal_wake_particle(EleProxy& ele, CoordProxy& orbit);
extern "C" void fortran_sr_transverse_wake_particle(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void sr_transverse_wake_particle(EleProxy& ele, CoordProxy& orbit);
extern "C" void fortran_sr_z_long_wake(
    void* ele /* 0D_NOT_type */,
    void* bunch /* 0D_NOT_type */,
    c_Real& z_ave /* 0D_NOT_real */);
void sr_z_long_wake(EleProxy& ele, BunchProxy& bunch, double z_ave);
extern "C" void fortran_srdt_calc(
    void* lat /* 0D_NOT_type */,
    void* srdt_sums /* 0D_NOT_type */,
    c_Int& order /* 0D_NOT_integer */,
    c_Int* n_slices_gen_opt /* 0D_NOT_integer */,
    c_Int* n_slices_sxt_opt /* 0D_NOT_integer */,
    void* per_ele_out /* 1D_ALLOC_type */);
SummationRdtProxy srdt_calc(
    LatProxy& lat,
    int order,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    optional_ref<SummationRdtProxyAlloc1D> per_ele_out = std::nullopt);

// Skipped unusable routine srdt_calc_with_cache:
// Variable inout sized array: cache(:,:,:) 3D_ALLOC_complex
extern "C" void fortran_srdt_lsq_solution(
    void* lat /* 0D_NOT_type */,
    void* var_indexes /* 1D_ALLOC_integer */,
    void* ls_soln /* 1D_ALLOC_real */,
    c_Int* n_slices_gen_opt /* 0D_NOT_integer */,
    c_Int* n_slices_sxt_opt /* 0D_NOT_integer */,
    c_Real* chrom_set_x_opt /* 0D_NOT_real */,
    c_Real* chrom_set_y_opt /* 0D_NOT_real */,
    c_RealArr weight_in /* 1D_NOT_real */);
RealAlloc1D srdt_lsq_solution(
    LatProxy& lat,
    IntAlloc1D& var_indexes,
    std::optional<int> n_slices_gen_opt = std::nullopt,
    std::optional<int> n_slices_sxt_opt = std::nullopt,
    std::optional<double> chrom_set_x_opt = std::nullopt,
    std::optional<double> chrom_set_y_opt = std::nullopt,
    std::optional<FixedArray1D<Real, 10>> weight_in = std::nullopt);

// Skipped unusable routine stack_file_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_start_branch_at(
    void* lat /* 0D_NOT_type */,
    c_Char ele_start /* 0D_NOT_character */,
    c_Bool& move_end_marker /* 0D_NOT_logical */,
    c_Bool& error /* 0D_NOT_logical */);
bool start_branch_at(
    LatProxy& lat,
    std::string ele_start,
    bool move_end_marker);
extern "C" bool fortran_stream_ele_end(
    c_Int& physical_end /* 0D_NOT_integer */,
    c_Int& ele_orientation /* 0D_NOT_integer */,
    c_Int& stream_end /* 0D_NOT_integer */);
void stream_ele_end(int physical_end, int ele_orientation, int stream_end);
extern "C" void fortran_string_attrib(
    c_Char attrib_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Char attrib_value /* 0D_NOT_character */);
std::string string_attrib(std::string attrib_name, EleProxy& ele);
extern "C" void fortran_strong_beam_sigma_calc(
    void* ele /* 0D_NOT_type */,
    c_Real& s_pos /* 0D_NOT_real */,
    c_RealArr sigma /* 1D_NOT_real */,
    c_Real& bbi_const /* 0D_NOT_real */,
    c_RealArr dsigma_ds /* 1D_NOT_real */);
struct StrongBeamSigmaCalc {
  FixedArray1D<Real, 2> sigma;
  double bbi_const;
  FixedArray1D<Real, 2> dsigma_ds;
};
StrongBeamSigmaCalc strong_beam_sigma_calc(EleProxy& ele, double s_pos);
extern "C" bool fortran_strong_beam_strength(
    void* ele /* 0D_NOT_type */,
    c_Real& strength /* 0D_NOT_real */);
void strong_beam_strength(EleProxy& ele, double strength);

// Skipped unusable routine strong_beam_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine summation_rdt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_curvature_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_displacement_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_displacement_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_surface_grid_displacement(
    void* ele /* 0D_NOT_type */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& z /* 0D_NOT_real */,
    c_RealArr dz_dxy /* 1D_NOT_real */,
    c_Bool* extend_grid /* 0D_NOT_logical */);
void surface_grid_displacement(
    EleProxy& ele,
    double x,
    double y,
    bool err_flag,
    double z,
    std::optional<FixedArray1D<Real, 2>> dz_dxy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt);

// Skipped unusable routine surface_h_misalign_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_h_misalign_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_segmented_pt_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine surface_segmented_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine switch_attrib_value_name:
// Variable-sized out character array: name_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_symp_lie_bmad(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    void* track /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_Bool* offset_ele /* 0D_NOT_logical */);
TrackProxy symp_lie_bmad(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<bool> offset_ele = std::nullopt);
extern "C" void fortran_t6_to_b123(
    c_RealArr t6 /* 2D_NOT_real */,
    c_RealArr abz_tunes /* 1D_NOT_real */,
    c_RealArr B1 /* 2D_NOT_real */,
    c_RealArr B2 /* 2D_NOT_real */,
    c_RealArr B3 /* 2D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct T6ToB123 {
  FixedArray2D<Real, 6, 6> B1;
  FixedArray2D<Real, 6, 6> B2;
  FixedArray2D<Real, 6, 6> B3;
  bool err_flag;
};
T6ToB123 t6_to_b123(
    FixedArray2D<Real, 6, 6> t6,
    FixedArray1D<Real, 3> abz_tunes);
extern "C" void fortran_taper_mag_strengths(
    void* lat /* 0D_NOT_type */,
    void* ref_lat /* 0D_NOT_type */,
    c_Char except /* 0D_NOT_character */,
    c_Bool* err_flag /* 0D_NOT_logical */);
void taper_mag_strengths(
    LatProxy& lat,
    optional_ref<LatProxy> ref_lat = std::nullopt,
    std::optional<std::string> except = std::nullopt,
    optional_ref<bool> err_flag = std::nullopt);
extern "C" void fortran_target_min_max_calc(
    c_RealArr r_corner1 /* 1D_NOT_real */,
    c_RealArr r_corner2 /* 1D_NOT_real */,
    c_Real& y_min /* 0D_NOT_real */,
    c_Real& y_max /* 0D_NOT_real */,
    c_Real& phi_min /* 0D_NOT_real */,
    c_Real& phi_max /* 0D_NOT_real */,
    c_Bool* initial /* 0D_NOT_logical */);
void target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double y_min,
    double y_max,
    double phi_min,
    double phi_max,
    std::optional<bool> initial = std::nullopt);

// Skipped unusable routine target_point_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_target_rot_mats(
    c_RealArr r_center /* 1D_NOT_real */,
    c_RealArr w_to_target /* 2D_NOT_real */,
    c_RealArr w_to_ele /* 2D_NOT_real */);
struct TargetRotMats {
  FixedArray2D<Real, 3, 3> w_to_target;
  FixedArray2D<Real, 3, 3> w_to_ele;
};
TargetRotMats target_rot_mats(FixedArray1D<Real, 3> r_center);
extern "C" void fortran_taylor_equal_taylor(
    void* taylor1 /* 0D_NOT_type */,
    void* taylor2 /* 0D_NOT_type */);
TaylorProxy taylor_equal_taylor(TaylorProxy& taylor2);
extern "C" void fortran_taylor_inverse(
    void* taylor_in /* 1D_ALLOC_type */,
    void* taylor_inv /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */);
struct TaylorInverse {
  TaylorProxyAlloc1D taylor_inv;
  bool err;
};
TaylorInverse taylor_inverse(TaylorProxyAlloc1D& taylor_in);

// Skipped unusable routine taylor_minus_taylor:
// Routine in configuration skip list

// Skipped unusable routine taylor_plus_taylor:
// Routine in configuration skip list
extern "C" void fortran_taylor_propagate1(
    void* orb_taylor /* 1D_ALLOC_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ref_in /* 0D_NOT_type */,
    void* spin_taylor /* 1D_ALLOC_type */);
bool taylor_propagate1(
    TaylorProxyAlloc1D& orb_taylor,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<CoordProxy> ref_in = std::nullopt,
    optional_ref<TaylorProxyAlloc1D> spin_taylor = std::nullopt);

// Skipped unusable routine taylor_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine taylor_term_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine taylor_to_genfield:
// Untranslated type: GenfieldProxy (0D_NOT_type)

// Skipped unusable routine taylor_to_mad_map:
// Untranslated type: MadEnergyProxy (0D_NOT_type)
// Untranslated type: MadMapProxy (0D_NOT_type)

// Skipped unusable routine taylor_to_real_8:
// Untranslated type: Real8Proxy (1D_NOT_type)
extern "C" void fortran_taylors_equal_taylors(
    void* taylor1 /* 1D_ALLOC_type */,
    void* taylor2 /* 1D_ALLOC_type */);
TaylorProxyAlloc1D taylors_equal_taylors(TaylorProxyAlloc1D& taylor2);
extern "C" void fortran_tilt_coords(
    c_Real& tilt_val /* 0D_NOT_real */,
    void* coord /* 1D_ALLOC_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void tilt_coords(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_tilt_coords_photon(
    c_Real& tilt_val /* 0D_NOT_real */,
    void* coord /* 1D_ALLOC_real */,
    c_RealArr w_mat /* 2D_NOT_real */);
void tilt_coords_photon(
    double tilt_val,
    RealAlloc1D& coord,
    std::optional<FixedArray2D<Real, 3, 3>> w_mat = std::nullopt);
extern "C" void fortran_tilt_mat6(
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Real& tilt /* 0D_NOT_real */);
void tilt_mat6(FixedArray2D<Real, 6, 6> mat6, double tilt);

// Skipped unusable routine time_runge_kutta_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine time_runge_kutta_periodic_kick_hook_def:
// Routine in configuration skip list
extern "C" void fortran_to_eta_reading(
    void* eta_actual /* 1D_ALLOC_real */,
    void* ele /* 0D_NOT_type */,
    c_Int& axis /* 0D_NOT_integer */,
    c_Bool& add_noise /* 0D_NOT_logical */,
    c_Real& reading /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct ToEtaReading {
  double reading;
  bool err;
};
ToEtaReading to_eta_reading(
    RealAlloc1D& eta_actual,
    EleProxy& ele,
    int axis,
    bool add_noise);
extern "C" void fortran_to_fieldmap_coords(
    void* ele /* 0D_NOT_type */,
    void* local_orb /* 0D_NOT_type */,
    c_Real& s_body /* 0D_NOT_real */,
    c_Int& ele_anchor_pt /* 0D_NOT_integer */,
    c_RealArr r0 /* 1D_NOT_real */,
    c_Bool& curved_ref_frame /* 0D_NOT_logical */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& cos_ang /* 0D_NOT_real */,
    c_Real& sin_ang /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void to_fieldmap_coords(
    EleProxy& ele,
    CoordProxy& local_orb,
    double s_body,
    int ele_anchor_pt,
    FixedArray1D<Real, 3> r0,
    bool curved_ref_frame,
    double x,
    double y,
    double z,
    double cos_ang,
    double sin_ang,
    bool err_flag);
extern "C" void fortran_to_orbit_reading(
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int& axis /* 0D_NOT_integer */,
    c_Bool& add_noise /* 0D_NOT_logical */,
    c_Real& reading /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct ToOrbitReading {
  double reading;
  bool err;
};
ToOrbitReading to_orbit_reading(
    CoordProxy& orb,
    EleProxy& ele,
    int axis,
    bool add_noise);
extern "C" void fortran_to_phase_and_coupling_reading(
    void* ele /* 0D_NOT_type */,
    c_Bool& add_noise /* 0D_NOT_logical */,
    void* reading /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
struct ToPhaseAndCouplingReading {
  BpmPhaseCouplingProxy reading;
  bool err;
};
ToPhaseAndCouplingReading to_phase_and_coupling_reading(
    EleProxy& ele,
    bool add_noise);
extern "C" bool fortran_to_photon_angle_coords(
    void* orb_in /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* orb_out /* 0D_NOT_type */);
CoordProxy to_photon_angle_coords(CoordProxy& orb_in, EleProxy& ele);
extern "C" void fortran_to_surface_coords(
    void* lab_orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* surface_orbit /* 0D_NOT_type */);
CoordProxy to_surface_coords(CoordProxy& lab_orbit, EleProxy& ele);
extern "C" void fortran_touschek_lifetime(
    void* mode /* 0D_NOT_type */,
    c_Real& Tl /* 0D_NOT_real */,
    void* lat /* 0D_NOT_type */);
double touschek_lifetime(NormalModesProxy& mode, LatProxy& lat);

// Skipped unusable routine touschek_lifetime_ele_by_ele:
// Untranslated type: MomentumApertureProxy (1D_ALLOC_type)

// Skipped unusable routine touschek_lifetime_with_aperture:
// Untranslated type: MomentumApertureProxy (1D_ALLOC_type)
extern "C" void fortran_touschek_rate1(
    void* mode /* 0D_NOT_type */,
    c_Real& rate /* 0D_NOT_real */,
    void* lat /* 0D_NOT_type */,
    c_Int* ix /* 0D_NOT_integer */,
    c_Real* s /* 0D_NOT_real */);
double touschek_rate1(
    NormalModesProxy& mode,
    LatProxy& lat,
    std::optional<int> ix = std::nullopt,
    std::optional<double> s = std::nullopt);
extern "C" void fortran_touschek_rate1_zap(
    void* mode /* 0D_NOT_type */,
    c_Real& rate /* 0D_NOT_real */,
    void* lat /* 0D_NOT_type */,
    c_Int* ix /* 0D_NOT_integer */,
    c_Real* s /* 0D_NOT_real */);
void touschek_rate1_zap(
    NormalModesProxy& mode,
    double rate,
    LatProxy& lat,
    optional_ref<int> ix = std::nullopt,
    optional_ref<double> s = std::nullopt);
extern "C" void fortran_track1(
    void* start_orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    void* track /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* ignore_radiation /* 0D_NOT_logical */,
    c_Bool* make_map1 /* 0D_NOT_logical */,
    c_Bool* init_to_edge /* 0D_NOT_logical */);
struct Track1 {
  CoordProxy end_orb;
  bool err_flag;
};
Track1 track1(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<TrackProxy> track = std::nullopt,
    std::optional<bool> ignore_radiation = std::nullopt,
    std::optional<bool> make_map1 = std::nullopt,
    std::optional<bool> init_to_edge = std::nullopt);
extern "C" void fortran_track1_beam(
    void* beam /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    void* centroid /* 1D_ALLOC_type */,
    c_Int* direction /* 0D_NOT_integer */);
bool track1_beam(
    BeamProxy& beam,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt);
extern "C" void fortran_track1_bmad(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
struct Track1Bmad {
  bool err_flag;
  TrackProxy track;
};
Track1Bmad track1_bmad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_bmad_photon(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool track1_bmad_photon(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_bunch(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    void* centroid /* 1D_ALLOC_type */,
    c_Int* direction /* 0D_NOT_integer */,
    void* bunch_track /* 0D_NOT_type */);
bool track1_bunch(
    BunchProxy& bunch,
    EleProxy& ele,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_csr(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* centroid /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Real* s_start /* 0D_NOT_real */,
    c_Real* s_end /* 0D_NOT_real */,
    void* bunch_track /* 0D_NOT_type */);
bool track1_bunch_csr(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_csr3d(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* centroid /* 1D_ALLOC_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Real* s_start /* 0D_NOT_real */,
    c_Real* s_end /* 0D_NOT_real */,
    void* bunch_track /* 0D_NOT_type */);
bool track1_bunch_csr3d(
    BunchProxy& bunch,
    EleProxy& ele,
    CoordProxyAlloc1D& centroid,
    std::optional<double> s_start = std::nullopt,
    std::optional<double> s_end = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_bunch_hom(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int* direction /* 0D_NOT_integer */,
    void* bunch_track /* 0D_NOT_type */);
void track1_bunch_hom(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);

// Skipped unusable routine track1_bunch_hook_def:
// Routine in configuration skip list
extern "C" void fortran_track1_bunch_space_charge(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* track_to_same_s /* 0D_NOT_logical */,
    void* bunch_track /* 0D_NOT_type */);
bool track1_bunch_space_charge(
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<bool> track_to_same_s = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track1_crystal(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_crystal(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);

// Skipped unusable routine track1_custom_def:
// Routine in configuration skip list
extern "C" void fortran_track1_diffraction_plate_or_mask(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_diffraction_plate_or_mask(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_high_energy_space_charge(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_high_energy_space_charge(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_lens(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_lens(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_linear(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */);
void track1_linear(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_lr_wake(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */);
void track1_lr_wake(BunchProxy& bunch, EleProxy& ele);
extern "C" void fortran_track1_mad(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */);
void track1_mad(CoordProxy& orbit, EleProxy& ele, LatParamProxy& param);
extern "C" void fortran_track1_mirror(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_mirror(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_mosaic_crystal(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_mosaic_crystal(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);
extern "C" void fortran_track1_multilayer_mirror(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_multilayer_mirror(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& orbit);

// Skipped unusable routine track1_postprocess_def:
// Routine in configuration skip list

// Skipped unusable routine track1_preprocess_def:
// Routine in configuration skip list
extern "C" void fortran_track1_radiation(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Int& edge /* 0D_NOT_integer */);
void track1_radiation(CoordProxy& orbit, EleProxy& ele, int edge);
extern "C" void fortran_track1_radiation_center(
    void* orbit /* 0D_NOT_type */,
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool* rad_damp /* 0D_NOT_logical */,
    c_Bool* rad_fluct /* 0D_NOT_logical */);
void track1_radiation_center(
    CoordProxy& orbit,
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> rad_damp = std::nullopt,
    std::optional<bool> rad_fluct = std::nullopt);
extern "C" void fortran_track1_runge_kutta(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
struct Track1RungeKutta {
  bool err_flag;
  TrackProxy track;
};
Track1RungeKutta track1_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_sample(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */);
void track1_sample(EleProxy& ele, LatParamProxy& param, CoordProxy& orbit);
extern "C" void fortran_track1_spin(
    void* start_orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    c_Bool* make_quaternion /* 0D_NOT_logical */);
struct Track1Spin {
  EleProxy ele;
  CoordProxy end_orb;
};
Track1Spin track1_spin(
    CoordProxy& start_orb,
    LatParamProxy& param,
    optional_ref<bool> make_quaternion = std::nullopt);

// Skipped unusable routine track1_spin_custom_def:
// Routine in configuration skip list
extern "C" void fortran_track1_spin_integration(
    void* start_orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */);
CoordProxy track1_spin_integration(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_spin_taylor(
    void* start_orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */);
CoordProxy track1_spin_taylor(
    CoordProxy& start_orb,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_sr_wake(
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */);
void track1_sr_wake(BunchProxy& bunch, EleProxy& ele);
extern "C" void fortran_track1_symp_lie_ptc(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* track /* 0D_NOT_type */);
TrackProxy track1_symp_lie_ptc(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track1_taylor(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* taylor /* 1D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track1_taylor(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<FixedArray1D<TaylorProxy, 6>> taylor = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track1_time_runge_kutta(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */,
    c_Real* t_end /* 0D_NOT_real */,
    c_Real* dt_step /* 0D_NOT_real */);
struct Track1TimeRungeKutta {
  bool err_flag;
  TrackProxy track;
};
Track1TimeRungeKutta track1_time_runge_kutta(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<double> t_end = std::nullopt,
    optional_ref<double> dt_step = std::nullopt);

// Skipped unusable routine track1_wake_hook_def:
// Routine in configuration skip list
extern "C" void fortran_track_a_beambeam(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    void* track /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
struct TrackABeambeam {
  TrackProxy track;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
TrackABeambeam track_a_beambeam(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_bend(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_bend(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_bend_photon(
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& length /* 0D_NOT_real */);
void track_a_bend_photon(CoordProxy& orb, EleProxy& ele, double length);
extern "C" void fortran_track_a_capillary(
    void* orb /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */);
void track_a_capillary(CoordProxy& orb, EleProxy& ele);
extern "C" void fortran_track_a_converter(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_converter(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_crab_cavity(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_crab_cavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_drift(
    void* orb /* 0D_NOT_type */,
    c_Real& length /* 0D_NOT_real */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */,
    c_Int* ele_orientation /* 0D_NOT_integer */,
    c_Bool* include_ref_motion /* 0D_NOT_logical */,
    c_Real* time /* 0D_NOT_real */);
void track_a_drift(
    CoordProxy& orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<int> ele_orientation = std::nullopt,
    std::optional<bool> include_ref_motion = std::nullopt,
    optional_ref<double> time = std::nullopt);
extern "C" void fortran_track_a_drift_photon(
    void* orb /* 0D_NOT_type */,
    c_Real& length /* 0D_NOT_real */,
    c_Bool& phase_relative_to_ref /* 0D_NOT_logical */);
void track_a_drift_photon(
    CoordProxy& orb,
    double length,
    bool phase_relative_to_ref);
extern "C" void fortran_track_a_foil(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_foil(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_gkicker(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_gkicker(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_lcavity(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_lcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_lcavity_old(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_lcavity_old(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_mask(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_mask(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_match(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool* err_flag /* 0D_NOT_logical */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_match(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_patch(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool* drift_to_exit /* 0D_NOT_logical */,
    c_Real& s_ent /* 0D_NOT_real */,
    c_Real& ds_ref /* 0D_NOT_real */,
    c_Bool* track_spin /* 0D_NOT_logical */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
struct TrackAPatch {
  double s_ent;
  double ds_ref;
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
};
TrackAPatch track_a_patch(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> track_spin = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_patch_photon(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool* drift_to_exit /* 0D_NOT_logical */,
    c_Bool* use_z_pos /* 0D_NOT_logical */);
void track_a_patch_photon(
    EleProxy& ele,
    CoordProxy& orbit,
    std::optional<bool> drift_to_exit = std::nullopt,
    std::optional<bool> use_z_pos = std::nullopt);
extern "C" void fortran_track_a_pickup(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool* err_flag /* 0D_NOT_logical */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_pickup(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    optional_ref<bool> err_flag = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_quadrupole(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_quadrupole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_rfcavity(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_rfcavity(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_sad_mult(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_sad_mult(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_sol_quad(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_sol_quad(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_thick_multipole(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
void track_a_thick_multipole(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_wiggler(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr mat6 /* 2D_NOT_real */,
    c_Bool* make_matrix /* 0D_NOT_logical */);
FixedArray2D<Real, 6, 6> track_a_wiggler(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param,
    std::optional<bool> make_matrix = std::nullopt);
extern "C" void fortran_track_a_zero_length_element(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* track /* 0D_NOT_type */);
struct TrackAZeroLengthElement {
  bool err_flag;
  TrackProxy track;
};
TrackAZeroLengthElement track_a_zero_length_element(
    CoordProxy& orbit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_track_all(
    void* lat /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Int& track_state /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* orbit0 /* 1D_ALLOC_type */,
    c_Bool* init_lost /* 0D_NOT_logical */);
struct TrackAll {
  int track_state;
  bool err_flag;
  CoordProxyAlloc1D orbit0;
};
TrackAll track_all(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> init_lost = std::nullopt);
extern "C" void fortran_track_beam(
    void* lat /* 0D_NOT_type */,
    void* beam /* 0D_NOT_type */,
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    void* centroid /* 1D_ALLOC_type */,
    c_Int* direction /* 0D_NOT_integer */,
    void* bunch_tracks /* 1D_ALLOC_type */);
bool track_beam(
    LatProxy& lat,
    BeamProxy& beam,
    optional_ref<EleProxy> ele1 = std::nullopt,
    optional_ref<EleProxy> ele2 = std::nullopt,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxyAlloc1D> bunch_tracks = std::nullopt);
extern "C" void fortran_track_bunch(
    void* lat /* 0D_NOT_type */,
    void* bunch /* 0D_NOT_type */,
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    void* centroid /* 1D_ALLOC_type */,
    c_Int* direction /* 0D_NOT_integer */,
    void* bunch_track /* 0D_NOT_type */);
bool track_bunch(
    LatProxy& lat,
    BunchProxy& bunch,
    optional_ref<EleProxy> ele1 = std::nullopt,
    optional_ref<EleProxy> ele2 = std::nullopt,
    optional_ref<CoordProxyAlloc1D> centroid = std::nullopt,
    std::optional<int> direction = std::nullopt,
    optional_ref<BunchTrackProxy> bunch_track = std::nullopt);
extern "C" void fortran_track_bunch_time(
    void* bunch /* 0D_NOT_type */,
    void* branch /* 0D_NOT_type */,
    c_Real& t_end /* 0D_NOT_real */,
    c_Real& s_end /* 0D_NOT_real */,
    void* dt_step /* 1D_ALLOC_real */,
    void* extra_field /* 1D_ALLOC_type */);
void track_bunch_time(
    BunchProxy& bunch,
    BranchProxy& branch,
    double t_end,
    double s_end,
    optional_ref<RealAlloc1D> dt_step = std::nullopt,
    optional_ref<EmFieldProxyAlloc1D> extra_field = std::nullopt);
extern "C" void fortran_track_bunch_to_s(
    void* bunch /* 0D_NOT_type */,
    c_Real& s /* 0D_NOT_real */,
    void* branch /* 0D_NOT_type */);
void track_bunch_to_s(BunchProxy& bunch, double s, BranchProxy& branch);
extern "C" void fortran_track_bunch_to_t(
    void* bunch /* 0D_NOT_type */,
    c_Real& t_target /* 0D_NOT_real */,
    void* branch /* 0D_NOT_type */);
void track_bunch_to_t(BunchProxy& bunch, double t_target, BranchProxy& branch);
extern "C" void fortran_track_complex_taylor(
    void* start_orb /* 1D_ALLOC_complex */,
    void* complex_taylor /* 1D_ALLOC_type */,
    void* end_orb /* 1D_ALLOC_complex */);
ComplexAlloc1D track_complex_taylor(
    ComplexAlloc1D& start_orb,
    ComplexTaylorProxyAlloc1D& complex_taylor);
extern "C" void fortran_track_from_s_to_s(
    void* lat /* 0D_NOT_type */,
    c_Real& s_start /* 0D_NOT_real */,
    c_Real& s_end /* 0D_NOT_real */,
    void* orbit_start /* 0D_NOT_type */,
    void* orbit_end /* 0D_NOT_type */,
    void* all_orb /* 1D_ALLOC_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Int& track_state /* 0D_NOT_integer */,
    c_Int* ix_ele_end /* 0D_NOT_integer */);
struct TrackFromSToS {
  CoordProxy orbit_end;
  CoordProxyAlloc1D all_orb;
  int track_state;
};
TrackFromSToS track_from_s_to_s(
    LatProxy& lat,
    double s_start,
    double s_end,
    CoordProxy& orbit_start,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ix_ele_end = std::nullopt);
extern "C" void fortran_track_many(
    void* lat /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */,
    c_Int& ix_start /* 0D_NOT_integer */,
    c_Int& ix_end /* 0D_NOT_integer */,
    c_Int& direction /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Int& track_state /* 0D_NOT_integer */);
int track_many(
    LatProxy& lat,
    CoordProxyAlloc1D& orbit,
    int ix_start,
    int ix_end,
    int direction,
    std::optional<int> ix_branch = std::nullopt);

// Skipped unusable routine track_many_hook_def:
// Routine in configuration skip list

// Skipped unusable routine track_point_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine track_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_track_to_surface(
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_RealArr w_surface /* 2D_NOT_real */);
FixedArray2D<Real, 3, 3> track_to_surface(
    EleProxy& ele,
    CoordProxy& orbit,
    LatParamProxy& param);
extern "C" void fortran_track_until_dead(
    void* start_orb /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */,
    void* end_orb /* 0D_NOT_type */,
    void* track /* 0D_NOT_type */);
struct TrackUntilDead {
  CoordProxy end_orb;
  TrackProxy track;
};
TrackUntilDead track_until_dead(CoordProxy& start_orb, LatProxy& lat);
extern "C" void fortran_tracking_rad_map_setup(
    void* ele /* 0D_NOT_type */,
    c_Real& tollerance /* 0D_NOT_real */,
    c_Int& ref_edge /* 0D_NOT_integer */,
    void* rad_map /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct TrackingRadMapSetup {
  RadMapProxy rad_map;
  bool err_flag;
};
TrackingRadMapSetup tracking_rad_map_setup(
    EleProxy& ele,
    double tollerance,
    int ref_edge);
extern "C" void fortran_transfer_ac_kick(
    void* ac_in /* 0D_PTR_type */,
    void* ac_out /* 0D_PTR_type */);
AcKickerProxy transfer_ac_kick(AcKickerProxy& ac_in);
extern "C" void fortran_transfer_branch(
    void* branch1 /* 0D_NOT_type */,
    void* branch2 /* 0D_NOT_type */);
BranchProxy transfer_branch(BranchProxy& branch1);
extern "C" void fortran_transfer_branch_parameters(
    void* branch_in /* 0D_NOT_type */,
    void* branch_out /* 0D_NOT_type */);
BranchProxy transfer_branch_parameters(BranchProxy& branch_in);
extern "C" void fortran_transfer_branches(
    void* branch1 /* 1D_ALLOC_type */,
    void* branch2 /* 1D_ALLOC_type */);
BranchProxyAlloc1D transfer_branches(BranchProxyAlloc1D& branch1);
extern "C" void fortran_transfer_ele(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool* nullify_pointers /* 0D_NOT_logical */);
EleProxy transfer_ele(
    EleProxy& ele1,
    std::optional<bool> nullify_pointers = std::nullopt);
extern "C" void fortran_transfer_ele_taylor(
    void* ele_in /* 0D_NOT_type */,
    void* ele_out /* 0D_NOT_type */,
    c_Int* taylor_order /* 0D_NOT_integer */);
EleProxy transfer_ele_taylor(
    EleProxy& ele_in,
    std::optional<int> taylor_order = std::nullopt);
extern "C" void fortran_transfer_eles(
    void* ele1 /* 1D_ALLOC_type */,
    void* ele2 /* 1D_ALLOC_type */);
EleProxyAlloc1D transfer_eles(EleProxyAlloc1D& ele1);
extern "C" void fortran_transfer_fieldmap(
    void* ele_in /* 0D_NOT_type */,
    void* ele_out /* 0D_NOT_type */,
    c_Int& who /* 0D_NOT_integer */);
EleProxy transfer_fieldmap(EleProxy& ele_in, int who);
extern "C" bool fortran_transfer_fixer_params(
    void* fixer /* 0D_NOT_type */,
    c_Bool& to_stored /* 0D_NOT_logical */,
    void* orbit /* 0D_NOT_type */,
    c_Char who /* 0D_NOT_character */,
    c_Bool& is_ok /* 0D_NOT_logical */);
bool transfer_fixer_params(
    EleProxy& fixer,
    bool to_stored,
    optional_ref<CoordProxy> orbit = std::nullopt,
    std::optional<std::string> who = std::nullopt);
extern "C" void fortran_transfer_lat(
    void* lat1 /* 0D_NOT_type */,
    void* lat2 /* 0D_NOT_type */);
LatProxy transfer_lat(LatProxy& lat1);
extern "C" void fortran_transfer_lat_parameters(
    void* lat_in /* 0D_NOT_type */,
    void* lat_out /* 0D_NOT_type */);
LatProxy transfer_lat_parameters(LatProxy& lat_in);
extern "C" void fortran_transfer_map_calc(
    void* lat /* 0D_NOT_type */,
    void* orb_map /* 1D_ALLOC_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Int* ix1 /* 0D_NOT_integer */,
    c_Int* ix2 /* 0D_NOT_integer */,
    void* ref_orb /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* one_turn /* 0D_NOT_logical */,
    c_Bool* unit_start /* 0D_NOT_logical */,
    c_Bool* concat_if_possible /* 0D_NOT_logical */,
    void* spin_map /* 1D_ALLOC_type */);
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
    void* lat /* 0D_NOT_type */,
    void* t_map /* 1D_ALLOC_type */,
    c_Real* s1 /* 0D_NOT_real */,
    c_Real* s2 /* 0D_NOT_real */,
    void* ref_orb_in /* 0D_NOT_type */,
    void* ref_orb_out /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* one_turn /* 0D_NOT_logical */,
    c_Bool* unit_start /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* concat_if_possible /* 0D_NOT_logical */,
    void* spin_map /* 1D_ALLOC_type */);
struct TransferMapFromSToS {
  CoordProxy ref_orb_out;
  bool err_flag;
};
TransferMapFromSToS transfer_map_from_s_to_s(
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
    void* twiss1 /* 0D_NOT_type */,
    void* twiss2 /* 0D_NOT_type */,
    c_RealArr mat /* 2D_NOT_real */);
FixedArray2D<Real, 2, 2> transfer_mat2_from_twiss(
    TwissProxy& twiss1,
    TwissProxy& twiss2);
extern "C" void fortran_transfer_mat_from_twiss(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_RealArr orb1 /* 1D_NOT_real */,
    c_RealArr orb2 /* 1D_NOT_real */,
    c_RealArr m /* 2D_NOT_real */);
FixedArray2D<Real, 6, 6> transfer_mat_from_twiss(
    EleProxy& ele1,
    EleProxy& ele2,
    FixedArray1D<Real, 6> orb1,
    FixedArray1D<Real, 6> orb2);
extern "C" void fortran_transfer_matrix_calc(
    void* lat /* 0D_NOT_type */,
    c_RealArr xfer_mat /* 2D_NOT_real */,
    c_RealArr xfer_vec /* 1D_NOT_real */,
    c_Int* ix1 /* 0D_NOT_integer */,
    c_Int* ix2 /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* one_turn /* 0D_NOT_logical */);
void transfer_matrix_calc(
    LatProxy& lat,
    FixedArray2D<Real, 6, 6> xfer_mat,
    std::optional<FixedArray1D<Real, 6>> xfer_vec = std::nullopt,
    std::optional<int> ix1 = std::nullopt,
    std::optional<int> ix2 = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> one_turn = std::nullopt);
extern "C" void fortran_transfer_twiss(
    void* ele_in /* 0D_NOT_type */,
    void* ele_out /* 0D_NOT_type */,
    c_Bool* reverse /* 0D_NOT_logical */);
EleProxy transfer_twiss(
    EleProxy& ele_in,
    std::optional<bool> reverse = std::nullopt);
extern "C" void fortran_transfer_wake(
    void* wake_in /* 0D_PTR_type */,
    void* wake_out /* 0D_PTR_type */);
WakeProxy transfer_wake(WakeProxy& wake_in);

// Skipped unusable routine transfer_wall3d:
// Routine in configuration skip list

// Skipped unusable routine transport_with_sr:
// Variable in sized array: M(:,:) 2D_NOT_real
// Variable in sized array: Bone(:,:) 2D_NOT_real
// Variable inout sized array: Yone(:,:) 2D_NOT_real

// Skipped unusable routine transport_with_sr_and_ibs:
// Variable in sized array: M(:,:) 2D_NOT_real
// Variable in sized array: Bone(:,:) 2D_NOT_real
// Variable in sized array: Yone(:,:) 2D_NOT_real
extern "C" void fortran_truncate_complex_taylor_to_order(
    void* complex_taylor_in /* 1D_ALLOC_type */,
    c_Int& order /* 0D_NOT_integer */,
    void* complex_taylor_out /* 1D_ALLOC_type */);
ComplexTaylorProxyAlloc1D truncate_complex_taylor_to_order(
    ComplexTaylorProxyAlloc1D& complex_taylor_in,
    int order);
extern "C" void fortran_twiss1_propagate(
    void* twiss1 /* 0D_NOT_type */,
    c_RealArr mat2 /* 2D_NOT_real */,
    c_Int& ele_key /* 0D_NOT_integer */,
    c_Real& length /* 0D_NOT_real */,
    void* twiss2 /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
struct Twiss1Propagate {
  TwissProxy twiss2;
  bool err;
};
Twiss1Propagate twiss1_propagate(
    TwissProxy& twiss1,
    FixedArray2D<Real, 2, 2> mat2,
    int ele_key,
    double length);
extern "C" void fortran_twiss3_at_start(
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_RealArr tune3 /* 1D_NOT_real */);
FixedArray1D<Real, 3> twiss3_at_start(
    LatProxy& lat,
    bool err_flag,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss3_from_twiss2(void* ele /* 0D_NOT_type */);
void twiss3_from_twiss2(EleProxy& ele);
extern "C" void fortran_twiss3_propagate1(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void twiss3_propagate1(EleProxy& ele1, EleProxy& ele2, bool err_flag);
extern "C" void fortran_twiss3_propagate_all(
    void* lat /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */);
void twiss3_propagate_all(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_twiss_and_track_at_s(
    void* lat /* 0D_NOT_type */,
    c_Real& s /* 0D_NOT_real */,
    void* ele_at_s /* 0D_NOT_type */,
    void* orb /* 1D_ALLOC_type */,
    void* orb_at_s /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* use_last /* 0D_NOT_logical */,
    c_Bool* compute_floor_coords /* 0D_NOT_logical */);
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
    void* branch /* 0D_NOT_type */,
    void* orbit_start /* 0D_NOT_type */,
    c_Real& s_end /* 0D_NOT_real */,
    void* orbit_end /* 0D_NOT_type */,
    void* ele_start /* 0D_NOT_type */,
    void* ele_end /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* compute_floor_coords /* 0D_NOT_logical */,
    c_Bool* compute_twiss /* 0D_NOT_logical */);
struct TwissAndTrackFromSToS {
  CoordProxy orbit_end;
  EleProxy ele_end;
  bool err;
};
TwissAndTrackFromSToS twiss_and_track_from_s_to_s(
    BranchProxy& branch,
    CoordProxy& orbit_start,
    double s_end,
    optional_ref<EleProxy> ele_start = std::nullopt,
    std::optional<bool> compute_floor_coords = std::nullopt,
    std::optional<bool> compute_twiss = std::nullopt);
extern "C" void fortran_twiss_and_track_intra_ele(
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& l_start /* 0D_NOT_real */,
    c_Real& l_end /* 0D_NOT_real */,
    c_Bool& track_upstream_end /* 0D_NOT_logical */,
    c_Bool& track_downstream_end /* 0D_NOT_logical */,
    void* orbit_start /* 0D_NOT_type */,
    void* orbit_end /* 0D_NOT_type */,
    void* ele_start /* 0D_NOT_type */,
    void* ele_end /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* compute_floor_coords /* 0D_NOT_logical */,
    c_Bool* compute_twiss /* 0D_NOT_logical */,
    c_Bool* reuse_ele_end /* 0D_NOT_logical */);
struct TwissAndTrackIntraEle {
  CoordProxy orbit_end;
  bool err;
};
TwissAndTrackIntraEle twiss_and_track_intra_ele(
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
    void* ele /* 0D_NOT_type */,
    void* start /* 0D_NOT_type */,
    void* end /* 0D_NOT_type */,
    void* average /* 0D_NOT_type */);
struct TwissAtElement {
  EleProxy start;
  EleProxy end;
  EleProxy average;
};
TwissAtElement twiss_at_element(EleProxy& ele);
extern "C" void fortran_twiss_at_start(
    void* lat /* 0D_NOT_type */,
    c_Int& status /* 0D_NOT_integer */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* type_out /* 0D_NOT_logical */);
int twiss_at_start(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> type_out = std::nullopt);

// Skipped unusable routine twiss_from_mat2:
// Variable inout sized array: mat_in(:,:) 2D_NOT_real

// Skipped unusable routine twiss_from_mat6:
// Variable in sized array: mat6(:,:) 2D_NOT_real
extern "C" void fortran_twiss_from_tracking(
    void* lat /* 0D_NOT_type */,
    void* ref_orb0 /* 0D_NOT_type */,
    c_Real& symp_err /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* d_orb /* 1D_ALLOC_real */);
struct TwissFromTracking {
  double symp_err;
  bool err_flag;
};
TwissFromTracking twiss_from_tracking(
    LatProxy& lat,
    CoordProxy& ref_orb0,
    optional_ref<RealAlloc1D> d_orb = std::nullopt);
extern "C" void fortran_twiss_propagate1(
    void* ele1 /* 0D_NOT_type */,
    void* ele2 /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* forward /* 0D_NOT_logical */);
bool twiss_propagate1(
    EleProxy& ele1,
    EleProxy& ele2,
    std::optional<bool> forward = std::nullopt);
extern "C" void fortran_twiss_propagate_all(
    void* lat /* 0D_NOT_type */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Int* ie_start /* 0D_NOT_integer */,
    c_Int* ie_end /* 0D_NOT_integer */);
bool twiss_propagate_all(
    LatProxy& lat,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<int> ie_start = std::nullopt,
    std::optional<int> ie_end = std::nullopt);

// Skipped unusable routine twiss_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_twiss_to_1_turn_mat(
    void* twiss /* 0D_NOT_type */,
    c_Real& phi /* 0D_NOT_real */,
    c_RealArr mat2 /* 2D_NOT_real */);
FixedArray2D<Real, 2, 2> twiss_to_1_turn_mat(TwissProxy& twiss, double phi);

// Skipped unusable routine type_complex_taylors:
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_type_coord(void* coord /* 0D_NOT_type */);
void type_coord(CoordProxy& coord);

// Skipped unusable routine type_ele:
// Variable-sized in character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_end_stuff:
// Variable-sized inout character array: li(:) 1D_ALLOC_character
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_expression_tree:
// Untranslated type: ExpressionTreeProxy (0D_NOT_type)

// Skipped unusable routine type_map:
// Untranslated type: Real8Proxy (1D_ALLOC_type)

// Skipped unusable routine type_map1:
// Untranslated type: Real8Proxy (1D_ALLOC_type)

// Skipped unusable routine type_ptc_fibre:
// Untranslated type: FibreProxy (0D_PTR_type)
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_ptc_internal_state:
// Untranslated type: InternalStateProxy (0D_NOT_type)
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_ptc_layout:
// Untranslated type: LayoutProxy (0D_NOT_type)

// Skipped unusable routine type_real_8_taylors:
// Untranslated type: Real8Proxy (1D_ALLOC_type)

// Skipped unusable routine type_taylors:
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine type_twiss:
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine universal_equal_universal:
// Untranslated type: UniversalTaylorProxy (0D_NOT_type)
// Untranslated type: UniversalTaylorProxy (0D_NOT_type)

// Skipped unusable routine universal_to_bmad_taylor:
// Untranslated type: UniversalTaylorProxy (0D_NOT_type)

// Skipped unusable routine unlink_fieldmap:
// Routine in configuration skip list

// Skipped unusable routine unlink_wall3d:
// Routine in configuration skip list
extern "C" void fortran_update_ele_from_fibre(void* ele /* 0D_NOT_type */);
void update_ele_from_fibre(EleProxy& ele);
extern "C" void fortran_update_fibre_from_ele(
    void* ele /* 0D_NOT_type */,
    c_Bool& survey_needed /* 0D_NOT_logical */);
bool update_fibre_from_ele(EleProxy& ele);
extern "C" void fortran_update_floor_angles(
    void* floor /* 0D_NOT_type */,
    void* floor0 /* 0D_NOT_type */);
void update_floor_angles(
    FloorPositionProxy& floor,
    optional_ref<FloorPositionProxy> floor0 = std::nullopt);
extern "C" bool fortran_valid_field_calc(
    void* ele /* 0D_NOT_type */,
    c_Int& field_calc /* 0D_NOT_integer */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void valid_field_calc(EleProxy& ele, int field_calc, bool is_valid);
extern "C" bool fortran_valid_fringe_type(
    void* ele /* 0D_NOT_type */,
    c_Int& fringe_type /* 0D_NOT_integer */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void valid_fringe_type(EleProxy& ele, int fringe_type, bool is_valid);
extern "C" bool fortran_valid_mat6_calc_method(
    void* ele /* 0D_NOT_type */,
    c_Int& species /* 0D_NOT_integer */,
    c_Int& mat6_calc_method /* 0D_NOT_integer */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void valid_mat6_calc_method(
    EleProxy& ele,
    int species,
    int mat6_calc_method,
    bool is_valid);
extern "C" bool fortran_valid_spin_tracking_method(
    void* ele /* 0D_NOT_type */,
    c_Int& spin_tracking_method /* 0D_NOT_integer */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void valid_spin_tracking_method(
    EleProxy& ele,
    int spin_tracking_method,
    bool is_valid);
extern "C" bool fortran_valid_tracking_method(
    void* ele /* 0D_NOT_type */,
    c_Int& species /* 0D_NOT_integer */,
    c_Int& tracking_method /* 0D_NOT_integer */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void valid_tracking_method(
    EleProxy& ele,
    int species,
    int tracking_method,
    bool is_valid);
extern "C" bool fortran_value_of_attribute(
    void* ele /* 0D_NOT_type */,
    c_Char attrib_name /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Real* err_value /* 0D_NOT_real */,
    c_Real& value /* 0D_NOT_real */);
bool value_of_attribute(
    EleProxy& ele,
    std::string attrib_name,
    std::optional<bool> err_print_flag,
    std::optional<double> err_value,
    double value);
extern "C" void fortran_value_to_line(
    c_Char line /* 0D_NOT_character */,
    c_Real& value /* 0D_NOT_real */,
    c_Char str /* 0D_NOT_character */,
    c_Char typ /* 0D_NOT_character */,
    c_Bool* ignore_if_zero /* 0D_NOT_logical */,
    c_Bool* use_comma /* 0D_NOT_logical */);
void value_to_line(
    std::string line,
    double value,
    std::string str,
    std::string typ,
    optional_ref<bool> ignore_if_zero = std::nullopt,
    optional_ref<bool> use_comma = std::nullopt);
extern "C" bool fortran_vec_to_polar(
    c_RealArr vec /* 1D_NOT_real */,
    c_Real* phase /* 0D_NOT_real */,
    void* polar /* 0D_NOT_type */);
void vec_to_polar(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    SpinPolarProxy& polar);
extern "C" bool fortran_vec_to_spinor(
    c_RealArr vec /* 1D_NOT_real */,
    c_Real* phase /* 0D_NOT_real */,
    c_ComplexArr spinor /* 1D_NOT_complex */);
void vec_to_spinor(
    FixedArray1D<Real, 3> vec,
    std::optional<double> phase,
    FixedArray1D<Complex, 2> spinor);
extern "C" bool fortran_verify_valid_name(
    c_Char name /* 0D_NOT_character */,
    c_Int& ix_name /* 0D_NOT_integer */,
    c_Bool* pure_name /* 0D_NOT_logical */,
    c_Bool* include_wild /* 0D_NOT_logical */,
    c_Bool& is_valid /* 0D_NOT_logical */);
bool verify_valid_name(
    std::string name,
    int ix_name,
    std::optional<bool> pure_name = std::nullopt,
    std::optional<bool> include_wild = std::nullopt);
extern "C" bool fortran_w_mat_for_bend_angle(
    c_Real& angle /* 0D_NOT_real */,
    c_Real& ref_tilt /* 0D_NOT_real */,
    c_RealArr r_vec /* 1D_NOT_real */,
    c_RealArr w_mat /* 2D_NOT_real */);
void w_mat_for_bend_angle(
    double angle,
    double ref_tilt,
    std::optional<FixedArray1D<Real, 3>> r_vec,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_tilt(
    c_Real& tilt /* 0D_NOT_real */,
    c_Bool* return_inverse /* 0D_NOT_logical */,
    c_RealArr w_mat /* 2D_NOT_real */);
void w_mat_for_tilt(
    double tilt,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_x_pitch(
    c_Real& x_pitch /* 0D_NOT_real */,
    c_Bool* return_inverse /* 0D_NOT_logical */,
    c_RealArr w_mat /* 2D_NOT_real */);
void w_mat_for_x_pitch(
    double x_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_for_y_pitch(
    c_Real& y_pitch /* 0D_NOT_real */,
    c_Bool* return_inverse /* 0D_NOT_logical */,
    c_RealArr w_mat /* 2D_NOT_real */);
void w_mat_for_y_pitch(
    double y_pitch,
    std::optional<bool> return_inverse,
    FixedArray2D<Real, 3, 3> w_mat);

// Skipped unusable routine wake_lr_mode_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wake_lr_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wake_sr_mode_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wake_sr_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wake_sr_z_long_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wake_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_wall3d_d_radius(
    void* position /* 1D_ALLOC_real */,
    void* ele /* 0D_NOT_type */,
    c_Int* ix_wall /* 0D_NOT_integer */,
    c_RealArr perp /* 1D_NOT_real */,
    c_Int& ix_section /* 0D_NOT_integer */,
    c_Bool& no_wall_here /* 0D_NOT_logical */,
    c_RealArr origin /* 1D_NOT_real */,
    c_Real& radius_wall /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Real& d_radius /* 0D_NOT_real */);
struct Wall3dDRadius {
  FixedArray1D<Real, 3> perp;
  int ix_section;
  bool no_wall_here;
  FixedArray1D<Real, 3> origin;
  double radius_wall;
  bool err_flag;
  double d_radius;
};
Wall3dDRadius wall3d_d_radius(
    RealAlloc1D& position,
    EleProxy& ele,
    std::optional<int> ix_wall = std::nullopt);
extern "C" void fortran_wall3d_initializer(
    void* wall3d /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
bool wall3d_initializer(Wall3dProxy& wall3d);
extern "C" void fortran_wall3d_section_initializer(
    void* section /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
bool wall3d_section_initializer(Wall3dSectionProxy& section);

// Skipped unusable routine wall3d_section_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wall3d_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" bool fortran_wall3d_to_position(
    void* orbit /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_RealArr position /* 1D_NOT_real */);
FixedArray1D<Real, 6> wall3d_to_position(CoordProxy& orbit, EleProxy& ele);

// Skipped unusable routine wall3d_vertex_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine wall_hit_handler_custom_def:
// Routine in configuration skip list

// Skipped unusable routine wiggler_modeling_common_struct_to_json:
// Routine module (bmad_json) in configuration skip list
extern "C" void fortran_word_to_value(
    c_Char word /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Real& value /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */);
void word_to_value(
    std::string word,
    LatProxy& lat,
    double value,
    bool err_flag,
    optional_ref<EleProxy> ele = std::nullopt);

// Skipped unusable routine write_2d:
// Variable inout sized array: grid(:,:) 2D_NOT_real
extern "C" void fortran_write_ascii_beam_file(
    c_Char file_name /* 0D_NOT_character */,
    void* beam /* 0D_NOT_type */,
    c_Bool* new_file /* 0D_NOT_logical */,
    c_Bool* alive_only /* 0D_NOT_logical */);
void write_ascii_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<bool> alive_only = std::nullopt);
extern "C" void fortran_write_astra_bend(
    c_Int& iu /* 0D_NOT_integer */,
    c_Real& strength /* 0D_NOT_real */,
    c_Int& id /* 0D_NOT_integer */,
    c_RealArr d1 /* 1D_NOT_real */,
    c_RealArr d2 /* 1D_NOT_real */,
    c_RealArr d3 /* 1D_NOT_real */,
    c_RealArr d4 /* 1D_NOT_real */);
void write_astra_bend(
    int iu,
    double strength,
    int id,
    FixedArray1D<Real, 2> d1,
    FixedArray1D<Real, 2> d2,
    FixedArray1D<Real, 2> d3,
    FixedArray1D<Real, 2> d4);

// Skipped unusable routine write_astra_ele:
// Untranslated type: StrIndexProxy (0D_NOT_type)
extern "C" void fortran_write_astra_field_grid_file(
    c_Int& astra_file_unit /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Real* dz /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteAstraFieldGridFile {
  double maxfield;
  bool err;
};
WriteAstraFieldGridFile write_astra_field_grid_file(
    int astra_file_unit,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);
extern "C" void fortran_write_astra_field_grid_file_3d(
    c_Char base_filename /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Real* dz /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteAstraFieldGridFile3d {
  double maxfield;
  bool err;
};
WriteAstraFieldGridFile3d write_astra_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);

// Skipped unusable routine write_astra_lattice_file:
// Untranslated type: AstraLatticeParamProxy (0D_NOT_type)
extern "C" void fortran_write_beam_file(
    c_Char file_name /* 0D_NOT_character */,
    void* beam /* 0D_NOT_type */,
    c_Bool* new_file /* 0D_NOT_logical */,
    c_Int* file_format /* 0D_NOT_integer */,
    void* lat /* 0D_NOT_type */,
    c_Bool* alive_only /* 0D_NOT_logical */);
void write_beam_file(
    std::string file_name,
    BeamProxy& beam,
    std::optional<bool> new_file = std::nullopt,
    std::optional<int> file_format = std::nullopt,
    optional_ref<LatProxy> lat = std::nullopt,
    std::optional<bool> alive_only = std::nullopt);
extern "C" void fortran_write_beam_floor_positions(
    c_Char file_name /* 0D_NOT_character */,
    void* beam /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool* new_file /* 0D_NOT_logical */);
void write_beam_floor_positions(
    std::string file_name,
    BeamProxy& beam,
    EleProxy& ele,
    std::optional<bool> new_file = std::nullopt);
extern "C" void fortran_write_binary_cartesian_map(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* cart_map /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void write_binary_cartesian_map(
    std::string file_name,
    EleProxy& ele,
    CartesianMapProxy& cart_map,
    bool err_flag);
extern "C" void fortran_write_binary_cylindrical_map(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* cl_map /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void write_binary_cylindrical_map(
    std::string file_name,
    EleProxy& ele,
    CylindricalMapProxy& cl_map,
    bool err_flag);
extern "C" void fortran_write_binary_grid_field(
    c_Char file_name /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    void* g_field /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void write_binary_grid_field(
    std::string file_name,
    EleProxy& ele,
    GridFieldProxy& g_field,
    bool err_flag);
extern "C" void fortran_write_blender_ele(
    c_Int& iu /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Bool* old_format /* 0D_NOT_logical */);
void write_blender_ele(
    int iu,
    EleProxy& ele,
    optional_ref<bool> old_format = std::nullopt);
extern "C" void fortran_write_blender_lat_layout(
    c_Char file_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */);
void write_blender_lat_layout(std::string file_name, LatProxy& lat);
extern "C" void fortran_write_bmad_lattice_file(
    c_Char bmad_file /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Int* output_form /* 0D_NOT_integer */,
    void* orbit0 /* 0D_NOT_type */);
bool write_bmad_lattice_file(
    std::string bmad_file,
    LatProxy& lat,
    std::optional<int> output_form = std::nullopt,
    optional_ref<CoordProxy> orbit0 = std::nullopt);

// Skipped unusable routine write_digested_bmad_file:
// Variable-sized in character array: file_names(:) 1D_ALLOC_character
// Untranslated type: ExtraParsingInfoProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine write_gpt_ele:
// Untranslated type: StrIndexProxy (0D_NOT_type)
extern "C" void fortran_write_gpt_field_grid_file_1d(
    c_Int& gpt_file_unit /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Real& ref_time /* 0D_NOT_real */,
    c_Real* dz /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteGptFieldGridFile1d {
  double maxfield;
  double ref_time;
  bool err;
};
WriteGptFieldGridFile1d write_gpt_field_grid_file_1d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);
extern "C" void fortran_write_gpt_field_grid_file_2d(
    c_Int& gpt_file_unit /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Real& ref_time /* 0D_NOT_real */,
    c_Real* dr /* 0D_NOT_real */,
    c_Real* dz /* 0D_NOT_real */,
    c_Real* r_max /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteGptFieldGridFile2d {
  double maxfield;
  double ref_time;
  bool err;
};
WriteGptFieldGridFile2d write_gpt_field_grid_file_2d(
    int gpt_file_unit,
    EleProxy& ele,
    std::optional<double> dr = std::nullopt,
    std::optional<double> dz = std::nullopt,
    std::optional<double> r_max = std::nullopt);
extern "C" void fortran_write_gpt_field_grid_file_3d(
    c_Char base_filename /* 0D_NOT_character */,
    void* ele /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Real& ref_time /* 0D_NOT_real */,
    c_Real* dz /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteGptFieldGridFile3d {
  double maxfield;
  double ref_time;
  bool err;
};
WriteGptFieldGridFile3d write_gpt_field_grid_file_3d(
    std::string base_filename,
    EleProxy& ele,
    std::optional<double> dz = std::nullopt);

// Skipped unusable routine write_gpt_lattice_file:
// Untranslated type: GptLatParamProxy (0D_NOT_type)
extern "C" void fortran_write_lat_line(
    c_Char line /* 0D_NOT_character */,
    c_Int& iu /* 0D_NOT_integer */,
    c_Bool& end_is_neigh /* 0D_NOT_logical */,
    c_Bool* do_split /* 0D_NOT_logical */,
    c_Bool* scibmad /* 0D_NOT_logical */);
void write_lat_line(
    std::string line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split = std::nullopt,
    std::optional<bool> scibmad = std::nullopt);
extern "C" void fortran_write_lattice_in_elegant_format(
    c_Char out_file_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* ref_orbit /* 1D_ALLOC_type */,
    c_Bool* use_matrix_model /* 0D_NOT_logical */,
    c_Bool* include_apertures /* 0D_NOT_logical */,
    c_Real* dr12_drift_max /* 0D_NOT_real */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err /* 0D_NOT_logical */);
bool write_lattice_in_elegant_format(
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<CoordProxyAlloc1D> ref_orbit = std::nullopt,
    std::optional<bool> use_matrix_model = std::nullopt,
    std::optional<bool> include_apertures = std::nullopt,
    std::optional<double> dr12_drift_max = std::nullopt,
    std::optional<int> ix_branch = std::nullopt);
extern "C" void fortran_write_lattice_in_foreign_format(
    c_Char out_type /* 0D_NOT_character */,
    c_Char out_file_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* ref_orbit /* 1D_ALLOC_type */,
    c_Bool* use_matrix_model /* 0D_NOT_logical */,
    c_Bool* include_apertures /* 0D_NOT_logical */,
    c_Real* dr12_drift_max /* 0D_NOT_real */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err /* 0D_NOT_logical */);
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
    c_Char out_type /* 0D_NOT_character */,
    c_Char out_file_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    void* ref_orbit /* 1D_ALLOC_type */,
    c_Bool* use_matrix_model /* 0D_NOT_logical */,
    c_Bool* include_apertures /* 0D_NOT_logical */,
    c_Real* dr12_drift_max /* 0D_NOT_real */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool& err /* 0D_NOT_logical */);
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
    c_Char out_file_name /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool* include_apertures /* 0D_NOT_logical */,
    c_Int* ix_branch /* 0D_NOT_integer */,
    c_Bool* err /* 0D_NOT_logical */);
void write_lattice_in_sad_format(
    std::string out_file_name,
    LatProxy& lat,
    optional_ref<bool> include_apertures = std::nullopt,
    optional_ref<int> ix_branch = std::nullopt,
    optional_ref<bool> err = std::nullopt);
extern "C" void fortran_write_lattice_in_scibmad(
    c_Char scibmad_file /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct WriteLatticeInScibmad {
  std::string scibmad_file;
  bool err_flag;
};
WriteLatticeInScibmad write_lattice_in_scibmad(LatProxy& lat);
extern "C" void fortran_write_line_element(
    c_Char line /* 0D_NOT_character */,
    c_Int& iu /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* lat /* 0D_NOT_type */);
void write_line_element(std::string line, int iu, EleProxy& ele, LatProxy& lat);
extern "C" void fortran_write_opal_field_grid_file(
    c_Int& opal_file_unit /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* param /* 0D_NOT_type */,
    c_Real& maxfield /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */);
struct WriteOpalFieldGridFile {
  double maxfield;
  bool err;
};
WriteOpalFieldGridFile write_opal_field_grid_file(
    int opal_file_unit,
    EleProxy& ele,
    LatParamProxy& param);
extern "C" void fortran_write_opal_lattice_file(
    c_Int& opal_file_unit /* 0D_NOT_integer */,
    void* lat /* 0D_NOT_type */,
    c_Bool& err /* 0D_NOT_logical */);
bool write_opal_lattice_file(int opal_file_unit, LatProxy& lat);
extern "C" void fortran_write_time_particle_distribution(
    c_Int& time_file_unit /* 0D_NOT_integer */,
    void* bunch /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Char style /* 0D_NOT_character */,
    void* branch /* 0D_NOT_type */,
    c_Char format /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
bool write_time_particle_distribution(
    int time_file_unit,
    BunchProxy& bunch,
    EleProxy& ele,
    std::optional<std::string> style = std::nullopt,
    optional_ref<BranchProxy> branch = std::nullopt,
    std::optional<std::string> format = std::nullopt);
extern "C" bool fortran_xlafun(
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void xlafun(double x, double y, double z, double res);
extern "C" bool fortran_xraylib_nist_compound(
    c_Char name /* 0D_NOT_character */,
    c_Int& indx /* 0D_NOT_integer */);
int xraylib_nist_compound(std::string name);

// Skipped unusable routine xrlcomplex_c_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine xsif_parser:
// Routine in configuration skip list

// Skipped unusable routine xy_disp_struct_to_json:
// Routine module (bmad_json) in configuration skip list

// Skipped unusable routine xyz_to_action:
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_ylafun(
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void ylafun(double x, double y, double z, double res);
extern "C" bool fortran_z_at_surface(
    void* ele /* 0D_NOT_type */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* extend_grid /* 0D_NOT_logical */,
    c_RealArr dz_dxy /* 1D_NOT_real */,
    c_Real& z /* 0D_NOT_real */);
struct ZAtSurface {
  bool err_flag;
  FixedArray1D<Real, 2> dz_dxy;
  double z;
};
ZAtSurface z_at_surface(
    EleProxy& ele,
    double x,
    double y,
    std::optional<bool> extend_grid = std::nullopt);
extern "C" void fortran_zero_ele_kicks(void* ele /* 0D_NOT_type */);
EleProxy zero_ele_kicks();
extern "C" void fortran_zero_ele_offsets(void* ele /* 0D_NOT_type */);
EleProxy zero_ele_offsets();
extern "C" void fortran_zero_lr_wakes_in_lat(void* lat /* 0D_NOT_type */);
void zero_lr_wakes_in_lat(LatProxy& lat);
extern "C" bool fortran_zlafun(
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& z /* 0D_NOT_real */,
    c_Real& res /* 0D_NOT_real */);
void zlafun(double x, double y, double z, double res);

// Skipped unusable routine zot_integrand:
// Untranslated type: CPtrProxy (0D_NOT_type)
} // namespace Bmad
